# -*- coding: utf-8 -*-

"""Main execution

This program searches for models that fit an specific motif.

"""

import time
import random
import copy
import json
import os
import cProfile
import pstats
import io
import numpy as np
from objects.organism_factory import OrganismFactory
from Bio import SeqIO

"""
Variable definition
"""
POPULATION_LENGTH = 0
DATASET_BASE_PATH_DIR = ""
RESULT_BASE_PATH_DIR = ""
POSITIVE_FILENAME = ""
NEGATIVE_FILENAME = ""
POPULATION_ORIGIN = ""
POPULATION_FILL_TYPE = ""
INPUT_FILENAME = ""
OUTPUT_FILENAME = ""
MAX_SEQUENCES_TO_FIT_POS = 0
MAX_SEQUENCES_TO_FIT_NEG = 0
MIN_ITERATIONS = 0
MIN_FITNESS = 0
RECOMBINATION_PROBABILITY = 0.0
THRESHOLD = 0.0

JSON_CONFIG_FILENAME = "config.json"
"""
Configuration of the object types
Populated by JSON read.
"""
configOrganism: dict = {}
configOrganismFactory: dict = {}
configConnector: dict = {}
configPssm: dict = {}


# list that holds the population
organism_population: list = []

# mean_nodes are the average number of nodes per organism in the population
# used to calculate organism complexity
mean_nodes: float = 0
# mean_fitness is the average fitness per organism in the population
# used to calculate organism complexity
mean_fitness: float = 0

# Initialize datasets
positive_dataset: list = []
negative_dataset: list = []


def main():
    """Main function for the motif seek
    """

    print("Loading parameters...")
    positive_dataset = read_fasta_file(
        DATASET_BASE_PATH_DIR + POSITIVE_FILENAME
    )
    negative_dataset = read_fasta_file(
        DATASET_BASE_PATH_DIR + NEGATIVE_FILENAME
    )

    mean_nodes = 0
    mean_recognizers = 0
    mean_fitness = 0
    """
    Generate initial population
    """
    # Instantiate organism Factory object with object configurations
    organism_factory = OrganismFactory(
        configOrganism, configOrganismFactory, configConnector, configPssm
    )
    # initialize organisms population
    organism_population = []

    # Generate population depending on origin and fill type.
    # Origin can be "random" or "file"(read a set of organisms from a file).
    if POPULATION_ORIGIN.lower() == "random":
        # For a random origin, we can generate #POPULATION_LENGTH organisms.

        for i in range(POPULATION_LENGTH):
            new_organism = organism_factory.get_organism()
            organism_population.append(new_organism)

            mean_nodes += new_organism.count_nodes()

	    # count recognizers
            new_organism.count_nodes()
            mean_recognizers += new_organism.num_recognizers

    elif POPULATION_ORIGIN.lower() == "file":
        # Set the file organisms and fill with random/same organisms
        # POPULATION_LENGTH must be >= len(fileOrganisms)
        file_organisms = organism_factory.import_organisms(INPUT_FILENAME)
        remaining_organisms = POPULATION_LENGTH - len(file_organisms)
        fill_organism_population = []

        if POPULATION_FILL_TYPE.lower() == "random":
            # FILL WITH RANDOM

            for i in range(remaining_organisms):
                new_organism = organism_factory.get_organism()
                fill_organism_population.append(new_organism)

        elif POPULATION_FILL_TYPE.lower() == "uniform":
            # FILL WITH SAME ORGANISMS IN FILE

            for i in range(remaining_organisms):
                new_organism = copy.deepcopy(
                    file_organisms[i % len(file_organisms)]
                )
                fill_organism_population.append(new_organism)
                new_organism.set_id(organism_factory.get_id())

        # join & calculate mean nodes
        organism_population = file_organisms + fill_organism_population

        for org in organism_population:
            mean_nodes += org.count_nodes()

	    # count recognizers
            org.count_nodes()
            mean_recognizers += org.num_recognizers

    else:
        raise Exception("Not a valid population origin, "
            + "check the configuration file.")
    

    # Convert node count into mean
    mean_nodes /= POPULATION_LENGTH
    mean_recognizers /= POPULATION_LENGTH
    print("len = {}".format(len(organism_population)))
    """
    Initialize iteration variables.
    """
    iterations = 0
    max_score = float("-inf")
    last_max_score = 0.0
    best_organism = (None, 0.0, 0, 0.0)
    max_organism = (None, 0.0, 0, 0.0)

    timeformat = "%Y-%m-%d--%H-%M-%S"
    print("Starting execution...")

    # Main loop, it iterates until organisms do not get a significant change
    # or MIN_ITERATIONS or MIN_FITNESS is reached.

    while not is_finished(
            END_WHILE_METHOD, iterations, max_score, last_max_score
    ):

        # Shuffle population & datasets
        # Organisms are shuffled for deterministic crowding selection
        # Datasets are shuffled for subsampling
        random.shuffle(organism_population)
        random.shuffle(negative_dataset)
        random.shuffle(positive_dataset)

        # Reset max_score
        last_max_score = max_score
        max_score = float("-inf")
        changed_best_score = False
        initial = time.time()

        a_fitness = []
        a_nodes = []

        # Deterministic crowding
        # Iterate over pairs of organisms
        for i in range(0, len(organism_population) - 1, 2):
            org1 = organism_population[i]
            org2 = organism_population[i + 1]

            # Cross parents to get children
            # Returns two children. Each child contains:
            #   - Child object itself
            #   - Similarity to organism 1
            #   - Similarity to organism 2
            #
            children = combine_organisms(org1, org2, organism_factory)

            child1 = children["child1"]["child"]
            child2 = children["child2"]["child"]

            # Mutate children
            child1.mutate(organism_factory)
            child2.mutate(organism_factory)

            # Match parent with its closest child for deterministic crowding
            # selection.
            # There are 2 possible combinations p_1-c1, p_2-c2 & p1-c2, p2-c1
            # We select a combination based on a sum of similarities in
            # combinations
            combination_1 = (
                children["child1"]["sim_org_1"]
                + children["child2"]["sim_org_2"]
            )  # Match the first parent to first child and second parent to
            # second child
            combination_2 = (
                children["child1"]["sim_org_2"]
                + children["child2"]["sim_org_1"]
            )  # Match the first parent to second child and second parent to
            # first child

            pair_children = []

            if combination_1 > combination_2:
                pair_children.append((org1, child1))
                pair_children.append((org2, child2))
            else:
                pair_children.append((org1, child2))
                pair_children.append((org2, child1))

            # Make two organisms compete
            # j index is used to re insert winning organism
            # into the population
            for j in range(len(pair_children)):

                first_organism = pair_children[j][0]  # Parent Organism
                second_organism = pair_children[j][1]  # Chid Organism
		
		
                # Boltzmannian fitness
                if FITNESS_FUNCTION == "boltzmannian":
                    performance1 = first_organism.get_boltz_fitness(positive_dataset[:MAX_SEQUENCES_TO_FIT_POS],
                                                                    negative_dataset[:MAX_SEQUENCES_TO_FIT_NEG],
                                                                    GENOME_LENGTH)
                    fitness1 = performance1["score"]
                    gini1 = performance1["avg_gini"]
                    
                    performance2 = second_organism.get_boltz_fitness(positive_dataset[:MAX_SEQUENCES_TO_FIT_POS],
                                                                     negative_dataset[:MAX_SEQUENCES_TO_FIT_NEG],
                                                                     GENOME_LENGTH)
                    fitness2 = performance2["score"]
                    gini2 = performance2["avg_gini"]
                    
                    fitness1 = round(fitness1, 8)
                    fitness2 = round(fitness2, 8)
                
                # Discriminative fitness
                elif FITNESS_FUNCTION == "discriminative":
                    positive_performance1 = first_organism.get_discriminative_fitness(positive_dataset[:MAX_SEQUENCES_TO_FIT_POS])
                    negative_performance1 = first_organism.get_discriminative_fitness(negative_dataset[:MAX_SEQUENCES_TO_FIT_NEG])
                    p_1 = positive_performance1["score"]
                    n_1 = negative_performance1["score"]
                    fitness1 =  p_1 - n_1
                    gini1 = positive_performance1["avg_gini"]
                    
                    positive_performance2 = second_organism.get_discriminative_fitness(positive_dataset[:MAX_SEQUENCES_TO_FIT_POS])
                    negative_performance2 = second_organism.get_discriminative_fitness(negative_dataset[:MAX_SEQUENCES_TO_FIT_NEG])
                    p_2 = positive_performance2["score"]
                    n_2 = negative_performance2["score"]
                    fitness2 =  p_2 - n_2
                    gini2 = positive_performance2["avg_gini"]
                
                else:
                    raise Exception("Not a valid fitness function name, "
                        + "check the configuration file.")


                if MAX_NODES != None:  # Upper_bound to complexity
                    
                    if first_organism.count_nodes() > MAX_NODES:
                        #print(first_organism.count_nodes(), "nodes")
                        fitness1 = -1000 * int(first_organism.count_nodes())
                    
                    if second_organism.count_nodes() > MAX_NODES:
                        #print(second_organism.count_nodes(), "nodes")
                        fitness2 = -1000 * int(second_organism.count_nodes())
                
                
                if INEQUALITY_PENALTY_METHOD=="avg_gini":
                    # INEQUALITY_PENALTY_PARAM acts as a penalty buffer
                    # The higher this parameter, the less important is the effect of the Gini penalty
                    # It's meant to work in the range [1, +inf)
                    effective_fitness_1 = fitness1 * (INEQUALITY_PENALTY_PARAM - gini1)
                    effective_fitness_2 = fitness2 * (INEQUALITY_PENALTY_PARAM - gini2)
                else:
                    effective_fitness_1 = fitness1
                    effective_fitness_2 = fitness2


                if (
                        effective_fitness_1 > effective_fitness_2
                ):  # The first organism wins
                    # Set it back to the population and save fitness
                    # for next iteration
                    organism_population[i + j] = first_organism
                    a_fitness.append(fitness1)
                    # If the parent wins, mean_nodes don't change
                    a_nodes.append(first_organism.count_nodes())

                    # Check if its the max score in that iteration
                    if effective_fitness_1 > max_score:
                        max_score = effective_fitness_1
                        max_organism = (
                            first_organism,
                            effective_fitness_1,
                            first_organism.count_nodes(),
                            gini1,
                        )

                    # Check if its the max score so far and if it is set it as
                    # best organism
                    if max_organism[1] > best_organism[1]:
                        # ID, EF, Nodes, Penalty applied
                        best_organism = max_organism
                        changed_best_score = True

                else:  # The second organism wins (child)
                    # Set it back to the population and save fitness for next
                    # iteration
                    organism_population[i + j] = second_organism
                    a_fitness.append(fitness2)
                    # If the child wins, update mean_nodes
                    # mean_nodes = ((meanNodes * POPULATION_LENGTH) +
                    # second_organism.count_nodes() -
                    # first_organism.count_nodes()) / POPULATION_LENGTH
                    a_nodes.append(second_organism.count_nodes())

                    # Pass tracking parameter from paretn to child
                    second_organism.set_is_tracked(first_organism.is_tracked)
                    if second_organism.is_tracked:
                        # Export it If its being tracked
                        print_ln(
                            "Evolution {}->{}".format(
                                first_organism._id, second_organism.ID
                            ),
                            RESULT_BASE_PATH_DIR + "evolution.txt",
                        )
                        filename = "tr{}_{}".format(
                            time.strftime(timeformat), second_organism._id
                        )
                        export_organism(
                            second_organism,
                            positive_dataset,
                            filename,
                            organism_factory,
                        )

                    # Check if its the max score in that iteration
                    if effective_fitness_2 > max_score:
                        max_score = effective_fitness_2
                        max_organism = (
                            second_organism,
                            effective_fitness_2,
                            second_organism.count_nodes(),
                            gini2,
                        )

                    # Check if its the max score so far and if it is set it as
                    # best organism
                    if effective_fitness_2 > best_organism[1]:
                        # ID, EF, Nodes, Penalty applied
                        best_organism = max_organism
                        changed_best_score = True

                # END FOR j

            # END FOR i

        # Mean fitness in the population
        mean_fitness = np.mean(a_fitness)
        # Standard deviation of fitness in the population
        standard_dev_fitness = np.std(a_fitness)
        # Inequality of fitness in the population (measured with the Gini coefficient)
        gini_fitness = gini_RSV(a_fitness)
        # Mean number of nodes per organism in the population
        mean_nodes = np.mean(a_nodes)

        # Show IDs of final array
        # print("-"*10)
        _m, _s = divmod((time.time() - initial), 60)
        _h, _m = divmod(_m, 60)
        s_time = "{}h:{}m:{:.2f}s".format(int(_h), int(_m), _s)
        print_ln(
            (
                "Iter: {} AF:{:.2f} SDF:{:.2f} GF:{:.2f} AN:{:.2f}"
                + " - MO: {} MF: {:.2f} MN: {} MP: {:.2f}"
                + " -  BO: {} BF: {:.2f} BN: {} BP: {:.2f} Time: {}"
            ).format(
                iterations,  # "Iter"
                mean_fitness,  # "AF"
                standard_dev_fitness,  # "SDF"
                gini_fitness,  # "GF"
                mean_nodes,  # "AN"
                max_organism[0]._id,  # "MO"
                max_organism[1],  # "MF" (fitness)
                max_organism[2],  # "MN" (nodes)
                max_organism[3],  # "MP" (penalty)
                best_organism[0]._id,  # "BO"
                best_organism[1],  # "BF" (fitness)
                best_organism[2],  # "BN" (nodes)
                best_organism[3],  # "BP" (penalty)
                s_time,  # Time
            ),
            RESULT_BASE_PATH_DIR + OUTPUT_FILENAME,
        )

        # Print against a random positive sequence
        random.shuffle(positive_dataset)
        print(max_organism[0].print_result(positive_dataset[0]))

        # Export organism if new best organism
        if changed_best_score:
            filename = "{}_{}".format(
                time.strftime(timeformat), best_organism[0]._id
            )
            export_organism(
                best_organism[0], positive_dataset, filename, organism_factory
            )
        # Periodic organism export
        if iterations % PERIODIC_EXPORT == 0:
            filename = "{}_{}".format(
                time.strftime(timeformat), max_organism[0]._id
            )
            export_organism(
                max_organism[0], positive_dataset, filename, organism_factory
            )

        # print("-"*10)
        iterations += 1
        # END WHILE

    # TODO: Maybe a good idea to export the full population after all
    # organism_factory.export_organisms(organism_population,
    #         RESULT_BASE_PATH_DIR+"final_population.json")


def is_finished(
        method: str, iterations: int, max_score: float, last_max_score: float
) -> bool:
    """Checks if main while loop is finished
    methods: 'Iterations', 'minScore', 'Threshold'

    Args:
        method: Name of the finishing method
        max_score: max score recorded on the current iteration
        last_max_score: max score recorded on the laset iteration
        iterations: Number of the current iteration

    Returns:
        True if program should finnish.
        False otherwise
    """

    if method.lower() == "iterations":
        return iterations >= MIN_ITERATIONS

    if method.lower() == "fitness":
        return max_score >= MIN_FITNESS

    if method.lower() == "threshold":
        return abs(last_max_score - max_score) <= THRESHOLD

    return True


def export_organism(
        organism, dataset: list, filename: str, factory: OrganismFactory
) -> None:
    """Exports a single organism in json format, visual format and its
    recognizers binding

    Args:
        organism (OrganismObject): organism to export
        dataset: Sequences to check the organism binding
        filename: Previous info to export filenames. Common in all filenames
        factory: Used to export in json format
    """

    organism_file = "{}{}_organism.txt".format(RESULT_BASE_PATH_DIR, filename)
    organism_file_json = "{}{}_organism.json".format(
        RESULT_BASE_PATH_DIR, filename
    )
    results_file = "{}{}_results.txt".format(RESULT_BASE_PATH_DIR, filename)

    organism.export(organism_file)
    organism.export_results(dataset, results_file)
    factory.export_organisms([organism], organism_file_json)


def combine_organisms(
        organism1, organism2, organism_factory: OrganismFactory
) -> dict:
    """Gets 2 organisms, and returns 2 children with format
    (child, similarity to parent 1, similarity to parent 2)

    Args:
        organism1 (OrganismObject): First organism for the crossover
        organism2 (OrganismObject): Second organism for the crossover
        organism_factory: factory used to set the children ids

    Retruns:
        A dictionary with 2 children:
        "child1": child derived from organism1
        "child2": child derived from organism1
    """
    # Save the number of nodes from the parents
    n_nodes_org_1 = organism1.count_nodes()
    n_nodes_org_2 = organism2.count_nodes()

    # Create the 2 childs and assign new IDs
    child1 = copy.deepcopy(organism1)
    child2 = copy.deepcopy(organism2)

    # Assign IDs to organisms and increase factory counter
    child1.set_id(organism_factory.get_id())
    child2.set_id(organism_factory.get_id())

    # Combine parents with probability p
    if random.random() < RECOMBINATION_PROBABILITY:

        # Select random nodes to swap from each child
        random_node_1 = random.randint(0, n_nodes_org_1 - 1)
        random_node_2 = random.randint(0, n_nodes_org_2 - 1)
        node1 = child1.get_node(random_node_1)
        node2 = child2.get_node(random_node_2)

        # Save the number of nodes taken from each  child
        n_nodes_from_org_1 = node1.count_nodes()
        n_nodes_from_org_2 = node2.count_nodes()

        # Get parents nodes of swapping nodes before swap
        parent_node_1 = child1.get_parent(node1._id)
        parent_node_2 = child2.get_parent(node2._id)

        # Swap nodes
        # Set nodes in oposite children
        # based on recipient parent node determine if incoming node goes to
        # left descendent or right descendent
        # if recipient node is root, subsitute with incoming
        if parent_node_1["is_root_node"]:
            # Its the root node of child 1
            child1.set_root_node(node2)
        else:
            if parent_node_1["is_left_side"]:
                # Child on left side
                parent_node_1["self"].set_node1(node2)
            else:
                # Child on right side
                parent_node_1["self"].set_node2(node2)

        if parent_node_2["is_root_node"]:
            # Its the root node of child 2
            child2.set_root_node(node1)
        else:
            if parent_node_2["is_left_side"]:
                # Child on left side
                parent_node_2["self"].set_node1(node1)
            else:
                # Child on right side
                parent_node_2["self"].set_node2(node1)

        n_nodes_child_1 = child1.count_nodes()
        n_nodes_child_2 = child2.count_nodes()

        # Reset children node IDs across the organism
        child1.reset_ids()
        child2.reset_ids()

        # dictionary with an organism and similarities to each parent
        # similatiries are computed as the number of nodes shared  between
        # each parent and child
        child_1_similarities = {
            "sim_org_1": (n_nodes_child_1 - n_nodes_from_org_2)
                         / n_nodes_child_1,
            "sim_org_2": n_nodes_from_org_2 / n_nodes_child_1,
            "child": child1,
        }

        child_2_similarities = {
            "sim_org_1": n_nodes_from_org_1 / n_nodes_child_2,
            "sim_org_2": (n_nodes_child_2 - n_nodes_from_org_1)
                         / n_nodes_child_2,
            "child": child2,
        }
    else:

        # If children are not recombined, return the same organisms and their
        # similarities
        child_1_similarities = {
            "sim_org_1": 1,  # Equal to organism 1
            "sim_org_2": 0,
            "child": child1,
        }

        child_2_similarities = {
            "sim_org_1": 0,
            "sim_org_2": 1,  # Equal to organism2
            "child": child2,
        }

    return {"child1": child_1_similarities, "child2": child_2_similarities}


def set_up():
    """Reads configuration file and sets up all program variables

    """

    # specify as global variable so it can be accesed in local
    # contexts outside setUp

    global END_WHILE_METHOD
    global POPULATION_LENGTH
    global DATASET_BASE_PATH_DIR
    global RESULT_BASE_PATH_DIR
    global POSITIVE_FILENAME
    global NEGATIVE_FILENAME
    global RESULT_PATH_PATH_DIR
    global MAX_SEQUENCES_TO_FIT_POS
    global MAX_SEQUENCES_TO_FIT_NEG
    global FITNESS_FUNCTION
    global GENOME_LENGTH
    global INEQUALITY_PENALTY_METHOD
    global INEQUALITY_PENALTY_PARAM
    global MIN_ITERATIONS
    global MIN_FITNESS
    global THRESHOLD
    global POPULATION_ORIGIN
    global POPULATION_FILL_TYPE
    global INPUT_FILENAME
    global OUTPUT_FILENAME
    global RECOMBINATION_PROBABILITY
    global PERIODIC_EXPORT
    global MAX_NODES
    global MIN_NODES

    # Config data
    global configOrganism
    global configOrganismFactory
    global configConnector
    global configPssm

    config = read_json_file(JSON_CONFIG_FILENAME)
    # Store config variables for main function
    POPULATION_LENGTH = config["main"]["POPULATION_LENGTH"]
    DATASET_BASE_PATH_DIR = config["main"]["DATASET_BASE_PATH_DIR"]
    RESULT_BASE_PATH_DIR = (
        config["main"]["RESULT_BASE_PATH_DIR"]
        + time.strftime("%Y%m%d%H%M%S")
        + "/"
    )
    POSITIVE_FILENAME = config["main"]["POSITIVE_FILENAME"]
    NEGATIVE_FILENAME = config["main"]["NEGATIVE_FILENAME"]
    MAX_SEQUENCES_TO_FIT_POS = config["main"]["MAX_SEQUENCES_TO_FIT_POS"]
    MAX_SEQUENCES_TO_FIT_NEG = config["main"]["MAX_SEQUENCES_TO_FIT_NEG"]
    FITNESS_FUNCTION = config["main"]["FITNESS_FUNCTION"]
    GENOME_LENGTH = config["main"]["GENOME_LENGTH"]
    INEQUALITY_PENALTY_METHOD = config["main"]["INEQUALITY_PENALTY_METHOD"]
    INEQUALITY_PENALTY_PARAM = config["main"]["INEQUALITY_PENALTY_PARAM"]
    MIN_ITERATIONS = config["main"]["MIN_ITERATIONS"]
    MIN_FITNESS = config["main"]["MIN_FITNESS"]
    THRESHOLD = config["main"]["THRESHOLD"]
    END_WHILE_METHOD = config["main"]["END_WHILE_METHOD"]
    POPULATION_ORIGIN = config["main"]["POPULATION_ORIGIN"]
    POPULATION_FILL_TYPE = config["main"]["POPULATION_FILL_TYPE"]
    INPUT_FILENAME = config["main"]["INPUT_FILENAME"]
    OUTPUT_FILENAME = config["main"]["OUTPUT_FILENAME"]
    RECOMBINATION_PROBABILITY = config["main"]["RECOMBINATION_PROBABILITY"]
    PERIODIC_EXPORT = config["main"]["PERIODIC_EXPORT"]
    MAX_NODES = config["organism"]["MAX_NODES"]
    MIN_NODES = config["organism"]["MIN_NODES"]

    # Create directory where the output and results will be stored
    os.mkdir(RESULT_BASE_PATH_DIR)

    # Store Config into variables to use later
    configOrganism = config["organism"]
    configOrganismFactory = config["organismFactory"]
    configConnector = config["connector"]
    configPssm = config["pssm"]

    # Throw config on a file
    parameters_path = RESULT_BASE_PATH_DIR + "parameters.txt"
    print_ln("-" * 50, parameters_path)
    print_ln(" " * 20 + "PARAMETERS", parameters_path)
    print_ln("-" * 50, parameters_path)

    print_config_json(config["main"], "Main Config", parameters_path)
    print_config_json(configOrganism, "Organism Config", parameters_path)
    print_config_json(
        configOrganismFactory, "Organism Factory Config", parameters_path
    )
    print_config_json(configConnector, "Connector Config", parameters_path)
    print_config_json(configPssm, "PSSM Config", parameters_path)

    print_ln("-" * 50, parameters_path)


def read_fasta_file(filename: str) -> list:
    """Reads a fasta file and returns an array of DNA sequences (strings)

    TODO: probably it can be useful to create our own Sequence object that
    creates the string and stores some properties from fasta format. Also
    we can adapt the current program to use Biopythons's Seq object.

    Args:
        filename: Name of the file that contains FASTA format sequences to read

    Returns:
        The set of sequences in string format

    """
    dataset = []

    fasta_sequences = SeqIO.parse(open(filename), "fasta")

    for fasta in fasta_sequences:
        dataset.append(str(fasta.seq).lower())

    return dataset


def read_json_file(filename: str) -> dict:
    """Reads a JSON file and returns a dictionary with the content

    Args:
        filename: Name of the json file to read

    Returns:
        Dictionary with the json file info

    """

    with open(filename) as json_content:

        return json.load(json_content)


def print_config_json(config: dict, name: str, path: str) -> None:
    """Print the config file on std out and send it to a file.
    It is useful so we can know which was the configuration on every run

    Args:
        config: Configuration file to print
        name: Title for the configuration file
        path: File to export the configuration info
    """
    print_ln("{}:".format(name), path)

    for key in config.keys():
        print_ln("{}: {}".format(key, config[key]), path)
    print_ln("\n", path)


def print_ln(string: str, name_file: str) -> None:
    """Shows the string on stdout and write it to a file
    (like the python's logging modules does)

    Args:
        string: Information to print on stdout and file
        name_file: path to the file to export the string
    """

    print(string)

    # Here we are sure file exists
    _f = open(name_file, "a+")
    _f.write(string + "\n")
    _f.close()


def gini_RSV(values_for_each_class):
    '''
    Gini coefficient, modified in order to be alble to deal with negative
    values as in "Inequality measures and the issue of negative incomes"
    (Raffinetti, Siletti, Vernizzi)

    Parameters
    ----------
    values_for_each_class : array-like object
        Values associated to each class.
        They don't need to be already sorted and/or normalized.
        They can also be negative.

    Returns
    -------
    giniRSV : float
        Ranges from 0 (perfect equality) to 1 (maximal inequality).

    '''
    
    N = len(values_for_each_class)
    
    numerator = 0
    for i in values_for_each_class:
        for j in values_for_each_class:
            numerator += abs(i - j)
    
    pos = 0  # sum over the positive values
    neg = 0  # sum over the negative values (in absolute value)
    for x in values_for_each_class:
        if x >= 0:
            pos += x
        else:
            neg += -x
    
    mu_RSV = (N - 1) * (pos + neg) / N**2  # modified mu parameter
    
    if mu_RSV == 0:
        # Manage two special cases (avoiding 0-division error):
        #   - when a single value is the input
        #   - when all the values in the input are 0
        # In both cases mu_RSV will be 0
        # No inequality is measurable, and therefore 0 is returned
        return 0
    denominator = 2 * N**2 * mu_RSV
    giniRSV = numerator / denominator
    
    return giniRSV


# Entry point to app execution
# It calculates the time, but could include other app stats

if __name__ == "__main__":

    INITIAL = time.time()
    # Reads configuration file and sets up all program variables
    set_up()

    # Profiling
    PROFILER = cProfile.Profile()
    PROFILER.enable()

    # Main function
    main()

    # Profiler output
    PROFILER.disable()
    STRING = io.StringIO()
    SORTBY = "cumulative"
    PS = pstats.Stats(PROFILER, stream=STRING).sort_stats(SORTBY)
    PS.print_stats()
    # print(STRING.getvalue())

    # Print final execution time and read parameters
    _M, _S = divmod((time.time() - INITIAL), 60)
    _H, _M = divmod(_M, 60)
    print_ln(
        "Time: {}h:{}m:{:.2f}s".format(int(_H), int(_M), _S),
        RESULT_BASE_PATH_DIR + "parameters.txt",
    )
