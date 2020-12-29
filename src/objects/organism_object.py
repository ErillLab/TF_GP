"""
Oragnism object
It allocates the full data structure
"""

import random
import numpy as np


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


class OrganismObject:
    """Organism object
    """

    def __init__(self, _id: int, conf: dict, max_pssm_length: int) -> None:
        """Organism constructor

        Args:
            _id: Organism identifier
            conf: Configuration from JSON file
            max_pssm_length: Maximum column size of the pssm recognizer
        """
        self._id = _id
        
        self.recognizers = []  # < - - - - - - - - - - NEW ATTRIBUTE
        self.connectors = []  # < - - - - - - - - - - NEW ATTRIBUTE
        
        #self.num_nodes = 0
        #self.num_recognizers = 0

        self.cumulative_fit_method = conf["CUMULATIVE_FIT_METHOD"]
        
        self.energy_threshold_method = conf["ENERGY_THRESHOLD_METHOD"]
        self.energy_threshold_value = conf["ENERGY_THRESHOLD_PARAM"]
        
        self.mutate_probability_substitute_pssm = conf[
            "MUTATE_PROBABILITY_SUBSTITUTE_PSSM"
        ]
        
        self.indel_method = conf[
            "INDEL_METHOD"
        ]
        self.mutate_probability_delete_recognizer = conf[
            "MUTATE_PROBABILITY_DELETE_RECOGNIZER"
        ]
        self.mutate_probability_insert_recognizer = conf[
            "MUTATE_PROBABILITY_INSERT_RECOGNIZER"
        ]
        self.mutate_probability_node_mutation = conf[
            "MUTATE_PROBABILITY_NODE_MUTATION"
        ]
        self.min_nodes = conf["MIN_NODES"]
        self.max_nodes = conf["MAX_NODES"]
        self.is_tracked = False

        self.is_automatic_placement_options = conf[
            "AUTOMATIC_PLACEMENT_OPTIONS"
        ]
        self.max_pssm_length = max_pssm_length


    # Setters an getters

    def set_is_tracked(self, new_tracked: bool):
        """Setter is_tracked

        Args:
            new_tracked: True if the organism should be tracked.
                         False otherwise.
        """
        self.is_tracked = new_tracked

    def mutate(self, org_factory) -> None:
        """Mutates an organism based on JSON configured probabilities

        Args:
            org_factory (OrganismFactory): Factory of organisms and node
                                           components
        """

        
        # Delete a recognizer (and one parent connector)
        if random.random() < self.mutate_probability_delete_recognizer:
            
            if self.indel_method == "blind":
    
                n_recognizers = self.count_recognizers()
                # Choose randomly the recognizer to be deleted
                recognizer_idx = random.randint(0, n_recognizers - 1)
                if recognizer_idx == 0:
                    # If the recognizer to delete is the first, the connector to
                    # delete is the one to the right
                    # ( same index: connector_idx = recognizer_idx = 0 )
                    connector_idx = recognizer_idx
                elif recognizer_idx == n_recognizers - 1:
                    # If the recognizer to delete is the last, the connector to
                    # delete is the one to the left
                    # ( different index: connector_idx = recognizer_idx - 1 )
                    connector_idx = recognizer_idx - 1
                else:
                    # if the recognizer to delete is in not a terminal recognizer
                    # of the chain, the parent connector to be deleted (left/riht)
                    # is chosen randomly
                    if random.random() < 0.5:
                        connector_idx = recognizer_idx
                    else:
                        connector_idx = recognizer_idx - 1
                
                new_recognizers = (self.recognizers[:recognizer_idx] +
                                   self.recognizers[recognizer_idx + 1:])
                new_connectors = (self.connectors[:connector_idx] +
                                   self.connectors[connector_idx + 1:])
                
                self.recognizers = new_recognizers
                self.connectors = new_connectors
            
            
            if self.indel_method == "intelligent":
                
                n_recognizers = self.count_recognizers()
                # Choose randomly the recognizer to be deleted
                recognizer_idx = random.randint(0, n_recognizers - 1)
                if recognizer_idx == 0:
                    # If the recognizer to delete is the first, the connector to
                    # delete is the one to the right
                    # ( same index: connector_idx = recognizer_idx = 0 )
                    connector_idx = recognizer_idx
                elif recognizer_idx == n_recognizers - 1:
                    # If the recognizer to delete is the last, the connector to
                    # delete is the one to the left
                    # ( different index: connector_idx = recognizer_idx - 1 )
                    connector_idx = recognizer_idx - 1
                else:
                    # if the recognizer to delete is in not a terminal recognizer
                    # of the chain, the parent connector to be deleted (left/riht)
                    # is chosen randomly
                    if random.random() < 0.5:
                        # Index of the parent connector that will be deleted
                        connector_idx = recognizer_idx
                        # Index of the parent connector that will be adjusted
                        connector_to_stretch = recognizer_idx - 1
                    else:
                        # Index of the parent connector that will be deleted
                        connector_idx = recognizer_idx - 1
                        # Index of the parent connector that will be adjusted
                        connector_to_stretch = recognizer_idx
                    
                    # Adjust parameters of the neighbour connector
                    ''' The parent connector that is not deleted is modified,
                    so that it can span the gap left by the deletion, witout
                    heavily affecting the placement of the nodes to the side of
                    the deletion point.
                    '''
                    # The average distance (mu) of the deleted connector is
                    # added to its mu
                    adj_mu = (self.connectors[connector_to_stretch]._mu +
                              self.connectors[connector_idx]._mu)
                    # Its standard deviation (sigma) becomes the square root of
                    # the sum of the squares of the sigmas (under assumption of
                    # independence between the two random variables)
                    adj_sigma = (self.connectors[connector_to_stretch]._sigma ** 2 +
                                self.connectors[connector_idx]._sigma ** 2)**(1/2)
                    # set new mu and new sigma
                    self.connectors[connector_to_stretch].set_mu(adj_mu)
                    self.connectors[connector_to_stretch].set_sigma(adj_sigma)
                
                new_recognizers = (self.recognizers[:recognizer_idx] +
                                   self.recognizers[recognizer_idx + 1:])
                new_connectors = (self.connectors[:connector_idx] +
                                   self.connectors[connector_idx + 1:])
                
                self.recognizers = new_recognizers
                self.connectors = new_connectors
        
        
        
        # Insert a recognizer (and one parent connector)
        if random.random() < self.mutate_probability_insert_recognizer:
            
            new_connector = org_factory.create_connector()
            new_recognizer = org_factory.create_pssm()
            
            if self.indel_method == "blind":
                
                n_recognizers = self.count_recognizers()
                # Choose randomly the recognizer next to which the insertion
                # is going to occur
                recognizer_idx = random.randint(0, n_recognizers - 1)
                # Choose randomly whether the insertion is going to be to the
                # left or to the right of the considered recognizer
                
                if random.random() < 0.5:  # Insertion occurs to the left
                    # First connector after insertion point and first
                    # recognizer after insertion point
                    connector_idx = recognizer_idx
                    
                else:  # Insertion occurs to the right
                    # First connector after insertion point
                    connector_idx = recognizer_idx
                    # First recognizer after insertion point
                    recognizer_idx += 1
                
                new_recognizers = (self.recognizers[:recognizer_idx] +
                                   [new_recognizer] +
                                   self.recognizers[recognizer_idx:])
                new_connectors = (self.connectors[:connector_idx] +
                                   [new_connector] +
                                   self.connectors[connector_idx:])
                
                self.recognizers = new_recognizers
                self.connectors = new_connectors
            
            
            if self.indel_method == "intelligent":
                
                n_recognizers = self.count_recognizers()
                # Choose randomly the recognizer next to which the insertion
                # is going to occur
                recognizer_idx = random.randint(0, n_recognizers - 1)
                # Choose randomly whether the insertion is going to be to the
                # left or to the right of the considered recognizer
                connector_to_compress = None
                if random.random() < 0.5:  # Insertion occurs to the left
                    # First connector after insertion point and first
                    # recognizer after insertion point
                    connector_idx = recognizer_idx
                    
                    if recognizer_idx != 0:
                        connector_to_compress = recognizer_idx - 1
                        # (No compression is required if insertion occurs to
                        # the left of the first recognizer of the chain)
                        
                else:  # Insertion occurs to the right
                    # First connector after insertion point
                    connector_idx = recognizer_idx
                    
                    if recognizer_idx != n_recognizers - 1:
                        connector_to_compress = recognizer_idx
                        # (No compression is required if insertion occurs to
                        # the right of the last recognizer of the chain)
                    
                    # First recognizer after insertion point
                    recognizer_idx += 1
                
                if connector_to_compress != None:
                    
                    # Adjust MUs
                    mu_old_connector = self.connectors[connector_to_compress]._mu
                    mu_new_connector = new_connector._mu
                    current_sum_mus = mu_old_connector + mu_new_connector
                    expected_sum_mus = mu_old_connector
                    
                    mu_scaling_factor =  expected_sum_mus / current_sum_mus
                    # Compress the neighbour connector
                    self.connectors[connector_to_compress].set_mu(
                        mu_old_connector * mu_scaling_factor
                    )
                    # Compress the new inserted connector
                    new_connector.set_mu(
                        mu_new_connector * mu_scaling_factor
                    )
                    
                    # Adjust SIGMAs
                    var_old_connector = self.connectors[connector_to_compress]._sigma ** 2
                    var_new_connector = new_connector._sigma**2
                    current_sum_variances = var_old_connector + var_new_connector
                    expected_sum_variances = var_old_connector
                    
                    var_scaling_factor =  expected_sum_variances / current_sum_variances
                    # Compress the neighbour connector
                    self.connectors[connector_to_compress].set_sigma(
                        np.sqrt(var_old_connector * var_scaling_factor)
                    )
                    # Compress the new inserted connector
                    new_connector.set_sigma(
                        np.sqrt(var_new_connector * var_scaling_factor)
                    )
                
                
                new_recognizers = (self.recognizers[:recognizer_idx] +
                                   [new_recognizer] +
                                   self.recognizers[recognizer_idx:])
                new_connectors = (self.connectors[:connector_idx] +
                                   [new_connector] +
                                   self.connectors[connector_idx:])
                
                self.recognizers = new_recognizers
                self.connectors = new_connectors
        
        '''
        # Mutate a random node
        if random.random() < self.mutate_probability_node_mutation:

            n_nodes = self.count_nodes()
            random_node_idx = random.randint(0, n_nodes - 1)
            if random_node_idx < self.count_recognizers():
                # mutate a recognizer
                self.recognizers[random_node_idx].mutate(org_factory)
            else:
                # mutate a connector
                connector_idx = random_node_idx - self.count_recognizers()
                self.connectors[connector_idx].mutate(org_factory)
        '''


    def get_seq_fitness(self, s_dna: str) -> dict:
        """Return the fitness of the organism for a given DNA sequence

        Args:
            s_dna: DNA sequence to analize

        Returns:
           score, blocked and blockers

        Description:
            This function implements the placement behavior for organism.
            The placement problem is defined as who to best position an
            organism on a sequence (i.e. how to maximize its fitness given
            the sequence).
            The implementation in this function follows the recursive 
            formulation of the organism, calling on the organisms root node
            "get_placement" function in order obtain the list of best possible
            placements for the organism. The best placement is chosen.
            
            If the root organism returns an empty list, because the placement
            procedure has been unable to identify a non-conflicting placement,
            the function sets the energy to a large negative value.
            
            The overall placement strategy, implemented via the recursive 
            calling to the get_placement function of connectors and PSSMs is 
            as follows:
                
                - PSSMs evaluate the sequence and return a list of candidate
                  placements, with a list of blocked positions for those 
                  placemenets
                - PSSMs connectors receive list of candidate positions for 
                  PSSMs, determine all viable (non-conflicting) combinations,
                  rank them taking into account their energy, and return them
                - General connectors receive list of candidate positions for 
                  lower connectors, determine viable combinations, rank them 
                  and return them
            
            Notes on the placement implementation:
                - This is a greedy placement strategy, in the sense that a 
                  connector only "sees" its subtree, so sister connectors
                  may independently choose placement options with conflicting
                  PSSM placements. The placement option number for connectors
                  hence plays a crucial role in determining whether the 
                  organism will be capable of identifying a correct placement.
                  This problem will be more severe as the depth of the organism
                  grows.
        """

        # Compute the number of nodes for automatic placement
        num_pssm = (self.count_nodes() - 1) / 2 + 1

        automatic_placement_options = int((self.max_pssm_length + 1) * num_pssm)
        
        # Set energy threshold method and value
        E_threshold_method = self.energy_threshold_method
        E_threshold_value = self.energy_threshold_value

        # call recursively to get the list of best placement options for 
        # the organism (root node)
        node_root = self.root_node.get_placement(
            s_dna,
            len(s_dna),
            automatic_placement_options,
            self.is_automatic_placement_options
                )

        # handle the case in which no viable placement options have been
        # identified
        if len(node_root) < 1:
            print("Too few placement options")
            return {
                "energy": -1000,
                "position": 0,
                "lock_vector": [],
                "recognizers_scores": []
                }

        # Apply lower bound to energy if required
        if E_threshold_method == "organism":
            if node_root[0]["energy"] < E_threshold_value:
                node_root[0]["energy"] = E_threshold_value
        
        # return score, blocks and blokcers and PSSMs scores in that sequence
        return node_root[0]
    
    def get_binding_energies(self, a_dna: list) -> list:
        """Return a list containing the binding energies on the sequences
        provided as input data

        Args:
            a_dna: list of dna sequences

        Returns:
            list of all the binding energies
        """
        scores = []
        ginis = []
        for s_dna in a_dna:
            sfit = self.get_seq_fitness(s_dna)
            energy = sfit["energy"]
            pssm_scores = sfit["recognizers_scores"]
            
            if len(pssm_scores) > 0:
                gini = gini_RSV(pssm_scores)  # Gini coefficient
                ginis.append(gini)
            
            scores.append(energy)

        # Compute the average Gini coefficient as the geometric mean
        if len(ginis) == 0:  # Case where no placement was available for any sequence
            avg_gini = 1  # maximum penalty is arbitrarily assigned
        else:
            avg_gini = np.prod(ginis) ** (1/len(ginis))  # geometric mean
        
        return {"energy_scores": scores, "avg_gini": avg_gini}
    
    def get_discriminative_fitness(self, a_dna: list) -> float:
        """Return the total Fitness for an array of DNA sequences and the
        fitness method

        Args:
            a_dna: list of dna sequences

        Returns:
            score assigned to the organism
        """

        scores = []
        ginis = []
        for s_dna in a_dna:
            sfit = self.get_seq_fitness(s_dna)
            energy = sfit["energy"]  # energy
            pssm_scores = sfit["recognizers_scores"]  # PSSMs scores
            
            if len(pssm_scores) > 0:
                gini = gini_RSV(pssm_scores)  # Gini coefficient
                ginis.append(gini)

            scores.append(energy)
        
        if self.cumulative_fit_method == "sum":
            # Compute fitness score as sum over the positive scores
            score = np.sum(scores)
        
        elif self.cumulative_fit_method == "mean":
            # Compute fitness score as average positive score
            score = np.mean(scores)
        
        # Compute the average Gini coefficient as the geometric mean
        if len(ginis) == 0:  # Case where no placement was available for any sequence
            avg_gini = 1  # maximum penalty is arbitrarily assigned
        else:
            avg_gini = np.prod(ginis) ** (1/len(ginis))  # geometric mean
        
        return {"score": score, "avg_gini": avg_gini}

    def get_boltz_fitness(self, pos_dataset: list, neg_dataset: list,
                          genome_length: int) -> float:
        """Returns the organism's fitness, defined as the probability that the regulator binds a
        positive sequence. All the binding energies are turned into probabilities according to a
        Boltzmannian distribution. The probability of binding a particular sequence, given the binding
        energy on that sequence, is p = e**binding_energy / Z
        where Z is the partition function.
        A high number of negative sequences is assumed to be present (emulating the environment of a
        regulator that needs to find its targets on an entire genome).
        A coefficient called neg_factor is computed, so that the value of Z can be as high as if there
        were as	many negative sequences as required to cover the entire genome.

        Args:
            pos_dataset: list of dna sequences in the positive dataset
            neg_dataset: list of dna sequences in the negative dataset
	    genome_length: integer representing the length of the genome

        Returns:
            fitness assigned to the organism
        """
        
        # Values onthe positive set
        pos_values = []
        ginis = []
        for s_dna in pos_dataset:
            sfit = self.get_seq_fitness(s_dna)
            boltz_exp = np.e**sfit["energy"]  # exp(energy)
            pssm_scores = sfit["recognizers_scores"]  # PSSMs scores
            if len(pssm_scores) > 0:
                gini = gini_RSV(pssm_scores)  # Gini coefficient
                ginis.append(gini)

            pos_values.append(boltz_exp)
        
        # Compute the average Gini coefficient as the geometric mean
        if len(ginis) == 0:  # Case where no placement was available for any sequence
            avg_gini = 1  # maximum penalty is arbitrarily assigned
        else:
            avg_gini = np.prod(ginis) ** (1/len(ginis))  # geometric mean
        
        # Values onthe negative set
        neg_values = []
        neg_lengths = []
        for s_dna in neg_dataset:
            sfit = self.get_seq_fitness(s_dna)
            boltz_exp = np.e**sfit["energy"]  # exp(energy)
            neg_values.append(boltz_exp)
            neg_lengths.append(len(s_dna))
        
        # Scaling factor, used to over-represent the negative scores, so that
        # it simulates a genome of specified length
        neg_factor = genome_length//sum(neg_lengths)
        
        # Partition function
        Z = sum(pos_values) + neg_factor * sum(neg_values)
        
        # Compute fitness score as a Boltzmannian probability
        boltz_fitness = sum(pos_values) / Z
        
        return {"score": boltz_fitness, "avg_gini": avg_gini}

    def count_nodes(self) -> int:
        """Returns the number of nodes of the organism

        Returns:
            Number of nodes of the organism
        """

        return 2 * len(self.recognizers) - 1
    
    def count_connectors(self) -> int:
        """Returns the number of connectors of the organism

        Returns:
            Number of connectors.
        """
        
        return len(self.connectors)

    def count_recognizers(self) -> int:
        """Returns the number of recognizers of the organism

        Returns:
            Number of recognizers.
        """
        
        return len(self.recognizers)
    
    def get_random_connector(self) -> int:
        """Returns the index of a random connector of the organism

        Returns:
            Integer between 0 and N-1 (both included), where N is the number of
            connectors the organism has.
        """
        
        num_connectors =  self.count_connectors()
        
        return random.randint(0, num_connectors - 1)

    def get_random_recognizer(self) -> int:
        """Returns the index of a random recognizer of the organism

        Returns:
            Integer between 0 and N-1 (both included), where N is the number of
            recognizers the organism has.
        """
        
        num_recognizers =  self.count_recognizers()
        
        return random.randint(0, num_recognizers - 1)
    
    def break_chain(self, connector_to_break, bond_to_keep):
        """Brakes an organism at the specified link and returns the resulting
        pair of chunks into a list. Each chunk is a dictionary with two keys:
        "recognizers" and "connectors".

        Parameters
        ----------
        connector_to_break : int
            Index of the connector where the chain will be broken.
        bond_to_keep : str
            If "left" the connector where the split occurs will stay linked to
            the left chunk, while its right bond will be broken (the opposite
            happens if its value is "right".

        Returns
        -------
        list
            A list with two elements, which are the two chunks of the splitted
            chain, both represented as a dictionary with two keys:
            "recognizers" and "connectors" (which point to lists of recognizers
            or connectors, respectively).

        """
        
        # Recognizers of left and right chunks
        L_recognizers = self.recognizers[:connector_to_break + 1]
        R_recognizers = self.recognizers[connector_to_break + 1:]
        
        # Connectors of left and right chunks
        if bond_to_keep=="left":
            L_connectors = self.connectors[:connector_to_break + 1]
            R_connectors = self.connectors[connector_to_break + 1:]
        elif bond_to_keep=="right":
            L_connectors = self.connectors[:connector_to_break]
            R_connectors = self.connectors[connector_to_break:]
        else:
            raise Exception('bond_to_keep needs to be "left" or "right".')
        
        L_chunk = {"recognizers": L_recognizers, "connectors": L_connectors}
        R_chunk = {"recognizers": R_recognizers, "connectors": R_connectors}
        
        return [L_chunk, R_chunk]
    
    def set_connectors(self, connectors_list):
        """Set the connectors of the organism to be those provided in the input
        list.
        """
        
        self.connectors = connectors_list
        
    def set_recognizers(self, recognizers_list):
        """Set the recognizers of the organism to be those provided in the
        input list.
        """
        
        self.recognizers = recognizers_list

    def print(self) -> None:
        """Prints the whole tree data structure
        """
        
        print("***** Organism {} *****".format(self._id))
        for i in range(len(self.recognizers) - 1):
            self.recognizers[i].print()
            self.connectors[i].print()
        self.recognizers[-1].print()

    def export(self, filename: str) -> None:
        """Exports the whole tree data structure

        Args:
            filename: Name of the file to export the organism
        """
        organism_file = open(filename, "w+")
        organism_file.write("***** Organism {} *****".format(self._id))
        
        for i in range(len(self.recognizers) - 1):
            self.recognizers[i].export(organism_file)
            self.connectors[i].export(organism_file)
        self.recognizers[-1].export(organism_file)

        organism_file.write("\n")
        organism_file.close()

    def export_results(self, a_dna: list, filename: str) -> None:
        """Exports all DNA sequences organism binding to a file

        Args:
            filename: Name of the file to export sequences
            a_dna: list fo sequences to export

        """

        # Sort the array, so its always shown in the same order
        # sorting is done by sequence, so first sequences start with "AAA.."
        a_dna.sort()

        results_file = open(filename, "w+")

        # for every DNA sequence
        for s_dna in a_dna:

            # call fitness evaluation for sequence
            sfit = self.get_seq_fitness(s_dna.lower())

            # write out the sequence
            results_file.write("\n{}\n".format(s_dna))

            # create an empy positions map
            map_positions = "-" * len(s_dna)

            # positions for PSSMs are in blocks and blocked lists, returned by
            # getSeqFitness. we zip them and then iterate over the zip to
            # print the PSSMs in their locations respective to the sequence
            positions = sfit["lock_vector"]
            for pssm in positions:
                # print _id, capped to the length of PSSM
                _p = round(pssm["position"])
                pssm_str = (str(pssm["id"]) * pssm["length"])[:pssm["length"]]

                # fill up map at correct positions
                map_positions = (
                    map_positions[:_p] + pssm_str + map_positions[_p + pssm["length"]:]
                )

            # write map to file for this sequence
            results_file.write(map_positions + "\n")
            # resultsFile.write(str(stuff) + "\n")

        results_file.close()

    def print_result(self, s_dna: str) -> str:
        """Prints the results of s_dna binding sites to stdout

        Args:
            s_dna: DNA sequence to export

        Returns:
            DNA sequence and binding sites of the organisms recognizer
        """

        s_dna = s_dna.lower()

        # call fitness evaluation for sequence
        sfit = self.get_seq_fitness(s_dna.lower())

        # create an empy positions map
        map_positions = "-" * len(s_dna)
        
        # positions for PSSMs are in blocked and blocked lists, returned by
        # getSeqFitness. we zip them and then iterate over the zip to
        # print the PSSMs in their locations respective to the sequence
        positions = sfit["lock_vector"]
        for pssm in positions:
            # print _id, capped to the length of PSSM
            _p = round(pssm["position"])
            pssm_str = (str(pssm["id"]) * pssm["length"])[:pssm["length"]]

            # fill up map at correct positions
            map_positions = (
                map_positions[:_p] + pssm_str + map_positions[_p + pssm["length"]:]
            )
            # handle two-digit positions, by alterning between digits

        # return map for this sequence
        return "{}\n{}".format(s_dna, map_positions)
