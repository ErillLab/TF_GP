"""Organism Factory creates organisms, connectors and pssms
"""

import random
import json
import numpy as np
from .organism_object import OrganismObject
from .connector_object import ConnectorObject
from .pssm_object import PssmObject


class OrganismFactory:
    """Factory
    """

    def __init__(self, conf_org, conf_org_fac, conf_con, conf_pssm) -> None:

        self._id = 0
        
        self.num_recognizers_lambda_param = conf_org_fac[
            "NUM_RECOGNIZERS_LAMBDA_PARAM"
        ]

        self.min_mu = conf_org_fac["MIN_MU"]
        self.max_mu = conf_org_fac["MAX_MU"]

        self.min_sigma = conf_org_fac["MIN_SIGMA"]
        self.max_sigma = conf_org_fac["MAX_SIGMA"]

        self.pwm_length = conf_org_fac["PWM_LENGTH"]

        self.pwm_probability_step = conf_org_fac[
            "PWM_PROBABILITY_STEP"
        ]  # It should be a BASE_PROBABILITY divisor Ex: 1, 2, 4, 5, 10, 25...
        self.pwm_probability_base = conf_org_fac["PWM_PROBABILITY_BASE"]
        self.pwm_probability_decimals = conf_org_fac[
            "PWM_PROBABILITY_DECIMALS"
        ]

        self.conf_org = conf_org
        self.conf_con = conf_con
        self.conf_pssm = conf_pssm
        
        # Add ENERGY_THRESHOLD_METHOD and ENERGY_THRESHOLD_PARAM from
        # conf_org to conf_con and conf_pssm
        self.conf_con["ENERGY_THRESHOLD_METHOD"] = conf_org["ENERGY_THRESHOLD_METHOD"]
        self.conf_pssm["ENERGY_THRESHOLD_METHOD"] = conf_org["ENERGY_THRESHOLD_METHOD"]
        
        self.conf_con["ENERGY_THRESHOLD_PARAM"] = conf_org["ENERGY_THRESHOLD_PARAM"]
        self.conf_pssm["ENERGY_THRESHOLD_PARAM"] = conf_org["ENERGY_THRESHOLD_PARAM"]
        

    def get_id(self) -> int:
        """Gives a new ID for an organism
        TODO: This should be a function so all the count of IDs, including
        assigned outside the class, keep consistency between all organisms

        Returns:
           a new non-repeated ID
        """
        self._id += 1
        return self._id

    def get_organism(self) -> OrganismObject:
        """It creates and returns a full organism datastructure

        Returns:
            A new organism based on JSON config file
        """

        new_organism = OrganismObject(self.get_id(), self.conf_org, self.conf_pssm["MAX_COLUMNS"])
        
        # The number of recognizers of the organism is randomly chosen from a
        # Poisson distribution with the lambda provided.
        # The Poisson distribution is truncated to integers larger than 0, so
        # that null organisms are avoided. The truncation is achieved by using
        # lambda - 1 instead of lambda (in this way the average number of
        # recognizers will be lower by one unit) and then shifting up the
        # values by one unit.
        number_of_recognizers = np.random.poisson(self.num_recognizers_lambda_param - 1)
        number_of_recognizers += 1
        
        for i in range(number_of_recognizers - 1):
            # New recognizer
            new_recognizer = self.create_pssm(self.pwm_length)
            new_organism.recognizers.append(new_recognizer)
            # New connector
            _mu = random.randint(self.min_mu, self.max_mu)
            _sigma = random.randint(self.min_sigma, self.max_sigma)
            new_connector = ConnectorObject(_mu, _sigma, self.conf_con)
            new_organism.connectors.append(new_connector)
        new_recognizer = self.create_pssm(self.pwm_length)  # last recognizer
        new_organism.recognizers.append(new_recognizer)

        return new_organism

    def create_connection(
            self, connection_probability: float
    ) -> ConnectorObject:
        """It returns a connector object with its nodes assigned depending on
        connectionProbability

        Args:
            connection_probability: Probability to generate a connector instead
                                    of a recognizer

        Returns:
            A new Connection with recognizers included
        """

        # Assign a random value to mu and sigma
        _mu = random.randint(self.min_mu, self.max_mu)
        _sigma = random.randint(self.min_sigma, self.max_sigma)

        # Create the new connection
        new_connection = ConnectorObject(_mu, _sigma, self.conf_con)

        # Set the connection node to a connector or PSSM object, based on a
        # random probability. If connector object is selected, probability of
        # getting another connector is reduced by PROBABILITY_REDUCED_FACTOR
        node1 = None
        if random.random() < connection_probability:
            node1 = self.create_connection(
                connection_probability * self.reducer_probability_factor
            )
        else:
            node1 = self.create_pssm(self.pwm_length)

        new_connection.set_node1(node1)

        node2 = None
        if random.random() < connection_probability:
            node2 = self.create_connection(
                connection_probability * self.reducer_probability_factor
            )
        else:
            node2 = self.create_pssm(self.pwm_length)

        new_connection.set_node2(node2)

        return new_connection

    def create_pssm(self, length: int) -> PssmObject:
        """It return a PSSM object with a specific length

        Args:
            min_length: minimum positions a recognizer can recognize
            max_length: maximum positions a recognizer can recognize

        Returns:
            A pssm with an initializated PWM
        """
        pwm = []
        # Generate as much as needed
        for _ in range(length):
            pwm.append(self.get_pwm_column())

        return PssmObject(np.array(pwm), self.conf_pssm)

    def get_pwm_column(self) -> dict:
        """Generates a single column of the pwm

        Returns:
            a random probability for each base [a, c, g, t]
        """

        initial_probability = (
            self.pwm_probability_base / self.pwm_probability_step
        )
        probabilities: list = []
        # Left probability is
        left_probability = initial_probability
        # Minimum and maximum number of probabilities to be generated
        min_probability = 0
        max_probability = 4
        # Number of decimals on the probability
        decimals = 2
        # Generate 4 random probabilities out of initial_probability, one for
        # each base

        # Add a probability while we have less than 3 and and total probability
        # is not 1
        while (
                left_probability > min_probability
                and len(probabilities) < max_probability - 1
        ):
            new_probability = random.randint(0, left_probability)
            probabilities.append(float(new_probability))
            left_probability -= new_probability
        # Add the last probability or fill with 0 probability
        if left_probability > 0:
            probabilities.append(initial_probability - sum(probabilities))
        else:
            while len(probabilities) < max_probability:
                probabilities.append(0.0)

        # Shuffle the array is needed so high probability is not always on
        # first positions
        random.shuffle(probabilities)

        # Transform probabilities array from integer
        # [0-(BASE_PROBABILITY / STEP)] to complementary float
        # probabilities [0.0-1.0]
        np_probabilities = (
            np.array(probabilities)
            * self.pwm_probability_step
            * (1 / self.pwm_probability_base)
        )
        probabilities = np_probabilities.tolist()

        # Return object with "decimals" decimals probability to each base
        return {
            "a": round(probabilities[0], decimals),
            "g": round(probabilities[1], decimals),
            "c": round(probabilities[2], decimals),
            "t": round(probabilities[3], decimals),
        }
    
    def get_children(self, parent1, parent2):
        
        index_1 = parent1.get_random_connector()
        index_2 = parent2.get_random_connector()
        
        child1 = OrganismObject(self.get_id(), self.conf_org, self.conf_pssm["MAX_COLUMNS"])
        child2 = OrganismObject(self.get_id(), self.conf_org, self.conf_pssm["MAX_COLUMNS"])
        
        if random.random() < 0.5:
            
            # First parent keeps the broken connector in the left chunk
            L1, R1 = parent1.break_chain(index_1, "left")

            if random.random() < 0.5:
                # Second parent keeps the broken connector in the left chunk
                L2, R2 = parent2.break_chain(index_2, "left")
            
                # Recombination: case 1 (left, left)
                
                # Child 1 is (L1 + R2)
                child1_reco = L1["recognizers"] + R2["recognizers"]
                child1_conn = L1["connectors"] + R2["connectors"]
                # Child 2 is (L2 + R1)
                child2_reco = L2["recognizers"] + R1["recognizers"]
                child2_conn = L2["connectors"] + R1["connectors"]
                
            else:
                # Second parent keeps the broken connector in the right chunk
                L2, R2 = parent2.break_chain(index_2, "right")
            
                # Recombination: case 2 (left, right)
                
                # Child 1 is (L1 + L2)
                child1_reco = L1["recognizers"] + L2["recognizers"]
                child1_conn = L1["connectors"] + L2["connectors"]
                # Child 2 is (R1 + R2)
                child2_reco = R1["recognizers"] + R2["recognizers"]
                child2_conn = R1["connectors"] + R2["connectors"]
                
        else:
            
            # First parent keeps the broken connector in the right chunk
            L1, R1 = parent1.break_chain(index_1, "right")
            
            if random.random() < 0.5:
                # Second parent keeps the broken connector in the left chunk
                L2, R2 = parent2.break_chain(index_2, "left")
            
                # Recombination: case 3 (right, left)
                
                # Child 1 is (L2 + L1)
                child1_reco = L2["recognizers"] + L1["recognizers"]
                child1_conn = L2["connectors"] + L1["connectors"]
                # Child 2 is (R2 + R1)
                child2_reco = R2["recognizers"] + R1["recognizers"]
                child2_conn = R2["connectors"] + R1["connectors"]
                
            
            else:
                # Second parent keeps the broken connector in the right chunk
                L2, R2 = parent2.break_chain(index_2, "right")
            
                # Recombination: case 4 (right, right)
                
                # Child 1 is (L1 + R2)
                child1_reco = L1["recognizers"] + R2["recognizers"]
                child1_conn = L1["connectors"] + R2["connectors"]
                # Child 2 is (L2 + R1)
                child2_reco = L2["recognizers"] + R1["recognizers"]
                child2_conn = L2["connectors"] + R1["connectors"]
        
        # Set child1 recognizers and connectors
        child1.set_recognizers(child1_reco)
        child1.set_connectors(child1_conn)
        
        # Set child2 recognizers and connectors
        child2.set_recognizers(child2_reco)
        child2.set_connectors(child2_conn)
        
        print("\n")
        child1.print()
        print("\n")
        child2.print()
        
        return [child1, child2]

    def import_organisms(self, file_name: str) -> list:
        """Import Organisms from file

        Args:
            file_name: Name of the file with the organisms to read as an input

        Returns:
            a list of organisms objects read from the file
        """
        organism_list = []

        with open(file_name) as json_file:
            organism_json = json.load(json_file)

        for organism in organism_json:

            new_organism = OrganismObject(
                self.get_id(), self.conf_org, self.conf_pssm["MAX_COLUMNS"]
            )
            
            new_org_recognizers = []  # The recognizers are collected here
            new_org_connectors = []  # The connectors are collected here
            
            for element in organism:
                if element["objectType"] == "pssm":
                    new_org_recognizers.append(self.import_pssm(element))
                elif element["objectType"] == "connector":
                    new_org_connectors.append(self.import_connector(element))
            
            # Set recognizers and connectors of the organism
            new_organism.set_recognizers(new_org_recognizers)
            new_organism.set_connectors(new_org_connectors)

            #if "isTracked" in organism.keys():  # !!!
            #    new_organism.set_is_tracked(organism["isTracked"])

            organism_list.append(new_organism)

        return organism_list

    def import_connector(self, connector: dict) -> ConnectorObject:
        """Import Connector from JSON object

        Args:
            connector: connector in dictionary format

        Returns:
            Connector object from given connector dictionary
        """
        new_connector = ConnectorObject(
            connector["mu"], connector["sigma"], self.conf_con
        )

        return new_connector

    def import_pssm(self, pssm: dict) -> PssmObject:
        """Import PSSM from JSON object

        Args:
            pssm: pssm recognizer in dictionary format

        Returns:
            PSSM Object from given  pssm dictionary

        """
        return PssmObject(np.array(pssm["pwm"]), self.conf_pssm)

    def export_organisms(self, a_organisms: list, filename: str) -> None:
        """Export a list of organisms to JSON format

        Args:
            a_organisms: list of organisms to export
            filename: name of the file to export all the organisms
        """
        
        list_json_organisms = []
        for o_organism in a_organisms:
            organism = []
            for i in range(o_organism.count_recognizers() - 1):
                organism.append(self.export_pssm(o_organism.recognizers[i]))
                organism.append(self.export_connector(o_organism.connectors[i]))
            organism.append(self.export_pssm(o_organism.recognizers[-1]))
            list_json_organisms.append(organism)
        
        with open(filename, "w+") as json_file:
            json.dump(list_json_organisms, json_file, indent=2)

    def export_connector(self, o_connector: ConnectorObject) -> dict:
        """Export connector object

        Args:
            o_connector: Connector to export

        Returns:
            Connector in dictionary format
        """
        connector = {}
        connector["objectType"] = "connector"
        connector["mu"] = o_connector._mu
        connector["sigma"] = o_connector._sigma

        return connector

    def export_pssm(self, o_pssm: PssmObject) -> dict:
        """Export PSSM object

        Args:
            o_pssm: PSSM object to export

        Returns:
            pssm in dictionary format

        """
        pssm = {}
        pssm["objectType"] = "pssm"
        pssm["pwm"] = o_pssm.pwm.tolist()
        return pssm
