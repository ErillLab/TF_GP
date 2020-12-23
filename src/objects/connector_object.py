"""C object
Connects two nodes at a specific distance

"""

# pylint: disable=E0402
# type: ignore
import random
import numpy as np
import math


def norm_pdf(x, mu, sigma):
    if sigma != 0:
        var = float(sigma)**2
        denom = (2*math.pi*var)**.5
        num = math.exp(-(float(x)-float(mu))**2/(2*var))
        p = num/denom
    else:
        # when sigma is 0
        if x == mu:
            p = 1
        else:
            p = 0
    return p


def norm_cdf(x, mu, sigma):
    # Cumulative distribution function for the normal distribution
    z = (x-mu)/abs(sigma)
    return (1.0 + math.erf(z / math.sqrt(2.0))) / 2.0



# pylint: disable=R0902
class ConnectorObject():
    """Connector Object is a node that connects two nodes
    """

    # pylint: disable=R0913
    def __init__(
            self,
            _mu: int,
            _sigma: int,
            config: dict,
    ):
        """Connector constructor gets mu, sigma and can get the two initial nodes.

        Args:
            _mu: Mean distance between node1 and node2
            _sigma: Variance between node 1 and node2
            config: Configurations loadad from config.json

        """
        super().__init__()
        self._mu = _mu  # Mean discance
        self._sigma = _sigma  # Variance between elements
        self.mutate_probability_sigma = config["MUTATE_PROBABILITY_SIGMA"]
        self.mutate_probability_mu = config["MUTATE_PROBABILITY_MU"]
        self.mutate_variance_sigma = config["MUTATE_VARIANCE_SIGMA"]
        self.mutate_variance_mu = config["MUTATE_VARIANCE_MU"]
        self.placement_options = config["PLACEMENT_OPTIONS"]
        self.energy_threshold_method = config["ENERGY_THRESHOLD_METHOD"]
        self.energy_threshold_value = config["ENERGY_THRESHOLD_PARAM"]

    # pylint: enable=R0913
    # Setters
    def set_mu(self, _mu: int) -> None:
        """Set mu variable

        Args:
            _mu: Mean distance between node1 and node2
        """
        self._mu = _mu

    def set_sigma(self, sigma: int) -> None:
        """Set sigma variable

        Args:
            sigma: Variance between node 1 and node2
        """
        self._sigma = sigma

    # pylint: disable=W0613
    def mutate(self, org_factory) -> None:
        """mutation for a connector

        Args:
            org_factory(organism_factory): Organism Facory
        """
        
        
        
        # LINEAR SIGMA MUTATION
        if random.random() < self.mutate_probability_sigma:
            # Alter sigma
            self._sigma += random.uniform(
                -self.mutate_variance_sigma, self.mutate_variance_sigma
            )
        
        
        '''
        # LOGARITHMIC SIGMA MUTATION
        if random.random() < self.mutate_probability_sigma:
            base = self.mutate_variance_sigma
            logb_sigma = np.log(self._sigma) / np.log(base)
            shift = random.uniform(-1, 1)
            # Apply a shift in the range (-1, 1) to the log-sigma
            logb_sigma += shift
            new_sigma = base**logb_sigma
            self._sigma = new_sigma
        '''

        if random.random() < self.mutate_probability_mu:
            # Alter mu
            self._mu += random.uniform(
                -self.mutate_variance_mu, self.mutate_variance_mu
            )

    # pylint: enable=W0613

    def print(self) -> None:
        """It prints the connector mu, sigma values and its children values in
           tree structure
        """
        print(" m: {} s: {}".format(self._mu, self._sigma))

    def export(self, export_file) -> None:
        """Exports Connector data to the given file

        Args:
            export_file (file): File to export the conector

        """
        export_file.write("\n m: {} s: {}".format(self._mu, self._sigma))

    # pylint: disable=R0201
    def is_connector(self) -> bool:
        """node is connector

        Returns:
            True because is a connector
        """
        return True

    def is_pssm(self) -> bool:
        """node is not a pssm

        Returns:
            False because is a connector
        """
        return False
