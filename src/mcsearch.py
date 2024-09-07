"""Module for performing a simple Monte Carlo search."""

import math
import random
from protein import Protein
from lattice import Lattice


class MCSearch:
    """
    Performs a MCsearch to optimize protein folding within a lattice.

    Attributes
    ----------
    current_lattice : Lattice
        The current lattice configuration containing the protein.
    nb_steps : int
        The number of Monte Carlo steps to be performed.
    temperature : float
        The temperature parameter used to determine the probability
        of accepting higher-energy configurations.
    """

    def __init__(self, lattice, nb_steps, temperature):
        """
        Initialize the MCSearch.

        Parameters
        ----------
        lattice : Lattice
            The lattice object containing the protein to optimize.
        nb_steps : int
            The total number of Monte Carlo steps to perform during the search.
        temperature : float
            The temperature used in the Metropolis criterion to accept
            or reject new configurations.
        """
        self.current_lattice = lattice
        self.nb_steps = nb_steps
        self.temperature = temperature

    def run(self):
        """
        Execute the Monte Carlo search to optimize the protein folding.

        The method runs through a set number of Monte Carlo steps,
        attempting random moves on random residues of the protein.
        The new lattice configuration is accepted if it has a lower
        or equal energy, or probabilistically if it has a higher energy
        based on the temperature.
        """
        for i in range(self.nb_steps):
            if i % 10 == 0:
                print(f"=============== Step {i} ===============")
            protein_length = self.current_lattice.protein.length
            residue_number = random.randint(1, protein_length)
            residue = self.current_lattice.protein.get_residue(residue_number)
            move = self.current_lattice.random_move()
            new_lattice = move(residue)
            new_lattice = self.current_lattice.pull_moves(residue)
            if new_lattice is not None:
                current_energy = self.current_lattice.compute_energy()
                new_energy = new_lattice.compute_energy()
                delta_energy = new_energy - current_energy
                if delta_energy <= 0:
                    self.current_lattice = new_lattice
                else:
                    random_prob = random.random()
                    if random_prob > math.exp(
                        -delta_energy / self.temperature
                    ):
                        self.current_lattice = new_lattice
        print(self.current_lattice)
        print(f"Energie = {self.current_lattice.compute_energy()}")


if __name__ == "__main__":
    AA_SEQUENCE = "MLSLGLLLLGLLQGVGKH"
    MCSEARCH = MCSearch(Lattice(Protein(AA_SEQUENCE)), 500, 180)
    MCSEARCH.run()
