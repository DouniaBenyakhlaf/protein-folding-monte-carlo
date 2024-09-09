"""Module for Replica exchange Monte Carlo search."""

import math
import random
from protein import Protein
from lattice import Lattice


class REMCSimulation:
    """
    Perform a Replica Exchange Monte Carlo (REMC) simulation.

    This class simulates the folding of a protein sequence using the
    Replica Exchange Monte Carlo method. The simulation maintains several
    replicas of the protein structure, each running at different temperature
    settings, allowing for efficient exploration of the conformational space.

    Attributes
    ----------
    aa_sequence : str
        The amino acid sequence of the protein being simulated.
    nb_steps : int
        The number of local Monte Carlo steps.
    tmp_min : int
        The minimum temperature value for the replicas.
    tmp_max : int
        The maximum temperature value for the replicas.
    nb_replicas : int
        The number of replicas in the simulation.
    optimal_energy : float
        The energy goal of the simulation.
    max_iter : int
        The maximum number of iterations allowed.
    replicas : dict
        A dictionary where each key is a temperature and the corresponding
        value is a replica of the lattice (protein structure) at that
        temperature.
    """

    def __init__(
        self,
        aa_sequence,
        nb_steps=500,
        tmp_min=160,
        tmp_max=220,
        nb_replicas=5,
        optimal_energy=-48,
        max_iter=50000,
    ):
        """
        Initialize a REMCSimulation instance.

        Parameters
        ----------
        aa_sequence : str
            The amino acid sequence of the protein to be simulated.
        nb_steps : int, optional
            Number of local Monte Carlo search steps per replica
            (default is 500).
        tmp_min : int, optional
            The minimum temperature for the simulation (default is 160).
        tmp_max : int, optional
            The maximum temperature for the simulation (default is 220).
        nb_replicas : int, optional
            The number of replicas to simulate, each with a unique temperature
            (default is 5).
        optimal_energy : float, optional
            The target energy value to be reached during the simulation
            (default is -48).
        max_iter : int, optional
            The maximum number of iterations allowed for the simulation
            (default is 50000).
        """
        self.aa_sequence = aa_sequence
        self.nb_steps = nb_steps
        self.tmp_min = tmp_min
        self.tmp_max = tmp_max
        self.nb_replicas = nb_replicas
        self.optimal_energy = optimal_energy
        self.max_iter = max_iter
        self.replicas = {}
        lattice = Lattice(Protein(self.aa_sequence))
        temperatures = self.linear_distribution_temperature(
            self.tmp_min, self.tmp_max
        )
        for tmp in temperatures:
            replica = lattice.copy()
            self.replicas[tmp] = replica

    def linear_distribution_temperature(self, tmp_min, tmp_max):
        """
        Generate a linearly distributed list of temperatures.

        This method creates a list of temperatures that are evenly
        spaced between the specified minimum and maximum temperature
        values. The number of temperatures generated corresponds to
        the number of replicas in the simulation.

        Parameters
        ----------
        tmp_min : int
            The minimum temperature value.
        tmp_max : int
            The maximum temperature value.

        Returns
        -------
        list of int
            A list of temperatures, linearly distributed between
            tmp_min and tmp_max, with a total length equal to
            self.nb_replicas.
        """
        temperatures = []
        for i in range(self.nb_replicas):
            tmp = int(
                tmp_min
                + ((i) / ((self.nb_replicas - 1)) * (tmp_max - tmp_min))
            )
            temperatures.append(tmp)
        return temperatures

    def swap_labels(self, temperature1, temperature2):
        """
        Swap the replicas associated with two temperatures.

        This method exchanges the replicas (protein structures) between two
        specified temperature values in the simulation. The replicas at
        temperature1 and temperature2 are swapped in the replicas dictionary.

        Parameters
        ----------
        temperature1 : int or float
            The first temperature whose associated replica will be swapped.
        temperature2 : int or float
            The second temperature whose associated replica will be swapped.
        """
        previous_replica1 = self.replicas[temperature1]
        self.replicas[temperature1] = self.replicas[temperature2]
        self.replicas[temperature2] = previous_replica1

    def mcsearch(self, current_lattice, temperature):
        """
        Perform a Monte Carlo search for a potential energy-minimizing move.

        This method applies a Monte Carlo search algorithm to the given lattice
        configuration. It randomly selects residues and attempts to apply a
        move (e.g., pull, crankshaft, or corner move) to minimize the energy
        of the system. The move is accepted based on the Metropolis criterion,
        which considers the energy difference and temperature.

        Parameters
        ----------
        current_lattice : Lattice
            The current lattice configuration on which the Monte Carlo
            search is performed.
        temperature : float
            The temperature value associated with the replica, used to
            compute the acceptance probability for energy-increasing moves.

        Returns
        -------
        Lattice
            The updated lattice configuration if a move is accepted,
            or the original lattice if no favorable move is found.
        """
        for _ in range(self.nb_steps):
            protein_length = current_lattice.protein.length
            residue_number = random.randint(1, protein_length)
            residue = current_lattice.protein.get_residue(residue_number)
            new_lattice = current_lattice.random_move(residue)
            if new_lattice is not None:
                current_energy = current_lattice.compute_energy()
                new_energy = new_lattice.compute_energy()
                delta_energy = new_energy - current_energy
                if delta_energy <= 0:
                    return new_lattice
                random_prob = random.random()
                if random_prob > math.exp(-delta_energy / temperature):
                    return new_lattice
        return current_lattice

    def run(self):
        """
        Execute the Replica Exchange Monte Carlo (REMC) simulation.

        This method runs the REMC simulation, performing Monte Carlo
        searches on each replica to minimize energy, and periodically
        attempting to exchange replicas between adjacent temperatures
        based on the Metropolis criterion.

        The process is repeated until the optimal energy is reached or
        the maximum number of iterations is exceeded. The energy of the
        system is printed periodically to track progress.
        """
        energy = 0
        offset = 0
        nb_iterations = 0
        while energy > self.optimal_energy and nb_iterations < self.max_iter:
            if nb_iterations % 500 == 0:
                print(
                    f"===============Iteration {nb_iterations}==============="
                )
            for temperature, replica in self.replicas.items():
                self.replicas[temperature] = self.mcsearch(
                    replica, temperature
                )
                new_energy = self.replicas[temperature].compute_energy()
                energy = min(energy, new_energy)
            i = offset
            while i + 1 < len(self.replicas):
                j = i + 1
                tmp_i = list(self.replicas.keys())[i]
                tmp_j = list(self.replicas.keys())[j]
                replica_i = self.replicas[tmp_i]
                replica_j = self.replicas[tmp_j]
                delta = (1 / tmp_j - 1 / tmp_i) * (
                    replica_i.compute_energy() - replica_j.compute_energy()
                )
                if delta <= 0:
                    self.swap_labels(tmp_i, tmp_j)
                else:
                    prob = random.random()
                    if prob <= math.exp(-delta):
                        self.swap_labels(tmp_i, tmp_j)
                i = i + 2
            offset = 1 - offset
            nb_iterations += 1
            if nb_iterations % 500 == 0:
                print(f"Energy_{nb_iterations} = {energy}")
        print(energy)


def expand_sequence(seq):
    """
    Expand a sequence of letters and numbers, replacing letters.

    This function processes a string where each letter is followed by a number
    representing the number of times the letter should be repeated. If no
    number follows a letter, it is repeated once by default. After expansion,
    the function replaces all occurrences of 'P' with 'R' and 'H' with 'G'.

    Parameters
    ----------
    seq : str
        The input string consisting of letters and optional numbers.

    Returns
    -------
    str
        A string where each letter is repeated according to the
        associated number, with 'P' replaced by 'R' and 'H'
        replaced by 'G'.

    Examples
    --------
    >>> expand_sequence("P2H3PH8")
    'RRRGGGGRRRRRRRR'

    >>> expand_sequence("P3H10PHP")
    'RRRGGGGGGGGGGRR'
    """
    result = ""
    i = 0
    while i < len(seq):
        if seq[i].isalpha():  # Vérifier si c'est une lettre
            letter = seq[i]
            i += 1
            num_str = ""
            while (
                i < len(seq) and seq[i].isdigit()
            ):  # Collecter le nombre associé
                num_str += seq[i]
                i += 1
            num = (
                int(num_str) if num_str else 1
            )  # Si pas de nombre, utiliser 1
            result += letter * num  # Répéter la lettre
        else:
            i += 1  # Ignorer les espaces
    return result.replace("P", "R").replace("H", "G")


if __name__ == "__main__":

    # Test sequence 1 : HPHPPHHPHPPHPHHPPHPH <=> "GRGRRGGRGRRGRGGRRGRG" # 9
    print(expand_sequence("P2H3PH8P3H10PHP3H12P4H6PH2PHP"))
    # AA_SEQUENCE = expand_sequence("P2H3PH8P3H10PHP3H12P4H6PH2PHP") # 32
    AA_SEQUENCE = expand_sequence(
        "P6HPH2P5H3PH5PH2P4H2P2H2PH5PH10PH2PH7P11H7P2HPH3P6HPH2"
    )
    REMC = REMCSimulation(AA_SEQUENCE)
    # print(REMC.replicas)
    # REMC.swap_labels(
    #     list(REMC.replicas.keys())[2], list(REMC.replicas.keys())[4]
    # )
    # print()
    # print(REMC.replicas)
    REMC.run()
