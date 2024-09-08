"""Module for Replica exchange Monte Carlo search."""

import math
import random
from protein import Protein
from lattice import Lattice


class REMCSimulation:

    def __init__(
        self,
        aa_sequence,
        nb_steps=500,
        tmp_min=160,
        tmp_max=220,
        nb_replicas=5,
        optimal_energy=-5,
    ):
        self.aa_sequence = aa_sequence
        self.nb_steps = nb_steps
        self.tmp_min = tmp_min
        self.tmp_max = tmp_max
        self.nb_replicas = nb_replicas
        self.optimal_energy = optimal_energy
        self.replicas = {}
        lattice = Lattice(Protein(self.aa_sequence))
        temperatures = self.linear_distribution_temperature(
            self.tmp_min, self.tmp_max
        )
        for tmp in temperatures:
            replica = lattice.copy()
            self.replicas[replica] = tmp

    def linear_distribution_temperature(self, tmp_min, tmp_max):
        temperatures = []
        for i in range(self.nb_replicas):
            tmp = int(
                tmp_min
                + ((i) / ((self.nb_replicas - 1)) * (tmp_max - tmp_min))
            )
            temperatures.append(tmp)
        return temperatures

    def swap_labels(self, replica1, replica2):
        previous_temp_rep1 = self.replicas[replica1]
        self.replicas[replica1] = self.replicas[replica2]
        self.replicas[replica2] = previous_temp_rep1

    def mcsearch(self, current_lattice):
        for i in range(self.nb_steps):
            protein_length = current_lattice.protein.length
            residue_number = random.randint(1, protein_length)
            residue = self.current_lattice.protein.get_residue(residue_number)
            move = self.current_lattice.random_move()
            new_lattice = move(residue)
            if new_lattice is not None:
                current_energy = self.current_lattice.compute_energy()
                new_energy = new_lattice.compute_energy()
                delta_energy = new_energy - current_energy
                if delta_energy <= 0:
                    return new_lattice
                else:
                    random_prob = random.random()
                    if random_prob > math.exp(
                        -delta_energy / self.temperature
                    ):
                        return new_lattice
        return current_lattice

    def run(self):
        pass


if __name__ == "__main__":
    AA_SEQUENCE = "MLSLGLLLLGLLQGVGKH"
    REMC = REMCSimulation(AA_SEQUENCE)
    print(REMC.replicas)
    REMC.swap_labels(
        list(REMC.replicas.keys())[2], list(REMC.replicas.keys())[4]
    )
    print()
    print(REMC.replicas)
