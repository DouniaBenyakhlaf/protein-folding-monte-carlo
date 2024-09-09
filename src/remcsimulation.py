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

    def swap_labels(self, temperature1, temperature2):
        previous_replica1 = self.replicas[temperature1]
        self.replicas[temperature1] = self.replicas[temperature2]
        self.replicas[temperature2] = previous_replica1

    def mcsearch(self, current_lattice, temperature):
        for i in range(self.nb_steps):
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
                else:
                    random_prob = random.random()
                    if random_prob > math.exp(-delta_energy / temperature):
                        return new_lattice
        return current_lattice

    def run(self):
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
                if new_energy < energy:
                    energy = new_energy
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
