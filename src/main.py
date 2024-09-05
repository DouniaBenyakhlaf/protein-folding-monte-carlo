import numpy as np

# from residue import * # already imported in protein class
from protein import *
from lattice import *

if __name__ == "__main__":

    # # =================== Test of Residue class ===================
    # print("=================== Test of Residue class ===================")
    # arginine = Residue("P", 3, 4, 7)
    # glycine = Residue("H", 2, 4, 8)
    # alanine = Residue("H", 52, 10, 8)
    # print(arginine)
    # print(glycine)
    # print(alanine)
    # print(arginine.is_connected(glycine))
    # print(arginine.is_connected(alanine))

    # # =================== Test of Protein class ===================
    # print("=================== Test of Protein class ===================")
    # protein_1 = Protein("AASASS")
    # print(protein_1)

    # # =================== Test of Lattice class ===================
    # print("=================== Test of Lattice class ===================")
    # lattice_1 = Lattice(protein_1)
    # print(lattice_1)
    # print(lattice_1.compute_energy())
    # lattice_1.possible_endmoves()

    # print("--------test of energy function--------")
    # residues = [
    #     Residue("P", 1, 2, 2),
    #     Residue("P", 2, 1, 2),
    #     Residue("P", 3, 1, 3),
    #     Residue("H", 4, 2, 3),
    #     Residue("H", 5, 2, 4),
    #     Residue("P", 6, 1, 4),
    #     Residue("P", 7, 1, 5),
    #     Residue("P", 8, 2, 5),
    #     Residue("P", 9, 3, 5),
    #     Residue("P", 10, 4, 5),
    #     Residue("P", 11, 4, 4),
    #     Residue("H", 12, 3, 4),
    #     Residue("H", 13, 3, 3),
    #     Residue("P", 14, 4, 3),
    #     Residue("P", 15, 4, 2),
    #     Residue("P", 16, 3, 2),
    # ]

    # grid = np.empty((16 * 2, 16 * 2), dtype=object)
    # grid[:] = None
    # for residue in residues:
    #     grid[residue.i_coord, residue.j_coord] = residue

    # protein_2 = Protein("SSSGGSSSSSSGGSSS")
    # for i in range(protein_2.length):
    #     protein_2.hp_sequence[i] = residues[i]

    # lattice_2 = Lattice(protein_2, grid)  # sequence "PPPHHPPPPPPHHPPP"
    # print(lattice_2)
    # print(lattice_2.compute_energy())
    # lattice_2.possible_endmoves()
    # lattice_2.possible_cornermoves(residues[4])

    # print("--------test of corner moves--------")
    # residues_v2 = [
    #     Residue("P", 1, 1, 1),
    #     Residue("P", 2, 1, 0),
    #     Residue("P", 3, 2, 0),
    # ]
    # grid2 = np.empty((3, 3), dtype=object)
    # grid2[:] = None
    # for residue in residues_v2:
    #     grid2[residue.i_coord, residue.j_coord] = residue

    # protein_3 = Protein("PPP")
    # for i in range(protein_3.length):
    #     protein_3.hp_sequence[i] = residues_v2[i]
    # lattice_3 = Lattice(protein_3, grid2)
    # print(lattice_3)
    # lattice_3.possible_cornermoves(residues_v2[1])
    # lattice_3.possible_cornermoves(residues_v2[0])
    # lattice_3.possible_crankshaftmoves(residues_v2[1])

    print("--------test of crankshaf moves--------")
    residues_v3 = [
        Residue("P", 1, 0, 1),
        Residue("P", 2, 1, 1),
        Residue("P", 3, 1, 0),
        Residue("P", 4, 2, 0),
        Residue("P", 5, 2, 1),
        Residue("P", 6, 3, 1),
    ]
    grid3 = np.empty((3 * 6, 3 * 6), dtype=object)
    grid3[:] = None
    for residue in residues_v3:
        grid3[residue.i_coord, residue.j_coord] = residue

    protein_4 = Protein("PPPPPP")
    for i in range(protein_4.length):
        protein_4.hp_sequence[i] = residues_v3[i]
    lattice_4 = Lattice(protein_4, grid3)
    print(lattice_4)
    lattice_4.possible_cornermoves(residues_v3[1])
    lattice_4.possible_crankshaftmoves(residues_v3[2])
    lattice_4.possible_crankshaftmoves(residues_v3[3])
