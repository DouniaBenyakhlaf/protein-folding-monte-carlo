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
    # lattice_1.end_moves()

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
    # print(lattice_2.end_moves(residues[0]))
    # print(lattice_2.corner_moves(residues[4]))

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
    # print(lattice_3.corner_moves(residues_v2[1]))
    # print(lattice_3.corner_moves(residues_v2[1]).corner_moves(residues_v2[1]))
    # print(lattice_3.crankshaft_moves(residues_v2[1]))

    # print("--------test of crankshaf moves--------")
    # residues_v3 = [
    #     Residue("P", 1, 0, 1),
    #     Residue("P", 2, 1, 1),
    #     Residue("P", 3, 1, 0),
    #     Residue("P", 4, 2, 0),
    #     Residue("P", 5, 2, 1),
    #     Residue("P", 6, 3, 1),
    # ]
    # grid3 = np.empty((3 * 6, 3 * 6), dtype=object)
    # grid3[:] = None
    # for residue in residues_v3:
    #     grid3[residue.i_coord, residue.j_coord] = residue

    # protein_4 = Protein("PPPPPP")
    # for i in range(protein_4.length):
    #     protein_4.hp_sequence[i] = residues_v3[i]
    # lattice_4 = Lattice(protein_4, grid3)
    # print(lattice_4)
    # print(lattice_4.end_moves(residues_v3[0]))
    # print(lattice_4.corner_moves(residues_v3[1]))
    # print(lattice_4.crankshaft_moves(residues_v3[2]))
    # print(lattice_4.crankshaft_moves(residues_v3[3]))

    print("--------test of pull moves--------")
    # test is_valid :
    # protein_5 = Protein("PPPPPP")
    # lattice_5 = Lattice(protein_5)
    # print(protein_5)
    # print(lattice_5)
    # print(protein_5.is_sequence_valid())
    # protein_5.hp_sequence[3].i_coord = 56
    # print(protein_5.is_sequence_valid())

    residues_v3 = [
        Residue("H", 1, 3, 0),
        Residue("H", 2, 2, 0),
        Residue("H", 3, 1, 0),
        Residue("H", 4, 1, 1),
        Residue("H", 5, 0, 1),
        Residue("H", 6, 0, 2),
        Residue("H", 7, 1, 2),
        Residue("H", 8, 2, 2),
        Residue("H", 9, 3, 2),
    ]
    grid3 = np.empty((3 * 9, 3 * 9), dtype=object)
    grid3[:] = None
    for residue in residues_v3:
        grid3[residue.i_coord, residue.j_coord] = residue

    protein_4 = Protein("PPPPPPPPP")
    for i in range(protein_4.length):
        protein_4.hp_sequence[i] = residues_v3[i]
    lattice_4 = Lattice(protein_4, grid3)
    print(lattice_4)
    new_lattic_4 = lattice_4.pull_moves(residues_v3[7])
    print(new_lattic_4)
    print(new_lattic_4.pull_moves(new_lattic_4.protein.get_residue(1)))

    # test copy()
    # lattice_6 = lattice_4.copy()
    # print(lattice_6)
