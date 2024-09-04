# from residue import * # already imported in protein class
from protein import *
from lattice import *

if __name__ == "__main__":

    # =================== Test of Residue class ===================
    print("=================== Test of Residue class ===================")
    arginine = Residue("P", 3, 4, 7)
    glycine = Residue("H", 2, 4, 8)
    alanine = Residue("H", 52, 10, 8)
    print(arginine)
    print(glycine)
    print(alanine)
    print(arginine.is_connected(glycine))
    print(arginine.is_connected(alanine))

    # =================== Test of Protein class ===================
    print("=================== Test of Protein class ===================")
    protein_1 = Protein("HHPHPP")
    print(protein_1)

    # =================== Test of Lattice class ===================
    print("=================== Test of Lattice class ===================")
    lattice_1 = Lattice(protein_1)
    print(lattice_1)
    print(lattice_1.compute_energy())
