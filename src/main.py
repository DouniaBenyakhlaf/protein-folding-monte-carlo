# from residue import * # already imported in protein class
from protein import *
from lattice import *

if __name__ == "__main__":

    # =================== Test of Residue class ===================
    print("=================== Test of Residue class ===================")
    arginine = Residue("P", 3, 4, 7)
    glycine = Residue("H", 2, 4, 8)
    alanine = Residue("H", 52, 10, 8)
    arginine.res_description()
    glycine.res_description()
    alanine.res_description()
    print(arginine.is_connected(glycine))
    print(arginine.is_connected(alanine))
    alanine.set_coordinate(11, 11)
    alanine.res_description()

    # =================== Test of Protein class ===================
    print("=================== Test of Protein class ===================")
    protein_1 = Protein("HHPHPP")
    protein_1.print_sequence()

    # =================== Test of Lattice class ===================
    print("=================== Test of Lattice class ===================")
    lattice_1 = Lattice(protein_1)
    lattice_1.print_lattice()
