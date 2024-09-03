from residue import *

if __name__ == "__main__":

    # =================== Test of Residue class ===================
    print("=================== Test of Residue class ===================")
    arginine = Residue("R", 3, "P", 4, 7)
    glycine = Residue("G", 2, "H", 4, 8)
    alanine = Residue("A", 52, "H", 10, 8)
    arginine.res_description()
    glycine.res_description()
    alanine.res_description()
    print(arginine.is_connected(glycine))
    print(arginine.is_connected(alanine))
    alanine.set_coordinate(11, 11)
    alanine.res_description()
