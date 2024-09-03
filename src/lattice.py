import numpy as np
from protein import *


class Lattice:
    """
    A class representing a 2D square lattice.

    Attributes
    ----------
    grid : numpy.ndarray
        An (n+2,n+2) numpy array where each cell contains a residue of a protein of size n.
    dim : int
        The dimension of the numpy array.
    protein : Protein
        The protein contained within the lattice.
    """

    def __init__(self, protein):
        """
        Initialize a Lattice instance with a protein.

        Parameters
        ----------
        protein : Protein
            The protein that will be contained within the lattice.
        """
        self.protein = protein
        # Simplification: protein.length + 2 to ensure that a residue always has 4 neighboring cells.
        self.dim = protein.get_length() + 2
        self.grid = np.empty((self.dim, self.dim), dtype=object)
        self.grid[:] = None
        # The protein is initially unfolded and placed horizontally in the lattice
        i = self.dim // 2  # the middle of the lattice
        for j in range(1, self.dim - 1):
            residue = self.protein.get_residue(j)
            residue.set_coordinate(j, i)
            self.grid[i, j] = residue

    def print_lattice(self):
        """
        Print a visual representation of the lattice grid.
        """
        # Rewrite all the print as seen in script course
        for i in range(self.dim):
            line = "|"
            for j in range(self.dim):
                residue = self.grid[i, j]
                if residue == None:
                    line += "  |"
                else:
                    line += residue.get_type()
                    line += str(residue.get_number())
                    line += "|"
            print(line)
