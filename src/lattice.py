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
        self.dim = protein.length
        if grid is None:
        self.grid = np.empty((self.dim, self.dim), dtype=object)
        self.grid[:] = None
        # The protein is initially unfolded and placed horizontally in the lattice
        i = self.dim // 2  # the middle of the lattice
        for j in range(self.dim):
            residue = self.protein.get_residue(
                j + 1
            )  # addition +1 because the aa numbers begins with 1
            residue.x_coord, residue.y_coord = j, i
            self.grid[i, j] = residue

    def __str__(self):
        """
        Return a string representation of the lattice grid.

        Returns
        -------
        str
            A formatted string depicting the lattice grid, with residues shown
            by their type and number, and empty cells represented by spaces.
        """
        description = ""
        for i in range(self.dim):
            line = "|"
            for j in range(self.dim):
                residue = self.grid[i, j]
                if residue is None:
                    line += "    |"
                else:
                    line += residue.type
                    line += str(residue.number)
                    line += "|"
            line += "\n"
            description += line
        return description

    def get_adjacents(self, pos_i, pos_j):
        """
        Retrieve the adjacent residues in the lattice.

        This method returns a list of residues adjacent to the given position
        in the lattice. It checks the four possible neighboring positions
        (up, down, left, right) and includes only the non-empty ones.

        Parameters
        ----------
        pos_i : int
            The row index of the residue in the lattice.
        pos_j : int
            The column index of the residue in the lattice.

        Returns
        -------
        list
            A list of Residue objects that are adjacent to the specified position.
            The list will be empty if there are no adjacent residues.
        """
        adjacents = []
        if pos_i + 1 < self.dim and self.grid[pos_i + 1, pos_j] is not None:
            adjacents.append(self.grid[pos_i + 1, pos_j])
        if pos_j + 1 < self.dim and self.grid[pos_i, pos_j + 1] is not None:
            adjacents.append(self.grid[pos_i, pos_j + 1])
        if pos_i - 1 > 0 and self.grid[pos_i - 1, pos_j] is not None:
            adjacents.append(self.grid[pos_i - 1, pos_j])
        if pos_j - 1 > 0 and self.grid[pos_i, pos_j - 1] is not None:
            adjacents.append(self.grid[pos_i, pos_j - 1])
        return adjacents

    def compute_energy(self):
        """
        Compute the energy of the protein configuration in the lattice.

        This method calculates the total energy of the protein by summing the
        interactions between adjacent hydrophobic (H) residues that are not
        sequentially connected in the protein sequence.

        The energy is decreased by 1 for each pair of adjacent hydrophobic
        residues that are not directly connected in the sequence.

        Returns
        -------
        int
            The computed energy of the protein configuration in the lattice
        """
        energy = 0
        for i in range(1, self.protein.length):
            residue1 = self.protein.get_residue(i)
            # not necessary to continue if the residue is not hydrophobic
            if residue1.type == "H":
                for j in range(i + 1, self.protein.length + 1):
                    residue2 = self.protein.get_residue(j)
                    # not necessary to continue if the neighbour is not hydrophobic
                    if (
                        residue2.type == "H"
                        and not residue1.is_connected(residue2)
                        and residue1.is_adjacent(residue2)
                    ):
                        energy -= 1
        return energy
