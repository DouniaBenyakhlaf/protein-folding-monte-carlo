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

    def __init__(self, protein, grid=None):
        """
        Initialize a Lattice instance with a protein.

        Parameters
        ----------
        protein : Protein
            The protein that will be contained within the lattice.
        """
        self.protein = protein
        if grid is None:
            self.dim = protein.length * 4  # dimension of the lattice
            print(self.dim)
            self.grid = np.empty((self.dim, self.dim), dtype=object)
            self.grid[:] = None
            # The protein is initially unfolded and placed horizontally in the lattice
            i = self.dim // 2  # the middle of the lattice
            j = self.dim // 2
            protein_position = 0
            while protein_position < self.protein.length:
                residue = self.protein.get_residue(
                    protein_position + 1
                )  # addition +1 because the aa numbers begins with 1
                residue.i_coord, residue.j_coord = i, j
                self.grid[i, j] = residue
                protein_position += 1
                j += 1
        else:
            self.grid = grid  # in order to test some conformation (energy function etc). Remove this part later.
            self.dim = grid.shape[0]

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
                    line += f"{residue.number:< 3d}"
                    line += "|"
            line += "\n"
            description += line
        return description

    def verify_dim(self, position):
        """
        Check if a position is within the bounds of the grid dimensions.

        This method verifies whether the given position is within the valid range of the grid dimensions.
        Specifically, it checks if both the x and y coordinates of the position are greater than 0 and less
        than the grid's dimension size.

        Parameters
        ----------
        position : tuple of int
            A tuple (x, y) representing the coordinates of the position to be checked.

        Returns
        -------
        bool
            True if the position is within the grid boundaries, otherwise False.
        """
        return (
            position[0] > 0
            and position[0] < self.dim
            and position[1] > 0
            and position[1] < self.dim
        )

    def get_adjacents_available_positions(self, residue):
        """
        Retrieve available adjacent positions in the lattice grid.

        This method checks the positions directly adjacent to the given
        coordinates `(pos_i, pos_j)` in the grid. It returns a list of
        coordinates where adjacent positions are empty (i.e., `None`).

        Parameters
        ----------
        pos_i : int
            The row index of the current position in the lattice.
        pos_j : int
            The column index of the current position in the lattice.

        Returns
        -------
        list of tuple of (int, int)
            A list of tuples, where each tuple represents the coordinates
            of an adjacent position that is available (i.e., empty).
        """
        pos_i, pos_j = residue.i_coord, residue.j_coord
        available_positions = []
        if pos_i + 1 < self.dim and self.grid[pos_i + 1, pos_j] is None:
            available_positions.append((pos_i + 1, pos_j))
        if pos_j + 1 < self.dim and self.grid[pos_i, pos_j + 1] is None:
            available_positions.append((pos_i, pos_j + 1))
        if pos_i - 1 >= 0 and self.grid[pos_i - 1, pos_j] is None:
            available_positions.append((pos_i - 1, pos_j))
        if pos_j - 1 >= 0 and self.grid[pos_i, pos_j - 1] is None:
            available_positions.append((pos_i, pos_j - 1))
        return available_positions

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

    def possible_endmoves(self):
        first_neighbour = self.protein.get_residue(2)
        last_neighbour = self.protein.get_residue(self.protein.length - 1)
        first_available_positions = self.get_adjacents_available_positions(
            first_neighbour
        )
        last_available_positions = self.get_adjacents_available_positions(
            last_neighbour
        )
        if len(first_available_positions) > 0:
            print(
                f"Available positions for the first residue : {first_available_positions}"
            )
        if len(last_available_positions) > 0:
            print(
                f"Available positions for the last residue : {last_available_positions}"
            )

    def possible_cornermoves(self, residue):
        if residue.number != self.protein.length and residue.number != 1:
            connected_neighbour_1 = self.protein.get_residue(
                (residue.number + 1)
            )
            connected_neighbour_2 = self.protein.get_residue(
                (residue.number - 1)
            )
            avail_pos_neighbour_1 = set(
                self.get_adjacents_available_positions(connected_neighbour_1)
            )
            avail_pos_neighbour_2 = set(
                self.get_adjacents_available_positions(connected_neighbour_2)
            )
            cornermove = list(
                avail_pos_neighbour_1.intersection(avail_pos_neighbour_2)
            )
            print(f"available cornermove position: {cornermove}")
        else:
            print("corner moves are not possible on end residues")

    def crankshaft_first_corner(self, residue):
        if self.protein.is_residue_in_ushaped_bend(residue, True):
            if (residue.number - 1) >= 1 and (
                residue.number + 2
            ) <= self.protein.length:
                res_i_minus_1 = self.protein.get_residue(residue.number - 1)
                res_i_plus_1 = self.protein.get_residue(residue.number + 1)
                res_i_plus_2 = self.protein.get_residue(residue.number + 2)
                sym_pos_i = residue.get_symetrical_position(res_i_minus_1)
                sym_pos_i_plus_1 = res_i_plus_1.get_symetrical_position(
                    res_i_plus_2
                )
                if (
                    self.verify_dim(sym_pos_i)
                    and self.grid[sym_pos_i[0], sym_pos_i[1]] is None
                ) and (
                    self.verify_dim(sym_pos_i_plus_1)
                    and self.grid[sym_pos_i_plus_1[0], sym_pos_i_plus_1[1]]
                    is None
                ):
                    print(
                        f"The position {sym_pos_i} is available for the residue {residue.number}"
                    )
                    print(
                        f"The position {sym_pos_i_plus_1} is available for the residue {residue.number + 1}"
                    )
                    return True
        return False

    def crankshaft_second_corner(self, residue):
        if self.protein.is_residue_in_ushaped_bend(residue, False):
            if (residue.number - 2) >= 1 and (
                residue.number + 1
            ) <= self.protein.length:
                res_i_minus_2 = self.protein.get_residue(residue.number - 2)
                res_i_minus_1 = self.protein.get_residue(residue.number - 1)
                res_i_plus_1 = self.protein.get_residue(residue.number + 1)
                sym_pos_i = residue.get_symetrical_position(res_i_plus_1)
                sym_pos_i_minus_1 = res_i_minus_1.get_symetrical_position(
                    res_i_minus_2
                )
                if (
                    self.verify_dim(sym_pos_i)
                    and self.grid[sym_pos_i[0], sym_pos_i[1]] is None
                ) and (
                    self.verify_dim(sym_pos_i_minus_1)
                    and self.grid[sym_pos_i_minus_1[0], sym_pos_i_minus_1[1]]
                    is None
                ):
                    print(
                        f"The position {sym_pos_i} is available for the residue {residue.number}"
                    )
                    print(
                        f"The position {sym_pos_i_minus_1} is available for the residue {residue.number-1}"
                    )
                    return True
        return False

    def possible_crankshaftmoves(self, residue):
        if self.crankshaft_first_corner(
            residue
        ) or self.crankshaft_second_corner(residue):
            print("Possible")
        else:
            print("Not possible")
