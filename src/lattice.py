"""Module for handling lattice."""

import numpy as np
from protein import Protein


class Lattice:
    """
    A class representing a 2D square lattice.

    Attributes
    ----------
    grid : numpy.ndarray
        An (n+2,n+2) numpy array where each cell
        contains a residue of a protein of size n.
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
            self.grid = np.empty((self.dim, self.dim), dtype=object)
            self.grid[:] = None
            # The protein is initially unfolded and placed horizontally
            # in the lattice
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
            # in order to test some conformation (energy function etc).
            # Remove this part later.
            self.grid = grid
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

        This method verifies whether the given position is within the valid
        range of the grid dimensions. Specifically, it checks if both the x
        and y coordinates of the position are greater than 0 and less than
        the grid's dimension size.

        Parameters
        ----------
        position : tuple of int
            A tuple (x, y) representing the coordinates of the position
            to be checked.

        Returns
        -------
        bool
            True if the position is within the grid boundaries,
            otherwise False.
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

    def move_residue(self, number_residue, new_position):
        """
        Move a residue to a new position on the grid.

        This method updates the grid by removing the residue from its current
        position and placing it at the new coordinates specified. It also
        updates the `i_coord` and `j_coord` attributes of the `residue`
        to reflect its new position.

        Parameters
        ----------
        number_residue : int
            The number of the residue to be repositioned.
        new_position : tuple of int
            A tuple (i, j) representing the new coordinates on the grid where
            the residue will be placed.
        """
        residue = self.protein.get_residue(number_residue)
        self.grid[residue.i_coord, residue.j_coord] = None
        residue.i_coord = new_position[0]
        residue.j_coord = new_position[1]
        self.grid[new_position[0], new_position[1]] = residue

    def copy(self):
        """
        Create a copy of the current lattice and its associated protein.

        This method creates a deep copy of the protein and lattice.
        It initializes a new protein object with the same amino acid
        sequence and then replicates the positions of each residue in
        the grid. The residues from the original protein are removed
        from the grid in the new lattice, and their coordinates are
        reassigned to match the original positions.
        Finally, the residues are placed in the appropriate locations
        on the new lattice grid.

        Returns
        -------
        Lattice
            A new Lattice object that is a copy of the current one
            with the protein residues positioned in the same way as
            in the original lattice.
        """
        new_protein = Protein(self.protein.aa_sequence)
        new_lattice = Lattice(new_protein)
        for i in range(self.protein.length):
            residue = self.protein.hp_sequence[i]
            new_lattice.move_residue(
                residue.number, (residue.i_coord, residue.j_coord)
            )
        return new_lattice

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
                    # not necessary to continue if the neighbour is
                    # not hydrophobic
                    if (
                        residue2.type == "H"
                        and not residue1.is_connected(residue2)
                        and residue1.is_adjacent(residue2)
                    ):
                        energy -= 1
        return energy

    def end_moves(self, residue):
        """
        Attempt to reposition a terminal residue.

        This method checks if the given residue is at one of the ends
        of the protein sequence. If it is the first or last residue,
        the function identifies the adjacent residue and retrieves
        its available positions on the grid. If there are available
        positions, it creates a new lattice configuration by moving
        the terminal residue to the first available position.

        Parameters
        ----------
        residue : Residue
            The residue object to be checked.

        Returns
        -------
        new_lattice : Lattice or None
            A new lattice configuration with the terminal residue
            moved to the first available position, or None if the
            residue is not terminal or no available positions exist.

        """
        neighbour = None
        if residue.number == 1:
            neighbour = self.protein.get_residue(2)
        elif residue.number == self.protein.length:
            neighbour = self.protein.get_residue(self.protein.length - 1)
        else:
            return None
        available_positions = self.get_adjacents_available_positions(neighbour)
        if len(available_positions) > 0:
            new_lattice = self.copy()
            new_lattice.move_residue(
                residue.number, available_positions[0]
            )  # first position by default
            return new_lattice
        return None

    def corner_moves(self, residue):
        """
        Attempt to reposition a non-terminal residue to a corner move.

        This method checks if the given residue is not a terminal one
        (i.e., neither the first nor the last in the protein sequence).
        If the residue is connected to two neighbors, the method
        identifies common available positions adjacent to both neighbors.
        If a corner move is found, the function returns a new lattice with
        the residue repositioned to the first available corner move.

        Parameters
        ----------
        residue : Residue
            The residue to check.

        Returns
        -------
        new_lattice : Lattice or None
            A new lattice with the residue moved to the first available
            corner position, or None if no corner move is possible or
            the residue is terminal.
        """
        if residue.number not in (1, self.protein.length):
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
            if len(cornermove) > 0:
                new_lattice = self.copy()
                new_lattice.move_residue(residue.number, cornermove[0])
                return new_lattice
        return None

    def crankshaft_first_corner(self, residue):
        """
        Reposition a residue in the first corner of a U-shaped bend.

        This method attempts to move both the residue and its adjacent neighbor
        to their symmetrical positions on the grid. The move is only valid if
        both new positions are within grid dimensions and unoccupied.

        Parameters
        ----------
        residue : Residue
            The residue object to be moved.

        Returns
        -------
        Lattice or None
            A new lattice configuration if the move is successful; otherwise,
            None if the move cannot be made.
        """
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
                and self.grid[sym_pos_i_plus_1[0], sym_pos_i_plus_1[1]] is None
            ):
                new_lattice = self.copy()
                new_lattice.move_residue(residue.number, sym_pos_i)
                new_lattice.move_residue(res_i_plus_1.number, sym_pos_i_plus_1)
                return new_lattice
        return None

    def crankshaft_second_corner(self, residue):
        """
        Reposition a residue in the second corner of a U-shaped bend.

        This method attempts to move both the residue and its adjacent neighbor
        to their symmetrical positions on the grid at the second corner of the
        bend. The move is only valid if both new positions are within grid
        dimensions and unoccupied.

        Parameters
        ----------
        residue : Residue
            The residue object to be moved.

        Returns
        -------
        Lattice or None
            A new lattice configuration if the move is successful; otherwise,
            None if the move cannot be made.
        """
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
                new_lattice = self.copy()
                new_lattice.move_residue(residue.number, sym_pos_i)
                new_lattice.move_residue(
                    res_i_minus_1.number, sym_pos_i_minus_1
                )
                return new_lattice
        return None

    def crankshaft_moves(self, residue):
        """
        Attempt a crankshaft move for a residue in a U-shaped bend.

        This method first determines if the residue is part of a U-shaped bend.
        If so, it checks whether the residue is located at the first or second
        corner of the bend. It then attempts to perform the appropriate
        crankshaft move (either the first or second corner move) and returns
        a new lattice configuration with the movement applied if successful.
        If the move is not possible, it returns None.

        Parameters
        ----------
        residue : Residue
            The residue object to be moved.

        Returns
        -------
        Lattice or None
            A new lattice configuration if a valid crankshaft move is made;
            otherwise, None if the move cannot be performed.
        """
        is_ushaped = self.protein.is_residue_in_ushaped_bend(residue)
        if not is_ushaped[0]:
            return None
        if is_ushaped[1] == 1:
            return self.crankshaft_first_corner(residue)
        return self.crankshaft_second_corner(residue)

    def successive_pulls(self, new_lattice, number_residue):
        """
        Perform a series of successive pull moves on the lattice.

        This method applies successive pull moves to adjust the
        positions of residues in the lattice to ensure that the
        sequence of residues becomes valid. It starts from the
        residue with a number number_residue - 2 and attempts
        to reposition each residue in sequence by moving it to
        the position of the residue two places ahead in the original
        lattice. This process continues until the sequence is valid
        or all residues from the starting point down to residue 1
        have been adjusted.

        Parameters
        ----------
        new_lattice : Lattice
            The lattice object where residues will be moved.
            This lattice is modified during the process to
            reflect the new positions of the residues.
        number_residue : int
            The starting residue number for the successive
            pulls. The function will start moving residues
            from number_residue - 2 down to residue 1.
        """
        current_number_residu = number_residue - 2
        while (
            current_number_residu >= 1
            and not new_lattice.protein.is_sequence_valid()
        ):
            new_coordinate_i = self.protein.get_residue(
                current_number_residu + 2
            ).i_coord
            new_coordinate_j = self.protein.get_residue(
                current_number_residu + 2
            ).j_coord
            new_lattice.move_residue(
                current_number_residu,
                (new_coordinate_i, new_coordinate_j),
            )
            current_number_residu -= 1

    @staticmethod
    def is_adjacent_position(pos1, pos2):
        """
        Check if two positions are adjacent on the grid.

        This method determines if two positions pos1 and pos2 are adjacent
        on the grid, meaning they share either the same row or column and are
        only one unit apart.

        Parameters
        ----------
        pos1 : tuple of int
            The (i, j) coordinates of the first position on the grid.
        pos2 : tuple of int
            The (i, j) coordinates of the second position on the grid.

        Returns
        -------
        bool
            True if the positions are adjacent (either horizontally
            or vertically), False otherwise.
        """
        return (abs(pos1[0] - pos2[0]) == 1 and pos1[1] == pos2[1]) or (
            abs(pos1[1] - pos2[1]) == 1 and pos1[0] == pos2[0]
        )

    @staticmethod
    def is_diagonal_position(pos1, pos2):
        """
        Check if two positions are diagonally adjacent on the grid.

        This method determines if two positions pos1 and pos2 are diagonally
        adjacent on the grid, meaning they differ by one unit both horizontally
        and vertically.

        Parameters
        ----------
        pos1 : tuple of int
            The (i, j) coordinates of the first position on the grid.
        pos2 : tuple of int
            The (i, j) coordinates of the second position on the grid.

        Returns
        -------
        bool
            True if the positions are diagonally adjacent, False otherwise.
        """
        return abs(pos1[0] - pos2[0]) == 1 and abs(pos1[1] - pos2[1]) == 1

    def get_position_l(self, residue):
        """
        Retrieve the position L.

        This method finds a position L on the grid that is diagonally
        adjacent to the given residue and adjacent to its neighbor
        (residue i+1). It iterates over all available adjacent positions
        of the neighbor and checks if any of them are diagonally adjacent
        to the current position of the residue. The first valid position
        is returned.

        Parameters
        ----------
        residue : Residue
            The residue object for which to find a diagonally
            adjacent position L.

        Returns
        -------
        position_l : tuple of int or None
            A tuple (i, j) representing the coordinates of
            the position L, or None if no suitable diagonal
            position is found.
        """
        position_l = None
        residu_plus_1 = self.protein.get_residue(residue.number + 1)
        adj_res_plus_1 = self.get_adjacents_available_positions(residu_plus_1)
        for elem in adj_res_plus_1:
            if Lattice.is_diagonal_position(
                (residue.i_coord, residue.j_coord), elem
            ):
                position_l = elem
                break  # stop at the first solution
        return position_l

    def get_position_c(self, residue, position_l):
        """
        Retrieve the position C.

        This method finds a position C on the grid that is adjacent to the
        current position of the specified residue and the given position L.
        It iterates over all available adjacent positions of the residue
        and checks if any of them are also adjacent to L. Once such a
        position is found, it is returned.

        Parameters
        ----------
        residue : Residue
            The residue object for which to find an adjacent position C.
        position_l : tuple of int
            The coordinates (i, j) of the position L that must be
            adjacent to C.

        Returns
        -------
        position_c : tuple of int or None
            A tuple (i, j) representing the coordinates of the position C,
            or None if no suitable position is found.
        """
        position_c = None
        adj_res = self.get_adjacents_available_positions(residue)
        for elem in adj_res:
            if Lattice.is_adjacent_position(elem, position_l):
                position_c = elem
                break
        return position_c

    def pull_moves(self, residue):
        """
        Attempt to perform a pull move on the given residue.

        This method attempts to perform a pull move by shifting the residue i
        to an adjacent position L and its preceding residue i-1
        (if applicable) to a new position C adjacent to both i and L.

        The pull move involves the following steps:
        1. Retrieve position L, which is adjacent to i+1 and diagonally
        adjacent to i. The check ensures that residue i+1 exists in the
        protein sequence.
        2. If L is available, retrieve position C, which is adjacent
        to both i and L. The case where C is occupied by i-1
        (the preceding residue) is not implemented, as it falls under
        a corner move case.
        3. If both C and L are accessible, perform the pull move by moving
        residue i to L. If residue i is the first one in the sequence, return
        the new lattice.
        4. If C is not occupied by i-1, move i to L and i-1 to C.
        5. If the sequence is invalid after the move, attempt successive pull
        moves on preceding residues until a valid sequence is obtained or the
        first residue is reached.

        Parameters
        ----------
        residue : Residue
            The residue object to be repositioned.

        Returns
        -------
        new_lattice : Lattice or None
            A new lattice configuration if a valid pull move is performed,
            or None if the move cannot be made.
        """
        if residue.number + 1 <= self.protein.length:
            position_l = self.get_position_l(residue)
            position_c = None
            if position_l is not None:
                position_c = self.get_position_c(residue, position_l)
            if position_c is not None:
                new_lattice = self.copy()
                new_lattice.move_residue(residue.number, position_l)
                if residue.number == 1:
                    return new_lattice
                new_lattice.move_residue(residue.number - 1, position_c)
                self.successive_pulls(new_lattice, residue.number)
                return new_lattice
        return None
