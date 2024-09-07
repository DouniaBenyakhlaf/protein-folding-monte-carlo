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

        This method creates a deep copy of the protein and lattice. It initializes a new protein
        object with the same amino acid sequence and then replicates the positions of each residue
        in the grid. The residues from the original protein are removed from the grid in the new lattice,
        and their coordinates are reassigned to match the original positions. Finally, the residues are
        placed in the appropriate locations on the new lattice grid.

        Returns
        -------
        Lattice
            A new Lattice object that is a copy of the current one with the protein residues
            positioned in the same way as in the original lattice.
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
                    # not necessary to continue if the neighbour is not hydrophobic
                    if (
                        residue2.type == "H"
                        and not residue1.is_connected(residue2)
                        and residue1.is_adjacent(residue2)
                    ):
                        energy -= 1
        return energy

    def end_moves(self, residue):
        """
        Attempt to reposition a terminal residue based on its neighbor's available positions.

        This method checks if the given residue is at one of the ends of the protein sequence.
        If it is the first or last residue, the function identifies the adjacent residue and
        retrieves its available positions on the grid. If there are available positions, it
        creates a new lattice configuration by moving the terminal residue to the first available
        position.

        Parameters
        ----------
        residue : Residue
            The residue object to be checked.

        Returns
        -------
        new_lattice : Lattice or None
            A new lattice configuration with the terminal residue moved to the first available
            position, or None if the residue is not terminal or no available positions exist.

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
        Attempt to reposition a non-terminal residue to a corner move, if available.

        This method checks if the given residue is not a terminal one (i.e., neither the
        first nor the last in the protein sequence). If the residue is connected to two
        neighbors, the method identifies common available positions adjacent to both
        neighbors. If a corner move is found, the function returns a new lattice with
        the residue repositioned to the first available corner move.

        Parameters
        ----------
        residue : Residue
            The residue to check.

        Returns
        -------
        new_lattice : Lattice or None
            A new lattice with the residue moved to the first available corner position,
            or None if no corner move is possible or the residue is terminal.
        """
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
            if len(cornermove) > 0:
                new_lattice = self.copy()
                new_lattice.move_residue(residue.number, cornermove[0])
                return new_lattice
        return None

    def crankshaft_first_corner(self, residue):
        """
        Attempt to reposition a residue (in the first corner
        of a U-shaped bend) and its neighbor.

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
        Attempt to reposition a residue (in the second corner of a U-shaped bend)
        and its neighbor.

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
        Attempt a crankshaft move for a residue based on its position in a
        U-shaped bend.

        This method first determines if the residue is part of a U-shaped bend.
        If so, it checks whether the residue is located at the first or second
        corner of the bend. It then attempts to perform the appropriate crankshaft
        move (either the first or second corner move) and returns a new lattice
        configuration with the movement applied if successful. If the move is not
        possible, it returns None.

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
        elif is_ushaped[1] == 1:
            return self.crankshaft_first_corner(residue)
        else:
            return self.crankshaft_second_corner(residue)

    # Essayer de simplifier cette fonction
    def possible_pullmoves(self, residue):
        # 1) recuperer L qui est adjacent a i+1 et diagonalement adjacent a i (verifier que i-1 existe)
        if residue.number + 1 <= self.protein.length:
            print("pas dernier residue")
            residu_plus_1 = self.protein.get_residue(residue.number + 1)
            adj_res_plus_1 = self.get_adjacents_available_positions(
                residu_plus_1
            )
            print(adj_res_plus_1)
            position_l = None
            position_c = None
            for elem in adj_res_plus_1:
                if (
                    abs(residue.i_coord - elem[0]) == 1
                    and abs(residue.j_coord - elem[1]) == 1
                ):
                    position_l = elem
                    break
            if position_l is not None:
                print("il existe une position L")
                # 2) recuperer C qui est adjacent à i et L
                # je n'implémente pas le cas ou C est occupe par i-1 car revient au corner moves
                # si je change d'avis il faut juste verifier si le residu i-1 est adjacent a L
                adj_res = self.get_adjacents_available_positions(residue)
                for elem in adj_res:
                    if (
                        abs(elem[0] - position_l[0]) == 1
                        and elem[1] == position_l[1]
                    ) or (
                        abs(elem[1] - position_l[1]) == 1
                        and elem[0] == position_l[0]
                    ):
                        position_c = elem
                        break
            # On peut faire le pull moves
            if position_c is not None:
                print("Il existe une position C")
                if residue.number == 1:
                    print("pas de i-1")
                    print(
                        f"New position of res {residue.number} = {position_l}\n"
                    )
                else:
                    # 3) si C et L sont accessible on peut faire le pull moves
                    # 4) si oui alors si C est occupe par le residu i-1, il faut juste bouger i dans L (= corner move)
                    # 5) si C pas occupe pas i-1, il faut bouger i dans L et i-1 dans C
                    # 6) si sequence pas valide alors tant que pas valide ou on atteint le premier residu on bouge res tmp a la position du res+2. On commence par res-2.
                    print(
                        f"New position of res {residue.number} = {position_l}\nNew position of res {residue.number - 1} = {position_c}"
                    )
                    current_number_residu = residue.number - 2
                    while (
                        current_number_residu >= 1
                        and not self.protein.is_sequence_valid()
                    ):
                        # ne rentre jamais ici car il faut appliquer les modifications cette fois-ci... Les étapes precedentes ont l'air de fonctionner
                        print("sequence pas correcte")
                        new_coordinate_i = self.protein.get_residue(
                            current_number_residu + 2
                        ).i_coord
                        new_coordinate_j = self.protein.get_residue(
                            current_number_residu + 2
                        ).j_coord
                        print(
                            f"New position of res {residue.number - current_number_residu} = {(new_coordinate_i, new_coordinate_j)}\n"
                        )
                        current_number_residu -= 1
