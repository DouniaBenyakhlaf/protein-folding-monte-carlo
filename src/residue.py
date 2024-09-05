class Residue:
    """
    A class used to model a residue.

    Attributes
    ----------
    name : str
        The code name of the residue (one letter).
    number : int
        The residue position number.
    type : str
        The residue type: hydrophobic (H) or polar (P).
    i_coord : int
        The i coordinate of the residue in the lattice.
    j_coord : int
        The j coordinate of the residue in the lattice.
    """

    def __init__(self, type, number, i_coord, j_coord):
        """
        Initializes a Residue instance with a type, number, and coordinates.

        Parameters
        ----------
        type : str
            The residue type : hydrophobic (H) or polar (P).
        number : int
            The residue number.
        i_coord : int
            The i coordinate of the residue in the lattice.
        j_coord : int
            The j coordinate of the residue in the lattice.
        """
        self.type = type
        self.number = number
        self.i_coord = i_coord
        self.j_coord = j_coord

    def __str__(self):
        """
        Return a string representation of the residue.

        Returns
        -------
        str
            A formatted string describing the residue, including its type, number,
            and lattice position.
        """
        return f"The residue {self.type}{self.number} is in position ({self.i_coord},{self.j_coord}) in the lattice."

    def is_connected(self, res):
        """
        Determines if the current residue is connected to another residue
        based on their position numbers.

        Parameters
        ----------
        res : Residue
            Another residue to check for connection.

        Returns
        -------
        bool
            True if the residues are connected (i.e., their position numbers
            differ by 1), False otherwise.
        """
        return abs(self.number - res.number) == 1

    def is_adjacent(self, residue2):
        """
        Determine if the current residue is adjacent to another residue.

        Two residues are considered adjacent if they are next to each other
        either horizontally or vertically in the lattice.

        Parameters
        ----------
        residue2 : Residue
            The other residue to compare with the current residue.

        Returns
        -------
        bool
            True if the current residue is adjacent to residue2, False otherwise.
        """
        return (
            abs(self.i_coord - residue2.i_coord) == 1
            and self.j_coord == residue2.j_coord
        ) or (
            abs(self.j_coord - residue2.j_coord) == 1
            and self.i_coord == residue2.i_coord
        )

    def get_symetrical_position(self, residue):
        """
        Calculate the symmetrical position of the current residue with respect to another residue.

        This method computes the symmetrical position of the current residue (self) relative to
        another residue. The symmetry is calculated in a 2D grid, based on the coordinates of the
        given residue.

        Parameters
        ----------
        residue : Residue
            The residue with respect to which the symmetrical position is calculated.

        Returns
        -------
        tuple of int
            A tuple (x, y) representing the symmetrical position of the current residue in the grid.
        """
        return (
            2 * residue.i_coord - self.i_coord,
            2 * residue.j_coord - self.j_coord,
        )
