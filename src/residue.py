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
