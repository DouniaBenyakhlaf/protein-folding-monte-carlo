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
    x_coord : int
        The x coordinate of the residue in the lattice.
    y_coord : int
        The y coordinate of the residue in the lattice.
    """

    def __init__(self, name, number, type, x_coord, y_coord):
        """
        Initializes a Residue instance with a name, number, type, and coordinates.

        Parameters
        ----------
        name : str
            The code name of the residue (one letter).
        number : int
            The residue number.
        type : str
            The residue type : hydrophobic (H) or polar (P).
        x_coord : int
            The x coordinate of the residue in the lattice.
        y_coord : int
            The y coordinate of the residue in the lattice.
        """
        self.name = name
        self.number = number
        self.type = type
        self.x_coord = x_coord
        self.y_coord = y_coord

    def res_description(self):
        """
        Prints a description of the residue, including its name, number, type,
        and position in the lattice.
        """
        print(
            f"The residue {self.name}{self.number} of type {self.type} is in position ({self.x_coord},{self.y_coord}) in the lattice."
        )

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
