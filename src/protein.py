from residue import *


class Protein:
    """
    A class used to model a protein.

    Attributes
    ----------
    sequence : list[Residue]
        A list of Residue objects representing the protein sequence.
    length : int
        The number of residues in the protein sequence.
    first_residue : Residue
        The first residue in the protein sequence.
    last_residue : Residue
        The last residue in the protein sequence.
    """

    def __init__(self, str_sequence):
        """
        Initialize a Protein instance with a string sequence.

        Parameters
        ----------
        str_sequence : str
            A string representing the sequence of the protein.
        """
        self.sequence = []
        protein_length = len(str_sequence)
        for position in range(protein_length):
            type_residue = str_sequence[position]
            new_residue = Residue(type_residue, position + 1, -1, -1)
            self.sequence.append(new_residue)
        self.length = protein_length
        self.first_residue = self.sequence[0]
        self.last_residue = self.sequence[-1]

    def print_sequence(self):
        """
        Prints a description the protein and each residue of the protein sequence.
        """
        print(f"Protein length = {self.length}")
        for residue in self.sequence:
            residue.res_description()
