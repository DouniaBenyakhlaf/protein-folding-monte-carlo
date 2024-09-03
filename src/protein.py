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

    # Dictionary containing the type of each residue
    RESIDUE_TYPE_DICT = {
        "A": "H",  # Alanine
        "R": "P",  # Arginine
        "N": "P",  # Asparagine
        "D": "P",  # Aspartic acid
        "C": "H",  # Cysteine
        "E": "P",  # Glutamic acid
        "Q": "P",  # Glutamine
        "G": "H",  # Glycine
        "H": "P",  # Histidine
        "I": "H",  # Isoleucine
        "L": "H",  # Leucine
        "K": "P",  # Lysine
        "M": "H",  # Methionine
        "F": "H",  # Phenylalanine
        "P": "H",  # Proline
        "S": "P",  # Serine
        "T": "P",  # Threonine
        "W": "H",  # Tryptophan
        "Y": "P",  # Tyrosine
        "V": "H",  # Valine
    }

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
            str_residue = str_sequence[position]
            type_residue = self.RESIDUE_TYPE_DICT[str_residue]
            new_residue = Residue(
                str_residue, position + 1, type_residue, -1, -1
            )
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
