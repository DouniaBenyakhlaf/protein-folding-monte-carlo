from residue import *


class Protein:
    """
    A class used to model a protein.

    Attributes
    ----------
    aa_sequence : str
        A string of the amino acid sequence (one letter code)
    hp_sequence : list[Residue]
        A list of Residue objects representing the protein sequence.
    length : int
        The number of residues in the protein sequence.
    first_residue : Residue
        The first residue in the protein sequence.
    last_residue : Residue
        The last residue in the protein sequence.
    """

    # Dictionary containing the type of each residue : H for hydrophobic and P for polar
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
            A string representing the amino acid sequence of the protein.
        """
        self.aa_sequence = str_sequence
        self.hp_sequence = []
        protein_length = len(str_sequence)
        for position in range(protein_length):
            type_residue = Protein.RESIDUE_TYPE_DICT[str_sequence[position]]
            new_residue = Residue(type_residue, position + 1, -1, -1)
            self.hp_sequence.append(new_residue)
        self.length = protein_length
        self.first_residue = self.hp_sequence[0]
        self.last_residue = self.hp_sequence[-1]

    def __str__(self):
        """
        Return a string representation of the protein.

        Returns
        -------
        str
            A formatted string describing the protein, including its length and
            the sequence of residues with their types and numbers.
        """
        description = f"Protein length = {self.length}\nAA_Sequence : {self.aa_sequence}\nHP_Sequence : "
        for residue in self.hp_sequence:
            description += f"{residue.type}{residue.number} "
        return description

    def get_residue(self, residue_number):
        """
        Retrieve a residue from the protein sequence by its number.

        Parameters
        ----------
        residue_number : int
            The position number of the residue in the sequence (1-based index).

        Returns
        -------
        Residue
            The Residue object at the specified position in the sequence.
        """
        return self.hp_sequence[residue_number - 1]

    def is_residue_in_ushaped_bend(self, residue, is_first_corner):
        """
        Check if a residue is part of a U-shaped bend in the chain.

        This method determines whether the given residue is part of a U-shaped bend by checking the
        adjacency of other residues. Depending on the value of `is_first_corner`, it either checks
        the adjacency between residue i-1 and i+2 or between residue i-2 and i+1.

        Parameters
        ----------
        residue : Residue
            The residue to be checked for being part of a U-shaped bend.
        is_first_corner : bool
            If True, checks if residue i-1 and i+2 are adjacent. If False, checks if residue i-2
            and i+1 are adjacent.

        Returns
        -------
        bool
            True if the residue is part of a U-shaped bend, otherwise False.
        """
        if is_first_corner:
            # residue i is part of a u-shaped bend in the chain if residue i-1 and i+2 are adjacents
            if (residue.number - 1) >= 1 and (
                residue.number + 2
            ) < self.length:
                return self.get_residue(residue.number - 1).is_adjacent(
                    self.get_residue(residue.number + 2)
                )
            else:
                return False
        else:
            # residue i is part of a u-shaped bend in the chain if residue i-2 and i+1 are adjacents
            if (residue.number - 2) >= 1 and (
                residue.number + 1
            ) < self.length:
                return self.get_residue(residue.number - 2).is_adjacent(
                    self.get_residue(residue.number + 1)
                )
            else:
                return False

    def is_sequence_valid(self):
        """
        Check if the HP sequence is valid based on adjacency of consecutive residues.

        This method verifies if each consecutive pair of residues in the HP sequence is adjacent.
        The sequence is considered valid if every residue is adjacent to the next one in the sequence.

        Returns
        -------
        bool
            True if all consecutive residues in the sequence are adjacent, otherwise False.
        """
        for i in range(self.length - 1):
            if not (self.hp_sequence[i].is_adjacent(self.hp_sequence[i + 1])):
                return False
        return True
