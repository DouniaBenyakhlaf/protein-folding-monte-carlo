"""Module for handling protein."""

import networkx as nx
import matplotlib.pyplot as plt
from residue import Residue


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

    # Dictionary containing the type of each residue :
    # H for hydrophobic and P for polar
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
        description = f"Protein length = {self.length}\n\
AA_Sequence : {self.aa_sequence}\nHP_Sequence : "
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

    def is_residue_in_ushaped_bend(self, residue):
        """
        Check if a residue is part of a U-shaped bend in the chain.

        This method determines whether the given residue is part of a
        U-shaped bend by checking the adjacency of other residues in
        the chain. The method verifies if the residue is in either the
        first or the second corner of a U-shaped bend based on the
        adjacency conditions.

        Parameters
        ----------
        residue : Residue
            The residue to be checked for being part of a U-shaped bend.

        Returns
        -------
        tuple of (bool, int)
            A tuple where the first element is a boolean indicating if
            the residue is part of a U-shaped bend, and the second element
            is an integer indicating the corner type:
            - 1 for the first corner
            - 2 for the second corner
            - 0 if the residue is not part of a U-shaped bend.
        """
        # residue i is part of a u-shaped bend in the chain if residue
        # i-1 and i+2 are adjacents (first corner)
        if (residue.number - 1) >= 1 and (residue.number + 2) < self.length:
            res_minus_1 = self.get_residue(residue.number - 1)
            res_plus_2 = self.get_residue(residue.number + 2)
            if res_minus_1.is_adjacent(res_plus_2):
                return (True, 1)
        # residue i is part of a u-shaped bend in the chain if residue
        # i-2 and i+1 are adjacents (second corner)
        if (residue.number - 2) >= 1 and (residue.number + 1) < self.length:
            res_minus_2 = self.get_residue(residue.number - 2)
            res_plus_1 = self.get_residue(residue.number + 1)
            if res_minus_2.is_adjacent(res_plus_1):
                return (True, 2)
        return (False, 0)

    def is_sequence_valid(self):
        """
        Check if the HP sequence is valid.

        This method verifies if each consecutive pair of residues in the
        HP sequence is adjacent. The sequence is considered valid if every
        residue is adjacent to the next one in the sequence.

        Returns
        -------
        bool
            True if all consecutive residues in the sequence are adjacent,
            otherwise False.
        """
        for i in range(self.length - 1):
            if not self.hp_sequence[i].is_adjacent(self.hp_sequence[i + 1]):
                return False
        return True

    def plot_protein_graph(self, filename):
        """
        Create and save a 2D graph of the protein residues to a file.

        This method uses the NetworkX and Matplotlib libraries to create a
        2D graph where each node represents a residue in the protein sequence.
        Connections between nodes indicate adjacent residues. The graph is
        saved to a file specified by the 'filename' parameter.

        Parameters
        ----------
        filename : str
            The path to the file where the graph will be saved.
        """
        graph = nx.Graph()

        # Add each residue as a node to the graph with its 2D coordinates.
        for residue in self.hp_sequence:
            graph.add_node(
                residue.number,
                pos=(residue.j_coord, -residue.i_coord),
                type=residue.type,
            )

        # Add edges between connected residues.
        for i in range(len(self.hp_sequence) - 1):
            if self.hp_sequence[i].is_connected(self.hp_sequence[i + 1]):
                graph.add_edge(
                    self.hp_sequence[i].number, self.hp_sequence[i + 1].number
                )

        # Get the positions of the nodes.
        pos = nx.get_node_attributes(graph, "pos")

        # Create node color mapping based on residue type.
        color_map = {
            residue.number: "red" if residue.type == "H" else "lightblue"
            for residue in self.hp_sequence
        }

        # Create labels for the nodes with the residue number and type.
        labels = {
            residue.number: f"{residue.number} ({residue.type})"
            for residue in self.hp_sequence
        }

        # Plot the graph.
        plt.figure(figsize=(10, 10))
        nx.draw(
            graph,
            pos,
            with_labels=True,
            labels=labels,
            node_color=[color_map[node] for node in graph.nodes],
            node_size=2800,
            font_weight="bold",
            font_size=14,
        )
        plt.title("2D Representation of the Protein Sequence")
        # Save the figure to the specified file.
        plt.savefig(filename)
        plt.close()  # Close the plot to free up memory

    def generate_pymol_script(
        self, name, state, filename="protein_representation.pml"
    ):
        """
        Generate a PyMOL script (.pml) for representing a protein.

        Uses pseudo-atoms and bonds for connecting adjacent residues.

        Parameters
        ----------
        protein : Protein
            An instance of the Protein class with residue coordinates.
        filename : str, optional
            The name of the file to save the PyMOL script to (default
            is "protein_representation.pml").
        """
        right = "w"
        if state > 1:
            right = "a"
        with open(filename, right, encoding="utf-8") as f:
            # Start a new object for the protein
            # f.write(f"frame {state}\n")
            # Write pseudo-atoms for each residue
            for residue in self.hp_sequence:
                res_name = f"residue_{residue.number}_{state}"
                x, y = (
                    residue.j_coord,
                    residue.i_coord,
                )  # PyMOL uses x, y, z coordinates
                f.write(
                    f'pseudoatom {name}, pos=[{x}, {y}, 0], name="{res_name}",\
resn="{residue.type}", chain="A", b={residue.number}, state={state}\n'
                )

            f.write(f"frame {state}\n")
            # Write bonds between consecutive residues
            for i in range(self.length - 1):
                res1 = f"residue_{self.hp_sequence[i].number}_{state}"
                res2 = f"residue_{self.hp_sequence[i + 1].number}_{state}"
                f.write(
                    f"bond ({name} and name {res1}), ({name} \
and name {res2})\n"
                )

            # Write commands to visualize the structure
            f.write(f"show sticks, {name}\n")
            f.write(f"zoom {name}\n")

            # Optional: Color residues based on type
            f.write("color red, (resn H)\n")
            f.write("color blue, (resn P)\n")
