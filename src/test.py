"""Test the replica exchange Monte Carlo algorithm for protein folding."""

import random
from protein import Protein
from remcsimulation import REMCSimulation


def verify_number():
    """
    Prompt the user for input and verify that it is a valid number.

    This function repeatedly prompts the user for input until a valid
    integer between 1 and 12 is provided. If the input is invalid, an
    appropriate error message is displayed.

    Returns
    -------
    int
        A valid number between 1 and 12.
    """
    while True:
        choice = input()
        try:
            number_choice = int(choice)
            if number_choice > 0 and number_choice <= 12:
                return number_choice
            else:
                print("Error: Please choose a valid number.")
        except ValueError:
            print("Error: Please choose a valid number.")


def verify_energy():
    """
    Prompt the user for input and verify that it is a valid integer.

    This function repeatedly prompts the user for input until a valid
    negative integer is provided. If the input is invalid, an appropriate
    error message is displayed.

    Returns
    -------
    int
        A valid negative integer.
    """
    while True:
        choice = input()
        try:
            number_choice = int(choice)
            if number_choice < 0:
                return number_choice
            else:
                print(
                    "Error: Please choose a valid energy value (must be negative)."
                )
        except ValueError:
            print("Error: Please choose a valid energy value.")


def verify_sequence(sequence):
    """
    Validate the protein sequence.

    Checks if each character in the provided sequence is a valid residue according
    to the residue type dictionary.

    Parameters
    ----------
        sequence : str
            The protein sequence to be validated.

    Returns
    -------
    tuple
        A tuple containing:
            - bool: True if the sequence is valid, False otherwise.
            - str: A message indicating the result of the validation.
    """
    for char in sequence:
        if char not in Protein.RESIDUE_TYPE_DICT:
            return False, f"Erreur : The residue '{char}' is not valid."
    return True, "The sequence is valid."


if __name__ == "__main__":

    random.seed(1)

    PROTEIN_1 = "GRGRRGGRGRRGRGGRRGRG"
    PROTEIN_2 = "GGRRGRRGRRGRRGRRGRRGRRGG"
    PROTEIN_3 = "RRGRRGGRRRRGGRRRRGGRRRRGG"
    PROTEIN_4 = "RRRGGRRGGRRRRRGGGGGGGRRGGRRRRGGRRGRR"
    PROTEIN_5 = "RRGRRGGRRGGRRRRRGGGGGGGGGGRRRRRRGGRRGGRRGRRGGGGG"
    PROTEIN_6 = "GGRGRGRGRGGGGRGRRRGRRRGRRRRGRRRGRRRGRGGGGRGRGRGRGG"
    PROTEIN_7 = "RRGGGRGGGGGGGGRRRGGGGGGGGGGRGRRRGGGGGGGGGGGGRRRRGGGGGGRGGRGR"
    PROTEIN_8 = (
        "GGGGGGGGGGGGRGRGRRGGRRGGRRGRRGGRRGGRRGRRGGRRGGRRGRGRGGGGGGGGGGGG"
    )
    PROTEIN_9 = "GGGGRRRRGGGGGGGGGGGGRRRRRRGGGGGGGGGGGGRRRGGGGGGGGGGGGRRRGGGGGGGGGGGGRRRGRRGGRRGGRRGRG"
    PROTEIN_10 = "RRRGGRRGGGGRRGGGRGGRGGRGGGGRRRRRRRRGGGGGGRRGGGGGGRRRRRRRRRGRGGRGGGGGGGGGGGRRGGGRGGRGRRGRGGGRRRRRRGGG"
    PROTEIN_11 = "RRRRRRGRGGRRRRRGGGRGGGGGRGGRRRRGGRRGGRGGGGGRGGGGGGGGGGRGGRGGGGGGGRRRRRRRRRRRGGGGGGGRRGRGGGRRRRRRGRGG"

    LIST_PROTEIN = [
        PROTEIN_1,
        PROTEIN_2,
        PROTEIN_3,
        PROTEIN_4,
        PROTEIN_5,
        PROTEIN_6,
        PROTEIN_7,
        PROTEIN_8,
        PROTEIN_9,
        PROTEIN_10,
        PROTEIN_11,
    ]
    OPTIM_ENERGIES = [-9, -9, -8, -14, -23, -21, -36, -42, -53, -50, -48]

    # Code to color text
    ORANGE = "\033[38;5;208m"
    PURPLE = "\033[95m"
    BLUE = "\033[34m"
    RED = "\033[31m"
    RESET = "\033[0m"

    print(
        f"{ORANGE}Welcome to the protein folding program using the Replica Exchange Monte Carlo algorithm.\n{RESET}"
    )

    print(
        f"{PURPLE}To test the algorithm on the following examples, choose a number:{RESET}"
    )
    for i in range(11):
        print(
            f"{BLUE}[{i+1}]{RESET} Protein {i+1} |  length = {len(LIST_PROTEIN[i])} | optimal_energy = {OPTIM_ENERGIES[i]}"
        )
    print(f"{BLUE}[12]{RESET} I want to test it on my own protein")
    print(
        f"\n{ORANGE}[INFO]{RESET} For the first 3 proteins that converge quickly, a PyMOL script is generated for some states of each replica in \
the results directory. The final state of each replica is also represented in a Networkx graph.\nFor other examples that do not converge after \
50,000 iterations, only the final state of each replica\nis represented in a PyMOL script and a Networkx graph, available in the results directory.\n"
    )
    NUMBER_CHOICE = verify_number()
    if NUMBER_CHOICE == 12:
        print(
            "Enter your protein sequence (Residues must be in one-letter code and uppercase.): "
        )
        SEQUENCE = input().upper()
        VALID, MESSAGE = verify_sequence(SEQUENCE)
        print(MESSAGE)
        print("Enter the optimal energy of your protein : ")
        ENERGY = verify_energy()
        REMC = REMCSimulation(SEQUENCE, optimal_energy=ENERGY, max_iter=50000)
        REMC.run()
        NUM_REPLICA = 0
        for replica in REMC.replicas.values():
            print(f"Energy replica {NUM_REPLICA} = {replica.compute_energy()}")
            FILE_NAME = f"replica_{NUM_REPLICA}"
            replica.protein.plot_protein_graph(f"../results/{FILE_NAME}")
            replica.protein.generate_pymol_script(
                FILE_NAME,
                1,
                filename=f"../results/{FILE_NAME}.pml",
            )
            NUM_REPLICA += 1
    else:
        print(f"You have chosen the protein {NUMBER_CHOICE}.")
        print(f"{RED}Start of the protein folding :{RESET}")
        REMC = REMCSimulation(
            LIST_PROTEIN[NUMBER_CHOICE - 1],
            optimal_energy=OPTIM_ENERGIES[NUMBER_CHOICE - 1],
            max_iter=50000,
        )
        if NUMBER_CHOICE <= 3:
            REMC.run(display=True)
            NUM_REPLICA = 0
            for replica in REMC.replicas.values():
                replica.protein.plot_protein_graph(
                    f"../results/replica_{NUM_REPLICA}"
                )
                print(
                    f"Energy replica {NUM_REPLICA} = {replica.compute_energy()}"
                )
                NUM_REPLICA += 1
        else:
            REMC.run()
            NUM_REPLICA = 0
            for replica in REMC.replicas.values():
                print(
                    f"Energy replica {NUM_REPLICA} = {replica.compute_energy()}"
                )
                FILE_NAME = f"replica_{NUM_REPLICA}"
                replica.protein.plot_protein_graph(f"../results/{FILE_NAME}")
                replica.protein.generate_pymol_script(
                    FILE_NAME,
                    1,
                    filename=f"../results/{FILE_NAME}.pml",
                )
                NUM_REPLICA += 1
        print(
            f"{PURPLE}Don’t forget to check the graphs/scripts generated in the results directory.{RESET}"
        )
