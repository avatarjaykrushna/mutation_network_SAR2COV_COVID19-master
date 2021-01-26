from utils import read_file, translation, check_for_neutral_mutations, possible_spike_protein_combinations, \
    spike_protein_mutations, mutation_network
import pprint


def main():
    genome = read_file("wuhan")
    pp = pprint.PrettyPrinter(indent=4)
    """Question 1"""
    # polypeptide = translation(genome)
    # neutral_mutations_count_per_codon = check_for_neutral_mutations(genome)
    # print(neutral_mutations_count_per_codon)

    """Question 3"""
    # print("** Question 3")
    # possible_spike_protein_combinations(amino_acids="YLNDT")

    """Question 4"""
    """    
    # non_neutral_one_mutation_away = spike_protein_mutations(genome)
    # print(non_neutral_one_mutation_away)
    """

    sar_2 = "TATTTGAACGATACC"

    # neutral_mutations_count_per_codon = check_for_neutral_mutations(genome)
    mutation_list = spike_protein_mutations(sar_2)
    print(mutation_list)

    mutation_network(amino_acid_seq="", n_mutations=15)


if __name__ == "__main__":
    main()
