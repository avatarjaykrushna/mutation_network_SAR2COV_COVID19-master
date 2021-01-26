import itertools
import matplotlib.pylab as plt
import collections

available_genome_sequences = ["wuhan", "philippines"]
codon_wheel = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": ">",
    "TAG": ">",
    "TGT": "C",
    "TGC": "C",
    "TGA": ">",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "<",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}
reverse_codon_wheel = {
    "F": ["ATC", "ATT"],
    "L": ["ATA", "ATG", "CTT", "CTC", "CTA", "CTG"],
    "S": ["ACT", "ACA", "ACG", "ACC", "AGC", "AGT"],
    "Y": ["TAT", "TAC"],
    ">": ["TAA", "TAG", "TGA"],
    "C": ["TGT", "TGC"],
    "W": ["TGG"],
    "P": ["CCT", "CCA", "CCC", "CCG"],
    "H": ["CAT", "CAC"],
    "Q": ["CAA", "CAG"],
    "R": ["CGT", "CGC", "CGA", "CGG", "AGG", "AGA"],
    "I": ["ATA", "ATC", "ATT"],
    "<": ["ATG"],
    "T": ["ACG", "ACA", "ACC", "ACT"],
    "N": ["AAC", "AAT"],
    "K": ["AAG", "AGA"],
    "V": ["GTA", "GTG", "GTC", "GTT"],
    "A": ["GCA", "GCG", "GCC", "GCT"],
    "D": ["GAC", "GAT"],
    "E": ["GAG", "GAA"],
    "G": ["GGA", "GGG", "GGC", "GGT"]
}


def read_file(file_id):
    if file_id in available_genome_sequences:
        gen_seq_file = open("data/sars_2_virus_" + file_id + "_genome_sequence.txt", "r")
    else:
        return None
    gen_seq_data = ""
    for line in gen_seq_file:
        gen_seq_data += line.replace(" ", "").rstrip()

    gen_seq_file.close()

    return gen_seq_data


def translation(genome):
    codon = ""
    polypeptide = ""

    for nucleotide in genome:
        if len(codon) < 3:
            codon += nucleotide.upper()
        else:
            polypeptide += codon_wheel.get(codon, "?")
            codon = nucleotide.upper()

    # ignore the remaining incomplete codon
    # polypeptide += codon

    return polypeptide


def get_mutation_list(codon):
    all_bases = ['A', 'C', 'G', 'T']
    mutation_list = []

    first_position_complement_list = [x for x in all_bases if x != codon[0]]
    second_position_complement_list = [x for x in all_bases if x != codon[1]]
    third_position_complement_list = [x for x in all_bases if x != codon[2]]

    mutation_list.append(codon[:2] + third_position_complement_list[0])
    mutation_list.append(codon[:2] + third_position_complement_list[1])
    mutation_list.append(codon[:2] + third_position_complement_list[2])

    mutation_list.append(codon[0] + second_position_complement_list[0] + codon[2])
    mutation_list.append(codon[0] + second_position_complement_list[1] + codon[2])
    mutation_list.append(codon[0] + second_position_complement_list[2] + codon[2])

    mutation_list.append(first_position_complement_list[0] + codon[1:])
    mutation_list.append(first_position_complement_list[1] + codon[1:])
    mutation_list.append(first_position_complement_list[2] + codon[1:])

    return mutation_list


def calculate_neutral_mutations(codon):
    neutral_count = 0
    single_mutations = get_mutation_list(codon)

    for mutation in single_mutations:
        if codon_wheel[mutation] == codon_wheel[codon]:
            neutral_count += 1

    return neutral_count


def calculate_non_neutral_mutations_with_mutated_codons(codon):
    neutral_count = 0
    mutated_codons = []
    single_mutations = get_mutation_list(codon)

    for mutation in single_mutations:
        if codon_wheel[mutation] != codon_wheel[codon]:
            neutral_count += 1
            mutated_codons.append(mutation)

    return neutral_count, mutated_codons


def calculate_non_neutral_mutations(codon):
    non_neutral_count = 0
    single_mutations = get_mutation_list(codon)

    for mutation in single_mutations:
        if codon_wheel[mutation] != codon_wheel[codon]:
            non_neutral_count += 1

    return non_neutral_count


def mutation_network(amino_acid_seq, n_mutations):
    if n_mutations is None:
        n_mutations = 0

    for codon in amino_acid_seq:
        pass


def check_for_neutral_mutations(genome):
    codon = ""
    amino_acid_sequence = {}

    for nucleotide in genome:
        if len(codon) < 3:
            codon += nucleotide.upper()
        else:
            neutral_mutations = calculate_neutral_mutations(codon)
            amino_acid_sequence[codon] = neutral_mutations
            codon = nucleotide.upper()

    x, y = zip(*sorted(amino_acid_sequence.items()))
    plt.plot(x, y)
    plt.xticks(x, x, rotation='vertical')
    plt.savefig("plots/neutral_amino_acid_count_per_codon.jpg")
    plt.show()

    return amino_acid_sequence


def possible_spike_protein_combinations(amino_acids):
    all_nucleotides = []

    for amino_acid in amino_acids:
        member_nucleotides = list(map(list, reverse_codon_wheel[amino_acid]))
        all_nucleotides.extend(list(itertools.chain(*member_nucleotides)))

    print("Total number of nucleotides present in the spike proteins:")
    print("Nucleotides: ", all_nucleotides)
    print("Count per nucleotide", collections.Counter(all_nucleotides))
    amino_acid_combinations = list(
        dict.fromkeys(["".join(x) for x in list(itertools.combinations(all_nucleotides, 3))]))
    print("Amino Acid from combinations", amino_acid_combinations)
    print("Length of above sequence", len(amino_acid_combinations))


def spike_protein_mutations(genome):
    codon = ""

    # for only one mutation
    no_of_sequecnces_one_mutation_away = 0
    for nucleotide in genome:
        if len(codon) < 3:
            codon += nucleotide.upper()
        else:
            non_neutral_mutations = calculate_non_neutral_mutations(codon)
            print(non_neutral_mutations)
            no_of_sequecnces_one_mutation_away = no_of_sequecnces_one_mutation_away + non_neutral_mutations
            codon = nucleotide.upper()

    # for 2-15 mutations away
    final = [no_of_sequecnces_one_mutation_away]
    codon = ""
    x = [1]
    for i in range(1, 5):
        possible_sequences = [genome]
        i += 1
        x.append(i)
        print("checking possibilities for ", i, " mutations away ")
        for _ in range(0, i):
            temp = []

            for sequence in possible_sequences:
                for nucleotide_index in range(0, len(sequence)):
                    # temp_seq=sequence
                    if len(codon) < 3:
                        codon += sequence[nucleotide_index].upper()
                    else:
                        non_neutral_mutations, mutated_codons = calculate_non_neutral_mutations_with_mutated_codons(
                            codon)
                        # no_of_sequecnces_one_mutation_away = no_of_sequecnces_one_mutation_away + non_neutral_mutations
                        codon = sequence[nucleotide_index].upper()
                        for new_codon in mutated_codons:
                            # temp_seq = sequence
                            # temp_seq[nucleotide_index-3] = new_codon[0]
                            # temp_seq[nucleotide_index - 2] = new_codon[1]
                            # temp_seq[nucleotide_index - 1] = new_codon[2]
                            temp.append(sequence[:nucleotide_index - 3] + codon + sequence[nucleotide_index:])
            possible_sequences = temp

        final.append(len(possible_sequences))
    # y = [log(y, 10) for y in final]
    # print(y)
    # plt.plot(x, y)
    # plt.xticks(x, x, rotation='vertical')
    # plt.savefig("plots/neutral_amino_acid_count_per_codon.jpg")
    # plt.show()

    return final

def question_4_log_values_for_non_neutral():
    #this method calculates the non neutral sequencec=s from 1 to 15 mutations away and writes log values in txt file
    
    f = open("logvalues.txt", "w+")

    possible_sequences=[]
    for i in range(1,16):
        print(i, )
        # comb = itertools.combinations([1,2,3,4,5,6,1,2,3,4,5,6,7,8,2], i)
        # possible_sequences.append(log(len(list(comb)*int(pow(3,i))),10))
        possible_sequences.append(log(int(pow(15,i)*pow(3, i)), 10))
        f.write("%f\n" % possible_sequences[-1])
    print(possible_sequences)

