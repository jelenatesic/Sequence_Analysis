# Assignment: DNA Sequence Analysis

# Background:
# DNA sequence analysis is a fundamental task in bioinformatics.
# In this assignment, i will work on a project to analyze DNA sequences and perform various bioinformatics tasks.

# __________________________________________________________________________________

# 1. Sequence Input: Write a Python program that prompts the user to input a DNA sequence.
# Ensure that the program validates the input to ensure it contains valid DNA nucleotides (A, T, C, G).
# If the input is not valid, the program should prompt the user to re-enter the sequence until a valid one is provided.

# 2. Sequence Information: Implement a function that takes a DNA sequence as input and displays the following information
# about the sequence:
# Length of the sequence
# Count of each nucleotide (A, T, C, G)
# GC content (percentage of G and C nucleotides)

# 3. Sequence Complement: Write a function that takes a DNA sequence as input and returns its complement sequence.
# The complement of a DNA sequence is formed by replacing each nucleotide with its complementary nucleotide
# (A with T, T with A, C with G, and G with C).

# 4. Sequence Translation: Implement a function that translates a DNA sequence into its corresponding protein sequence.
# Use the standard genetic code to map codons (triplets of nucleotides) to amino acids. The translation should start
# from the first codon (ATG) and continue until a stop codon is encountered.

# 5. Sequence Alignment: Develop a function that performs pairwise alignment of two DNA sequences using the Needleman-Wunsch
# algorithm or any other alignment algorithm of your choice. The function should display the aligned sequences, including
# gaps if necessary, along with the alignment score.


# _______________________________________ functions ___________________________________________________________________

# *** Sequence input
def sequence_input():

    """
    The putpose of this function is to take imported sequence, check
    if it is in the right format, and returne it.

    :return sequence in list
    """
    m = 1
    while m:
        m = 0
        seq = input()

        seq_list = list(seq)

        for i in seq_list:
            if i == '\n':
                seq_list.remove(i)



        for i in range(0, len(seq_list)):
            if seq_list[i].lower() not in ["a", "g", "c", "t"]:
                print("not a nucleotide", seq_list[i])
                m = 1
                break

    # print('Your seq is in correct forme: ', seq)
    return seq_list


# *** Sequence information
def seq_info(seq_list):

    """
    The purpose of this function is to take a sequence and give some basic information
    about it:
    Length of the sequence
    Count of each nucleotide (A, T, C, G)
    GC content (percentage of G and C nucleotides)

    :param seq_list:
    :return:
    """

    a = 0
    g = 0
    c = 0
    t = 0

    for i in range(0,len(seq_list)):
        if seq_list[i].lower() == 'a':
            a += 1
        elif seq_list[i].lower() == 'g':
            g += 1
        elif seq_list[i].lower() == 'c':
            c += 1
        else:
            t += 1
    print("Length of your sequence is: ", len(seq_list))

    print(f"In your sequence, there are: \n{a} of A, \n{g} of G, \n{c} of C, \n{t} of T")
    print(f"Percentage of G in sequence is {np.round(g/len(seq_list)*100, 2)}%")
    print(f"Percentage of C in sequence is {np.round(c/len(seq_list)*100, 2)}%")


# *** Complement Sequence
def compliment_seq(seq_list):

    """
    The purpose of this function is to take a DNA sequence and return its complement
    sequence.

    :param seq_list:
    :return: complement sequence
    """

    seq_compl_list = []
    for i in seq_list:
        if i.lower() == 'a':
            seq_compl_list.append('t')
        elif i.lower() == 'g':
            seq_compl_list.append('c')
        elif i.lower() == 'c':
            seq_compl_list.append('g')
        else:
            seq_compl_list.append('a')

    seq_compl = "".join(seq_compl_list).upper()
    f"For the imported sequence {''.join(seq_list)}, the complement is {seq_compl}"
    return seq_compl



# 4. Sequence translation function
def transtating_dna(seq_list):
    """
    This function takes the DNA sequence and gives back the translated
    sequence, using codon maps

    :param seq_list:
    :return: translated list
    """
    M = ["ATG"]   # Methionine
    F = ["TTT", "TTC"]   # Phenylalanine
    L = ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"] # Leucine
    I = ["ATT", "ATC", "ATA"]  # Isoleucine
    V = ["GTT", "GTC", "GTA", "GTG"]   # Valine
    S = ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"] # Serine
    P = ["CCT", "CCC", "CCA", "CCG"]   # Proline
    T = ["ACT", "ACC", "ACA", "ACG"]   # Threonine
    A = ["GCT", "GCC", "GCA", "GCG"]   # Alanine
    Y = ["TAT", "TAC"]   # Tyrosine
    STOP = ["TAA", "TAG", "TGA"]
    H = ["CAT", "CAC"]  # Histidine
    Q = ["CAA", "CAG"]  # Glutamine
    N = ["AAT", "AAC"]  # Asparagine
    K = ["AAA", "AAG"]  # Lysine
    D = ["GAT", "GAC"]  # Aspartate
    E = ["GAA", "GAG"]  # Glutamate
    C = ["TGT", "TGC"]  # Cysteine
    W = ["TGG"]     # Tryptophan
    R = ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"]  # Arganine
    G = ["GGT", "GGC", "GGA", "GGG"]   # Glycine

    proteine_list = {"M": M, "F": F, "L": L, "I":I,"V": V,"S": S,"P": P,"T": T,"A": A,"Y": Y,"H": H,"Q": Q,
                     "N": N,"K": K,"D": D,"E": E,"C": C,"W": W,"R": R,"G": G,"*": STOP}


    # ***Finding start codon M
    for i in range(0,len(seq_list)-3):
        if [f"{seq_list[i]}{seq_list[i+1]}{seq_list[i+2]}"] == M:
            start_index = [i, i+1, i+2]
            break


    # print(start_index)

    # ***Making new list with triplets from start codon
    seq_list_triplets_all = []

    for i in range(start_index[0], (len(seq_list)), 3):
        seq_list_triplets_all.append(f"{seq_list[i]}{seq_list[i+1]}{seq_list[i+2]}")

    # print(seq_list_triplets_all)

    # ***Finding stop codon
    # seq_list_triplets = []
    # for i in seq_list_triplets_all:
    #     if i not in STOP:
    #         seq_list_triplets.append(i)
    #     else:
    #         seq_list_triplets.append(i)
    #         break
    #
    # print(seq_list_triplets)


    # ***Translating
    seq_translated = []

    for triplet in seq_list_triplets_all:
        for proteine_name, proteine in proteine_list.items():
            if triplet in proteine:
                seq_translated.append(proteine_name)

    # print("".join(seq_translated))

    return "".join(seq_translated)


# *** Sequence Alignment
def sequence_alignment():

    """
    The purpose of this function is to do a Needleman-Wunsch algorithm
    for Global alignment, returning the best alignment for the two
    inported sequences.

    :return: alignmen_1 and alignment_2
    """

    seq1 = input('seq1: ')
    seq2 = input('seq2: ')

    # Use these values to calculate scores
    import numpy as np


    gap_penalty = -1
    match_award = 1
    mismatch_penalty = -1

    # Defining empty matrix
    s1_len = len(seq1)
    s1_len
    s2_len = len(seq2)
    s2_len

    # Defining list
    seq1_list = list(seq1)
    seq1_list
    seq2_list = list(seq2)
    seq2_list

    # *** Making zero matrixes
    main_matrix = np.zeros(shape= (s2_len+1, s1_len+1))
    main_matrix
    match_checker_matrix = np.zeros((s2_len, s1_len))
    match_checker_matrix

    # *** Filling the match checker matrix accroding to match or missmatch
    for i in range(s2_len):
        for j in range(s1_len):
            if seq1_list[j] == seq2_list[i]:
                match_checker_matrix[i][j] = match_award
            else:
                match_checker_matrix[i][j] = mismatch_penalty

    match_checker_matrix


    # *** Filling out the first row and column of the matrix using the gap penalty.
    for i in range(0, s2_len):
        main_matrix[i+1][0] = main_matrix[i][0] + gap_penalty

    for j in range(0, s1_len):
        main_matrix[0][j + 1] = main_matrix[0][j] + gap_penalty

    main_matrix

    # **** Fill out the rest of the matrix using a nested loop, for every cell calculate the following:
    # match = matrix[i-1, j-1] + match_score(seq1[j-1], seq2[i-1])
    # delete = matrix[i - 1][j] + gap_penalty
    # insert = matrix[i][j - 1] + gap_penalty


    for i in range(1, main_matrix.shape[0]): # red
        for j in range(1, main_matrix.shape[1]): # kolona

            F_diag = main_matrix[i-1][j-1] + match_checker_matrix[i-1][j-1]

            F_horizontal = main_matrix[i][j-1] + gap_penalty

            F_vertical = main_matrix[i-1][j] + gap_penalty

            # Finding the max value
            Fij = max(F_diag, F_vertical, F_horizontal)

            # Filling the matrix
            main_matrix[i][j] = Fij


    main_matrix

    # ***** Traceback

    aligned_1 = ""
    aligned_2 = ""

    # Length of sequences
    ti = len(seq2_list) # for j
    tj = len(seq1_list) # for i

    while(ti > 0 and tj > 0):

        # When they are aligned
        if(ti >0 and tj > 0 and main_matrix[ti][tj] == main_matrix[ti-1][tj-1] + match_checker_matrix[ti-1][tj-1]):

            aligned_1 = seq1_list[tj-1] + aligned_1
            aligned_2 = seq2_list [ti-1] + aligned_2

            ti -= 1
            tj -= 1

        # When there is vertical movement, along i-axis, gap in seq_1
        elif (ti > 0 and main_matrix[ti][tj] == main_matrix[ti-1][tj] + gap_penalty): # vertikalno pomeranje

            aligned_1 = "_" + aligned_1
            aligned_2 = seq2_list[ti-1] + aligned_2 # gap is in the seq which did not change during move

            ti -= 1

        # When there is horizontal movement along j-axis, gap in seq_2
        else:
            aligned_1 = seq1_list[tj-1] + aligned_1
            aligned_2 = "_" + aligned_2

            tj -= 1


    return np.array([[aligned_1], [aligned_2]])

# ________________________________________ testiing _______________________________________________

# 1. Inputing seq
sequence_list_1 = sequence_input()

# 2. Sequence info
seq_info(sequence_list_1)

# 2. Complement sequence
compliiment_sequence_1 = compliment_seq(sequence_list_1)
print(f"complement for the seq \n {''.join(sequence_list_1)} is \n {compliiment_sequence_1}")

# 4. Translated seq
tranlated_seq = transtating_dna(sequence_list_1)
tranlated_seq

# 5. Sequence alignment
# seq_1 = ACTGAT
# seq_2 = ACATAGCT
sequence_alignment()


