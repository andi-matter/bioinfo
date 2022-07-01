import numpy as np
import sys
import pandas as pd

#fasta_sequence_file = "sequences.fasta"  # REPLACE FOR ARGV
#blosum_table = "BLOSUM62.txt"  # REPLACE FOR ARGV
fasta_sequence_file = sys.argv[1]
blosum_table = sys.argv[2]


# ------------------------------------------------------------------------------ FILE READING


if fasta_sequence_file[-6:] != ".fasta" or len(fasta_sequence_file) < 6:
    fasta_sequence_file += ".fasta"

## Read sequences/patterns from fasta file
with open(fasta_sequence_file, 'r') as i:
    lines = i.read().splitlines()

sequence_list = []
sequence_meta_list = []
sequence_length_list = []
current_sequence = []
string_index = -1
# count = 0
for line in open(fasta_sequence_file, 'r'):
    if line[0] == ">":
        if len(current_sequence) > 0:
            sequence_list.append(''.join(current_sequence))
        current_sequence = []
        sequence_length_list.append(0)
        string_index += 1
        sequence_meta_list.append(line[1:].strip())
        continue
    else:
        line = line.strip()
        current_sequence.append(line)
        sequence_length_list[string_index] += len(line)
sequence_list.append(''.join(current_sequence))

## Read scoring table
with open(blosum_table, 'r') as i:
    lines = i.read().splitlines()

score_matrix = []
alphabet = []
alphabet_bool = False

count = 0
for line in lines:
  count += 1
  if line[0] == "#" or line[0] == "":
      continue
  if not alphabet_bool:
    alphabet = line.replace(' ', '')
    alphabet_bool = True
    continue
  row = ''.join(line.replace(' ', ',').replace(',,', ',')[2:-1])
  row = np.fromstring(row, dtype=int, sep=',')
  score_matrix.append(row)
score_matrix = np.array(score_matrix)

# ------------------------------------------------------------------------------ PAIRWISE ALIGNMENT/SIMILARITY COST


def cost_pairwise(seq_1_index, seq_2_index):
    long_ind = seq_1_index
    short_ind = seq_2_index

    seq_short = np.char.array(list(sequence_list[short_ind]))
    seq_long = np.char.array(list(sequence_list[long_ind]))

    len_short = sequence_length_list[short_ind]
    len_long = sequence_length_list[long_ind]

    gap_penalty_pos = alphabet.index("*")

    # initialize alignment matrix
    alignment_matrix = np.zeros((len_short + 1, len_long + 1), dtype="int")  # alignment scores

    # create row /col header with respective gap costs
    init_long = np.zeros(len_long + 1)
    init_short = np.zeros(len_short + 1)

    for k in range(1, len_short + 1):
        init_short[k] = score_matrix[alphabet.index(seq_short[k - 1]), gap_penalty_pos]
    for j in range(1, len_long + 1):
        init_long[j] = score_matrix[alphabet.index(seq_long[j - 1]), gap_penalty_pos]

    alignment_matrix[0] = init_long
    alignment_matrix[:, 0] = init_short

    abs_max = 0  # keeping track of global maximum

    # print("Calculating alignment matrix...", end="")
    for k in range(len_short):
        for j in range(len_long):
            k_alph = alphabet.index(seq_short[k])
            j_alph = alphabet.index(seq_long[j])
            options = [alignment_matrix[k + 1, j] + score_matrix[j_alph, gap_penalty_pos],
                       alignment_matrix[k, j + 1] + score_matrix[k_alph, gap_penalty_pos],
                       alignment_matrix[k, j] + score_matrix[k_alph, j_alph]]
            # print(options)
            ind = np.argmax(
                options)  # which option was picked - if theres two equal maxima i just let numpy decide :) (think it takes the first one)
            loc_max = options[ind]  # what maximum actually was
            alignment_matrix[k + 1, j + 1] = loc_max

    cost = alignment_matrix[-1][-1]
    return cost


# ------------------------------------------------------------------------------ SIMILARITY MATRIX

# upper triangular matrix
def sim_matrix(sequence_list_1):
  len_list = len(sequence_list_1)
  matrix = np.zeros((len_list, len_list))
  for seq_1_index in range(1, len_list):
    for seq_2_index in range(seq_1_index):
      matrix[seq_2_index, seq_1_index] = cost_pairwise(seq_1_index, seq_2_index)

  return matrix

#OUTPUT SIMILARITY MATRIX
matrix = sim_matrix(sequence_list)
df = pd.DataFrame(matrix, columns=sequence_meta_list, index=sequence_meta_list)
print(df)

df.to_csv(r'output_guidetree.txt', sep=' ')

# ------------------------------------------------------------------------------ HIERARCHICAL CLUSTERING

def hierarchical_clustering3(matrix_1, sequence_list_1):
    len_list = len(sequence_list_1)
    merge_matrix = np.zeros((2 * len_list - 2, 2 * len_list - 2))  # there are len-1 merges but last is uninteresting

    name_indices = np.arange(2 * len_list - 2, dtype=int)
    #print("name_indices", name_indices)
    names = []
    for i in range(2 * len_list - 2):
        names.append([str((i + 1) * (i + 1 < len_list + 1))])
    #print("names", names)
    branches = []

    merge_matrix[0:matrix_1.shape[0], 0:matrix_1.shape[1]] += matrix_1  # fill old sim matrix into new
    merge_matrix += merge_matrix.T  # complete to not only have upper triangle

    for col_ind in range(len_list, 2 * len_list - 2):
        # print("#########################################")
        # print("#########################################")
        # print("col_ind", col_ind)
        # print(merge_matrix)
        # find max value indices
        max_row, max_col = np.array(np.where(np.max(merge_matrix) == merge_matrix))[:, 0]
        # print("max_row, max_col" ,max_row, max_col)
        # print("col", merge_matrix[:, max_col])
        # print("row", merge_matrix[:, max_row])

        # calculate new averages
        new_merge_matrix = np.copy(merge_matrix)
        new_merge_matrix[:, col_ind] = 0.5 * (merge_matrix[:, max_row] + merge_matrix[:, max_col])
        new_merge_matrix[col_ind, :] = 0.5 * (merge_matrix[max_row, :] + merge_matrix[max_col, :])
        merge_matrix = new_merge_matrix
        # print(merge_matrix)

        # zero all the values from indices that get merged
        merge_matrix[max_row, :] = 0
        merge_matrix[:, max_col] = 0
        merge_matrix[max_col, :] = 0
        merge_matrix[:, max_row] = 0
        # print(merge_matrix)

        # save what is merged
        bigger = 0
        smaller = 0
        if max_row < max_col:
            bigger = max_col
            smaller = max_row
        else:
            bigger = max_row
            smaller = max_col

        smaller_ind = np.where(name_indices == smaller)[0]
        # print("smaller", smaller_ind)
        # print(smaller_ind)
        bigger_ind = np.where(name_indices == bigger)[0]
        # print("bigger", bigger)
        # print(bigger_ind)
        name_small = names[smaller]
        name_big = names[bigger]

        current_merge = ""
        for name in [name_small, name_big]:
            current_merge += ", "
            if len(name) > 1:
                for j in range(len(name)):
                    current_merge += (name[j] + "+")
                current_merge = current_merge[:-1]
            else:
                current_merge += name[0]

        current_merge = "(" + current_merge[2:] + ")"
        print(current_merge)

        new_name = sorted(name_small + name_big)
        # print("new_name", new_name)
        names[col_ind] = new_name
        names[smaller] = ["0"]
        names[bigger] = ["0"]
        # print("names", names)
        current_branches = []
        for name in names:
            if name != ["0"]:
                current_branches += name
        # print("current_branches", current_branches)
        branches.append(current_branches)

        current_branches = sorted(current_branches)
        output = ""
        for i in range(len(current_branches)):
            branch = current_branches[i]
            output += ", "
            if len(branch) > 1:
                for j in range(len(branch)):
                    output += (branch[j] + "+")
                output = output[:-1]
            else:
                output += branch[0]
        output = "(" + output[2:] + ")"
        
        
        
        print(output)

#OUTPUT NEW MERGES AND CURRENT BRANCHES
hierarchical_clustering3(matrix, matrix[0])