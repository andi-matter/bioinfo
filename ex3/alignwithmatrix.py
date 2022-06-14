# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 15:35:20 2022

@author: uni
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import numpy as np
import sys
import time
# 

# fasta_sequence_file = "sequence_pair.fasta" # REPLACE FOR ARGV
# scoring_matrix_file = "BLOSUM62.txt"


fasta_sequence_file = sys.argv[1]
scoring_matrix_file = sys.argv[2]

#------------------------------------------------------------------------------ FILE READING

# START = time.time()

if fasta_sequence_file[-6:] != ".fasta" or len(fasta_sequence_file) < 6:
        fasta_sequence_file += ".fasta"

if scoring_matrix_file[-4:] != ".txt" or len(scoring_matrix_file) < 4:
        scoring_matrix_file += ".txt"


## Read sequences/patterns from fasta file 
with open(fasta_sequence_file, 'r') as i:
    lines = i.read().splitlines()

sequence_list = []
sequence_length_list = []
current_sequence = []
string_index = -1
# count = 0
for line in open(fasta_sequence_file, 'r'):
    # line = line.strip()
    # count += 1
    # progress = count % 10000
    # if 0 == progress:
        # print("Its still going... " + str(count))
    if line[0] == ">":
        if len(current_sequence) > 0:
            sequence_list.append(''.join(current_sequence))
        current_sequence = []
        # sequence_list.append("")
        sequence_length_list.append(0)
        string_index += 1
        continue
    else:
        line = line.strip()
        current_sequence.append(line)
        # sequence_list[string_index] += line
        sequence_length_list[string_index] += len(line)
sequence_list.append(''.join(current_sequence))




#------------- READ SCORING MATRIX jumping through many, many hoops

scoring_file = open(scoring_matrix_file, "r")
scoring_lines = scoring_file.readlines()

skip_header = 0

for line in scoring_lines:
    if line[0] == "#":
        skip_header += 1
        continue
    
    # print(line)
    
    skip_header += 1
    alphabet = [c for c in "".join(line.split())]
    ascii_alphabet = [ord(c) for c in "".join(line.split())]
    break


len_alph = len(alphabet)


scoring_matrix = np.asarray(np.genfromtxt(scoring_matrix_file, comments="#", names=alphabet, dtype=int, usecols=np.arange(1, len_alph+1), skip_header=skip_header))

scoring_matrix = scoring_matrix.view(int).reshape(scoring_matrix.shape + (-1, ))

gap_penalty = scoring_matrix[0, -1]


#-------------------- assign sequence entries index of respective alphabet entry in scoring matrix

num_sequences = []

for k in range(len(sequence_list)):
    num_sequ = np.array([alphabet.index(c) for c in sequence_list[k]])
    num_sequences.append(num_sequ)
    
# print(num_sequences)


long_ind = np.argmax(sequence_length_list)
short_ind = 1 - long_ind

seq_short = num_sequences[short_ind]
seq_long = num_sequences[long_ind]

len_short = sequence_length_list[short_ind]
len_long = sequence_length_list[long_ind]



#------------------------------------------------------------------------------ ALIGNMENT MATRIX

## initialising matrix
alignment_matrix = np.zeros((len_short+1, len_long+1), dtype="int")# alignment scores

init_long = gap_penalty * np.arange(0, len_long + 1)
init_short = gap_penalty * np.arange(0, len_short + 1)

alignment_matrix[0] = init_long
alignment_matrix[:, 0] = init_short

# print(alignment_matrix)



# print("Calculating alignment matrix...", end="")
for k in range(len_short):
    perc = int(k/len_short*100)
    
    char_short = seq_short[k]
    
    
    # if perc%10 == 5:
        # print(f"{perc}...", end ="")
    for j in range(len_long):
        
        char_long = seq_long[j]
        
        # print(char_short, char_long)
        
        
        # print(k,j, "cost",  cost_matrix[k, j])
        options = [alignment_matrix[k + 1, j] + gap_penalty,
                   alignment_matrix[k, j + 1] + gap_penalty,
                   alignment_matrix[k, j] + scoring_matrix[char_short, char_long]]
                    # definition a la lecture using scoring matrix entries
        # print(options)
              
        alignment_matrix[k+1, j+1] = max(options)
        


print(alignment_matrix[-1, -1])

