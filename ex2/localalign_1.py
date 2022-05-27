# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import sys
import time

fasta_sequence_file = "KRAS_human_vs_mouse.fasta" # REPLACE FOR ARGV
# fasta_sequence_file = sys.argv[1]

#------------------------------------------------------------------------------ FILE READING

START = time.time()

if fasta_sequence_file[-6:] != ".fasta" or len(fasta_sequence_file) < 6:
        fasta_sequence_file += ".fasta"


## Read sequences/patterns from fasta file (this is basically the same as last time cause
## i could not for the life of me find the method that one dude mentioned)
with open(fasta_sequence_file, 'r') as i:
    lines = i.read().splitlines()

sequence_list = []
sequence_length_list = []
string_index = -1
count = 0
for line in open(fasta_sequence_file, 'r'):
    # line = line.strip()
    count += 1
    progress = count % 10000
    if 0 == progress:
        print("Its still going... " + str(count))
    if line[0] == ">":
        sequence_list.append("")
        sequence_length_list.append(0)
        string_index += 1
        continue
    else:
        line = line.strip()
        sequence_list[string_index] += line
        sequence_length_list[string_index] += len(line)


print("File read-in complete.") # maybe delete this b4 handin


#------------------------------------------------------------------------------ COST MATRIX


## cost matrix:
    ## unit cost model is M = 1, I/R/D = -1
    ## ones everywhere where same letter, -1 everywhere else
    
    
## defines some test strings, comment to use file stuff (see above)
# str1 = "AATCGGGATCGATAGCTACGATAGTGTACAGATACTCGCATAGCTA"
# str2 =        "GTCTGATACTACG"      
# sequence_list = [str1, str2]
# sequence_length_list = [len(str1), len(str2)]
## end test string stuff    

long_ind = np.argmax(sequence_length_list)
short_ind = 1 - long_ind

seq_short = np.char.array(list(sequence_list[short_ind]))
seq_long = np.char.array(list(sequence_list[long_ind]))

len_short = sequence_length_list[short_ind]
len_long = sequence_length_list[long_ind]

# print(len_short, len_long)


''' slow af cost matrix
#this is really slow!

start = time.time()

similarity_matrix = -1 * np.ones((len_short, len_long)) # init matrix filled with -1

for k in range(len_short):
    for j in range(len_long):
        if seq_short[k] == seq_long[j]: similarity_matrix[k, j] = 1
        
stop = time.time()

print(stop - start)
'''

## pretty fast cost matrix! could probs be faster with byte calcs but im lazy and dumb
cost_matrix = -np.ones((len_short, len_long))
cost_matrix[seq_short[:,None] == seq_long] = 1

# print(cost_matrix)


#------------------------------------------------------------------------------ ALIGNMENT MATRIX

## initialising of a bunch of matrices i need
alignment_matrix = np.zeros((len_short+1, len_long+1), dtype="int")# alignment scores
vector_matrix = np.zeros((len_short+1, len_long+1), dtype="int")    # directional info for max()
graphic_matrix = np.empty((len_short+1, len_long+1), dtype="str")  # or dont need but were useful

symbols = ["--", "||", "\\", ".."]

maxima = [] # gathering all global maxima in alignment matrix
curr_max = 0  # keeping track of global maximum

print("Calculating alignment matrix...", end="")
for k in range(len_short):
    perc = int(k/len_short*100)
    if perc%10 == 5:
        print(f"{perc}...", end ="")
    for j in range(len_long):
        # print(k,j, "cost",  cost_matrix[k, j])
        options = [alignment_matrix[k + 1, j] - 1,
                   alignment_matrix[k, j + 1] - 1,
                   alignment_matrix[k, j] + cost_matrix[k, j],
                   0] # definition a la lexture with unit cost for I/D/R
        # print(options)
        ind = np.argmax(options) # which option was picked - if theres two equal maxima i just let numpy decide :) (think it takes the first one)
        loc_max = options[ind]   # what maximum actually was
        alignment_matrix[k+1, j+1] = loc_max
        
        # update maximum tracking stuff
        if loc_max > curr_max:
            # print("greater than")
            maxima = [[k, j]]
            curr_max = loc_max
            
        elif loc_max == curr_max:
            # print("equal")
            maxima.append([k, j])
            curr_max = loc_max
        
        # print(alignment_matrix)
        
        
        vector_matrix[k+1, j+1] = ind # write down which direction we went to next cell
        graphic_matrix[k+1, j+1] = symbols[ind]
        
alignment_matrix[0, :] = 0 # reset first column and row to zeros (not sure why this is necessary rn but there was a reason earlier)
alignment_matrix[:, 0] = 0

 ##printf bugfixing
# print(alignment_matrix)
# print()
# print(vector_matrix)
# print()
# print(graphic_matrix)
# print()
# print(maxima)
# print()

# for k in range(len_short):
#     for j in range(len_long):
#         print(graphic_matrix[k,j], "\t", end="")
#     print()


TWEEN = time.time()
delta = TWEEN - START
print(f"Alignment matrix done after {delta/60:.4f} minutes. \n")

print(int(curr_max)) # KEEP

# ----------------------------------------------------------------------------- BACKWARD PATHFINDING


directions = [[0, -1], [-1, 0], [-1, -1], [0,0]] # cell walking directions depending on which option was chosen in alignment matrix

# maximum length of total alignment 
max_len = len_short + len_long
# print(max_len)

print()
print("Finding paths...")

count_max = 0

for k, j in maxima: # each maximum gets one alignment path
    print("path ", count_max, ":", end=" ")
    core_seq_short = [] # initialise/reset core sequences from last path
    core_seq_long = []
      
    core_seq_short.append(seq_short[k]) # the same index is NOT ACCIDENTAL!
    core_seq_long.append(seq_long[j]) # note down character at maximum
    
    score = alignment_matrix[k+1, j+1] # this is the maximum, i could also just write score = curr_max now that i think about it
    
    store_k = k # we need these for later
    store_j = j
    
    # print(core_seq_short)
    
    
    while score != 0:   # end path at 0 cell
        perc = 100 - int(k/(len_short)*100)
        if perc%15 == 0:
            print(f"{perc}...", end ="")
        direct = directions[vector_matrix[k, j]] # find out which direction to go from this cell (backward)

        k = k  + direct[0] # go the direction
        j = j  + direct[1]
        
        l = max(k,j)
        
        seq_short_next = seq_short[k] # get next character
        seq_long_next = seq_long[j]
        
        
        if direct[0] == 0: seq_short_next = "_" # unless direction was to not move, then insert blank
        if direct[1] == 0: seq_long_next = "_"  
        
        core_seq_short.append(seq_short_next) # write down character again INDEX INTENTIONAL
        core_seq_long.append(seq_long_next) 
              
        score = alignment_matrix[k, j] # see next score
    
    # print()
    # print("str1 ", str1)
    # print("str2 ", str2)
    # print()
    
    l = max(j, k)
    
    stars_after = max(len_short - store_k, len_long - store_j) - 1
  #  stars_before = max( )
      
    core_seq_long.append("*" * l)
    core_seq_short.append("*" * l)

    core_seq_long.reverse()
    core_seq_short.reverse()
    
    core_seq_long.append("*" * stars_after)
    core_seq_short.append("*" * stars_after)

    print()
    print(''.join(core_seq_long)) # KEEP join array to one string and print it
    print(''.join(core_seq_short)) # currently in the order they have in the pdf
    print()


STOP = time.time()
delta = STOP - START

print(f"\n Total calculation time {delta:.4f} seconds. That's {delta/60:.4f} minutes. Oof.")



