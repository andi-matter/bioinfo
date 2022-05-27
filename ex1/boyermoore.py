# Read console argument (here the file name / directory)
import sys
import numpy as np
fasta_sequence_file = sys.argv[1]
fasta_patterns_file = sys.argv[2]

if fasta_sequence_file[-6:] != ".fasta" or len(fasta_sequence_file) < 6:
        fasta_sequence_file += ".fasta"

if fasta_patterns_file[-6:] != ".fasta" or len(fasta_patterns_file) < 6:
        fasta_patterns_file += ".fasta"

# Bad character rule (BCR)
def preprocessing_bc(pattern_list, pattern_length_list, alphabet, len_alphabet):
    # Initialize table for every pattern
    shift_table_list = []
    for pattern_ind in range(len(pattern_list)):
        pattern = pattern_list[pattern_ind]
        len_pattern = pattern_length_list[pattern_ind]
        # Shift-Tabelle:
        # Zeilenindex <-> Buchstabe
        # Spalteindex <-> Position im Pattern
        shift_table_row = np.arange(1, len_pattern + 1)
        shift_table = []
        for row in range(len_alphabet):
            shift_table.append(shift_table_row)
        shift_table = np.array(shift_table)
        # print(shift_table)
        # Loop over every possible letter
        for letter_ind in range(len_alphabet):
            letter = alphabet[letter_ind]
            # For each given letter loop from right to left through every table column/ pattern position
            for position in range(len_pattern - 1, -1, -1):
                # Keep track of how far until the next letter occurence
                count = 0
                if letter != pattern[position]:
                    # Loop from every (start) position to left end or until match and continue outer loop
                    for new_position in range(position - 1, -1, -1):
                        count += 1
                        if letter == pattern[
                            new_position]:  # (This exact table entry will never be requested since bcr is only acivated for mismatch)
                            # Save count only if the letter occurs again
                            shift_table[letter_ind, position] = count
                            break

        shift_table_list.append(shift_table)
    return shift_table_list

# Good suffix rule (GSR)
def preprocessing_gs(pattern_list, pattern_length_list, alphabet, len_alphabet):
    # For every possible suffix and each letter of the alphabet the first reoccurence (from left to right) position
    #  of that suffix with the respective preceding letter will listed

    # Initialize table for all patterns
    suffix_table_list = []
    prefix_suffix_table_list = []
    for pattern_ind in range(len(pattern_list)):
        pattern = pattern_list[pattern_ind]
        len_pattern = pattern_length_list[pattern_ind]

        # Table for one pattern
        #   Row <-> Letter of alphabet
        #   Col <-> Suffix
        #   Entry <-> Shift to first reoccurence of suffix where letter of row precedes
        suffix_table = np.zeros([len_alphabet, int((len_pattern - 2) / 2) - 1], dtype=np.int8)

        # All possible suffixes are examined
        # 1st suffix -> i=0 with len=2, 2nd -> i=1 with len=3, ... -> this is unique, i.e. the suffix str must not be saved
        for suffix_ind in range(int((len_pattern - 2) / 2) - 1):  # max suffix len = (n-2)/2, min len = 2

            suffix_ind += 2  # Now index and length are equal
            suffix_len = suffix_ind
            suffix = pattern[-suffix_ind:]  # Current suffix
            # print("Suffix", suffix)
            suffix_positions = []
            suffix_following_letters = []

            # Start searching for duplicate suffixes next to original
            for substring_index in range(len_pattern - 1 - suffix_len, suffix_len - 2,
                                         -1):  # substring_index = index of last substring character
                substring = pattern[substring_index - suffix_len + 1:substring_index + 1]

                # print("substring", substring)

                # Save position for first match with each preceding letter from alphabet
                if suffix == substring:
                    pre_letter = pattern[substring_index - suffix_len]
                    pre_letter_ind = alphabet.index(pre_letter)
                    if 0 == suffix_table[pre_letter_ind, suffix_ind - 2]:
                        # actually dont save position but its difference to save the neceassary shift
                        suffix_table[pre_letter_ind, suffix_ind - 2] = len_pattern - 1 - substring_index

        suffix_table_list.append(suffix_table)

    return suffix_table_list

# Boyer Moore algorithm with BCR and GSR
def boyermoore(pattern_list, pattern_length_list, alphabet, len_alphabet, sequence, len_sequence):

    # Create shift table for bad character rule for each pattern
    shift_table_bc_list = preprocessing_bc(pattern_list, pattern_length_list,
                                           alphabet, len_alphabet)
    # Shift-Tabelle:
    # Zeilenindex <-> Buchstabe
    # Spalteindex <-> Position im Pattern
    #  Entry <-> Shift to "good character"
    # print(shift_table_bc_list)

    # Create "shift" table for good suffix rule for each pattern
    shift_table_gs_list = preprocessing_gs(pattern_list, pattern_length_list,
                                           alphabet, len_alphabet)
    #   Row <-> Letter of alphabet
    #   Col <-> Suffix
    #   Entry <-> Shift to first reoccurence of suffix where letter of row precedes

    match_positions_list = []
    bcr_list = []  # keep track of bcr count
    # Loop over all patterns
    for pattern_ind in range(len(pattern_list)):
        #print("Pattern no.", pattern_ind + 1, "/", len(pattern_list))
        pattern = pattern_list[pattern_ind]

        shift_table_bc = shift_table_bc_list[pattern_ind]

        shift_table_gs = shift_table_gs_list[pattern_ind]

        len_pattern = pattern_length_list[pattern_ind]
        max_suffix_len = int((len_pattern - 2) / 2)

        match_positions = []
        bcr_count = 0
        match_count = 0
        seq_position = 0  # Position in sequence refers to where the left most letter of the pattern is compared to
        # Check along entire sequence
        last_seq_position = len_sequence - len_pattern
        while seq_position <= last_seq_position:
            len_match = 0
            # Loop over pattern backwards until complete match/ single mismtach with seq substring
            for letter_pos_pattern in range(len_pattern - 1, -1, -1):
                # Sinlge letter match
                if pattern[letter_pos_pattern] == sequence[seq_position + letter_pos_pattern]:
                    len_match += 1
                    # Whole substring match!!!
                    if len_match == len_pattern:
                        match_count += 1
                        match_positions.append(seq_position + 1)  # +1 to convert from index to position
                        shift = 1  # for a large alphabet this could be improved (e.g. shift to next pattern occurence of right most (1st) matching letter)
                        seq_position += shift
                        
                # Mismatch: find out how far to shift pattern to the right
                else:
                    # Bad character rule (bcr)
                    mismatch_letter = sequence[seq_position + letter_pos_pattern]
                    mismatch_letter_ind = alphabet.index(mismatch_letter)
                    shift_bc = shift_table_bc[mismatch_letter_ind, letter_pos_pattern]

                    # Good suffix rule (gsr)
                    shift_gs = 0
                    if 2 <= len_match <= max_suffix_len:
                        shift_gs = shift_table_gs[mismatch_letter_ind, len_match - 2]
                    # shift_gs = 0

                    # Is bcr or gsr (and which more) useful?
                    if shift_bc > 1 and shift_bc >= shift_gs:
                        bcr_count += 1
                    seq_position += max(shift_bc, shift_gs, 1)
                    break

        match_positions_list.append(match_positions)
        bcr_list.append(bcr_count)

    return match_positions_list, bcr_list


#print("Search pattern apperances in:", fasta_sequences_file)
#print("Use patterns from:", fasta_patterns_file)

# Read sequences/patterns from fasta file
with open(fasta_patterns_file, 'r') as i:
    lines = i.read().splitlines()

pattern_list = []
pattern_length_list = []
string_index = -1
count = 0
for line in lines:
    count += 1
    progress = count % 10000
    if 0 == progress:
        print("Its still going... " + str(count))
    if line[0] == ">":
        pattern_list.append("")
        pattern_length_list.append(0)
        string_index += 1
        continue
    else:
        pattern_list[string_index] += line
        pattern_length_list[string_index] += len(line)



sequence = ""
with open(fasta_sequence_file) as infile:
    for line in infile:
        sequence += line.strip()

# strip first comment line
with open(fasta_sequence_file) as f:
    first_line = f.readline().strip()
len_first = len(first_line)
sequence = sequence[len_first:]
len_sequence = len(sequence)


alphabet = ["a", "c", "g", "t", "n"]
len_alphabet = len(alphabet)
match_positions_list, bcr_list = boyermoore(pattern_list, pattern_length_list, alphabet, len_alphabet, sequence, len_sequence)

# Print results
for pattern_ind in range(len(pattern_list)):
    output = ""
    output += str(bcr_list[pattern_ind]) + " / "
    output += str(len(match_positions_list[pattern_ind])) + " /"
    match_positions = match_positions_list[pattern_ind]
    for occurence in range(min(10, len(match_positions_list[pattern_ind]))):
        output += " " + str(match_positions[occurence]) + " /"
    print(output[:-1])

sys.exit()
