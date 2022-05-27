import sys


if __name__=="__main__":
    # Read console argument (here the file name / directory)
    fasta_file_name = sys.argv[1]

    if fasta_file_name[-6:] != ".fasta" or len(fasta_file_name)<6:
        fasta_file_name += ".fasta"

    #print(fasta_file_name)
    # Open fasta file and read line by line

    with open(fasta_file_name, 'r') as i:
        lines = i.read().splitlines()

    string_list = []
    string_length_list = []
    string_index = -1
    for line in lines:
        if line[0] == ">":
            string_list.append("")
            string_length_list.append(0)
            string_index += 1
            continue
        else:
            string_list[string_index] += line
            string_length_list[string_index] += len(line)

    for len in string_length_list:
        print(str(len))


    """output = ""
    for len in string_length_list:
        output += str(len) + " | """

    #print()

    sys.exit()
    #input("Press any key to continue.")
