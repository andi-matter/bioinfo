# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 10:59:19 2022

@author: erns_ae
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import sys

import time




graph_file = "graph_short.txt"
# graph_file = "graph_mid.txt"
# graph_file = "graph_7.txt"
# graph_file = "9606.protein.physical.links.v11.5.txt"

# graph_file = sys.argv[1]

if graph_file[-4:] != ".txt":
    graph_file += ".txt"

start = time.time()

gtable = pd.read_csv(graph_file, sep=" ", header=None, usecols=[0,1])

print("read time", time.time() - start)

# print(gtable)

id_dict = {}

unique_indices = 0

adjacency_lists = []

degrees = []

start = time.time()

for index, row in gtable.iterrows():
    vert_A, vert_B = row
    
    try:
        ind_A = id_dict[vert_A][0]
    except:
        id_dict[vert_A] = [unique_indices, 0]
        ind_A = unique_indices
        adjacency_lists.append([])
        degrees.append(0)
        unique_indices += 1
        
    try:
        ind_B = id_dict[vert_B][0]
    except:
        id_dict[vert_B] = [unique_indices, 0]
        ind_B = unique_indices
        adjacency_lists.append([])
        degrees.append(0)
        unique_indices += 1
        
    #print(ind_A, ind_B)
        
    adjacency_lists[ind_A].append(ind_B)
    adjacency_lists[ind_B].append(ind_A)
    
    degrees[ind_A] += 1
    degrees[ind_B] += 1
        

print("loop time", time.time() - start)

#print(hash_lists)

#degrees = pd.DataFrame.from_dict(id_dict, orient="index")[1].to_numpy()

#print(degrees) 
#print(unique_indices)

degree_binned, bin_edges = np.histogram(degrees, bins=np.arange(1, max(degrees)+2))


# PUT IN FOR HANDIN
# for k in range(len(degree_binned)):
    # print(degree_binned[-k-1], end=" ")


start = time.time()

plot = False

if plot:
    adjacency_file = open("adjacency.txt", 'w+')



for list_ in adjacency_lists:
    list_.sort()
    if plot:
        adjacency_file.write(str(list_))
        adjacency_file.write("\n")
    
if plot:
    adjacency_file.close()

# print(adjacency_lists)


print("sort time", time.time() - start)

#print(np.sort(np.array(adjacency_lists)))


#print(rev_deg_counts)
# print(rev_deg)



if plot:
    
    rev_deg_counts = degree_binned[::-1]
    rev_deg = bin_edges[-2::-1]

    
    np.savetxt("degrees_unbinned.txt", degrees, fmt="%i")
    
    plt.scatter(rev_deg, rev_deg_counts)
    plt.show()
    plt.hist(degrees)