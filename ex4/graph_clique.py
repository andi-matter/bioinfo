# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 14:41:31 2022

@author: erns_ae
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import sys

import time

# graph_file = "Weighted_pairwise-gene-clusters.txt"
# graph_file = "graph_eng.txt"
# graph_file = "graph_short.txt"
# graph_file = "graph_mid.txt"
# graph_file = "graph_7.txt"
# graph_file = "9606.protein.physical.links.v11.5.txt"

graph_file = sys.argv[1]

if graph_file[-4:] != ".txt":
    graph_file += ".txt"

start = time.time()

gtable = pd.read_csv(graph_file, sep='\s+', header=None, usecols=[0,1])

# print("read time", time.time() - start)


#----------------
id_dict = {} # look up index given id
index_dict = {} # look up id given index
unique_indices = 0 # count unique node names
adjacency_lists = [] 
degrees = [] # dont actually need this here, could use for optimising

start = time.time()

for index, row in gtable.iterrows():
    vert_A, vert_B = row
    
    try:
        ind_A = id_dict[vert_A][0] # find index for protein A
    except:
        id_dict[vert_A] = [unique_indices, 0] # if new protein, assing new index
        ind_A = unique_indices
        index_dict[ind_A] = vert_A
        adjacency_lists.append([]) # new protein, needs new adjacency list
        degrees.append(0) # and degree counter
        unique_indices += 1 # new protein found
        
    try:
        ind_B = id_dict[vert_B][0] # find index for protein B
    except:
        id_dict[vert_B] = [unique_indices, 0] # if new protein, assing new index
        ind_B = unique_indices
        index_dict[ind_B] = vert_B
        adjacency_lists.append([]) # new protein, needs new adjacency list
        degrees.append(0) # and degree counter
        unique_indices += 1 # new protein found
        
    # now just register that A and B are connected
    adjacency_lists[ind_A].append(ind_B)
    adjacency_lists[ind_B].append(ind_A)
    
    degrees[ind_A] += 1
    degrees[ind_B] += 1       

# print("loop time", time.time() - start)

#------------
start = time.time()

plot = False # postprocessing data save keeping for timings sake

if plot:
    adjacency_file = open("adjacency.txt", 'w+')

for list_ in adjacency_lists:
    list_.sort() # sort adjacency lists (not super sure if the numpy functions actually are faster if sorted but we can hope)
    if plot:
        adjacency_file.write(str(list_))
        adjacency_file.write("\n")
    
if plot:
    adjacency_file.close()

# print("sort time", time.time() - start)
#------------------



start = time.time()

## >  <
# ignore this, my work laptop doesnt have these symbols >:[

# cliques of size 2
clique_2 = []

for k in range(unique_indices):   
    for j in range(k, unique_indices):
        if j in adjacency_lists[k]:
            clique_2.append([k,j])

         

i = 2

condition = True

S_lists = [clique_2] # first bunch of cliques


stop = False
while not stop:
    
    i += 1     
    
    # S_i := ∅;
    S_i = []      
    S_i1 = S_lists[i-1-2] # this is S_i-1

    # for j := 1 to |S_i-1|   
    for j in range(1, len(S_i1)):
        
        # for k := j+1 to |S_i-1|
        for k in range(j+1, len(S_i1)):
            
            # T := S_i-1[j] ∩ S_i-1[k];
            T = np.array(S_i1[j])[np.in1d(S_i1[j], S_i1[k], assume_unique=True)]
            
            if len(T) == i-2:
                N = np.union1d(S_i1[j], S_i1[k])
                notN = np.setxor1d(S_i1[j], S_i1[k])
                
                # check if N is clique
                # N is clique if all elements connected
                # N made of two sub-cliques
                # just have to check if not-overlapping elements are connected to each other
                if notN[0] in adjacency_lists[notN[1]]:                   
                    S_i.append(N)
                    
    # delete doubled cliques from list
    S_i = np.unique(S_i, axis=0)
    
    stop =  len(S_i) == 0 # stop if no new cliques were found    
    if stop:
        break
    
    S_lists.append(S_i) # register new cliques
    
print(len(S_lists[-1][0]))
for clique in S_lists[-1]:
    #print(clique)
    print([index_dict[k] for k in clique])