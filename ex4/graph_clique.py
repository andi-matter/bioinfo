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

graph_file = "Weighted_pairwise-gene-clusters.txt"
# graph_file = "graph_eng.txt"
# graph_file = "graph_short.txt"
# graph_file = "graph_mid.txt"
# graph_file = "graph_7.txt"
# graph_file = "9606.protein.physical.links.v11.5.txt"

# graph_file = sys.argv[1]

if graph_file[-4:] != ".txt":
    graph_file += ".txt"

start = time.time()

gtable = pd.read_csv(graph_file, sep='\s+', header=None, usecols=[0,1])

print("read time", time.time() - start)


#----------------
id_dict = {} # look up index given id
index_dict = {} # look up id given index
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
        index_dict[ind_A] = vert_A
        adjacency_lists.append([])
        degrees.append(0)
        unique_indices += 1
        
    try:
        ind_B = id_dict[vert_B][0]
    except:
        id_dict[vert_B] = [unique_indices, 0]
        ind_B = unique_indices
        index_dict[ind_B] = vert_B
        adjacency_lists.append([])
        degrees.append(0)
        unique_indices += 1
        
    #print(ind_A, ind_B)
        
    adjacency_lists[ind_A].append(ind_B)
    adjacency_lists[ind_B].append(ind_A)
    
    degrees[ind_A] += 1
    degrees[ind_B] += 1       

print("loop time", time.time() - start)

#------------
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

print("sort time", time.time() - start)
#------------------

# cliques of size 2

start = time.time()

## >  <
# ignore this, my work notebook doesnt have these symbols >:[

clique_2 = []

for k in range(unique_indices):
    
    for j in range(k, unique_indices):
        if j in adjacency_lists[k]:
            clique_2.append([k,j])
        # else: print("NO")
        
# print(len(clique_2))
# print(unique_indices)    

# print(adjacency_lists)   

i = 2

condition = True

S_lists = [clique_2]

max_deg = max(degrees)
print("max deg", max_deg)

stop = False



while not stop:
    
    i += 1
    
    if i == 10:
        stop = True
    
    
    # print(i) 
    
    # S_i := ∅;
    S_i = []
    
    
    
    S_i1 = S_lists[i-1-2]
    # print("S_i1", S_i1)
    
    if len(S_i1[-1]) > max_deg + 1:
        break # cant have cliques of bigger than max number of connections on one node + 1
    
    # for j := 1 to |S_i-1|
    
    last_N = []
    
    
    for j in range(1, len(S_i1)):
        
        # for k := j+1 to |S_i-1|
        for k in range(j+1, len(S_i1)):
            
            # print(S_i1)
            # print(S_i1[j])
            # print(S_i1[k])
            
            # print(np.in1d(S_i1[j], S_i1[k], assume_unique=True))
            
            # T := S_i-1[j] ∩ S_i-1[k];
            T = np.array(S_i1[j])[np.in1d(S_i1[j], S_i1[k], assume_unique=True)]

            # print(T)
            
            if len(T) == i-2:
                # print(S_i1[j], S_i1[k])
                N = np.union1d(S_i1[j], S_i1[k])
                notN = np.setxor1d(S_i1[j], S_i1[k])
                # print(len(N), notN)
                
                # if notN[0] == 435:
                # print("N ", N)
                # print("NOT N", notN)
                # print("adj", adjacency_lists[notN[1]])
                
                # check if N is clique
                # N is clique if all elements connected
                # N made of two sub-cliques
                # just have to check if not-overlapping elements are connected to each other
                if notN[0] in adjacency_lists[notN[1]]:
                    # print("S_i here" , S_i)
                    
                    already_in = False
                    
                    for last_N in S_i:
                        # print("ALL", np.all(np.isin(N, last_N)))
                        if np.all(np.isin(N, last_N)): 
                            already_in = True
                    if not already_in:
                        S_i.append(N)
                                       
            # print("\n j break")
        # print("\n k break \n") 
        
    # print("\n i break", S_i)
    
    stop =  len(S_i) == 0
    
    if stop:
        break
    
    S_lists.append(S_i)
    
    
# for list_ in S_lists:
    # print(list_)
   
print("cliq time", time.time() - start)
    
print("Maximum cliques: ")
for clique in S_lists[-1]:
    print(clique)
    print([index_dict[k] for k in clique])