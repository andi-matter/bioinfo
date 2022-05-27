# -*- coding: utf-8 -*-
"""
Created on Tue May  3 18:56:18 2022

@author: uni
"""

## task 1, fastaread

import io
import sys

def fastaread(filename):
    
    #print(filename[-6:])
    
    if filename[-6:] != ".fasta" or len(filename) < 6:
        filename += ".fasta"
    
    f = io.open(filename, "r", encoding="utf-8")
    #print("EHRE")

    #content = f.read()
    #print(content)
    count = 0
    segments = 0
    started = False
    sequence = ""
    
    for line in f:
        
        stripped = line.strip()
        
        if len(stripped) == 0:
            continue
        
        if stripped[0] == ">":
            started = True   
            if segments > 0:
                print(count)
                #print(sequence)
                #print("--------")
            
            segments += 1
            
            count = 0
            sequence = ""
                
            #print(stripped)
            
        else:
            if not started: continue
            count += len(stripped)
            sequence += stripped
            
    print(count)
    #print(sequence)
    
    f.close()
    
    return
    
if __name__ == "__main__":
    filename = "pattern.fasta" #sys.argv[1]
        
    fastaread(filename)