#!/usr/env/bin/ python3
"""
[split_inter_intra_dca.py]

Written for protein-protein DCA score files

Takes in 3 column text file, and length of your first protein,
splits the data into two files, one for inter and one for intra protein contacts
"""

### IMPORTS
from sys import argv
import numpy as np

### FUNCTIONS
def checkbounds(respairarray, protlen1, protlen2):
    """Takes in a 3 column np array of residue pairs, and the 
       protein lengths, to return only interprotein pairs"""
    interpairslist = []
    intrapairslist = []
    for res1, res2, val in respairarray:
        if (res1 <= protlen1 and res2 >= protlen1+1) or (res1 >= protlen1+1 and res2 <= protlen1):
            interpairslist.append([res1, res2, val])
        else:
            intrapairslist.append([res1, res2, val])
    interpairsarray = np.asarray(interpairslist)
    intrapairsarray = np.asarray(intrapairslist)
    return interpairsarray, intrapairsarray
    
def main():
    filename = argv[1]
    protein_length_1 = int(argv[2])
    protein_length_2= int(argv[3])
    A = np.loadtxt(filename, dtype=float, delimiter=" ")
    InterPairs, IntraPairs = checkbounds(A, protein_length_1, protein_length_2)
    np.savetxt('interprotscores.dat', InterPairs, delimiter=" ")
    np.savetxt('intraprotscores.dat', IntraPairs, delimiter=" ")

### MAIN
if __name__ == "__main__":
    main()
