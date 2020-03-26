#!/usr/env/bin/ python3
"""
[match_respairs.py]

Written for 3-column protein-protein DCA score files

Finds matching respairs (ids) and returns their coev vals
"""

### IMPORTS
from sys import argv
import numpy as np

### FUNCTIONS
def checkrespairs(respairarray1, respairarray2):
    """Takes in two np arrays of 3 columns, format: res1, res2, dca_value
       matches residue pairs and returns 4 column matrix with appended dca_vals """
    common_respairslist = []

    respairs1 = respairarray1[respairarray1[:,0].argsort()] #sort arrays by first columns
    respairs2 = respairarray2[respairarray2[:,0].argsort()]
    
    for res1, res2, val in respairs1:
        section = respairs2[np.where(respairs2[:,0] == res1)]
        for resi1, resi2, value in section:
            if res1 == resi1:
                if res2 == resi2:
                    common_respairslist.append([res1, res2, val, value])
                    break


    commonpairsarray = np.asarray(common_respairslist)

    return commonpairsarray

def main():
    respairsfile1 = argv[1]
    respairsfile2 = argv[2]
    A = np.loadtxt(respairsfile1, dtype=float, delimiter=" ")
    B = np.loadtxt(respairsfile2, dtype=float, delimiter=" ")
    CommonPairsArray = checkrespairs(A, B)
    print(CommonPairsArray)

    np.savetxt('common_respairs.dat', CommonPairsArray, delimiter=" ")

### MAIN
if __name__ == "__main__":
    main()
