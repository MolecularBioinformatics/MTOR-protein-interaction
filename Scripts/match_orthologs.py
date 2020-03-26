#!/env/usr/bin python3
"""
YH ~ 27.2.2020
[match_orthologs.py]

Takes two MSAs, identifies common organisms, returns
one sequence (longest/most complete) per organism in 
separate MSAs
"""

### IMPORTS
from sys import argv
import fastautils as fa 
from collections import OrderedDict as OD

### FUNCTIONS
def parse_id_for_org(listofids, idtype="uniprot"):
    """Takes a list of identifiers, returns
       list of organisms of origin"""
    orglist = []
    if idtype == "uniprot":
        for elem in listofids:
            shortid = elem.split()[0]
            orglist.append(shortid[-6:])
    return orglist

def find_common_orgs(dictofmsa1, dictofmsa2):
    """Takes in two msas to find common organisms
       returns two msa dicts with only common orgs"""
    commonorgs = []
    orgs1 = list(dictofmsa1.keys())
    orgs2 = list(dictofmsa2.keys())
    orglist1 = parse_id_for_org(orgs1)
    orglist2 = parse_id_for_org(orgs2)
    commonorgs = list(set(orglist1).intersection(orglist2))
    if not commonorgs:
        print('No common organisms between your MSAs')
    return commonorgs, orgs1, orgs2

def select_seqs(listofids, dictofmsa, listofcommonorgs):
    """Takes in an msa (dict), list of orgs, and a common list 
       returns revised msa with only one best choice seq per org in
       common orgs""" 
    newdict = OD()
    for org in listofcommonorgs:
        seqlist = []
        seqidlist = []
        for seqid in listofids:
            shortid = seqid.split()[0]
            if shortid[-6:] == org:
                seqlist.append(dictofmsa[seqid])
                seqidlist.append(seqid)
        bestseq = max(seqlist, key=len)
        bestseqid = seqlist.index(max(seqlist, key=len))
        newdict[seqidlist[bestseqid]] = dictofmsa[seqidlist[bestseqid]]
    return newdict

def main():
    msafile1 = argv[1]
    msafile2 = argv[2]
    msa_dict1 = fa.parse_fasta(open(msafile1, "r"))
    msa_dict2 = fa.parse_fasta(open(msafile2, "r"))
    common_orgs, org_list1, org_list2 = find_common_orgs(msa_dict1, msa_dict2)
    newmsa1 = select_seqs(org_list1, msa_dict1, common_orgs)
    newmsa2 = select_seqs(org_list2, msa_dict2, common_orgs)

    fa.write_fasta(newmsa1, f"{msafile1}_new.fa")
    fa.write_fasta(newmsa2, f"{msafile2}_new.fa")

### MAIN

if __name__ == "__main__":
    main()
