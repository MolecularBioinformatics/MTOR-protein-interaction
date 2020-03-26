"""
Module containing functions for dealing with Bio.PDB objects, create contact
matrices, and map contact matrices to multiple sequence alignments
"""

import numpy as np

from collections import OrderedDict
from Bio import SeqIO
import warnings

global THREE_TO_ONE
THREE_TO_ONE = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N',
                'ASP': 'D', 'CYS': 'C', 'GLU': 'E',
                'GLN': 'Q', 'GLY': 'G', 'HIS': 'H',
                'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
                'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                'SER': 'S', 'THR': 'T', 'TRP': 'W',
                'TYR': 'Y', 'VAL': 'V'}

global LEGAL_AAS
LEGAL_AAS = ['ALA', 'ARG', 'ASN', 'ASP',
             'CYS', 'GLU', 'GLN', 'GLY',
             'HIS', 'ILE', 'LEU', 'LYS',
             'MET', 'PHE', 'PRO', 'SER',
             'THR', 'TRP', 'TYR', 'VAL']


#############################
# Bio.PDB object processing #
#############################

def remove_water(chain):
    """
    Given a protein chain, remove waters forming part of the sequence.

    Arguments
    ---------
    chain:  Bio.PDB.Chain.Chain object

    Returns
    -------
    filtered_chain: list of Bio.PDB.Residue.Residue, with no "water residues"
    """
    filtered_chain = []
    for res in list(chain):
        if res.__dict__['resname'] == 'HOH':
            pass
        else:
            filtered_chain.append(res)
    return filtered_chain


def remove_non_legal(chain, legal_aas=LEGAL_AAS):
    """
    Given a protein chain, remove non-legal amino acids forming part of the
    sequence.
    This is intended to remove non-amino acids from the chain.

    Arguments
    ---------
    chain:  Bio.PDB.Chain.Chain object

    Returns
    -------
    filtered_chain: list of Bio.PDB.Residue.Residue
    """
    filtered_chain = []
    for res in list(chain):
        if res.__dict__['resname'] in legal_aas:
            filtered_chain.append(res)
    return filtered_chain


def three_letter_to_one_letter(chain, code=THREE_TO_ONE):
    """
    Translate' a chain of amino acids in three letter code to one letter code.
    Note that the returned object contains no structural information: it is
    only the sequence of the protein.

    Arguments
    ---------
    chain: Bio.PDB.Chain.Chain object or list of Bio.PDB.Residue.Residue
    code:  dict, conversion table

    Returns
    -------
    translated_chain: string
    """
    translated_chain = []
    for res in list(chain):
        try:
            translated_chain.append(THREE_TO_ONE[res.__dict__['resname']])
        except KeyError:
            warnings.warn(f'Unknown amino acid encountered: {res}, skipping', RuntimeWarning)
    return ''.join(translated_chain)


###########################
# Contact matrix building #
###########################

def calc_residue_distance(residue_one, residue_two):
    """
    Given two residues, calculate the Euclidean distance between them.
    The distance is measured between the beta carbons (or, in the case of
    glycine, with respect to the alpha carbon).

    Arguments
    ---------
    residue_one: Bio.PDB.Residue.Residue object
    residue_two: Bio.PDB.Residue.Residue object

    Returns
    -------
    euclid_dist:  float, Euclidean distance between the two residues
    """
    is_one_glycine = (residue_one.__dict__['resname'] == 'GLY')
    is_two_glycine = (residue_two.__dict__['resname'] == 'GLY')

    testCA_one = "CA" in residue_one.__dict__['child_dict']
    testCA_two = "CA" in residue_two.__dict__['child_dict']

    testCB_one = "CB" in residue_one.__dict__['child_dict']
    testCB_two = "CB" in residue_two.__dict__['child_dict']

    problemresidue = False

    if is_one_glycine and is_two_glycine and testCA_one and testCA_two:
        diff_vector = residue_one["CA"].coord - residue_two["CA"].coord
        #print('{}  {}'.format(residue_one.__dict__['_id'][1], residue_two.__dict__['_id'][1])) 
    elif is_one_glycine and not is_two_glycine and testCB_two and testCA_one:
        diff_vector = residue_one["CA"].coord - residue_two["CB"].coord
        #print('{}  {}'.format(residue_one.__dict__['_id'][1], residue_two.__dict__['_id'][1])) 
    elif not is_one_glycine and is_two_glycine and testCB_one and testCA_two:
        diff_vector = residue_one["CB"].coord - residue_two["CA"].coord
        #print('{}  {}'.format(residue_one.__dict__['_id'][1], residue_two.__dict__['_id'][1])) 
    elif testCB_one and testCB_two:
        diff_vector = residue_one["CB"].coord - residue_two["CB"].coord
        #print('{}  {}'.format(residue_one.__dict__['_id'][1], residue_two.__dict__['_id'][1])) 
    else:
        problemresidue = True
        diff_vector = 100 # needs to be high enough to not count, set at 100
        print('{}  {}'.format(residue_one.__dict__['_id'][1], residue_two.__dict__['_id'][1])) 
    
    #if problemresidue:
    #    print('=========== Missing Carbon Problem ===========')


    euclid_dist = np.sqrt(np.sum(diff_vector * diff_vector))
    return euclid_dist


def calc_dist_matrix(chain_one, chain_two):
    """
    Given two proteins chains from a PDB file, calculate a matrix containing
    all pairwise Euclidean distances between residues from the two different
    chains.

    Arguments
    ---------
    chain_one, chain_two:  Bio.PDB.Chain.Chain object, or
                           list of Bio.PDB.Residue.Residue

    Returns
    -------
    dist_mtx:   array-like, matrix containing pairwise Euclidean distances
                between residues, of dimensions (chain_one, chain_two)
    """
    dist_mtx = np.zeros((len(chain_one), len(chain_two)))
    for i, residue_one in enumerate(chain_one):
        for j, residue_two in enumerate(chain_two):
            dist_mtx[i,j] = calc_residue_distance(residue_one, residue_two)

    return dist_mtx

def return_residue_ids(chain_one, chain_two):
    """
    Given two protein chains from a PDB file, return the residue ID numbers
    as given in the PDB file corresponding to both chains, 
    useful for later visualisation if original indices are needed.

    Arguments
    ---------
    chain_one, chain_two: Bio.PDB.Chain.Chain object, or
                          list of Bio.PDB.Residue.Residue

    Returns
    -------
    resi_ids_ch1/2: dictionaries of the numeric indices paired with PDB 
    indices
    """
    resi_ids_ch1 = {}
    resi_ids_ch2 = {}
    for i, residue_one in enumerate(chain_one):
        resi_ids_ch1[i] = residue_one.__dict__['_id'][1]
    for j, residue_two in enumerate(chain_two):
        resi_ids_ch2[j] = residue_two.__dict__['_id'][1]

    return resi_ids_ch1, resi_ids_ch2
   

def get_contact_matrix(dist_mtx, threshold):
    """
    Discretize distance matrix according to a given threshold.
    This is done with the purpose of creating a contact map.

    Arguments
    ---------
    dist_mtx:   array-like, matrix containing pairwise Euclidean distances
                between residues, of dimensions (chain_one, chain_two)
    Returns
    -------
    contact_mtx: array-like, Boolean matrix; entries with a True value
                 represent a contact
    """
    return dist_mtx <= threshold

##################
# Matrix merging #
##################

def array_union(*args):
    """
    Given an arbitrary number of binary 2-dimensional arrays, find the
    union array.
    """

    # Check that all arrays are of the same dimensions
    for array in args[1:]:
        assert args[0].shape == array.shape
    # Fill the array
    union_ar = np.zeros_like(args[0][0])
    for array in args:
        for idx, val in np.ndenumerate(array):
            if val == 1:
                idx = (idx[1], idx[2])
                union_ar[idx] = 1
    return union_ar


###################################
# Map contact matrix to alignment #
###################################


def map_contact_mtx_to_alignment(contact_mtx, orig_seq_one, orig_seq_two,
                                 aln_seq_one, aln_seq_two):
    """
    Function to map a contact matrix derived from PDB files to aligned protein
    sequences. Put another way, given a contact matrix of dimensions
    corresponding to the lengths of the proteins in the PDB file, this
    function returns a contact matrix adjusted to the size of the same proteins
    once they have been aligned, with the contacts in the corresponding
    positions.

    Arguments
    ---------
    contact_mtx:    array-like, contact matrix derived from the PDB file, of
                    dimensions (orig_seq_one, orig_seq_two)
    orig_seq_one:   str, protein sequence extracted from the PDB file
    orig_seq_two:   str, protein sequence extracted from the PDB file
    aln_seq_one:    str, protein sequence orig_seq_one once it's been aligned
    aln_seq_two:    str, protein sequence orig_seq_two once it's been aligned

    Return
    ------
    mapped_contact_mtx: array-like, contact matrix of dimensions
                        (aln_seq_one, aln_seq_two)
    """
    # Map sequences extracted from PDB file to aligned sequences
    seq_map_one = map_seq_to_alignment(orig_seq_one, aln_seq_one)
    seq_map_two = map_seq_to_alignment(orig_seq_two, aln_seq_two)
    # Get from the contact matrix which residues are in contact
    # for both sequences
    mapped_coords = map_contacts_to_alignment(
        contact_mtx, seq_map_one, seq_map_two)
    # Build a contact matrix adjusted to the dimensions of the aligned proteins
    mapped_contact_mtx = reconstruct_contact_mtx(
        mapped_coords, aln_seq_one, aln_seq_two)

    return mapped_contact_mtx


def map_seq_to_alignment(orig_seq, aln_seq):
    """
    Function to identify which positions in a sequence correspond to positions
    in the same sequence once it is aligned.

    Arguments
    ---------
    orig_seq: str, unaligned sequence
    aln_seq:  str, aligned sequence

    Returns
    -------
    seq_map: OrderedDict, in the format
             {(index in original sequence, index in aligned sequence)}

    """

    # Check that the two sequences are indeed the same,
    # and that no modifications have taken place during alignment
    # (e.g. MAFFT might delete positions when adding sequences to an alignment
    # under certain parameters)
    assert orig_seq == aln_seq.replace('-','')

    seq_map = OrderedDict()
    # Search each residue in the original sequence in the aligned sequence
    for i, res_one in enumerate(orig_seq):
        for j, res_two in enumerate(aln_seq):
            # They have to be the same residue and not have been mapped before
            if res_one == res_two and j not in seq_map.values():
                seq_map[i] = j
                break
    keys = seq_map.keys()
    vals = seq_map.values()

    # Check that the list of keys is sequential
    assert list(keys) == list(range(min(keys), max(keys) + 1))
    # Check that there are no duplicates in the values
    assert sorted(list(set(vals))) == list(vals)
    return seq_map


def get_contacts_from_matrix(contact_mtx):
    """
    Retrieve from the contact matrix the indexes at which there are contacts
    for the two sequences.

    Arguments
    ---------
    contact_mtx: array-like, contact matrix

    Returns
    -------
    contact_coords: list, contains tuples of indexes (i,j) indicating the
                    positions in the contact matrix where there are contacts
    """
    contact_coords = []
    for idx, val in np.ndenumerate(contact_mtx):
        if val == 1:
            contact_coords.append(idx)

    return contact_coords


def map_contacts_to_alignment(contact_mtx, seq_map_one, seq_map_two):
    """
    Given a contact matrix, find what positions in the aligned sequences
    correspond to the contacts.

    Arguments
    ---------
    contact_mtx: array-like, contact matrix derived from the PDB file, of
                 dimensions (orig_seq_one, orig_seq_two)
    seq_map:     OrderedDict, in the format
                 {(index in original sequence, index in aligned sequence)}

    Returns
    -------
    mapped_coords: list of tuples, where each tuple (i,j) indicates a pair
                   of positions that make contact
    """
    contact_coords = get_contacts_from_matrix(contact_mtx)
    mapped_coords = []
    for contact in contact_coords:
        try:
            pos_one = seq_map_one[contact[0]]
            pos_two = seq_map_two[contact[1]]
            mapped_coords.append((pos_one, pos_two))
        except KeyError:
            warnings.warn(RuntimeWarning, f"""Contact {contact} not present in alignments:
            	possible if dealing with sequence fragments. Check your input.""")
    return mapped_coords


def reconstruct_contact_mtx(mapped_coords, aln_seq_one, aln_seq_two):
    """
    Given indexes for the contacts corresponding to positions in the multiple
    sequence alignments, build a contact matrix of appropiate dimensions.

    Arguments
    ---------
    mapped_coords: list of tuples, where each tuple (i,j) indicates a pair
                   of positions that make contact
    aln_seq_one:    str, protein sequence orig_seq_one once it's been aligned
    aln_seq_two:    str, protein sequence orig_seq_two once it's been aligned

    Returns
    -------
    mapped_contact_mtx: array-like, contact matrix of dimensions
                        (aln_seq_one, aln_seq_two)
    """
    mapped_contact_mtx = np.zeros((len(aln_seq_one), len(aln_seq_two)))
    for contact in mapped_coords:
        mapped_contact_mtx[contact] = 1
    return mapped_contact_mtx
