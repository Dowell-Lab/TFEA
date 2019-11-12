#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''This file contains scripts to generate motif annotation files from a motif 
    names.tsv file and a bed annotation file that has ';' delimited gene names
    as 4th column
'''
#==============================================================================
__author__ = 'Jonathan D. Rubin and Rutendo F. Sigauke'
__credits__ = ['Jonathan D. Rubin', 'Rutendo F. Sigauke', 'Jacob T. Stanley',
                'Robin D. Dowell']
__maintainer__ = 'Jonathan D. Rubin'
__email__ = 'Jonathan.Rubin@colorado.edu'

#Imports
#==============================================================================
import numpy as np

#==============================================================================
def insert_single_motif(sequences=None, sequence_n=None, motif_database=None, 
                        motif=None, rank_pdf=None, rank_inserts=None, 
                        distance_pdf=None, distance_inserts=None, 
                        motif_number=None, seed=None):

    #Get PSSM of motif from motif_database file
    PSSM = get_PSSM(motif_database=motif_database, motif=motif)
            
    #Get ranks and positions to insert motif
    if rank_inserts is None:
        # rank_inserts = get_rank_inserts(sequence_n=sequence_n, 
        #                                 rank_pdf=rank_pdf, seed=seed)
        rank_inserts = get_fixednumber_rank_inserts(sequence_n=sequence_n, 
                                                    rank_pdf=rank_pdf, 
                                                    motif_number=motif_number, 
                                                    seed=seed)
    if distance_inserts is None:
        distance_inserts = get_distance_inserts(sequence_n=sequence_n, 
                                                distance_pdf=distance_pdf, 
                                                rank_insert_size=len(rank_inserts), 
                                                seed=seed)
    
    #Insert motifs at given ranks with given distances
    sequences = insert_motif(PSSM=PSSM, rank_inserts=rank_inserts, 
                            distance_inserts=distance_inserts, 
                            sequences=sequences)

    return sequences

#==============================================================================
def get_PSSM(motif_database=None, motif=None):
    '''
    Obtain a pssm model from a meme formatted database file. Warning: If there are multiple
    motif matches, this function will return the last match in the database.
    
    Parameters
    ----------
    motif_database : string
        full path to a MEME formatted (.meme) file
    motif : string
        the name of a motif that exactly matches a motif in motif_database
        
    Returns
    -------
    PSSM : list or array
        a list of lists where each corresponding to position then alphabet probability
    '''
    motif_hit = False
    PSSM = list()
    with open(motif_database,'r') as F:
        for line in F:
            if 'MOTIF' in line:
                if motif in line:
                    motif_hit = True
                else:
                    motif_hit = False
            elif motif_hit and 'URL' not in line and 'letter-probability' not in line and line != '\n':
                acgt_probabilities = [float(x) for x in line.strip('\n').split()]
                total_prob = sum(acgt_probabilities)
                acgt_probabilities = [x/total_prob for x in acgt_probabilities] #Convert to probabilities
                PSSM.append(acgt_probabilities)
            
    return PSSM

#==============================================================================
def get_rank_inserts(sequence_n=None, rank_pdf=None, seed=None):
    if seed is None:
        rank_inserts = list(np.random.choice(sequence_n, size=(sequence_n), 
                        replace=True, p=rank_pdf))
    else:
        np.random.seed(seed)
        rank_inserts = list(np.random.choice(sequence_n, size=(sequence_n), 
                        replace=True, p=rank_pdf))

    return rank_inserts

#==============================================================================
def get_distance_inserts(sequence_n=None, distance_pdf=None, 
                            rank_insert_size=None, seed=None):
    if seed is None:
        distance_inserts = list(np.random.choice(len(distance_pdf), 
                            size=(rank_insert_size), p=distance_pdf))
    else:
        np.random.seed(seed)
        distance_inserts = list(np.random.choice(len(distance_pdf), 
                            size=(rank_insert_size), p=distance_pdf))

    return distance_inserts

#==============================================================================
def get_fixednumber_rank_inserts(sequence_n=None, rank_pdf=None, 
                                    motif_number=None, seed=None):
    if seed is None:
        rank_inserts = list(np.random.choice(sequence_n, size=(motif_number), 
                            replace=False, p=rank_pdf))
    else:
        np.random.seed(seed)
        rank_inserts = list(np.random.choice(sequence_n, size=(motif_number), 
                            replace=False, p=rank_pdf))
    
    return rank_inserts

#==============================================================================
def insert_motif(PSSM=None, alphabet=['A','C','G','T'], rank_inserts=None, 
                    distance_inserts=None, sequences=None):
    '''
    This function inserts a motif PSSM in a list of sequences based on 
    rank_inserts and distance_inserts
    
    Parameters
    ----------
    PSSM : list or array
        a list of lists containing a motif PSSM
        
    alphabet : list or array
        the alphabet corresponding to the PSSM
        
    rank_inserts : list or array
        a list of indices determining which regions will get an inserted motif
        
    distance_inserts : list or array
        a list of distances that corresponds to each rank_inserts region
        
    sequences : list or array
        a list of sequences
        
    Returns
    -------
    sequences : list or array
        a list of sequences modified to contain the input PSSM at the locations 
        specified
    '''
    motif_length = len(PSSM)
    modified_sequences = [''] * len(sequences)
    for i in range(len(sequences)):
        sequence_to_modify = sequences[i]
        if i in rank_inserts:
            distance = distance_inserts[rank_inserts.index(i)]
            sequence_length = len(sequence_to_modify)
            start = distance-int(motif_length/2)
            if start < 0:
                start = 0
            stop = start+motif_length
            if stop > sequence_length:
                start = sequence_length-motif_length
                stop = sequence_length

            
            # print(f"Motif inserted: {start}-{stop}")
            # print(f'Previous Sequence {sequence_to_modify[start:stop]}')

            for j in range(start,stop):
                sequence_to_modify = (sequence_to_modify[:j] 
                                        + alphabet[np.random.choice(len(alphabet), p=PSSM[j-start])] #Pull from PSSM
#                                         + alphabet[PSSM[j-start].index(max(PSSM[j-start]))] #Consensus
                                        + sequence_to_modify[j + 1:])
            
            # print(f'Modified Sequence {sequence_to_modify[start:stop]}')
        
        modified_sequences[i] = sequence_to_modify
        
    return modified_sequences