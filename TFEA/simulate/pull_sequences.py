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

#Main Function
#==============================================================================
def run(outputpath, fasta_file=None, sequence_n=None):
    sequence_pull = get_sequence_indices(fasta=fasta_file, 
                                            sequence_n=sequence_n)

    sequences = pull_sequences(fasta=fasta_file, 
                                sequence_pull=sequence_pull, 
                                sequence_n=sequence_n)

    write_fasta(sequences=sequences, outputpath=outputpath)

#Auxiliary Functions
#==============================================================================
def count_fasta(fasta):
    '''
    This function counts how many regions are in a fasta file
    
    Parameters
    ----------
    fasta : string
        full path to a fasta (.fa) file
        
    Returns
    -------
    count : int
        number of regions in the fasta file
    '''
    
    with open(fasta, 'r') as F:
        count = sum([1 for line in F if '>' in line])
        
    return count

#==============================================================================
def get_sequence_indices(fasta=None, sequence_n=None, seed=None):
    '''
    Get the indices of sequences to pull from a fasta file
    
    Parameters
    ----------
    fasta : string
        full path to a fasta file containing random sequences
        
    sequence_n : int
        the total number of sequences that motifs will be inserted into
        
    Returns
    -------
    sequence_pull : list or array
        a list of indices corresponding to sequences to pull from the inputted fasta file
    '''
    total_sequences = count_fasta(fasta)
    if seed is not None:
        sequence_pull = np.random.choice(total_sequences, size=int(sequence_n), 
                                            replace=False)
    else:
        np.random.seed(seed)
        sequence_pull = np.random.choice(total_sequences, size=int(sequence_n), 
                                            replace=False)
    
    return sequence_pull

#==============================================================================
def pull_sequences(fasta=None, sequence_pull=None, sequence_n=None):
    '''
    Get sequences from a fasta file
    
    Parameters
    ----------
    fasta : string
        full path to a fasta file containing random sequences
        
    sequence_pull : list or array
        a list of indices corresponding to which sequences to pull from fasta file
        
    Returns
    -------
    sequences : list or array
        a list of sequences pulled from the fasta file
    '''
    sequences = [''] * sequence_n
    counter = 0
    sequence_index = -1
    pull_sequence = False
    with open(fasta,'r') as F:
        for line in F:
            if '>' in line:
                pull_sequence = False
                if counter in sequence_pull:
                    pull_sequence = True
                    sequence_index += 1
                counter += 1
            elif pull_sequence:
                sequences[sequence_index] += line.strip('\n')
                
    return sequences

#==============================================================================
def write_fasta(sequences=None, outputpath=None):
    '''
    Writes a fasta file given a list of sequences. Each fasta region will be 
        named arbitrarily.
    
    Parameters
    ----------
    sequences : list or array
        a list of strings containing sequences to write into a fasta file
        
    outputpath : str
        a string pointing to a full path to an output file (inlcudes filename)
        
    Returns
    -------
    None
    '''
    sequence_counter = 1
    with open(outputpath,'w') as outfile:
        for sequence in sequences:
            outfile.write('>' + str(sequence_counter) + '\n')
            outfile.write(sequence + '\n')
            sequence_counter += 1