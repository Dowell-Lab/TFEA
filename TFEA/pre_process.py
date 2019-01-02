#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''This file contains functions related to pre-processing input into
    a format that can be used with enrichment analysis functions
'''
#==============================================================================
__author__ = 'Jonathan D. Rubin and Rutendo F. Sigauke'
__credits__ = ['Jonathan D. Rubin', 'Rutendo F. Sigauke', 'Jacob T. Stanley',
                'Robin D. Dowell']
__maintainer__ = 'Jonathan D. Rubin'
__email__ = 'Jonathan.Rubin@colorado.edu'
__version__ = '4.0'
#==============================================================================
#Imports
#==============================================================================
import os
import pickle
import subprocess
#==============================================================================
#Functions
#==============================================================================
def run(input_file=None, input_type=None, scanner=None, background=None, 
        output=None, motifs=None):
    if input_type == 'bed+bam':
        print("Under construction")
    if input_type == 'counts':
        print("Under construction")
    if input_type == 'fasta':
        ranked_list = fasta_input(fasta=input_file, scanner=scanner, 
                                    background=background, 
                                    motifs=motifs)
                                    
        return ranked_list
    if input_type == 'ranked_list':
        print("Under construction")

    
#==============================================================================
#Functions shared between input types
#==============================================================================
def fimo(motif=None, bg_file=None, fasta_file=None, motifdatabase=None,
        thresh=0.0001):
    '''This function runs fimo on a given fastafile for a single motif in a 
        provided motif database. The output is cut and sorted to convert into 
        a sorted bed file

        Parameters
        ----------
        tempdir : string
            full path to temp directory in output directory (created by TFEA)

        motifdatabase : string
            full path to a motif database file in meme format

        bgfile : string
            full path to a markov background model

        motif : string
            the name of a motif that matches a motif within motifdatabase

        fastafile : string
            full path to a fasta file that fimo will perform motif scanning on
            
        Returns
        -------
        fimo_out : string
            full path to where fimo output which is stored within the tempdir 
            directory.
    '''
    command = ["fimo",  "--skip-matched-sequence",  "--text", 
                "--verbosity", "1", 
                "--thresh", str(thresh),
                "--bgfile", bg_file, 
                "--motif", motif,
                motifdatabase,
                fasta_file]
    fimo_out = subprocess.check_output(command).decode('utf-8')

    return fimo_out

def make_bg_file(fastafile=None, order='2', outputpath=None):
    '''This function runs meme's fasta-get-markov function that generates a 
        background markov file (for use with fimo) from a fasta file.

        Parameters
        ----------
        tempdir : string
            full path to temp directory in output directory (created by TFEA)

        fastafile : string
            full path to fasta file that will be used to generate the markov
            background model file

        order : string
            an integer formatted as a string where a user may specify what 
            order markov model they would like (default='0')
            
        Returns
        -------
        None
    '''
    markov_background = outputpath
    with open(markov_background) as outfile:
        command = ["fasta-get-markov", "-m", order, fastafile]
        subprocess.call(command, stdout=outfile)

    return markov_background
#==============================================================================
#Functions associated with input being of type 'fasta'
#==============================================================================
def fasta_input(fasta=None, scanner=None, background=None, motifs=None):
    motifs = motif_function(motifs)

    fimo_out = fimo(motif=motif, bg_file=bg_file, ranked_fasta_file=fasta, 
                    tempdir=None, motifdatabase=None, thresh=0.0001)
    
    ranked_list = other_function(fimo_out)

    return ranked_list
#==============================================================================

#==============================================================================

#==============================================================================

#==============================================================================
def test(a, arg1=None):
    print(arg1, a)

if __name__ == "__main__":
    import multiprocessing
    from functools import partial

    pool = multiprocessing.Pool(3)
    a = range(3)
    pool.map(partial(test, arg1=1), a)
