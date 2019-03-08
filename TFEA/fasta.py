#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''This module takes as input a bed file (must be ranked if downstream TFEA 
    analysis desired) and returns a fasta file.
'''

#==============================================================================
__author__ = 'Jonathan D. Rubin and Rutendo F. Sigauke'
__credits__ = ['Jonathan D. Rubin', 'Rutendo F. Sigauke', 'Jacob T. Stanley',
                'Robin D. Dowell']
__maintainer__ = 'Jonathan D. Rubin'
__email__ = 'Jonathan.Rubin@colorado.edu'

#Imports
#==============================================================================
import os
import subprocess
from pathlib import Path

import config
from exceptions import FileEmptyError, SubprocessError

#Main Script
#==============================================================================
def main(ranked_file=None, md_bedfile1=None, md_bedfile2=None, md=config.MD, 
            genomefasta=config.GENOMEFASTA, tempdir=Path(config.TEMPDIR)):
    '''The main script of the FASTA module. Takes as input a bed file (must
        be ranked if downstream TFEA desired) and converts it to a fasta file.
        Can take input from the RANK module.

    Parameters
    ----------
    ranked_file : str
        Full path to a ranked bed file
    md_bedfile1 : str
        Full path to a bed file corresponding to a single condition. Only 
        required if md score analysis desired
    md_bedfile2 : str
        Full path to a bed file corresponding to a single condition. Only 
        required if md score analysis desired
    md : boolean
        Whether md score analysis is desired. If True, requires bed files for
        each condition. These can be generated in the COMBINE module.
    genomefasta : str
        Full path to a fasta file for desired genome
    tempdir : str
        Full path to a directory where files will be saved

    Returns
    -------
    fasta_file : str
        Full path to the output fasta file
    md_fasta1 : str
        Full path to the output fasta file for use with md score analysis 
        corresponding to a single condition
    md_fasta2 : str
        Full path to the output fasta file for use with md score analysis 
        corresponding to a single condition

    Raises
    ------
    FileEmptyError
        If any resulting file is empty
    '''
    fasta_file = getfasta(bedfile=ranked_file, genomefasta=genomefasta, 
                            tempdir=tempdir, outname='fasta_file.fa')
    if os.stat(fasta_file).st_size == 0:
        raise FileEmptyError("Error in FASTA module. Resulting fasta file is empty.")

    if md:
        md_fasta1 = getfasta(bedfile=md_bedfile1, genomefasta=genomefasta, 
                                tempdir=tempdir, outname='md_bedfile1.fa')
        md_fasta2 = getfasta(bedfile=md_bedfile2, genomefasta=genomefasta, 
                                tempdir=tempdir, outname='md_bedfile2.fa')
        if os.stat(md_fasta1).st_size == 0 or os.stat(md_fasta2).st_size == 0:
            raise FileEmptyError("Error in FASTA module. Resulting md fasta file is empty.")
        
        return fasta_file, md_fasta1, md_fasta2
    
    return fasta_file

#Functions
#==============================================================================
def getfasta(bedfile=None, genomefasta=None, tempdir=None, outname=None):
    '''Converts a bed file to a fasta file using bedtools. Outputs into the 
        tempdir directory created by TFEA.

    Parameters
    ----------
    bedfile : string
        full path to a bed file

    genomefasta : string
        full path to a fasta file for the genome of interest
        
    tempdir : string
        full path to temp directory in output directory (created by TFEA)

    Returns
    -------
    ranked_file_fasta : string 
        full path to a fasta file containing the inputted bed file regions in 
        fasta format 
    '''
    fasta_file = tempdir / outname
    getfasta_command = ["bedtools", "getfasta", "-name",
                        "-fi", genomefasta, 
                        "-bed", bedfile,
                        "-fo", fasta_file]
    
    try:
        subprocess.run(getfasta_command, check=True, 
                        stdout=subprocess.PIPE, 
                        stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        raise SubprocessError(e.stderr.decode())

    return fasta_file

#==============================================================================
def fasta_linecount(fastafile=None):
    linecount = 0
    with open(fastafile) as F:
        for line in F:
            if line[0] == '>':
                linecount += 1
    
    return linecount
            
#==============================================================================