#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''This file contains a list of functions associated with counting reads that 
    map to inputted regions
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
from exceptions import FileEmptyError

#Main Script
#==============================================================================
def main(bedfile=None, bam1=config.vars.BAM1, bam2=config.vars.BAM2, 
            tempdir=config.vars.TEMPDIR, label1=config.vars.LABEL1, 
            label2=config.vars.LABEL2):
    '''This is the main script of the COUNT module which takes as input a
        bedfile and bam files and counts reads over the bed regions

    Parameters
    ----------
    bedfile : str
        Full path to a bedfile. Can be generated with the COMBINE module or
        defined by a user
    bam1 : list
        A list of strings specifying full paths to bam files associated with 
        a single condition (replicates)
    bam2 : list
        A list of strings specifying full paths to bam files associated with 
        a single condition (replicates)
    tempdir : str
        Full path to a directory where files will be saved
    label1 : str
        An informative label for the condition corresponding to files within
        bam1. Used to add a header line to the count file for downstream
        DE-Seq analysis.
    label2 : str
        An informative label for the condition corresponding to files within
        bam2. Used to add a header line to the count file for downstream
        DE-Seq analysis.

    Returns
    -------
    count_file : str
        Full path to a count file that contains the counts over regions within
        inputted bed files for all inputted bam files

    Raises
    ------
    FileEmptyError
        If any resulting file is empty
    '''
    count_file = count_reads(bedfile=bedfile, bam1=bam1, bam2=bam2, 
                            tempdir=tempdir, label1=label1, label2=label2)
    sample_number = (len(bam1)+len(bam2))
    millions_mapped = sum_reads(count_file=count_file, 
                                sample_number=sample_number)
    
    if os.stat(count_file).st_size == 0:
        raise FileEmptyError("Error in COUNT module. Resulting file is empty.")
        
    return count_file

#Functions
#==============================================================================
def count_reads(bedfile=None, bam1=None, bam2=None, tempdir=None, label1=None, 
                label2=None):
    '''Counts reads across regions in a given bed file using bam files inputted
        by a user

    Parameters
    ----------
    bedfile : string
        full path to a bed file containing full regions of interest which will 
        be counted using bedtools multicov

    bam1 : list or array
        a list of full paths to bam files pertaining to a single condition 
        (i.e. replicates of a single treatment)

    bam2 : list or array
        a list of full paths to bam files pertaining to a single condition 
        (i.e. replicates of a single treatment)

    tempdir : string
        full path to temp directory in output directory (created by TFEA)

    label1 : string
        the name of the treatment or condition corresponding to bam1 list

    label2 : string
        the name of the treatment or condition corresponding to bam2 list

    Returns
    -------
    None
    '''
    #This os.system call runs bedtools multicov to count reads in all specified
    #BAMs for given regions in BED
    count_file = tempdir / "count_file.bed"
    multicov_command = ["bedtools", "multicov", 
                        "-bams"] + bam1+bam2 + ["-bed", bedfile]
    with open(count_file, 'w') as outfile:
        subprocess.run(multicov_command, stdout=outfile, check=True)

    #This section adds a header to the count_file and reformats it to remove 
    #excess information and add a column with the region for later use
    count_file_header = tempdir / "count_file.header.bed"
    outfile = open(count_file_header, 'w')
    outfile.write("chrom\tstart\tstop\tregion\t" 
                    + '\t'.join([label1]*len(bam1)) + "\t" 
                    + '\t'.join([label2]*len(bam2)) + "\n")

    with open(count_file) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            chrom,start,stop = line[:3]
            counts = line[-(len(bam1)+len(bam2)):]
            outfile.write('\t'.join([chrom,start,stop]) + "\t" 
                            + chrom + ":" + start + "-" + stop + "\t"
                            + '\t'.join(counts) + "\n")
    outfile.close()
    return count_file_header

#==============================================================================
def sum_reads(count_file=None, sample_number=None):
    '''This function calculates millions mapped reads to regions of interest
        to be used later for normalization of meta plot
    
    Parameters
    ----------
    count_file : string
        The full path to a count file containing reads mapping to regions of
        interest
    
    sample_number : int
        The total number of samples

    Returns
    -------
    millions_mapped : list
        A list of ints that corresponds to the number of millions mapped per
        sample in the order that appears in the count_file (by column)
    '''
    millions_mapped = [0.0]*sample_number
    with open(count_file) as F:
        F.readline()
        for line in F:
            line = line.strip('\n').split('\t')
            [x+float(y) for x,y in zip(millions_mapped,line[-sample_number:])]

    return millions_mapped