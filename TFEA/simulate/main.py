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
import sys
import argparse
from pathlib import Path

from TFEA import scanner

#Run function
#==============================================================================
def run():
    parser = parse_arguments().parse_args()
    output = Path(parser.output)

    if parser.random_folder:
        random_folder = Path(parser.random_folder)
    elif parser.random_fasta:
        from TFEA.simulate import pull_sequences
        pull_sequences.run(parser.random_fasta, simulations=parser.simulations,
                            sequences=parser.sequences, output=output)
        random_folder = output
        output = output / 'simulations.txt'
    else:
        print("This part still under construction")
        sys.exit()

    for fasta_file in random_folder.glob('/*.fa'):
        scanner.main(use_config=False, fasta_file=fasta_file, scanner='fimo', 
                        md=False, largewindow=1500, smallwindow=150, 
                        fimo_background='largewindow', 
                        genomefasta=parser.genomefasta, tempdir=output, 
                        fimo_motifs=parser.motifs, singlemotif=False, 
                        fimo_thresh=parser.fimo_thresh,debug=False, mdd=False, 
                        jobid=0, cpus=parser.cpus)
        



#Argument Parsing
#==============================================================================
def parse_arguments():
    '''Parse user arguments
    '''
    parser = argparse.ArgumentParser(description=("Simulate data to determine "
                                                "false positive rate"))

    parser.add_argument('--output', '-o', help=("Full path to output file."))
    parser.add_argument('--motifs', '-m', help=("Full path to a .meme formatted "
                        "databse with TF motifs"))
    parser.add_argument('--fimo_thresh', help=("P-value threshold for "
                                    "calling FIMO motif hits. Default: 1e-6"))
    parser.add_argument('--genomefasta', help=("Genomic fasta file"))
    parser.add_argument('--simulations', '-n', help=("Number of simulations "
                        "to run"))
    parser.add_argument('--sequences', '-s', help=("Number of sequences to use "
                        "per simulation"))
    parser.add_argument('--random_fasta', '-r', help=("Full path to a fasta file "
                        "containing random sequences from which to pull [optional]"), 
                        default=False)
    parser.add_argument('--random_folder', '-f', help=("Full path to a folder "
                        "containing fasta files [optional]"), 
                        default=False)
    parser.add_argument('--cpus', help=("Number of cpus to use"))
    parser.add_argument('--sbatch', help=("Use to submit an sbatch job"), 
                        metavar="EMAIL")
    parser.add_argument('--mem', help=("Amount of memory to use (use only with "
                        "sbatch argument)"))

                        
    return parser

