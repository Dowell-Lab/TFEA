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
import shutil
import subprocess
from pathlib import Path

import numpy as np
from scipy import stats


from TFEA.simulate import pull_sequences
from TFEA.simulate import motif_insert

#Run function
#==============================================================================
def run():
    parser = parse_arguments()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    parser = parser.parse_args()
    output = Path(parser.output)

    #Get random sequences used to insert motifs from a fasta file
    if parser.random_fasta:
        sequence_pull = pull_sequences.get_sequence_indices(fasta=parser.random_fasta, 
                                                            sequence_n=parser.sequence_n, 
                                                            seed=parser.seed)
        sequences = pull_sequences.pull_sequences(fasta=parser.random_fasta, 
                                                    sequence_pull=sequence_pull, 
                                                    sequence_n=parser.sequence_n)
    else:
        print("This part still under construction")
        sys.exit()

    #Determine which motifs will be inserted into fasta file
    include_motifs = parser.include_motifs.split(',') if parser.include_motifs is not None else None
    exclude_motifs = parser.exclude_motifs.split(',') if parser.exclude_motifs is not None else None
    if exclude_motifs is not None and include_motifs is None:
        include_motifs = [m for m in get_motifs(parser.motifs) if m not in exclude_motifs]

    

    #For each motif insert it at desired locations based on insert_parameters
    if include_motifs is not None:
        #Create insert parameters used to determine where to insert motifs
        insert_parameters = zip(parser.distance_mu.split(','), 
                                parser.distance_sigma.split(','), 
                                parser.rank_range.split(','), 
                                parser.motif_number.split(','))
        #Loop through motifs to insert and insert the motif in the randomly pulled
        # sequences
        window_size = parser.largewindow*2
        for motif in include_motifs:
            for d_mu, d_sigma, rank_range, motif_number in insert_parameters:
                if d_sigma != 'uniform':
                    distance_pdf = list(stats.norm.pdf(np.arange(window_size+1), 
                                        loc=parser.largewindow - int(d_mu), 
                                        scale=int(d_sigma)))
                    distance_pdf_sum = sum(distance_pdf)
                    if distance_pdf_sum < 1:
                        distance_pdf = [x+((1-distance_pdf_sum)/len(distance_pdf)) for x in distance_pdf]
                elif d_sigma == 'uniform':
                    distance_pdf = [1.0/float(window_size+1) for _ in range(window_size+1)]
                rank_pdf = [0] * parser.sequence_n
                start, stop = [int(x) for x in rank_range.split('-')]
                if stop-start != 0:
                    equal_probability = 1/(stop - start)
                    rank_pdf[start:stop] = [equal_probability for _ in range(stop-start)]
                    rank_pdf[start] = rank_pdf[start] + 1.0 - sum(rank_pdf[start:stop])
                    sequences = motif_insert.insert_single_motif(sequences=sequences, 
                                                    sequence_n=parser.sequence_n, 
                                                    motif_database=parser.motifs, 
                                                    motif=motif, 
                                                    rank_pdf=rank_pdf, 
                                                    distance_pdf=distance_pdf, 
                                                    motif_number=int(motif_number),
                                                    seed=parser.seed)
    write_fasta(sequences=sequences, outputpath=output) #Write output

    #Optionally, split fasta for use with MDD-scores
    if parser.mdd != False:
        split_fasta(fasta_file=output, 
                    output_file1=output.parent / (output.stem + '.mdd1.fa'), 
                    output_file2=output.parent / (output.stem + '.mdd2.fa'), 
                    percent=float(parser.mdd))


#Argument Parsing
#==============================================================================
def parse_arguments():
    '''Parse user arguments
    '''
    parser = argparse.ArgumentParser(description=("Simulate data to determine "
                                                "false positive rate"))

    parser.add_argument('--output', '-o', help=("Full path to output fasta file."))
    parser.add_argument('--motifs', '-m', help=("Full path to a .meme formatted "
                                                "databse with TF motifs"))
    parser.add_argument('--include_motifs', '-i', help=("Comma-delimited list "
                                                        "of motifs to "
                                                        "include from the .meme "
                                                        "database"))
    parser.add_argument('--exclude_motifs', '-e', help=("Comma-delimited list "
                                                        "of motifs to "
                                                        "exclude from the "
                                                        ".meme database"))
    parser.add_argument('--random_fasta', '-f', help=("Full path to a fasta file "
                                                        "containing random "
                                                        "sequences from which "
                                                        "to pull [optional]"), 
                        default=False)
    parser.add_argument('--genomefasta', help=("Genomic fasta file. Required "
                                                "if random_fasta not provided"))
    parser.add_argument('--sequence_n', '-s', help=("Number of sequences to use "
                                                    "per simulation. "
                                                    "Default: 10000"), 
                        type=int)
    parser.add_argument('--seed', help=("Seed for random state. Default: Time"), 
                        type=int)           
    parser.add_argument('--largewindow', help=("Largewindow. Default: 1500"), 
                        default=1500, type=int)
    parser.add_argument('--smallwindow', help=("Smallwindow. Default: 150"), 
                        default=150, type=int)
    parser.add_argument('--distance_mu', '-dm', help=("Average distance to add "
                        "a motif from middle. Can be multiple comma-separated "
                        "averages (Ex: '0,5') but must correspond to the number of "
                        "sigmas, rank ranges, and motif numbers."))
    parser.add_argument('--distance_sigma', '-ds', help=("Variance of the "
                        "distance to add a motif from mu. Can be multiple "
                        "comma-separated sigmas (Ex: '150,300') but must correspond "
                        "to the number of mus, rank ranges, and motif numbers."
                        "If 'uniform' specified, a uniform distribution will be"
                        " used spanning the entire region."))
    parser.add_argument('--rank_range', '-rr', help=("Rank range in which to "
                        "add a motif. Can be multiple comma-separated ranges "
                        "(Ex: '0-100,500-10000') but must correspond to the number "
                        "of mus, sigmas, and motif numbers."))
    parser.add_argument('--motif_number', '-mn', help=("Number of motifs to "
                        "insert at given mu, sigma, and range. Can be multiple "
                        "comma-separated numbers (Ex: '100,100') but must "
                        "correspond to the number of mu, sigmas, and rank "
                        "ranges. Cannot be more than the number of regions "
                        "within ranked ranges."))
    parser.add_argument('--mdd', help=("Determines percentage at which to "
                        "split fasta file for MDD analysis. Default: False"), 
                        default=False)
                        
    return parser

#==============================================================================
def get_motifs(motif_database):
    motif_list = list()
    with open(motif_database,'r') as F:
        for line in F:
            if 'MOTIF' in line:
                motif_list.append(line.strip().split()[-1])
    return motif_list

#==============================================================================
def write_fasta(sequences=None, outputpath=None):
    '''
    Writes a fasta file given a list of sequences. Each fasta region will be named arbitrarily.
    
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

#==============================================================================
def split_fasta(fasta_file=None, output_file1=None, output_file2=None, 
                    percent=None):
    line_cutoff = sum([1 for line in Path(fasta_file).read_text().split('\n')])*percent - 1
    with open(output_file1, 'w') as outfile1:
        with open(output_file2, 'w') as outfile2:
            with open(fasta_file) as F:
                i = 0
                for line in F:
                    if i < line_cutoff:
                        outfile2.write(line)
                    else:
                        outfile1.write(line)
                        
                    i += 1