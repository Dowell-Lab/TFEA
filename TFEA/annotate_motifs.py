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

#Main Function
#==============================================================================
def main():
    '''Main executable script
    '''
    parser = parse_arguments()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = vars(parser.parse_args())

    #Parse gene_names file, create dictionary that ties genes to motifs
    motif_gene = {}
    with open(inputs['gene_names']) as F:
        for line in F:
            motif, gene = line.strip('\n').split('\t')
            if gene not in motif_gene:
                motif_gene[gene] = [motif]
            else:
                motif_gene[gene].append(motif)

    #Parse annotation file, tie motifs to genomic coordinates
    motif_annotation = {}
    with open(inputs['annotation']) as F:
        for line in F:
            linelist = line.strip('\n').split('\t')
            names = linelist[3].split(';')
            for n in names:
                if n in motif_gene:
                    motifs = motif_gene[n]
                    for motif in motifs:
                        motif_annotation[motif] = ['\t'.join(linelist[:3]), linelist[5]]

    #Write results to output
    with open(inputs['output'], 'w') as outfile:
        for motif in motif_annotation:
            outfile.write('\t'.join([motif_annotation[motif][0], motif, '0', motif_annotation[motif][1]]) + '\n')

#Secondary Functions
#==============================================================================
def parse_arguments():
    '''Parse user arguments
    '''
    parser = argparse.ArgumentParser(description=("Generate motif annotation "
                                                    "files"))
    parser.add_argument('--output', '-o', help=("Full path to output file."))
    parser.add_argument('--gene_names', '-g', help=("Full path to a file that "
                        "ties motif names to gene names."))
    parser.add_argument('--annotation', '-a', help=("Full path to a bed file "
                        "that has ';' delimited gene names as 4th column."))
    return parser
            
#Independent script functionality
#==============================================================================
if __name__ == "__main__":
    main()
