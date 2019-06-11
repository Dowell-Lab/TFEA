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

from TFEA import scanner
from TFEA import enrichment
from TFEA import exceptions

#Run function
#==============================================================================
def run():
    parser = parse_arguments()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    parser = parser.parse_args()
    output = Path(parser.output)
    parent_folder = output.parent
    if parser.sbatch:
        script = Path(__file__).absolute().parent / 'simulate.sbatch'
        args = sys.argv
        args[args.index('--sbatch')+1] = 'False'
        try:
            subprocess.run(["sbatch", 
                                "--error=" + (parent_folder / "%x.err").as_posix(), 
                                "--output=" + (parent_folder / "%x.out").as_posix(), 
                                "--mail-user=" + parser.sbatch, 
                                "--export=cmd=" + ' '.join(args), 
                                "--ntasks=" + str(parser.cpus),
                                "--mem=" + str(parser.mem),
                                script], stderr=subprocess.PIPE, 
                                stdout=subprocess.PIPE, check=True)
        except subprocess.CalledProcessError as e:
            raise exceptions.SubprocessError(e.stderr.decode())

    if parser.random_folder:
        random_folder = Path(parser.random_folder)
    elif parser.random_fasta:
        from TFEA.simulate import pull_sequences
        pull_sequences.run(parser.random_fasta, simulations=parser.simulations,
                            sequences=parser.sequences, output=output)
        random_folder = parent_folder
    else:
        print("This part still under construction")
        sys.exit()

    motif_counts = {}
    simulations = 0
    for fasta_file in random_folder.glob('/*.fa'):
        motif_distances, _, _, _, _ = scanner.main(use_config=False, 
                        fasta_file=fasta_file, scanner='fimo', 
                        md=False, largewindow=parser.largewindow, 
                        smallwindow=parser.smallwindow, 
                        fimo_background='largewindow', 
                        genomefasta=parser.genomefasta, tempdir=output, 
                        fimo_motifs=parser.motifs, singlemotif=False, 
                        fimo_thresh=parser.fimo_thresh, debug=False, mdd=False, 
                        jobid=0, cpus=parser.cpus)
        results, _, _ = enrichment.main(use_config=False, 
                        motif_distances=motif_distances, 
                        enrichment='auc', permutations=1000, debug=False, 
                        largewindow=parser.largewindow, 
                        smallwindow=parser.smallwindow, 
                        md=False, mdd=False, 
                        cpus=parser.cpus, jobid=0, suppress_plots=True)
        for result in results:
            motif = result[0]
            if motif not in motif_counts:
                motif_counts[motif] = 0
            else:
                motif_counts[motif] += 1
        
        simulations += 1
    
    with open(output, 'w') as outfile:
        user_input = vars(parser)
        for key in user_input:
            outfile.write('#' + key + '=' + str(user_input[key]) + '\n')
        outfile.write('#Simulations=' + str(simulations) + '\n')
        outfile.write('#Motif\tFalsePositiveNumber\n')
        for motif in motif_counts:
            outfile.write(motif + '\t' + motif_counts[motif] + '\n')
    
    if not parser.keep_fasta:
        for fasta_file in random_folder.glob('/*.fa'):
            fasta_file.unlink()


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
    parser.add_argument('--padjcutoff', help=("P-Adjusted cutoff to be "
                        "considered as signficant. Default: 0.01"), default=0.01)
    parser.add_argument('--random_folder', '-f', help=("Full path to a folder "
                        "containing fasta files [optional]"), 
                        default=False)
    parser.add_argument('--keep_fasta', help=("Turn on to keep fasta files"), 
                        action='store_true')
    parser.add_argument('--cpus', help=("Number of cpus to use"))
    parser.add_argument('--sbatch', help=("Use to submit an sbatch job"), 
                        metavar="EMAIL", default=False)
    parser.add_argument('--mem', help=("Amount of memory to use (use only with "
                        "sbatch argument)"))
    parser.add_argument('--largewindow', help=("Largewindow. Default: 1500"), 
                        default=1500)
    parser.add_argument('--smallwindow', help=("Smallwindow. Default: 150"), 
                        default=150)
                        
    return parser

