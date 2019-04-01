#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''This file contains functions related to parsing user arguments
'''
#==============================================================================
__author__ = 'Jonathan D. Rubin and Rutendo F. Sigauke'
__credits__ = ['Jonathan D. Rubin', 'Rutendo F. Sigauke', 'Jacob T. Stanley',
                'Robin D. Dowell']
__maintainer__ = 'Jonathan D. Rubin'
__email__ = 'Jonathan.Rubin@colorado.edu'
__version__ = '4.0'

#Imports
#==============================================================================
import argparse

#Functions
#==============================================================================
def read_arguments():
    '''This function consolidates all argument parsing into a single place

    Parameters
    ----------
    None

    Returns
    -------
    parser : object
        a parser object from argparse built-in python library
    '''
    parser = argparse.ArgumentParser(description=("Transcription Factor Enrichment "
                                                    "Analysis (TFEA)"))

    # The main inputs into TFEA to run the complete pipeline
    inputs = parser.add_argument_group('Main Inputs', 
                                        'Inputs required for full pipeline')
    inputs.add_argument('--output', '-o', help=("Full path to output directory. "
                        "If it exists, overwrite its contents."), required=True)
    inputs.add_argument('--bed1', nargs='?', help=("Bed files associated with "
                        "condition 1"), type=file)
    inputs.add_argument('--bed2', nargs='?', help=("Bed files associated with "
                        "condition 2"), type=file)
    inputs.add_argument('--bam1', nargs='?', help=("Sorted bam files "
                        "associated with condition 1"), type=file)
    inputs.add_argument('--bam2', nargs='?', help=("Sorted bam files "
                        "associated with condition 2"), type=file)
    inputs.add_argument('--label1', help=("An informative label for "
                        "condition 1"), type=str)
    inputs.add_argument('--label2', help=("An informative label for "
                        "condition 2"), type=str)
    inputs.add_argument('--config','-c', help=("A configuration file that a "
                        "user may use instead of specifying flags. Command "
                        "line flags will overwrite options within the config "
                        "file. See examples in the config_files folder."), 
                        type=file)
    inputs.add_argument('--sbatch','-s', help=("Submits an sbatch (slurm) job. "
                        "If specified, input an e-mail address."), type=str)
    inputs.add_argument('--test','-t', action='store_true', 
                        help="Performs unit testing")
    inputs.set_defaults(config=False, sbatch=False, test=False)

    # Processed inputs if a user desires to run TFEA from a specific point
    processed_inputs = parser.add_argument_group('Processed Inputs', 
                                        'Input options for pre-processed data')
    processed_inputs.add_argument('--combined_file', help=("A single bed file "
                                            "combining regions of interest."), 
                                            type=file)
    processed_inputs.add_argument('--ranked_file', help=("A bed file containing each "
                                            "regions rank as the 4th column."), 
                                            type=file)
    processed_inputs.add_argument('--fasta_file', help=("A fasta file containing "
                                            "sequences to analyzed ranked by "
                                            "the user."), type=file)
    

    # Secondary analysis inputs
    secondary_inputs = parser.add_argument_group('Secondary Analysis Inputs', 
                                                'Input options for performing '
                                                'MD-Score and Differential '
                                                'MD-Score analysis')
    secondary_inputs.add_argument('--md', help=("Switch that controls whether "
                                    "to perform MD analysis."), type=bool)
    secondary_inputs.add_argument('--mdd', help=("Switch that controls whether "
                                    "to perform differential MD analysis."), 
                                    type=bool)
    secondary_inputs.add_argument('--md_bedfile1', help=("A bed file for MD-Score "
                                    "analysis associated with condition 1."), type=file)
    secondary_inputs.add_argument('--md_bedfile2', help=("A bed file for MD-Score "
                                    "analysis associated with condition 2."), type=file)
    secondary_inputs.add_argument('--mdd_bedfile1', help=("A bed file for "
                                    "Differential MD-Score analysis associated "
                                    "with condition 1."), type=file)
    secondary_inputs.add_argument('--mdd_bedfile2', help=("A bed file for "
                                    "Differential MD-Score analysis associated "
                                    "with condition 2."), type=file)
    secondary_inputs.add_argument('--md_fasta1', help=("A fasta file for MD-Score "
                                    "analysis associated with condition 1."), type=file)
    secondary_inputs.add_argument('--md_fasta2', help=("A fasta file for MD-Score "
                                    "analysis associated with condition 2."), type=file)
    secondary_inputs.add_argument('--mdd_fasta1', help=("A fasta file for "
                                    "Differential MD-Score analysis associated "
                                    "with condition 1."), type=file)
    secondary_inputs.add_argument('--mdd_fasta2', help=("A fasta file for "
                                    "Differential MD-Score analysis associated "
                                    "with condition 2."), type=file)
    secondary_inputs.set_defaults(md=False, mdd=False)

    # Module switches
    module_switches = parser.add_argument_group('Module Switches', 
                                            'Switches for different modules')
    module_switches.add_argument('--combine', help=("Method for combining input "
                                    "bed files"), choices=[False, 'merge all',
                                    'tfit clean','tfit remove small',
                                    'intersect/merge'])
    module_switches.add_argument('--count', help=("Swtich that determines "
                                    "whether counting is performed."), 
                                    choices=[True, False])
    module_switches.add_argument('--rank', help=("Method for ranking combined "
                                    "bed file"), choices=[False, 'deseq'])
    module_switches.add_argument('--fasta', help=("Swtich that determines "
                                    "whether converting to fasta is performed."),
                                    choices=[True, False])
    module_switches.add_argument('--scanner', help=("Method for scanning fasta "
                                    "files for motifs"), choices=[False, 'fimo', 
                                    'genome hits'])
    module_switches.add_argument('--enrichment', help=("Method for calculating "
                                    "enrichment"), choices=[False, 'auc', 
                                    'auc_bgcorrect'])
    module_switches.add_argument('--debug', help=("Print memory usage to stderr"), 
                                    choices=[True, False])
    module_switches.set_defaults(combine=True, count=True, rank='deseq', 
                                    fasta=True, scanner='fimo', 
                                    enrichment='auc', debug=False)
    
    # Fasta Options
    fasta_options = parser.add_argument_group('Fasta Options', ('Options for "
                                        "performing conversion of bed to fasta"))
    fasta_options.add_argument('--genomefasta', help=("Genomic fasta file"), 
                                    type=file)

    # Scanner Options
    scanner_options = parser.add_argument_group('Scanner Options', 
                                        'Options for performing motif scanning')
    scanner_options.add_argument('--fimo_thresh', help=("P-value threshold for "
                                    "calling FIMO motif hits"))
    scanner_options.add_argument('--fimo_motifs', help=("Full path to a .meme "
                                    "formatted motif databse file. Some "
                                    "databases included in motif_databases "
                                    "folder."))
    scanner_options.add_argument('--fimo_background', help=("Options for "
                                    "choosing mononucleotide background "
                                    "distribution to use with FIMO. "
                                    "{'largewindow', 'smallwindow', int, file}"))
    scanner_options.add_argument('--genomehits', help=("A folder containing "
                                    "bed files with pre-calculated motif hits "
                                    "to a genome. For use with 'genome hits' "
                                    "scanner option."), type=file)
    scanner_options.add_argument('--singlemotif', help=("Option to run "
                                    "analysis on a subset of motifs within "
                                    "specified motif database or genome hits. "
                                    "Can be a single motif or a comma-separated "
                                    "list of motifs."))
    scanner_options.set_defaults(fimo_thres=0.000001, 
                                    fimo_background='largewindow', 
                                    singlemotif=False)

    # Enrichment Options
    enrichment_options = parser.add_argument_group('Enrichment Options', 
                                'Options for performing enrichment analysis')
    enrichment_options.add_argument('--permutations', help=("Number of "
                                        "permutations to perfrom for "
                                        "calculating p-value."), type=int)
    enrichment_options.add_argument('--largewindow', help=("The size (bp) of a "
                                        "large window around input regions "
                                        "that captures background."), type=int)
    enrichment_options.add_argument('--smallwindow', help=("The size (bp) of a "
                                        "small window arount input regions "
                                        "that captures signal."), type=int)
    enrichment_options.set_defaults(permutations=1000, largewindow=1500, 
                                        smallwindow=150)
    
    # Output Options
    output_options = parser.add_argument_group('Output Options', 
                                                'Options for the output.')
    output_options.add_argument('--dpi', help=("Resolution of output figures."), 
                                    type=int)
    output_options.add_argument('--padjcutoff', help=("A p-adjusted cutoff "
                                    "value that determines some plotting output."))
    output_options.add_argument('--pvalcutoff', help=("A p-value cutoff "
                                    "value that determines some plotting output."))
    output_options.add_argument('--textOnly', '-t', help="Suppress html output.", 
                                    action='store_true')
    output_options.set_defaults(dpi=None, padjcutoff=0.001, pvalcutoff=0.01, 
                                textOnly=False)

    return parser

def validate_parser(parser):
    boolean = True
    return boolean