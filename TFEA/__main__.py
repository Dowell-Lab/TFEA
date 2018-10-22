#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''This is the main script that runs TFEA. It imports config, create_html, 
    independent_functions, and dependent_functions modules contained within
    this package
'''
#==============================================================================
__author__ = 'Jonathan D. Rubin and Rutendo F. Sigauke'
__credits__ = ['Jonathan D. Rubin', 'Rutendo F. Sigauke', 'Jacob T. Stanley',
                'Robin D. Dowell']
__maintainer__ = 'Jonathan D. Rubin'
__email__ = 'Jonathan.Rubin@colorado.edu'
__version__ = '3.0'
#==============================================================================
#PRIMARY IMPORTS
#==============================================================================
import os
import sys
import time
import argparse
import configparser
import preconfig_functions
#==============================================================================
#ARGUMENT PARSING
#==============================================================================
#TFEA source directory
srcdirectory = os.path.dirname(os.path.realpath(__file__))

#argparse to add arguments to this python package
parser = argparse.ArgumentParser(description='Transcription Factor Enrichment \
                                    Analysis (TFEA) takes as input a \
                                    configuration file (.ini) and outputs \
                                    a folder containing TFEA results.',
                                usage='TFEA --config CONFIG.ini [--sbatch \
                                    email@address.com]')

parser.add_argument('--config','-c',metavar='',help='REQUIRED. A \
                        configuration file containing .ini suffix \
                        (ex. config.ini). See example in the examples folder.')

parser.add_argument('--sbatch','-s',default=False,metavar='',help='OPTIONAL. \
                        Submits an sbatch job. If specified, input an e-mail \
                        address.')

#Display help message when no args are passed.
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

#If user provided arguments, then parse them
sbatch = parser.parse_args().sbatch
configfile = parser.parse_args().config
config_object = configparser.ConfigParser(
                            interpolation=configparser.ExtendedInterpolation())
config_object.read(configfile)
#==============================================================================
#CREATING DIRECTORIES
#==============================================================================
#If user specifies the --sbatch flag, then we first create the output 
#directories then run the sbatch script with the 'SUBMITTED' command submitted 
#to the --sbatch flag so we know not to remake output directories. If --sbatch 
#flag not specified, simply make output directories and continue.
if sbatch == False:
    output, tempdir, figuredir, e_and_o = preconfig_functions.make_out_directories(
                                                dirs=True, 
                                                config_object=config_object)
elif str(sbatch) == 'SUBMITTED':
    output, tempdir, figuredir, e_and_o = preconfig_functions.make_out_directories(
                                                dirs=False, 
                                                config_object=config_object)
else:
    output, tempdir, figuredir, e_and_o = preconfig_functions.make_out_directories(
                                                dirs=True, 
                                                config_object=config_object)
    scriptdir = os.path.join(os.path.dirname(srcdirectory), 'scripts')
    script = os.path.join(scriptdir, 'run_main.sbatch')
    email = str(sbatch)
    os.system("sbatch --error=" + e_and_o + "/%x.err --output=" + e_and_o 
                + "/%x.out --mail-user=" + email + " --export=cmd='" 
                + srcdirectory+" --config " + configfile 
                + " --sbatch SUBMITTED' " + script)

    sys.exit(("TFEA has been submitted using an sbatch script, use qstat to "
            "check its progress."))

#Run the config_parser script which will create variables for all folders and 
#paths to use throughout TFEA
preconfig_functions.parse_config(srcdirectory=srcdirectory, 
                                    config_object=config_object,
                                    output=output,tempdir=tempdir,
                                    figuredir=figuredir)

#Verify config file to make sure user has inputted all necessary variables
# config_dict = independent_functions.verify_config_object(config=config)
preconfig_functions.verify_config_file()
#==============================================================================
#SECONDARY IMPORTS
#==============================================================================
from multiprocessing import Pool
import multiprocessing as mp
import config
import independent_functions
import dependent_functions
#==============================================================================
#MAIN SCRIPT
#==============================================================================
#This module takes the input list of BED files, concatenates them, and then 
#merges them via bedtools.
COMBINEtime = time.time()
if config.COMBINE:
    # bedfile = independent_functions.merge_bed(beds=config.BED1+config.BED2, 
    #                                             tempdir=tempdir)
    # bedfile = independent_functions.tfit_clean_merge(
    #                           beds=config.BED1+config.BED2, tempdir=tempdir)
    bedfile = independent_functions.intersect_merge_bed(bed1=config.BED1,
                                                        bed2=config.BED2, 
                                                        tempdir=tempdir)
else:
    bedfile = config.BED1 + config.BED2
    bedfile = bedfile[0]
COMBINEtime = time.time()-COMBINEtime 

#This module counts reads from all Bam files in BAM1 and BAM2 and creates 
#count_file with this info.
COUNTtime = time.time()
if config.COUNT:
    print "Counting reads in regions..."
    count_file = independent_functions.count_reads(bedfile=bedfile, 
                                                    bam1=config.BAM1, 
                                                    bam2=config.BAM2, 
                                                    tempdir=tempdir, 
                                                    label1=config.LABEL1, 
                                                    label2=config.LABEL2)
    print "done"
elif config.DESEQ:
    #If you don't want to perform multibamcov but still want to perform
    # DE-Seq, user must provide a count_file
    count_file = config.COUNT_FILE
COUNTtime = time.time()-COUNTtime

#This module runs DESeq on the count_file produced above - this is used to rank
#inputted regions
DESEQtime = time.time()
if config.DESEQ:
    print "Running DESeq..."
    deseq_file = dependent_functions.deseq_run(bam1=config.BAM1, 
                                                bam2=config.BAM2, 
                                                tempdir=tempdir, 
                                                count_file=count_file,
                                                label1=config.LABEL1, 
                                                label2=config.LABEL2,
                                                figuredir=figuredir)

    ranked_center_file = independent_functions.rank_deseqfile(
                                                deseq_file=deseq_file,
                                                tempdir=tempdir)
    print "done"
elif config.CALCULATE:
    #If user does not want to perform DE-Seq but still wants TFEA to caluclate
    # TF enrichment, user must provide a ranked_center_file
    ranked_center_file = config.RANKED_CENTER_FILE
DESEQtime = time.time()-DESEQtime


#This is the bulk of the analysis of this package, it performs:
#   1. GC distribution across regions for plotting
#   2. Millions mapped calculation for meta eRNA plots
#   3. Motif distance calculation to the center of inputted regions
#   4. Enrichment score calculation via AUC method
#   5. Random shuffle simulation and recalculation of enrichment score
#   6. Plotting and generation of html report
CALCULATEtime = time.time()
if config.CALCULATE:
    dependent_functions.calculate(tempdir=tempdir, 
                                    outputdir=output, 
                                    bam1=config.BAM1, 
                                    bam2=config.BAM2, 
                                    label1=config.LABEL1,
                                    label2=config.LABEL2,
                                    motif_hits=config.MOTIF_HITS, 
                                    singlemotif=config.SINGLEMOTIF, 
                                    ranked_center_file=ranked_center_file, 
                                    COMBINEtime=COMBINEtime, 
                                    COUNTtime=COUNTtime, DESEQtime=DESEQtime, 
                                    CALCULATEtime=CALCULATEtime, 
                                    fimo=config.FIMO, 
                                    plot=config.PLOT, 
                                    pool=config.POOL,
                                    padj_cutoff=config.PADJCUTOFF, 
                                    logos=config.LOGOS, 
                                    figuredir=figuredir,
                                    largewindow=config.LARGEWINDOW, 
                                    smallwindow=config.SMALLWINDOW,
                                    motifdatabase=config.MOTIFDATABASE,
                                    genomefasta=config.GENOMEFASTA)


#Here we simply remove large bed files that are produced within this package. 
#This option can be turned on/off within the config file.
if not config.TEMP:
    print "Removing temporary bed files..."
    files_to_keep = ['count_file.bed', 'DESeq.R', 'DESeq.res.txt', 
                    'DESeq.Rout', 'ranked_file.bed', 'markov_background.txt', 
                    'ranked_file.center.bed', 'ranked_file.center.sorted.bed']
    for file1 in os.listdir(tempdir):
        if file1 not in files_to_keep:
             os.remove(os.path.join(tempdir, file1))

print "done"
