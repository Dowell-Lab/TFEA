#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''This is the main script that runs TFEA. It imports config, create_html, 
    independent_functions, and dependent_functions modules contained within
    this package
'''
#==============================================================================
__author__ = 'Jonathan D. Rubin and Rutendo Sigauke'
__credits__ = ['Jonathan D. Rubin', 'Rutendo Sigauke', 'Jacob Stanley',
                'Robin Dowell']
__maintainer__ = 'Jonathan D. Rubin'
__email__ = 'Jonathan.Rubin@colorado.edu'
__version__ = '3.0'
#==============================================================================
#MAIN IMPORTS
#==============================================================================
import os
import sys
import time
import argparse
import configparser
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
#==============================================================================
#SECONDARY IMPORTS
#==============================================================================
from multiprocessing import Pool
import multiprocessing as mp
import independent_functions
import dependent_functions
#==============================================================================
#MAIN SCRIPT
#==============================================================================
#If user provided arguments, then parse them
sbatch = parser.parse_args().sbatch
configfile = parser.parse_args().config
config = configparser.ConfigParser(
                            interpolation=configparser.ExtendedInterpolation())
config.read(configfile)

#If user specifies the --sbatch flag, then we first create the output 
#directories then run the sbatch script with the 'SUBMITTED' command submitted 
#to the --sbatch flag so we know not to remake output directories. If --sbatch 
#flag not specified, simply make output directories and continue.
if sbatch == False:
    output,tempdir,figuredir,e_and_o = independent_functions.make_out_directories(
                                                                True,config)
elif str(sbatch) == 'SUBMITTED':
    output,tempdir,figuredir,e_and_o = independent_functions.make_out_directories(
                                                                False,config)
else:
    output,tempdir,figuredir,e_and_o = independent_functions.make_out_directories(
                                                                True,config)
    scriptdir = os.path.join(os.path.dirname(srcdirectory), 'scripts')
    script = os.path.join(scriptdir, 'run_main.sbatch')
    email = str(sbatch)
    os.system("sbatch --error=" + e_and_o + "%x.err --output=" + e_and_o 
                + "%x.out --mail-user="+email + " --export=cmd='" 
                + srcdirectory+" --config " +configfile 
                + " --sbatch SUBMITTED' " + script)

    sys.exit("TFEA has been submitted using an sbatch script, use qstat to \
            check its progress.")


#Run the config_parser script which will create variables for all folders and 
#paths to use throughout TFEA
independent_functions.parse_config(srcdirectory=srcdirectory,config=config,
                                    output=output,tempdir=tempdir,
                                    figuredir=figuredir)

#Verify config file to make sure user has inputted all necessary variables
independent_functions.verify_config()

#Import config file once it's created
import config

#This module takes the input list of BED files, concatenates them, and then 
#merges them via bedtools.
COMBINEtime = time.time()
if config.COMBINE:
    bedfile = independent_functions.merge_bed()
else:
    bedfile = config.BEDS[0]
COMBINEtime = time.time()-COMBINEtime 

#This module counts reads from all Bam files in BAM1 and BAM2 and creates 
#count_file with this info.
COUNTtime = time.time()
if config.COUNT:
    print "Counting reads in regions..."
    independent_functions.count_reads(bedfile=bedfile)
    print "done"
COUNTtime = time.time()-COUNTtime

#This module runs DESeq on the count_file produced above - this is used to rank
#inputted regions
DESEQtime = time.time()
if config.DESEQ:
    print "Running DESeq..."
    dependent_functions.deseq_run()
    independent_functions.rank_deseqfile()
    print "done"
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
    dependent_functions. calculate(COMBINEtime=COMBINEtime, 
                                    COUNTtime=COUNTtime, DESEQtime=DESEQtime, 
                                    CALCULATEtime=CALCULATEtime)


#Here we simply remove large bed files that are produced within this package. 
#This option can be turned on/off within the config file.
if not config.TEMP:
    print "Removing temporary bed files..."
    os.system("rm " + tempdir + '*.sorted.distance.bed')
    os.system("rm " + tempdir + '*.fa')

print "done"
