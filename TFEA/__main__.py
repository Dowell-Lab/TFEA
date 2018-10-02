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

#If user provided arguments, then parse them
sbatch = parser.parse_args().sbatch
configfile = parser.parse_args().config
config_object = configparser.ConfigParser(
                            interpolation=configparser.ExtendedInterpolation())
config_object.read(configfile)
#==============================================================================
#Functions
#==============================================================================
def make_out_directories(dirs=False, config_object=None):
    '''Creates output directories in a user-specified location where all TFEA 
        outputs will go.

    Parameters
    ----------
    dirs : boolean
        determines whether output folders will be created or not (default: 
        False)
        
    config : dict
        a configparser object that contains variables within the config file

    Returns
    -------
    output : string 
        full path to the parent output directory

    tempdir : string
        full path to the temporary directory where files are stored

    figuredir : string
        full path to the directory containing figures and plots
        
    e_and_o : string
        full path to the directory that stores stdout and stderr files
    '''
    #Output directory
    output = config_object['DATA']['OUTPUT'].strip("'")
    label1 = config_object['DATA']['LABEL1'].strip("'")
    label2 = config_object['DATA']['LABEL2'].strip("'")
    outfoldername = 'TFEA_'+label1+'-'+label2+'_'
    if dirs:
        if not os.path.isdir(os.path.join(output, outfoldername + '0')):
            output = os.path.join(output, outfoldername + '0')
            os.makedirs(output)
        else:
            outputfolders = list()
            for folder in os.listdir(output):
                if outfoldername in folder:
                    outputfolders.append(int(folder.split('_')[-1]))
            output = os.path.join(output, outfoldername + str(max(outputfolders)+1))
            os.makedirs(output)
    else:
        outputfolders = list()
        for folder in os.listdir(output):
            if outfoldername in folder:
                outputfolders.append(int(folder.split('_')[-1]))
        output = os.path.join(output, outfoldername + str(max(outputfolders)))


    #Temporary files will go in this directory
    tempdir = os.path.join(output, 'temp_files')
    if dirs:
        if not os.path.isdir(tempdir):
            os.makedirs(tempdir)

    #Error and out files will go in this directory
    e_and_o = os.path.join(output, 'e_and_o')
    if dirs:
        if not os.path.isdir(e_and_o):
            os.makedirs(e_and_o)


    #Directory where plots used in html file will be stored.
    figuredir = os.path.join(output, 'plots')
    if dirs:
        if not os.path.isdir(figuredir):
            os.makedirs(figuredir)

    return output,tempdir,figuredir,e_and_o
#==============================================================================
#CREATING DIRECTORIES
#==============================================================================
#If user specifies the --sbatch flag, then we first create the output 
#directories then run the sbatch script with the 'SUBMITTED' command submitted 
#to the --sbatch flag so we know not to remake output directories. If --sbatch 
#flag not specified, simply make output directories and continue.
if sbatch == False:
    output, tempdir, figuredir, e_and_o = make_out_directories(dirs=True, 
                                                config_object=config_object)
elif str(sbatch) == 'SUBMITTED':
    output, tempdir, figuredir, e_and_o = make_out_directories(dirs=False, 
                                                config_object=config_object)
else:
    output, tempdir, figuredir, e_and_o = make_out_directories(dirs=True, 
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
#Run the config_parser script which will create variables for all folders and 
#paths to use throughout TFEA
independent_functions.parse_config(srcdirectory=srcdirectory, 
                                    config_object=config_object,
                                    output=output,tempdir=tempdir,
                                    figuredir=figuredir)

#Verify config file to make sure user has inputted all necessary variables
# config_dict = independent_functions.verify_config_object(config=config)
independent_functions.verify_config_file()

#Import config file once it's created
import config

#This module takes the input list of BED files, concatenates them, and then 
#merges them via bedtools.
COMBINEtime = time.time()
if config.COMBINE:
    # bedfile = independent_functions.merge_bed(beds=config.BEDS, 
    #                                             tempdir=tempdir)
    # bedfile = independent_functions.tfit_clean_merge(beds=config.BEDS, 
    #                                             tempdir=tempdir)
    bedfile = independent_functions.intersect_merge_bed(bed1=config.BEDS[0:2],
                                                        bed2=config.BEDS[2:4], 
                                                        tempdir=tempdir)
else:
    bedfile = config.BEDS[0]
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
else:
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
else:
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
    os.system("rm " + os.path.join(tempdir, '*.sorted.distance.bed'))
    os.system("rm " + os.path.join(tempdir, '*.fa'))

print "done"
