#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''This is the main script that runs when the TFEA package is run.
'''
#==============================================================================
__author__ = 'Jonathan D. Rubin and Rutendo F. Sigauke'
__credits__ = ['Jonathan D. Rubin', 'Rutendo F. Sigauke', 'Jacob T. Stanley',
                'Robin D. Dowell']
__maintainer__ = 'Jonathan D. Rubin'
__email__ = 'Jonathan.Rubin@colorado.edu'
__version__ = '4.0'

#PRIMARY IMPORTS
#==============================================================================
import os
import sys
import time
import shutil
import pathlib
import argparse
import subprocess
import configparser

import preconfig_functions
from exceptions import SubprocessError

#ARGUMENT PARSING
#==============================================================================
'''We begin by parsing user arguments. TFEA requires a config file (.ini) 
    containing arguments needed to properly run TFEA. Additionally, a user may 
    specify the --sbatch flag and provide an email address to submit an sbatch 
    job to the currently logged in compute cluster. Finally, if a user simply 
    wants to test to make sure TFEA is working properly, they can specify the 
    --test flag. This testing also checks that the config file is formatted 
    properly and contains all required arguments. A test.ini config file is 
    provided within the test directory of TFEA.
'''
#TFEA source directory
srcdirectory = pathlib.Path(__file__).absolute().parent

#argparse to add arguments to this python package
parser = argparse.ArgumentParser(description=("Transcription Factor Enrichment "
                                    "Analysis (TFEA) takes as input a "
                                    "configuration file (.ini) and outputs "
                                    "a folder containing TFEA results."),
                                usage=("TFEA --config CONFIG.ini [--sbatch "
                                    "email@address.com]"))

parser.add_argument('--config','-c', default=False, metavar='',help=("REQUIRED. "
                        "A configuration file containing .ini suffix "
                        "(ex. config.ini). See example in the examples folder."))

parser.add_argument('--sbatch','-s', default=False, metavar='',help=("OPTIONAL. "
                        "Submits an sbatch job. If specified, input an e-mail "
                        "address."))

parser.add_argument('--test','-t', default=False, metavar='', action='store_const',
                        const=True, help=("OPTIONAL. "
                        "Tests TFEA to make sure it is running properly. Also "
                        "verifies that the config file has all appropriate "
                        "arguments."))

#Display help message when no args are passed.
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

#If user provided arguments, then parse them
sbatch = parser.parse_args().sbatch
configfile = parser.parse_args().config
test = parser.parse_args().test
if configfile != False:
    config_object = configparser.ConfigParser(
                            interpolation=configparser.ExtendedInterpolation())
    config_object.read(configfile)

#TEST module
#==============================================================================
'''This module runs test cases on TFEA module. Makes sure everything is installed
    and running properly. It is recommended to run this before use.
'''
if test:
    subprocess.call(["python3", "-m", 
                    "unittest", "-v", "-f", 
                    srcdirectory / 'test' / 'test_basic.py'])
    sys.exit()

#CREATING DIRECTORIES
#==============================================================================
'''TFEA creates temporary files (can be saved with the TEMP switch within the
    config file), plots, standard error and out files, and a .txt or 
    .html output. These are all stored within a parent output directory
    specified within the config file and subdirectories for each class of file.. 
'''
#If user specifies the --sbatch flag, then we first create the output 
#directories then run the sbatch script with the 'SUBMITTED' command submitted 
#to the --sbatch flag so we know not to remake output directories. If --sbatch 
#flag not specified, simply make output directories and continue.
if sbatch == False: #No flag
    outputdir, tempdir, figuredir, e_and_o = preconfig_functions.make_out_directories(
                                                dirs=True, 
                                                config_object=config_object)
    try:
        shutil.copy2(configfile, outputdir) #Copy config file to output
    except shutil.SameFileError:
        pass
elif str(sbatch) == 'SUBMITTED': #Internal flag
    outputdir, tempdir, figuredir, e_and_o = preconfig_functions.make_out_directories(
                                                dirs=False, 
                                                config_object=config_object)
else: #--sbatch specified
    outputdir, tempdir, figuredir, e_and_o = preconfig_functions.make_out_directories(
                                                dirs=True, 
                                                config_object=config_object)
    try:
        shutil.copy2(configfile, outputdir) #Copy config file to output
    except shutil.SameFileError:
        pass
    scriptdir = srcdirectory.parent / 'cluster_scripts'
    script = scriptdir / 'run_main.sbatch'
    email = str(sbatch)
    error_file = e_and_o / 'TFEA.err'
    # error_file.unlink()
    try:
        sbatch_out = subprocess.run(["sbatch", 
                            "--error=" + (e_and_o / "%x.err").as_posix(), 
                            "--output=" + (e_and_o / "%x.out").as_posix(), 
                            "--mail-user=" + email, 
                            "--export=cmd=" + srcdirectory.as_posix() 
                            + " --config " + configfile + " --sbatch SUBMITTED", 
                            script], stderr=subprocess.PIPE, 
                            stdout=subprocess.PIPE, check=True)
    except subprocess.CalledProcessError as e:
        raise SubprocessError(e.stderr.decode())
    
    (tempdir / 'jobid.txt').write_text(sbatch_out.stdout.decode().split()[-1])
    print(("TFEA has been submitted using an sbatch script. \nIt can be "
            "monitored using:\ntail -f " + error_file.as_posix()))
    sys.exit()

#VERIFICATION OF CONFIG FILE
#==============================================================================
'''This module will run but is still awaiting TFEA inputs to be finalized
'''
#This second verification ensures that at time of sbatch submission, the config
#file is correct (in case a second instance of TFEA is running)
preconfig_functions.parse_config(srcdirectory=srcdirectory, 
                                    config_object=config_object,
                                    outputdir=outputdir,tempdir=tempdir,
                                    figuredir=figuredir)

#Verify config file to make sure user has inputted all necessary variables
# config_dict = independent_functions.verify_config_object(config=config)
preconfig_functions.verify_config_file()

#==============================================================================
#MAIN SCRIPT
#==============================================================================
print("TFEA start: ", file=sys.stderr)

#SECONDARY IMPORTS
#==============================================================================
import time
import datetime
import resource
from multiprocessing import Pool
import multiprocessing as mp

import config
import multiprocess
from exceptions import InputError

#Print multiprocessing information to stderr
if config.DEBUG:
    mp.log_to_stderr()
    multiprocess.current_mem_usage()

#COMBINE module
#==============================================================================
'''This module is a pre-processing step where a user may specify how to handle
    multiple bed file inputs. The goal is to arrive at a single bed file to
    input into subsequent modules.
'''
COMBINEtime = time.time()
if config.COMBINE != False:
    print("Combining Regions...", end=' ', flush=True, file=sys.stderr)
    import combine
    if config.MD:
        bedfile, md_bedfile1, md_bedfile2 = combine.main() #Main COMBINE function
    else:
        bedfile = combine.main() #Main COMBINE without MD specified
        md_bedfile1 = None
        md_bedfile2 = None
    
    COMBINEtime = time.time()-COMBINEtime
    print("done in: " + str(datetime.timedelta(seconds=int(COMBINEtime))), file=sys.stderr)
else: #If COMBINE switch turned off, then check config file for necessary files
    if config.COUNT:
        bedfile = config.COMBINE_FILE
    if config.MD and config.FASTA:
        try:
            md_bedfile1 = config.MD_BEDFILE1
            md_bedfile2 = config.MD_BEDFILE2
        except:
            raise InputError("MD module specified without COMBINE module but no MD bed files inputted.")
    COMBINEtime = time.time()-COMBINEtime
if config.DEBUG:
    multiprocess.current_mem_usage()

#COUNT module
#==============================================================================
'''This module counts reads over a bed file
'''
COUNTtime = time.time()
if config.COUNT:
    print("Counting reads in regions...", end=' ', flush=True, file=sys.stderr)
    import count
    count_file = count.main(bedfile=bedfile) #Main COUNT function
    COUNTtime = time.time()-COUNTtime
    print("done in: " + str(datetime.timedelta(seconds=int(COUNTtime))), file=sys.stderr)
elif config.RANK == 'deseq':
    try:
        count_file = config.COUNT_FILE
    except:
        raise InputError("DE-Seq analysis specified but COUNT module set to False and COUNT_FILE is missing.")
    COUNTtime = time.time()-COUNTtime
else:
    COUNTtime = time.time()-COUNTtime
if config.DEBUG:
    multiprocess.current_mem_usage()

#RANK module
#==============================================================================
'''This module decides how to rank regions within the bed files. If genome
    hits specified then the ranked output will only contain the center of each
    region (since we will perform bedtools closest later)
'''
RANKtime = time.time()
if config.RANK != False:
    print("Ranking regions...", end=' ', flush=True, file=sys.stderr)
    import rank
    ranked_file = rank.main(count_file=count_file)
    if config.MDD:
        mdd_bedfile1, mdd_bedfile2 = rank.create_mdd_files(
                                                    ranked_file=ranked_file,
                                                    tempdir=config.TEMPDIR,
                                                    percent=0.2)
    RANKtime = time.time()-RANKtime
    print("done in: " + str(datetime.timedelta(seconds=int(RANKtime))), file=sys.stderr)
elif config.FASTA:
    try:
        ranked_file = config.RANKED_FILE
        mdd_bedfile1 = config.MDD_BEDFILE1
        mdd_bedfile2 = config.MDD_BEDFILE2
    except:
        raise InputError("FASTA conversion specified but RANK module set to False and RANKED_FILE is missing.")
    RANKtime = time.time()-RANKtime
else:
    ranked_file = None
    mdd_bedfile1 = None
    mdd_bedfile2 = None
    RANKtime = time.time()-RANKtime
if config.DEBUG:
    multiprocess.current_mem_usage()

#FASTA module
#==============================================================================
'''This module converts bed files to fasta, only used if scanning on
    the fly with fimo or homer
'''
FASTAtime = time.time()
if config.FASTA and config.SCANNER != 'genome hits':
    print("Generating fasta file...", end=' ', flush=True, file=sys.stderr)
    import fasta
    fasta_file = fasta.main(ranked_file=ranked_file, outname='ranked_file.fa')
    if config.MD:
        md_fasta1 = fasta.main(ranked_file=md_bedfile1, outname='md1_fasta.fa')
        md_fasta2 = fasta.main(ranked_file=md_bedfile2, outname='md2_fasta.fa')
    if config.MDD:
        mdd_fasta1 = fasta.main(ranked_file=mdd_bedfile1, outname='mdd1_fasta.fa')
        mdd_fasta2 = fasta.main(ranked_file=mdd_bedfile2, outname='mdd2_fasta.fa')
    FASTAtime = time.time()-FASTAtime
    print("done in: " + str(datetime.timedelta(seconds=int(RANKtime))), file=sys.stderr)
elif config.SCANNER != 'genome hits':
    if config.SCANNER != False:
        try:
            fasta_file = config.FASTA_FILE
        except:
            raise InputError("SCANNER module specified but FASTA set to False and FASTA_FILE not provided.")
    if config.MD:
        try:
            md_fasta1 = config.MD_FASTA1
            md_fasta2 = config.MD_FASTA2
        except:
            raise InputError("MD module specified without FASTA module but MD_FASTA1 and MD_FASTA2 not provided.")
    FASTAtime = time.time()-FASTAtime
else:
    fasta_file = None
    md_fasta1 = None
    md_fasta2 = None
    mdd_fasta1 = None
    mdd_fasta2 = None
    FASTAtime = time.time()-FASTAtime
if config.DEBUG:
    multiprocess.current_mem_usage()

#SCANNER module
#==============================================================================
'''This module returns motif distances to regions of interest. This is
    accomplished either by scanning regions on the fly using fimo or homer, or 
    by running bedtools closest on region centers compared to a database of
    motif hits across the genome.
'''
SCANNERtime = time.time()
if config.SCANNER != False:
    print("Scanning regions using " + config.SCANNER + "...", end=' ', flush=True, file=sys.stderr)
    import scanner
    motif_distances = scanner.main(fasta_file=fasta_file,  
                                        ranked_file=ranked_file)
    if config.MD:
        # motif_distances, md_distances1, md_distances2 = scanner.main(
        #                                                 fasta_file=fasta_file, 
        #                                                 md_fasta1=md_fasta1, 
        #                                                 md_fasta2=md_fasta2, 
        #                                                 md_bedfile1=md_bedfile1,
        #                                                 md_bedfile2=md_bedfile2,
        #                                                 ranked_file=ranked_file)
        md_distances1 = scanner.main(fasta_file=md_fasta1,  
                                        ranked_file=md_bedfile1)
        md_distances2 = scanner.main(fasta_file=md_fasta2,  
                                        ranked_file=md_bedfile2)
    if config.MDD:
        mdd_distances1 = scanner.main(fasta_file=mdd_fasta1,  
                                        ranked_file=mdd_bedfile1)
        mdd_distances2 = scanner.main(fasta_file=mdd_fasta2,  
                                        ranked_file=mdd_bedfile2)
    SCANNERtime = time.time()-SCANNERtime
    print("done in: " + str(datetime.timedelta(seconds=int(SCANNERtime))), file=sys.stderr)
else:
    SCANNERtime = time.time()-SCANNERtime
if config.DEBUG:
    multiprocess.current_mem_usage()
    
#ENRICHMENT module
#==============================================================================
'''Where the bulk of TFEA analysis occurs. Some components of plotting module 
    are contained within this enrichment module
'''
ENRICHMENTtime = time.time()
if config.ENRICHMENT != False:
    print("Calculating enrichment...", end=' ', flush=True, file=sys.stderr)
    import enrichment
    results = enrichment.tfea(motif_distances=motif_distances)
    if config.MD:
        md_results = enrichment.md(md_distances1=md_distances1, 
                                    md_distances2=md_distances2)
    if config.MDD:
        mdd_results = enrichment.md(md_distances1=mdd_distances1, 
                                    md_distances2=mdd_distances2)
    ENRICHMENTtime = time.time()-ENRICHMENTtime
    print("done in: " + str(datetime.timedelta(seconds=int(SCANNERtime))), file=sys.stderr)
else:
    ENRICHMENTtime = time.time()-ENRICHMENTtime
if config.DEBUG:
    multiprocess.current_mem_usage()

#OUTPUT module
#==============================================================================
'''A module to write output to either a txt or html file
'''
OUTPUTtime = time.time()
if config.OUTPUT_TYPE != False:
    print("Creating output...", end=' ', flush=True, file=sys.stderr)
    import output
    output.txt_output(results=results, outname='results.txt', 
                        sortindex=[-2,-1])
    if config.MD:
        header = ['TF', 'MD-Score', 'Events', 'p-val']
        output.txt_output(results=md_results, outname='md_results.txt', 
                            header=header, sortindex=[-1])
    if config.MDD:
        header = ['TF', 'MD-Score', 'Events', 'p-val']
        output.txt_output(results=mdd_results, outname='mdd_results.txt', 
                            header=header, sortindex=[-1])
    OUTPUTtime = time.time()-OUTPUTtime
    if config.OUTPUT_TYPE == 'plots':
        output.plot_output()
    if config.OUTPUT_TYPE == 'html':
        output.summary_html_output(config_object=config_object, 
                                    outputdir=outputdir)
        module_list = [('COMBINE', config.COMBINE, COMBINEtime), 
                        ('COUNT', config.COUNT, COUNTtime), 
                        ('RANK', config.RANK, RANKtime), 
                        ('FASTA', config.FASTA, FASTAtime), 
                        ('SCANNER', config.SCANNER, SCANNERtime), 
                        ('ENRICHMENT', config.ENRICHMENT, ENRICHMENTtime), 
                        ('OUTPUT', config.OUTPUT_TYPE, OUTPUTtime)]
        output.html_output(results=results, module_list=module_list)
        
    print("done in: " + str(datetime.timedelta(seconds=int(OUTPUTtime))), file=sys.stderr)
    if config.DEBUG:
        multiprocess.current_mem_usage()
else:
    raise InputError("Warning: OUTPUT_TYPE either not specified or not supported") 

print("TFEA done.", file=sys.stderr)

#REMOVE TEMP FILES
#==============================================================================
# if not config.TEMP:
#     print("Removing temporary bed files...")
#     files_to_keep = ['count_file.bed', 'DESeq.R', 'DESeq.res.txt', 
#                     'DESeq.Rout', 'ranked_file.bed', 'markov_background.txt']
#     for file1 in os.listdir(tempdir):
#         if file1 not in files_to_keep:
#              os.remove(os.path.join(tempdir, file1))

#     print("\t", "done")

#==============================================================================
#END OF MAIN SCRIPT
#==============================================================================