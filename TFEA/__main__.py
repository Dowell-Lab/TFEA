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
import argparse
import subprocess
import configparser

import preconfig_functions

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
srcdirectory = os.path.dirname(os.path.realpath(__file__))

#argparse to add arguments to this python package
parser = argparse.ArgumentParser(description=("Transcription Factor Enrichment "
                                    "Analysis (TFEA) takes as input a "
                                    "configuration file (.ini) and outputs "
                                    "a folder containing TFEA results."),
                                usage=("TFEA --config CONFIG.ini [--sbatch "
                                    "email@address.com]"))

parser.add_argument('--config','-c',default=False, metavar='',help=("REQUIRED. "
                        "A configuration file containing .ini suffix "
                        "(ex. config.ini). See example in the examples folder."))

parser.add_argument('--sbatch','-s',default=False,metavar='',help=("OPTIONAL. "
                        "Submits an sbatch job. If specified, input an e-mail "
                        "address."))

parser.add_argument('--test','-t', default=False, metavar='', help=("OPTIONAL. "
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
config_object = configparser.ConfigParser(
                            interpolation=configparser.ExtendedInterpolation())
config_object.read(configfile)

#CREATING DIRECTORIES
#==============================================================================
'''TFEA creates temporary files (can be saved with the TEMP switch within the
    config file), plots, standard error and out files, and a .txt or 
    .html output. These are all stored within a parent output directory
    specified within the config file and subdirectories for each class of file.
    If the specified output directory exists, TFEA will remove that directory 
    and all its contents and create new directories. 
'''
#If user specifies the --sbatch flag, then we first create the output 
#directories then run the sbatch script with the 'SUBMITTED' command submitted 
#to the --sbatch flag so we know not to remake output directories. If --sbatch 
#flag not specified, simply make output directories and continue.
if sbatch == False:
    outputdir, tempdir, figuredir, e_and_o = preconfig_functions.make_out_directories(
                                                dirs=True, 
                                                config_object=config_object)
    try:
        shutil.copy2(configfile, outputdir)
    except shutil.SameFileError:
        pass
elif str(sbatch) == 'SUBMITTED':
    outputdir, tempdir, figuredir, e_and_o = preconfig_functions.make_out_directories(
                                                dirs=False, 
                                                config_object=config_object)
else:
    outputdir, tempdir, figuredir, e_and_o = preconfig_functions.make_out_directories(
                                                dirs=True, 
                                                config_object=config_object)
    try:
        shutil.copy2(configfile, outputdir)
    except shutil.SameFileError:
        pass
    scriptdir = os.path.join(os.path.dirname(srcdirectory), 'scripts')
    script = os.path.join(scriptdir, 'run_main.sbatch')
    email = str(sbatch)
    os.system("sbatch --error=" + e_and_o + "/%x.err --output=" + e_and_o 
                + "/%x.out --mail-user=" + email + " --export=cmd='" 
                + srcdirectory+" --config " + configfile 
                + " --sbatch SUBMITTED' " + script)

    sys.exit(("TFEA has been submitted using an sbatch script, use qstat to "
            "check its progress."))

#TEST module
#==============================================================================
'''If a user
'''
if test:
    preconfig_functions.parse_config(srcdirectory=srcdirectory, 
                                    config_object=config_object,
                                    outputdir=outputdir,tempdir=tempdir,
                                    figuredir=figuredir)

    preconfig_functions.verify_config_file()
    sys.exit()

#VERIFICATION OF CONFIG FILE
#==============================================================================
'''
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

#SECONDARY IMPORTS
#==============================================================================
from multiprocessing import Pool
import multiprocessing as mp
import config
'''
'''
#Calculate how many cpus can be used for later parallelization steps
cpus = mp.cpu_count()
if cpus > 64:
    cpus = 64

#COMBINE module
#==============================================================================
'''
'''
COMBINEtime = time.time()
if config.COMBINE != False:
    print("Combining Regions...", end=' ', flush=True)
    import combine_functions
    if config.COMBINE == 'merge all':
        bedfile = combine_functions.merge_bed(beds=config.BED1+config.BED2, 
                                                tempdir=tempdir)
        if config.MD:
            md_bedfile1, md_bedfile2 = combine_functions.combine_md(
                                                    bed1=config.BED1, 
                                                    bed2=config.BED2, 
                                                    tempdir=tempdir, 
                                                    method=config.COMBINE, 
                                                    size_cut=200)
    elif config.COMBINE == 'tfit clean':
        bedfile = combine_functions.tfit_clean(
                                beds=config.BED1+config.BED2, tempdir=tempdir)
        if config.MD:
            md_bedfile1, md_bedfile2 = combine_functions.combine_md(
                                                    bed1=config.BED1, 
                                                    bed2=config.BED2, 
                                                    tempdir=tempdir, 
                                                    method=config.COMBINE, 
                                                    size_cut=200)
    elif config.COMBINE == 'tfit remove small':
        bedfile = combine_functions.tfit_remove_small(
                                beds=config.BED1+config.BED2, tempdir=tempdir)
        if config.MD:
            md_bedfile1, md_bedfile2 = combine_functions.combine_md(
                                                    bed1=config.BED1, 
                                                    bed2=config.BED2, 
                                                    tempdir=tempdir, 
                                                    method=config.COMBINE, 
                                                    size_cut=200)
    elif config.COMBINE == 'intersect/merge':
        bedfile = combine_functions.intersect_merge_bed(bed1=config.BED1,
                                                        bed2=config.BED2, 
                                                        tempdir=tempdir)
        if config.MD:
            md_bedfile1, md_bedfile2 = combine_functions.combine_md(
                                                    bed1=config.BED1, 
                                                    bed2=config.BED2, 
                                                    tempdir=tempdir, 
                                                    method=config.COMBINE, 
                                                    size_cut=200, 
                                                    largewindow=config.LARGEWINDOW)
    print("done")
else:
    if config.COUNT:
        bedfile = config.COMBINE_FILE
    if config.MD and config.FASTA:
        try:
            md_bedfile1 = config.MD_BEDFILE1
            md_bedfile2 = config.MD_BEDFILE2
        except:
            sys.exit(1, "MD module specified without COMBINE module but no MD bed files inputted.")
COMBINEtime = time.time()-COMBINEtime

#COUNT module
#==============================================================================
'''
'''
COUNTtime = time.time()
if config.COUNT:
    print("Counting reads in regions...", end=' ', flush=True)
    import count_functions
    
    count_file = count_functions.count_reads(bedfile=bedfile, 
                                                bam1=config.BAM1, 
                                                bam2=config.BAM2, 
                                                tempdir=tempdir, 
                                                label1=config.LABEL1, 
                                                label2=config.LABEL2)
    sample_number = (len(config.BAM1)+len(config.BAM2))
    millions_mapped = count_functions.sum_reads(count_file=count_file, 
                                                sample_number=sample_number)
    print("done")
elif config.DESEQ:
    #If you don't want to perform multibamcov but still want to perform
    # DE-Seq, user must provide a count_file
    count_file = config.COUNT_FILE
COUNTtime = time.time()-COUNTtime

#DESEQ module
#==============================================================================
'''
'''
DESEQtime = time.time()
if config.DESEQ:
    print("Running DESeq...", end=' ', flush=True)
    import deseq_functions
    
    deseq_functions.write_deseq_script(bam1=config.BAM1, bam2=config.BAM2, 
                                                tempdir=tempdir, 
                                                count_file=count_file,
                                                label1=config.LABEL1, 
                                                label2=config.LABEL2)
    with open(os.path.join(tempdir, 'DESeq.Rout'), 'w') as stdout:
        subprocess.call("R < " + os.path.join(tempdir, "DESeq.R") + " --no-save", 
                        shell=True, stdout=stdout)
    deseq_file = os.path.join(tempdir, 'DESeq.res.txt')

    # deseq_functions.plot_deseq_MA(deseq_file=deseq_file, label1=config.LABEL1, 
    #                                 label2=config.LABEL2, figuredir=figuredir)
    print("done")
elif config.RANK != False:
    #If user does not want to perform DE-Seq but still wants TFEA to caluclate
    # TF enrichment, user must provide a ranked_center_file
    deseq_file = config.DESEQ_FILE
DESEQtime = time.time()-DESEQtime

#RANK module
#==============================================================================
'''
'''
RANKtime = time.time()
if config.RANK != False:
    if config.RANK == 'deseq':
        print("Ranking regions...", end=' ', flush=True)
        import rank_functions
        if config.SCANNER == 'genome hits':
            ranked_file = rank_functions.deseq_center(deseq_file=deseq_file,
                                                    tempdir=tempdir)
        else:
            ranked_file = rank_functions.deseq(deseq_file=deseq_file, 
                                                tempdir=tempdir, 
                                                largewindow=config.LARGEWINDOW)
        print("done")
elif config.FASTA:
    ranked_file = config.RANKED_FILE
RANKtime = time.time()-RANKtime

#FASTA module
#==============================================================================
'''
'''
FASTAtime = time.time()
if config.FASTA:
    print("Generating fasta file...", end=' ', flush=True)
    import fasta_functions
    fasta_file = fasta_functions.getfasta(bedfile=ranked_file, 
                                            genomefasta=config.GENOMEFASTA, 
                                            tempdir=tempdir, 
                                            outname='fasta_file.fa')
    if config.MD:
        md_fasta1 = fasta_functions.getfasta(bedfile=md_bedfile1, 
                                                genomefasta=config.GENOMEFASTA, 
                                                tempdir=tempdir, 
                                                outname='md_bedfile1.fa')
        md_fasta2 = fasta_functions.getfasta(bedfile=md_bedfile2, 
                                                genomefasta=config.GENOMEFASTA, 
                                                tempdir=tempdir, 
                                                outname='md_bedfile2.fa')
    
    print("done")
else:
    if config.SCANNER != False:
        fasta_file = config.FASTA_FILE
    if config.MD:
        try:
            md_fasta1 = config.MD_FASTA1
            md_fasta2 = config.MD_FASTA2
        except:
            sys.exit(1, "MD module specified without FASTA module but no MD fasta files inputted.")
FASTAtime = time.time()-FASTAtime

#SCANNER module
#==============================================================================
'''
'''
SCANNERtime = time.time()
if config.SCANNER != False:
    print("Scanning regions using " + config.SCANNER + "...", end=' ', flush=True)
    import scanner_functions

    #FIMO
    if config.SCANNER == 'fimo':
        if config.FIMO_BACKGROUND == 'largewindow':
            background_file = scanner_functions.fimo_background_file(
                                window=int(config.LARGEWINDOW), 
                                tempdir=tempdir, bedfile=ranked_file, 
                                genomefasta=config.GENOMEFASTA, order='1')
        elif config.FIMO_BACKGROUND == 'smallwindow':
            background_file = scanner_functions.fimo_background_file(
                                window=int(config.SMALLWINDOW), 
                                tempdir=tempdir, bedfile=ranked_file, 
                                genomefasta=config.GENOMEFASTA, order='1')
        elif type(config.FIMO_BACKGROUND) == int:
            background_file = scanner_functions.fimo_background_file(
                                window=config.FIMO_BACKGROUND, 
                                tempdir=tempdir, bedfile=ranked_file, 
                                genomefasta=config.GENOMEFASTA, order='1')
        elif type(config.FIMO_BACKGROUND) == str:
            background_file = config.FIMO_BACKGROUND
        else:
            background_file = None

        if config.SINGLEMOTIF != False:
            motif_list = [config.SINGLEMOTIF]
        else:
            motif_list = scanner_functions.fimo_motif_names(
                                            motifdatabase=config.FIMO_MOTIFS)

        fimo_keywords = dict(bg_file=background_file, fasta_file=fasta_file, 
                            tempdir=tempdir, motifdatabase=config.FIMO_MOTIFS, 
                            thresh=config.FIMO_THRESH, 
                            largewindow=config.LARGEWINDOW)
        p = Pool(cpus)
        motif_distances = list()
        for motif in motif_list:
            motif_distances.append(p.apply_async(scanner_functions.fimo, 
                                            (motif,), fimo_keywords).get())
        if config.MD:
            fimo_keywords = dict(bg_file=background_file, fasta_file=md_fasta1, 
                            tempdir=tempdir, motifdatabase=config.FIMO_MOTIFS, 
                            thresh=config.FIMO_THRESH, 
                            largewindow=config.LARGEWINDOW)
            p = Pool(cpus)
            md_distances1 = list()
            for motif in motif_list:
                md_distances1.append(p.apply_async(scanner_functions.fimo, 
                                                (motif,), fimo_keywords).get())
            
            fimo_keywords = dict(bg_file=background_file, fasta_file=md_fasta2, 
                            tempdir=tempdir, motifdatabase=config.FIMO_MOTIFS, 
                            thresh=config.FIMO_THRESH, 
                            largewindow=config.LARGEWINDOW)
            p = Pool(cpus)
            md_distances2 = list()
            for motif in motif_list:
                md_distances2.append(p.apply_async(scanner_functions.fimo, 
                                                (motif,), fimo_keywords).get())
            


    #HOMER
    elif config.SCANNER == 'homer':
        print("Homer scanning is incomplete")

    #GENOME HITS
    elif config.SCANNER == 'genome hits':
        if config.SINGLEMOTIF == False:
            motif_list = os.listdir(config.GENOMEHITS)
        else:
            motif_list = os.path.join(config.GENOMEHITS, 
                                    config.SINGLEMOTIF)

        bedtools_distance_keywords = dict(
                                    genomehits=config.GENOMEHITS, 
                                    ranked_center_file=ranked_file, 
                                    tempdir=tempdir, 
                                    distance_cutoff=config.LARGEWINDOW)

        p = Pool(cpus)
        motif_distances = list()
        for motif in motif_list:
            motif_distances.append(p.apply_async(
                                        scanner_functions.bedtools_closest,
                                        (motif,), 
                                        bedtools_distance_keywords).get())

    print("done")
SCANNERtime = time.time()-SCANNERtime
    
#ENRICHMENT module
#==============================================================================
'''Where the bulk of TFEA analysis occurs. Some components of plotting module 
    are contained within this enrichment module
'''
ENRICHMENTtime = time.time()
if config.ENRICHMENT != False:
    print("Calculating enrichment...", end=' ', flush=True)
    import enrichment_functions
    if config.ENRICHMENT == 'auc':
        auc_keywords = dict(output_type=config.OUTPUT_TYPE, 
                            permutations=config.PERMUTATIONS)
        p = Pool(cpus)
        results = list()
        for distances in motif_distances:
            results.append(p.apply_async(enrichment_functions.area_under_curve, 
                            (distances,), auc_keywords).get())
        enrichment_functions.padj_bonferroni(results)
    elif config.ENRICHMENT == 'anderson-darling':
        p = Pool(cpus)
        results = list()
        for distances in motif_distances:
            results.append(p.apply_async(enrichment_functions.anderson_darling, 
                            (distances,)).get())


    if config.MD:
        md_keywords = dict(smallwindow=config.SMALLWINDOW)
        p = Pool(cpus)
        md_results = list()
        for distances in zip(sorted(md_distances1),sorted(md_distances2)):
            md_results.append(p.apply_async(enrichment_functions.md_score, 
                        (distances,), md_keywords).get())
        md_results = enrichment_functions.md_score_p(md_results)

    print("done")
ENRICHMENTtime = time.time()-ENRICHMENTtime

#OUTPUT module
#==============================================================================
'''
'''
OUTPUTtime = time.time()
if config.OUTPUT_TYPE != False:
    print("Creating output...", end=' ', flush=True)
    import output_functions
    output_functions.txt_output(results=results, outputdir=outputdir, 
                                outname='results.txt', header=None, 
                                sortindex=-1)
    if config.MD:
        header = ['TF', 'MD-Score', 'Events', 'p-val']
        output_functions.txt_output(results=md_results, outputdir=outputdir, 
                                outname='md_results.txt', header=header, 
                                sortindex=-1)
    if config.OUTPUT_TYPE == 'html':
        output_functions.summary_html_output(config_object=config_object, 
                                                outputdir=outputdir)
        OUTPUTtime = time.time()-OUTPUTtime
        module_list = [('COMBINE', config.COMBINE, COMBINEtime), 
                        ('COUNT', config.COUNT, COUNTtime), 
                        ('DESEQ', config.DESEQ, DESEQtime), 
                        ('RANK', config.RANK, RANKtime), 
                        ('FASTA', config.FASTA, FASTAtime), 
                        ('SCANNER', config.SCANNER, SCANNERtime), 
                        ('ENRICHMENT', config.ENRICHMENT, ENRICHMENTtime), 
                        ('OUTPUT', config.OUTPUT_TYPE, OUTPUTtime)]
        output_functions.html_output(results=results, module_list=module_list)
        
    print("done")
else:
    print("Warning: OUTPUT_TYPE either not specified or not supported") 

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




# #Create a meta plot using input regions and input bam files
# #==============================================================================
# if config.METAPLOT:
#     print("Creating meta plots...")
#     meta_profile_dict = dependent_functions.meta_plot(
#                                     ranked_center_file=ranked_center_file, 
#                                     largewindow=config.LARGEWINDOW, 
#                                     bam1=config.BAM1, bam2=config.BAM2, 
#                                     tempdir=tempdir, cpus=cpus, 
#                                     figuredir=figuredir, label1=config.LABEL1, 
#                                     label2=config.LABEL2, dpi=config.DPI,
#                                     millions_mapped=millions_mapped)
#     print("done")
# else:
#     meta_profile_dict={}
# #Either perform motif scanning on the fly or use genome-wide motif hits 
# # provided by user
# #==============================================================================

# #GENOMEWIDEHITS refers to precomputed motif hits accross the genome, if a user
# #specifies this option, they must populate MOTIF_GENOMEWIDE_HITS with a full
# #path to a folder containing bed files where each bed file is the genomewide
# #hits of a particular motif
# if config.GENOMEWIDEHITS:
#     if config.SINGLEMOTIF == False:
#         #If all motifs should be analyzed, then the list of motif files is
#         #simply all files within MOTIF_GENOMEWIDE_HITS
#         list_motif_files = os.listdir(config.MOTIF_GENOMEWIDE_HITS)
#     else:
#         #If a single motif should be analyzed then the list of motif files is
#         #simply the full path to that bed file
#         list_motif_files = os.path.join(config.MOTIF_GENOMEWIDE_HITS, 
#                                 config.SINGLEMOTIF)

#     #Build the list of arguments to compute motif to region distances using 
#     #bedtools
#     bedtools_distance_args = [(motif_file, config.MOTIF_GENOMEWIDE_HITS, 
#                             ranked_center_file) for motif_file in list_motif_files]
#     p = Pool(cpus)

#     #Execute bedtools_distance function
#     motif_list = p.map(dependent_functions.bedtools_distance, 
#                                                         bedtools_distance_args)
# elif fastafile:
#     ranked_small_regions_fasta_file = independent_functions.smallfasta_from_fasta(
#         fastafile=fastafile, 
#         outname=os.path.join(tempdir, 'ranked_smallregions.fa'))

#     bg_file = independent_functions.fasta_markov(
#                                     fastafile=ranked_smallregions_fasta_file,
#                                     tempdir=tempdir)
# else:
#     #If a user desires TFEA to perform motif scanning, SMALLWINDOW will be used
#     #to obtain background ACGT content, and actual scanning will be performed
#     #over LARGEWINDOW (mostly for display purposes)
#     ranked_smallregions_file = independent_functions.get_regions(
#                                         tempdir=tempdir, 
#                                         ranked_center_file=ranked_center_file,
#                                         window=config.SMALLWINDOW,
#                                         outname='ranked_smallregions.bed')

#     ranked_largeregions_file = independent_functions.get_regions(
#                                         tempdir=tempdir, 
#                                         ranked_center_file=ranked_center_file,
#                                         window=config.LARGEWINDOW,
#                                         outname='ranked_largeregions.bed')

#     ranked_smallregions_fasta_file = independent_functions.getfasta(
#                                             genomefasta=config.GENOMEFASTA, 
#                                             tempdir=tempdir,
#                                             bedfile=ranked_smallregions_file,
#                                             outname='ranked_smallregions.fa')

#     ranked_largeregions_fasta_file = independent_functions.getfasta(
#                                             genomefasta=config.GENOMEFASTA, 
#                                             tempdir=tempdir,
#                                             bedfile=ranked_largeregions_file,
#                                             outname='ranked_largeregions.fa')
    
#     bg_file = independent_functions.fasta_markov(
#                                     fastafile=ranked_smallregions_fasta_file,
#                                     tempdir=tempdir)

# if config.FIMO:
#     if config.SINGLEMOTIF == False:
#         motif_list = independent_functions.get_motif_names(
#                                             motifdatabase=config.MOTIFDATABASE)
#     else:
#         motif_list = [config.SINGLEMOTIF]

#     #Perform FIMO motif scanning using parallelization
#     fimo_args = [(motif, config.LARGEWINDOW, tempdir, config.MOTIFDATABASE, 
#         figuredir, ranked_center_file, bg_file, 
#         ranked_largeregions_fasta_file) for motif in motif_list]
#     p = Pool(cpus)
#     motif_list = p.map(dependent_functions.fimo_distance, fimo_args)

# elif config.HOMER:
#     homer_out = independent_functions.homer(tempdir=tempdir, 
#                         srcdirectory=os.path.dirname(srcdirectory), 
#                         fastafile=ranked_largeregions_fasta_file, cpus=cpus,
#                         motif_file=config.HOMER_MOTIF_FILE)
    
#     motif_list = independent_functions.homer_parse(largewindow=config.LARGEWINDOW, 
#                                         tempdir=tempdir, homer_out=homer_out,
#                                         ranked_center_file=ranked_center_file)
# #==============================================================================

# #This is the bulk of the analysis of this package, it performs:
# #   1. GC distribution across regions for plotting
# #   2. Millions mapped calculation for meta eRNA plots
# #   3. Motif distance calculation to the center of inputted regions
# #   4. Enrichment score calculation via AUC method
# #   5. Random shuffle simulation and recalculation of enrichment score
# #   6. Plotting and generation of html report
# #==============================================================================
# CALCULATEtime = time.time()
# if config.CALCULATE:
#     dependent_functions.calculate(tempdir=tempdir, 
#                                     outputdir=output, 
#                                     bam1=config.BAM1, 
#                                     bam2=config.BAM2, 
#                                     label1=config.LABEL1,
#                                     label2=config.LABEL2,
#                                     motif_list=motif_list, 
#                                     singlemotif=config.SINGLEMOTIF, 
#                                     ranked_center_file=ranked_center_file, 
#                                     COMBINEtime=COMBINEtime, 
#                                     COUNTtime=COUNTtime, DESEQtime=DESEQtime, 
#                                     CALCULATEtime=CALCULATEtime, 
#                                     fimo=config.FIMO, 
#                                     plot=config.PLOTALL, 
#                                     padj_cutoff=config.PADJCUTOFF, 
#                                     logos=config.LOGOS, 
#                                     figuredir=figuredir,
#                                     largewindow=config.LARGEWINDOW, 
#                                     smallwindow=config.SMALLWINDOW,
#                                     motifdatabase=config.MOTIFDATABASE,
#                                     genomefasta=config.GENOMEFASTA, 
#                                     meta_profile_dict=meta_profile_dict,
#                                     millions_mapped=millions_mapped)
# #==============================================================================

# #Here we simply remove large bed files that are produced within this package. 
# #This option can be turned on/off within the config file.
# #==============================================================================
# if not config.TEMP:
#     print("Removing temporary bed files...")
#     files_to_keep = ['count_file.bed', 'DESeq.R', 'DESeq.res.txt', 
#                     'DESeq.Rout', 'ranked_file.bed', 'markov_background.txt', 
#                     'ranked_file.center.bed', 'ranked_file.center.sorted.bed']
#     for file1 in os.listdir(tempdir):
#         if file1 not in files_to_keep:
#              os.remove(os.path.join(tempdir, file1))

# print("done")
# #==============================================================================