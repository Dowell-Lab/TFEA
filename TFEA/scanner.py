#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''This module outputs a list of motif distances to inputted regions for all
    desired TF motifs. This ouptut can be used in the ENRICHMENT module to 
    calculate motif enrichment within input regions.
'''

#==============================================================================
__author__ = 'Jonathan D. Rubin'
__credits__ = ['Jonathan D. Rubin', 'Rutendo F. Sigauke', 'Jacob T. Stanley',
                'Robin D. Dowell']
__maintainer__ = 'Jonathan D. Rubin'
__email__ = 'Jonathan.Rubin@colorado.edu'

#Imports
#==============================================================================
import os
import sys
import time
import datetime
import subprocess
import traceback
from pathlib import Path

import numpy as np
from pybedtools import BedTool
from pybedtools import featurefuncs

from TFEA import multiprocess
from TFEA import exceptions

#Main Script
#==============================================================================
def main(use_config=True, fasta_file=False, md_fasta1=False, md_fasta2=False, 
            ranked_file=None, md_bedfile1=None, md_bedfile2=None, 
            scanner=None, md=None, largewindow=None, smallwindow=None, 
            genomehits=None, fimo_background=None, genomefasta=None, 
            tempdir=None, fimo_motifs=None, singlemotif=None, fimo_thresh=None,
            debug=None, mdd=None, jobid=None, cpus=None):
    '''This is the main script of the SCANNER module. It returns motif distances
        to regions of interest by either scanning fasta files on the fly using
        fimo or homer or by using bedtools closest on a center bed file and 
        a database of bed files corresponding to motif hits across the genome

    Parameters
    ----------
    use_config : boolean
        Whether to use a config module to assign variables.
    fasta_file : str
        Full path to a fasta file
    md_fasta1 : str
        Full path to a fasta file corresponding to a single condition. Only 
        required if md score analysis desired
    md_fasta2 : str
        Full path to a fasta file corresponding to a single condition. Only 
        required if md score analysis desired
    ranked_file : str
        Full path to a ranked bed file used in calculating background for 
        fimo scanning. Only necessary if fimo scanning desired
    scanner : str
        Scanning method desired
    md : boolean
        Whether md score analysis is desired. If True, requires bed files for
        each condition. These can be generated in the COMBINE module.
    largewindow : int
        Half-length of total window size to use when defining cutoffs for 
        how far out to measure motif distances
    smallwindow : int
        Half-length of window size to use when defining cutoffs for significant
        motif hits
    genomehits : str
        Full path to a folder containing bed files of motif hits across the 
        genome
    fimo_background : int, str, or boolean
        Defines whether to use a background file when performing fimo motif
        scanning. A user can specify any int for window size, smallwindow, 
        largewindow, or False if not desired.
    genomefasta : str
        Full path to a fasta file for desired genome
    tempdir : str
        Full path to a directory where files will be saved
    fimo_motifs : str
        Full path to a .meme formatted motif database
    singlemotif : str or boolean
        Whether to perform scanning on only a subset of motifs. A user can
        specify a single motif or a ',' separated list of motifs.
    fimo_thresh : str
        A float formatted as a string to be used when calling fimo to specify
        the p-value cutoff threshold
    debug : boolean
        Whether to print debug statements specifically within the multiprocess
        module

    Returns
    -------
    motif_distances : list of lists
        A list containing a list for each motif scanned. For each motif, the 
        list begins with the motif name as a string and is followed by int
        values corresponding to the motif distance for each region (ranked). A
        '.' value means the motif was not within the given region
    md_distances1 : list of lists
        A list containing a list for each motif scanned. For each motif, the 
        list begins with the motif name as a string and is followed by int
        values corresponding to the motif distance for each region (ranked). A
        '.' value means the motif was not within the given region
    md_distances2 : list of lists
        A list containing a list for each motif scanned. For each motif, the 
        list begins with the motif name as a string and is followed by int
        values corresponding to the motif distance for each region (ranked). A
        '.' value means the motif was not within the given region

    Raises
    ------
    InputError
        If an unknown scanner option is specified
    '''
    start_time = time.time()
    if use_config:
        from TFEA import config
        fasta_file = config.vars['FASTA_FILE']
        md_fasta1 = config.vars['MD_FASTA1']
        md_fasta2 = config.vars['MD_FASTA2']
        mdd_fasta1 = config.vars['MDD_FASTA1']
        mdd_fasta2 = config.vars['MDD_FASTA2']
        ranked_file = config.vars['RANKED_FILE']
        md_bedfile1 = config.vars['MD_BEDFILE1']
        md_bedfile2 = config.vars['MD_BEDFILE2']
        mdd_bedfile1 = config.vars['MDD_BEDFILE1']
        mdd_bedfile2 = config.vars['MDD_BEDFILE2']
        scanner = config.vars['SCANNER']
        md = config.vars['MD']
        largewindow = config.vars['LARGEWINDOW']
        smallwindow = config.vars['SMALLWINDOW']
        genomehits = config.vars['GENOMEHITS']
        fimo_background = config.vars['FIMO_BACKGROUND']
        genomefasta = config.vars['GENOMEFASTA']
        tempdir = config.vars['TEMPDIR']
        fimo_motifs = config.vars['FIMO_MOTIFS']
        singlemotif = config.vars['SINGLEMOTIF']
        fimo_thresh = config.vars['FIMO_THRESH']
        debug = config.vars['DEBUG']
        mdd = config.vars['MDD']
        mdd_pval = config.vars['MDD_PVAL']
        mdd_percent = config.vars['MDD_PERCENT']
        pvals = config.vars['PVALS']
        cpus = config.vars['CPUS']
        jobid = config.vars['JOBID']

    print("Scanning regions using " + scanner + "...", flush=True, file=sys.stderr)

    motif_distances = None
    md_distances1 = None
    md_distances2 = None
    mdd_distances1 = None
    mdd_distances2 = None

    if not fasta_file and scanner != 'genome hits':
        fasta_file = getfasta(bedfile=ranked_file, genomefasta=genomefasta, 
                                tempdir=tempdir, outname='ranked_file.fa')
        if os.stat(fasta_file).st_size == 0:
            raise exceptions.FileEmptyError("Error in SCANNER module. Converting RANKED_FILE to fasta failed.")
    if md:
        if not md_fasta1:
            md_fasta1 = getfasta(bedfile=md_bedfile1, genomefasta=genomefasta, 
                                tempdir=tempdir, outname='md1_fasta.fa')
        if not md_fasta2:
            md_fasta2 = getfasta(bedfile=md_bedfile2, genomefasta=genomefasta, 
                                tempdir=tempdir, outname='md2_fasta.fa')
        if os.stat(md_fasta1).st_size == 0 or os.stat(md_fasta2).st_size == 0:
            raise exceptions.FileEmptyError("Error in SCANNER module. Converting MD bedfiles to fasta failed.")
    if mdd:
        if not mdd_fasta1:
            mdd_fasta1 = getfasta(bedfile=mdd_bedfile1, genomefasta=genomefasta, 
                                tempdir=tempdir, outname='mdd1_fasta.fa')
        if not mdd_fasta2:
            mdd_fasta2 = getfasta(bedfile=mdd_bedfile2, genomefasta=genomefasta, 
                                tempdir=tempdir, outname='mdd2_fasta.fa')
        if os.stat(mdd_fasta1).st_size == 0 or os.stat(mdd_fasta2).st_size == 0:
            raise exceptions.FileEmptyError("Error in SCANNER module. Converting MDD bedfiles to fasta failed.")

    #FIMO
    if scanner == 'fimo':
        #Get background file, if none desired set to 'None'
        if fasta_file and fimo_background:
            background_file = fasta_markov(tempdir=tempdir, fastafile=fasta_file, order='1')
        elif fimo_background == 'largewindow':
            background_file = fimo_background_file(
                                window=int(largewindow), 
                                tempdir=tempdir, bedfile=ranked_file, 
                                genomefasta=genomefasta, order='1')
        elif fimo_background == 'smallwindow':
            background_file = fimo_background_file(
                                window=int(smallwindow), 
                                tempdir=tempdir, bedfile=ranked_file, 
                                genomefasta=genomefasta, order='1')
        elif type(fimo_background) == int:
            background_file = fimo_background_file(
                                window=fimo_background, 
                                tempdir=tempdir, bedfile=ranked_file, 
                                genomefasta=genomefasta, order='1')
        elif type(fimo_background) == str:
            background_file = fimo_background
        else:
            background_file = None

        #Get motifs to scan through
        if singlemotif != False:
            motif_list = singlemotif.split(',')
        else:
            motif_list = fimo_motif_names(motifdatabase=fimo_motifs)

        #Perform fimo on desired motifs
        print("\tTFEA:", file=sys.stderr)
        fimo_keywords = dict(bg_file=background_file, fasta_file=fasta_file, 
                            tempdir=tempdir, motifdatabase=fimo_motifs, 
                            thresh=fimo_thresh, 
                            largewindow=largewindow)

        motif_distances = multiprocess.main(function=fimo, args=motif_list, 
                                            kwargs=fimo_keywords, debug=debug, 
                                            jobid=jobid, cpus=cpus)

        #FIMO for md score fasta files
        if md:
            print("\tMD:", file=sys.stderr)
            fimo_keywords = dict(bg_file=background_file, fasta_file=md_fasta1, 
                            tempdir=tempdir, motifdatabase=fimo_motifs, 
                            thresh=fimo_thresh, 
                            largewindow=largewindow)
            md_distances1 = multiprocess.main(function=fimo, args=motif_list, 
                                                kwargs=fimo_keywords, 
                                                debug=debug, jobid=jobid, 
                                                cpus=cpus)
            
            fimo_keywords = dict(bg_file=background_file, fasta_file=md_fasta2, 
                            tempdir=tempdir, motifdatabase=fimo_motifs, 
                            thresh=fimo_thresh, 
                            largewindow=largewindow)
            md_distances2 = multiprocess.main(function=fimo, args=motif_list, 
                                                kwargs=fimo_keywords, 
                                                debug=debug, jobid=jobid,
                                                cpus=cpus)
            
            if use_config:
                config.vars['MD_DISTANCES1'] = md_distances1
                config.vars['MD_DISTANCES2'] = md_distances2
        
        if mdd:
            print("\tMDD:", file=sys.stderr)
            print(f'\t Completed: 0/{len(motif_distances)} ', end=' ', file=sys.stderr)
            fimo_keywords = dict(bg_file=background_file, fasta_file=mdd_fasta1, 
                            tempdir=tempdir, motifdatabase=fimo_motifs, 
                            thresh=fimo_thresh, 
                            largewindow=largewindow)
            mdd_distances1 = multiprocess.main(function=fimo, args=motif_list, 
                                                kwargs=fimo_keywords, 
                                                debug=debug, jobid=jobid,
                                                cpus=cpus)
            
            fimo_keywords = dict(bg_file=background_file, fasta_file=mdd_fasta2, 
                            tempdir=tempdir, motifdatabase=fimo_motifs, 
                            thresh=fimo_thresh, 
                            largewindow=largewindow)
            mdd_distances2 = multiprocess.main(function=fimo, args=motif_list, 
                                                kwargs=fimo_keywords, 
                                                debug=debug, jobid=jobid,
                                                cpus=cpus)
            # mdd_distances1 = []
            # mdd_distances2 = []
            # mdd_sorted_indices = np.argsort(pvals)
            # for i, single_motif_distances in enumerate(motif_distances, 1):
            #     motif = single_motif_distances[0]
            #     mdd_distances = single_motif_distances[1:]
            #     print("pval len:", len(pvals), file=sys.stderr)
            #     print("mdd_indices len:", len(mdd_sorted_indices), file=sys.stderr)
            #     print("mdd_dist len:", len(mdd_distances), file=sys.stderr)
            #     mdd_sorted_distances = [mdd_distances[i] for i in mdd_sorted_indices]
            #     if mdd_percent != False:
            #         cutoff = int(len(mdd_sorted_distances)*mdd_percent)
            #         mdd_distances2.append([motif] + mdd_sorted_distances[:cutoff])
            #         mdd_distances1.append([motif] + mdd_sorted_distances[cutoff:])
            #     else:
            #         sorted_pvals = [pvals[i] for i in mdd_sorted_indices]
            #         cutoff = int(len([p for p in sorted_pvals if p < mdd_pval]))
            #         mdd_distances2.append([motif] + mdd_sorted_distances[:cutoff])
            #         mdd_distances1.append([motif] + mdd_sorted_distances[cutoff:])
            #     # print(f'\r\t Completed: {i}/{len(motif_distances)} ', end=' ', flush=True, file=sys.stderr)
            if use_config:
                config.vars['MDD_DISTANCES1'] = mdd_distances1
                config.vars['MDD_DISTANCES2'] = mdd_distances2
            
    #HOMER
    elif scanner== 'homer':
        raise exceptions.InputError("Homer scanning is not supported at this time.")

    #GENOME HITS
    elif scanner == 'genome hits':
        #Get motifs to analyze
        if singlemotif == False:
            motif_list = os.listdir(genomehits)
        else:
            motif_list = [os.path.join(genomehits, motif) for motif in singlemotif.split(',')]

        #Perform bedtools closest to get distances
        ranked_file = get_center(bedfile=ranked_file, outname=ranked_file)
        print("\tTFEA:", file=sys.stderr)
        bedtools_distance_keywords = dict(genomehits=genomehits, 
                                            ranked_center_file=ranked_file, 
                                            tempdir=tempdir, 
                                            distance_cutoff=largewindow, 
                                            rank_index=3)
        
        motif_distances = multiprocess.main(function=bedtools_closest, 
                                            args=motif_list, 
                                            kwargs=bedtools_distance_keywords, 
                                            debug=debug, jobid=jobid,
                                            cpus=cpus)

        #GENOME HITS for md score bed files
        if md:
            print("\tMD:", file=sys.stderr)
            md_bedfile1 = get_center(bedfile=md_bedfile1, outname=md_bedfile1)
            bedtools_distance_keywords = dict(genomehits=genomehits, 
                                                ranked_center_file=md_bedfile1, 
                                                tempdir=tempdir, 
                                                distance_cutoff=largewindow)

            md_distances1 = multiprocess.main(function=bedtools_closest, 
                                            args=motif_list, 
                                            kwargs=bedtools_distance_keywords, 
                                            debug=debug, jobid=jobid,
                                            cpus=cpus)

            md_bedfile2 = get_center(bedfile=md_bedfile2, outname=md_bedfile2)
            bedtools_distance_keywords = dict(genomehits=genomehits, 
                                                ranked_center_file=md_bedfile2, 
                                                tempdir=tempdir, 
                                                distance_cutoff=largewindow)

            md_distances2 = multiprocess.main(function=bedtools_closest, 
                                            args=motif_list, 
                                            kwargs=bedtools_distance_keywords, 
                                            debug=debug, jobid=jobid,
                                            cpus=cpus)
            if use_config:
                config.vars['MD_DISTANCES1'] = md_distances1
                config.vars['MD_DISTANCES2'] = md_distances2
        if mdd:
            print("\tMDD:", file=sys.stderr)
            print(f'\t Completed: 0/{len(motif_distances)} ', end=' ', file=sys.stderr)
            mdd_bedfile1 = get_center(bedfile=mdd_bedfile1, outname=mdd_bedfile1)
            bedtools_distance_keywords = dict(genomehits=genomehits, 
                                                ranked_center_file=mdd_bedfile1, 
                                                tempdir=tempdir, 
                                                distance_cutoff=largewindow)

            mdd_distances1 = multiprocess.main(function=bedtools_closest, 
                                            args=motif_list, 
                                            kwargs=bedtools_distance_keywords, 
                                            debug=debug, jobid=jobid,
                                            cpus=cpus)

            mdd_bedfile2 = get_center(bedfile=mdd_bedfile2, outname=mdd_bedfile2)
            bedtools_distance_keywords = dict(genomehits=genomehits, 
                                                ranked_center_file=mdd_bedfile2, 
                                                tempdir=tempdir, 
                                                distance_cutoff=largewindow)

            mdd_distances2 = multiprocess.main(function=bedtools_closest, 
                                            args=motif_list, 
                                            kwargs=bedtools_distance_keywords, 
                                            debug=debug, jobid=jobid,
                                            cpus=cpus)
            # mdd_distances1 = []
            # mdd_distances2 = []
            # mdd_sorted_indices = np.argsort(pvals)
            # for i, single_motif_distances in enumerate(motif_distances, 1):
            #     motif = single_motif_distances[0]
            #     mdd_distances = single_motif_distances[1:]
            #     mdd_sorted_distances = [mdd_distances[i] for i in mdd_sorted_indices]
            #     if mdd_percent != False:
            #         cutoff = int(len(mdd_sorted_distances)*mdd_percent)
            #         mdd_distances2.append([motif] + mdd_sorted_distances[:cutoff])
            #         mdd_distances1.append([motif] + mdd_sorted_distances[cutoff:])
            #     else:
            #         sorted_pvals = [pvals[i] for i in mdd_sorted_indices]
            #         cutoff = int(len([p for p in sorted_pvals if p < mdd_pval]))
            #         mdd_distances2.append([motif] + mdd_sorted_distances[:cutoff])
            #         mdd_distances1.append([motif] + mdd_sorted_distances[cutoff:])
            #    # print(f'\r\t Completed: {i}/{len(motif_distances)} ', end=' ', flush=True, file=sys.stderr)
            if use_config:
                config.vars['MDD_DISTANCES1'] = mdd_distances1
                config.vars['MDD_DISTANCES2'] = mdd_distances2
    else:
        raise exceptions.InputError("SCANNER option not recognized.")

    if use_config:
        config.vars['MOTIF_DISTANCES'] = motif_distances

    total_time = time.time() - start_time
    if use_config:
        config.vars['SCANNERtime'] = total_time

    #Remove large fasta files from output folder
    # if fasta_file:
    #     fasta_file.unlink()
    # if md_fasta1:
    #     md_fasta1.unlink()
    # if md_fasta2:
    #     md_fasta2.unlink()
    # if mdd_fasta1:
    #     mdd_fasta1.unlink()
    # if mdd_fasta2:
    #     mdd_fasta2.unlink()

    print("done in: " + str(datetime.timedelta(seconds=int(total_time))), file=sys.stderr)

    if debug:
        multiprocess.current_mem_usage(jobid)

    return motif_distances, md_distances1, md_distances2, mdd_distances1, mdd_distances2

#Functions
#==============================================================================
def getfasta(bedfile=None, genomefasta=None, tempdir=None, outname=None):
    '''Converts a bed file to a fasta file using bedtools. Outputs into the 
        tempdir directory created by TFEA.

    Parameters
    ----------
    bedfile : string
        full path to a bed file

    genomefasta : string
        full path to a fasta file for the genome of interest
        
    tempdir : string
        full path to temp directory in output directory (created by TFEA)

    Returns
    -------
    ranked_file_fasta : string 
        full path to a fasta file containing the inputted bed file regions in 
        fasta format 
    '''
    fasta_file = tempdir / outname
    #pybedtools implementation (incomplete)
    # pybed = BedTool(bedfile).sequence(fi=genomefasta).saveas(fasta_file)

    getfasta_command = ["bedtools", "getfasta",
                        "-fi", genomefasta, 
                        "-bed", bedfile,
                        "-fo", fasta_file]
    
    try:
        subprocess.run(getfasta_command, check=True, 
                        stdout=subprocess.PIPE, 
                        stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        raise exceptions.SubprocessError(e.stderr.decode())

    return fasta_file

#==============================================================================
def fasta_linecount(fastafile=None):
    linecount = 0
    with open(fastafile) as F:
        for line in F:
            if line[0] == '>':
                linecount += 1
    
    return linecount

#==============================================================================
def fasta_names(fastafile=None):
    names = list()
    with open(fastafile) as F:
        for line in F:
            if line[0] == '>':
                names.append(line[1:].strip('\n'))
    
    return names
            
#==============================================================================
def fimo_background_file(window=None, tempdir=None, bedfile=None, 
                            genomefasta=None, order=None):
    background_bed_file = tempdir / "fimo_background.bed"
    outfile = open(background_bed_file,'w')
    with open(bedfile) as F:
        F.readline()
        for line in F:
            chrom, start, stop, name = line.strip('\n').split('\t')
            center = (int(start)+int(stop))/2
            start = int(center - window)
            stop = int(center + window)
            outfile.write('\t'.join([chrom, str(start), str(stop), name]) + '\n')
    outfile.close()
    
    getfasta(bedfile=background_bed_file, 
                                genomefasta=genomefasta, tempdir=tempdir, 
                                outname='fimo_background.fa')

    fasta_file = tempdir / 'fimo_background.fa'
    background_file = fasta_markov(tempdir=tempdir, fastafile=fasta_file, 
                                    order=order)
    
    return background_file

#==============================================================================
def fasta_markov(tempdir=None, fastafile=None, order=None):
    '''This function runs meme's fasta-get-markov function that generates a 
        background markov file (for use with fimo) from a fasta file.

    Parameters
    ----------
    tempdir : string
        full path to temp directory in output directory (created by TFEA)

    fastafile : string
        full path to fasta file that will be used to generate the markov
        background model file

    order : string
        an integer formatted as a string where a user may specify what order
        markov model they would like (default='0')
        
    Returns
    -------
    None
    '''
    markov_background = tempdir / "markov_background.txt"
    try:
        with open(markov_background, 'w') as output:
            subprocess.run(["fasta-get-markov", "-m", order, fastafile],
                            stdout=output, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        raise exceptions.SubprocessError(e.stderr.decode())
        
    return markov_background

#==============================================================================
def fimo(motif, bg_file=None, fasta_file=None, tempdir=None, 
        motifdatabase=None, thresh=None, largewindow=None):
    '''This function runs fimo on a given fastafile for a single motif in a 
        provided motif database. The output is cut and sorted to convert into 
        a sorted bed file

    Parameters
    ----------
    tempdir : string
        full path to temp directory in output directory (created by TFEA)

    motifdatabase : string
        full path to a motif database file in meme format

    bgfile : string
        full path to a markov background model

    motif : string
        the name of a motif that matches a motif within motifdatabase

    fastafile : string
        full path to a fasta file that fimo will perform motif scanning on
        
    Returns
    -------
    fimo_out : string
        full path to where fimo output which is stored within the tempdir 
        directory.
    '''
    fimo_out = tempdir / (motif+'.fimo.bed')
    if bg_file is not None:
        command = ("fimo", "--skip-matched-sequence", 
                    "--verbosity", "1", 
                    "--thresh", str(thresh),
                    "--bgfile", bg_file,
                    "--motif", motif, 
                    motifdatabase, fasta_file)
    else:
        command = ("fimo", "--skip-matched-sequence", 
                    "--verbosity", "1", 
                    "--thresh",str(thresh),
                    "--motif", motif, 
                    motifdatabase, fasta_file)

    try:
        fimo_out = subprocess.check_output(command, stderr=subprocess.PIPE).decode('UTF-8')
    except subprocess.CalledProcessError as e:
        raise exceptions.SubprocessError(e.stderr.decode())

    # fasta_count = fasta_linecount(fastafile=fasta_file)
    names = fasta_names(fastafile=fasta_file)
    distances = fimo_parse_stdout(fimo_stdout=fimo_out, 
                                    largewindow=largewindow, 
                                    names=names)
                                    # linecount=fasta_count)

    del fimo_out

    return [motif] + distances

#==============================================================================
def fimo_motif_names(motifdatabase=None):
    '''Extracts motif names from a MEME formatted motif database

    Parameters
    ----------
    motifdatabase : string
        full path to a meme formatted file with motifs to be used in TFEA
        
    Returns
    -------
    motif_list : list
        a list of motif names to be analyzed in TFEA
    '''
    motif_list = list()
    with open(motifdatabase) as F:
        for line in F:
            if 'MOTIF' in line:
                line = line.strip('\n').split()
                motif_name = line[-1]
                motif_list.append(motif_name)

    return motif_list

#==============================================================================
def fimo_parse(fimo_file=None, largewindow=None, retain='distance', 
                linecount=None):
    '''Parses a fimo output file and writes into a new file that is formatted
        in a way that can be parsed within existing TFEA functions

    Parameters
    ----------
    largewindow : float
        the size of the larger window to perform TFEA. Specified by user in
        config file

    tempdir : string
        full path to temp directory in output directory (created by TFEA)

    fimo_file : string
        full path to the fimo output file to be parsed by this function

    motif_file : string
        the name of the motif being parsed, this function will create a file
        using this motif_file string
        
    Returns
    -------
    outname : string
        the full path to the file to be used by other TFEA functions
    '''
    d = dict()
    with open(fimo_file) as F:
        header = F.readline().strip('\n').split('\t')
        if len(header) > 1:
            start_index = header.index('start')
            stop_index = header.index('stop')
            name_index = header.index('sequence_name')
            score_index = header.index('score')
            for line in F:
                line = line.strip('\n').split('\t')
                rank = line[name_index].split(',')[-1]
                rank = int(rank)
                start = line[start_index]
                stop = line[stop_index]
                distance = ((int(start)+int(stop))/2)-int(largewindow)
                score = line[score_index]
                if rank not in d:
                    d[rank] = [rank, score, distance]
                elif retain == 'score':
                    prev_score = float(d[rank][-2])
                    if prev_score < float(score):
                        d[rank] = [rank, score, distance]
                elif retain == 'distance':
                    prev_distance = float(d[rank][-1])
                    if prev_distance < float(distance):
                        d[rank] = [rank, score, distance]
    distances = list()
    for rank in range(1, linecount+1):
        if rank in d:
            distances.append(d[rank][-1])
        else:
            distances.append('.')

    return distances

#==============================================================================
def fimo_parse_stdout(fimo_stdout=None, largewindow=None, retain='score', 
                        names=None):
    '''
    '''
    d = dict()
    lines = fimo_stdout.split('\n')
    header = lines[0].split('\t')
    max_score = 0
    if len(header) > 1:
        start_index = header.index('start')
        stop_index = header.index('stop')
        name_index = header.index('sequence_name')
        score_index = header.index('score')
    for item in lines[1:-1]: #To remove header and empty last line in output
        line_list = item.split('\t')
        id = line_list[name_index]
        start = line_list[start_index]
        stop = line_list[stop_index]
        distance = ((int(start)+int(stop))/2)-int(largewindow)
        score = float(line_list[score_index])
        if score > max_score:
            max_score = score
        if id not in d:
            d[id] = [score, distance]
        elif retain == 'score':
            prev_score = d[id][0]
            if prev_score < score:
                d[id] = [score, distance]
        elif retain == 'distance':
            prev_distance = d[id][1]
            if prev_distance < distance:
                d[id] = [score, distance]
        elif retain == 'avg':
            d[id].append(score)
            d[id].append(distance)

    if retain == 'avg':
        for id in d:
            norm_scores = [d[id][i]/max_score for i in range(len(d[id])) if i%2 == 0]
            distances = [d[id][i] for i in range(len(d[id])) if i%2 != 0]
            avg_distance = [distances[i] for i in range(len(distances))]


    # distances = np.empty(len(names), dtype=object)
    distances = [] 
    for name in names:
        if name in d:
            distances.append(d[name][1])
        else:
            distances.append('.')

    return distances

#==============================================================================
def bedtools_closest(motif, genomehits=None, ranked_center_file=None, 
                        tempdir=None, distance_cutoff=None, rank_index=None):
    '''Calculates nearest motif hit from a bed file. TFEA provides this 
        function with a bed file containing the center of the inputted regions.

    Parameters
    ----------
    TFresults : list of lists
        contains calculated enrichment scores for all TFs of interest specified
        by the user
        
    Returns
    -------
    motif_distance_bed_sorted : string
        full path to where the sorted motif distance file was outputted
    '''
    try:
        motif_path = genomehits / motif
        if os.stat(motif_path).st_size == 0:
            return [motif] + ['.' for i in range(os.stat(ranked_center_file).st_size)]

        command = ("bedtools", "closest", "-D", "ref", "-t", "first", "-a", 
                    ranked_center_file, "-b", motif_path)
        closest_out = tempdir / (motif + '.closest.bed')

        # import sys
        # print(' '.join([str(c) for c in command]) + ' > ' + closest_out.as_posix(), file=sys.stderr)

        try:
            closest_out.write_bytes(subprocess.check_output(command, stderr=subprocess.PIPE))
        except subprocess.CalledProcessError as e:
            raise exceptions.SubprocessError(e.stderr.decode())
        
        
        distances = list()
        ranks = list()
        with open(closest_out) as F:
            for line in F:
                linelist = line.strip('\n').split('\t')
                distance = int(linelist[-1])
                if rank_index is not None:
                    rank = int(linelist[rank_index].split(',')[-1])
                    ranks.append(rank)
                if abs(distance) <= distance_cutoff:
                    distances.append(distance)
                else:
                    distances.append('.')
        if rank_index is not None:
            distances = [x for i,x in sorted(zip(ranks,distances))]
        closest_out.unlink()

    except Exception as e:
        # This prints the type, value, and stack trace of the
        # current exception being handled.
        print(traceback.print_exc())
        raise e

    return [motif.strip('.bed')] + distances

#==============================================================================
def get_center(bedfile=None, outname=None):
    '''This function takes a bed file and returns the center for all regions
        in the bed file.

    Parameters
    ----------
    bedfile : str
        full path to a bed file
    outname : str
        full path to the output centered bed file

    Returns
    -------
    outname : str
        full path to the output centered bed file
    '''
    outbed = BedTool(str(bedfile)).each(featurefuncs.center).sort()
    outbed.saveas(outname)

    return outname