#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''This file contains a list of functions associated with outputting a ranked
    bed file based on some metric.
'''

#==============================================================================
__author__ = 'Jonathan D. Rubin and Rutendo F. Sigauke'
__credits__ = ['Jonathan D. Rubin', 'Rutendo F. Sigauke', 'Jacob T. Stanley',
                'Robin D. Dowell']
__maintainer__ = 'Jonathan D. Rubin'
__email__ = 'Jonathan.Rubin@colorado.edu'

#Imports
#==============================================================================
import os
import subprocess
import traceback

import fasta_functions


#Functions
#==============================================================================
def fimo_background_file(window=None, tempdir=None, bedfile=None, 
                            genomefasta=None, order=None):
    background_bed_file = os.path.join(tempdir, "fimo_background.bed")
    outfile = open(background_bed_file,'w')
    with open(bedfile) as F:
        for line in F:
            chrom, start, stop, name = line.strip('\n').split('\t')
            center = (int(start)+int(stop))/2
            start = int(center - window)
            stop = int(center + window)
            outfile.write('\t'.join([chrom, str(start), str(stop), name]) + '\n')
    outfile.close()
    
    fasta_functions.getfasta(bedfile=background_bed_file, 
                                genomefasta=genomefasta, tempdir=tempdir, 
                                outname='fimo_background.fa')

    fasta_file = os.path.join(tempdir, 'fimo_background.fa')
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
    markov_background = os.path.join(tempdir, "markov_background.txt")
    os.system("fasta-get-markov -m " + order + " " + fastafile +
              " > " + markov_background)

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
    import fasta_functions
    fimo_out = os.path.join(tempdir, motif + '.fimo.bed')
    if bg_file != None:
        command = ("fimo", "--skip-matched-sequence", 
                    "--verbosity", "1", 
                    "--thresh",str(thresh),
                    "--bgfile", bg_file,
                    "--motif", motif, 
                    motifdatabase, fasta_file)
    else:
        # command = ("fimo --skip-matched-sequence --verbosity 1 --thresh " 
        #             + str(thresh) 
        #             + " --motif " + motif + " " + motifdatabase 
        #             + " " + fasta_file
        #             + " > " + fimo_out)
        command = ("fimo", "--skip-matched-sequence", 
                    "--verbosity", "1", 
                    "--thresh",str(thresh),
                    "--motif", motif, 
                    motifdatabase, fasta_file)
    with open(os.devnull, 'w') as devnull:
        # subprocess.call(command, shell=True, stderr=devnull)
        fimo_out = subprocess.check_output(command, stderr=devnull).decode('UTF-8')

    fasta_linecount = fasta_functions.fasta_linecount(fastafile=fasta_file)
    # distances = fimo_parse(fimo_file=fimo_out, largewindow=largewindow, 
    #                         linecount=fasta_linecount)
    distances = fimo_parse_stdout(fimo_stdout=fimo_out, 
                                    largewindow=largewindow, 
                                    linecount=fasta_linecount)

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
def fimo_parse(fimo_file=None, largewindow=None, retain='score', 
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
                        linecount=None):
    '''
    '''
    d = dict()
    lines = fimo_stdout.split('\n')
    header = lines[0].split('\t')
    if len(header) > 1:
        start_index = header.index('start')
        stop_index = header.index('stop')
        name_index = header.index('sequence_name')
        score_index = header.index('score')
    for item in lines[1:-1]: #To remove header and empty last line in output
        line_list = item.split('\t')
        rank = line_list[name_index].split(',')[-1]
        rank = int(float(rank))
        start = line_list[start_index]
        stop = line_list[stop_index]
        distance = ((int(start)+int(stop))/2)-int(largewindow)
        score = line_list[score_index]

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
def bedtools_closest(motif, genomehits=None, ranked_center_file=None, 
                        tempdir=None, distance_cutoff=None):
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
        motif_path = os.path.join(genomehits, motif)
        if os.stat(motif_path).st_size == 0:
            return [motif] + ['.' for i in range(os.stat(ranked_center_file).st_size)]
        # bedtools_closest_out = os.path.join(tempdir, motif + '.closestBed.bed')

        # command = ("bedtools closest -D ref -t first -a " 
        #             + ranked_center_file + " -b " + motif_path + " > "
        #             + bedtools_closest_out)

        command = ("bedtools closest -D ref -t first -a " 
                    + ranked_center_file + " -b " + motif_path)

        closest_out = subprocess.check_output(command, shell=True).decode('UTF-8')
        
        distances = list()
        # with open(bedtools_closest_out) as F:
        #     for line in F:
        for line in closest_out.split('\n')[1:-1]:
            distance = int(line.strip('\n').split('\t')[-1])
            if distance <= distance_cutoff:
                distances.append(distance)
            else:
                distances.append('.')

    except Exception as e:
        # This prints the type, value, and stack trace of the
        # current exception being handled.
        print(traceback.print_exc())
        raise e

    return [motif] + distances

#TESTS
#==============================================================================
if __name__ == "__main__":
    import scanner_functions
    motif='P53_HUMAN.H11MO.0.A'
    fasta_file = './test/random_10.fa'
    motifdatabase = './motif_databases/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme'
    tempdir='./test'
    thresh=0.0001
    fimo_out = os.path.join(tempdir, motif + '.fimo.bed')

    results = scanner_functions.fimo(motif, fasta_file=fasta_file, 
                                    tempdir=tempdir, motifdatabase=motifdatabase, 
                                    thresh=thresh, largewindow=1500)

    print(results)