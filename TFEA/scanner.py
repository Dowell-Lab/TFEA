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
import subprocess
import traceback
from pathlib import Path

import config
import fasta
import multiprocess
from exceptions import FileEmptyError, InputError, OutputError, SubprocessError

#Main Script
#==============================================================================
def main(fasta_file=None, md_fasta1=None, md_fasta2=None, ranked_file=None, 
        md_bedfile1=None, md_bedfile2=None, 
        scanner=config.SCANNER, md=config.MD, largewindow=config.LARGEWINDOW,
        smallwindow=config.SMALLWINDOW, genomehits=Path(config.GENOMEHITS), 
        fimo_background=config.FIMO_BACKGROUND, genomefasta=Path(config.GENOMEFASTA), 
        tempdir=Path(config.TEMPDIR), fimo_motifs=Path(config.FIMO_MOTIFS), 
        singlemotif=config.SINGLEMOTIF, fimo_thresh=config.FIMO_THRESH,
        debug=config.DEBUG):
    '''This is the main script of the SCANNER module. It returns motif distances
        to regions of interest by either scanning fasta files on the fly using
        fimo or homer or by using bedtools closest on a center bed file and 
        a database of bed files corresponding to motif hits across the genome

    Parameters
    ----------
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
    #FIMO
    if scanner == 'fimo':
        #Get background file, if none desired set to 'None'
        if fimo_background == 'largewindow':
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
        fimo_keywords = dict(bg_file=background_file, fasta_file=fasta_file, 
                            tempdir=tempdir, motifdatabase=fimo_motifs, 
                            thresh=fimo_thresh, 
                            largewindow=largewindow)

        motif_distances = multiprocess.main(function=fimo, args=motif_list, 
                                            kwargs=fimo_keywords, debug=debug)

        #FIMO for md score fasta files
        if md:
            fimo_keywords = dict(bg_file=background_file, fasta_file=md_fasta1, 
                            tempdir=tempdir, motifdatabase=fimo_motifs, 
                            thresh=fimo_thresh, 
                            largewindow=largewindow)
            md_distances1 = multiprocess.main(function=fimo, args=motif_list, 
                                                kwargs=fimo_keywords, 
                                                debug=debug)
            
            fimo_keywords = dict(bg_file=background_file, fasta_file=md_fasta2, 
                            tempdir=tempdir, motifdatabase=fimo_motifs, 
                            thresh=fimo_thresh, 
                            largewindow=largewindow)
            md_distances2 = multiprocess.main(function=fimo, args=motif_list, 
                                                kwargs=fimo_keywords, 
                                                debug=debug)

            return motif_distances, md_distances1, md_distances2
            
    #HOMER
    elif scanner== 'homer':
        raise InputError("Homer scanning is not supported at this time.")

    #GENOME HITS
    elif scanner == 'genome hits':
        #Get motifs to analyze
        if singlemotif == False:
            motif_list = os.listdir(genomehits)
        else:
            motif_list = [os.path.join(genomehits, motif) for motif in singlemotif.split(',')]

        #Perform bedtools closest to get distances
        bedtools_distance_keywords = dict(genomehits=genomehits, 
                                            ranked_center_file=ranked_file, 
                                            tempdir=tempdir, 
                                            distance_cutoff=largewindow)
        
        motif_distances = multiprocess.main(function=bedtools_closest, 
                                            args=motif_list, 
                                            kwargs=bedtools_distance_keywords, 
                                            debug=debug)

        #GENOME HITS for md score bed files
        if md:
            bedtools_distance_keywords = dict(genomehits=genomehits, 
                                                ranked_center_file=md_bedfile1, 
                                                tempdir=tempdir, 
                                                distance_cutoff=largewindow)

            md_distances1 = multiprocess.main(function=bedtools_closest, 
                                            args=motif_list, 
                                            kwargs=bedtools_distance_keywords, 
                                            debug=debug)
            
            bedtools_distance_keywords = dict(genomehits=genomehits, 
                                                ranked_center_file=md_bedfile2, 
                                                tempdir=tempdir, 
                                                distance_cutoff=largewindow)

            md_distances2 = multiprocess.main(function=bedtools_closest, 
                                            args=motif_list, 
                                            kwargs=bedtools_distance_keywords, 
                                            debug=debug)

            return motif_distances, md_distances1, md_distances2
    else:
        raise InputError("SCANNER option not recognized.")
            
    return motif_distances

#Functions
#==============================================================================
def fimo_background_file(window=None, tempdir=None, bedfile=None, 
                            genomefasta=None, order=None):
    background_bed_file = tempdir / "fimo_background.bed"
    outfile = open(background_bed_file,'w')
    with open(bedfile) as F:
        for line in F:
            chrom, start, stop, name = line.strip('\n').split('\t')
            center = (int(start)+int(stop))/2
            start = int(center - window)
            stop = int(center + window)
            outfile.write('\t'.join([chrom, str(start), str(stop), name]) + '\n')
    outfile.close()
    
    fasta.getfasta(bedfile=background_bed_file, 
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
        raise SubprocessError(e.stderr.decode())
        
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
    if bg_file != None:
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
        raise SubprocessError(e.stderr.decode())

    fasta_linecount = fasta.fasta_linecount(fastafile=fasta_file)
    # distances = fimo_parse(fimo_file=fimo_out, largewindow=largewindow, 
    #                         linecount=fasta_linecount)
    distances = fimo_parse_stdout(fimo_stdout=fimo_out, 
                                    largewindow=largewindow, 
                                    linecount=fasta_linecount)

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
        motif_path = genomehits / motif
        if os.stat(motif_path).st_size == 0:
            return [motif] + ['.' for i in range(os.stat(ranked_center_file).st_size)]
        # bedtools_closest_out = os.path.join(tempdir, motif + '.closestBed.bed')

        # command = ("bedtools closest -D ref -t first -a " 
        #             + ranked_center_file + " -b " + motif_path + " > "
        #             + bedtools_closest_out)

        command = ("bedtools", "closest", "-D", "ref", "-t", "first", "-a", 
                    ranked_center_file, "-b", motif_path)
        closest_out = tempdir / (motif + 'closest.bed')
        try:
            # closest_out = subprocess.check_output(command, stderr=subprocess.PIPE).decode('UTF-8')
            closest_out.write_bytes(subprocess.check_output(command, stderr=subprocess.PIPE))
        except subprocess.CalledProcessError as e:
            raise SubprocessError(e.stderr.decode())
        
        
        distances = list()
        with open(closest_out) as F:
            for line in F:
        # for line in closest_out.split('\n')[1:-1]:
                distance = int(line.strip('\n').split('\t')[-1])
                if distance <= distance_cutoff:
                    distances.append(distance)
                else:
                    distances.append('.')

        closest_out.unlink()

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