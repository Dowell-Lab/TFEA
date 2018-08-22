__author__ = 'Jonathan Rubin'

import os
import math
import combine_bed, config

#Takes a single ranked_file and returns bed file with the centers for each region and a bed file with the distance to the closest motif for each centered region
def run(ranked_center_file,MOTIF_PATH):
    #Get closest motif hit to center base of each region
    # os.system("bedtools closest -d -t first -a " + ranked_center_file.split('.bed')[0] + ".sorted.bed -b " + MOTIF_PATH + " > " + ranked_center_file.split('.bed')[0] + ".sorted.distance.bed")
    os.system("bedtools closest -D ref -t first -a " + ranked_center_file.split('.bed')[0] + ".sorted.bed -b " + MOTIF_PATH + " > " + '/' + '/'.join(ranked_center_file.split('/')[:-1]) + '/' + MOTIF_PATH.split('/')[-1] + ".sorted.distance.bed")

    return '/' + '/'.join(ranked_center_file.split('/')[:-1]) + '/' + MOTIF_PATH.split('/')[-1] + ".sorted.distance.bed"

def fimo_distance(fimo_file,MOTIF_FILE):
    d = dict()
    with open(fimo_file) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            pval,fc,rank = line[0].split(',')
            start,stop = line[1:2]
            distance = ((int(start)+int(stop))/2)-config.LARGEWINDOW
            score = line[5]
            if rank not in d:
                d[rank] = [rank,pval,fc,score,distance]
            else:
                prev_distance = math.fabs(d[rank][-1])
                if prev_distance > math.fabs(distance):
                    d[rank] = [rank,pval,fc,score,distance]

    outname = config.FILEDIR+MOTIF_FILE+'.sorted.distance.bed'
    outfile = open(outname,'w')
    for i in range(len(d)):
        rank = str(i)
        outfile.write('\t'.join(d[rank])+'\n')
    outfile.close()

    return outname

#Function to calculate GC content for all regions - creates panel underneath enrichment score figure in html output
def get_gc(ranked_file,window=config.SMALLWINDOW):

    '''This function calculates gc content over all eRNAs. It uses the SMALLWINDOW variable within
    the config file instead of the whole region. It performs a running average of window size = total_regions/1000

    Parameters                                                                                                                                                                           
    ----------                                                                                                                                                                           
    ranked_bed_file : tab delimited bedfile with regions ranked by some metric of 
        differential transcription                                                                                                                                               
                                                                                                                                       
    Returns                                                                                                                                                                             
    -------                                                                                                                                                                             
    AUC : float                                                                                                                                                                   
        the AUC for a given tf

    p-value : float                                                                                                                                         
                 
    Raises
    ------
    ValueError
        when tf does not have any hits

    '''

    outfile = open(config.FILEDIR+"ranked_file.windowed.bed",'w')
    with open(config.FILEDIR + "ranked_file.bed") as F:
        for line in F:
            line = line.strip('\n').split('\t')
            chrom,start,stop = line[:3]
            center = (int(start)+int(stop))/2
            newstart = center - window
            newstop = center + window
            outfile.write(chrom + '\t' + str(newstart) + '\t' + str(newstop) + '\t' + '\t'.join(line[3:]) + '\n')
    outfile.close()

    ranked_file_windowed_fasta = combine_bed.get_fasta(config.FILEDIR+"ranked_file.windowed.bed")

    gc_array = []
    with open(ranked_file_windowed_fasta) as F:
        for line in F:
            if '>' not in line:
                line = line.strip('\n')
                gc_array.append(convert_sequence_to_array(line))



def convert_sequence_to_array(sequence):

    '''This function converts a DNA sequence (ACGT alphabet) to an array, collapsing GCs
    to 1's and ATs to 0's

    Parameters                                                                                                                                                                           
    ----------                                                                                                                                                                           
    sequence : a string containing ACGT character                                                                                                                                               
                                                                                                                                       
    Returns                                                                                                                                                                             
    -------                                                                                                                                                                             
    array : list
        of floats corresponding to 1 or 0 depending on whether the 
        input sequence was GC or AT at a specific site       

    Raises
    ------
    ValueError
        when a character is not ACTG                                                                                                                                                  

    '''
    array = []
    for character in sequence:
        character = character.upper()
        if character == 'G' or character == 'C':
            array.append(1.0)
        elif character == 'A' or character == 'T':
            array.append(0.0)
        else:
            print "Warning: Character not recognized"
            array.append(0.0)

    return array 







