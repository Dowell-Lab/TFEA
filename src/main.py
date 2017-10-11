__author__ = 'Jonathan Rubin'

import os
from config import BED,BAM1,BAM2,SINGLEMOTIF,DATABASE,GENOME,MEMEDB,DESEQFILE
import count_reads
import rank_regions
import motif_distance
import ES_calculator

def run():
    #Booleans that allow you to run specific parts of this python package
    count = True
    rank = True
    distance = True
    calculate = True

    #Choose what type of ranking metric to be used to rank regions of interest:
    #Options: "log2fc", "deseqfile"
    rank_metric = "deseqfile"

    #Home directory
    homedir = os.path.dirname(os.path.realpath(__file__))

    #Directory where all temp files will be stored
    filedir = parent_dir(homedir) + '/files/'

    #Path to count file. Can be changed if using your own count file. Generated in count_reads module
    count_file = filedir + "count_file.bed"

    #Path to ranked fasta file. Generated in rank_regions module
    ranked_file = filedir + "ranked_file.fasta"

    #This module counts reads from all Bam files in BAM1 and BAM2 and creates count_file with this info.
    if count:
        print "Counting reads in regions..."
        count_reads.run(BED,BAM1,BAM2,filedir)
        print "done"

    #This module ranks regions in BED based on some specified metric
    if rank:
        print "Ranking regions based on " + rank_metric + "..."
        if rank_metric == "log2fc":
            rank_regions.log2fc(count_file,filedir,GENOME,BAM1,BAM2)
        if rank_metric == "deseqfile":
            rank_regions.deseqfile(DESEQFILE,filedir)
        print "done"

    #Scans ranked BED regions for motifs of interest and records them in distance file
    if distance:
        print "Finding motif hits in regions..."
        motif_distance.run(ranked_file,filedir,MEMEDB,DATABASE,SINGLEMOTIF)
        print "done"

    #Calculates an Enrichment Score and a Normalized Enrichment Score for all specified motifs
    if calculate:
        ES_calculator.run()


#Return parent directory
def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])
    
    return newdir