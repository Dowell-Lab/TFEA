__author__ = 'Jonathan Rubin'

import os
from config import BED,BAM1,BAM2,SINGLEMOTIF,DATABASE,GENOME,MEMEDB,MOTIF_HITS,DESEQFILE
import count_reads
import rank_regions
import motif_distance
import ES_calculator

def run():
    #Booleans that allow you to run specific parts of this python package
    count = False
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

    figuredir = parent_dir(homedir) + '/figures/'

    #Path to count file. Can be changed if using your own count file. Generated in count_reads module
    count_file = filedir + "count_file.bed"

    #Path to ranked file. Can be changed if using your own ranked file. Generated in rank_regions module
    ranked_file = filedir + "ranked_file.bed"

    ranked_center_distance_file = filedir + "ranked_file.center.sorted.distance.bed"

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
            rank_regions.deseqfile(DESEQFILE,GENOME,filedir,MOTIF_HITS,SINGLEMOTIF)
        print "done"

    #Scans ranked BED regions for motifs of interest and records them in distance file
    if distance:
        fimo = False
        if fimo:
            print "Warning: This part of this package is incomplete"
            motif_distance.runfimo(ranked_file,filedir,MEMEDB,DATABASE)
        else:
            print "Finding motif hits in regions..."
            if SINGLEMOTIF == False:
                outfile1 = open(filedir + 'results.tmp', 'w')
                outfile1.write('TF-Motif\tES\tNES\tp-value\n')
                ES = list()
                for MOTIF_FILE in os.listdir(MOTIF_HITS):
                    total_hits = motif_distance.run(ranked_file,filedir,MOTIF_HITS+MOTIF_FILE)
                    if calculate:
                        results = (MOTIF_FILE,ES_calculator.run(ranked_center_distance_file,figuredir,filedir,total_hits))
                        outfile1.write('\t'.join([str(val) for val in results]) +  '\n')
                        ES.append(results)
                outfile = open(filedir + 'results.txt', 'w')
                outfile.write('TF-Motif\tES\tNES\tp-value\n')
                for val in sorted(ES, key=lambda x: x[3]):
                    outfile.write('\t'.join([str(i) for i in val]) +  '\n')

            else:
                total_hits = motif_distance.run(ranked_file,filedir,MOTIF_HITS+SINGLEMOTIF)
                ES_calculator.run_deseq(ranked_center_distance_file,figuredir)
        print "done"


#Return parent directory
def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])
    
    return newdir