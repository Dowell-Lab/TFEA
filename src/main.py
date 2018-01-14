__author__ = 'Jonathan Rubin'

import os
import time
from config import COMBINE,COUNT,RANK,DISTANCE,CALCULATE,BEDS,BAM1,BAM2,SINGLEMOTIF,DATABASE,GENOME,MEMEDB,MOTIF_HITS,DESEQFILE
import combine_bed
import count_reads
import rank_regions
import DESeq
import motif_distance
import ES_calculator
import create_html

def run():
    #Choose what type of ranking metric to be used to rank regions of interest:
    #Options: "log2fc", "deseqfile"
    rank_metric = "deseqfile"

    #Home directory
    homedir = os.path.dirname(os.path.realpath(__file__))

    #Directory where all temp files will be stored
    filedir = parent_dir(homedir) + '/files/'

    figuredir = parent_dir(homedir) + '/output/plots/'

    #Output directory
    output = parent_dir(homedir) + '/output/'

    #Path to count file. Can be changed if using your own count file. Generated in count_reads module
    count_file = filedir + "count_file.bed"

    #Path to ranked file. Can be changed if using your own ranked file. Generated in rank_regions module
    ranked_file = filedir + "ranked_file.bed"

    ranked_center_distance_file = filedir + "ranked_file.center.sorted.distance.bed"

    #This module takes the input list of BED files, concatenates them, and then merges them (bedtools).
    if COMBINE:
        BED = combine_bed.run(BEDS,filedir)
    else:
        BED = BEDS[0]

    #This module counts reads from all Bam files in BAM1 and BAM2 and creates count_file with this info.
    if COUNT:
        print "Counting reads in regions..."
        count_reads.run(BED,BAM1,BAM2,filedir)
        print "done"

    #This module ranks regions in BED based on some specified metric
    if RANK:
        print "Ranking regions based on " + rank_metric + "..."
        if rank_metric == "log2fc":
            rank_regions.log2fc(count_file,filedir,GENOME,BAM1,BAM2)
        if rank_metric == "deseqfile":
            rank_regions.deseqfile(DESEQFILE,GENOME,filedir,MOTIF_HITS,SINGLEMOTIF)
        if rank_metric == "deseq":
            DESEQFILE2 = DESeq.run(LABEL1,LABEL2,BAM1,BAM2,scriptdir,filedir)
            rank_regions.deseqfile(DESEQFILE2,GENOME,filedir,MOTIF_HITS,SINGLEMOTIF)
        print "done"

    #Scans ranked BED regions for motifs of interest and records them in distance file
    if DISTANCE:
        fimo = False
        if fimo:
            print "Warning: This part of this package is incomplete"
            motif_distance.runfimo(ranked_file,filedir,MEMEDB,DATABASE)
        else:
            print "Finding motif hits in regions..."
            if SINGLEMOTIF == False:
                TFresults = list()
                NESlist = list()
                for MOTIF_FILE in os.listdir(MOTIF_HITS):
                    a = time.time()
                    total_hits = motif_distance.run(ranked_file,filedir,MOTIF_HITS+MOTIF_FILE)
                    print "Motif distances done in: ",time.time()-a,"s"

                    #This module is where the bulk of the analysis is done. The functions below calculate ES,NES,p-value,FDR for each TF motif in
                    #the HOCOMOCO database.
                    if CALCULATE:
                        a = time.time()
                        results = ES_calculator.run(MOTIF_FILE,ranked_center_distance_file,figuredir,filedir,total_hits)
                        TFresults.append(results)
                        NESlist.append(results[1])
                        print results[0],results[1],results[2],results[3]
                        print "ES calculation done in: ", time.time()-a,"s"
                TFresults = ES_calculator.FDR(TFresults,NESlist)
                outfile = open(output + 'results.txt', 'w')
                outfile.write('TF-Motif\tES\tNES\tP-value\tFDR\n')
                for val in sorted(TFresults, key=lambda x: x[3]):
                    outfile.write('\t'.join([str(val[i]) for i in range(len(val)) if i!=4]) +  '\n')
                create_html.run(TFresults,output)

            #Note if you set the SINGLEMOTIF variable to a specific TF, this program will be unable to accurately determine a NES or FDR for the given motif.
            #This will return only
            else:
                total_hits = motif_distance.run(ranked_file,filedir,MOTIF_HITS+SINGLEMOTIF)
                results = ES_calculator.run(ranked_center_distance_file,figuredir,filedir,total_hits)
                print SINGLEMOTIF, results[0], results[1], results[2], results[3]
        print "done"


#Return parent directory
def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])
    
    return newdir