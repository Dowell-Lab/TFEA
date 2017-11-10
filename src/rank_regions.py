__author__ = 'Jonathan Rubin'

import math
import os
import subprocess
from pybedtools import BedTool

def log2fc(count_file,filedir,GENOME,BAM1,BAM2):
    #This loop calculates log2fc based on the average normalized counts for all replicates and appends that and chrom,start,stop to a list called ranks
    condition1 = len(BAM1)
    condition2 = len(BAM2)
    ranks = list()
    with open(count_file) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            chrom,start,stop = line[:3]
            cond1avg = sum([float(x) for x in line[3:3+condition1]])/float(condition1)
            cond2avg = sum([float(y) for y in line[3+condition1:3+condition1+condition2]])/float(condition2)
            log2fc = math.log(cond1avg/cond2avg,2)
            ranks.append((chrom,start,stop,str(log2fc)))

    #Creates an out file called ranked_file.bed which just has chrom,start,stop for all regions sorted by log2fc
    outfile = open(filedir + "ranked_file.bed",'w')
    for region in sorted(ranks, key=lambda x: x[3], reverse=True):
        outfile.write('\t'.join(region) + '\n')


def deseqfile_up_down(DESEQFILE,GENOME,filedir,MOTIF_HITS,SINGLEMOTIF):
    rankup = list()
    rankdown = list()
    with open(DESEQFILE) as F:
        F.readline()
        for line in F:
            try:
                line = line.strip('\n').split('\t')
                pval = format(float(line[-2]),'.12f')
                chrom,start,stop = line[1].split(',')
                chrom = chrom.strip('"')
                stop = stop.strip('"')
                center = (int(start)+int(stop))/2
                fc = float(line[5])
                if fc > 1:
                    rankup.append((chrom,str(center),str(center+1),pval))
                else:
                    rankdown.append((chrom,str(center),str(center+1),pval))
            except ValueError:
                pass


    if SINGLEMOTIF == False:
        print "This module isn't finished yet"
    else:
        for MOTIF_FILE in os.listdir(MOTIF_HITS):
            if MOTIF_FILE == SINGLEMOTIF:
                a = BedTool(rankup)
                b = BedTool(MOTIF_HITS + MOTIF_FILE)
                c = a.closest(b,d=True)
                d = BedTool(rankdown)
                e = d.closest(b,d=True)



    outfile = open(filedir + "ranked_file_up.bed",'w')
    for region in sorted(c, key=lambda x: x[3]):
        outfile.write('\t'.join(region[:4]) + '\t' + region[-1] + '\n')

    outfile = open(filedir + "ranked_file_down.bed",'w')
    for region in sorted(e, key=lambda x: x[3]):
        outfile.write('\t'.join(region[:4]) + '\t' + region[-1] + '\n')


def deseqfile(DESEQFILE,GENOME,filedir,MOTIF_HITS,SINGLEMOTIF):
    #Parse a deseq file and obtain the exact middle of each region (for motif distance calc later) and pvalue (to rank)
    ranked = list()
    with open(DESEQFILE) as F:
        F.readline()
        for line in F:
            try:
                line = line.strip('\n').split('\t')
                pval = format(float(line[-2]),'.12f')
                chrom,start,stop = line[1].split(',')
                chrom = chrom.strip('"')
                stop = stop.strip('"')
                ranked.append((chrom,start,stop,pval))
            except ValueError:
                pass

    #Save ranked regions in a bed file (pvalue included)
    outfile = open(filedir + "ranked_file.bed",'w')
    for region in sorted(c, key=lambda x: x[3]):
        outfile.write('\t'.join(region) + '\n')



