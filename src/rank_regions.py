__author__ = 'Jonathan Rubin'

import math
import os
import subprocess
# from pybedtools import BedTool

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


def deseqfile(DESEQFILE,filedir):
    #Parse a deseq file and obtain the exact middle of each region (for motif distance calc later) and pvalue (to rank)
    ranked = list()
    up = list()
    down = list()
    with open(DESEQFILE) as F:
        F.readline()
        for line in F:
            line = line.strip('\n').split('\t')
            pval = format(float(line[-2]),'.12f')
            region = line[1].split(':')
            chrom = region[0]
            coordinates = region[1].split('-')
            start = coordinates[0]
            stop = coordinates[1]
            chrom = chrom.strip('"')
            stop = stop.strip('"')
            fc = float(line[4])
            if fc < 1:
                down.append((chrom,start,stop,pval,str(fc)))
            else:
                up.append((chrom,start,stop,pval,str(fc)))
            ranked.append((chrom,start,stop,pval))

    #Save ranked regions in a bed file (pvalue included)
    outfile = open(filedir + "ranked_file.bed",'w')
    # r=1
    # for region in sorted(ranked, key=lambda x: x[3]):
    #     outfile.write('\t'.join(region) + '\t' + str(r) + '\n')
    #     r += 1
    r=1
    for region in sorted(up, key=lambda x: x[3]):
        outfile.write('\t'.join(region) + '\t' + str(r) + '\n')
        r += 1
    for region in sorted(down, key=lambda x: x[3], reverse=True):
        outfile.write('\t'.join(region) + '\t' + str(r) + '\n')
        r += 1

    outfile.close()



