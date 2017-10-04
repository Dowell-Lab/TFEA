__author__ = 'Jonathan Rubin'

import math
import os

def log2fc(count_file,filedir,BAM1,BAM2):
    condition1 = float(len(BAM1))
    condition2 = float(len(BAM2))
    ranks = list()
    with open(count_file) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            chrom,start,stop = line[:3]
            cond1avg = sum([float(x) for x in line[3:3+condition1]])/condition1
            cond2avg = sum([float(y) for y in line[3+condition1:3+condition1+condition2]])/condition2
            log2fc = math.log(cond1avg/cond2avg,2)
            ranks.append((chrom,start,stop,str(log2fc)))


    outfile = open(filedir + "ranked_file.bed",'w')
    for region in sorted(ranks, key=lambda x: x[3], reverse=True):
        outfile.write('\t'.join(region[:-1]) + '\n')

    os.system("bedtools getfasta -fi " + genome + " -bed " + filedir + "ranked_file.bed > " + filedir + "ranked_file.fasta")