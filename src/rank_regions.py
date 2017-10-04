__author__ = 'Jonathan Rubin'

import math
import os

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
        outfile.write('\t'.join(region[:-1]) + '\n')

    #This os.system call uses bedtools to convert the ranked_file.bed into fasta format (ranked_file.fasta)
    # os.system("bedtools getfasta -fi " + GENOME + " -bed " + filedir + "ranked_file.bed -fo " + filedir + "ranked_file.fasta")
    os.system("bedtools getfasta -fi /Users/joru1876/scratch_backup/hg19_reference_files/hg19_all.fa -bed /Users/joru1876/scratch_backup/TFEA/files/ranked_file.bed -fo /Users/joru1876/scratch_backup/TFEA/files/ranked_file.fasta")
    print "bedtools getfasta -fi " + GENOME + " -bed " + filedir + "ranked_file.bed -fo " + filedir + "ranked_file.fasta"