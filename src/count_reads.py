__author__ = 'Jonathan Rubin'

import os

def run(BED,BAM1,BAM2,filedir):
    os.system("bedtools multicov -bams " + " ".join(BAM1) + " " + " ".join(BAM2) + " -bed " + BED + " > " + filedir + "count_file.bed")
    read_totals = [0]*(len(BAM1)+len(BAM2))
    line_storage = list()
    with open(filder + "count_file.bed") as F:
        for line in F:
            line = line.strip('\n').split('\t')
            chrom,start,stop = line[:3]
            counts = [float(x)+1.0 for x in line[3:]]
            read_totals = [a+b for a,b in zip(read_totals,counts)]
            line_storage.append((chrom,start,stop,counts))

    normalization_factors = [total_reads[0]/a for a in total_reads]
    outfile = open("count_file.bed",'w')
    for line in line_storage:
        chrom,start,stop,counts = line
        counts = [str(c*n) for c,n in zip(counts,normalization_factors)]
        outfile.write('\t'.join([chrom,start,stop]) + '\t' + '\t'.join(counts) + '\n')