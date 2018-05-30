__author__ = 'Jonathan Rubin'

import os
from config import BAM1, BAM2, LABEL1, LABEL2, FILEDIR

def run(BED):
    #This os.system call runs bedtools multicov to count reads in all specified BAMs for given regions in BED
    os.system("bedtools multicov -bams " + " ".join(BAM1) + " " + " ".join(BAM2) + " -bed " + BED + " > " + FILEDIR + "count_file.bed")

    #This section adds a header to the count_file and reformats it to remove excess information and add a column with the region for later use
    outfile = open(FILEDIR + "count_file.header.bed",'w')
    outfile.write("chrom\tstart\tstop\tregion\t" + '\t'.join([LABEL1]*len(BAM1)) + "\t" + '\t'.join([LABEL2]*len(BAM2)) + "\n")
    with open(FILEDIR + "count_file.bed") as F:
        for line in F:
            line = line.strip('\n').split('\t')
            chrom,start,stop = line[:3]
            counts = line[-(len(BAM1)+len(BAM2)):]
            outfile.write('\t'.join([chrom,start,stop]) + "\t"+chrom+":"+start+"-"+stop+"\t" + '\t'.join(counts) + "\n")