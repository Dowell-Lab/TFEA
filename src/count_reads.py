__author__ = 'Jonathan Rubin'

import os

def run(BED,BAM1,BAM2,LABEL1,LABEL2,filedir):
    #This os.system call runs bedtools multicov to count reads in all specified BAMs for given regions in BED
    os.system("bedtools multicov -bams " + " ".join(BAM1) + " " + " ".join(BAM2) + " -bed " + BED + " > " + filedir + "count_file.bed")

    #This section adds a header to the count_file and reformats it to remove excess information and add a column with the region for later use
    outfile = open(filedir + "count_file.header.bed",'w')
    outfile.write("chrom\tstart\tstop\tregion\t" + '\t'.join([LABEL1]*len(BAM1)) + "\t" + '\t'.join([LABEL2]*len(BAM2)) + "\n")
    with open(filedir + "count_file.bed") as F:
        for line in F:
            line = line.strip('\n').split('\t')
            chrom,start,stop = line[:3]
            counts = line[-len(BAM1)+len(BAM2):]
            outfile.write('\t'.join([chrom,start,stop]) + "\t"+chrom+":"+start+"-"+stop+"\t" + '\t'.join(counts) + "\n")
        

    # #This for loop goes through the count file to record the total number of mapped reads for each BAM file and saves all lines in line_storage to input into an outfile
    # total_reads = [0]*(len(BAM1)+len(BAM2))
    # line_storage = list()
    # with open(filedir + "count_file.bed") as F:
    #     for line in F:
    #         line = line.strip('\n').split('\t')
    #         chrom,start,stop = line[:3]
    #         counts = [float(x)+1.0 for x in line[3:]]
    #         total_reads = [a+b for a,b in zip(total_reads,counts)]
    #         line_storage.append((chrom,start,stop,counts))

    # #This loop goes through line_storage lines normalizing counts relative to the total counts for the first BAM file then writes all results to count_file.bed
    # try:
    #     normalization_factors = [total_reads[0]/a for a in total_reads]
    # except ZeroDivisionError:
    #     print "One of your BAM files has 0 reads mapping to your regions of interest"
    # outfile = open(filedir + "count_file.bed",'w')
    # for line in line_storage:
    #     chrom,start,stop,counts = line
    #     counts = [str(c*n) for c,n in zip(counts,normalization_factors)]
    #     outfile.write('\t'.join([chrom,start,stop]) + '\t' + '\t'.join(counts) + '\n')