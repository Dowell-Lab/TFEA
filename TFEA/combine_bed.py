__author__ = 'Jonathan Rubin'

import os
from config import LARGEWINDOW

def run(BEDS,filedir):
    os.system("cat " + " ".join(BEDS) + " > " + filedir + "combined_input.bed")
    os.system("sort -k1,1 -k2,2n " + filedir + "combined_input.bed > " + filedir + "combined_input.sorted.bed")
    os.system("bedtools merge -i " + filedir + "combined_input.sorted.bed > " + filedir + "combined_input.merge.bed")
    #os.system("bedtools intersect -v -a " + filedir + "combined_input.merge1.bed -b /scratch/Users/rusi2317/projects/rotation/bin/TFEA/files/hg19_TSS_sorted.bed > " + filedir + "combined_input.merge.bed") ##to remove TSS regions
    ## os.system("rm " + filedir + "combined_input.bed")
    return filedir + "combined_input.merge.bed"

#5/23/18: This function returns a fasta file from a bed input
def getfasta(genomefasta,bedfile,filedir):
    os.system("bedtools getfasta -name -fi "+genomefasta+" -bed "+bedfile+" > combined_input.merge.fa")

    return filedir + "combined_input.merge.fa"

#5/23/18: This function outputs a zero order markov background model file (is text file right?) for use with fimo later on
def get_bgfile(fastafile,filedir):
    os.system("fasta-get-markov "+fastafile+" "+filedir+"markov_background.txt")

    return filedir + "markov_background.txt"

def get_regions(ranked_center_file,filedir):
    outfile = open(filedir+'ranked_file.fullregions.bed','w')
    with open(ranked_center_file) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            chrom,start,stop = line[:3]
            start = str(int(start)-LARGEWINDOW)
            stop = str(int(stop)+LARGEWINDOW)
            pval,fc,rank = line[3:]
            name = ','.join([pval,fc,rank])
            outfile.write('\t'.join([chrom,start,stop,name]) + '\n')

    return filedir + 'ranked_file.fullregions.bed'




##def run(BEDS,filedir):
##    os.system("sort -k1,1 -k2,2n " + BEDS[0] + " > " + filedir + "bed1_sort.bed")
##    os.system("sort -k1,1 -k2,2n " + BEDS[1] + " > " + filedir + "bed2_sort.bed")
##    os.system("bedtools intersect -a " + filedir + "bed1_sort.bed -b " + filedir + "bed2_sort.bed > " + filedir + "combined_input.bed")
##    os.system("sort -k1,1 -k2,2n " + filedir + "combined_input.bed > " + filedir + "combined_input_sorted.bed")
##    os.system("bedtools merge -i " + filedir + "combined_input_sorted.bed > " + filedir + "combined_input.merge.bed")
    
##    return filedir + "combined_input.merge.bed"


#def run(BEDS,filedir):
#    os.system("sort -k1,1 -k2,2n " + BEDS[0] + " > " + filedir + "bed1_sort.bed")
#    os.system("sort -k1,1 -k2,2n " + BEDS[1] + " > " + filedir + "bed2_sort.bed")
#    os.system("bedtools intersect -v -a " + filedir + "bed1_sort.bed -b " + filedir + "bed2_sort.bed > " + filedir + "combined_inputa.bed")
#    os.system("bedtools intersect -v -b " + filedir + "bed1_sort.bed -a " + filedir + "bed2_sort.bed > " + filedir + "combined_inputb.bed")
#    os.system("cat " + filedir + "combined_inputa.bed " + filedir + "combined_inputb.bed > " + filedir + "combined_input.bed")
#    os.system("sort -k1,1 -k2,2n " + filedir + "combined_input.bed > " + filedir + "combined_input_sorted.bed")
#    os.system("bedtools merge -i " + filedir + "combined_input_sorted.bed > " + filedir + "combined_input.merge.bed")
    
#    return filedir + "combined_input.merge.bed"

#def run(BEDS,filedir):
#    os.system("sort -k1,1 -k2,2n " + BEDS[0] + " > " + filedir + "bed1_sort.bed")
#    os.system("sort -k1,1 -k2,2n " + BEDS[1] + " > " + filedir + "bed2_sort.bed")
#    os.system("bedtools intersect -v -a " + filedir + "bed1_sort.bed -b " + filedir + "bed2_sort.bed > " + filedir + "combined_inputa.bed")
#    os.system("bedtools intersect -v -b " + filedir + "bed1_sort.bed -a " + filedir + "bed2_sort.bed > " + filedir + "combined_inputb.bed")
#    os.system("cat " + filedir + "combined_inputa.bed " + filedir + "combined_inputb.bed > " + filedir + "combined_input.bed")
#    os.system("sort -k1,1 -k2,2n " + filedir + "combined_input.bed > " + filedir + "combined_input_sorted1.bed")
#    os.system("bedtools intersect -v -a " + filedir + "combined_input_sorted1.bed -b /scratch/Users/rusi2317/projects/rotation/bin/TFEA/files/hg19_TSS_sorted.bed > " + filedir + "combined_input_sorted.bed")
#    os.system("bedtools merge -i " + filedir + "combined_input_sorted.bed > " + filedir + "combined_input.merge.bed")
    
#    return filedir + "combined_input.merge.bed"
