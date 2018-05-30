__author__ = 'Jonathan Rubin'

import os
from config import FILEDIR, BEDS, GENOMEFASTA, LARGEWINDOW, RANKED_CENTER_FILE

def run():
    os.system("cat " + " ".join(BEDS) + " > " + FILEDIR + "combined_input.bed")
    os.system("sort -k1,1 -k2,2n " + FILEDIR + "combined_input.bed > " + FILEDIR + "combined_input.sorted.bed")
    os.system("bedtools merge -i " + FILEDIR + "combined_input.sorted.bed > " + FILEDIR + "combined_input.merge.bed")
    #os.system("bedtools intersect -v -a " + filedir + "combined_input.merge1.bed -b /scratch/Users/rusi2317/projects/rotation/bin/TFEA/files/hg19_TSS_sorted.bed > " + filedir + "combined_input.merge.bed") ##to remove TSS regions
    ## os.system("rm " + filedir + "combined_input.bed")
    return FILEDIR + "combined_input.merge.bed"

#5/23/18: This function returns a fasta file from a bed input
def getfasta(bedfile):
    os.system("bedtools getfasta -name -fi "+GENOMEFASTA+" -bed "+bedfile+" > combined_input.merge.fa")

    return FILEDIR + "combined_input.merge.fa"

#5/23/18: This function outputs a zero order markov background model file (is text file right?) for use with fimo later on
def get_bgfile(fastafile):
    os.system("fasta-get-markov "+fastafile+" "+FILEDIR+"markov_background.txt")

    return FILEDIR + "markov_background.txt"

def get_regions():
    outfile = open(FILEDIR+'ranked_file.fullregions.bed','w')
    with open(RANKED_CENTER_FILE) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            chrom,start,stop = line[:3]
            start = str(int(start)-LARGEWINDOW)
            stop = str(int(stop)+LARGEWINDOW)
            pval,fc,rank = line[3:]
            name = ','.join([pval,fc,rank])
            outfile.write('\t'.join([chrom,start,stop,name]) + '\n')

    return FILEDIR + 'ranked_file.fullregions.bed'




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
