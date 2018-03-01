__author__ = 'Jonathan Rubin'

import math
import os

def deseqfile(DESEQFILE,filedir):
    #Parse a deseq file and obtain the exact middle of each region (for motif distance calc later) and pvalue (to rank)
    ranked = list()
    up = list()
    down = list()
    with open(DESEQFILE) as F:
        F.readline()
        for line in F:
            line = line.strip('\n').split('\t')
            try:
                pval = format(float(line[-2]),'.12f')
            except ValueError:
                pval = format(1.0,'.12f')
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

    #Get center base for each region
    outfile = open(filedir+"ranked_file.center.bed",'w')
    with open(filedir + "ranked_file.bed") as F:
        for line in F:
            line = line.strip('\n').split('\t')
            chrom,start,stop = line[:3]
            center = (int(start)+int(stop))/2
            outfile.write(chrom + '\t' + str(center) + '\t' + str(center+1) + '\t' + '\t'.join(line[3:]) + '\n')
    outfile.close()


    os.system("sort -k1,1 -k2,2n " + filedir+"ranked_file.center.bed" + " > " + filedir + "ranked_file.center.sorted.bed")
    os.system("rm " + filedir + "combined_input.bed")
    os.system("rm " + filedir + "combined_input.merge.bed")
    os.system("rm " + filedir + "combined_input.sorted.bed")
    os.system("rm " + filedir + "count_file.bed")
    os.system("rm " + filedir + "count_file.header.bed")
    os.system("rm " + filedir + "ranked_file.bed")
    os.system("rm " + filedir + "ranked_file.center.bed")



