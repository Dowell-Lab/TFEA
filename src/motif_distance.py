__author__ = 'Jonathan Rubin'

import os

def runfimo(ranked_file,filedir,MEMEDB,DATABASE,SINGLEMOTIF):
    if os.path.isdir(filedir + "fimo_out/"):
        os.system("rm -r " + filedir + "fimo_out/")
    if SINGLEMOTIF == False:
        os.system("fimo -o " + filedir + "fimo_out/ " + MEMEDB + DATABASE + " " + ranked_file)
    else:
        os.system("fimo -o " + filedir + "fimo_out/ " + "--motif " + SINGLEMOTIF + " " + MEMEDB + DATABASE + " " + ranked_file)

def run(ranked_file,filedir,MOTIF_HITS,SINGLEMOTIF):
    outfile = open(filedir+"ranked_file.center.bed")
    with open(ranked_file) as F:
        for line in F:
            chrom,start,stop = line.strip('\n').split('\t')
            center = (int(start)+int(stop))/2
            outfile.write(chrom + '\t' + str(center) + '\t' + str(center+1) + '\n')

    if SINGLEMOTIF == False:
        print "This module isn't finished yet"
    else:
        for MOTIF_FILE in os.listdir(MOTIF_HITS):
            if MOTIF_FILE == SINGLEMOTIF:
                os.system("bedtools closest -d -a " + filedir + "ranked_file.center.bed -b " + MOTIF_HITS + MOTIF_FILE + " > " + filedir + "ranked_file.center.distance.bed")