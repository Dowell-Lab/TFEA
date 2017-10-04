__author__ = 'Jonathan Rubin'

import os

def run(BED,BAM1,BAM2,filedir):
    os.system("bedtools multicov -bams " + " ".join(BAM1) + " " + " ".join(BAM2) + " -bed " + BED + " > " + filedir + "count_file.bed")
    print "bedtools multicov -bams " + " ".join(BAM1) + " " + " ".join(BAM2) + " -bed " + BED + " > " + filedir + "count_file.bed"
