__author__ = 'Joanthan Rubin'

import os

def run(BEDS,filedir):
    os.system("cat " + " ".join(BEDS) + " > " + filedir + "combined_input.bed")
    os.system("sort -k1,1 -k2,2n " + filedir + "combined_input.bed > " + filedir + "combined_input.sorted.bed")
    os.system("bedtools merge -i " + filedir + "combined_input.sorted.bed > " + filedir + "combined_input.merge.bed")
    # os.system("rm " + filedir + "combined_input.bed")
    return filedir + "combined_input.merge.bed"