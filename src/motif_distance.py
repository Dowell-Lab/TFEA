__author__ = 'Jonathan Rubin'

import os

def run(ranked_file,filedir,MEMEDB,DATABASE,SINGLEMOTIF):
    if SINGLEMOTIF == False:
        os.system("fimo -o " + filedir + "fimo_out/ " + MEMEDB + DATABASE + " " + ranked_file)
    else:
        os.system("fimo -o " + filedir + "fimo_out/" + " --motif " + SINGLEMOTIF + " " + MEMEDB + DATABASE + " " + ranked_file)
