__author__ = 'Jonathan Rubin'

import os

def run(ranked_file,filedir,MEME,DATABASE,SINGLEMOTIF):
    if SINGLEMOTIF == False:
        os.system("fimo -o " + filedir + "fimo_out/ " + MEME + '/db/fasta_databases/' + DATABASE + " " + ranked_file)
    else:
        os.system("fimo -o " + filedir + "fimo_out/" + " --motif " + SINGLEMOTIF + " " + MEME + '/db/fasta_databases/' + DATABASE + " " + ranked_file)
