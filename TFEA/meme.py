__author__ = 'Jonathan Rubin'

import os
from config import FILEDIR, MOTIFDATABASE, FIGUREDIR

#5/23/18:This function runs fimo on a given fastafile for a single motif in a provided motif database. The output is cut and sorted to convert into a sorted bed file
def fimo(bgfile,motif,fastafile):
    os.system("fimo --text --bgfile "+bgfile+" --motif "+motif+" "+MOTIFDATABASE+" "+fastafile+" > "+FILEDIR+motif+".txt")
    os.system("cut -d$'\\t' -f -2 "+FILEDIR+motif+".txt | sort -k1,1 -k2,2n > " + FILEDIR+motif+".sorted.bed")

#5/23/18: This function creates meme logos for use in the output html
def meme2images(motif):
    os.system("meme2images -png -rc --motif "+motif+" "+MOTIFDATABASE+" "+FIGUREDIR)

#5/29/18: This function runs meme's fasta-get-markov function that generates a background markov file (for use with fimo) from a fasta file. This function is also within combine_bed.py
def fasta_markov(fastafile,order='0'):
    os.system("fasta-get-markov -m " + order + " " + fastafile + " > " + FILEDIR + "Markov_Background.txt")