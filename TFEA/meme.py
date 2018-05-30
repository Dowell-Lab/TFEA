__author__ = 'Jonathan Rubin'

import os

#5/23/18:This function runs fimo on a given fastafile for a single motif in a provided motif database. The output is cut and sorted to convert into a sorted bed file
def fimo(bgfile,motif,motifdatabase,fastafile,filedir):
    os.system("fimo --text --bgfile "+bgfile+" --motif "+motif+" "+motifdatabase+" "+fastafile+" > "+filedir+motif+".txt")
    os.system("cut -d$'\\t' -f -2 | sort -k1,1 -k2,2n > " + filedir+motif+".sorted.bed")

#5/23/18: This function creates meme logos for use in the output html
def meme2images(outputdir,motif,motifdatabase):
    os.system("meme2images -png -rc --motif "+motif+" "+motifdatabase+" "+outputdir)

#5/29/18: This function runs meme's fasta-get-markov function that generates a background markov file (for use with fimo) from a fasta file. This function is also within combine_bed.py
def fasta_markov(order='0',fastafile,filedir):
    os.system("fasta-get-markov -m " + order + " " + fastafile + " > " + fildir + "Markov_Background.txt")