__author__ = 'Jonathan Rubin'

import os

#5/23/18:This function runs fimo on a given fastafile for a single motif in a provided motif database. The output is cut and sorted to convert into a sorted bed file
def fimo(bgfile,motif,motifdatabase,fastafile,filedir):
    os.system("fimo --text --bgfile "+bgfile+" --motif "+motif+" "+motifdatabase+" "+fastafile+" > "+filedir+motif+".txt")
    os.system("cut -d$'\\t' -f -2 | sort -k1,1 -k2,2n > " + filedir+motif+".sorted.bed")

#5/23/18: This function creates meme logos for use in the output html
def meme2images(outputdir,motif,motifdatabase):
    os.system("meme2images -png -rc --motif "+motif+" "+motifdatabase+" "+outputdir)