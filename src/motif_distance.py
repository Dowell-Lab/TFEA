__author__ = 'Jonathan Rubin'

import os

def run(ranked_file,filedir):
    os.system("fimo -o " + filedir + "fimo_out/ " + ranked_file)
