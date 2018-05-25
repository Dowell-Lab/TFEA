__author__ = 'Jonathan Rubin'

import os
import main
import config

def run(homedir,configfile,script,email):
    output,filedir,figuredir,e_and_o = main.make_out_directories(True)
    print "sbatch --error=" + e_and_o + "%x.err --output=" + e_and_o + "%x.out --mail-user="+email+" --export=src="+homedir+",config=" +configfile+ " " + script
    os.system("sbatch --error=" + e_and_o + "%x.err --output=" + e_and_o + "%x.out --mail-user="+email+" --export=src="+homedir+",config=" +configfile+ " " + script)

#Return parent directory
def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])
    
    return newdir
