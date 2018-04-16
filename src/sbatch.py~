__author__ = 'Jonathan Rubin'

import os
import main
import config

def run(script):
    output,filedir,figuredir,e_and_o = main.make_out_directories(True)
    os.system("sbatch --error=" + e_and_o + "%x.err --output=" + e_and_o + "%x.out --export=output="+output+",filedir=" +filedir+",figuredir="
                +figuredir+",e_and_o=" + e_and_o + " " + script)


#Return parent directory
def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])
    
    return newdir

if __name__ == "__main__":
    homedir = os.path.dirname(os.path.realpath(__file__))
    scriptdir = parent_dir(homedir) + '/scripts/'
    script = scriptdir + 'run_main.sbatch'
    run(script)