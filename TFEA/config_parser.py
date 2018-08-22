__author__ = 'Jonathan Rubin'

import configparser
import os

def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])

    return newdir


def run(srcdirectory,config,output,filedir,figuredir):
    outfile = open(srcdirectory+'config.py','w')
    for key in config:
        for item in config[key]:
            outfile.write(item.upper()+'='+config[key][item]+'\n')

    homedir = os.path.dirname(os.path.realpath(__file__))

    outfile.write('OUTPUTDIR="'+output+'"\n')
    outfile.write('FILEDIR="'+filedir+'"\n')
    outfile.write('FIGUREDIR="'+figuredir+'"\n')
    outfile.write('GC_ARRAY=list()')

    #Path to count file. Can be changed if using your own count file. Generated in count_reads module
    count_file = filedir + "count_file.header.bed"
    outfile.write('COUNT_FILE="'+count_file+'"\n')

    #Path to DESeq file. Can be changed if using your own DESeq file. Generated in DESeq module
    deseq_file = filedir + "DESeq.res.txt"
    outfile.write('DESEQ_FILE="'+deseq_file+'"\n')

    #Path to ranked file. Can be changed if using your own ranked file. Generated in rank_regions module
    ranked_file = filedir + "ranked_file.bed"
    outfile.write('RANKED_FILE="'+ranked_file+'"\n')

    #Path to ranked centered file. Just a bed file with single basepair coordinates for the exact middle of each bed region
    ranked_center_file = filedir + "ranked_file.center.bed"
    outfile.write('RANKED_CENTER_FILE="'+ranked_center_file+'"\n')

    #Path to the centered ranked file with measures of distance to the motif
    ranked_center_distance_file = filedir + "ranked_file.center.sorted.distance.bed"
    outfile.write('RANKED_CENTER_DISTANCE_FILE="'+ranked_center_distance_file+'"\n')

    #Path to a directory full of motif logos for all TFs in the HOCOMOCO database (v10)
    logos = parent_dir(homedir) + '/human_logo/'
    #logos = srcdirectory + 'human_logo/'
    outfile.write('LOGOS="'+logos+'"\n')

    #Path to mouse directory with motif logos in HOCOMOCO v10
    ##logos = parent_dir(homedir) + '/mouse_logo/'

    outfile.close()

