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
    outfile.write('GC_ARRAY=list()\n')

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

#This function verifies that the config file contains all necessary variables
def verify():
    import config
    try:
        if type(config.COMBINE) != bool:
            raise TypeError('COMBINE variable must be a boolean. This switch determines whether TFEA merges bed files within BED input.')
    except NameError:
        raise NameError('COMBINE variable not found in config.ini file. This switch determines whether TFEA merges bed files within BED input.')

    try:
        if type(config.COUNT) != bool:
            raise TypeError('COUNT variable must be a boolean. This switch determines whether TFEA performs read counting over BED regions.')
    except NameError:
        raise NameError('COUNT variable not found in config.ini file. This switch determines whether TFEA performs read counting over BED regions.')

    try:
        if type(config.DESEQ) != bool:
            raise TypeError('DESEQ variable must be a boolean. This switch determines whether TFEA performs DE-Seq analysis on counted BED regions.')
    except NameError:
        raise NameError('DESEQ variable not found in config.ini file. This switch determines whether TFEA performs DE-Seq analysis on counted BED regions.')

    try:
        if type(config.CALCULATE) != bool:
            raise TypeError('CALCULATE variable must be a boolean. This switch determines whether TFEA performs its standard enrichment score calculation and plotting.')
    except NameError:
        raise NameError('CALCULATE variable not found in config.ini file. This switch determines whether TFEA performs its standard enrichment score calculation and plotting.')

    try:
        if type(config.POOL) != bool:
            raise TypeError('POOL variable must be a boolean. This switch determines whether TFEA runs the analysis in parallel using the multiprocessing library in python.')
    except NameError:
        raise NameError('POOL variable not found in config.ini file. This switch determines whether TFEA runs the analysis in parallel using the multiprocessing library in python.')

    try:
        if type(config.SINGLEMOTIF) != bool or type(config.SINGLEMOTIF) != str:
            raise TypeError('SINGLEMOTIF variable must be a boolean or string. This switch determines whether TFEA performs its analysis on a single motif or all. If not False, set to a string matching a motif name.')
    except NameError:
        raise NameError('SINGLEMOTIF variable not found in config.ini file. This switch determines whether TFEA performs its analysis on a single motif or all. If not False, set to a string matching a motif name.')

    try:
        if type(config.FIMO) != bool:
            raise TypeError('FIMO variable must be a boolean. This switch determines whether TFEA uses FIMO to get motif hits or whether a database of motif hit calls (bed format) is used.')
    except NameError:
        raise NameError('FIMO variable not found in config.ini file. This switch determines whether TFEA uses FIMO to get motif hits or whether a database of motif hit calls (bed format) is used.')

    try:
        if type(config.TEMP) != bool:
            raise TypeError('TEMP variable must be a boolean. This switch determines whether TFEA saves large temporary files. If True, temporary files will be stored in the temp_files directory. Warning: There will be many large files.')
    except NameError:
        raise NameError('TEMP variable not found in config.ini file. This switch determines whether TFEA saves large temporary files. If True, temporary files will be stored in the temp_files directory. Warning: There will be many large files.')

    try:
        if type(config.PLOT) != bool:
            raise TypeError('PLOT variable must be a boolean. This switch determines whether TFEA outputs plots for all motifs provided or just significant ones. Warning: Setting this to True will slow down TFEA and create large output folders.')
    except NameError:
        raise NameError('PLOT variable not found in config.ini file. This switch determines whether TFEA outputs plots for all motifs provided or just significant ones. Warning: Setting this to True will slow down TFEA and create large output folders.')

    try:
        if type(config.OUTPUT) != str:
            raise TypeError('OUTPUT variable must be a string. Determines where TFEA stores output.')
    except NameError:
        raise NameError('OUTPUT variable not found in config.ini file. Determines where TFEA stores output.')

    try:
        if type(config.BEDS) != list:
            raise TypeError('BED variable must be a list. Input a list of regions (bed-format) to perform analysis over. If merging not desired, simply input a single BED file within a python list.')
    except NameError:
        raise NameError('BED variable not found in config.ini file. Input a list of regions (bed-format) to perform analysis over. If merging not desired, simply input a single BED file within a python list.')

    try:
        if type(config.BAM1) != list:
            raise TypeError('BAM1 variable must be a list. Input a list of BAM files to obtain read depth and coverage over regions of interest. One or multiple bams can be specified but they must be within a python list.')
    except NameError:
        raise NameError('BAM1 variable not found in config.ini file. Input a list of BAM files to obtain read depth and coverage over regions of interest. One or multiple bams can be specified but they must be within a python list.')

    try:
        if type(config.LABEL1) != str:
            raise TypeError('LABEL1 variable must be a string. Define a treatment or condition to label BAM1 files. Used in plotting and output folder naming.')
    except NameError:
        raise NameError('LABEL1 variable not found in config.ini file. Define a treatment or condition to label BAM1 files. Used in plotting and output folder naming.')

    try:
        if type(config.BAM2) != list:
            raise TypeError('BAM2 variable must be a string.')
    except NameError:
        raise NameError('BAM2 variable not found in config.ini file.')

    try:
        if type(config.LABEL2) != str:
            raise TypeError('LABEL2 variable must be a string.')
    except NameError:
        raise NameError('LABEL2 variable not found in config.ini file.')

    try:
        if type(config.PADJCUTOFF) != float:
            raise TypeError('PADJCUTOFF variable must be a string.')
    except NameError:
        raise NameError('PADJCUTOFF variable not found in config.ini file.')

    try:
        if type(config.LARGEWINDOW) != float:
            raise TypeError('LARGEWINDOW variable must be a string.')
    except NameError:
        raise NameError('LARGEWINDOW variable not found in config.ini file.')

    try:
        if type(config.SMALLWINDOW) != float:
            raise TypeError('SMALLWINDOW variable must be a string.')
    except NameError:
        raise NameError('SMALLWINDOW variable not found in config.ini file.')

    try:
        if type(config.MOTIF_HITS) != float:
            raise TypeError('MOTIF_HITS variable must be a string.')
    except NameError:
        raise NameError('MOTIF_HITS variable not found in config.ini file.')

    try:
        if type(config.GENOMEFASTA) != float:
            raise TypeError('GENOMEFASTA variable must be a string.')
    except NameError:
        raise NameError('GENOMEFASTA variable not found in config.ini file.')

    try:
        if type(config.MOTIFDATABASE) != float:
            raise TypeError('MOTIFDATABASE variable must be a string.')
    except NameError:
        raise NameError('MOTIFDATABASE variable not found in config.ini file.')

    print "Config file verified, all inputs present and correct type."




