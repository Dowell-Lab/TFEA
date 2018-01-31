__author__ = 'Jonathan Rubin'

import os
import time
import sys
import combine_bed
import count_reads
import rank_regions
import DESeq
import motif_distance
import ES_calculator
import create_html
import config

def make_out_directories(dirs):
    #Output directory
    output = config.OUTPUT
    if not os.path.isdir(output + 'TFEA_output-0/'):
        output = output + 'TFEA_output-0/'
        if dirs:
            os.makedirs(output)
    else:
        outputfolders = list()
        for folder in os.listdir(output):
            if 'TFEA_output' in folder:
                outputfolders.append(int(folder.split('-')[1]))
        output = output + 'TFEA_output-' + str(max(outputfolders)+1) + '/'
        if dirs:
            os.makedirs(output)


    #Temporary files will go in this directory
    filedir = output + 'temp_files/'
    if dirs:
        if not os.path.isdir(filedir):
            os.makedirs(filedir)

    #Error and out files will go in this directory
    e_and_o = output + 'e_and_o/'
    if dirs:
        if not os.path.isdir(e_and_o):
            os.makedirs(e_and_o)


    #Directory where plots used in html file will be stored.
    figuredir = output + 'plots/'
    if dirs:
        if not os.path.isdir(figuredir):
            os.makedirs(figuredir)

    return output,filedir,figuredir,e_and_o

def run():
    #Home directory
    homedir = os.path.dirname(os.path.realpath(__file__))

    directories = sys.argv[1:]
    print len(directories)
    if len(directories) == 0:
        output,filedir,figuredir,e_and_o = make_out_directories(True)
    else:
        output,filedir,figuredir,e_and_o = directories

    #Directory where all temp files will be stored
    # filedir = parent_dir(homedir) + '/files/'

    #Path to count file. Can be changed if using your own count file. Generated in count_reads module
    count_file = filedir + "count_file.header.bed"

    #Path to DESeq file. Can be changed if using your own DESeq file. Generated in DESeq module
    deseq_file = output + "DESeq.res.txt"

    #Path to ranked file. Can be changed if using your own ranked file. Generated in rank_regions module
    ranked_file = filedir + "ranked_file.bed"

    #Path to ranked centered file. Just a bed file with single basepair coordinates for the exact middle of each bed region
    ranked_center_file = filedir + "ranked_file.center.bed"

    #Path to the centered ranked file with measures of distance to the motif
    ranked_center_distance_file = filedir + "ranked_file.center.sorted.distance.bed"

    #Path to a directory full of motif logos for all TFs in the HOCOMOCO database (v10)
    logos = parent_dir(homedir) + '/logo/'

    #This module takes the input list of BED files, concatenates them, and then merges them via bedtools.
    COMBINEtime = time.time()
    if config.COMBINE:
        BED = combine_bed.run(config.BEDS,filedir)
    else:
        BED = config.BEDS[0]
    COMBINEtime = time.time()-COMBINEtime 

    #This module counts reads from all Bam files in BAM1 and BAM2 and creates count_file with this info.
    COUNTtime = time.time()
    if config.COUNT:
        print "Counting reads in regions..."
        count_reads.run(BED,config.BAM1,config.BAM2,config.LABEL1,config.LABEL2,filedir)
        print "done"
    COUNTtime = time.time()-COUNTtime

    #This module runs DESeq on specified 
    DESEQtime = time.time()
    if config.DESEQ:
        print "Running DESeq..."
        DESeq.run(config.LABEL1,config.LABEL2,config.BAM1,config.BAM2,output,count_file)
        rank_regions.deseqfile(deseq_file,filedir)
        print "done"
    DESEQtime = time.time()-DESEQtime

    #Scans ranked BED regions for motifs of interest and records them in distance file
    if config.CALCULATE:
        print "Finding motif hits in regions..."
        if config.SINGLEMOTIF == False:
            TFresults = list()
            NESlist = list()
            CALCULATEtime = 0.0
            for MOTIF_FILE in os.listdir(config.MOTIF_HITS):
                a = time.time()
                motif_distance.run(ranked_center_file,config.MOTIF_HITS+MOTIF_FILE)

                #This module is where the bulk of the analysis is done. The functions below calculate ES,NES,p-value,FDR for each TF motif in
                #the HOCOMOCO database.
                results = ES_calculator.run(MOTIF_FILE,ranked_center_distance_file,figuredir,logos)
                TFresults.append(results)
                CALCULATEtime += time.time()-a
                print MOTIF_FILE + " calculation done in: " + str(CALCULATEtime) + "s"

            TFresults = ES_calculator.FDR(TFresults,figuredir)
            create_html.run(TFresults,output,COMBINEtime,COUNTtime,DESEQtime,CALCULATEtime)

        #Note if you set the SINGLEMOTIF variable to a specific TF, this program will be unable to accurately determine an FDR for the given motif.
        else:
            motif_distance.run(ranked_file,filedir,config.MOTIF_HITS+SINGLEMOTIF)
            results = ES_calculator.run(config.SINGLEMOTIF,ranked_center_distance_file,figuredir,logos)
            create_html.single_motif(results,output)
        print "done"


#Return parent directory
def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])
    
    return newdir