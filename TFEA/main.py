__author__ = 'Jonathan Rubin'

import os
import time
import sys
import argparse
import configparser
import config_parser
from multiprocessing import Pool
import multiprocessing as mp

def run():
    #Home directory, gets the full path (no '/' at the end) to the folder containing this script
    homedir = os.path.dirname(os.path.realpath(__file__))

    parser = argparse.ArgumentParser(description='Transcription Factor Enrichment Analysis (TFEA) takes as input a configuration file (.ini) and outputs a folder containing TFEA results.',usage='TFEA --config CONFIG.ini [--sbatch email@address.com]')
    parser.add_argument('--config',help='REQUIRED. A configuration file containing .ini suffix (ex. config.ini). See example in the examples folder.')
    parser.add_argument('--sbatch',default=False,help='OPTIONAL. Submits an sbatch job. If specified, input an e-mail address.')
    if len(sys.argv)==1:
        # display help message when no args are passed.
        parser.print_help()
        sys.exit(1)
    sbatch = parser.parse_args().sbatch
    configfile = parser.parse_args().config
    config_parser.run(homedir+'/',str(configfile))

    if sbatch == False:
        output,filedir,figuredir,e_and_o = make_out_directories(True)
    elif str(sbatch) == 'SUBMITTED':
        output,filedir,figuredir,e_and_o = make_out_directories(False)
    else:
        output,filedir,figuredir,e_and_o = make_out_directories(True)
        scriptdir = parent_dir(homedir) + '/scripts/'
        script = scriptdir + 'run_main.sbatch'
        email = str(sbatch)
        os.system("sbatch --error=" + e_and_o + "%x.err --output=" + e_and_o + "%x.out --mail-user="+email+" --export=src="+homedir+",config=" +configfile+ " " + script)
        sys.exit("TFEA has been submitted using an sbatch script, use qstat to check its progress.")


    import combine_bed
    import count_reads
    import rank_regions
    import DESeq
    import motif_distance
    import ES_calculator
    import create_html
    import config
    import meta_eRNA

    #Path to count file. Can be changed if using your own count file. Generated in count_reads module
    count_file = filedir + "count_file.header.bed"

    #Path to DESeq file. Can be changed if using your own DESeq file. Generated in DESeq module
    deseq_file = filedir + "DESeq.res.txt"

    #Path to ranked file. Can be changed if using your own ranked file. Generated in rank_regions module
    ranked_file = filedir + "ranked_file.bed"

    #Path to ranked centered file. Just a bed file with single basepair coordinates for the exact middle of each bed region
    ranked_center_file = filedir + "ranked_file.center.bed"

    #Path to the centered ranked file with measures of distance to the motif
    ranked_center_distance_file = filedir + "ranked_file.center.sorted.distance.bed"

    #Path to a directory full of motif logos for all TFs in the HOCOMOCO database (v10)
    logos = parent_dir(homedir) + '/human_logo/'

    #Path to mouse directory with motif logos in HOCOMOCO v10
    ##logos = parent_dir(homedir) + '/mouse_logo/'


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
        DESeq.run(config.LABEL1,config.LABEL2,config.BAM1,config.BAM2,filedir,count_file)
        rank_regions.deseqfile(deseq_file,filedir)
        print "done"
    DESEQtime = time.time()-DESEQtime

    #Scans ranked BED regions for motifs of interest and records them in distance file
    if config.CALCULATE:
        # p = Pool(mp.cpu_count())
        # millions_mapped = p.map(meta_eRNA.get_millions_mapped_pool,config.BAM1+config.BAM2)
        # millions_mapped = meta_eRNA.get_millions_mapped(config.BAM1+config.BAM2)

        p = Pool(64)
        args = [(x,filedir) for x in config.BAM1+config.BAM2]
        millions_mapped = p.map(meta_eRNA.samtools_flagstat,args)

        # millions_mapped = list()
        # for bam in config.BAM1+config.BAM2:
        #     millions_mapped.append(meta_eRNA.samtools_flagstat(bam,filedir))

        print millions_mapped
        print "Finding motif hits in regions..."
        if config.SINGLEMOTIF == False:
            TFresults = list()
            NESlist = list()
            CALCULATEtime = 0.0
            if config.POOL:
                a = time.time()
                args = [(x,ranked_center_distance_file,ranked_center_file,figuredir,millions_mapped,logos) for x in os.listdir(config.MOTIF_HITS)]
                cpus = mp.cpu_count()
                if cpus > 64:
                    cpus = 64
                p = Pool(cpus)
                TFresults = p.map(ES_calculator.run,args)
                CALCULATEtime += time.time() - a
                create_html.createTFtext(TFresults,output)
            else:
                for MOTIF_FILE in os.listdir(config.MOTIF_HITS):
                    a = time.time()
                    # motif_distance.run(ranked_center_file,config.MOTIF_HITS+MOTIF_FILE)

                    #This module is where the bulk of the analysis is done. The functions below calculate ES,NES,p-value,FDR for each TF motif in
                    #the HOCOMOCO database.
                    results = ES_calculator.run((MOTIF_FILE,ranked_center_distance_file,ranked_center_file,figuredir,millions_mapped,logos,filedir))
                    if results != "no hits":
                        TFresults.append(results)
                        NESlist.append(results[2])
                        CALCULATEtime += time.time()-a
                        print MOTIF_FILE + " calculation done in: " + str(CALCULATEtime) + "s"
                    else:
                        print "No motifs within specified window for: ", MOTIF_FILE
            TFresults = [x for x in TFresults if x != "no hits"]
            TFresults = sorted(TFresults, key=lambda x: x[3])
            TFresults = ES_calculator.FDR(TFresults,NESlist,figuredir)
            create_html.run(TFresults,output,COMBINEtime,COUNTtime,DESEQtime,CALCULATEtime)

        #Note if you set the SINGLEMOTIF variable to a specific TF, this program will be unable to accurately determine an FDR for the given motif.
        else:
            # motif_distance.run(ranked_center_file,config.MOTIF_HITS+config.SINGLEMOTIF)
            results = ES_calculator.run((config.SINGLEMOTIF,ranked_center_distance_file,ranked_center_file,figuredir,millions_mapped,logos,filedir))
            create_html.single_motif(results,output)
    print "done"


def make_out_directories(dirs):
    import config
    #Output directory
    output = config.OUTPUT
    if dirs:
        if not os.path.isdir(output + 'TFEA_output-0/'):
            output = output + 'TFEA_output-0/'
            os.makedirs(output)
        else:
            outputfolders = list()
            for folder in os.listdir(output):
                if 'TFEA_output' in folder:
                    outputfolders.append(int(folder.split('-')[1]))
            output = output + 'TFEA_output-' + str(max(outputfolders)+1) + '/'
            os.makedirs(output)
    else:
        outputfolders = list()
        for folder in os.listdir(output):
            if 'TFEA_output' in folder:
                outputfolders.append(int(folder.split('-')[1]))
        output = output + 'TFEA_output-' + str(max(outputfolders)) + '/'


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

#Return parent directory
def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])
    
    return newdir
