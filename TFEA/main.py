__author__ = 'Jonathan Rubin'

import os,time,sys,argparse,configparser, config_parser
from multiprocessing import Pool
import multiprocessing as mp

def run():
    #Home directory, gets the full path (no '/' at the end) to the folder containing this script
    homedir = os.path.dirname(os.path.realpath(__file__))

    #argparse to add arguments to this python package
    parser = argparse.ArgumentParser(description='Transcription Factor Enrichment Analysis (TFEA) takes as input a configuration file (.ini) and outputs a folder containing TFEA results.',usage='TFEA --config CONFIG.ini [--sbatch email@address.com]')
    parser.add_argument('--config','-c',metavar='',help='REQUIRED. A configuration file containing .ini suffix (ex. config.ini). See example in the examples folder.')
    parser.add_argument('--sbatch','-s',default=False,metavar='',help='OPTIONAL. Submits an sbatch job. If specified, input an e-mail address.')

    #Display help message when no args are passed.
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    #If user provided arguments, then parse them
    sbatch = parser.parse_args().sbatch
    configfile = parser.parse_args().config
    temp = parser.parse_args().temp
    plot = parser.parse_args().plot
    config = configparser.ConfigParser(interpolation = configparser.ExtendedInterpolation())
    config.read(configfile)

    #Run the config_parser script which will create variables for all folders and paths to use throughout TFEA
    config_parser.run(homedir+'/',config,output,filedir,figuredir)

    #Verify config file to make sure user has inputted all necessary variables
    # config_parser.verify()

    #If user specifies the --sbatch flag, then we first create the output directories then run the sbatch script with the 'SUBMITTED' command submitted to the 
    #--sbatch flag so we know not to remake output directories. If --sbatch flag not specified, simply make output directories and continue.
    if sbatch == False:
        output,filedir,figuredir,e_and_o = make_out_directories(True,config)
    elif str(sbatch) == 'SUBMITTED':
        output,filedir,figuredir,e_and_o = make_out_directories(False,config)
    else:
        output,filedir,figuredir,e_and_o = make_out_directories(True,config)
        scriptdir = parent_dir(homedir) + '/scripts/'
        script = scriptdir + 'run_main.sbatch'
        email = str(sbatch)
        # os.system("sbatch --error=" + e_and_o + "%x.err --output=" + e_and_o + "%x.out --mail-user="+email+" --export=src="+homedir+",config=" +configfile+ " " + script)
        os.system("sbatch --error=" + e_and_o + "%x.err --output=" + e_and_o + "%x.out --mail-user="+email+" --export=cmd='"+homedir+" --config " +configfile+ " --sbatch SUBMITTED' " + script)
        sys.exit("TFEA has been submitted using an sbatch script, use qstat to check its progress.")


    #Import scripts from this package
    import config, combine_bed, count_reads, rank_regions, DESeq, motif_distance, ES_calculator, create_html, meta_eRNA


    #This module takes the input list of BED files, concatenates them, and then merges them via bedtools.
    COMBINEtime = time.time()
    if config.COMBINE:
        BED = combine_bed.run()
    else:
        BED = config.BEDS[0]
    COMBINEtime = time.time()-COMBINEtime 

    #This module counts reads from all Bam files in BAM1 and BAM2 and creates count_file with this info.
    COUNTtime = time.time()
    if config.COUNT:
        print "Counting reads in regions..."
        count_reads.run(BED)
        print "done"
    COUNTtime = time.time()-COUNTtime

    #This module runs DESeq on the count_file produced above - this is used to rank inputted regions
    DESEQtime = time.time()
    if config.DESEQ:
        print "Running DESeq..."
        DESeq.run()
        rank_regions.deseqfile()
        print "done"
    DESEQtime = time.time()-DESEQtime


    #This is the bulk of the analysis of this package, it performs:
    #   1. GC distribution across regions for plotting
    #   2. Millions mapped calculation for meta eRNA plots
    #   3. Motif distance calculation to the center of inputted regions
    #   4. Enrichment score calculation via AUC method
    #   5. Random shuffle simulation and recalculation of enrichment score
    #   6. Plotting and generation of html report

    CALCULATEtime = time.time()
    if config.CALCULATE:
        print "Calculating GC content of regions..."
        #This line gets an array of GC values for all inputted regions
        motif_distance.get_gc(config.RANKED_FILE)

        print "done\nCalculating millions mapped reads for bam files..."
        #Here we determine how many cpus to use for parallelization
        cpus = mp.cpu_count()
        if cpus > 64:
            cpus = 64

        #Here we calculate millions mapped reads for use with the metaeRNA module
        p = Pool(cpus)
        args = [(x) for x in config.BAM1+config.BAM2]
        millions_mapped = p.map(meta_eRNA.samtools_flagstat,args)

        print "done\nFinding motif hits in regions..."
        if config.SINGLEMOTIF == False:
            TFresults = list()
            if config.POOL:
                a = time.time()
                args = [(x,millions_mapped) for x in os.listdir(config.MOTIF_HITS)]
                p = Pool(cpus)
                TFresults = p.map(ES_calculator.run, args)
            else:
                for MOTIF_FILE in os.listdir(config.MOTIF_HITS):
                    results = ES_calculator.run((MOTIF_FILE,millions_mapped))
                    if results != "no hits":
                        TFresults.append(results)
                    else:
                        print "No motifs within specified window for: ", MOTIF_FILE
            CALCULATEtime = time.time()-CALCULATEtime
            create_html.createTFtext(TFresults, output)
            TFresults = ES_calculator.PADJ(TFresults)
            create_html.run(TFresults,COMBINEtime,COUNTtime,DESEQtime,CALCULATEtime)

        #Note if you set the SINGLEMOTIF variable to a specific TF, this program will be unable to determine an PADJ for the given motif.
        else:
            results = ES_calculator.run((config.SINGLEMOTIF,millions_mapped))
            create_html.single_motif(results)

        print "done"

    #Here we simply remove large bed files that are produced within this package. This option can be turned off by specifying the --temp flag
    if not config.TEMP:
        print "Removing temporary bed files..."
        os.system("rm " + filedir + '*.sorted.distance.bed')
        os.system("rm " + filedir + '*.fa')

    print "done"

#This function detects whether there are existing output folders that match the current output folder and if so sequentially adds an integer to the end
#of the output folder name so it doesn't overwrite
def make_out_directories(dirs,config):
    #Output directory
    output = config['DATA']['OUTPUT'].strip("'")
    label1 = config['DATA']['LABEL1'].strip("'")
    label2 = config['DATA']['LABEL2'].strip("'")
    outfoldername = 'TFEA_'+label1+'-'+label2+'_'
    if dirs:
        if not os.path.isdir(output + outfoldername + '0/'):
            output = output + outfoldername + '0/'
            os.makedirs(output)
        else:
            outputfolders = list()
            for folder in os.listdir(output):
                if outfoldername in folder:
                    outputfolders.append(int(folder.split('_')[-1]))
            output = output + outfoldername + str(max(outputfolders)+1) + '/'
            os.makedirs(output)
    else:
        outputfolders = list()
        for folder in os.listdir(output):
            if outfoldername in folder:
                outputfolders.append(int(folder.split('_')[-1]))
        output = output + outfoldername + str(max(outputfolders)) + '/' 


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
