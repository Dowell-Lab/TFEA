#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''This file contains a list of independent functions that do not call other 
    functions
'''
#==============================================================================
__author__ = 'Jonathan D. Rubin and Rutendo Sigauke'
__credits__ = ['Jonathan D. Rubin', 'Rutendo Sigauke', 'Jacob Stanley',
                'Robin Dowell']
__maintainer__ = 'Jonathan D. Rubin'
__email__ = 'Jonathan.Rubin@colorado.edu'
#==============================================================================
import matplotlib
matplotlib.use('Agg')
import os
import math
import matplotlib.pyplot as plt
import numpy as np
import config
#==============================================================================
#Functions
#==============================================================================
def make_out_directories(dirs=False,config=None):
    '''Creates output directories in a user-specified location where all TFEA 
        outputs will go.

    Parameters
    ----------
    dirs : boolean
        determines whether output folders will be created or not (default: 
        False)
        
    config : dict
        a configparser object that contains variables within the config file

    Returns
    -------
    output : string 
        full path to the parent output directory

    tempdir : string
        full path to the temporary directory where files are stored

    figuredir : string
        full path to the directory containing figures and plots
        
    e_and_o : string
        full path to the directory that stores stdout and stderr files
    '''
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
    tempdir = output + 'temp_files/'
    if dirs:
        if not os.path.isdir(tempdir):
            os.makedirs(tempdir)

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

    return output,tempdir,figuredir,e_and_o
#==============================================================================

#==============================================================================
def parse_config(srcdirectory='',config='',output='',tempdir='',figuredir=''):
    '''Creates the config.py file which is used in many aspects of TFEA. Within
        this config.py file, it writes all variables provided in the config
        parameter and also writes output, tempdir, and figuredir full paths

    Parameters
    ----------
    srcdirectory : string
        full path to TFEA source directory
        
    config : dict
        a configparser object that contains variables within the config file

    output : string 
        full path to the parent output directory

    tempdir : string
        full path to the temporary directory where files are stored

    figuredir : string
        full path to the directory containing figures and plots

    Returns
    -------
    None
    '''
    with open(srcdirectory+'config.py','w') as outfile:
        for key in config:
            for item in config[key]:
                outfile.write(item.upper()+'='+config[key][item]+'\n')

        outfile.write('OUTPUTDIR="'+output+'"\n')
        outfile.write('TEMPDIR="'+tempdir+'"\n')
        outfile.write('FIGUREDIR="'+figuredir+'"\n')

        #Path to a directory full of motif logos for all TFs in the HOCOMOCO database (v10)
        logos = os.path.join(os.path.dirname(srcdirectory),'human_logo')
        #logos = srcdirectory + 'human_logo/'
        outfile.write('LOGOS="'+logos+'"\n')

        #Path to mouse directory with motif logos in HOCOMOCO v10
        ##logos = parent_dir(homedir) + '/mouse_logo/'
#==============================================================================

#==============================================================================
def verify_config():
    '''Verifies that all necessary variables are present within the inputted 
        config file and that they are the correct variable types.
    
    Parameters
    ----------
    None

    Returns
    -------
    None

    Raises
    ------
    TypeError
        When variable type does not match required type

    NameError
        When variable is not found within config
    '''
    try:
        if type(config.COMBINE) != bool:
            raise TypeError('COMBINE variable must be a boolean. This switch \
                            determines whether TFEA merges bed files within \
                            BED input.')
    except NameError:
        raise NameError('COMBINE variable not found in config.ini file. This \
                            switch determines whether TFEA merges bed files \
                            within BED input.')

    try:
        if type(config.COUNT) != bool:
            raise TypeError('COUNT variable must be a boolean. This switch \
                                determines whether TFEA performs read \
                                counting over BED regions.')
    except NameError:
        raise NameError('COUNT variable not found in config.ini file. This \
                            switch determines whether TFEA performs read \
                            counting over BED regions.')

    try:
        if type(config.DESEQ) != bool:
            raise TypeError('DESEQ variable must be a boolean. This switch \
                                determines whether TFEA performs DE-Seq \
                                analysis on counted BED regions.')
    except NameError:
        raise NameError('DESEQ variable not found in config.ini file. This \
                            switch determines whether TFEA performs DE-Seq \
                            analysis on counted BED regions.')

    try:
        if type(config.CALCULATE) != bool:
            raise TypeError('CALCULATE variable must be a boolean. This \
                                switch determines whether TFEA performs its \
                                standard enrichment score calculation and \
                                plotting.')
    except NameError:
        raise NameError('CALCULATE variable not found in config.ini file. \
                            This switch determines whether TFEA performs its \
                            standard enrichment score calculation and \
                            plotting.')

    try:
        if type(config.POOL) != bool:
            raise TypeError('POOL variable must be a boolean. This switch \
                                determines whether TFEA runs the analysis \
                                in parallel using the multiprocessing library \
                                in python.')
    except NameError:
        raise NameError('POOL variable not found in config.ini file. This \
                            switch determines whether TFEA runs the analysis \
                            in parallel using the multiprocessing library in \
                            python.')

    try:
        if type(config.SINGLEMOTIF) != bool or type(config.SINGLEMOTIF) != str:
            raise TypeError('SINGLEMOTIF variable must be a boolean or \
                                string. This switch determines whether TFEA \
                                performs its analysis on a single motif or \
                                all. If not False, set to a string matching a \
                                motif name.')
    except NameError:
        raise NameError('SINGLEMOTIF variable not found in config.ini file. \
                            This switch determines whether TFEA performs its \
                            analysis on a single motif or all. If not False, \
                            set to a string matching a motif name.')

    try:
        if type(config.FIMO) != bool:
            raise TypeError('FIMO variable must be a boolean. This switch \
                                determines whether TFEA uses FIMO to get \
                                motif hits or whether a database of motif hit \
                                calls (bed format) is used.')
    except NameError:
        raise NameError('FIMO variable not found in config.ini file. This \
                            switch determines whether TFEA uses FIMO to get \
                            motif hits or whether a database of motif hit \
                            calls (bed format) is used.')

    try:
        if type(config.TEMP) != bool:
            raise TypeError('TEMP variable must be a boolean. This switch \
                                determines whether TFEA saves large temporary \
                                files. If True, temporary files will be \
                                stored in the temp_files directory. Warning: \
                                There will be many large files.')
    except NameError:
        raise NameError('TEMP variable not found in config.ini file. This \
                            switch determines whether TFEA saves large \
                            temporary files. If True, temporary files will be \
                            stored in the temp_files directory. Warning: \
                            There will be many large files.')

    try:
        if type(config.PLOT) != bool:
            raise TypeError('PLOT variable must be a boolean. This switch \
                                determines whether TFEA outputs plots for \
                                all motifs provided or just significant ones. \
                                Warning: Setting this to True will slow down \
                                TFEA and create large output folders.')
    except NameError:
        raise NameError('PLOT variable not found in config.ini file. This \
                            switch determines whether TFEA outputs plots for \
                            all motifs provided or just significant ones. \
                            Warning: Setting this to True will slow down TFEA \
                            and create large output folders.')

    try:
        if type(config.OUTPUT) != str:
            raise TypeError('OUTPUT variable must be a string. Determines \
                                where TFEA stores output.')
    except NameError:
        raise NameError('OUTPUT variable not found in config.ini file. \
                            Determines where TFEA stores output.')

    try:
        if type(config.BEDS) != list:
            raise TypeError('BED variable must be a list. Input a list of \
                                regions (bed-format) to perform analysis \
                                over. If merging not desired, simply input a \
                                single BED file within a python list.')
    except NameError:
        raise NameError('BED variable not found in config.ini file. Input a \
                            list of regions (bed-format) to perform analysis \
                            over. If merging not desired, simply input a \
                            single BED file within a python list.')

    try:
        if type(config.BAM1) != list:
            raise TypeError('BAM1 variable must be a list. Input a list of \
                                BAM files to obtain read depth and coverage \
                                over regions of interest. One or multiple \
                                bams can be specified but they must be within \
                                a python list.')
    except NameError:
        raise NameError('BAM1 variable not found in config.ini file. Input a \
                            list of BAM files to obtain read depth and \
                            coverage over regions of interest. One or \
                            multiple bams can be specified but they must be \
                            within a python list.')

    try:
        if type(config.LABEL1) != str:
            raise TypeError('LABEL1 variable must be a string. Define a \
                                treatment or condition to label BAM1 files. \
                                Used in plotting and output folder naming.')
    except NameError:
        raise NameError('LABEL1 variable not found in config.ini file. Define \
                            a treatment or condition to label BAM1 files. \
                            Used in plotting and output folder naming.')

    try:
        if type(config.BAM2) != list:
            raise TypeError('BAM2 variable must be a list. Input a list of \
                                BAM files to obtain read depth and coverage \
                                over regions of interest. One or multiple \
                                bams can be specified but they must be within \
                                a python list.')
    except NameError:
        raise NameError('BAM2 variable not found in config.ini file. Input a \
                            list of BAM files to obtain read depth and \
                            coverage over regions of interest. One or \
                            multiple bams can be specified but they must be \
                            within a python list.')

    try:
        if type(config.LABEL2) != str:
            raise TypeError('LABEL2 variable must be a string. Define a \
                                treatment or condition to label BAM1 files. \
                                Used in plotting and output folder naming.')
    except NameError:
        raise NameError('LABEL2 variable not found in config.ini file. Define \
                            a treatment or condition to label BAM1 files. \
                            Used in plotting and output folder naming.')

    try:
        if type(config.PADJCUTOFF) != float:
            raise TypeError('PADJCUTOFF variable must be a float. Provide a \
                                p-adjusted cutoff value to determine which \
                                TFs are plotted and called as significant.')
    except NameError:
        raise NameError('PADJCUTOFF variable not found in config.ini file. \
                                Provide a p-adjusted cutoff value to \
                                determine which TFs are plotted and called as \
                                significant.')

    try:
        if type(config.LARGEWINDOW) != float:
            raise TypeError('LARGEWINDOW variable must be a float. Provide a \
                                larger window to be used by TFEA for plotting \
                                and GC-content calculation.')
    except NameError:
        raise NameError('LARGEWINDOW variable not found in config.ini file. \
                                Provide a larger window to be used by TFEA \
                                for plotting and GC-content calculation.')

    try:
        if type(config.SMALLWINDOW) != float:
            raise TypeError('SMALLWINDOW variable must be a float. Provide a \
                                smaller window to be used by TFEA for \
                                plotting.')
    except NameError:
        raise NameError('SMALLWINDOW variable not found in config.ini file. \
                            Provide a smaller window to be used by TFEA for \
                            plotting.')

    try:
        if type(config.MOTIF_HITS) != float:
            raise TypeError('MOTIF_HITS variable must be a string. Provide \
                                the full path to a folder containing bed \
                                files with motif hits across the genome. \
                                These can be generated using TFEAs compile \
                                module.')
    except NameError:
        raise NameError('MOTIF_HITS variable not found in config.ini file. \
                            Provide the full path to a folder containing bed \
                            files with motif hits across the genome. These \
                            can be generated using TFEAs compile module.')

    try:
        if type(config.GENOMEFASTA) != float:
            raise TypeError('GENOMEFASTA variable must be a string. Provide \
                                the full path to a fasta file containing a \
                                genome of interest to perform motif scannning \
                                over.')
    except NameError:
        raise NameError('GENOMEFASTA variable not found in config.ini file. \
                            Provide the full path to a fasta file containing \
                            a genome of interest to perform motif scannning \
                            over.')

    try:
        if type(config.MOTIFDATABASE) != float:
            raise TypeError('MOTIFDATABASE variable must be a string. Provide \
                                the full path to a database of motifs in meme \
                                format.')
    except NameError:
        raise NameError('MOTIFDATABASE variable not found in config.ini file. \
                            Provide the full path to a database of motifs in \
                            meme format.')

    print "Config file verified, all inputs present and correct type."
#==============================================================================

#==============================================================================
def sbatch_submit(srcdirectory='',configpath='',script='',email='',
                    config=None):
    '''Submits an sbatch job using the configuration provided to the script 
        variable that runs TFEA with options specified within a give config 
        file

    Parameters
    ----------
    srcdirectory : string
        full path to TFEA source directory
        
    configpath : string
        full path to the config file

    script : string
        full path to the sbatch script containing configuration options

    email : string
        user e-mail to receive sbatch job notifications

    config : dict
        a configparser object

    Returns
    -------
    None
    '''
    #First make output directories, only save path to the e_and_o folder
    _, _, _, e_and_o = make_out_directories(dirs=True, 
                                                                config=config)

    #Submit the sbatch job
    os.system("sbatch --error=" + e_and_o + "%x.err --output=" + e_and_o 
                + "%x.out --mail-user="+email+" --export=src="+srcdirectory
                + ",config=" +configpath+ " " + script)
#==============================================================================

#==============================================================================
def combine_bed(beds=config.BEDS,tempdir=config.TEMPDIR):
    '''Concatenates, sorts, and merges (bedtools) a list of bed files. Outputs 
        into the tempdir directory created by TFEA

    Parameters
    ----------
    beds : list or array
        full paths to bed files (strings)
        
    tempdir : string
        full path to tempdir directory in output directory (created by TFEA)

    Returns
    -------
    combined_input_merged_bed : string 
        full path to a bed file containing the merged regions inputted by the 
        user 
    '''
    os.system("cat " + " ".join(beds) + " > " + tempdir + "combined_input.bed")

    os.system("sort -k1,1 -k2,2n " + tempdir + "combined_input.bed > " 
                + tempdir + "combined_input.sorted.bed")

    os.system("bedtools merge -i " + tempdir + "combined_input.sorted.bed > " 
                + tempdir + "combined_input.merge.bed")

    combined_input_merged_bed = tempdir + "combined_input.merge.bed"

    return combined_input_merged_bed
#==============================================================================

#==============================================================================
def getfasta(bedfile='', genomefasta=config.GENOMEFASTA, 
                tempdir=config.TEMPDIR):
    '''Converts a bed file to a fasta file using bedtools. Outputs into the 
        tempdir directory created by TFEA.

    Parameters
    ----------
    bedfile : string
        full path to a bed file

    genomefasta : string
        full path to a fasta file for the genome of interest
        
    tempdir : string
        full path to temp directory in output directory (created by TFEA)

    Returns
    -------
    ranked_file_fasta : string 
        full path to a fasta file containing the inputted bed file regions in 
        fasta format 
    '''
    os.system("bedtools getfasta -name -fi "+genomefasta+" -bed "+bedfile
                + " -fo " + tempdir + "ranked_file.fullregions.fa")

    ranked_file_fasta = tempdir + "ranked_file.fullregions.fa"

    return ranked_file_fasta
#==============================================================================

#==============================================================================
def get_bgfile(fastafile='', tempdir=config.TEMPDIR):
    '''Obtains a zero order markov background model (used in FIMO) from a fasta
        file. Outputs into the tempdir directory created by TFEA.

    Parameters
    ----------
    bedfile : string
        full path to a bed file

    genomefasta : string
        full path to a fasta file for the genome of interest
        
    tempdir : string
        full path to temp directory in output directory (created by TFEA)

    Returns
    -------
    markov_background : string 
        full path to a meme-formatted zero order markov background file 
        (txt file)
    '''
    os.system("fasta-get-markov " + fastafile + " " + tempdir
                + "markov_background.txt")

    markov_background = tempdir + "markov_background.txt"

    return markov_background
#==============================================================================

#==============================================================================
def get_regions(tempdir=config.TEMPDIR, 
                ranked_center_file=config.RANKED_CENTER_FILE,
                largewindow=config.LARGEWINDOW):
    '''Takes in a bed file that contains regions centered on the user-inputted 
        bed files and outputs a 'full regions' bed file which simply adds a 
        user-defined window to either side of these centered regions.

    Parameters
    ----------
    tempdir : string
        full path to temp directory in output directory (created by TFEA)

    ranked_center_file : string
        full path to a bed file with centered coordinates
        
    largewindow : float
        half the desired window size

    Returns
    -------
    ranked_full_regions : string 
        full path to a bed file with full regions to be compared with TFEA
    '''
    outfile = open(tempdir+'ranked_file.fullregions.bed','w')
    with open(ranked_center_file) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            chrom,start,stop = line[:3]
            start = str(int(start)-int(largewindow))
            stop = str(int(stop)+int(largewindow))
            pval,fc,rank = line[3:]
            name = ','.join([rank,pval,fc])
            outfile.write('\t'.join([chrom,start,stop,name]) + '\n')

    ranked_full_regions = tempdir + 'ranked_file.fullregions.bed'
    return ranked_full_regions
#==============================================================================

#==============================================================================
def count_reads(bedfile='', bam1=config.BAM1, bam2=config.BAM2, 
                tempdir=config.TEMPDIR, label1=config.LABEL1, 
                label2=config.LABEL2):
    '''Counts reads across regions in a given bed file using bam files inputted
        by a user

    Parameters
    ----------
    bedfile : string
        full path to a bed file containing full regions of interest which will 
        be counted using bedtools multicov

    bam1 : list or array
        a list of full paths to bam files pertaining to a single condition 
        (i.e. replicates of a single treatment)

    bam2 : list or array
        a list of full paths to bam files pertaining to a single condition 
        (i.e. replicates of a single treatment)

    tempdir : string
        full path to temp directory in output directory (created by TFEA)

    label1 : string
        the name of the treatment or condition corresponding to bam1 list

    label2 : string
        the name of the treatment or condition corresponding to bam2 list

    Returns
    -------
    None
    '''
    #This os.system call runs bedtools multicov to count reads in all specified
    #BAMs for given regions in BED
    os.system("bedtools multicov -bams " + " ".join(bam1) + " " 
                + " ".join(bam2) + " -bed " + bedfile + " > " + tempdir 
                + "count_file.bed")

    #This section adds a header to the count_file and reformats it to remove 
    #excess information and add a column with the region for later use
    outfile = open(tempdir + "count_file.header.bed",'w')
    outfile.write("chrom\tstart\tstop\tregion\t" 
                    + '\t'.join([label1]*len(bam1)) + "\t" 
                    + '\t'.join([label2]*len(bam2)) + "\n")

    with open(tempdir + "count_file.bed") as F:
        for line in F:
            line = line.strip('\n').split('\t')
            chrom,start,stop = line[:3]
            counts = line[-(len(bam1)+len(bam2)):]
            outfile.write('\t'.join([chrom,start,stop]) + "\t" 
                            + chrom + ":" + start + "-" + stop + "\t"
                            + '\t'.join(counts) + "\n")
#==============================================================================

#==============================================================================
def write_deseq_script(bam1=config.BAM1, bam2=config.BAM2, 
                        tempdir=config.TEMPDIR, count_file=config.COUNT_FILE,
                        label1=config.LABEL1, label2=config.LABEL2):
    '''Writes an R script within the tempdir directory in TFEA output to run 
        either DE-Seq or DE-Seq2 depending on the number of user-inputted 
        replicates.

    Parameters
    ----------
    bedfile : string
        full path to a bed file containing full regions of interest which will
        be counted using bedtools multicov

    bam1 : list or array
        a list of full paths to bam files pertaining to a single condition 
        (i.e. replicates of a single treatment)

    bam2 : list or array
        a list of full paths to bam files pertaining to a single condition 
        (i.e. replicates of a single treatment)

    tempdir : string
        full path to temp directory in output directory (created by TFEA)

    label1 : string
        the name of the treatment or condition corresponding to bam1 list

    label2 : string
        the name of the treatment or condition corresponding to bam2 list

    Returns
    -------
    None
    '''
    #If more than 1 replicate, use DE-Seq2
    if (len(bam1) > 1 and len(bam2) > 1):
        with open(tempdir + 'DESeq.R','w') as outfile:
            outfile.write('''#!/usr/bin/env Rscript
sink("'''+tempdir+'''DESeq.Rout")
library("DESeq2")
'data <- read.delim("'''+count_file+'''", sep="\t", header=TRUE)
countsTable <- subset(data, \
                select=c('''\
                +', '.join([str(i) for i in range(5,5+len(bam1)+len(bam2))])\
                +'''))

rownames(countsTable) <- data$region
conds <- as.data.frame(c(''' + ', '.join(['"'+label1+'"']*len(bam1)) \
                        + ', ' + ', '.join(['"'+label2+'"']*len(bam2)) \
                        + '''))

colnames(conds) <- c("treatment")
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = countsTable, \
                                            colData = conds, \
                                            design = ~ treatment)

dds <- DESeq(ddsFullCountTable)
res1 <- results(dds,alpha = 0.05, contrast=c("treatment","'''+label2+'''",\
                                            "'''+label1+'''"))

resShrink <- lfcShrink(dds, res = res1, contrast = c("treatment","'''+label2\
                                                    +'''","'''+label1+'''"))

resShrink$fc <- 2^(resShrink$log2FoldChange)
res <- resShrink[c(1:3,7,4:6)]
write.table(res, file = "'''+tempdir+'''DESeq.res.txt", append = FALSE, \
            sep= "\t" )
sink()''')
    #else, there must only be 1 repliceate, use DE-Seq
    else:
        with open(tempdir + 'DESeq.R','w') as outfile:
            outfile.write('''#!/usr/bin/env Rscript
sink("'''+tempdir+'''DESeq.Rout")
library("DESeq")
data <- read.delim("'''+count_file+'''", sep="\t", header=TRUE)
countsTable <- subset(data, \
                select=c('''\
                +', '.join([str(i) for i in range(5,5+len(bam1)+len(bam2))])\
                +'''))

rownames(countsTable) <- data$region
conds <- c(''' + ', '.join(['"'+label1+'"']*len(bam1)) + ', ' \
                + ', '.join(['"'+label2+'"']*len(bam2)) + ''')
cds <- newCountDataSet( countsTable, conds )
cds <- estimateSizeFactors( cds )
sizeFactors(cds)
cds <- estimateDispersions( cds ,method="blind", sharingMode="fit-only")
res <- nbinomTest( cds, "'''+label1+'''", "'''+label2+'''" )
rownames(res) <- res$id
write.table(res, file="'''+tempdir+'''DESeq.res.txt", append=FALSE, sep="\t" )
sink()''')
#==============================================================================

#==============================================================================
def plot_deseq_MA(deseq_file='', label1=config.LABEL1, label2=config.LABEL2,
                    figuredir=config.FIGUREDIR):
    '''Plots the DE-Seq MA-plot using the full regions of interest and saves it
    to the figuredir directory created in TFEA output folder

    Parameters
    ----------
    deseqfile : string
        full path to the deseq file (specifically .res.txt)

    label1 : string
        the name of the treatment or condition corresponding to bam1 list

    label2 : string
        the name of the treatment or condition corresponding to bam2 list

    figuredir : string
        full path to figure directory in output directory (created by TFEA)

    Returns
    -------
    None
    '''
    up_x = list()
    up_y = list()
    up_p = list()
    dn_x = list()
    dn_y = list()
    dn_p = list()
    with open(deseq_file,'r') as F:
        header = F.readline().strip('\n').split('\t')
        basemean_index = header.index('"baseMean"')
        log2fc_index = header.index('"log2FoldChange"')
        for line in F:
            line = line.strip('\n').split('\t')
            try:
                log2fc = float(line[log2fc_index+1])
                basemean = math.log(float(line[basemean_index+1]),10)
                pval = float(line[-2])
                if log2fc > 0:
                    up_x.append(basemean)
                    up_y.append(log2fc)
                    up_p.append(pval)
                else:
                    dn_x.append(basemean)
                    dn_y.append(log2fc)
                    dn_p.append(pval)
            except:
                pass

    x = [x for _,x in sorted(zip(up_p,up_x))] \
        + [x for _,x in sorted(zip(dn_p,dn_x),reverse=True)]

    y = [y for _,y in sorted(zip(up_p,up_y))] \
        + [y for _,y in sorted(zip(dn_p,dn_y),reverse=True)]

    c = plt.cm.RdYlGn(np.linspace(0, 1, len(x)))
    #Creates an MA-Plot of the region expression
    F = plt.figure(figsize=(7,6))
    ax = plt.subplot(111)
    plt.scatter(x=x,y=y,color=c,edgecolor='')
    ax.set_title("DE-Seq MA-Plot",fontsize=14)
    ax.set_ylabel("Log2 Fold-Change ("+label2+"/"+label1+")",fontsize=14)
    ax.set_xlabel("Log10 Average Expression",fontsize=14)
    ax.tick_params(axis='y', which='both', left='off', right='off', 
                    labelleft='on')
    ax.tick_params(axis='x', which='both', bottom='off', top='off', 
                    labelbottom='on')
    plt.savefig(figuredir + 'DESEQ_MA_Plot.png',bbox_inches='tight')
#==============================================================================

#==============================================================================
def permutations(distances=list(), permutations=1000):
    '''Generates permutations of the distances and calculates AUC for each 
        permutation.

    Parameters
    ----------
    distances : list or array
        normalized distances 
        
    permutations : int
        number of times to permute (default=1000)
        
    Returns
    -------
    es_permute : list 
        list of AUC calculated for permutations 
       
    '''
    es_permute = []
    triangle_area = 0.5*(len(distances))
    for i in range(permutations):
        random_distances = np.random.permutation(distances)
        cum_distances = np.cumsum(random_distances)
        es = np.trapz(cum_distances)
        auc = es - triangle_area
        es_permute.append(auc)

    return es_permute
#==============================================================================

#==============================================================================
def padj_bonferroni(TFresults=list()):
    '''This function iterates through TFEA results, removes TFs that returned 
        "no hits" and calculates a p-adj using the Bonferroni Correction for 
        each TF motif appending it to the given TFresults array

    Parameters
    ----------
    TFresults : list of lists
        contains calculated enrichment scores for all TFs of interest specified
        by the user
        
    Returns
    -------
    TFresults : list of lists
        same as input with an additional p-adjusted value appended to each TF
    '''
    TFresults = [x for x in TFresults if x != "no hits"]
    for i in range(len(TFresults)):
        PVAL = TFresults[i][-1]
        #Using Bonferroni Correction
        PADJ = 1 if PVAL*len(TFresults) > 1 else PVAL*len(TFresults)
        TFresults[i].append(PADJ)

    return TFresults
#==============================================================================

#==============================================================================
def motif_distance_bedtools_closest(ranked_center_file='',
                            motif_path=''):
    '''Calculates nearest motif hit from a bed file. TFEA provides this 
        function with a bed file containing the center of the inputted regions.

    Parameters
    ----------
    TFresults : list of lists
        contains calculated enrichment scores for all TFs of interest specified
        by the user
        
    Returns
    -------
    motif_distance_bed_sorted : string
        full path to where the sorted motif distance file was outputted
    '''
    os.system("bedtools closest -D ref -t first -a " 
                + ranked_center_file.split('.bed')[0] + ".sorted.bed -b " 
                + motif_path + " > " + '/' 
                + '/'.join(ranked_center_file.split('/')[:-1]) + '/' 
                + motif_path.split('/')[-1] + ".sorted.distance.bed")

    motif_distance_bed_sorted = '/' \
                                + '/'.join(ranked_center_file.split('/')[:-1])\
                                + '/' + motif_path.split('/')[-1]\
                                + ".sorted.distance.bed"

    return motif_distance_bed_sorted
#==============================================================================

#==============================================================================
def fimo(tempdir=config.TEMPDIR, motifdatabase=config.MOTIFDATABASE, bgfile='',
            motif='', fastafile=''):
    '''This function runs fimo on a given fastafile for a single motif in a 
        provided motif database. The output is cut and sorted to convert into 
        a sorted bed file

    Parameters
    ----------
    tempdir : string
        full path to temp directory in output directory (created by TFEA)

    motifdatabase : string
        full path to a motif database file in meme format

    bgfile : string
        full path to a markov background model

    motif : string
        the name of a motif that matches a motif within motifdatabase

    fastafile : string
        full path to a fasta file that fimo will perform motif scanning on
        
    Returns
    -------
    fimo_out : string
        full path to where fimo output which is stored within the tempdir 
        directory.
    '''
    os.system("fimo --text --bgfile "+bgfile+" --motif "+motif+" " +
              motifdatabase+" "+fastafile+" > "+tempdir+motif+".txt")
    os.system("cut -d$'\\t' -f 3- "+tempdir+motif +
              ".txt > " + tempdir+motif+".bed")

    fimo_out = tempdir+motif+".bed"

    return fimo_out
#==============================================================================

#==============================================================================
def meme2images(motifdatabase=config.MOTIFDATABASE, figuredir=config.FIGUREDIR,
                motif=''):
    '''This function creates meme logos for use in the output html

    Parameters
    ----------
    motifdatabase : string
        full path to a motif database file in meme format

    motif : string
        the name of a motif that matches a motif within motifdatabase

    figuredir : string
        full path to figure directory in output directory (created by TFEA)
        
    Returns
    -------
    None
    '''
    os.system("meme2images -png -rc --motif " + motif + " " + motifdatabase
              + " " + figuredir)
#==============================================================================

#==============================================================================
def fasta_markov(tempdir=config.TEMPDIR, fastafile='', order='0'):
    '''This function runs meme's fasta-get-markov function that generates a 
        background markov file (for use with fimo) from a fasta file.
    Parameters
    ----------
    tempdir : string
        full path to temp directory in output directory (created by TFEA)

    fastafile : string
        full path to fasta file that will be used to generate the markov
        background model file

    order : string
        an integer formatted as a string where a user may specify what order
        markov model they would like (default='0')
        
    Returns
    -------
    None
    '''
    os.system("fasta-get-markov -m " + order + " " + fastafile +
              " > " + tempdir + "Markov_Background.txt")
#==============================================================================

#==============================================================================
def fimo_distance(largewindow=config.LARGEWINDOW, tempdir=config.TEMPDIR, 
                    fimo_file='', motif_file=''):
    '''Parses a fimo output file and writes into a new file that is formatted
        in a way that can be parsed within existing TFEA functions
    Parameters
    ----------
    largewindow : float
        the size of the larger window to perform TFEA. Specified by user in
        config file

    tempdir : string
        full path to temp directory in output directory (created by TFEA)

    fimo_file : string
        full path to the fimo output file to be parsed by this function

    motif_file : string
        the name of the motif being parsed, this function will create a file
        using this motif_file string
        
    Returns
    -------
    outname : string
        the full path to the file to be used by other TFEA functions
    '''
    d = dict()
    with open(fimo_file) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            pval, fc, rank = line[0].split(',')
            start, stop = line[1:2]
            distance = ((int(start)+int(stop))/2)-largewindow
            score = line[5]
            if rank not in d:
                d[rank] = [rank, pval, fc, score, distance]
            else:
                prev_distance = math.fabs(d[rank][-1])
                if prev_distance > math.fabs(distance):
                    d[rank] = [rank, pval, fc, score, distance]

    outname = tempdir+motif_file+'.sorted.distance.bed'
    outfile = open(outname, 'w')
    for i in range(len(d)):
        rank = str(i)
        outfile.write('\t'.join(d[rank])+'\n')
    outfile.close()

    return outname
#==============================================================================

#==============================================================================
def convert_sequence_to_array(sequence=''):
    '''This function converts a DNA sequence (ACGT alphabet) to an array, collapsing GCs
    to 1's and ATs to 0's

    Parameters
    ----------
    sequence : string
        a string containing ACGT character
        
    Returns
    -------
    array : list
        a list of floats corresponding to 1.0 or 0.0 depending on whether the 
        input sequence was GC or AT at a specific site

    Raises
    ------
    Warning : str
        when a character is not ACTG. Value is given 0.0
    '''
    array = []
    for character in sequence:
        character = character.upper()
        if character == 'G' or character == 'C':
            array.append(1.0)
        elif character == 'A' or character == 'T':
            array.append(0.0)
        else:
            print "Warning: Character not recognized"
            array.append(0.0)

    return array 
#==============================================================================

#==============================================================================
def rank_deseqfile(deseq_file=config.DESEQ_FILE, tempdir=config.TEMPDIR):
    '''This function parses a DE-seq output file and creates a new file with 
        the center of each region ranked by p-value 
    
    Parameters
    ----------
    deseq_file : string
        full path to a DE-Seq output file
    
    tempdir : string
        full path to the tempdir directory in the output directory (created by 
        TFEA)
        
    Returns
    -------
    None
    '''
    up = list()
    down = list()
    with open(deseq_file) as F:
        header = F.readline().strip('\n').split('\t')
        fc_index = [i for i in range(len(header)) if header[i]=='"fc"' or header[i]=='"foldChange"'][0]
        for line in F:
            line = line.strip('\n').split('\t')
            if line[fc_index] != 'NA':
                try:
                    pval = format(float(line[-2]),'.12f')
                except ValueError:
                    pval = format(1.0,'.12f')
                region = line[0].split(':')
                chrom = region[0]
                coordinates = region[1].split('-')
                start = coordinates[0]
                stop = coordinates[1]
                chrom = chrom.strip('"')
                stop = stop.strip('"')
                fc = float(line[fc_index+1])
                if fc < 1:
                    down.append((chrom,start,stop,pval,str(fc)))
                else:
                    up.append((chrom,start,stop,pval,str(fc)))

    #Save ranked regions in a bed file (pvalue included)
    outfile = open(tempdir + "ranked_file.bed",'w')
    r=1
    for region in sorted(up, key=lambda x: x[3]):
        outfile.write('\t'.join(region) + '\t' + str(r) + '\n')
        r += 1
    for region in sorted(down, key=lambda x: x[3], reverse=True):
        outfile.write('\t'.join(region) + '\t' + str(r) + '\n')
        r += 1
    outfile.close()

    #Get center base for each region
    outfile = open(tempdir+"ranked_file.center.bed",'w')
    with open(tempdir + "ranked_file.bed") as F:
        for line in F:
            line = line.strip('\n').split('\t')
            chrom,start,stop = line[:3]
            center = (int(start)+int(stop))/2
            outfile.write(chrom + '\t' + str(center) + '\t' + str(center+1) 
                            + '\t' + '\t'.join(line[3:]) + '\n')
    outfile.close()


    os.system("sort -k1,1 -k2,2n " + tempdir+"ranked_file.center.bed" + " > " 
                + tempdir + "ranked_file.center.sorted.bed")
#==============================================================================

#==============================================================================
def samtools_flagstat(args, tempdir=config.TEMPDIR):
    '''Performs samtools flagstat on a bam file. Then parses the samtools 
        flagstat output and returns a the millions mapped reads for the given
        bam file.

    Parameters
    ----------
    args : string
        full paths to a bam files

    tempdir : string
        full path to the tempdir directory within the output directory (created
        by TFEA)
        
    Returns
    -------
    millions_mapped : float
        millions mapped reads for the given bam file
    '''
    bam = args
    filename = bam.split('/')[-1]
    os.system("samtools flagstat " + bam + " > " + tempdir + filename 
                + ".flagstat")
    with open(tempdir+filename+".flagstat") as F:
        lines = F.readlines()
        millions_mapped = float(lines[4].strip('\n').split(' ')[0])/1000000.0

        return millions_mapped