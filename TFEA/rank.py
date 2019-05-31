#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''This module takes as input a count file produced from the COUNT module or 
    user-defined and ranks regions based on some metric. This script outputs 
    files necessary to use with the SCANNER module
'''

#==============================================================================
__author__ = 'Jonathan D. Rubin and Rutendo F. Sigauke'
__credits__ = ['Jonathan D. Rubin', 'Rutendo F. Sigauke', 'Jacob T. Stanley',
                'Robin D. Dowell']
__maintainer__ = 'Jonathan D. Rubin'
__email__ = 'Jonathan.Rubin@colorado.edu'

#Imports
#==============================================================================
import os
import sys
import time
import datetime
import subprocess
from pathlib import Path

from pybedtools import BedTool
import HTSeq as hts
import numpy as np

from TFEA import exceptions
from TFEA import multiprocess

#Main Script
#==============================================================================
def main(use_config=True, combined_file=None, rank=None, scanner=None, 
            bam1=None, bam2=None, tempdir=None, label1=None, label2=None, 
            largewindow=None, mdd=None, mdd_bedfile1=None, mdd_bedfile2=None, 
            debug=False, jobid=None):
    '''This is the main script of the RANK module which takes as input a
        count file and bam files and ranks the regions within the count file
        according to a user specified 

    Parameters
    ----------
    use_config : boolean
        Whether to use a config module to assign variables.
    count_file : str
        Full path to a count file. Can be generated with the COUNT module or
        defined by a user
    scanner : str
        Type of scanning to be used in the SCANNER module. If set to 
        'genome hits', then return the center of regions, otherwise return
        the full window using largewindow as the half-window size.
    bam1 : list
        A list of strings specifying full paths to bam files associated with 
        a single condition (replicates)
    bam2 : list
        A list of strings specifying full paths to bam files associated with 
        a single condition (replicates)
    tempdir : str
        Full path to a directory where files will be saved
    label1 : str
        An informative label for the condition corresponding to files within
        bam1. Used to add a header line to the count file for downstream
        DE-Seq analysis.
    label2 : str
        An informative label for the condition corresponding to files within
        bam2. Used to add a header line to the count file for downstream
        DE-Seq analysis.
    largewindow : int
        Half-length of window size to use when generating md-score related
        bed files
    mdd : boolean
        a switch which determines whether to create bed files for mdd analysis

    Returns
    -------
    ranked_file : str
        Full path to a ranked file that contains regions ranked based on 
        DE-Seq p-value
    mdd_bedfile1 : str
        Full path to a bed file with regions for mdd analysis
    mdd_bedfile2 : str
        Full path to a bed file with regions for mdd analysis

    Raises
    ------
    FileEmptyError
        If any resulting file is empty
    '''
    start_time = time.time()
    if use_config:
        from TFEA import config
        combined_file=config.vars.COMBINED_FILE
        rank=config.vars.RANK
        scanner=config.vars.SCANNER
        bam1=config.vars.BAM1
        bam2=config.vars.BAM2
        tempdir=config.vars.TEMPDIR
        label1=config.vars.LABEL1
        label2=config.vars.LABEL2
        largewindow=config.vars.LARGEWINDOW
        mdd=config.vars.MDD 
        mdd_bedfile1=config.vars.MDD_BEDFILE1
        mdd_bedfile2=config.vars.MDD_BEDFILE2
        mdd_pval=config.vars.MDD_PVAL
        mdd_percent=config.vars.MDD_PERCENT
        debug = config.vars.DEBUG
        jobid = config.vars.JOBID
    print("Ranking regions...", end=' ', flush=True, file=sys.stderr)

    #Begin by counting reads from bam files over the combined_file produced
    # by the combine module
    count_file = count_reads(bedfile=combined_file, bam1=bam1, bam2=bam2, 
                            tempdir=tempdir, label1=label1, label2=label2)
    sample_number = (len(bam1)+len(bam2))
    millions_mapped = sum_reads(count_file=count_file, 
                                sample_number=sample_number)
    
    
    if os.stat(count_file).st_size == 0:
        raise exceptions.FileEmptyError("Error in RANK module. Counting failed.")

    if rank == 'deseq' or rank == 'fc':
        ranked_file, pvals, fcs = deseq(bam1=bam1, bam2=bam2, tempdir=tempdir, 
                                    count_file=count_file, label1=label1, 
                                    label2=label2, largewindow=largewindow, 
                                    rank=rank)
        q1regions, q2regions, q3regions, q4regions = quartile_split(ranked_file)
        meta_profile = meta_profile_quartiles(q1regions, q2regions, q3regions, q4regions, 
                            bam1=bam1, bam2=bam2, millions_mapped=millions_mapped, 
                            largewindow=largewindow)
    else:
        raise exceptions.InputError("RANK option not recognized.")
    if os.stat(ranked_file).st_size == 0:
        raise exceptions.FileEmptyError("Error in RANK module. DE-Seq running or parsing failed.")

    if use_config:
        config.vars.RANKED_FILE = ranked_file
        config.vars.PVALS = pvals
        config.vars.FCS = fcs
        config.vars.META_PROFILE = meta_profile

    if mdd and (not mdd_bedfile1 or not mdd_bedfile2):
        mdd_bedfile1, mdd_bedfile2 = create_mdd_files(ranked_file=ranked_file,
                                                        tempdir=tempdir,
                                                        percent=mdd_percent,
                                                        pval_cut=mdd_pval)
        if os.stat(mdd_bedfile1).st_size == 0 or os.stat(mdd_bedfile2).st_size == 0:
            raise exceptions.FileEmptyError("Error in RANK module. MDD bed file creation failed.")
        if use_config:
            if not config.vars.MDD_BEDFILE1:
                config.vars.MDD_BEDFILE1 = mdd_bedfile1
            if not config.vars.MDD_BEDFILE2:
                config.vars.MDD_BEDFILE2 = mdd_bedfile2
    
    total_time = time.time() - start_time
    if use_config:
        config.vars.RANKtime = total_time
    print("done in: " + str(datetime.timedelta(seconds=int(total_time))), file=sys.stderr)

    if debug:
        multiprocess.current_mem_usage(jobid)

#Functions
#==============================================================================
def count_reads(bedfile=None, bam1=None, bam2=None, tempdir=None, label1=None, 
                label2=None):
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
    count_file = tempdir / "count_file.bed"

    #pybedtools implementation (incomplete)
    # pybed_count = BedTool(bedfile.as_posix()).multi_bam_coverage(bams=bam1+bam2)
    # pybed_count.saveas(count_file, trackline=("#chrom\tstart\tstop\tregion\t" 
    #                                             + '\t'.join([label1]*len(bam1)) + "\t" 
    #                                             + '\t'.join([label2]*len(bam2)) + "\n"))

    #Bedtools implementation
    multicov_command = ["bedtools", "multicov", 
                        "-bams"] + bam1+bam2 + ["-bed", bedfile]
    with open(count_file, 'w') as outfile:
        subprocess.run(multicov_command, stdout=outfile, check=True)

    # This section adds a header to the count_file and reformats it to remove 
    # excess information and add a column with the region for later use
    count_file_header = tempdir / "count_file.header.bed"
    outfile = open(count_file_header, 'w')
    outfile.write("#chrom\tstart\tstop\tregion\t" 
                    + '\t'.join([label1]*len(bam1)) + "\t" 
                    + '\t'.join([label2]*len(bam2)) + "\n")

    with open(count_file) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            chrom,start,stop = line[:3]
            counts = line[-(len(bam1)+len(bam2)):]
            outfile.write('\t'.join([chrom,start,stop]) + "\t" 
                            + chrom + ":" + start + "-" + stop + "\t"
                            + '\t'.join(counts) + "\n")
    outfile.close()
    return count_file_header

#==============================================================================
def sum_reads(count_file=None, sample_number=None):
    '''This function calculates millions mapped reads to regions of interest
        to be used later for normalization of meta plot
    
    Parameters
    ----------
    count_file : string
        The full path to a count file containing reads mapping to regions of
        interest
    
    sample_number : int
        The total number of samples

    Returns
    -------
    millions_mapped : list
        A list of ints that corresponds to the number of millions mapped per
        sample in the order that appears in the count_file (by column)
    '''
    millions_mapped = [0.0]*sample_number
    with open(count_file) as F:
        F.readline()
        for line in F:
            line = line.strip('\n').split('\t')
            [x+float(y) for x,y in zip(millions_mapped,line[-1*sample_number:])]

    return millions_mapped

#==============================================================================
def write_deseq_script(bam1=None, bam2=None, tempdir=None, count_file=None, 
                        label1=None, label2=None):
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
        Rfile = open(tempdir / 'DESeq.R','w')
        Rfile.write('''library("DESeq2")
data <- read.delim("'''+count_file.as_posix()+'''", sep="\t", header=TRUE)
countsTable <- subset(data, select=c('''
                +', '.join([str(i) for i in range(5,5+len(bam1)+len(bam2))])
                +'''))

rownames(countsTable) <- data$region
conds <- as.data.frame(c(''' 
                        + ', '.join(['"'+label1+'"']*len(bam1)) 
                        + ', ' 
                        + ', '.join(['"'+label2+'"']*len(bam2)) 
                        + '''))

colnames(conds) <- c("treatment")
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = countsTable,
                                            colData = conds,
                                            design = ~ treatment)

dds <- DESeq(ddsFullCountTable)
res <- results(dds, alpha = 0.05, contrast=c("treatment", "'''+label2+'''",
                                                            "'''+label1+'''"))
res$fc <- 2^(res$log2FoldChange)
res <- res[c(1:3,7,4:6)]

write.table(res, file = "'''    + (tempdir / 'DESeq.res.txt').as_posix() 
                                + '''", append = FALSE, sep= "\t" )
sink()''')
    else:
        Rfile = open(tempdir /  'DESeq.R','w')
        Rfile.write('''library("DESeq")
data <- read.delim("'''+count_file.as_posix()+'''", sep="\t", header=TRUE)

countsTable <- subset(data, select=c('''
            +', '.join([str(i) for i in range(5,5+len(bam1)+len(bam2))])
            +'''))

rownames(countsTable) <- data$region
conds <- c('''  + ', '.join(['"'+label1+'"']*len(bam1)) 
                + ', ' 
                + ', '.join(['"'+label2+'"']*len(bam2)) 
                + ''')

cds <- newCountDataSet( countsTable, conds )
cds <- estimateSizeFactors( cds )
sizeFactors(cds)                                                               
cds <- estimateDispersions( cds ,method="blind", \
                        sharingMode="fit-only")

res <- nbinomTest( cds, "'''+label1+'''", "'''+label2+'''" )
rownames(res) <- res$id                      
write.table(res, file = "'''    + os.path.join(tempdir,'DESeq.res.txt') 
                                + '''", append = FALSE, sep= "\t" )''')
    Rfile.close()

#==============================================================================
def deseq(bam1=None, bam2=None, tempdir=None, count_file=None, label1=None, 
            label2=None, largewindow=None, rank=None):
    #Write the DE-Seq R script
    write_deseq_script(bam1=bam1, bam2=bam2, tempdir=tempdir, 
                        count_file=count_file, label1=label1, label2=label2)

    #Execute the DE-Seq R script
    # with open(tempdir / 'DESeq.Rout', 'w') as stdout:
    deseqR = tempdir / "DESeq.R"
    deseqout = tempdir / 'DESeq.Rout'
    deseq_file = tempdir / 'DESeq.res.txt'
    with open(deseqout, 'w') as output:
        exitcode = subprocess.run(["Rscript", "--vanilla", deseqR], stdout=output, 
                            stderr=output)
    if exitcode.returncode != 0:
        errormessage = deseqout.read_text()
        if 'Error' in errormessage:
            printmessage = errormessage[errormessage.index('Error'):]
        raise exceptions.SubprocessError(printmessage)

    ranked_file = deseq_parse(deseq_file=deseq_file, tempdir=tempdir, 
                                largewindow=largewindow, rank=rank)

    return ranked_file

#==============================================================================
def deseq_parse(deseq_file=None, tempdir=None, largewindow=None, rank=None):
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
    ranked_center_file : string
        full path to a bed file that contains the center of regions of interest
        ranked via DE-Seq p-value
    '''
    #Parse the DE-Seq output file
    up = list()
    down = list()
    with open(deseq_file) as F:
        header = F.readline().strip('\n').split('\t')
        fc_index = [i for i in range(len(header)) 
                    if header[i]=='"fc"' or header[i]=='"foldChange"'][0]
        for line in F:
            line = line.strip('\n').split('\t')
            if line[fc_index+1] != 'NA':
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
                center = int((int(start)+int(stop))/2)
                start = center - int(largewindow)
                stop = center + int(largewindow)
                fc = float(line[fc_index+1])
                if fc < 1:
                    down.append([chrom, str(start), str(stop), str(fc), str(pval)])
                else:
                    up.append([chrom, str(start), str(stop), str(fc), str(pval)])

    #Save ranked regions in a bed file (pvalue included)
    ranked_file = tempdir / "ranked_file.bed"
    pvals = list()
    fcs = list()
    with open(ranked_file,'w') as outfile:
        outfile.write('\t'.join(['#chrom', 'start', 'stop', 'rank, fc, p-value']) 
                        + '\n')
        r=1
        if rank == 'deseq':
            for region in sorted(up, key=lambda x: x[4]):
                outfile.write('\t'.join(region[:3]) 
                            + '\t' + ','.join(region[3:]+[str(r)]) 
                            + '\n')
                pvals.append(float(region[-1]))
                fcs.append(float(region[4]))
                r += 1
            for region in sorted(down, key=lambda x: x[4], reverse=True):
                outfile.write('\t'.join(region[:3]) 
                            + '\t' + ','.join(region[3:]+[str(r)]) 
                            + '\n')
                pvals.append(float(region[-1]))
                fcs.append(float(region[4]))
                r += 1
        elif rank == 'fc':
            for region in sorted(up+down, key=lambda x: x[3], reverse=True):
                outfile.write('\t'.join(region[:3]) 
                            + '\t' + ','.join(region[3:]+[str(r)]) 
                            + '\n')
                pvals.append(float(region[-1]))
                fcs.append(float(region[4]))
                r += 1

    return ranked_file, pvals, fcs

#==============================================================================
def create_mdd_files(ranked_file=None, percent=False, pval_cut=False, tempdir=None):
    '''This function creates bed files to be used with motif displacement
        differential (mdd) calculation. It separates a 'ranked_file' which
        contains regions ranked by differential expression via DE-Seq

    Parameters
    ----------
    ranked_file : str
        a file containing regions ranked on p-value obtained from DE-Seq
    percent : float
        if not 'None', percentage cutoff used to separate bed files
    pval : float
        if not 'None', pval cutoff value used to separate bed files
    tempdir : str
        full path to a temporary directory to store files

    Returns
    -------
    mdd_bedfile1 : str
        full path to a file containing bed regions associated with 
        non-differentially transcribed regions
    mdd_bedfile2 : str
        full path to a file containing bed regions associated with 
        differentially transcribed regions
    '''
    mdd_bedfile1 = tempdir / 'mdd_bedfile1.bed'
    mdd_bedfile2 = tempdir / 'mdd_bedfile2.bed'
    regions = list()
    with open(ranked_file) as F:
        for line in F:
            if line[0] != '#':
                chrom, start, stop, rank_values = line.strip('\n').split('\t')
                fc, pval, rank = rank_values.split(',')
                fc = float(fc)
                pval = float(pval)
                rank = int(rank)
                regions.append((chrom, start, stop, fc, pval, rank))
    regions = sorted(regions, key=lambda x: x[-2])
    if percent != False:
        index = int(len(regions) * percent)
        mdd1_regions = regions[index:]
        mdd2_regions = regions[:index]
    elif pval_cut != False:
        mdd1_regions = [region for region in regions if region[-2] >= pval_cut]
        mdd2_regions = [region for region in regions if region[-2] < pval_cut]
    else:
        raise exceptions.InputError("No cutoff value specified for creating mdd bed files")

    with open(mdd_bedfile1, 'w') as ofile1:
        ofile1.write('\t'.join(['#chrom', 'start', 'stop', 'fc,pval,oldrank', 'rank']) + '\n')
        rank1 = 1
        for region1 in mdd1_regions:
            ofile1.write('\t'.join([x for x in region1[:3]]) 
                        + '\t' + ','.join([str(x) for x in region1[3:-1]] + [str(rank1)])
                        + '\n')
            rank1 += 1

    with open(mdd_bedfile2, 'w') as ofile2:
        ofile2.write('\t'.join(['#chrom', 'start', 'stop', 'fc,pval,oldrank', 'rank']) + '\n')
        rank2 = 1
        for region2 in mdd2_regions:
            ofile2.write('\t'.join([x for x in region2[:3]]) 
                        + '\t' + ','.join([str(x) for x in region2[3:-1]] + [str(rank2)])
                        + '\n')
            rank2 += 1

    return mdd_bedfile1, mdd_bedfile2

#==============================================================================
def quartile_split(ranked_file):
    '''Takes a ranked_file and outputs regions separated by quartiles
    '''
    regions = list()
    with open(ranked_file) as F:
        F.readline()
        for line in F:
            linelist = line.strip('\n').split('\t')
            regions.append(linelist[:3])
    
    q1 = int(round(len(regions)*.25))
    q2 = int(round(len(regions)*.5))
    q3 = int(round(len(regions)*.75))

    q1regions = regions[:q1]
    q2regions = regions[q1:q2]
    q3regions = regions[q2:q3]
    q4regions = regions[q3:]

    return q1regions, q2regions, q3regions, q4regions
#==============================================================================
def meta_profile_quartiles(q1regions, q2regions, q3regions, q4regions, 
                            bam1=None, bam2=None, millions_mapped=None, 
                            largewindow=None):
    '''This function creates a metaprofile from 4 regions and stores them in
        a dictionary. 
    '''
    q1posprofile1, q1negprofile1, q1posprofile2, q1negprofile2 = meta_profile(
                                        regionlist=q1regions, 
                                        millions_mapped=millions_mapped, 
                                        largewindow=largewindow, bam1=bam1, 
                                        bam2=bam2)
    
    q2posprofile1, q2negprofile1, q2posprofile2, q2negprofile2 = meta_profile(
                                        regionlist=q2regions, 
                                        millions_mapped=millions_mapped, 
                                        largewindow=largewindow, bam1=bam1, 
                                        bam2=bam2)

    q3posprofile1, q3negprofile1, q3posprofile2, q3negprofile2 = meta_profile(
                                        regionlist=q3regions, 
                                        millions_mapped=millions_mapped, 
                                        largewindow=largewindow, bam1=bam1, 
                                        bam2=bam2)

    q4posprofile1, q4negprofile1, q4posprofile2, q4negprofile2 = meta_profile(
                                        regionlist=q4regions, 
                                        millions_mapped=millions_mapped, 
                                        largewindow=largewindow, bam1=bam1, 
                                        bam2=bam2)

    meta_profile_dict = {'q1posprofile1': q1posprofile1, 
                        'q1negprofile1': q1negprofile1,
                        'q1posprofile2': q1posprofile2, 
                        'q1negprofile2': q1negprofile2, 
                        'q2posprofile1': q2posprofile1, 
                        'q2negprofile1': q2negprofile1, 
                        'q2posprofile2': q2posprofile2, 
                        'q2negprofile2': q2negprofile2,
                        'q3posprofile1': q3posprofile1, 
                        'q3negprofile1': q3negprofile1, 
                        'q3posprofile2': q3posprofile2, 
                        'q3negprofile2': q3negprofile2,
                        'q4posprofile1': q4posprofile1, 
                        'q4negprofile1': q4negprofile1, 
                        'q4posprofile2': q4posprofile2, 
                        'q4negprofile2': q4negprofile2}

    return meta_profile_dict

#==============================================================================
def meta_profile(regionlist=None, region_file=None, millions_mapped=None, 
                largewindow=None, bam1=None, bam2=None):
    '''This function returns average profiles for given regions of interest.
        A user may input either a list of regions or a bed file
    Parameters
    ----------
    regionlist : list
        a list of regions of interest. Format [(chrom, start, stop), (), ...]
    
    region_file : string
        full path to a bed file containing regions of interest.
    millions_mapped : list
        a list of floats corresponding to millions mapped reads for inputted
        bam files. These must be in order corresponding to the order of bam1
        files followed by bam2 files.
    
    largewindow : float
        the window with which to compute profiles for
    bam1 : list
        a list of full paths to bam files corresponding to a condition or 
        treatment
    
    bam2 : list
        a list of full paths to bam files corresponding to a condition or 
        treatment
        
    Returns
    -------
    posprofile1 : list
        profile for the positive strand corresponding to condition 1. Each 
        value in the list is a bp
    negprofile1 : list
        profile for the negative strand corresponding to condition 1. Each 
        value in the list is a bp
    posprofile2 : list
        profile for the positive strand corresponding to condition 2. Each 
        value in the list is a bp
    negprofile2 : list
        profile for the negative strand corresponding to condition 2. Each 
        value in the list is a bp
    Raises
    ------
    '''
    regions=list()
    if regionlist != None and region_file == None:
        for chrom, start, stop in regionlist:
            regions.append(
                hts.GenomicInterval(chrom, int(start), int(stop), '.'))
    elif region_file != None and regionlist == None:
        with open(region_file) as F:
            for line in F:
                line = line.strip('\n').split('\t')
                chrom, start, stop = line[:3]
                regions.append(
                    hts.GenomicInterval(chrom, int(start), int(stop), '.'))
    elif region_file == None and regionlist == None:
        raise NameError("Must input either a regionlist or region_file.")
    else:
        raise NameError("Only input one of regionlist or region_file variables.")
    

    hts_bam1 = list()
    hts_bam2 = list()
    for bam in bam1:
        hts_bam1.append(hts.BAM_Reader(bam))
    for bam in bam2:
        hts_bam2.append(hts.BAM_Reader(bam))

    if len(hts_bam1) == 0 or len(hts_bam2) == 0:
        raise ValueError("One of bam1 or bam2 variables is empty.")
    # if len(millions_mapped) == 0 or len(millions_mapped) != len(hts_bam1) + len(hts_bam2):
    #     raise ValueError(("Millions_mapped variable is empty or does not match "
    #                         "length of bam variables combined."))


    posprofile1 = np.zeros(2*int(largewindow))  
    negprofile1 = np.zeros(2*int(largewindow))
    posprofile2 = np.zeros(2*int(largewindow))   
    negprofile2 = np.zeros(2*int(largewindow))
    rep1number = float(len(hts_bam1))
    rep2number = float(len(hts_bam2))
    mil_map1 = 0.0
    mil_map2 = 0.0
    for window in regions:
        avgposprofile1 = np.zeros(2*int(largewindow))
        avgnegprofile1 = np.zeros(2*int(largewindow))
        i = 0
        for sortedbamfile in hts_bam1:
            # mil_map = float(millions_mapped[i])
            i += 1
            tempposprofile = np.zeros(2*int(largewindow))
            tempnegprofile = np.zeros(2*int(largewindow))
            for almnt in sortedbamfile[ window ]:
                mil_map1 += 1.0
                if almnt.iv.strand == '+':
                    start_in_window = almnt.iv.start - window.start
                    end_in_window = almnt.iv.end - window.end \
                                    + 2*int(largewindow)
                    start_in_window = max( start_in_window, 0 )
                    end_in_window = min( end_in_window, 2*int(largewindow) )
                    tempposprofile[ start_in_window : end_in_window ] += 1.0
                if almnt.iv.strand == '-':
                    start_in_window = almnt.iv.start - window.start
                    end_in_window   = almnt.iv.end - window.end \
                                        + 2*int(largewindow)
                    start_in_window = max( start_in_window, 0 )
                    end_in_window = min( end_in_window, 2*int(largewindow) )
                    tempnegprofile[ start_in_window : end_in_window ] += -1.0
            # pos_sum = np.sum(tempposprofile)
            # neg_sum = np.sum(tempnegprofile)
            # if pos_sum != 0:
            #     tempposprofile = [x/pos_sum for x in tempposprofile]
            # if neg_sum != 0:
            #     tempnegprofile = [-(x/neg_sum) for x in tempnegprofile]
            avgposprofile1 = [x+y for x,y in zip(avgposprofile1, tempposprofile)]
            avgnegprofile1 = [x+y for x,y in zip(avgnegprofile1, tempnegprofile)]
        # avgposprofile1 = [x/rep1number/mil_map for x in avgposprofile1]
        # avgnegprofile1 = [x/rep1number/mil_map for x in avgnegprofile1]
        avgposprofile1 = [x/rep1number for x in avgposprofile1]
        avgnegprofile1 = [x/rep1number for x in avgnegprofile1]
        posprofile1 = [x+y for x,y in zip(posprofile1,avgposprofile1)]
        negprofile1 = [x+y for x,y in zip(negprofile1, avgnegprofile1)]

        avgposprofile2 = np.zeros(2*int(largewindow))
        avgnegprofile2 = np.zeros(2*int(largewindow))
        i = len(hts_bam1)
        for sortedbamfile in hts_bam2:
            # mil_map = float(millions_mapped[i])
            i += 1
            tempposprofile = np.zeros(2*int(largewindow))
            tempnegprofile = np.zeros(2*int(largewindow))
            for almnt in sortedbamfile[ window ]:
                if almnt.iv.strand == '+':
                    mil_map2 += 1.0
                    start_in_window = almnt.iv.start - window.start
                    end_in_window   = almnt.iv.end - window.end \
                                        + 2*int(largewindow)
                    start_in_window = max( start_in_window, 0 )
                    end_in_window = min( end_in_window, 2*int(largewindow) )
                    tempposprofile[ start_in_window : end_in_window ] += 1.0
                if almnt.iv.strand == '-':
                    start_in_window = almnt.iv.start - window.start
                    end_in_window   = almnt.iv.end - window.end \
                                        + 2*int(largewindow)
                    start_in_window = max( start_in_window, 0 )
                    end_in_window = min( end_in_window, 2*int(largewindow) )
                    tempnegprofile[ start_in_window : end_in_window ] += -1.0
            # pos_sum = np.sum(tempposprofile)
            # neg_sum = np.sum(tempnegprofile)
            # if pos_sum != 0:
            #     tempposprofile = [x/pos_sum for x in tempposprofile]
            # if neg_sum != 0:
            #     tempnegprofile = [-(x/neg_sum) for x in tempnegprofile]
            avgposprofile2 = [x+y for x,y in zip(avgposprofile2,tempposprofile)]
            avgnegprofile2 = [x+y for x,y in zip(avgnegprofile2, tempnegprofile)]
        # avgposprofile2 = [x/rep2number/mil_map for x in avgposprofile2]
        # avgnegprofile2 = [x/rep2number/mil_map for x in avgnegprofile2]
        avgposprofile2 = [x/rep2number for x in avgposprofile2]
        avgnegprofile2 = [x/rep2number for x in avgnegprofile2]
        posprofile2 = [x+y for x,y in zip(posprofile2,avgposprofile2)]
        negprofile2 = [x+y for x,y in zip(negprofile2, avgnegprofile2)]
    
    mil_map1 = mil_map1/rep1number
    mil_map2 = mil_map2/rep2number

    posprofile1 = [x/mil_map1 for x in posprofile1]
    negprofile1 = [x/mil_map1 for x in negprofile1]
    posprofile2 = [x/mil_map2 for x in posprofile2]
    negprofile2 = [x/mil_map2 for x in negprofile2]


    return posprofile1, negprofile1, posprofile2, negprofile2