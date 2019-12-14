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
import warnings
import ujson
from pathlib import Path
from multiprocessing import Manager

from pybedtools import BedTool
import HTSeq as hts
import numpy as np

from TFEA import exceptions
from TFEA import multiprocess
from TFEA import plot

#Main Script
#==============================================================================
def main(use_config=True, combined_file=None, rank=None, scanner=None, 
            bam1=None, bam2=None, tempdir=None, label1=None, label2=None, 
            largewindow=None, mdd=False, mdd_bedfile1=False, mdd_bedfile2=False, 
            motif_annotations=False, debug=False, jobid=None, figuredir=None, 
            output_type=None, basemean_cut=None, plot_format=None):
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
        figuredir=config.vars['FIGUREDIR']
        combined_file=config.vars['COMBINED_FILE']
        rank=config.vars['RANK']
        scanner=config.vars['SCANNER']
        bam1=config.vars['BAM1']
        bam2=config.vars['BAM2']
        bg1=config.vars['BG1']
        bg2=config.vars['BG2']
        tempdir=config.vars['TEMPDIR']
        label1=config.vars['LABEL1']
        label2=config.vars['LABEL2']
        largewindow=config.vars['LARGEWINDOW']
        mdd=config.vars['MDD']
        mdd_bedfile1=config.vars['MDD_BEDFILE1']
        mdd_bedfile2=config.vars['MDD_BEDFILE2']
        mdd_pval=config.vars['MDD_PVAL']
        mdd_percent=config.vars['MDD_PERCENT']
        motif_annotations = config.vars['MOTIF_ANNOTATIONS']
        debug = config.vars['DEBUG']
        jobid = config.vars['JOBID']
        output_type = config.vars['OUTPUT_TYPE']
        basemean_cut = config.vars['BASEMEAN_CUT']
        plot_format = config.vars['PLOT_FORMAT']
        meta_profile_dict = {}
        metaprofile = config.vars['METAPROFILE']
    print("Ranking regions...", flush=True, file=sys.stderr)

    #Begin by counting reads from bam files over the combined_file produced
    # by the combine module
    if bam1 and bam2:
        count_file = count_reads(bedfile=combined_file, bam1=bam1, bam2=bam2, 
                            tempdir=tempdir, label1=label1, label2=label2)
        sample_number = len(bam1+bam2)
    elif bg1 and bg2:
        count_file = count_reads_bedtools(bedfile=combined_file, bg1=bg1, 
                                        bg2=bg2, tempdir=tempdir, label1=label1, 
                                        label2=label2)
        sample_number = len(bg1+bg2)
    millions_mapped = sum_reads(count_file=count_file, sample_number=sample_number)
    if motif_annotations:
        if bam1 and bam2:
            motif_fpkm = motif_count_reads(bedfile=motif_annotations, 
                                            bam1=bam1, bam2=bam2, 
                                            tempdir=tempdir, 
                                            label1=label1, label2=label2, 
                                            millions_mapped=millions_mapped)
        elif bg1 and bg2:
            motif_fpkm = motif_count_reads_bg(bedfile=motif_annotations, 
                                            bg1=bg1, bg2=bg2, 
                                            tempdir=tempdir, 
                                            label1=label1, label2=label2, 
                                            millions_mapped=millions_mapped)
        if use_config:
            config.vars['MOTIF_FPKM'] = motif_fpkm
        
    
    if os.stat(count_file).st_size == 0:
        raise exceptions.FileEmptyError("Error in RANK module. Counting failed.")

    if rank == 'deseq' or rank == 'fc':
        if bam1 and bam2:
            ranked_file, pvals, fcs = deseq(bam1=bam1, bam2=bam2, tempdir=tempdir, 
                                    count_file=count_file, label1=label1, 
                                    label2=label2, largewindow=largewindow, 
                                    rank=rank, figuredir=figuredir, 
                                    basemean_cut=basemean_cut, plot_format=plot_format)
        elif bg1 and bg2:
            ranked_file, pvals, fcs = deseq(bam1=bg1, bam2=bg2, tempdir=tempdir, 
                                    count_file=count_file, label1=label1, 
                                    label2=label2, largewindow=largewindow, 
                                    rank=rank, figuredir=figuredir, 
                                    basemean_cut=basemean_cut, plot_format=plot_format)
        if output_type == 'html' and metaprofile:
            print("\tGenerating Meta-Profile per Quartile:", file=sys.stderr)
            q1regions, q2regions, q3regions, q4regions = quartile_split(ranked_file)
            meta_profile_dict = meta_profile_quartiles(q1regions, q2regions, q3regions, q4regions, 
                                bam1=bam1, bam2=bam2, bg1=bg1, bg2=bg2, 
                                largewindow=largewindow, 
                                millions_mapped=millions_mapped, 
                                tempdir=tempdir)
        else:
            meta_profile_dict = False
    else:
        raise exceptions.InputError("RANK option not recognized.")
    if os.stat(ranked_file).st_size == 0:
        raise exceptions.FileEmptyError("Error in RANK module. DE-Seq running or parsing failed.")

    if use_config:
        config.vars['RANKED_FILE'] = ranked_file
        config.vars['PVALS'] = pvals
        config.vars['FCS'] = fcs
        config.vars['META_PROFILE'] = meta_profile_dict

    if mdd and (not mdd_bedfile1 or not mdd_bedfile2):
        mdd_bedfile1, mdd_bedfile2 = create_mdd_files(ranked_file=ranked_file,
                                                        tempdir=tempdir,
                                                        percent=mdd_percent,
                                                        pval_cut=mdd_pval)
        if os.stat(mdd_bedfile1).st_size == 0 or os.stat(mdd_bedfile2).st_size == 0:
            raise exceptions.FileEmptyError("Error in RANK module. MDD bed file creation failed.")
        if use_config:
            if not config.vars['MDD_BEDFILE1']:
                config.vars['MDD_BEDFILE1'] = mdd_bedfile1
            if not config.vars['MDD_BEDFILE2']:
                config.vars['MDD_BEDFILE2'] = mdd_bedfile2
    
    total_time = time.time() - start_time
    if use_config:
        config.vars['RANKtime'] = total_time
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
    count_file_header = tempdir / "count_file.header.bed"

    #pybedtools implementation (incomplete)
    # pybed_count = BedTool(bedfile.as_posix()).multi_bam_coverage(bams=bam1+bam2)
    # pybed_count.saveas(count_file, trackline=("#chrom\tstart\tstop\tregion\t" 
    #                                             + '\t'.join([label1]*len(bam1)) + "\t" 
    #                                             + '\t'.join([label2]*len(bam2)) + "\n"))
    
    for bamfile in bam1+bam2:
        samtools_index_command = ["samtools", "index", bamfile]
        subprocess.run(samtools_index_command)

    #Bedtools implementation
    multicov_command = ["bedtools", "multicov", 
                        "-bams"] + bam1+bam2 + ["-bed", bedfile]
    with open(count_file, 'w') as outfile:
        subprocess.run(multicov_command, stdout=outfile, check=True)

    # This section adds a header to the count_file and reformats it to remove 
    # excess information and add a column with the region for later use
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
def count_reads_bedtools(bedfile=None, bg1=None, bg2=None, tempdir=None, label1=None, 
                            label2=None):
    '''Counts reads across regions in a given bed file using bam files inputted
        by a user

    Parameters
    ----------
    bedfile : string
        full path to a bed file containing full regions of interest which will 
        be counted using bedtools multicov

    bg1 : list or array
        a list of full paths to bam files pertaining to a single condition 
        (i.e. replicates of a single treatment)

    bg2 : list or array
        a list of full paths to bam files pertaining to a single condition 
        (i.e. replicates of a single treatment)

    tempdir : string
        full path to temp directory in output directory (created by TFEA)

    label1 : string
        the name of the treatment or condition corresponding to bg1 list

    label2 : string
        the name of the treatment or condition corresponding to bg2 list

    Returns
    -------
    None
    '''
    #This os.system call runs bedtools multicov to count reads in all specified
    #BAMs for given regions in BED
    count_file = tempdir / "count_file.bed"
    count_file_header = tempdir / "count_file.header.bed"

    #pybedtools implementation (incomplete)
    # pybed_count = BedTool(bedfile.as_posix()).multi_bam_coverage(bams=bg1+bg2)
    # pybed_count.saveas(count_file, trackline=("#chrom\tstart\tstop\tregion\t" 
    #                                             + '\t'.join([label1]*len(bg1)) + "\t" 
    #                                             + '\t'.join([label2]*len(bg2)) + "\n"))

    #Bedtools implementation
    multicov_command = ["bedtools", "map", 
                        "-a", bedfile, ] + bg1+bg2 + ["-bed", bedfile]
    multicov_command = ""
    bedfile = str(bedfile)
    for i,bgfile in enumerate(bg1+bg2):
        bgfile = f"<(awk -F'\\t' -v OFS='\\t' 'function abs(x) {{return ((x < 0.0) ? -x : x)}} {{print $1,$2,$3,abs($4)}}' {bgfile})"
        if i == 0:
            multicov_command = "bedtools map -null 0 -o sum -c 4 -a "+bedfile+" -b "+bgfile+" | "
        else:
            multicov_command = multicov_command + "bedtools map -null 0 -o sum -c 4 -a stdin -b "+bgfile
    with open(count_file, 'w') as outfile:
        try:
            subprocess.run(['/bin/bash','-c', multicov_command],
                            stdout=outfile, stderr=subprocess.PIPE, check=True)
        except subprocess.CalledProcessError as e:
            raise exceptions.SubprocessError(e.stderr.decode())

    # This section adds a header to the count_file and reformats it to remove 
    # excess information and add a column with the region for later use
    outfile = open(count_file_header, 'w')
    outfile.write("#chrom\tstart\tstop\tregion\t" 
                    + '\t'.join([label1]*len(bg1)) + "\t" 
                    + '\t'.join([label2]*len(bg2)) + "\n")

    with open(count_file) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            chrom,start,stop = line[:3]
            counts = line[-(len(bg1)+len(bg2)):]
            outfile.write('\t'.join([chrom,start,stop]) + "\t" 
                            + chrom + ":" + start + "-" + stop + "\t"
                            + '\t'.join(counts) + "\n")
    outfile.close()
    return count_file_header

#==============================================================================
def motif_count_reads(bedfile=None, bam1=None, bam2=None, tempdir=None, 
                            label1=None, label2=None, millions_mapped=None):
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
    count_file = tempdir / "motif_counts.bed"
    count_file_header = tempdir / "motif_counts.header.fpkm.bed"

    #pybedtools implementation (incomplete)
    # pybed_count = BedTool(bedfile.as_posix()).multi_bam_coverage(bams=bam1+bam2)
    # pybed_count.saveas(count_file, trackline=("#chrom\tstart\tstop\tregion\t" 
    #                                             + '\t'.join([label1]*len(bam1)) + "\t" 
    #                                             + '\t'.join([label2]*len(bam2)) + "\n"))

    #Bedtools implementation
    multicov_command = ["bedtools", "multicov", "-s", 
                        "-bams"] + bam1+bam2 + ["-bed", bedfile]
    with open(count_file, 'w') as outfile:
        subprocess.run(multicov_command, stdout=outfile, check=True)

    # This section adds a header to the count_file and reformats it to remove 
    # excess information and add a column with the region for later use
    outfile = open(count_file_header, 'w')
    outfile.write("#chrom\tstart\tstop\tregion\t" 
                    + '\t'.join([label1]*len(bam1)) + "\t" 
                    + '\t'.join([label2]*len(bam2)) + "\n")
    motif_fpkm = {}
    with open(count_file) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            chrom,start,stop,motif = line[:4]
            length = float(stop) - float(start)
            fpkm = [float(c)/(m/1000000.0)/(length/1000.0) for c,m in zip(line[-(len(bam1)+len(bam2)):], millions_mapped)]
            motif_fpkm[motif] = np.mean(fpkm)
            outfile.write('\t'.join([chrom,start,stop]) + "\t" 
                            + motif + "\t"
                            + '\t'.join([str(f) for f in fpkm]) + "\n")
    outfile.close()
    return motif_fpkm

#==============================================================================
def motif_count_reads_bg(bedfile=None, bg1=None, bg2=None, tempdir=None, 
                            label1=None, label2=None, millions_mapped=None):
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
    count_file = tempdir / "motif_counts.bed"
    count_file_header = tempdir / "motif_counts.header.fpkm.bed"

    #pybedtools implementation (incomplete)
    # pybed_count = BedTool(bedfile.as_posix()).multi_bam_coverage(bams=bam1+bam2)
    # pybed_count.saveas(count_file, trackline=("#chrom\tstart\tstop\tregion\t" 
    #                                             + '\t'.join([label1]*len(bam1)) + "\t" 
    #                                             + '\t'.join([label2]*len(bam2)) + "\n"))

    #Bedtools implementation
    multicov_command = ["bedtools", "map", 
                        "-a", bedfile, ] + bg1+bg2 + ["-bed", bedfile]
    multicov_command = ""
    bedfile = str(bedfile)
    for i,bgfile in enumerate(bg1+bg2):
        bgfile = f"<(awk -F'\\t' -v OFS='\\t' 'function abs(x) {{return ((x < 0.0) ? -x : x)}} {{print $1,$2,$3,abs($4)}}' {bgfile})"
        if i == 0:
            multicov_command = "bedtools map -null 0 -c 4 -a "+bedfile+" -b "+bgfile+" | "
        else:
            multicov_command = multicov_command + "bedtools map -null 0 -c 4 -a stdin -b "+bgfile
    with open(count_file, 'w') as outfile:
        try:
            subprocess.run(['/bin/bash','-c', multicov_command],
                            stdout=outfile, stderr=subprocess.PIPE, check=True)
        except subprocess.CalledProcessError as e:
            raise exceptions.SubprocessError(e.stderr.decode())

    # This section adds a header to the count_file and reformats it to remove 
    # excess information and add a column with the region for later use
    outfile = open(count_file_header, 'w')
    outfile.write("#chrom\tstart\tstop\tregion\t" 
                    + '\t'.join([label1]*len(bg1)) + "\t" 
                    + '\t'.join([label2]*len(bg2)) + "\n")
    motif_fpkm = {}
    with open(count_file) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            chrom,start,stop,motif = line[:4]
            length = float(stop) - float(start)
            fpkm = [float(c)/(m/1000000.0)/(length/1000.0) for c,m in zip(line[-(len(bg1)+len(bg2)):], millions_mapped)]
            motif_fpkm[motif] = np.mean(fpkm)
            outfile.write('\t'.join([chrom,start,stop]) + "\t" 
                            + motif + "\t"
                            + '\t'.join([str(f) for f in fpkm]) + "\n")
    outfile.close()
    return motif_fpkm

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
            millions_mapped = [x+float(y) for x,y in zip(millions_mapped,line[-1*sample_number:])]

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
print("Size Factors")
print(sizeFactors(dds))
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
print("Size Factors")
print(sizeFactors(cds))
cds <- estimateDispersions( cds ,method="blind", \
                        sharingMode="fit-only")

res <- nbinomTest( cds, "'''+label1+'''", "'''+label2+'''" )
rownames(res) <- res$id                      
write.table(res, file = "'''    + os.path.join(tempdir,'DESeq.res.txt') 
                                + '''", append = FALSE, sep= "\t" )''')
    Rfile.close()

#==============================================================================
def deseq(bam1=None, bam2=None, tempdir=None, count_file=None, label1=None, 
            label2=None, largewindow=None, rank=None, figuredir=None, 
            basemean_cut=None, plot_format=None):
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

    
    plot.plot_deseq_MA(deseq_file=deseq_file, label1=label1, label2=label2, 
                        figuredir=figuredir, basemean_cut=basemean_cut, 
                        plot_format=plot_format)

    ranked_file, pvals, fcs = deseq_parse(deseq_file=deseq_file, tempdir=tempdir, 
                                largewindow=largewindow, rank=rank, 
                                basemean_cut=basemean_cut)

    return ranked_file, pvals, fcs

#==============================================================================
def deseq_parse(deseq_file=None, tempdir=None, largewindow=None, rank=None, 
                basemean_cut=0):
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
        basemean_index = [i for i in range(len(header)) if header[i]=='"baseMean"'][0]
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
                start = start if start > 0 else 0
                stop = center + int(largewindow)
                fc = float(line[fc_index+1])
                basemean = float(line[basemean_index+1])
                if basemean > basemean_cut:
                    if fc < 1:
                        down.append([chrom, str(start), str(stop), str(fc), str(pval)])
                    else:
                        up.append([chrom, str(start), str(stop), str(fc), str(pval)])

    #Save ranked regions in a bed file (pvalue included)
    ranked_file = tempdir / "ranked_file.bed"
    pvals = list()
    fcs = list()
    with open(ranked_file,'w') as outfile:
        outfile.write('\t'.join(['#chrom', 'start', 'stop', 'fc,p-value,rank']) 
                        + '\n')
        r=1
        if rank == 'deseq':
            for region in sorted(up, key=lambda x: x[4]):
                outfile.write('\t'.join(region[:3]) 
                            + '\t' + ','.join(region[3:]+[str(r)]) 
                            + '\n')
                pvals.append(float(region[-1]))
                fcs.append(float(region[3]))
                r += 1
            for region in sorted(down, key=lambda x: x[4], reverse=True):
                outfile.write('\t'.join(region[:3]) 
                            + '\t' + ','.join(region[3:]+[str(r)]) 
                            + '\n')
                pvals.append(float(region[-1]))
                fcs.append(float(region[3]))
                r += 1
        elif rank == 'fc':
            for region in sorted(up+down, key=lambda x: x[3], reverse=True):
                outfile.write('\t'.join(region[:3]) 
                            + '\t' + ','.join(region[3:]+[str(r)]) 
                            + '\n')
                pvals.append(float(region[-1]))
                fcs.append(float(region[3]))
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
    
    q1 = int(round(np.percentile(np.arange(1, len(regions),1), 25)))
    q2 = int(round(np.percentile(np.arange(1, len(regions),1), 50)))
    q3 = int(round(np.percentile(np.arange(1, len(regions),1), 75)))

    q1regions = regions[:q1]
    q2regions = regions[q1:q2]
    q3regions = regions[q2:q3]
    q4regions = regions[q3:]

    return q1regions, q2regions, q3regions, q4regions
#==============================================================================
def meta_profile_quartiles(q1regions, q2regions, q3regions, q4regions, 
                            bam1=None, bam2=None, bg1=None, bg2=None, 
                            largewindow=None,
                            tempdir=None, millions_mapped=None):
    '''This function creates a metaprofile from 4 regions and stores them in
        a dictionary. 
    '''
    if len(q1regions + q2regions + q3regions + q4regions)*largewindow < 6e7 and bam1 and bam2:
        q1regions = ['q1'] + q1regions
        q2regions = ['q2'] + q2regions
        q3regions = ['q3'] + q3regions
        q4regions = ['q4'] + q4regions
        kwargs = dict(largewindow=largewindow, bam1=bam1, bam2=bam2)
        meta_profile_tuples = multiprocess.main(function=meta_profile, 
                        args=[q1regions, q2regions, q3regions, q4regions], 
                        kwargs=kwargs)
    elif bam1 and bam2:
        print("\tToo many regions to parallelize over. Performing computation in serial.", 
                file=sys.stderr)
        q1regions = ['q1'] + q1regions
        q2regions = ['q2'] + q2regions
        q3regions = ['q3'] + q3regions
        q4regions = ['q4'] + q4regions
        meta_profile_tuples = np.empty(4,dtype=list)
        meta_profile_tuples[0] = meta_profile(regionlist=q1regions, 
                                        largewindow=largewindow, bam1=bam1, 
                                        bam2=bam2)

        meta_profile_tuples[1] = meta_profile(regionlist=q2regions, 
                                        largewindow=largewindow, bam1=bam1, 
                                        bam2=bam2)

        meta_profile_tuples[2] = meta_profile(regionlist=q3regions, 
                                        largewindow=largewindow, bam1=bam1, 
                                        bam2=bam2)
    
        meta_profile_tuples[3] = meta_profile(regionlist=q4regions, 
                                        largewindow=largewindow, bam1=bam1, 
                                        bam2=bam2)
    elif bg1 and bg2:
        regionlist = q1regions+q2regions+q3regions+q4regions
        meta_profile_tuples = meta_profile_bg(regionlist=regionlist, 
                                                largewindow=largewindow, 
                                                bg1=bg1, bg2=bg2)

    meta_profile_dict = {}
    if bam1 and bam2:
        mil_map1 = sum(millions_mapped[:len(bam1)])/len(bam1)/1e6
        mil_map2 = sum(millions_mapped[-len(bam2):])/len(bam2)/1e6
    elif bg1 and bg2:
        mil_map1 = sum(millions_mapped[:len(bg1)])/len(bg1)/1e6
        mil_map2 = sum(millions_mapped[-len(bg2):])/len(bg2)/1e6
    for profile_list in meta_profile_tuples:
        for key, profile in profile_list:
            if key[-1] == '1':
                profile = [[y/mil_map1 for y in x] for x in profile]
            elif key[-1] == '2':
                profile = [[y/mil_map2 for y in x] for x in profile]
            meta_profile_dict[key] = profile
    if tempdir is None:
        return meta_profile_dict
    else:
        meta_profile_folder = tempdir / 'meta_profile'
        meta_profile_folder.mkdir(exist_ok=True)
        
        for key in meta_profile_dict:
            for i, profile in enumerate(meta_profile_dict[key]):
                profile_file = meta_profile_folder / (f'{key}_{i}')
                profile_file.write_text(ujson.dumps(profile))
                
        return meta_profile_folder

        #JSON single file
        # meta_profile_file = tempdir / 'meta_profile.json'
        # meta_profile_file.write_text(ujson.dumps(meta_profile_dict))

        # return meta_profile_file

        #Pickling
        # with open(meta_profile_file, 'wb') as f:
        #     pickle.dump(meta_profile_dict, f)

        # #Pure Python
        # with open(meta_profile_file, 'w') as outfile:
        #     for key in meta_profile_dict:
        #         outfile.write('#' + key + '\n')
        #         outfile.write('\n'.join([','.join([str(y) for y in x]) for x in profile]) + '\n')
        

#==============================================================================
def meta_profile(regionlist=None, largewindow=None, bam1=None, bam2=None):
    '''This function returns average profiles for given regions of interest.
        A user may input either a list of regions or a bed file
    Parameters
    ----------
    regionlist : list
        a list of regions of interest. Format: [(chrom, start, stop), (), ...]
        First value is special 'key_prefix' which simply keeps track of which
        quartile we're in.
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
    key_prefix = regionlist[0]
    regions=list()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for chrom, start, stop in regionlist[1:]:
            regions.append(hts.GenomicInterval(chrom, int(start), int(stop), '.'))

    hts_bam1 = list()
    hts_bam2 = list()
    for bam in bam1:
        hts_bam1.append(hts.BAM_Reader(bam))
    for bam in bam2:
        hts_bam2.append(hts.BAM_Reader(bam))

    if len(hts_bam1) == 0 or len(hts_bam2) == 0:
        raise ValueError("One of bam1 or bam2 variables is empty.")

    # posprofile1 = np.zeros(2*int(largewindow))  
    # negprofile1 = np.zeros(2*int(largewindow))
    # posprofile2 = np.zeros(2*int(largewindow))   
    # negprofile2 = np.zeros(2*int(largewindow))
    # posprofile1 = []
    # negprofile1 = []
    # posprofile2 = []   
    # negprofile2 = []
    posprofile1 = np.empty((len(regions), 2*int(largewindow)))
    negprofile1 = np.empty((len(regions), 2*int(largewindow)))
    posprofile2 = np.empty((len(regions), 2*int(largewindow)))
    negprofile2 = np.empty((len(regions), 2*int(largewindow)))
    rep1number = float(len(hts_bam1))
    rep2number = float(len(hts_bam2))
    for i, window in enumerate(regions):
        avgposprofile1 = np.zeros(2*int(largewindow))
        avgnegprofile1 = np.zeros(2*int(largewindow))
        # i = 0
        for sortedbamfile in hts_bam1:
            # i += 1
            tempposprofile = np.zeros(2*int(largewindow))
            tempnegprofile = np.zeros(2*int(largewindow))
            for almnt in sortedbamfile[ window ]:
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
            # avgposprofile1 = [x+y for x,y in zip(avgposprofile1, tempposprofile)]
            # avgnegprofile1 = [x+y for x,y in zip(avgnegprofile1, tempnegprofile)]

            avgposprofile1 = np.add(avgposprofile1, tempposprofile)
            avgnegprofile1 = np.add(avgnegprofile1, tempnegprofile)
        # avgposprofile1 = [x/rep1number for x in avgposprofile1]
        # avgnegprofile1 = [x/rep1number for x in avgnegprofile1]
        avgposprofile1 = avgposprofile1/rep1number
        avgnegprofile1 = avgnegprofile1/rep1number
        # posprofile1 = [x+y for x,y in zip(posprofile1, avgposprofile1)]
        # negprofile1 = [x+y for x,y in zip(negprofile1, avgnegprofile1)]
        # posprofile1.append(avgposprofile1)
        # negprofile1.append(avgnegprofile1)
        posprofile1[i] = avgposprofile1
        negprofile1[i] = avgnegprofile1

        avgposprofile2 = np.zeros(2*int(largewindow))
        avgnegprofile2 = np.zeros(2*int(largewindow))
        # i = len(hts_bam1)
        for sortedbamfile in hts_bam2:
            # i += 1
            tempposprofile = np.zeros(2*int(largewindow))
            tempnegprofile = np.zeros(2*int(largewindow))
            for almnt in sortedbamfile[ window ]:
                if almnt.iv.strand == '+':
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
            # avgposprofile2 = [x+y for x,y in zip(avgposprofile2, tempposprofile)]
            # avgnegprofile2 = [x+y for x,y in zip(avgnegprofile2, tempnegprofile)]
            avgposprofile2 = np.add(avgposprofile2, tempposprofile)
            avgnegprofile2 = np.add(avgnegprofile2, tempnegprofile)
        # avgposprofile2 = [x/rep2number for x in avgposprofile2]
        # avgnegprofile2 = [x/rep2number for x in avgnegprofile2]
        avgposprofile2 = avgposprofile2/rep2number
        avgnegprofile2 = avgnegprofile2/rep2number
        # posprofile2 = [x+y for x,y in zip(posprofile2, avgposprofile2)]
        # negprofile2 = [x+y for x,y in zip(negprofile2, avgnegprofile2)]
        # posprofile2.append(avgposprofile2)
        # negprofile2.append(avgnegprofile2)
        posprofile2[i] = avgposprofile2
        negprofile2[i] = avgnegprofile2
    
    # mil_map1 = mil_map1/rep1number/1e6
    # mil_map2 = mil_map2/rep2number/1e6

    # posprofile1 = [[x/mil_map1 for x in posprofile1]]
    # negprofile1 = [[x/mil_map1 for x in negprofile1]]
    # posprofile2 = [[x/mil_map2 for x in posprofile2]]
    # negprofile2 = [[x/mil_map2 for x in negprofile2]]


    return (key_prefix + 'posprofile1', posprofile1), (key_prefix + 'negprofile1', negprofile1), (key_prefix + 'posprofile2', posprofile2), (key_prefix + 'negprofile2', negprofile2)

#==============================================================================
def meta_profile_bg(regionlist=None, largewindow=None, bg1=None, bg2=None):
    '''This function returns average profiles for given regions of interest.
        A user may input either a list of regions or a bed file
    Parameters
    ----------
    regionlist : list
        a list of regions of interest. Format: [(chrom, start, stop), (), ...].
        First value is special 'key_prefix' which simply keeps track of which
        quartile we're in.
    largewindow : float
        the window with which to compute profiles for
    bg1 : list
        a list of full paths to bedgraph files corresponding to a condition or 
        treatment
    bg2 : list
        a list of full paths to bedgraph files corresponding to a condition or 
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
    from ncls import NCLS
    number_of_regions = len(regionlist)
    total_window = 2*int(largewindow)
    
    posprofile1 = np.empty((number_of_regions, total_window))
    negprofile1 = np.empty((number_of_regions, total_window))
    posprofile2 = np.empty((number_of_regions, total_window))
    negprofile2 = np.empty((number_of_regions, total_window))

    #Create a dictionary formatted {'chrom': rank: (start, stop, rank)}
    region_dict = {}
    for i, (chrom, start, stop) in enumerate(regionlist):
        start = int(start)
        stop = int(stop)
        if chrom not in region_dict:
            region_dict[chrom] = {}
        region_dict[chrom][i] = (start, stop, i)

    current_chrom = None
    id_pointer = {}
    for bgfile in bg1:
        with open(bgfile, 'r') as F:
            for i, line in enumerate(F):
                if '#' not in line[0]:
                    chrom, start, stop, value = line.strip().split('\t')
                    start = int(start)
                    stop = int(stop)
                    value = int(value)
                    if current_chrom == None: #If starting, initialize
                        current_chrom = chrom
                        id_pointer[i] = (chrom, start, stop, value)
                        bedgraph_regions = [(start, stop, i)]
                    elif chrom != current_chrom: #If new chrom transition:
                        #1. Perform intersect computation
                        reg_starts, reg_stops, reg_ids = zip(*[region_dict[chrom][i] for i in region_dict[chrom]])
                        bg_starts, bg_stops, bg_ids = zip(*bedgraph_regions)
                        region_ncls = NCLS(np.array(reg_starts), 
                                            np.array(reg_stops), 
                                            np.array(reg_ids))
                        all_overlap_ids = region_ncls.all_overlaps_both(np.array(bg_starts), 
                                                            np.array(bg_stops), 
                                                            np.array(bg_ids))
                        #2. Add profile 
                        for i, j in np.column_stack((all_overlap_ids[0], all_overlap_ids[1])):
                            chrom, bg_start, bg_stop, bg_value = id_pointer[i]
                            reg_start, reg_stop, _ = region_dict[chrom][j]
                            j_start = bg_start - reg_start if bg_start - reg_start > 0 else 0
                            j_stop = bg_stop - reg_start if bg_stop - reg_start < total_window else total_window
                            if bg_value > 0:
                                posprofile1[j][j_start:j_stop] += bg_value
                            else:
                                negprofile1[j][j_start:j_stop] += bg_value
                        
                        #3. Initialize new round
                        current_chrom = chrom
                        id_pointer = {}
                        id_pointer[i] = (chrom, start, stop, value)
                        bedgraph_regions = [(start, stop, i)]
                    else: #If in same chormosome, continue appending to list
                        bedgraph_regions.append((start,stop,i))
                        id_pointer[i] = (chrom, start, stop, value)
        #Reached end of file. Perform computation and profiling 
        #1. Perform intersect computation and profiling
        reg_starts, reg_stops, reg_ids = zip(*[region_dict[chrom][i] for i in region_dict[chrom]])
        bg_starts, bg_stops, bg_ids = zip(*bedgraph_regions)
        region_ncls = NCLS(np.array(reg_starts), 
                            np.array(reg_stops), 
                            np.array(reg_ids))
        all_overlap_ids = region_ncls.all_overlaps_both(np.array(bg_starts), 
                                                        np.array(bg_stops), 
                                                        np.array(bg_ids))
        
        for i, j in np.column_stack((all_overlap_ids[0], all_overlap_ids[1])):
            chrom, bg_start, bg_stop, bg_value = id_pointer[i]
            reg_start, reg_stop, _ = region_dict[chrom][j]
            j_start = bg_start - reg_start if bg_start - reg_start > 0 else 0
            j_stop = bg_stop - reg_start if bg_stop - reg_start < total_window else total_window
            if bg_value > 0:
                posprofile1[j][j_start:j_stop] += bg_value
            else:
                negprofile1[j][j_start:j_stop] += bg_value

    for bgfile in bg2:
        with open(bgfile, 'r') as F:
            for i, line in enumerate(F):
                if '#' not in line[0]:
                    chrom, start, stop, value = line.strip().split('\t')
                    start = int(start)
                    stop = int(stop)
                    value = int(value)
                    if current_chrom == None: #If starting, initialize
                        current_chrom = chrom
                        id_pointer[i] = (chrom, start, stop, value)
                        bedgraph_regions = [(start, stop, i)]
                    elif chrom != current_chrom: #If new chrom transition:
                        #1. Perform intersect computation
                        reg_starts, reg_stops, reg_ids = zip(*[region_dict[chrom][i] for i in region_dict[chrom]])
                        bg_starts, bg_stops, bg_ids = zip(*bedgraph_regions)
                        region_ncls = NCLS(np.array(reg_starts), 
                                            np.array(reg_stops), 
                                            np.array(reg_ids))
                        all_overlap_ids = region_ncls.all_overlaps_both(np.array(bg_starts), 
                                                            np.array(bg_stops), 
                                                            np.array(bg_ids))
                        #2. Add profile 
                        for i, j in np.column_stack((all_overlap_ids[0], all_overlap_ids[1])):
                            chrom, bg_start, bg_stop, bg_value = id_pointer[i]
                            reg_start, reg_stop, _ = region_dict[chrom][j]
                            j_start = bg_start - reg_start if bg_start - reg_start > 0 else 0
                            j_stop = bg_stop - reg_start if bg_stop - reg_start < total_window else total_window
                            if bg_value > 0:
                                posprofile2[j][j_start:j_stop] += bg_value
                            else:
                                negprofile2[j][j_start:j_stop] += bg_value
                        
                        #3. Initialize new round
                        current_chrom = chrom
                        id_pointer = {}
                        id_pointer[i] = (chrom, start, stop, value)
                        bedgraph_regions = [(start, stop, i)]
                    else: #If in same chormosome, continue appending to list
                        bedgraph_regions.append((start,stop,i))
                        id_pointer[i] = (chrom, start, stop, value)
        #Reached end of file. Perform computation and profiling 
        #1. Perform intersect computation and profiling
        reg_starts, reg_stops, reg_ids = zip(*[region_dict[chrom][i] for i in region_dict[chrom]])
        bg_starts, bg_stops, bg_ids = zip(*bedgraph_regions)
        region_ncls = NCLS(np.array(reg_starts), 
                            np.array(reg_stops), 
                            np.array(reg_ids))
        all_overlap_ids = region_ncls.all_overlaps_both(np.array(bg_starts), 
                                                        np.array(bg_stops), 
                                                        np.array(bg_ids))
        
        for i, j in np.column_stack((all_overlap_ids[0], all_overlap_ids[1])):
            chrom, bg_start, bg_stop, bg_value = id_pointer[i]
            reg_start, reg_stop, _ = region_dict[chrom][j]
            j_start = bg_start - reg_start if bg_start - reg_start > 0 else 0
            j_stop = bg_stop - reg_start if bg_stop - reg_start < total_window else total_window
            if bg_value > 0:
                posprofile2[j][j_start:j_stop] += bg_value
            else:
                negprofile2[j][j_start:j_stop] += bg_value
    
    q1 = int(round(np.percentile(np.arange(1, len(regionlist),1), 25)))
    q2 = int(round(np.percentile(np.arange(1, len(regionlist),1), 50)))
    q3 = int(round(np.percentile(np.arange(1, len(regionlist),1), 75)))

    q1posprofile1 = posprofile1[:q1]
    q1negprofile1 = negprofile1[:q1]
    q1posprofile2 = posprofile2[:q1]
    q1negprofile2 = negprofile2[:q1]

    q2posprofile1 = posprofile1[q1:q2]
    q2negprofile1 = negprofile1[q1:q2]
    q2posprofile2 = posprofile2[q1:q2]
    q2negprofile2 = negprofile2[q1:q2]

    q3posprofile1 = posprofile1[q2:q3]
    q3negprofile1 = negprofile1[q2:q3]
    q3posprofile2 = posprofile2[q2:q3]
    q3negprofile2 = negprofile2[q2:q3]

    q4posprofile1 = posprofile1[q3:]
    q4negprofile1 = negprofile1[q3:]
    q4posprofile2 = posprofile2[q3:]
    q4negprofile2 = negprofile2[q3:]

    return np.array([[('q1' + 'posprofile1', q1posprofile1), 
                    ('q1' + 'negprofile1', q1negprofile1), 
                    ('q1' + 'posprofile2', q1posprofile2), 
                    ('q1' + 'negprofile2', q1negprofile2)],
                    [('q2' + 'posprofile1', q2posprofile1), 
                    ('q2' + 'negprofile1', q2negprofile1), 
                    ('q2' + 'posprofile2', q2posprofile2), 
                    ('q2' + 'negprofile2', q2negprofile2)],
                    [('q3' + 'posprofile1', q3posprofile1), 
                    ('q3' + 'negprofile1', q3negprofile1), 
                    ('q3' + 'posprofile2', q3posprofile2), 
                    ('q3' + 'negprofile2', q3negprofile2)],
                    [('q4' + 'posprofile1', q4posprofile1), 
                    ('q4' + 'negprofile1', q4negprofile1), 
                    ('q4' + 'posprofile2', q4posprofile2), 
                    ('q4' + 'negprofile2', q4negprofile2)]])
