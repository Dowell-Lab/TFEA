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
import subprocess
from pathlib import Path

import config
from exceptions import FileEmptyError, InputError, SubprocessError

#Main Script
#==============================================================================
def main(count_file=None, rank=config.RANK, scanner=config.SCANNER, 
            bam1=config.BAM1, bam2=config.BAM2, tempdir=Path(config.TEMPDIR), 
            label1=config.LABEL1, label2=config.LABEL2, 
            largewindow=config.LARGEWINDOW):
    '''This is the main script of the RANK module which takes as input a
        count file and bam files and ranks the regions within the count file
        according to a user specified 

    Parameters
    ----------
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

    Returns
    -------
    ranked_file : str
        Full path to a count file that contains the counts over regions within
        inputted bed files for all inputted bam files

    Raises
    ------
    FileEmptyError
        If any resulting file is empty
    '''
    if rank == 'deseq':
        if scanner == 'genome hits':
            ranked_file = deseq(bam1=bam1, bam2=bam2, tempdir=tempdir, 
                                    count_file=count_file, label1=label1, 
                                    label2=label2, largewindow=largewindow, 
                                    center=True)
        else:
            ranked_file = deseq(bam1=bam1, bam2=bam2, tempdir=tempdir, 
                                count_file=count_file, label1=label1, 
                                label2=label2, largewindow=largewindow, 
                                center=False)
    else:
        raise InputError("RANK option not recognized.")
    if os.stat(ranked_file).st_size == 0:
        raise FileEmptyError("Error in RANK module. Resulting file is empty.")

    return ranked_file

#Functions
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
res1 <- results(dds, alpha = 0.05, contrast=c("treatment", "'''+label2+'''",
                                                            "'''+label1+'''"))

resShrink <- lfcShrink(dds, res = res1, contrast = c("treatment", "'''+label2+'''",
                                                                "'''+label1+'''"))

resShrink$fc <- 2^(resShrink$log2FoldChange)
res <- resShrink[c(1:3,7,4:6)]
write.table(res, file = "''' + (tempdir / 'DESeq.res.txt').as_posix() + '''", append = FALSE, sep= "\t" )
sink()''')
    else:
        Rfile = open(tempdir /  'DESeq.R','w')
        Rfile.write('sink("' + tempdir /  'DESeq.Rout' + '")\n')
        Rfile.write('library("DESeq")\n')
        Rfile.write('data <- read.delim("'+count_file+'", sep="\t", \
                        header=TRUE)\n')

        Rfile.write('countsTable <- subset(data, select=c('
            +', '.join([str(i) for i in range(5,5+len(bam1)+len(bam2))])
            +'))\n')

        Rfile.write('rownames(countsTable) <- data$region\n')
        Rfile.write('conds <- c(' + ', '.join(['"'+label1+'"']*len(bam1)) 
                        + ', ' 
                        + ', '.join(['"'+label2+'"']*len(bam2)) 
                        + ')\n')

        Rfile.write('cds <- newCountDataSet( countsTable, conds )\n')
        Rfile.write('cds <- estimateSizeFactors( cds )\n')
        Rfile.write('sizeFactors(cds)\n')                                                               
        Rfile.write('cds <- estimateDispersions( cds ,method="blind", \
                        sharingMode="fit-only")\n')

        Rfile.write('res <- nbinomTest( cds, "'+label1+'", "'+label2+'" )\n')
        Rfile.write('rownames(res) <- res$id\n')                      
        Rfile.write('write.table(res, file = "'
                        + os.path.join(tempdir,'DESeq.res.txt') 
                        + '", append = FALSE, sep= "\t" )\n')

        # Rfile.write('sink()')
    Rfile.close()

#==============================================================================
def deseq(bam1=None, bam2=None, tempdir=None, count_file=None, label1=None, 
            label2=None, largewindow=None, center=False):
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
        raise SubprocessError(printmessage)

    if center:
        ranked_file = deseq_parse_center(deseq_file=deseq_file, tempdir=tempdir)
    else:
        ranked_file = deseq_parse(deseq_file=deseq_file, tempdir=tempdir, 
                                    largewindow=largewindow)

    return ranked_file

#==============================================================================
def deseq_parse_center(deseq_file=None, tempdir=None):
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
                fc = float(line[fc_index+1])
                if fc < 1:
                    down.append((chrom,start,stop,pval,str(fc)))
                else:
                    up.append((chrom,start,stop,pval,str(fc)))

    #Save ranked regions in a bed file (pvalue included)
    outfile = open(tempdir / "ranked_file.bed",'w')
    r=1
    for region in sorted(up, key=lambda x: x[3]):
        outfile.write('\t'.join(region) + '\t' + str(r) + '\n')
        r += 1
    for region in sorted(down, key=lambda x: x[3], reverse=True):
        outfile.write('\t'.join(region) + '\t' + str(r) + '\n')
        r += 1
    outfile.close()

    #Get center base for each region
    outfile = open(tempdir / "ranked_file.center.bed",'w')
    with open(tempdir / "ranked_file.bed") as F:
        for line in F:
            line = line.strip('\n').split('\t')
            chrom,start,stop = line[:3]
            center = int((int(start)+int(stop))/2)
            outfile.write(chrom + '\t' + str(center) + '\t' + str(center+1) 
                            + '\t' + '\t'.join(line[3:]) + '\n')
    outfile.close()

    with open(tempdir / "ranked_file.center.sorted.bed", 'w') as output:
        subprocess.run(["bedtools", "sort", 
                        "-i", tempdir / "ranked_file.center.bed"], 
                        stdout=output, check=True)

    ranked_center_file = tempdir / "ranked_file.center.sorted.bed"

    return ranked_center_file

#==============================================================================
def deseq_parse(deseq_file=None, tempdir=None, largewindow=None):
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
    with open(ranked_file,'w') as outfile:
        # outfile.write('\t'.join(['chrom', 'start', 'stop', 'rank, fc, p-value']) 
        #                 + '\n')
        r=1
        for region in sorted(up, key=lambda x: x[4]):
            outfile.write('\t'.join(region[:3]) 
                        + '\t' + ','.join(region[3:]+[str(r)]) 
                        + '\n')
            r += 1
        for region in sorted(down, key=lambda x: x[4], reverse=True):
            outfile.write('\t'.join(region[:3]) 
                        + '\t' + ','.join(region[3:]+[str(r)]) 
                        + '\n')
            r += 1

    return ranked_file