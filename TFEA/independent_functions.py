#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'Jonathan D. Rubin and Rutendo Sigauke'
__credits__ = ['Jonathan D. Rubin', 'Rutendo Sigauke', 'Jacob Stanley',
                'Robin Dowell']
__maintainer__ = 'Jonathan D. Rubin'
__email__ = 'Jonathan.Rubin@colorado.edu'

'''This file contains a list of independent functions that do not call other functions
'''

import matplotlib
matplotlib.use('Agg')
import os
import math
import matplotlib.pyplot as plt
import numpy as np

import config

#===================================================================================
def combine_bed(beds=config.BEDS,tempdir=config.TEMPDIR):
    '''Concatenates, sorts, and merges (bedtools) a list of bed files. Outputs into the
        tempdir directory created by TFEA

    Parameters
    ----------
    beds : list or array
        full paths to bed files (strings)
        
    tempdir : string
        full path to tempdir directory in output directory (created by TFEA)

    Returns
    -------
    combined_input_merged_bed : string 
        full path to a bed file containing the merged regions inputted by the user 
    '''
    os.system("cat " + " ".join(beds) + " > " + tempdir + "combined_input.bed")

    os.system("sort -k1,1 -k2,2n " + tempdir + "combined_input.bed > " + tempdir \
                + "combined_input.sorted.bed")

    os.system("bedtools merge -i " + tempdir + "combined_input.sorted.bed > " \
                + tempdir + "combined_input.merge.bed")

    combined_input_merged_bed = tempdir + "combined_input.merge.bed"

    return combined_input_merged_bed
#===================================================================================

#===================================================================================
def getfasta(bedfile=bedfile, genomefasta=config.GENOMEFASTA, tempdir=config.TEMPDIR):
    '''Converts a bed file to a fasta file using bedtools. Outputs into the tempdir 
        directory created by TFEA.

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
        full path to a fasta file containing the inputted bed file regions in fasta
        format 
    '''
    os.system("bedtools getfasta -name -fi "+genomefasta+" -bed "+bedfile+" -fo " \
                + tempdir + "ranked_file.fullregions.fa")

    ranked_file_fasta = tempdir + "ranked_file.fullregions.fa"

    return ranked_file_fasta
#===================================================================================

#===================================================================================
def get_bgfile(fastafile=fastafile, tempdir=config.TEMPDIR):
    '''Obtains a zero order markov background model (used in FIMO) from a fasta file. 
        Outputs into the tempdir directory created by TFEA.

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
        full path to a meme-formatted zero order markov background file (txt file)
    '''
    os.system("fasta-get-markov "+fastafile+" "+tempdir+"markov_background.txt")

    markov_background = tempdir + "markov_background.txt"

    return markov_background
#===================================================================================

#===================================================================================
def get_regions(tempdir=config.TEMPDIR, ranked_center_file=config.RANKED_CENTER_FILE,
                largewindow=config.LARGEWINDOW):
    '''Takes in a bed file that contains regions centered on the user-inputted bed files
        and outputs a 'full regions' bed file which simply adds a user-defined window
        to either side of these centered regions.

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
#===================================================================================

#===================================================================================
def count_reads(bedfile=bedfile, bam1=config.BAM1, bam2=config.BAM2, 
                tempdir=config.TEMPDIR, label1=config.LABEL1, label2=config.LABEL2):
    '''Counts reads across regions in a given bed file using bam files inputted by a
        user

    Parameters
    ----------
    bedfile : string
        full path to a bed file containing full regions of interest which will be
        counted using bedtools multicov

    bam1 : list or array
        a list of full paths to bam files pertaining to a single condition (i.e.
        replicates of a single treatment)

    bam2 : list or array
        a list of full paths to bam files pertaining to a single condition (i.e.
        replicates of a single treatment)

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
    #This os.system call runs bedtools multicov to count reads in all specified BAMs
    #for given regions in BED
    os.system("bedtools multicov -bams " + " ".join(bam1) + " " + " ".join(bam2) \
                + " -bed " + BED + " > " + tempdir + "count_file.bed")

    #This section adds a header to the count_file and reformats it to remove excess 
    #information and add a column with the region for later use
    outfile = open(tempdir + "count_file.header.bed",'w')
    outfile.write("chrom\tstart\tstop\tregion\t" + '\t'.join([label1]*len(bam1)) \
                    + "\t" + '\t'.join([label2]*len(bam2)) + "\n")
    with open(tempdir + "count_file.bed") as F:
        for line in F:
            line = line.strip('\n').split('\t')
            chrom,start,stop = line[:3]
            counts = line[-(len(bam1)+len(bam2)):]
            outfile.write('\t'.join([chrom,start,stop]) + "\t"+chrom+":"+start+"-"\
                            +stop+"\t" + '\t'.join(counts) + "\n")
#===================================================================================

#===================================================================================
def write_deseq_script(bam1=config.BAM1, bam2=config.BAM2, tempdir=config.TEMPDIR,
                        count_file=config.COUNT_FILE, label1=config.LABEL1,
                        label2=config.LABEL2):
    '''Writes an R script within the tempdir directory in TFEA output to run either 
    DE-Seq or DE-Seq2 depending on the number of user-inputted replicates.

    Parameters
    ----------
    bedfile : string
        full path to a bed file containing full regions of interest which will be
        counted using bedtools multicov

    bam1 : list or array
        a list of full paths to bam files pertaining to a single condition (i.e.
        replicates of a single treatment)

    bam2 : list or array
        a list of full paths to bam files pertaining to a single condition (i.e.
        replicates of a single treatment)

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
        outfile = open(tempdir + 'DESeq.R','w')
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
write.table(res, file = "'''+tempdir+'''DESeq.res.txt", append = FALSE, sep= "\t" )
sink()''')
#===================================================================================

#===================================================================================
def plot_deseq_MA(deseq_file=deseq_file, label1=config.LABEL1, label2=config.LABEL2,
                    figuredir=config.FIGUREDIR):
    '''Plots the DE-Seq MA-plot using the full regions of interest and saves it to the 
    figuredir directory created in TFEA output folder

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

    x = [x for _,x in sorted(zip(up_p,up_x))] 
        + [x for _,x in sorted(zip(dn_p,dn_x),reverse=True)]

    y = [y for _,y in sorted(zip(up_p,up_y))]
        + [y for _,y in sorted(zip(dn_p,dn_y),reverse=True)]

    c = plt.cm.RdYlGn(np.linspace(0,1,len(x)))
    #Creates an MA-Plot of the region expression
    F = plt.figure(figsize=(7,6))
    ax = plt.subplot(111)
    plt.scatter(x=x,y=y,color=c,edgecolor='')
    ax.set_title("DE-Seq MA-Plot",fontsize=14)
    ax.set_ylabel("Log2 Fold-Change ("+label2+"/"+label1+")",fontsize=14)
    ax.set_xlabel("Log10 Average Expression",fontsize=14)
    ax.tick_params(axis='y', which='both', left='off', right='off', labelleft='on')
    ax.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
    plt.savefig(figuredir + 'DESEQ_MA_Plot.png',bbox_inches='tight')
#===================================================================================

#===================================================================================
def permutations(distances=distances, permutations=1000):
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
#===================================================================================

#===================================================================================
def padj_bonferroni(TFresults=TFresults):
    '''This function iterates through TFEA results, removes TFs that returned 
        "no hits" and calculates a p-adj using the Bonferroni Correction for each 
        TF motif appending it to the given TFresults array

    Parameters
    ----------
    TFresults : list of lists
        contains calculated enrichment scores for all TFs of interest specified by
        the user
        
    Returns
    -------
    TFresults : list of lists
        same as input with an additional p-adjusted value appended to each TF
    '''
    TFresults = [x for x in TFresults if x != "no hits"]
    for i in range(len(TFresults)):                
        #Using Bonferroni Correction
        PADJ = 1 if PVAL*len(TFresults) > 1 else PVAL*len(TFresults)
        TFresults[i].append(PADJ)

    return TFresults
#===================================================================================

#===================================================================================
def motif_distance_bedtools_closest(ranked_center_file=ranked_center_file,
                            motif_path=MOTIF_PATH):
    '''Calculates nearest motif hit from a bed file. TFEA provides this function
        with a bed file containing the center of the inputted regions.

    Parameters
    ----------
    TFresults : list of lists
        contains calculated enrichment scores for all TFs of interest specified by
        the user
        
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

    motif_distance_bed_sorted = '/' + '/'.join(ranked_center_file.split('/')[:-1]) \
                                + '/' + motif_path.split('/')[-1] \
                                + ".sorted.distance.bed"

    return motif_distance_bed_sorted
#===================================================================================

#===================================================================================



