#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''This file contains a list of independent functions that do not call other 
    functions
'''
#==============================================================================
__author__ = 'Jonathan D. Rubin and Rutendo F. Sigauke'
__credits__ = ['Jonathan D. Rubin', 'Rutendo F. Sigauke', 'Jacob T. Stanley',
                'Robin D. Dowell']
__maintainer__ = 'Jonathan D. Rubin'
__email__ = 'Jonathan.Rubin@colorado.edu'
#==============================================================================
import matplotlib
matplotlib.use('Agg')
import os
import math
import datetime
import subprocess
import traceback
import matplotlib.pyplot as plt
# import seaborn as sns
import matplotlib.cm as cm
from matplotlib import gridspec
from scipy.stats import norm 
import numpy as np
import config
#==============================================================================
#Functions
#==============================================================================

#==============================================================================
def merge_bed(beds=None, tempdir=None):
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
    combined_input_merged_bed = os.path.join(tempdir, 
                                                "combined_input.merge.bed")

    os.system("cat " + " ".join(beds) 
                + " | bedtools sort -i stdin | bedtools merge -i stdin > " 
                + combined_input_merged_bed)

    return combined_input_merged_bed
#==============================================================================

#==============================================================================
def tfit_clean_merge(beds=None, tempdir=None, size_cut=500):
    '''Takes in a list of bed files outputted by Tfit and for each region in
        each bed file splits the region based on its genomic size and the 
        size_cut variable. The 'large regions' get merged via bedtools, the 
        'small regions' get instersected via bedtools, then expanded to be
        size_cut in length, then merged with the 'large regions'

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
    #First define a large_regions file and a small_regions list that will
    #store small regions files
    large_regions = os.path.join(tempdir, 'large_regions.bed')
    small_regions_list = list()
    large_regions_file = open(large_regions,'w')

    #Loop through the bed files and accordingly either append regions to 
    # the large regions file or to individual small regions files (append these
    # to the list)
    for bedfile in beds:
        small_regions = os.path.join(tempdir, 
            bedfile.split('/')[-1].split('.bed')[0] + '.small_regions.bed')
        small_regions_list.append(small_regions)
        small_regions_file = open(small_regions,'w')
        with open(bedfile) as F:
            for line in F:
                if '#' not in line:
                    chrom,start,stop = line.strip('\n').split('\t')[:3]
                    start = int(start)
                    stop = int(stop)
                    if stop-start >= size_cut:
                        large_regions_file.write('\t'.join([chrom,str(start),
                                                            str(stop)]) + '\n')
                    else:
                        small_regions_file.write('\t'.join([chrom,str(start),
                                                            str(stop)]) + '\n')
        small_regions_file.close()

    large_regions_file.close()

    #Perform bedtools intersect on the small regions
    command = ("bedtools intersect -a " + small_regions_list[0] + " -b " 
                + small_regions_list[1])
    for bedfile in small_regions_list[2:]:
        command = command + " | bedtools intersect -a stdin -b " + bedfile
    command = (command + " > " 
                + os.path.join(tempdir, "small_regions.intersect.bed"))
    # print command
    os.system(command)

    #Expand the intersected regions so they are size_cut in length
    small_regions_expanded = os.path.join(tempdir, 
                                        "small_regions.intersect.expanded.bed")
    small_regions_expanded_outfile = open(small_regions_expanded,'w')
    with open(os.path.join(tempdir, "small_regions.intersect.bed")) as F:
        for line in F:
            chrom,start,stop = line.strip('\n').split('\t')[:3]
            start = int(start)
            stop = int(stop)
            center = (start+stop)/2
            new_start = center - size_cut/2
            new_stop = center + size_cut/2
            small_regions_expanded_outfile.write('\t'.join([chrom, 
                                                            str(new_start), 
                                                            str(new_stop)])
                                                                + '\n')
    small_regions_expanded_outfile.close()

    #Define the output file
    combined_input_merged_bed = os.path.join(tempdir, 
                                                "combined_input.merge.bed")

    #Concatenate the large and small regions, sort, and merge
    os.system("cat " + large_regions + " " + small_regions_expanded 
                + " | bedtools sort -i stdin | bedtools merge -i stdin > "
                +  combined_input_merged_bed)

    return combined_input_merged_bed
#==============================================================================

#==============================================================================
def get_motif_hits(motif_file=None, header=False):
    '''Counts number of lines in a bed file. Meant to be used with motif hit
        files.

    Parameters
    ----------
    motif_file : string
        full path to a file to be counted

    Returns
    -------
    line_count : int 
        number of lines within the provided motif_file 
    '''
    line_count = subprocess.check_output(["wc", "-l", motif_file]).encode("utf-8")
    line_count = int(line_count.split()[0])
    if header:
        line_count = line_count-1
    return line_count
#==============================================================================

#==============================================================================
def intersect_merge_bed(bed1=None, bed2=None, tempdir=None):
    '''Takes in two lists of bed files, each containing replicates for one 
        condition. Intersects replicates using bedtools then merges intersected
        regions.

    Parameters
    ----------
    bed1 : list or array
        full paths to bed files for condition 1 (strings)
    
    bed2 : list or array
        full paths to bed files for condition 2 (strings)
        
    tempdir : string
        full path to tempdir directory in output directory (created by TFEA)

    Returns
    -------
    combined_input_merged_bed : string 
        full path to a bed file containing the merged regions inputted by the 
        user 
    '''
    #Define the output file
    combined_input_merged_bed = os.path.join(tempdir, 
                                                "combined_input.merge.bed")

    if len(bed1) > 1:
        #Build command to perform bedtools intersect on condition1 beds
        intersect1 = ("bedtools intersect -a " + bed1[0] + " -b " + bed1[1])
        for bedfile in bed1[2:]:
            intersect1 = (intersect1 + " | bedtools intersect -a stdin -b " 
                        + bedfile)
    else:
        intersect1 = "cat " + bed1[0]

    if len(bed2) > 1:
        #Build command to perform bedtools intersect on condition2 beds
        intersect2 = ("bedtools intersect -a " + bed2[0] + " -b " + bed2[1])
        for bedfile in bed2[2:]:
            intersect2 = (intersect2 + " | bedtools intersect -a stdin -b " 
                        + bedfile)
    else:
        intersect2 = "cat " + bed2[0]

    #Build full command which pipes both intersect commands into cat, then 
    # sorts and merges this resulting bed file
    command = ("cat <(" + intersect1 + ") <(" + intersect2 
                + ") | bedtools sort -i stdin | bedtools merge -i stdin > " 
                + combined_input_merged_bed)
    
    #Need to use subprocess here because this command is bash not sh
    subprocess.call(['bash', '-c', command])

    return combined_input_merged_bed
#==============================================================================

#==============================================================================
def getfasta(bedfile=None, genomefasta=None, tempdir=None, outname=None):
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
    os.system("bedtools getfasta -name -fi " + genomefasta 
                + " -bed " + bedfile
                + " -fo " + os.path.join(tempdir, outname))

    ranked_file_fasta = os.path.join(tempdir, outname)

    return ranked_file_fasta
#==============================================================================

#==============================================================================
def smallfasta_from_fasta(fastafile=None, outname=None, 
                            smallwindow=config.SMALLWINDOW,
                            largewindow=config.LARGEWINDOW):
    with open(outname,'w') as outfile:
        with open(fastafile) as F:
            for line in F:
                if '>' in line:
                    outfile.write(line)
                else:
                    line = line.strip('\n')
                    start = int(largewindow-smallwindow)
                    stop = int(largewindow+smallwindow)
                    outfile.write(line[start:stop] + '\n')

    return outname
#==============================================================================

#==============================================================================
def get_bgfile(fastafile=None, tempdir=None):
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
    markov_background = os.path.join(tempdir, "markov_background.txt")
    os.system("fasta-get-markov " + fastafile + " " + markov_background)

    return markov_background
#==============================================================================

#==============================================================================
def get_regions(tempdir=None, ranked_center_file=None, window=None, 
                outname=None):
    '''Takes in a bed file that contains regions centered on the user-inputted 
        bed files and outputs a 'full regions' bed file which simply adds a 
        user-defined window to either side of these centered regions.

    Parameters
    ----------
    tempdir : string
        full path to temp directory in output directory (created by TFEA)

    ranked_center_file : string
        full path to a bed file with centered coordinates
        
    window : float
        half the desired window size

    Returns
    -------
    ranked_full_regions : string 
        full path to a bed file with full regions to be compared with TFEA
    '''
    outfile = open(os.path.join(tempdir, outname),'w')
    with open(ranked_center_file) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            chrom,start,stop = line[:3]
            start = str(int(start)-int(window))
            stop = str(int(stop)+int(window))
            rank, pval, fc = line[3:]
            name = ','.join([rank, pval, fc])
            outfile.write('\t'.join([chrom,start,stop,name]) + '\n')

    ranked_full_regions = os.path.join(tempdir, outname)
    return ranked_full_regions
#==============================================================================

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
    os.system("bedtools multicov -bams " + " ".join(bam1) + " " 
                + " ".join(bam2) + " -bed " + bedfile + " > " 
                + os.path.join(tempdir, "count_file.bed"))

    #This section adds a header to the count_file and reformats it to remove 
    #excess information and add a column with the region for later use
    count_file = os.path.join(tempdir, "count_file.header.bed")
    outfile = open(count_file,'w')
    outfile.write("chrom\tstart\tstop\tregion\t" 
                    + '\t'.join([label1]*len(bam1)) + "\t" 
                    + '\t'.join([label2]*len(bam2)) + "\n")

    with open(os.path.join(tempdir, "count_file.bed")) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            chrom,start,stop = line[:3]
            counts = line[-(len(bam1)+len(bam2)):]
            outfile.write('\t'.join([chrom,start,stop]) + "\t" 
                            + chrom + ":" + start + "-" + stop + "\t"
                            + '\t'.join(counts) + "\n")

    return count_file
#==============================================================================

#==============================================================================
def sum_reads(count_file=None, bam1=config.BAM1, bam2=config.BAM2):
    '''
    '''
    sample_number = (len(bam1)+len(bam2))
    millions_mapped = [0.0]*sample_number
    with open(count_file) as F:
        F.readline()
        for line in F:
            line = line.strip('\n').split('\t')
            [x+float(y) for x,y in zip(millions_mapped,line[-sample_number:])]

    return millions_mapped
#==============================================================================

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
        outfile = open(os.path.join(tempdir, 'DESeq.R'),'w')
        outfile.write('sink("' + os.path.join(tempdir,'DESeq.Rout') + '")\n')
        outfile.write('library("DESeq2")\n')
        outfile.write('data <- read.delim("'+count_file+'", sep="\t", \
                        header=TRUE)\n')
        outfile.write('countsTable <- subset(data, select=c('
                +', '.join([str(i) for i in range(5,5+len(bam1)+len(bam2))])
                +'))\n')

        outfile.write('rownames(countsTable) <- data$region\n')
        outfile.write('conds <- as.data.frame(c(' 
                        + ', '.join(['"'+label1+'"']*len(bam1)) 
                        + ', ' 
                        + ', '.join(['"'+label2+'"']*len(bam2)) 
                        + '))\n')

        outfile.write('colnames(conds) <- c("treatment")\n')
        outfile.write('ddsFullCountTable <- DESeqDataSetFromMatrix(\
                                                    countData = countsTable, \
                                                    colData = conds, \
                                                    design = ~ treatment)\n')

        outfile.write('dds <- DESeq(ddsFullCountTable)\n')
        outfile.write('res1 <- results(dds,alpha = 0.05, \
                                        contrast=c("treatment",\
                                                        "'+label2+'",\
                                                        "'+label1+'"))\
                                                        \n')

        outfile.write('resShrink <- lfcShrink(dds, res = res1, \
                                                contrast = c("treatment",\
                                                "'+label2+'",\
                                                "'+label1+'"))\n')

        outfile.write('resShrink$fc <- 2^(resShrink$log2FoldChange)\n')
        outfile.write('res <- resShrink[c(1:3,7,4:6)]\n')
        outfile.write('write.table(res, file = "'
                        + os.path.join(tempdir,'DESeq.res.txt') 
                        + '", append = FALSE, sep= "\t" )\n')
        outfile.write('sink()')
    else:
        outfile = open(os.path.join(tempdir, 'DESeq.R'),'w')
        outfile.write('sink("'+os.path.join(tempdir, 'DESeq.Rout') + '")\n')
        outfile.write('library("DESeq")\n')
        outfile.write('data <- read.delim("'+count_file+'", sep="\t", \
                        header=TRUE)\n')

        outfile.write('countsTable <- subset(data, select=c('
            +', '.join([str(i) for i in range(5,5+len(bam1)+len(bam2))])
            +'))\n')

        outfile.write('rownames(countsTable) <- data$region\n')
        outfile.write('conds <- c(' + ', '.join(['"'+label1+'"']*len(bam1)) 
                        + ', ' 
                        + ', '.join(['"'+label2+'"']*len(bam2)) 
                        + ')\n')

        outfile.write('cds <- newCountDataSet( countsTable, conds )\n')
        outfile.write('cds <- estimateSizeFactors( cds )\n')
        outfile.write('sizeFactors(cds)\n')                                                               
        outfile.write('cds <- estimateDispersions( cds ,method="blind", \
                        sharingMode="fit-only")\n')

        outfile.write('res <- nbinomTest( cds, "'+label1+'", "'+label2+'" )\n')
        outfile.write('rownames(res) <- res$id\n')                      
        outfile.write('write.table(res, file = "'
                        + os.path.join(tempdir,'DESeq.res.txt') 
                        + '", append = FALSE, sep= "\t" )\n')

        outfile.write('sink()')
    outfile.close()
#==============================================================================

#==============================================================================
def plot_deseq_MA(deseq_file=None, label1=None, label2=None, figuredir=None):
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
    plt.savefig(os.path.join(figuredir, 'DESEQ_MA_Plot.png'),
                bbox_inches='tight')
#==============================================================================

#==============================================================================
def permutations_auc(distances=None, trend=None, permutations=1000):
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
    triangle_area = np.trapz(trend)
    for i in range(permutations):
        random_distances = np.random.permutation(distances)
        cum_distances = np.cumsum(random_distances)
        es = np.trapz(cum_distances)
        auc = es - triangle_area
        es_permute.append(auc)

    return es_permute
#==============================================================================

#==============================================================================
def permutations_youden_rank(distances=None, trend=None, permutations=1000):
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
    for i in range(permutations):
        random_distances = np.random.permutation(distances)
        cumscore = np.cumsum(random_distances)
        youden = cumscore - trend
        youden_max = max(cumscore - trend)
        rank_max = np.where(youden == youden_max)[0][0]/float(len(youden))
        es_permute.append(rank_max)

    return es_permute
#==============================================================================

#==============================================================================
def pvalue_global_youden_rank(TFresults=None, permutations=1000):
    '''Calculates a p-value based on a distribution of youden_ranks

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
    all_youden = [motif[1] for motif in TFresults]
    ##significance calculator                                                                                                                                                            
    mu = np.mean(all_youden)
    sigma = np.std(all_youden)
    for i in range(len(TFresults)):
        actual_youden_rank = TFresults[i][1]
        p = min(norm.cdf(actual_youden_rank,mu,sigma), 
                            1-norm.cdf(actual_youden_rank,mu,sigma))
        TFresults[i].append(p)
#==============================================================================

#==============================================================================
def padj_bonferroni(TFresults=None):
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
        PVAL = TFresults[i][-2]
        #Using Bonferroni Correction
        PADJ = 1 if PVAL*len(TFresults) > 1 else PVAL*len(TFresults)
        TFresults[i].append(PADJ)

    return TFresults
#==============================================================================

#==============================================================================
def motif_distance_bedtools_closest(ranked_center_file=None, motif_path=None):
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
    motif_distance_bed_sorted = '/'.join(ranked_center_file.split('/')[:-1])\
                                + '/' + motif_path.split('/')[-1]\
                                + ".sorted.distance.bed"

    os.system("bedtools closest -D ref -t first -a " 
                + ranked_center_file + " -b " 
                + motif_path + " > " + motif_distance_bed_sorted)
    

    return motif_distance_bed_sorted
#==============================================================================

#==============================================================================
defget_motif_names (motifdatabase=None):
    '''Extracts motif names from a MEME formatted motif database

    Parameters
    ----------
    motifdatabase : string
        full path to a meme formatted file with motifs to be used in TFEA
        
    Returns
    -------
    motif_list : list
        a list of motif names to be analyzed in TFEA
    '''
    motif_list = list()
    with open(motifdatabase) as F:
        for line in F:
            if 'MOTIF' in line:
                line = line.strip('\n').split()
                motif_name = line[-1]
                motif_list.append(motif_name+'.bed')

    return motif_list
#==============================================================================

#==============================================================================
def fimo(motif=None, bg_file=None, ranked_fasta_file=None, tempdir=None, 
        motifdatabase=None, thresh=config.FIMO_THRESH):
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
    fimo_out = os.path.join(tempdir, motif)
    command = ("fimo --skip-matched-sequence --verbosity 1 --thresh " + str(thresh) 
                + " --bgfile " + bg_file 
                + " --motif " + motif.strip('.bed') + " " + motifdatabase 
                + " " + ranked_fasta_file
                + " > " + fimo_out)
    with open(os.devnull, 'w') as devnull:
        subprocess.call(command, shell=True, stderr=devnull)

    return fimo_out
#==============================================================================

#==============================================================================
def homer(tempdir=None, srcdirectory=None, fastafile=None, motif_file=None, 
            cpus=None):
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
    homer_out = os.path.join(tempdir,"homer_out.txt")
    command = ("homer2 find -p " + str(cpus) + " -i " + fastafile 
        + " -m " 
        + motif_file
        + " > " + homer_out)
    # print command
    os.system(command)
    # with open(os.devnull, 'w') as devnull:
    #     subprocess.call(command, shell=True, stderr=devnull)

    return homer_out
#==============================================================================

#==============================================================================
def homer_parse(largewindow=None, tempdir=None, homer_out=None,
                ranked_center_file=None):
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
    motif_names = list()
    motif_files = list()
    with open(homer_out) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            pval, fc, rank = line[0].split(',')
            distance = line[1]
            score = line[-1]
            motif_name = line[-3].split('/')[0]
            if motif_name not in motif_names:
                motif_names.append(motif_name)
                motif_files.append(open(os.path.join(tempdir,motif_name+'.txt'),'w'))
            motif_files[motif_names.index(motif_name)].write('\t'.join([score, 
                                                        '', '', pval, fc, rank, 
                                                        str(distance)]) + '\n')
    
    for file1 in motif_files:
        file1.close()

    motif_list = list()
    for motif_name in motif_names:
        temp_lines = dict()
        hits = 0
        motif_name = os.path.join(tempdir,motif_name+'.txt')
        with open(motif_name) as F:
            for line in F:
                line = line.strip('\n').split('\t')
                score = line[0]
                pval, fc, rank, distance = line[3:]
                rank = int(rank)
                temp_lines[rank] = (rank, pval, fc, score, distance)
                hits+=1
        motif_list.append((motif_name, hits))
        outname = os.path.join(tempdir,motif_name+'.sorted.distance.bed')
        outfile = open(outname, 'w')

        with open(ranked_center_file, 'r') as F:
            for line in F:
                line = line.strip('\n').split('\t')
                pval, fc, rank = line[3:]
                if rank in temp_lines:
                    rank, pval, fc, score, distance = temp_lines[rank]
                    outfile.write('\t'.join([score, '', '', pval, fc, rank, 
                                    str(distance)]) + '\n')
                else:
                    outfile.write('\t'.join(
                                        [str(0.0), '', '', pval, fc, rank, 
                                        str(int(largewindow**2))]) + '\n')
        outfile.close()

    return motif_list
#==============================================================================

#==============================================================================
def meme2images(motifdatabase=None, figuredir=None, motif=None):
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
    command = ("meme2images -rc --motif " + motif.strip('.bed') + " " + motifdatabase
              + " " + figuredir)
    # os.system(command)
    with open(os.devnull) as devnull:
        subprocess.call(command, shell=True, stderr=devnull)
#==============================================================================

#==============================================================================
def fasta_markov(tempdir=None, fastafile=None, order='2'):
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
    markov_background = os.path.join(tempdir, "markov_background.txt")
    os.system("fasta-get-markov -m " + order + " " + fastafile +
              " > " + markov_background)

    return markov_background
#==============================================================================

#==============================================================================
def fimo_parse(largewindow=None, tempdir=None, fimo_file=None,
                ranked_center_file=None):
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
        header = F.readline().strip('\n').split('\t')
        if len(header) > 1:
            try:
                start_index = header.index('start')
                stop_index = header.index('stop')
                name_index = header.index('sequence_name')
                score_index = header.index('score')
            except:
                pass
            for line in F:
                line = line.strip('\n').split('\t')
                pval, fc, rank = line[name_index].split(',')
                start = line[start_index]
                stop = line[stop_index]
                distance = ((int(start)+int(stop))/2)-int(largewindow)
                score = line[score_index]
                if rank not in d:
                    d[rank] = [rank, pval, fc, score, distance]
                else:
                    prev_score = float(d[rank][-2])
                    if prev_score < float(score):
                        d[rank] = [rank, pval, fc, score, distance]
    
    outname = fimo_file.strip('.bed')+'.sorted.distance.bed'
    outfile = open(outname, 'w')
    with open(ranked_center_file, 'r') as F:
        for line in F:
            line = line.strip('\n').split('\t')
            pval, fc, rank = line[3:]
            if rank in d:
                rank, pval, fc, score, distance = d[rank]
                outfile.write('\t'.join([score, '', '', pval, fc, rank, 
                                str(distance)]) + '\n')
            else:
                outfile.write('\t'.join(
                                    [str(0.0), '', '', pval, fc, rank, 
                                    str(int(largewindow)*10)]) + '\n')
    outfile.close()

    return outname
#==============================================================================

#==============================================================================
def convert_sequence_to_array(sequence=None):
    '''This function converts a DNA sequence (ACGT alphabet) to an array, 
    collapsing GCs to 1's and ATs to 0's

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
            array.append(0.0)

    return array 
#==============================================================================

#==============================================================================
def rank_deseqfile(deseq_file=None, tempdir=None):
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
    outfile = open(os.path.join(tempdir, "ranked_file.bed"),'w')
    r=1
    for region in sorted(up, key=lambda x: x[3]):
        outfile.write('\t'.join(region) + '\t' + str(r) + '\n')
        r += 1
    for region in sorted(down, key=lambda x: x[3], reverse=True):
        outfile.write('\t'.join(region) + '\t' + str(r) + '\n')
        r += 1
    outfile.close()

    #Get center base for each region
    outfile = open(os.path.join(tempdir, "ranked_file.center.bed"),'w')
    with open(os.path.join(tempdir, "ranked_file.bed")) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            chrom,start,stop = line[:3]
            center = (int(start)+int(stop))/2
            outfile.write(chrom + '\t' + str(center) + '\t' + str(center+1) 
                            + '\t' + '\t'.join(line[3:]) + '\n')
    outfile.close()


    os.system("sort -k1,1 -k2,2n " 
                + os.path.join(tempdir, "ranked_file.center.bed") 
                + " > " 
                + os.path.join(tempdir, "ranked_file.center.sorted.bed"))

    ranked_center_file = os.path.join(tempdir, "ranked_file.center.sorted.bed")

    return ranked_center_file
#==============================================================================

#==============================================================================
def samtools_flagstat(args):
    '''Performs samtools flagstat on a bam file. Then parses the samtools 
        flagstat output and returns a the millions mapped reads for the given
        bam file.

    Parameters
    ----------
    args : tuple
        contains arguments for this function:
        bam : string
            full paths to a bam files
        tempdir : string
            full path to the tempdir directory within the output directory 
            (created by TFEA)
        
    Returns
    -------
    millions_mapped : float
        millions mapped reads for the given bam file
    '''
    bam, tempdir = args
    filename = bam.split('/')[-1]
    os.system("samtools flagstat " + bam + " > " 
                + os.path.join(tempdir, filename + ".flagstat"))
    with open(os.path.join(tempdir,filename+".flagstat")) as F:
        lines = F.readlines()
        millions_mapped = float(lines[4].strip('\n').split(' ')[0])/1000000.0

        return millions_mapped
#==============================================================================

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
    import HTSeq as hts
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

    print millions_mapped


    posprofile1 = np.zeros(2*int(largewindow))  
    negprofile1 = np.zeros(2*int(largewindow))
    posprofile2 = np.zeros(2*int(largewindow))   
    negprofile2 = np.zeros(2*int(largewindow))
    rep1number = float(len(hts_bam1))
    rep2number = float(len(hts_bam2))
    for window in regions:
        avgposprofile1 = np.zeros(2*int(largewindow))
        avgnegprofile1 = np.zeros(2*int(largewindow))
        i = 0
        for sortedbamfile in hts_bam1:
            mil_map = float(millions_mapped[i])/1000000.0
            i += 1
            tempposprofile = np.zeros(2*int(largewindow))
            tempnegprofile = np.zeros(2*int(largewindow))
            for almnt in sortedbamfile[ window ]:
                if almnt.iv.strand == '+':
                    start_in_window = almnt.iv.start - window.start
                    end_in_window = almnt.iv.end - window.end \
                                    + 2*int(largewindow)
                    start_in_window = max( start_in_window, 0 )
                    end_in_window = min( end_in_window, 2*int(largewindow) )
                    tempposprofile[ start_in_window : end_in_window ] += 1.0/mil_map
                if almnt.iv.strand == '-':
                    start_in_window = almnt.iv.start - window.start
                    end_in_window   = almnt.iv.end - window.end \
                                        + 2*int(largewindow)
                    start_in_window = max( start_in_window, 0 )
                    end_in_window = min( end_in_window, 2*int(largewindow) )
                    tempnegprofile[ start_in_window : end_in_window ] += -1.0/mil_map
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
            mil_map = float(millions_mapped[i])/1000000.0
            i += 1
            tempposprofile = np.zeros(2*int(largewindow))
            tempnegprofile = np.zeros(2*int(largewindow))
            for almnt in sortedbamfile[ window ]:
                if almnt.iv.strand == '+':
                    start_in_window = almnt.iv.start - window.start
                    end_in_window   = almnt.iv.end - window.end \
                                        + 2*int(largewindow)
                    start_in_window = max( start_in_window, 0 )
                    end_in_window = min( end_in_window, 2*int(largewindow) )
                    tempposprofile[ start_in_window : end_in_window ] += 1.0/mil_map
                if almnt.iv.strand == '-':
                    start_in_window = almnt.iv.start - window.start
                    end_in_window   = almnt.iv.end - window.end \
                                        + 2*int(largewindow)
                    start_in_window = max( start_in_window, 0 )
                    end_in_window = min( end_in_window, 2*int(largewindow) )
                    tempnegprofile[ start_in_window : end_in_window ] += -1.0/mil_map
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
    
    # mil_map1 = mil_map1/rep1number
    # mil_map2 = mil_map2/rep2number

    # posprofile1 = [x/mil_map1 for x in posprofile1]
    # negprofile1 = [x/mil_map1 for x in negprofile1]
    # posprofile2 = [x/mil_map2 for x in posprofile2]
    # negprofile2 = [x/mil_map2 for x in negprofile2]


    return posprofile1, negprofile1, posprofile2, negprofile2
#==============================================================================

#==============================================================================
def meta_profile2(regionlist=None, region_file=None, millions_mapped=None, 
                largewindow=None, bam1=None, bam2=None):

    import HTSeq
    bamfile = HTSeq.BAM_Reader( "SRR001432_head.bam" )
    gtffile = HTSeq.GFF_Reader( "Homo_sapiens.GRCh37.56_chrom1.gtf" )
    halfwinwidth = largewindow
    fragmentsize = 200

    tsspos = HTSeq.GenomicArrayOfSets( "auto", stranded=False )
    for feature in gtffile:
        if feature.type == "exon" and feature.attr["exon_number"] == "1":
            p = feature.iv.start_d_as_pos
            window = HTSeq.GenomicInterval( p.chrom, p.pos - halfwinwidth, p.pos + halfwinwidth, "." )
            tsspos[ window ] += p

    profile = np.zeros( 2*halfwinwidth, dtype="i" )
    for almnt in bamfile:
        if almnt.aligned:
            almnt.iv.length = fragmentsize
            s = set()
            for step_iv, step_set in tsspos[ almnt.iv ].steps():
                s |= step_set
            for p in s:
                if p.strand == "+":
                    start_in_window = almnt.iv.start - p.pos + halfwinwidth
                    end_in_window   = almnt.iv.end   - p.pos + halfwinwidth
                else:
                    start_in_window = p.pos + halfwinwidth - almnt.iv.end
                    end_in_window   = p.pos + halfwinwidth - almnt.iv.start
                start_in_window = max( start_in_window, 0 )
                end_in_window = min( end_in_window, 2*halfwinwidth )
                profile[ start_in_window : end_in_window ] += 1

    import HTSeq as hts
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
            mil_map = float(millions_mapped[i])
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
        avgposprofile1 = [x/rep1number/mil_map for x in avgposprofile1]
        avgnegprofile1 = [x/rep1number/mil_map for x in avgnegprofile1]
        # avgposprofile1 = [x/rep1number for x in avgposprofile1]
        # avgnegprofile1 = [x/rep1number for x in avgnegprofile1]
        posprofile1 = [x+y for x,y in zip(posprofile1,avgposprofile1)]
        negprofile1 = [x+y for x,y in zip(negprofile1, avgnegprofile1)]

        avgposprofile2 = np.zeros(2*int(largewindow))
        avgnegprofile2 = np.zeros(2*int(largewindow))
        i = len(hts_bam1)
        for sortedbamfile in hts_bam2:
            mil_map = float(millions_mapped[i])
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
        avgposprofile2 = [x/rep2number/mil_map for x in avgposprofile2]
        avgnegprofile2 = [x/rep2number/mil_map for x in avgnegprofile2]
        # avgposprofile2 = [x/rep2number for x in avgposprofile2]
        # avgnegprofile2 = [x/rep2number for x in avgnegprofile2]
        posprofile2 = [x+y for x,y in zip(posprofile2,avgposprofile2)]
        negprofile2 = [x+y for x,y in zip(negprofile2, avgnegprofile2)]
    
    # mil_map1 = mil_map1/rep1number
    # mil_map2 = mil_map2/rep2number

    # posprofile1 = [x/mil_map1 for x in posprofile1]
    # negprofile1 = [x/mil_map1 for x in negprofile1]
    # posprofile2 = [x/mil_map2 for x in posprofile2]
    # negprofile2 = [x/mil_map2 for x in negprofile2]


    return posprofile1, negprofile1, posprofile2, negprofile2
#==============================================================================

#==============================================================================
def enrichment_plot(largewindow=None, smallwindow=None, figuredir=None,
                    cumscore=None, sorted_distances=None, logpval=None, 
                    updistancehist=None, downdistancehist=None, 
                    gc_array=None, motif_file=None, dpi=None, save=True, 
                    score=None, q1_distances=None, q2_distances=None, 
                    q3_distances=None, q4_distances=None, 
                    meta_profile_dict=None, label1=None, label2=None):
    '''This function plots the TFEA enrichment plot.

    Parameters
    ----------
    largewindow : float
        a user specified larger window used for plotting purposes and to do
        some calculations regarding the user-provided regions of interest

    smallwindow : float
        a user specified smaller window used for plotting purposes and to do
        some calculations regarding the user-provided regions of interest

    figuredir : string
        the full path to the figuredir within the ouptut directory containing
        all figures and plots

    cumscore : list
        the cumulative score of the running sum as we walk through the ranked
        regions

    sorted_distances : list or array
        a sorted (based on rank) list of motif distances. Negative corresponds
        to upstream of region

    logpval : list or array
        a way to visualize the ranking of the regions of interest. It is
        the log10 of the p-value with the sign (positive or negative) based on
        whether the fold change of the region is over 1 or less than 1.

    updistancehist : list or array
        the first quartile of ranked regions. These are presumably higher in
        condition1

    downdistancehist : list or array
        the fourth quartile of ranked regions. These are presumably higher in
        condition2

    gc_array : list
        an array of gc richness of regions of interest. It is recommended that
        this array be no larger than 1000 bins.

    motif_file : string
        the name of the motif thats associated with all the input data. Used 
        for figure naming purposes.

    Returns
    -------
    None
    '''
    import config
    from scipy.stats import gaussian_kde
    dpi = config.DPI
    try:
        #Begin plotting section
        len_cumscore = float(len(cumscore))
        F = plt.figure(figsize=(15.5,12))
        # xvals = range(0, int(len_cumscore))
        xvals = np.linspace(start=0, stop=1, num=len_cumscore)
        limits = [0, 1]

        #With GC-Content
        # gs = gridspec.GridSpec(4, 1, height_ratios=[2, 2, 1, 1])

        #Without GC-Content
        # gs = gridspec.GridSpec(3, 1, height_ratios=[3, 1, 2])

        outer_gs = gridspec.GridSpec(2, 1, height_ratios=[2,1])
        enrichment_gs = gridspec.GridSpecFromSubplotSpec(4, 1, 
                                            subplot_spec=outer_gs[0], 
                                            height_ratios=[4, 1, 4, 2], 
                                            hspace=.1)
        lineplot = plt.subplot(enrichment_gs[0])
        barplot = plt.subplot(enrichment_gs[1])
        scatterplot = plt.subplot(enrichment_gs[2])
        fillplot = plt.subplot(enrichment_gs[3])

        #This is the enrichment score plot (i.e. line plot)
        # ax0 = plt.subplot(gs[0])
        lineplot.plot(xvals,cumscore,color='green')
        lineplot.plot([0, 1],[0, 1], '--', alpha=0.75)
        lineplot.set_title(motif_file.split('.bed')[0] + ' Enrichment Plot',
                        fontsize=14)
        lineplot.set_ylabel('Enrichment Score (ES)', fontsize=10)
        lineplot.tick_params(axis='y', which='both', left='on', right='off', 
                        labelleft='on')
        lineplot.tick_params(axis='x', which='both', bottom='off', top='off', 
                        labelbottom='off')
        lineplot.set_ylim([0,1])
        lineplot.set_xlim(limits)

        #This is the barplot right below the enrichment score line plot
        x = xvals
        scatterplot.scatter(x, sorted_distances, edgecolor="", color="black", 
                        s=10, alpha=0.25)
        scatterplot.tick_params(axis='y', which='both', left='off', right='off', 
                        labelleft='on') 
        scatterplot.tick_params(axis='x', which='both', bottom='off', top='off', 
                        labelbottom='off')
        plt.yticks([-int(largewindow),0,int(largewindow)],
                    [str(-int(largewindow)/1000.0),'0',\
                    str(int(largewindow)/1000.0)])
        scatterplot.set_xlim(limits)
        scatterplot.set_ylim([-int(largewindow),int(largewindow)])
        scatterplot.set_ylabel('Distance (kb)', fontsize=10)

        norm    = matplotlib.colors.Normalize(vmin=min(score), vmax=max(score))
        cmap    = cm.Greys
        m       = cm.ScalarMappable(norm=norm, cmap=cmap)
        colors  = [m.to_rgba(c) for c in score] 
        barplot.bar(x,[1 for l in x], edgecolor="", color=colors)
        barplot.tick_params(axis='y', which='both', left='off', right='off', 
                        labelleft='off') 
        barplot.tick_params(axis='x', which='both', bottom='off', top='off', 
                        labelbottom='off')
        barplot.set_xlim(limits)
        barplot.set_ylim([0,1])
        barplot.set_ylabel('Score', fontsize=10)

        #This is the rank metric fill plot
        fillplot.fill_between(xvals,0,logpval,facecolor='grey',edgecolor="")
        fillplot.tick_params(axis='y', which='both', left='on', right='off', 
                        labelleft='on')
        fillplot.tick_params(axis='x', which='both', bottom='off', top='off', 
                        labelbottom='on')
        ylim = math.fabs(max([x for x in logpval if -500 < x < 500],key=abs))
        fillplot.set_ylim([-ylim,ylim])
        fillplot.yaxis.set_ticks([int(-ylim),0,int(ylim)])
        fillplot.set_xlim(limits)
        fillplot.set_xlabel('Relative Rank (n='+str(int(len_cumscore))+')', 
                        fontsize=14)
        fillplot.set_ylabel('Rank Metric',fontsize=10)
        # try:
        #     ax2.axvline(len(updistancehist)+1,color='green',alpha=0.25)
        # except ValueError:
        #     pass
        # try:
        #     ax2.axvline(len(xvals) - len(downdistancehist), color='purple', 
        #                 alpha=0.25)
        # except ValueError:
        #     pass



        #This is the GC content plot
        # ax3 = plt.subplot(gs[2])
        # ax3.set_xlim(limits)
        # # plt.imshow(gc_array, cmap='hot', interpolation='nearest')
        # sns.heatmap(gc_array, cbar=False, xticklabels='auto',
        #             yticklabels='auto')

        # plt.yticks([0,int(largewindow),int(largewindow*2)],
        #             [str(-int(largewindow)/1000.0),'0',\
        #             str(int(largewindow)/1000.0)])

        # ax3.tick_params(axis='y', which='both', left='on', right='off', 
        #                 labelleft='on')

        # ax3.tick_params(axis='x', which='both', bottom='off', top='off', 
        #                 labelbottom='off')

        # ax3.set_ylabel('GC content per kb',fontsize=10)

        meta_gs = gridspec.GridSpecFromSubplotSpec(3, 4, 
                            subplot_spec=outer_gs[1], hspace=0.1, wspace=0.1)

        #Initiate meta plots
        ax3 = plt.subplot(meta_gs[:2, 0])
        ax4 = plt.subplot(meta_gs[:2, 1])
        ax5 = plt.subplot(meta_gs[:2, 2])
        ax6 = plt.subplot(meta_gs[:2, 3])

        #Initiate heatmaps
        ax7 = plt.subplot(meta_gs[2, 0])
        ax8 = plt.subplot(meta_gs[2, 1])
        ax9 = plt.subplot(meta_gs[2, 2])
        ax10 = plt.subplot(meta_gs[2, 3])

        if config.METAPLOT:
            xvals = range(-int(largewindow),int(largewindow))
            q1posprofile1 = meta_profile_dict['q1posprofile1']
            q1negprofile1 = meta_profile_dict['q1negprofile1']
            q1posprofile2 = meta_profile_dict['q1posprofile2']
            q1negprofile2 = meta_profile_dict['q1negprofile2']
            q2posprofile1 = meta_profile_dict['q2posprofile1']
            q2negprofile1 = meta_profile_dict['q2negprofile1']
            q2posprofile2 = meta_profile_dict['q2posprofile2']
            q2negprofile2 = meta_profile_dict['q2negprofile2']
            q3posprofile1 = meta_profile_dict['q3posprofile1']
            q3negprofile1 = meta_profile_dict['q3negprofile1']
            q3posprofile2 = meta_profile_dict['q3posprofile2']
            q3negprofile2 = meta_profile_dict['q3negprofile2']
            q4posprofile1 = meta_profile_dict['q4posprofile1']
            q4negprofile1 = meta_profile_dict['q4negprofile1']
            q4posprofile2 = meta_profile_dict['q4posprofile2']
            q4negprofile2 = meta_profile_dict['q4negprofile2']
            ylim = [min(q1negprofile1+q1negprofile2+q2negprofile1+q2negprofile2
                        +q3negprofile1+q3negprofile2+q4negprofile1+q4negprofile2),
                    max(q1posprofile1+q1posprofile2+q2posprofile1+q2posprofile2
                        +q3posprofile1+q3posprofile2+q4posprofile1+q4posprofile2)]

            # First quartile plot
            line1, = ax3.plot(xvals,q1posprofile1,color='blue',label=label1)
            ax3.plot(xvals,q1negprofile1,color='blue')
            line2, = ax3.plot(xvals,q1posprofile2,color='red',label=label2)
            ax3.plot(xvals,q1negprofile2,color='red')
            ax3.legend(loc=2,fontsize='small')
            ax3.set_title('Q1',fontsize=14)
            ax3.tick_params(axis='y', which='both', left='off', right='off', 
                            labelleft='on')
            ax3.tick_params(axis='x', which='both', bottom='off', top='off', 
                            labelbottom='off')
            ax3.set_ylabel('Reads per Millions Mapped',fontsize=10)
            plt.yticks([-int(largewindow),0,int(largewindow)],
                        [str(-int(largewindow)/1000.0),'0',\
                        str(int(largewindow)/1000.0)])
            ax3.set_ylim(ylim)


            # Second quartile plot
            line1, = ax4.plot(xvals,q2posprofile1,color='blue',label=label1)
            ax4.plot(xvals,q2negprofile1,color='blue')
            line2, = ax4.plot(xvals,q2posprofile2,color='red',label=label2)
            ax4.plot(xvals,q2negprofile2,color='red')
            # ax2.legend(loc=1)
            ax4.set_title('Q2',fontsize=14)
            ax4.tick_params(axis='y', which='both', left='off', right='off', 
                            labelleft='off')
            ax4.tick_params(axis='x', which='both', bottom='off', top='off', 
                            labelbottom='off')
            ax4.set_ylim(ylim)

            # Third quartile plot
            line1, = ax5.plot(xvals,q3posprofile1,color='blue',label=label1)
            ax5.plot(xvals,q3negprofile1,color='blue')
            line2, = ax5.plot(xvals,q3posprofile2,color='red',label=label2)
            ax5.plot(xvals,q3negprofile2,color='red')
            # ax3.legend(loc=1)
            ax5.set_title('Q3',fontsize=14)
            ax5.tick_params(axis='y', which='both', left='off', right='off', 
                            labelleft='off')
            ax5.tick_params(axis='x', which='both', bottom='off', top='off', 
                            labelbottom='off')
            ax5.set_ylim(ylim)

            # Fourth quartile plot
            line1, = ax6.plot(xvals,q4posprofile1,color='blue',label=label1)
            ax6.plot(xvals,q4negprofile1,color='blue')
            line2, = ax6.plot(xvals,q4posprofile2,color='red',label=label2)
            ax6.plot(xvals,q4negprofile2,color='red')
            # ax4.legend(loc=1)
            ax6.set_title('Q4',fontsize=14)
            ax6.tick_params(axis='y', which='both', left='off', right='off', 
                            labelleft='off')
            ax6.tick_params(axis='x', which='both', bottom='off', top='off', 
                            labelbottom='off')
            ax6.set_ylim(ylim)


        #Distance distribution heatmaps
        bins = 100
        xlim = [int(-largewindow), int(largewindow)]
        counts, edges = np.histogram(q1_distances, bins=bins)
        edges = (edges[1:]+edges[:-1])/2.0
        norm    = matplotlib.colors.Normalize(vmin=min(counts), 
                                                vmax=max(counts))
        cmap    = cm.YlOrRd
        m       = cm.ScalarMappable(norm=norm, cmap=cmap)
        colors  = [m.to_rgba(c) for c in counts] 
        
        ax7.bar(edges,np.ones((len(edges),)), color=colors, 
                    width=(edges[-1]-edges[0])/len(edges), edgecolor=colors)
        ax7.set_ylim([0,1])
        ax7.set_xlim(xlim)
        ax7.tick_params(axis='y', which='both', left='off', right='off', 
                        labelleft='off') 
        ax7.tick_params(axis='x', which='both', bottom='off', top='off', 
                        labelbottom='on')
        ax7.set_xlabel('Motif to Region Center Distance (bp)')


        counts,edges = np.histogram(q2_distances, bins=bins)
        edges        = (edges[1:]+edges[:-1])/2. 
        norm    = matplotlib.colors.Normalize(vmin=min(counts), 
                                                vmax=max(counts))
        cmap    = cm.YlOrRd
        m       = cm.ScalarMappable(norm=norm, cmap=cmap)
        colors  = [m.to_rgba(c) for c in counts] 
        
        ax8.bar(edges,np.ones((len(edges),)), color=colors, 
                    width=(edges[-1]-edges[0])/len(edges), edgecolor=colors)
        ax8.set_ylim([0,1])
        ax8.set_xlim(xlim)
        ax8.tick_params(axis='y', which='both', left='off', right='off', 
                        labelleft='off') 
        ax8.tick_params(axis='x', which='both', bottom='off', top='off', 
                        labelbottom='on')
        ax8.set_xlabel('Motif to Region Center Distance (bp)')


        counts,edges = np.histogram(q3_distances, bins=bins)
        edges        = (edges[1:]+edges[:-1])/2. 
        norm    = matplotlib.colors.Normalize(vmin=min(counts), 
                                                vmax=max(counts))
        cmap    = cm.YlOrRd
        m       = cm.ScalarMappable(norm=norm, cmap=cmap)
        colors  = [m.to_rgba(c) for c in counts] 
        
        ax9.bar(edges,np.ones((len(edges),)), color=colors, 
                    width=(edges[-1]-edges[0])/len(edges), edgecolor=colors)
        ax9.set_ylim([0,1])
        ax9.set_xlim(xlim)
        ax9.tick_params(axis='y', which='both', left='off', right='off', 
                        labelleft='off') 
        ax9.tick_params(axis='x', which='both', bottom='off', top='off', 
                         labelbottom='on')
        ax9.set_xlabel('Motif to Region Center Distance (bp)')


        counts,edges = np.histogram(q4_distances, bins=bins)
        edges        = (edges[1:]+edges[:-1])/2. 
        norm    = matplotlib.colors.Normalize(vmin=min(counts), 
                                                vmax=max(counts))
        cmap    = cm.YlOrRd
        m       = cm.ScalarMappable(norm=norm, cmap=cmap)
        colors  = [m.to_rgba(c) for c in counts] 
        
        ax10.bar(edges,np.ones((len(edges),)), color=colors, 
                    width=(edges[-1]-edges[0])/len(edges), edgecolor=colors)
        ax10.set_ylim([0,1])
        ax10.set_xlim(xlim)
        ax10.tick_params(axis='y', which='both', left='off', right='off', 
                        labelleft='off') 
        ax10.tick_params(axis='x', which='both', bottom='off', top='off', 
                        labelbottom='on')
        ax10.set_xlabel('Motif Distance (kb)')

        counts,edges = np.histogram(q4_distances, bins=bins)
        edges        = (edges[1:]+edges[:-1])/2. 

        #Makes sure that axis labels are properly spaced
        plt.tight_layout()

        if save:
            plt.savefig(os.path.join(figuredir, motif_file
                    + '_enrichment_plot.png'), dpi=dpi, bbox_inches='tight')
            plt.close()
        else:
            plt.show()

    except Exception as e:
        # This prints the type, value, and stack trace of the
        # current exception being handled.
        print traceback.print_exc()
        raise e
#==============================================================================

#==============================================================================
def distance_heatmap_plot(figuredir=None, motif_file=None, q1_distances=None, 
                        q2_distances=None, q3_distances=None, 
                        q4_distances=None, bins=None, dpi=None):
    '''
    '''
    F = plt.figure(figsize=(15.5,1))

    # ax0 = plt.subplot(141)
    # norm    = matplotlib.colors.Normalize(vmin=min(score), vmax=max(score))
    # cmap    = cm.Greys
    # m       = cm.ScalarMappable(norm=norm, cmap=cmap)
    # colors  = [m.to_rgba(c) for c in score] 
    # ax0.bar(x,[1 for l in x], edgecolor="", color=colors)
    # ax0.tick_params(axis='y', which='both', left='off', right='off', 
    #                 labelleft='off') 
    # ax0.tick_params(axis='x', which='both', bottom='off', top='off', 
    #                 labelbottom='off')
    # ax0.set_xlim(limits)
    # ax0.set_ylim([0,1])
    # ax0.set_ylabel('Score', fontsize=10)

    counts,edges = np.histogram(q1_distances, bins=bins)
    edges        = (edges[1:]+edges[:-1])/2. 
    # plt.bar(edges, counts, width=(edges[-1]-edges[0])/bins)

    norm    = matplotlib.colors.Normalize(vmin=min(counts), vmax=max(counts))
    cmap    = cm.YlOrRd
    m       = cm.ScalarMappable(norm=norm, cmap=cmap)
    colors  = [m.to_rgba(c) for c in counts] 
    
    ax0 = F.add_subplot(141)
    ax0.bar(edges,np.ones((len(edges),)), color=colors, 
                    width=(edges[-1]-edges[0])/len(edges) , edgecolor=colors)
    ax0.set_ylim([0,1])
    # ax0.set_xlim(-1, 1)
    ax0.tick_params(axis='y', which='both', left='off', right='off', 
                    labelleft='off') 
    ax0.tick_params(axis='x', which='both', bottom='off', top='off', 
                    labelbottom='on')
    ax0.set_xlabel('Motif Distance (kb)')


    counts,edges = np.histogram(q2_distances, bins=bins)
    edges        = (edges[1:]+edges[:-1])/2. 
    # plt.bar(edges, counts, width=(edges[-1]-edges[0])/bins)

    norm    = matplotlib.colors.Normalize(vmin=min(counts), vmax=max(counts))
    cmap    = cm.YlOrRd
    m       = cm.ScalarMappable(norm=norm, cmap=cmap)
    colors  = [m.to_rgba(c) for c in counts] 
    
    ax1 = F.add_subplot(141)
    ax1.bar(edges,np.ones((len(edges),)), color=colors, 
                    width=(edges[-1]-edges[0])/len(edges) , edgecolor=colors)
    ax1.set_ylim([0,1])
    ax1.set_xlim(-1, 1)
    ax1.tick_params(axis='y', which='both', left='off', right='off', 
                    labelleft='off') 
    ax1.tick_params(axis='x', which='both', bottom='off', top='off', 
                    labelbottom='on')
    ax1.set_xlabel('Motif Distance (kb)')


    counts,edges = np.histogram(q3_distances, bins=bins)
    edges        = (edges[1:]+edges[:-1])/2. 
    # plt.bar(edges, counts, width=(edges[-1]-edges[0])/bins)

    norm    = matplotlib.colors.Normalize(vmin=min(counts), vmax=max(counts))
    cmap    = cm.YlOrRd
    m       = cm.ScalarMappable(norm=norm, cmap=cmap)
    colors  = [m.to_rgba(c) for c in counts] 
    
    ax2 = F.add_subplot(142)
    ax2.bar(edges,np.ones((len(edges),)), color=colors, 
                    width=(edges[-1]-edges[0])/len(edges) , edgecolor=colors)
    ax2.set_ylim([0,1])
    ax2.set_xlim(-1, 1)
    ax2.tick_params(axis='y', which='both', left='off', right='off', 
                    labelleft='off') 
    ax2.tick_params(axis='x', which='both', bottom='off', top='off', 
                    labelbottom='on')
    ax2.set_xlabel('Motif Distance (kb)')


    counts,edges = np.histogram(q3_distances, bins=bins)
    edges        = (edges[1:]+edges[:-1])/2. 
    # plt.bar(edges, counts, width=(edges[-1]-edges[0])/bins)

    norm    = matplotlib.colors.Normalize(vmin=min(counts), vmax=max(counts))
    cmap    = cm.YlOrRd
    m       = cm.ScalarMappable(norm=norm, cmap=cmap)
    colors  = [m.to_rgba(c) for c in counts] 
    
    ax3 = F.add_subplot(143)
    ax3.bar(edges,np.ones((len(edges),)), color=colors, 
                    width=(edges[-1]-edges[0])/len(edges) , edgecolor=colors)
    ax3.set_ylim([0,1])
    ax3.set_xlim(-1, 1)
    ax3.tick_params(axis='y', which='both', left='off', right='off', 
                    labelleft='off') 
    ax3.tick_params(axis='x', which='both', bottom='off', top='off', 
                    labelbottom='on')
    ax3.set_xlabel('Motif Distance (kb)')

    counts,edges = np.histogram(q4_distances, bins=bins)
    edges        = (edges[1:]+edges[:-1])/2. 
    # plt.bar(edges, counts, width=(edges[-1]-edges[0])/bins)

    norm    = matplotlib.colors.Normalize(vmin=min(counts), vmax=max(counts))
    cmap    = cm.YlOrRd
    m       = cm.ScalarMappable(norm=norm, cmap=cmap)
    colors  = [m.to_rgba(c) for c in counts] 
    
    ax4 = F.add_subplot(144)
    ax4.bar(edges,np.ones((len(edges),)), color=colors, 
                    width=(edges[-1]-edges[0])/len(edges) , edgecolor=colors)
    ax4.set_ylim([0,1])
    ax4.set_xlim(-1, 1)
    ax4.tick_params(axis='y', which='both', left='off', right='off', 
                    labelleft='off') 
    ax4.tick_params(axis='x', which='both', bottom='off', top='off', 
                    labelbottom='on')
    ax4.set_xlabel('Motif Distance (kb)')


    plt.savefig(os.path.join(figuredir, motif_file 
                                        + '_distance_distribution.png'), 
                dpi=None, bbox_inches='tight')
    plt.close()
#==============================================================================

#==============================================================================
def simulation_plot(figuredir=None, simES=None, actualES=None,
                        motif_file=None, dpi=None):
    '''This function plots the simulated 'enrichment' scores against the
        observed 'enrichment' score

    Parameters
    ----------
    figuredir : string
        the full path to the figuredir within the ouptut directory containing
        all figures and plots

    simES : list
        a list of simulated 'enrichment' scores calculated by randomizing
        region rank

    actualES : float
        the observed 'enchrichment' score. Can be calculated in many different
        ways..

    motif_file : string
        the name of the motif thats associated with all the input data. Used 
        for figure naming purposes.

    Returns
    -------
    None
    '''
    F = plt.figure(figsize=(7,6))
    ax2 = plt.subplot(111)
    maximum = max(simES)
    minimum = min(simES)
    ax2.hist(simES,bins=100)
    width = (maximum-minimum)/100.0
    rect = ax2.bar(actualES,ax2.get_ylim()[1],color='red',width=width*2)[0]
    height = rect.get_height()
    ax2.text(rect.get_x() + rect.get_width()/2., 1.05*height, 
                'Observed ES', ha='center', va='bottom')

    ax2.set_xlim([min(minimum,actualES)-(width*40), \
                max(maximum,actualES)+(width*40)])

    ax2.set_ylim([0,(1.05*height)+5])
    ax2.tick_params(axis='y', which='both', left='off', right='off', 
                    labelleft='on')

    ax2.tick_params(axis='x', which='both', bottom='off', top='off', 
                    labelbottom='on')

    plt.title('Distribution of Simulated Enrichment Scores',fontsize=14)
    ax2.set_ylabel('Number of Simulations',fontsize=14)
    ax2.set_xlabel('Area Under the Curve (AUC)',fontsize=14)
    plt.savefig(os.path.join(figuredir, motif_file 
                + '_simulation_plot.png'),dpi=dpi, bbox_inches='tight')

    plt.cla()
#==============================================================================

#==============================================================================
def moustache_plot(figuredir=None, ESlist=None, PADJlist=None, sigx=None, 
                    sigy=None, dpi=None):

    '''This function plots a moustache plot for all motifs. In the x-axis, 
        all observed 'enrichment' scores are plotted against the adjusted
        pvalue of each motif

    Parameters
    ----------
    figuredir : string
        the full path to the figuredir within the ouptut directory containing
        all figures and plots

    ESlist : list or array
        a list of 'enrichment' scores to be plotted on the x-axis

    PADJlist : list or array
        a list of p-adjusted values to be plotted on the y-axis

    sigx : list or array
        a list of significant motifs to be colored red

    sigy : list or array
        a list of significant motifs to be colored red

    Returns
    -------
    None
    '''
    import config
    import math
    dpi = config.DPI
    PADJlist = [-math.log(x,10) if x > 0 else 0 for x in PADJlist]
    sigy = [-math.log(x,10) if x > 0 else 0 for x in sigy]
    max_val = max(PADJlist)
    PADJlist = [x if x != 0 else max_val for x in PADJlist]
    sigy = [x if x != 0 else max_val for x in sigy]
    plt.figure(figsize=(7,6))
    ax = plt.subplot(111)
    ax.scatter(ESlist,PADJlist,color='black',edgecolor='')
    ax.scatter(sigx,sigy,color='red',edgecolor='')
    ax.set_title("TFEA Moustache Plot",fontsize=14)
    ax.set_xlabel("Area Under the Curve (AUC)",fontsize=14)
    ax.set_ylabel("P-value (PADJ)",fontsize=14)
    xlimit = math.fabs(max(ESlist,key=abs))
    ylimit = math.fabs(max(PADJlist,key=abs))
    ax.set_xlim([-xlimit,xlimit])
    ax.set_ylim([0,ylimit])
    ax.tick_params(axis='y', which='both', left='off', right='off', 
                    labelleft='on')

    ax.tick_params(axis='x', which='both', bottom='off', top='off', 
                    labelbottom='on')

    plt.savefig(os.path.join(figuredir, 'TFEA_Results_Moustache_Plot.png'), 
                    dpi=dpi, bbox_inches='tight')
    plt.cla()
#==============================================================================

#==============================================================================
def pval_histogram_plot(figuredir=None, PVALlist=None):
    '''This function plots a histogram of p-values for all motifs

    Parameters
    ----------
    figuredir : string
        the full path to the figuredir within the ouptut directory containing
        all figures and plots

    PVALlist : list or array
        a list of p-values corresponding to the observed 'enrichment' score
        compared to the distribution of simulated 'enrichment' scores.
    
    Returns
    -------
    None
    '''
    plt.figure(figsize=(7,6))
    ax = plt.subplot(111)
    binwidth = 1.0/100.0
    # print PVALlist
    ax.hist(PVALlist,bins=np.arange(0,0.5+binwidth,binwidth),color='green')
    ax.set_title("TFEA P-value Histogram",fontsize=14)
    ax.set_xlabel("P-value",fontsize=14)
    ax.set_ylabel("Count",fontsize=14)
    ax.tick_params(axis='y', which='both', left='off', right='off', 
                    labelleft='on')

    ax.tick_params(axis='x', which='both', bottom='off', top='off', 
                    labelbottom='on')

    plt.savefig(os.path.join(figuredir, 'TFEA_Pval_Histogram.png'),
                    bbox_inches='tight')
    plt.cla()
#==============================================================================

#==============================================================================
def MA_plot(figuredir=None, label1=None, label2=None, POSlist=None, 
            ESlist=None, MAx=None, MAy=None):
    '''This function plots an 'MA' plot with the 'enrichment' score on the 
        y-axis and the number of hits within the largewindow in the x-axis

    Parameters
    ----------
    figuredir : string
        the full path to the figuredir within the ouptut directory containing
        all figures and plots

    label1 : string
        a label that corresponds to condition1

    label2 : string
        a label that corresponds to condition2

    POSlist : list or array
        a list of 'positive' hits for each motif defined as being within a 
        largewindow

    ESlist : list or array
        a list of 'enrichment' scores for each motif

    MAx : list or array
        a list of x-values corresponding to significant motifs to be colored 
        red

    MAx : list or array
        a list of y-values corresponding to significant motifs to be colored 
        red

    Returns
    -------
    None
    '''
    plt.figure(figsize=(7,6))
    ax = plt.subplot(111)
    ax.scatter(POSlist,ESlist,color='black',edgecolor='')
    ax.scatter(MAx,MAy,color='red',edgecolor='')
    ax.set_title("TFEA MA-Plot",fontsize=14)
    ax.set_ylabel("Area Under the Curve (AUC)", fontsize=14)

    ax.set_xlabel("Motif Hits Log10",fontsize=14)
    ax.tick_params(axis='y', which='both', left='off', right='off', 
                    labelleft='on')

    ax.tick_params(axis='x', which='both', bottom='off', top='off', 
                    labelbottom='on')

    plt.savefig(os.path.join(figuredir, 'TFEA_NES_MA_Plot.png'),
                bbox_inches='tight')
    plt.cla()
#==============================================================================

#==============================================================================
def meta_eRNA_quartiles(figuredir=None, label1=None, label2=None, 
                        q1label1pos=None, q1label1neg=None, q1label2pos=None,
                        q1label2neg=None, q2label1pos=None, q2label1neg=None,
                        q2label2pos=None, q2label2neg=None, q3label1pos=None, 
                        q3label1neg=None, q3label2pos=None, q3label2neg=None, 
                        q4label1pos=None, q4label1neg=None, q4label2pos=None, 
                        q4label2neg=None, largewindow=None, dpi=None):
    '''This function creates a plot with meta profiles for inputted regions
        of interest. It creates a separate plot for quartiles 1, 2/3, and 4.
        Additionally, under each meta plot it also produces a heatmap histogram
        of motif hits.

    Parameters
    ----------
    figuredir : string
        the full path to the figuredir within the ouptut directory containing
        all figures and plots

    label1 : string
        a label that corresponds to condition1

    label2 : string
        a label that corresponds to condition2

    POSlist : list or array
        a list of 'positive' hits for each motif defined as being within a 
        largewindow

    ESlist : list or array
        a list of 'enrichment' scores for each motif

    MAx : list or array
        a list of x-values corresponding to significant motifs to be colored 
        red

    MAx : list or array
        a list of y-values corresponding to significant motifs to be colored 
        red

    Returns
    -------
    None
    '''
    F = plt.figure(figsize=(15.5,3))
    
    xvals = range(-int(largewindow),int(largewindow))
    ylim = [min(q1label1neg+q1label2neg+q2label1neg+q2label2neg+q3label1neg
                +q3label2neg+q4label1neg+q4label2neg),
            max(q1label1pos+q1label2pos+q2label1pos+q2label2pos+q3label1pos
            +q3label2pos+q4label1pos+q4label2pos)]

    # First quartile plot
    ax0 = plt.subplot(141)
    line1, = ax0.plot(xvals,q1label1pos,color='blue',label=label1)
    ax0.plot(xvals,q1label1neg,color='blue')
    line2, = ax0.plot(xvals,q1label2pos,color='red',label=label2)
    ax0.plot(xvals,q1label2neg,color='red')
    ax0.legend(loc=2,fontsize='small')
    ax0.set_title('Q1',fontsize=14)
    ax0.tick_params(axis='y', which='both', left='off', right='off', labelleft='on')
    ax0.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
    ax0.set_ylabel('Reads per Millions Mapped',fontsize=14)
    ax0.set_xlabel('Distance to eRNA Origin (bp)')
    ax0.set_ylim(ylim)


    # Second quartile plot
    ax2 = plt.subplot(142)
    line1, = ax2.plot(xvals,q2label1pos,color='blue',label=label1)
    ax2.plot(xvals,q2label1neg,color='blue')
    line2, = ax2.plot(xvals,q2label2pos,color='red',label=label2)
    ax2.plot(xvals,q2label2neg,color='red')
    # ax2.legend(loc=1)
    ax2.set_title('Q2',fontsize=14)
    ax2.tick_params(axis='y', which='both', left='off', right='off', labelleft='off')
    ax2.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
    ax2.set_xlabel('Distance to eRNA Origin (bp)')
    ax2.set_ylim(ylim)

    # Second quartile plot
    ax3 = plt.subplot(143)
    line1, = ax3.plot(xvals,q3label1pos,color='blue',label=label1)
    ax3.plot(xvals,q3label1neg,color='blue')
    line2, = ax3.plot(xvals,q3label2pos,color='red',label=label2)
    ax3.plot(xvals,q3label2neg,color='red')
    # ax3.legend(loc=1)
    ax3.set_title('Q3',fontsize=14)
    ax3.tick_params(axis='y', which='both', left='off', right='off', labelleft='off')
    ax3.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
    ax3.set_xlabel('Distance to eRNA Origin (bp)')
    ax3.set_ylim(ylim)

    # Third quartile plot
    ax4 = plt.subplot(144)
    line1, = ax4.plot(xvals,q4label1pos,color='blue',label=label1)
    ax4.plot(xvals,q4label1neg,color='blue')
    line2, = ax4.plot(xvals,q4label2pos,color='red',label=label2)
    ax4.plot(xvals,q4label2neg,color='red')
    # ax4.legend(loc=1)
    ax4.set_title('Q4',fontsize=14)
    ax4.tick_params(axis='y', which='both', left='off', right='off', labelleft='off')
    ax4.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
    ax4.set_xlabel('Distance to eRNA Origin (bp)')
    ax4.set_ylim(ylim)

    plt.savefig(os.path.join(figuredir, 'meta_plot.png'), dpi=dpi, 
        bbox_inches='tight')
    plt.cla()
#==============================================================================

#==============================================================================
def create_text_output(outputdir=None, TFresults=None):
    '''Creates a .txt output of results from TFEA

    Parameters
    ----------
    outputdir : string
        the full path to TFEA output directory (created by TFEA)

    TFresults : list or array
        a list of lists contining 'enrichment' scores, normalized 'enrichment' 
        scores, p-value, p-adj, and number of hits for each individual motif

    Returns
    -------
    None
    '''
    TFresults = sorted(TFresults, key=lambda x: x[3])
    outfile = open(os.path.join(outputdir, 'results.txt'), 'w')
    outfile.write('TF-Motif\tES\tNES\tP-value\tPADJ\tHITS\n')
    for val in TFresults:
        outfile.write('\t'.join([str(val[i]) for i in range(len(val))]) + '\n')
    outfile.close()
#==============================================================================

#==============================================================================
def create_summary_html():
    '''Creates the main html output and also individual html outputs for each
        motif
    
    Parameters
    ----------
    outputdir : string
        the full path to the output directory created by TFEA

    beds : list or array
        a list of full paths to bed files to be considered as regions of 
        interest

    label1 : string
        an informative label describing sample corresponding to condition1

    label2 : string
        an informative label describing sample corresponding to condition2

    bam1 : list or array
        a list of full paths to bam files corresponding to condition1

    bam2 : list or array
        a list of full paths to bam files corresponding to condition2

    singlemotif : boolean or string
        either False if all motifs should be considered in TFEA or the name of
        a specific motif to be analyzed

    motif_hits : string
        the full path to a directory containing motif hits across the genome

    output : string
        the full path to a user-specified output directory. TFEA will create
        a new folder within this directory - this is called outputdir

    padj_cutoff : float
        the cutoff value for determining significance

    plot : boolean
        a switch that controls whether all motifs are plotted or just 
        significant ones defined by the p-adj cutoff

    combine : boolean
        a switch that determines whether bed files within the beds variable
        get combined and merged using bedtools

    count : boolean
        a switch that controls whether reads are counted over the regions of
        interest

    deseq : boolean
        a switch that controls whether DE-Seq is performed on the inputted
        regions that have been counted over

    calculate : boolean
        a switch that determines whether the TFEA calculation is performed

    TFresults : list or array
        a list of lists contining 'enrichment' scores, normalized 'enrichment' 
        scores, p-value, p-adj, and number of hits for each individual motif

    COMBINEtime : float
        the time it took to combine and merge the bed files using bedtools

    COUNTtime : float
        the time it took to count reads over regions of interest

    DESEQtime : float
        the time it took to perform DE-Seq using the counts file

    CALCULATEtime : float
        the time it took to perform TFEA

    Returns
    -------
    None
    '''
    #Using a config file
    import config
    outputdir = config.OUTPUTDIR
    bed1 = config.BED1
    bed2 = config.BED2
    label1 = config.LABEL1
    label2 = config.LABEL2
    bam1 = config.BAM1
    bam2 = config.BAM2
    singlemotif = config.SINGLEMOTIF
    motif_hits = config.MOTIF_GENOMEWIDE_HITS
    output = config.OUTPUT
    padj_cutoff = config.PADJCUTOFF
    smallwindow = config.SMALLWINDOW
    largewindow = config.LARGEWINDOW
    plot = config.PLOTALL
    combine = config.COMBINE
    count = config.COUNT
    deseq = config.DESEQ
    calculate = config.CALCULATE
    fimo = config.FIMO
    homer = config.HOMER
    temp = config.TEMP
    logos = config.LOGOS
    motifdatabase = config.MOTIFDATABASE
    genomefasta = config.GENOMEFASTA

    #summary.html contains all user-defined variables, and also information 
    #about module used
    outfile = open(os.path.join(outputdir,'summary.html'),'w')
    outfile.write("""<!DOCTYPE html>
            <html>
            <head>
            <title>Variables Used</title>
            </head>
            <body>
                <h1>Variables Used</h1>
                <p>OUTPUT = """+output+"""
                <p>BED1 = """+str(bed1)+"""</p>
                <p>BED2 = """+str(bed2)+"""</p>
                <p>LABEL1 = """+label1+"""</p>
                <p>LABEL2 = """+label2+"""</p>
                <p>BAM1 = """+str(bam1)+"""</p>
                <p>BAM2 = """+str(bam2)+"""</p>
                <p>COMBINE = """+str(combine)+"""</p>
                <p>COUNT = """+str(count)+"""</p>
                <p>DESEQ = """+str(deseq)+"""</p>
                <p>CALCULATE = """+str(calculate)+"""</p>
                <p>SINGLEMOTIF = """+str(singlemotif)+"""</p>
                <p>PLOT = """+str(plot)+"""</p>
                <p>FIMO = """+str(fimo)+"""</p>
                <p>HOMER = """+str(homer)+"""</p>
                <p>TEMP = """+str(temp)+"""</p>
                <p>PADJCUTOFF = """+str(padj_cutoff)+"""</p>
                <p>SMALLWINDOW = """+str(smallwindow)+"""</p>
                <p>LARGEWINDOW = """+str(largewindow)+"""</p>
                <p>MOTIF_HITS = """+str(motif_hits)+"""</p>
                <p>LOGOS = """+str(logos)+"""</p>
                <p>MOTIFDATABASE = """+str(motifdatabase)+"""</p>
                <p>GENOMEFASTA = """+str(genomefasta)+"""</p>
            </body>""")
    outfile.close()
#==============================================================================

#==============================================================================
def create_motif_result_html(TFresults=None):
    '''Creates the main html output and also individual html outputs for each
        motif
    
    Parameters
    ----------
    outputdir : string
        the full path to the output directory created by TFEA

    beds : list or array
        a list of full paths to bed files to be considered as regions of 
        interest

    label1 : string
        an informative label describing sample corresponding to condition1

    label2 : string
        an informative label describing sample corresponding to condition2

    bam1 : list or array
        a list of full paths to bam files corresponding to condition1

    bam2 : list or array
        a list of full paths to bam files corresponding to condition2

    singlemotif : boolean or string
        either False if all motifs should be considered in TFEA or the name of
        a specific motif to be analyzed

    motif_hits : string
        the full path to a directory containing motif hits across the genome

    output : string
        the full path to a user-specified output directory. TFEA will create
        a new folder within this directory - this is called outputdir

    padj_cutoff : float
        the cutoff value for determining significance

    plot : boolean
        a switch that controls whether all motifs are plotted or just 
        significant ones defined by the p-adj cutoff

    combine : boolean
        a switch that determines whether bed files within the beds variable
        get combined and merged using bedtools

    count : boolean
        a switch that controls whether reads are counted over the regions of
        interest

    deseq : boolean
        a switch that controls whether DE-Seq is performed on the inputted
        regions that have been counted over

    calculate : boolean
        a switch that determines whether the TFEA calculation is performed

    TFresults : list or array
        a list of lists contining 'enrichment' scores, normalized 'enrichment' 
        scores, p-value, p-adj, and number of hits for each individual motif

    COMBINEtime : float
        the time it took to combine and merge the bed files using bedtools

    COUNTtime : float
        the time it took to count reads over regions of interest

    DESEQtime : float
        the time it took to perform DE-Seq using the counts file

    CALCULATEtime : float
        the time it took to perform TFEA

    Returns
    -------
    None
    '''
    #Using a config file
    import config
    outputdir = config.OUTPUTDIR
    padj_cutoff = config.PADJCUTOFF
    singlemotif = config.SINGLEMOTIF
    plot = config.PLOTALL

    #For each TF motif with an PADJ value less than a cutoff, an html file is 
    #created to be used in results.html
    positivelist = [x[0] for x in TFresults 
                    if x[2] > 0 and (plot or x[-1] < padj_cutoff)]
    negativelist = [x[0] for x in TFresults 
                    if x[2] < 0 and (plot or x[-1] < padj_cutoff)]

    for i in range(len(TFresults)):
        MOTIF_FILE,ES,NES,PVAL,POS,PADJ = TFresults[i] 
        if plot or PADJ < padj_cutoff or singlemotif:
            if ES > 0:
                try:
                    NEXT_MOTIF = positivelist[positivelist.index(MOTIF_FILE)+1]
                except IndexError:
                    NEXT_MOTIF = positivelist[0]
                try:
                    PREV_MOTIF = positivelist[positivelist.index(MOTIF_FILE)-1]
                except IndexError:
                    PREV_MOTIF = positivelist[len(positivelist)]
            else:
                try:
                    NEXT_MOTIF = negativelist[negativelist.index(MOTIF_FILE)+1]
                except IndexError:
                    NEXT_MOTIF = negativelist[0]
                try:
                    PREV_MOTIF = negativelist[negativelist.index(MOTIF_FILE)-1]
                except IndexError:
                    PREV_MOTIF = negativelist[len(negativelist)]
            direct_logo = MOTIF_FILE.strip('HO_') + "_direct.png"
            reverse_logo = MOTIF_FILE.strip('HO_') + "_revcomp.png"
            outfile = open(os.path.join(outputdir, 'plots', MOTIF_FILE 
                            + '.results.html'),'w')
            outfile.write("""<!DOCTYPE html>
    <html>
    <head>
    <title>"""+MOTIF_FILE+""" Results</title>
    <style>
        table {
            font-family: arial, sans-serif;
            border-collapse: collapse;
            width: 100%;
        }

        .row {
        display: flex; /* equal height of the children */
        width: 100%;
        padding-bottom: 50px
        }

        img {
            max-width: 100%;
            max-height: 100%;
        }

        td, th {
            border: 1px solid #dddddd;
            text-align: left;
            padding: 8px;
        }

        tr:nth-child(even) {
            background-color: #dddddd;
        }
    </style>
    </head>
    <body style="width:1300px; margin:0 auto;">
        <div>
            <div style="float:left">
                <a href="./"""+PREV_MOTIF+""".results.html">PREV</a>
            </div>
            <div style="float:right">
                <a href="./"""+NEXT_MOTIF+""".results.html">NEXT</a>
            </div>
            <div style="text-align:center">
                <a href="../results.html">ALL</a>
        </div>
        <div class="row">
        </div>
            <h1>"""+MOTIF_FILE+""" Results</h1>
        <div>
            <div style="float: middle; width: 1300px; padding-bottom:25px; \
                padding-top:25px">
                <table> 
                    <tr>
                        <th>TF Motif</th>
                        <th>AUC</th> 
                        <th>P-value</th>
                        <th>PADJ</th>
                        <th>HITS</th>
                    </tr>
                    <tr>
                        <td>"""+MOTIF_FILE+"""</td>
                        <td>"""+str("%.2f" % ES)+"""</td>
                        <td>"""+str("%.4g" % PVAL)+"""</td>
                        <td>"""+str("%.4g" % PADJ)+"""</td>
                        <td>"""+str(POS)+"""</td>
                    </tr>
                </table>
            </div>
        </div>
        <div>
            <div style="float: left; width 1250px; padding-bottom:50px; \
                padding-top:50px">
                <img src="./"""+MOTIF_FILE+"""_enrichment_plot.png" \
                    alt="Enrichment Plot">
            </div>
        </div>
        <div class="row">
            <div style="float: right; width: 600px">
                <p>Forward:</p>
                <img src="./"""+direct_logo+"""" \
                    alt="Forward Logo">
                <p></p>
                <p>Reverse:</p>
                <img src="./"""+reverse_logo+"""" \
                    alt="Reverse Logo">
            </div>
            <div style="float:left; width: 600px">
                <img src="./"""+MOTIF_FILE+"""_simulation_plot.png" \
                    alt="Simulation Plot">
            </div>
        </div>

    </body>
    </html>""")
            outfile.close()
            PREV_MOTIF = MOTIF_FILE
#==============================================================================

#==============================================================================
def create_main_results_html(TFresults=None, COMBINEtime=None, COUNTtime=None, 
                                DESEQtime=None, CALCULATEtime=None):
    '''Creates the main html output and also individual html outputs for each
        motif
    
    Parameters
    ----------
    outputdir : string
        the full path to the output directory created by TFEA

    beds : list or array
        a list of full paths to bed files to be considered as regions of 
        interest

    label1 : string
        an informative label describing sample corresponding to condition1

    label2 : string
        an informative label describing sample corresponding to condition2

    bam1 : list or array
        a list of full paths to bam files corresponding to condition1

    bam2 : list or array
        a list of full paths to bam files corresponding to condition2

    singlemotif : boolean or string
        either False if all motifs should be considered in TFEA or the name of
        a specific motif to be analyzed

    motif_hits : string
        the full path to a directory containing motif hits across the genome

    output : string
        the full path to a user-specified output directory. TFEA will create
        a new folder within this directory - this is called outputdir

    padj_cutoff : float
        the cutoff value for determining significance

    plot : boolean
        a switch that controls whether all motifs are plotted or just 
        significant ones defined by the p-adj cutoff

    combine : boolean
        a switch that determines whether bed files within the beds variable
        get combined and merged using bedtools

    count : boolean
        a switch that controls whether reads are counted over the regions of
        interest

    deseq : boolean
        a switch that controls whether DE-Seq is performed on the inputted
        regions that have been counted over

    calculate : boolean
        a switch that determines whether the TFEA calculation is performed

    TFresults : list or array
        a list of lists contining 'enrichment' scores, normalized 'enrichment' 
        scores, p-value, p-adj, and number of hits for each individual motif

    COMBINEtime : float
        the time it took to combine and merge the bed files using bedtools

    COUNTtime : float
        the time it took to count reads over regions of interest

    DESEQtime : float
        the time it took to perform DE-Seq using the counts file

    CALCULATEtime : float
        the time it took to perform TFEA

    Returns
    -------
    None
    '''
    #Using a config file
    import config
    outputdir = config.OUTPUTDIR
    label1 = config.LABEL1
    label2 = config.LABEL2
    padj_cutoff = config.PADJCUTOFF
    plot = config.PLOTALL
    combine = config.COMBINE
    count = config.COUNT
    deseq = config.DESEQ
    calculate = config.CALCULATE

    outfile = open(os.path.join(outputdir, 'results.html'),'w')
    outfile.write("""<!DOCTYPE html>
    <html>
    <head>
    <title>TFEA Results """ + label1 + """ vs. """ + label2 +"""</title>
    <style>
        table {
            font-family: arial, sans-serif;
            border-collapse: collapse;
            width: 100%;
        }

        .row {
        display: flex; /* equal height of the children */
        width: 100%;
        padding-bottom: 50px
        }

        img {
            max-width: 100%;
            max-height: 100%;
        }

        td, th {
            border: 1px solid #dddddd;
            text-align: left;
            padding: 8px;
        }

        tr:nth-child(even) {
            background-color: #dddddd;
        }
    </style>
    </head>
    <body style="width:1300px; margin:0 auto;">

        <h1>TFEA Results """ + label1 + """ vs. """ + label2 + """</h1>
        <div class="row">
            <div style="float: left; width: 45%">
                <img src="./plots/TFEA_NES_MA_Plot.png" alt="NES MA-Plot">
            </div>
            <div style="float: left; width: 45%">
                <img src="./plots/DESEQ_MA_Plot.png" alt="DE-Seq MA-Plot">
            </div>
        </div>
        <div class="row">
            <div style="float: left; width:45%">
                <img src="./plots/TFEA_Pval_Histogram.png" alt="TFEA P-Value \
                    Histogram">
            </div>
            <div id="Summary of Variables Used" style="float: right; \
                width: 45%">
                <p><a href="./Summary.html">Full Summary of Variables Used\
                </a></p>
                <p><b>PADJ < """ + str(padj_cutoff) + """</b></p>
                <table>
                    <tr>
                        <th>Module</th>
                        <th>Switch</th>
                        <th>Time (hh:mm:ss)</th>
                    </tr>
                    <tr>
                        <td>COMBINE</td>
                        <td>"""+str(combine)+"""</td>
                        <td>"""+str(datetime.timedelta(
                                    seconds=int(COMBINEtime)))
                        +"""</td>
                    </tr>
                    <tr>
                        <td>COUNT</td>
                        <td>"""+str(count)+"""</td>
                        <td>"""+str(datetime.timedelta(seconds=int(COUNTtime)))
                        +"""</td>
                    </tr>
                    <tr>
                        <td>DESEQ</td>
                        <td>"""+str(deseq)+"""</td>
                        <td>"""+str(datetime.timedelta(seconds=int(DESEQtime)))
                        +"""</td>
                    </tr>
                    <tr>
                        <td>CALCULATE</td>
                        <td>"""+str(calculate)+"""</td>
                        <td>"""+str(datetime.timedelta(
                                    seconds=int(CALCULATEtime)))
                        +"""</td>
                    </tr>
                    <tr>
                        <td><b>TOTAL</b></td>
                        <td> </td>
                        <td>"""+str(datetime.timedelta(
                            seconds=int(COMBINEtime)
                                    +int(COUNTtime)
                                    +int(DESEQtime)
                                    +int(CALCULATEtime)))
                            +"""</td>
                    </tr>
                </table>   
            </div>
        </div>
        <div>
            <div id="Positive Enrichment Score" style="float: left; width:45%">
                <h1>Positive Enrichment Score</h1>
                <table> 
                    <tr>
                        <th>TF Motif</th>
                        <th>AUC</th>
                        <th>P-value</th>
                        <th>PADJ</th>
                        <th>HITS</th>
                    </tr>
                """)

    for MOTIF_FILE,ES,NES,PVAL,POS,PADJ in TFresults:
        if NES > 0:
            if PADJ < padj_cutoff:
                outfile.write("""
            <tr style="color: red;">
                <td><a href="./plots/"""+MOTIF_FILE+""".results.html">"""
                    +MOTIF_FILE+"""</td>
                <td>"""+str("%.4f" % ES)+"""</td>
                <td>"""+str("%.3g" % PVAL)+"""</td>
                <td>"""+str("%.3g" % PADJ)+"""</td>
                <td>"""+str(POS)+"""</td>
            </tr>
                    """)
            elif plot:
                outfile.write("""
            <tr>
                <td><a href="./plots/"""+MOTIF_FILE+""".results.html">"""
                    +MOTIF_FILE+"""</td>
                <td>"""+str("%.4f" % ES)+"""</td>
                <td>"""+str("%.3g" % PVAL)+"""</td>
                <td>"""+str("%.3g" % PADJ)+"""</td>
                <td>"""+str(POS)+"""</td>
            </tr>
                    """)

            else:
                outfile.write("""
            <tr>
                <td>"""+MOTIF_FILE+"""</td>
                <td>"""+str("%.4f" % ES)+"""</td>
                <td>"""+str("%.3g" % PVAL)+"""</td>
                <td>"""+str("%.3g" % PADJ)+"""</td>
                <td>"""+str(POS)+"""</td>
            </tr>
                    """)


    outfile.write("""            
        </table>
    </div>

    <div id="Negative Enrichment Score" style="float: right; width: 45%">
        <h1>Negative Enrichment Score</h1>
        <table> 
            <tr>
                <th>TF Motif</th>
                <th>AUC</th> 
                <th>P-value</th>
                <th>PADJ</th>
                <th>HITS</th>
            </tr>
                """)

    for MOTIF_FILE,ES,NES,PVAL,POS,PADJ in TFresults:
        if NES < 0:
            if PADJ < padj_cutoff:
                outfile.write("""
            <tr style="color: red;">
                <td><a href="./plots/"""+MOTIF_FILE+""".results.html">"""
                    +MOTIF_FILE+"""</td>
                <td>"""+str("%.4f" % ES)+"""</td>
                <td>"""+str("%.3g" % PVAL)+"""</td>
                <td>"""+str("%.3g" % PADJ)+"""</td>
                <td>"""+str(POS)+"""</td>
            </tr>
                    """)
            elif plot:
                outfile.write("""
            <tr>
                <td><a href="./plots/"""+MOTIF_FILE+""".results.html">"""
                    +MOTIF_FILE+"""</td>
                <td>"""+str("%.4f" % ES)+"""</td>
                <td>"""+str("%.3g" % PVAL)+"""</td>
                <td>"""+str("%.3g" % PADJ)+"""</td>
                <td>"""+str(POS)+"""</td>
            </tr>
                    """)
            else:
                outfile.write("""
            <tr>
                <td>"""+MOTIF_FILE+"""</td>
                <td>"""+str("%.4f" % ES)+"""</td>
                <td>"""+str("%.3g" % PVAL)+"""</td>
                <td>"""+str("%.3g" % PADJ)+"""</td>
                <td>"""+str(POS)+"""</td>
            </tr>
                    """)

    outfile.write("""        
            </table>
        </div>
        </div>

    </body>
    </html>""")

    outfile.close()
#==============================================================================

#==============================================================================
def get_TFresults_from_txt(results_file=None):
    '''This function parses a results_file into a TFresults array

    Parameters
    ----------
    results_file : string
        file containing TFEA results
        
    Returns
    -------
    TFresults : array
        results for all motifs within the results_file
                 
    Raises
    ------
    None
    '''
    TFresults = list()
    with open(results_file) as F:
        F.readline()
        for line in F:
            line = line.strip('\n').split('\t')
            TFresults.append(line)
    return TFresults
#==============================================================================

#==============================================================================
if __name__ == "__main__":
    # results_file = "/Users/jonathanrubin/Google_Drive/Colorado University/Jonathan/TFEA_outputs/IRIS/new_outputs/TFEA_30_IFN-30_IFN_CA_0/results.txt"
    # TFresults = get_TFresults_from_txt(results_file=results_file)
    # ESlist = [float(x[1]) for x in TFresults]
    # PADJlist = [float(x[-3]) for x in TFresults]
    # padj_cutoff = 0.001
    # sigx = [x for x, p in zip(ESlist, PADJlist) if p < padj_cutoff]
    # sigy = [p for x, p in zip(ESlist, PADJlist) if p < padj_cutoff]
    # figuredir = "/Users/jonathanrubin/Google_Drive/Colorado University/Jonathan/TFEA_outputs/IRIS/new_outputs/TFEA_30_IFN-30_IFN_CA_0/plots/"
    # import config
    # config.DPI = 1200
    # moustache_plot(figuredir=figuredir, ESlist=ESlist, PADJlist=PADJlist, 
    #                 sigx=sigx, sigy=sigy)

    #==================METAeRNA TEST===================
    import matplotlib.pyplot as plt
    posprofile1, negprofile1, posprofile2, negprofile2 = meta_profile(regionlist=[('chr17', '67603195','67603800')], 
                                                            millions_mapped=[54307224,158978147], 
                                                            largewindow=1500.0, 
                                                            bam1=['../../TFEA_outputs/SRR1105736.fastqbowtie2.sorted.bam'], 
                                                            bam2=['../../TFEA_outputs/SRR1105738.fastqbowtie2.sorted.bam'])
    print sum(posprofile1), sum(negprofile1), sum(posprofile2), sum(negprofile2)

    F = plt.figure(figsize=(15.5,6))
    ax0 = plt.subplot(111)
    xvals = range(-int(1500),int(1500))
    line1, = ax0.plot(xvals,posprofile1,color='blue',label='DMSO')
    ax0.plot(xvals,negprofile1,color='blue')
    line2, = ax0.plot(xvals,posprofile2,color='red',label='Nutlin')
    ax0.plot(xvals,negprofile2,color='red')
    ax0.legend(loc=1)
    ax0.set_title('Meta Plot over all regions',fontsize=14)
    ax0.tick_params(axis='y', which='both', left='off', right='off', labelleft='on')
    ax0.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
    ax0.set_ylabel('Normalized Read Coverage',fontsize=14)
    ax0.set_xlabel('Distance to eRNA Origin (bp)')
    plt.savefig('./test.png',bbox_inches='tight')
    plt.show()
    # plt.cla()

    #==================Enrichment Plot TEST=================
    # import math
    # sorted_distances = np.arange(0,150,10000)
    # cumscore = np.cumsum(sorted_distances)
    # logpval = np.arange(0.0001,1,10000)
    # logpval = [math.exp(-x) for x in logpval]
    # enrichment_plot(largewindow=150.0, smallwindow=15.0, figuredir='../',
    #                 cumscore=cumscore, sorted_distances=sorted_distances, logpval=logpval, 
    #                 updistancehist=None, downdistancehist=None, 
    #                 gc_array=None, dpi=None, save=False)