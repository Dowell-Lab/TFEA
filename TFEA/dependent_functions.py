#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''This file contains a list of dependent functions that rely on other 
    functions to work properly
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
import sys
import math
import numpy as np
from multiprocessing import Pool
import multiprocessing as mp
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.cm as cm
from matplotlib import gridspec
from scipy.stats import norm 
import time
import config
import independent_functions
#==============================================================================
#Functions
#==============================================================================
def deseq_run(bam1=list(), bam2=list(), tempdir=str(), count_file=str(),
                label1=str(), label2=str(), figuredir=str()):
    '''Writes the DE-Seq script, runs the DE-Seq script, and plots the DE-Seq 
        MA plot

    Parameters
    ----------
    tempdir : string
        full path to temp directory in output directory (created by TFEA)

    Returns
    ----------
    deseq_file : string
        full path to the outputted deseq file
    '''
    independent_functions.write_deseq_script(bam1=bam1, bam2=bam2, 
                                                tempdir=tempdir, 
                                                count_file=count_file,
                                                label1=label1, label2=label2)
    os.system("R < " + os.path.join(tempdir, "DESeq.R") + " --no-save")
    deseq_file = os.path.join(tempdir, 'DESeq.res.txt')
    independent_functions.plot_deseq_MA(deseq_file=deseq_file,
                                            label1=label1, label2=label2, 
                                            figuredir=figuredir)
    
    return deseq_file
#==============================================================================

#==============================================================================
def calculate_es_auc(args):
    '''This function calculates the AUC for any TF based on the TF motif hits 
        relative to the bidirectionals. The calculated AUC is used as a proxy 
        for the enrichemnt of the TFs.

    Parameters
    ----------
    args : tuple 
        contains two arguments that are unpacked within this function:
            motif_file : string
                the name of a motif bed file contained within MOTIF_HITS

            millions_mapped : list or array
                contiains a list of floats that corresponds to the millions 
                mapped reads for each bam file

            fimo : boolean
                a switch that controls whether motif hits are generated on the 
                fly using fimo

            ranked_center_file : string
                the full path to a bed file containing the center of regions 
                sorted by DE-Seq p-value

            motif_hits : string
                the full path to a directory containing motif hits across the 
                genome

            plot : boolean
                a switch that controls whether TFEA generates plots for all 
                motifs or just significant ones
            
            padj_cutoff : float
                the cutoff value for calling a motif as significant

            logos : string
                the full path to a directory containing meme logos for motifs

            figuredir : string
                the full path to the figure directory within the output 
                directory where figures and plots are stored

            largewindow : float
                a user specified larger window used for plotting purposes and 
                to do some calculations regarding the user-provided regions of 
                interest

            smallwindow : float
                a user specified smaller window used for plotting purposes and 
                to do some calculations regarding the user-provided regions of 
                interest

    Returns
    -------
    AUC : float
        the area under the curve (AUC) for a given tf

    p-value : float
        theoretical p-value calculated by comparing the observed value  to the
        distribution of all simulations (default:1000)
                 
    Raises
    ------
    ValueError
        when tf does not have any hits
    '''
    motif_file = args[0]
    millions_mapped = args[1]
    gc_array = args[2]
    fimo  = args[3]
    ranked_center_file = args[4]
    motif_hits = args[5]
    plot = args[6] 
    padj_cutoff = args[7]
    logos = args[8]
    figuredir = args[9]
    largewindow = args[10]
    smallwindow = args[11]
    genomefasta = args[12]
    tempdir = args[13]
    motifdatabase = args[14]

    if fimo:
        ranked_center_distance_file = fimo_distance(
                                        ranked_center_file=ranked_center_file,
                                        motif_file=motif_file,
                                        genomefasta=genomefasta, 
                                        largewindow=largewindow, 
                                        tempdir=tempdir, 
                                        figuredir=figuredir,
                                        motifdatabase=motifdatabase)
    else:
        ranked_center_distance_file = independent_functions.motif_distance_bedtools_closest(
                                        ranked_center_file=ranked_center_file,
                                        motif_path=motif_hits+motif_file)

    distances = []
    ranks = []
    pvals= []
    pos = 0
    fc = []

    with open(ranked_center_distance_file) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            distance = int(line[-1])
            rank = int(line[5])
            pvals.append(float(line[3]))
            fc.append(float(line[4]))
            ranks.append(rank)
            distances.append(distance)
            if distance < largewindow:
                pos+=1

    #sort distances based on the ranks from TF bed file
    #and calculate the absolute distance
    sorted_distances = [x for _,x in sorted(zip(ranks, distances))]
    distances_abs = [abs(x) for x in sorted_distances] 

    #filter any TFs/files without and hits
    if len(distances_abs) == 0:
        return "no hits"

    #Get -exp() of distance and get cumulative scores
    score = [math.exp(-x) for x in distances_abs] 
    total = float(sum(score))
    normalized_score = [x/total for x in score]
    cumscore = np.cumsum(normalized_score)

    #The AUC is the relative to the "random" line
    actualES = np.trapz(cumscore) - (0.5*len(cumscore))

    #Calculate random AUC
    simES = independent_functions.permutations_auc(distances=normalized_score)

    ##significance calculator                                                                                                                                                            
    mu = np.mean(simES)
    NES = actualES/abs(mu)
    sigma = np.std(simES)


    if actualES > 0:
        p = 1-norm.cdf(actualES,mu,sigma)
    else:
        p = norm.cdf(actualES,mu,sigma)

    plot_individual_graphs(plot=plot, padj_cutoff=padj_cutoff,
                            figuredir=figuredir, logos=logos, 
                            largewindow=largewindow, 
                            smallwindow=smallwindow,
                            distances_abs=distances_abs, 
                            sorted_distances=sorted_distances,ranks=ranks,
                            pvals=pvals, fc=fc, cumscore=cumscore, 
                            motif_file=motif_file, p=p, simES=simES, 
                            actualES=actualES, gc_array=gc_array)

    return [motif_file.split('.bed')[0],actualES,NES,p,pos]
#==============================================================================

#==============================================================================
def calculate_es_youden_rank(args):
    '''This function calculates the AUC for any TF based on the TF motif hits 
        relative to the bidirectionals. The calculated AUC is used as a proxy 
        for the enrichemnt of the TFs.

    Parameters
    ----------
    args : tuple 
        contains two arguments that are unpacked within this function:
            motif_file : string
                the name of a motif bed file contained within MOTIF_HITS

            millions_mapped : list or array
                contiains a list of floats that corresponds to the millions 
                mapped reads for each bam file

            fimo : boolean
                a switch that controls whether motif hits are generated on the 
                fly using fimo

            ranked_center_file : string
                the full path to a bed file containing the center of regions 
                sorted by DE-Seq p-value

            motif_hits : string
                the full path to a directory containing motif hits across the 
                genome

            plot : boolean
                a switch that controls whether TFEA generates plots for all 
                motifs or just significant ones
            
            padj_cutoff : float
                the cutoff value for calling a motif as significant

            logos : string
                the full path to a directory containing meme logos for motifs

            figuredir : string
                the full path to the figure directory within the output 
                directory where figures and plots are stored

            largewindow : float
                a user specified larger window used for plotting purposes and 
                to do some calculations regarding the user-provided regions of 
                interest

            smallwindow : float
                a user specified smaller window used for plotting purposes and 
                to do some calculations regarding the user-provided regions of 
                interest

    Returns
    -------
    AUC : float
        the area under the curve (AUC) for a given tf

    p-value : float
        theoretical p-value calculated by comparing the observed value  to the
        distribution of all simulations (default:1000)
                 
    Raises
    ------
    ValueError
        when tf does not have any hits
    '''
    motif_file = args[0]
    millions_mapped = args[1]
    gc_array = args[2]
    fimo  = args[3]
    ranked_center_file = args[4]
    motif_hits = args[5]
    plot = args[6] 
    padj_cutoff = args[7]
    logos = args[8]
    figuredir = args[9]
    largewindow = args[10]
    smallwindow = args[11]
    genomefasta = args[12]
    tempdir = args[13]
    motifdatabase = args[14]

    if fimo:
        ranked_center_distance_file = fimo_distance(
                                        ranked_center_file=ranked_center_file,
                                        motif_file=motif_file,
                                        genomefasta=genomefasta, 
                                        largewindow=largewindow, 
                                        tempdir=tempdir, 
                                        figuredir=figuredir,
                                        motifdatabase=motifdatabase)
    else:
        ranked_center_distance_file = independent_functions.motif_distance_bedtools_closest(
                                        ranked_center_file=ranked_center_file,
                                        motif_path=motif_hits+motif_file)

    distances = []
    ranks = []
    pvals= []
    pos = 0
    fc = []

    with open(ranked_center_distance_file) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            distance = int(line[-1])
            rank = int(line[5])
            pvals.append(float(line[3]))
            fc.append(float(line[4]))
            ranks.append(rank)
            distances.append(distance)
            if distance < largewindow:
                pos+=1

    #sort distances based on the ranks from TF bed file
    #and calculate the absolute distance
    rank_number = len(ranks)
    ranks = [float(rank)/rank_number for rank in ranks]
    sorted_distances = [x for _,x in sorted(zip(ranks, distances))]
    distances_abs = [abs(x) for x in sorted_distances] 

    #filter any TFs/files without and hits
    if len(distances_abs) == 0:
        return "no hits"

    #Get -exp() of distance and get cumulative scores
    score = [math.exp(-x) for x in distances_abs] 
    total = float(sum(score))
    normalized_score = [x/total for x in score]
    cumscore = np.cumsum(normalized_score)

    trend = np.arange(0,1,1/len(ranks))

    #The AUC is the relative to the "random" line
    youden = cumscore - trend
    youden_max = max(cumscore - trend)
    rank_max = np.where(youden == youden_max)[0][0]/len(youden)

    #Calculate random AUC
    simES = independent_functions.permutations_youden_rank(
                                                    distances=normalized_score, 
                                                    trend=trend)

    ##significance calculator                                                                                                                                                            
    mu = np.mean(simES)
    NES = youden/abs(mu)
    sigma = np.std(simES)


    if youden > 0:
        p = 1-norm.cdf(rank_max,mu,sigma)
    else:
        p = norm.cdf(rank_max,mu,sigma)

    plot_individual_graphs(plot=plot, padj_cutoff=padj_cutoff,
                            figuredir=figuredir, logos=logos, 
                            largewindow=largewindow, 
                            smallwindow=smallwindow,
                            distances_abs=distances_abs, 
                            sorted_distances=sorted_distances,ranks=ranks,
                            pvals=pvals, fc=fc, cumscore=cumscore, 
                            motif_file=motif_file, p=p, simES=simES, 
                            actualES=rank_max, gc_array=gc_array)

    return [motif_file.split('.bed')[0], rank_max, NES, p, pos]
#==============================================================================

#==============================================================================
def plot_individual_graphs(plot=bool(), padj_cutoff=float(),
                            figuredir=str(), logos=str(), 
                            largewindow=float(), 
                            smallwindow=float(),
                            distances_abs=list(), sorted_distances=list(),
                            ranks=list(),pvals=list(),fc=list(), 
                            cumscore=list(), motif_file='',p=float(),
                            simES=list(), actualES=float(), gc_array=list()):
    '''This function plots all TFEA related graphs for an individual motif

    Parameters
    ----------
    plot : boolean
        switch to determine whether all motifs are plotted or just significant
        ones
    
    padj_cutoff : float
        the cutoff value that determines signficance. This is technically
        evaluated against the non-adjusted p-value since p-adj has not been
        calculated yet at this point
    
    figuredir : string
        the full path to the figuredir within the ouptut directory containing
        all figures and plots

    logos : string
        full path to a directory containing meme logos. These are copied to the
        output directory to be displayed in the html results.

    largewindow : float
        a user specified larger window used for plotting purposes and to do
        some calculations regarding the user-provided regions of interest

    smallwindow : float
        a user specified smaller window used for plotting purposes and to do
        some calculations regarding the user-provided regions of interest

    distances_abs : list or array
        a sorted (based on rank) list of absolute motif to region center 
        distances

    sorted_distances : list or array
        a sorted (based on rank) list of motif distances. Negative corresponds
        to upstream of region

    ranks : list or array
        a list of ranks that correspond to each region of interest. The ranking
        comes from DE-Seq p-value

    pvals : list or array
        a list of DE-Seq p-values for all regions

    p : float
        the p-value of the given motif based on simulations

    simES : list
        a list of simulated 'enrichment' scores calculated by randomizing
        region rank

    actualES : float
        the observed 'enchrichment' score. Can be calculated in many different
        ways..

    gc_array : list
        an array of gc richness of regions of interest. It is recommended that
        this array be no larger than 1000 bins.

    Returns
    -------
    None
    '''
    #Only plot things if user selects to plot all or if the pvalue is less than
    #the cutoff
    if plot or p < padj_cutoff:

        #Filter distances into quartiles for plotting purposes
        q1 = round(np.percentile(np.arange(1, len(distances_abs),1), 25))
        q3 = round(np.percentile(np.arange(1, len(distances_abs),1), 75))
        updistancehist = distances_abs[0:int(q1)]
        middledistancehist =  distances_abs[int(q1):int(q3)]
        downdistancehist = distances_abs[int(q3):len(distances_abs)]

        
        #Get log pval to plot for rank metric
        sorted_pval = [x for _,x in sorted(zip(ranks, pvals))]
        sorted_fc = [x for _,x in sorted(zip(ranks, fc))]
        try:
            logpval = [math.log(x,10) if y < 1 else -math.log(x,10) \
                        for x,y in zip(sorted_pval,sorted_fc)]
        except ValueError:
            logpval = sorted_pval

        #Get motif logos from logo directory                              
        if 'HO_' in motif_file:
            os.system("scp " + logos 
                        + motif_file.split('.bed')[0].split('HO_')[1] 
                        + "_direct.png " + figuredir)

            os.system("scp " + logos 
                        + motif_file.split('.bed')[0].split('HO_')[1] 
                        + "_revcomp.png " + figuredir)
        else:
            os.system("scp " + logos + motif_file.split('.bed')[0] 
                        + "_direct.png " + figuredir)

            os.system("scp " + logos + motif_file.split('.bed')[0] 
                        + "_revcomp.png " + figuredir)

        #Plot the enrichment plot
        independent_functions.enrichment_plot(largewindow=largewindow,
                                            smallwindow=smallwindow,
                                            figuredir=figuredir,
                                            cumscore=cumscore, 
                                            sorted_distances=sorted_distances, 
                                            logpval=logpval, 
                                            updistancehist=updistancehist, 
                                            downdistancehist=downdistancehist, 
                                            gc_array=gc_array,
                                            motif_file=motif_file)

        #Plot the simulation plot
        independent_functions.simulation_plot(figuredir=figuredir, simES=simES,
                                            actualES=actualES,
                                            motif_file=motif_file)

        #Plot the distance distribution histograms plot
        independent_functions.distance_distribution_plot(
                                        largewindow=largewindow,
                                        smallwindow=smallwindow,
                                        figuredir=figuredir,
                                        updistancehist=updistancehist, 
                                        middledistancehist=middledistancehist,
                                        downdistancehist=downdistancehist, 
                                        motif_file=motif_file)
#==============================================================================

#==============================================================================
def plot_global_graphs(padj_cutoff=float(), label1=str(), label2=str(), 
                        figuredir=str(), TFresults=list()):
    '''This function plots graphs that are displayed on the main results.html 
        filethat correspond to results relating to all analyzed TFs.

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

    ESlist = [i[1] for i in TFresults]
    NESlist = [i[2] for i in TFresults]
    PVALlist = [i[3] for i in TFresults]
    POSlist = [math.log(i[4],10) if i[4] != 0 else 0 for i in TFresults]
    PADJlist = [i[5] for i in TFresults]

    sigx = [x for x, p in zip(ESlist, PADJlist) if p < padj_cutoff]
    sigy = [p for x, p in zip(ESlist, PADJlist) if p < padj_cutoff]

    MAy = sigx
    MAx = [x for x, p in zip(POSlist, PADJlist) if p < padj_cutoff]

    #Creates a moustache plot of the global PADJs vs. ESs                                                                                                                                                                                       
    independent_functions.moustache_plot(figuredir=figuredir, ESlist=ESlist,
                                            PADJlist=PADJlist, sigx=sigx, 
                                            sigy=sigy)

    #Creates a histogram of p-values                                                                                                                                                                                                            
    independent_functions.pval_histogram_plot(figuredir=figuredir, 
                                                PVALlist=PVALlist)

    #Creates an MA-plot with NES on Y-axis and positive hits on X-axis                                                                                                                                                                                                      
    independent_functions.MA_plot(POSlist=POSlist, ESlist=ESlist, MAx=MAx, 
                                    MAy=MAy, figuredir=figuredir, label1=label1,
                                    label2=label2)
#==============================================================================

#==============================================================================
def get_gc_array(tempdir=str(), ranked_file=str(), genomefasta=str(), 
                    window=int(), bins=1000):
    '''This function calculates gc content over all eRNAs. It uses the 
    LARGEWINDOW variable within the config file instead of the whole region. 
    It performs a running average of window size = total_regions/bins

    Parameters
    ----------
    ranked_bed_file : str
        full path to a tab delimited bedfile with regions ranked by some metric
        of differential transcription
    window : int
        the size of the window (in bp) that you wish to calculate gc content 
        for
    bins : int 
        the number of equal sized bins to compute for the gc array
        
    Returns
    -------
    final_array : numpy array
        a list of gc content averaged to have bins equal to bins specified in
        paramters
    '''
    #First, create a bed file with the correct coordinates centered on the 
    #given regions with the specified window size on either side
    outfile = open(os.path.join(tempdir, "ranked_file.windowed.bed"),'w')
    with open(os.path.join(tempdir, "ranked_file.bed")) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            chrom,start,stop = line[:3]
            center = (int(start)+int(stop))/2
            newstart = center - window
            newstop = center + window
            outfile.write(chrom + '\t' + str(newstart) + '\t' + str(newstop) 
                            + '\t' + '\t'.join(line[3:]) + '\n')
    outfile.close()

    #Convert the bed file created above into a fasta file using meme (contained
    #within the combine_bed file for some reason)
    ranked_file_windowed_fasta = independent_functions.getfasta(
                                    genomefasta=genomefasta, 
                                    tempdir=tempdir,
                                    bedfile=os.path.join(tempdir, 
                                                "ranked_file.windowed.bed"))

    #Create a gc_array which simply contains all sequences in fasta file c
    #ollapsed into 1.0 for G/C and 0.0 for A/T
    gc_array = []
    with open(ranked_file_windowed_fasta) as F:
        for line in F:
            if '>' not in line:
                line = line.strip('\n')
                gc_content = independent_functions.convert_sequence_to_array(
                                                                sequence=line)
                gc_array.append(gc_content)

    #The length of each bin is equal to the total positions over the number of 
    #desired bins
    binwidth = len(gc_array)/bins

    #Collapse the gc_array into an array containing the correct number of bins 
    #(specified by user)
    final_array = []
    ##First, step through the gc_array with binwidth step size (i)
    for i in range(0,len(gc_array),binwidth):
        ##Initialize a position_average list that will store the mean value for
        #each position (along the window)
        position_average = []
        ##Now we step through the total window size (window*2) position by 
        #position
        for k in range(window*2):
            ##new_array stores for each position in the window, a binwidth 
            #amount of data points to be averaged
            new_array = []
            ##Now, if we are not at the end of the gc_array, we will step 
            #through for each position, a binwidth amount of values and store 
            #them in new_array
            if i+binwidth < len(gc_array):
                for j in range(i,i+binwidth):
                    new_array.append(gc_array[j][k])
            else:
                for j in range(i,len(gc_array)):
                    new_array.append(gc_array[j][k])
            ##Finally, we simply append the average of the new_array into the 
            #position_average list
            position_average.append(np.mean(new_array))
        ##And this poition_average list is appended to the actual GC_ARRAY 
        #within the config file for later use
        final_array.append(position_average)

    final_array = np.array(final_array).transpose()

    return final_array
#==============================================================================

#==============================================================================
def calculate(tempdir=str(), outputdir=str(), ranked_center_file=str(),
                bam1=list(), bam2=list(), config_dict=dict(), 
                motif_hits=str(), singlemotif=str(), 
                COMBINEtime=float(), COUNTtime=float(), DESEQtime=float(), 
                CALCULATEtime=float(), fimo=bool(), plot=bool(), 
                padj_cutoff=float(), logos=str(), figuredir=str(),
                largewindow=float(), smallwindow=float(), genomefasta=str(),
                motifdatabase=str(), pool=bool(), label1=str(), label2=str()):
    '''This function performs the bulk of transcription factor enrichment
        analysis (TFEA), it does the following:
            1. GC distribution across regions for plotting
            2. Millions mapped calculation for meta eRNA plots
            3. Motif distance calculation to the center of inputted regions
            4. Enrichment score calculation via AUC method
            5. Random shuffle simulation and recalculation of enrichment score
            6. Plotting and generation of html report

    Parameters
    ----------
    tempdir : string
        the full path to the tempdir directory within the outputdir created 
        by TFEA

    outputdir : string
        the full path to the output directory created by TFEA

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

    '''
    print "Calculating GC content of regions..."
    #This line gets an array of GC values for all inputted regions
    gc_array = get_gc_array(
                        ranked_file=os.path.join(tempdir, "ranked_file.bed"),
                        window=int(largewindow), tempdir=tempdir, 
                        genomefasta=genomefasta)

    print "done\nCalculating millions mapped reads for bam files..."
    #Here we determine how many cpus to use for parallelization
    cpus = mp.cpu_count()
    if cpus > 64:
        cpus = 64

    #Here we calculate millions mapped reads for use with the metaeRNA module
    p = Pool(cpus)
    args = [(x,tempdir) for x in bam1+bam2]
    millions_mapped = p.map(independent_functions.samtools_flagstat,args)

    print "done\nFinding motif hits in regions..."
    if singlemotif == False:
        TFresults = list()
        if pool:
            args = [(x, millions_mapped, gc_array, fimo, ranked_center_file,
                        motif_hits, plot, padj_cutoff, logos, figuredir,
                        largewindow, smallwindow, genomefasta, tempdir, 
                        motifdatabase) for x in os.listdir(motif_hits)]
            p = Pool(cpus)
            TFresults = p.map(calculate_es_auc, args)
        else:
            for motif_file in os.listdir(motif_hits):
                results = calculate_es_auc((motif_file, millions_mapped, 
                        gc_array, fimo, ranked_center_file, motif_hits, plot, 
                        padj_cutoff, logos, figuredir, largewindow, 
                        smallwindow, genomefasta, tempdir, motifdatabase))
                if results != "no hits":
                    TFresults.append(results)
                else:
                    print "No motifs within specified window for: ", motif_file


        CALCULATEtime = time.time()-CALCULATEtime
        TFresults = independent_functions.padj_bonferroni(TFresults=TFresults)
        plot_global_graphs(padj_cutoff=padj_cutoff, label1=label1, 
                            label2=label2, figuredir=figuredir, 
                            TFresults=TFresults)
        independent_functions.create_text_output(TFresults=TFresults, 
                                                    outputdir=outputdir)
        independent_functions.create_html_output(TFresults=TFresults, 
                                                    config_dict=config_dict, 
                                                    COMBINEtime=COMBINEtime, 
                                                    COUNTtime=COUNTtime, 
                                                    DESEQtime=DESEQtime, 
                                                    CALCULATEtime=CALCULATEtime)

    #Note if you set the SINGLEMOTIF variable to a specific TF, this program 
    #will be unable to determine an PADJ for the given motif.
    else:
        results = calculate_es_auc((singlemotif, millions_mapped, gc_array, 
                        fimo, ranked_center_file, motif_hits, plot, 
                        padj_cutoff, logos, figuredir, largewindow, 
                        smallwindow, genomefasta, tempdir, motifdatabase))
        independent_functions.create_single_motif_html(results=results, 
                                                        outputdir=outputdir)

    print "done"
#==============================================================================

#==============================================================================
def fimo_distance(motif_file=str(),genomefasta=str(), largewindow=float(),
                tempdir=str(), motifdatabase=str(), figuredir=str(),
                ranked_center_file=str()):
    ranked_fullregions_file = independent_functions.get_regions(
                                        tempdir=tempdir, 
                                        ranked_center_file=ranked_center_file,
                                        largewindow=largewindow)

    ranked_fasta_file = independent_functions.getfasta(
                                            genomefasta=genomefasta, 
                                            tempdir=tempdir,
                                            bedfile=ranked_fullregions_file)

    bg_file = independent_functions.get_bgfile(fastafile=ranked_fasta_file,
                                                tempdir=tempdir)

    fimo_file = independent_functions.fimo(bgfile=bg_file, 
                                            motif=motif_file,
                                            fastafile=ranked_fasta_file,
                                            tempdir=tempdir, 
                                            motifdatabase=motifdatabase)

    independent_functions.meme2images(motif=motif_file, 
                                        motifdatabase=motifdatabase, 
                                        figuredir=figuredir)

    ranked_center_distance_file = independent_functions.fimo_parse(
                                                        fimo_file=fimo_file,
                                                        motif_file=motif_file,
                                                        largewindow=largewindow,
                                                        tempdir=tempdir)

    return ranked_center_distance_file
#==============================================================================

#==============================================================================