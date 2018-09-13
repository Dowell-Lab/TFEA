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
def deseq_run(tempdir=config.TEMPDIR):
    '''Writes the DE-Seq script, runs the DE-Seq script, and plots the DE-Seq 
        MA plot

    Parameters
    ----------
    tempdir : string
        full path to temp directory in output directory (created by TFEA)

    Returns
    ----------
    None
    '''
    independent_functions.write_deseq_script()
    os.system("R < " + tempdir + "DESeq.R --no-save")
    independent_functions.plot_deseq_MA(tempdir+'DESeq.res.txt')
#==============================================================================

#==============================================================================
def calculate_es_auc(args, fimo=config.FIMO, 
                        ranked_center_file=config.RANKED_CENTER_FILE,
                        motif_hits=config.MOTIF_HITS, plot=config.PLOT,
                        padj_cutoff=config.PADJCUTOFF, logos=config.LOGOS,
                        figuredir=config.FIGUREDIR,
                        largewindow=config.LARGEWINDOW,
                        smallwindow=config.SMALLWINDOW):
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
        a switch that controls whether motif hits are generated on the fly
        using fimo

    ranked_center_file : string
        the full path to a bed file containing the center of regions sorted by
        DE-Seq p-value

    motif_hits : string
        the full path to a directory containing motif hits across the genome

    plot : boolean
        a switch that controls whether TFEA generates plots for all motifs or
        just significant ones
    
    padj_cutoff : float
        the cutoff value for calling a motif as significant

    logos : string
        the full path to a directory containing meme logos for motifs

    figuredir : string
        the full path to the figure directory within the output directory where
        figures and plots are stored

    largewindow : float
        a user specified larger window used for plotting purposes and to do
        some calculations regarding the user-provided regions of interest

    smallwindow : float
        a user specified smaller window used for plotting purposes and to do
        some calculations regarding the user-provided regions of interest

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
    motif_file, millions_mapped, gc_array = args
    if fimo:
        ranked_fullregions_file = independent_functions.get_regions()
        ranked_fasta_file = independent_functions.getfasta(
                                bedfile=ranked_fullregions_file)

        background_file = independent_functions.get_bgfile(
                                fastafile=ranked_fasta_file)

        fimo_file = independent_functions.fimo(background_file, motif_file,
                                                ranked_fasta_file)

        independent_functions.meme2images(motif_file)
        ranked_center_distance_file = independent_functions.fimo_distance(
                                                                    fimo_file,
                                                                    motif_file)
    else:
        ranked_center_distance_file = \
            independent_functions.motif_distance_bedtools_closest(
                                                        ranked_center_file,
                                                        motif_hits+motif_file)

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
    simES = independent_functions.permutations(normalized_score)

    ##significance calculator                                                                                                                                                            
    mu = np.mean(simES)
    NES = actualES/abs(mu)
    sigma = np.std(simES)


    if actualES > 0:
        p = 1-norm.cdf(actualES,mu,sigma)
    else:
        p = norm.cdf(actualES,mu,sigma)

    plot_individual_graphs(distances_abs=distances_abs, 
                            sorted_distances=sorted_distances,ranks=ranks,
                            pvals=pvals, fc=fc, cumscore=cumscore, 
                            motif_file=motif_file, p=p, simES=simES, 
                            actualES=actualES, gc_array=gc_array)

    return [motif_file.split('.bed')[0],actualES,NES,p,pos]
#==============================================================================

#==============================================================================
def plot_individual_graphs(plot=config.PLOT, padj_cutoff=config.PADJCUTOFF,
                            figuredir=config.FIGUREDIR, logos=config.LOGOS, 
                            largewindow=config.LARGEWINDOW, 
                            smallwindow=config.SMALLWINDOW,
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
        independent_functions.enrichment_plot(cumscore=cumscore, 
                                            sorted_distances=sorted_distances, 
                                            logpval=logpval, 
                                            updistancehist=updistancehist, 
                                            downdistancehist=downdistancehist, 
                                            gc_array=gc_array,
                                            motif_file=motif_file)

        #Plot the simulation plot
        independent_functions.simulation_plot(simES=simES,actualES=actualES,
                                            motif_file=motif_file)

        #Plot the distance distribution histograms plot
        independent_functions.distance_distribution_plot(
                                        updistancehist=updistancehist, 
                                        middledistancehist=middledistancehist,
                                        downdistancehist=downdistancehist, 
                                        motif_file=motif_file)
#==============================================================================

#==============================================================================
def plot_global_graphs(padj_cutoff=config.PADJCUTOFF, label1=config.LABEL1,
                        label2=config.LABEL2, figuredir=config.FIGUREDIR, 
                        TFresults=list()):
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
    TFresults = independent_functions.padj_bonferroni(TFresults=TFresults)

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
    independent_functions.moustache_plot(ESlist=ESlist,PADJlist=PADJlist,
                                            sigx=sigx,sigy=sigy)

    #Creates a histogram of p-values                                                                                                                                                                                                            
    independent_functions.pval_histogram_plot(PVALlist=PVALlist)

    #Creates an MA-plot with NES on Y-axis and positive hits on X-axis                                                                                                                                                                                                      
    independent_functions.MA_plot(POSlist=POSlist, ESlist=ESlist, MAx=MAx, 
                                    MAy=MAy)
#==============================================================================

#==============================================================================
def get_gc_array(tempdir=config.TEMPDIR, ranked_file='',
                    window=int(config.LARGEWINDOW),bins=1000):
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
    outfile = open(tempdir+"ranked_file.windowed.bed",'w')
    with open(tempdir + "ranked_file.bed") as F:
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
    ranked_file_windowed_fasta = independent_functions.getfasta(tempdir
                                                +"ranked_file.windowed.bed")

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
def calculate(tempdir=config.TEMPDIR, outputdir=config.OUTPUTDIR, 
                bam1=config.BAM1, bam2=config.BAM2, 
                motif_hits=config.MOTIF_HITS, singlemotif=config.SINGLEMOTIF, 
                COMBINEtime=float(), COUNTtime=float(), DESEQtime=float(), 
                CALCULATEtime=float()):
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
                        ranked_file=os.path.join(tempdir, "ranked_file.bed"))

    print "done\nCalculating millions mapped reads for bam files..."
    #Here we determine how many cpus to use for parallelization
    cpus = mp.cpu_count()
    if cpus > 64:
        cpus = 64

    #Here we calculate millions mapped reads for use with the metaeRNA module
    p = Pool(cpus)
    args = [(x) for x in bam1+bam2]
    millions_mapped = p.map(independent_functions.samtools_flagstat,args)

    print "done\nFinding motif hits in regions..."
    if singlemotif == False:
        TFresults = list()
        if config.POOL:
            a = time.time()
            args = [(x,millions_mapped,gc_array) for x in os.listdir(motif_hits)]
            p = Pool(cpus)
            TFresults = p.map(calculate_es_auc, args)
        else:
            for MOTIF_FILE in os.listdir(motif_hits):
                results = calculate_es_auc((MOTIF_FILE, millions_mapped,
                                                gc_array))
                if results != "no hits":
                    TFresults.append(results)
                else:
                    print "No motifs within specified window for: ", MOTIF_FILE
        CALCULATEtime = time.time()-CALCULATEtime
        TFresults = independent_functions.padj_bonferroni(TFresults=TFresults)
        independent_functions.create_text_output(TFresults=TFresults)
        independent_functions.create_html_output(TFresults=TFresults, 
                                                COMBINEtime=COMBINEtime, 
                                                COUNTtime=COUNTtime, 
                                                DESEQtime=DESEQtime, 
                                                CALCULATEtime=CALCULATEtime)

    #Note if you set the SINGLEMOTIF variable to a specific TF, this program 
    #will be unable to determine an PADJ for the given motif.
    else:
        results = calculate_es_auc((singlemotif, millions_mapped, gc_array))
        independent_functions.create_single_motif_html(results=results)

    print "done"
#==============================================================================

#==============================================================================
