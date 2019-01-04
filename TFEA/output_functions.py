#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''This file contains a list of functions associated with outputting a ranked
    bed file based on some metric.
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
import numpy as np

#Functions
#==============================================================================
def plot_individual_graphs(plot=None, padj_cutoff=None,
                            figuredir=None, logos=None, 
                            largewindow=None, score=None, 
                            smallwindow=None,
                            distances_abs=None, sorted_distances=None,
                            ranks=None, pvals=None, fc=None, 
                            cumscore=None, motif_file=None, p=None,
                            simES=None, actualES=None, gc_array=None,
                            meta_profile_dict=None, label1=None, label2=None):
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
        if config.HOMER:
            pass
        else:
            #Get motif logos from logo directory
            os.system("scp '" + logos 
                        + motif_file.strip('.sorted.distance.bed').strip('HO_') 
                        + "_direct.png' " + figuredir)

            os.system("scp '" + logos 
                        + motif_file.strip('.sorted.distance.bed').strip('HO_') 
                        + "_revcomp.png' " + figuredir)
        

        #Filter distances into quartiles for plotting purposes
        q1 = int(round(np.percentile(np.arange(1, len(sorted_distances),1), 25)))
        q2 = int(round(np.percentile(np.arange(1, len(sorted_distances),1), 50)))
        q3 = int(round(np.percentile(np.arange(1, len(sorted_distances),1), 75)))
        q1_distances = [x for x in sorted_distances[:q1] if x <= largewindow]
        q2_distances = [x for x in sorted_distances[q1:q2] if x <= largewindow]
        q3_distances = [x for x in sorted_distances[q2:q3] if x <= largewindow]
        q4_distances = [x for x in sorted_distances[q3:] if x <= largewindow]

        
        #Get log pval to plot for rank metric
        sorted_pval = [x for _,x in sorted(zip(ranks, pvals))]
        sorted_fc = [x for _,x in sorted(zip(ranks, fc))]
        logpval = list()
        for x,y in zip(sorted_pval,sorted_fc):
            try:
                if y < 1:
                    logpval.append(math.log(x,10))
                else:
                    logpval.append(-math.log(x,10))
            except:
                logpval.append(0.0)

        #Plot the enrichment plot
        independent_functions.enrichment_plot(largewindow=largewindow,
                                            smallwindow=smallwindow,
                                            figuredir=figuredir,
                                            cumscore=cumscore, 
                                            sorted_distances=sorted_distances, 
                                            logpval=logpval,  
                                            gc_array=gc_array, score=score, 
                                            motif_file=motif_file, 
                                            q1_distances=q1_distances, 
                                            q2_distances=q2_distances, 
                                            q3_distances=q3_distances, 
                                            q4_distances=q4_distances,
                                            meta_profile_dict=meta_profile_dict,
                                            dpi=None, label1=label1, 
                                            label2=label2)

        #Plot the simulation plot
        independent_functions.simulation_plot(figuredir=figuredir, simES=simES,
                                            actualES=actualES,
                                            motif_file=motif_file)

        # #Plot the distance distribution histograms plot
        # independent_functions.meta_plot(figuredir=figuredir, 
        #                                     meta_profile_dict=meta_profile_dict,
        #                                     motif_file=motif_file, 
        #                                     q1_distances=q1_distances, 
        #                                     q2_distances=q2_distances, 
        #                                     q3_distances=q3_distances, 
        #                                     q4_distances=q4_distances, 
        #                                     bins=100, 
        #                                     dpi=None)

#==============================================================================
def plot_global_graphs(padj_cutoff=None, label1=None, label2=None, 
                        figuredir=None, TFresults=None):
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
    POSlist = [i[4] for i in TFresults]
    PADJlist = [i[5] for i in TFresults]

    POSlist = [math.log(x,10) if x > 0 else 0.0 for x in POSlist]

    sigx = [x for x, p in zip(ESlist, PADJlist) if p < padj_cutoff]
    sigy = [p for x, p in zip(ESlist, PADJlist) if p < padj_cutoff]

    MAy = sigx
    MAx = [x for x, p in zip(POSlist, PADJlist) if p < padj_cutoff]

    #Creates a moustache plot of the global PADJs vs. ESs                                                                                                                                                                                       
    independent_functions.moustache_plot(figuredir=figuredir, ESlist=ESlist,
                                            PADJlist=PVALlist, sigx=sigx, 
                                            sigy=sigy)

    #Creates a histogram of p-values                                                                                                                                                                                                            
    independent_functions.pval_histogram_plot(figuredir=figuredir, 
                                                PVALlist=PVALlist)

    #Creates an MA-plot with NES on Y-axis and positive hits on X-axis                                                                                                                                                                                                      
    independent_functions.MA_plot(POSlist=POSlist, ESlist=ESlist, MAx=MAx, 
                                    MAy=MAy, figuredir=figuredir, label1=label1,
                                    label2=label2)