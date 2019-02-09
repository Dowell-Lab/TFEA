#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''This file contains functions that plot results of enrichment analysis
'''
#==============================================================================
__author__ = 'Jonathan D. Rubin and Rutendo F. Sigauke'
__credits__ = ['Jonathan D. Rubin', 'Rutendo F. Sigauke', 'Jacob T. Stanley',
                'Robin D. Dowell']
__maintainer__ = 'Jonathan D. Rubin'
__email__ = 'Jonathan.Rubin@colorado.edu'
__version__ = '4.0'

#Imports
#==============================================================================
import os
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


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
    import config
    #Only plot things if user selects to plot all or if the pvalue is less than
    #the cutoff
    if plot or p < padj_cutoff:
        if config.SCANNER == 'homer':
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
        enrichment_plot(largewindow=largewindow,
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
        simulation_plot(figuredir=figuredir, simES=simES, actualES=actualES,
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
    moustache_plot(figuredir=figuredir, ESlist=ESlist,
                                            PADJlist=PVALlist, sigx=sigx, 
                                            sigy=sigy)

    #Creates a histogram of p-values                                                                                                                                                                                                            
    pval_histogram_plot(figuredir=figuredir, PVALlist=PVALlist)

    #Creates an MA-plot with NES on Y-axis and positive hits on X-axis                                                                                                                                                                                                      
    MA_plot(POSlist=POSlist, ESlist=ESlist, MAx=MAx, 
                                    MAy=MAy, figuredir=figuredir, label1=label1,
                                    label2=label2)

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
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib.gridspec as gridspec
    import traceback

    import config
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
        figure_title = motif_file.split('.bed')[0] + ' Enrichment Plot'
        lineplot(title=figure_title, ax=lineplot, xvals=xvals, yvals=cumscore, 
                    xlimits=limits)
        scatterplot(ax=scatterplot, xvals=xvals, yvals=sorted_distances, 
                    xlimits=None, largewindow=largewindow)
        barplot(ax=barplot, xvals=xvals, colorarray=score, xlimits=limits)
        fillplot(ax=fillplot, xvals=xvals, yvals=logpval, xlimits=limits, 
                ylimits=None)

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
        print(traceback.print_exc())
        raise e

#==============================================================================
def lineplot(title=None, ax=None, xvals=None, yvals=None, xlimits=None):
    #This is the enrichment score plot (i.e. line plot)
    # ax0 = plt.subplot(gs[0])
    ax.plot(xvals,yvals,color='green')
    ax.plot([0, 1],[0, 1], '--', alpha=0.75)
    ax.set_title(title, fontsize=14)
    ax.set_ylabel('Enrichment Score (ES)', fontsize=10)
    ax.tick_params(axis='y', which='both', left='on', right='off', 
                    labelleft='on')
    ax.tick_params(axis='x', which='both', bottom='off', top='off', 
                    labelbottom='off')
    ax.set_ylim([0,1])
    ax.set_xlim(xlimits)

#==============================================================================
def scatterplot(ax=None, xvals=None, yvals=None, xlimits=None, largewindow=None):
    import matplotlib.pyplot as plt
    #This is the barplot right below the enrichment score line plot
    ax.scatter(xvals, yvals, edgecolor="", color="black", 
                    s=10, alpha=0.25)
    ax.tick_params(axis='y', which='both', left='off', right='off', 
                    labelleft='on') 
    ax.tick_params(axis='x', which='both', bottom='off', top='off', 
                    labelbottom='off')
    plt.yticks([-int(largewindow),0,int(largewindow)],
                [str(-int(largewindow)/1000.0),'0',\
                str(int(largewindow)/1000.0)])
    ax.set_xlim(xlimits)
    ax.set_ylim([-int(largewindow),int(largewindow)])
    ax.set_ylabel('Distance (kb)', fontsize=10)

#==============================================================================
def barplot(ax=None, xvals=None, colorarray=None, xlimits=None):
    import matplotlib
    import matplotlib.cm as cm
    norm    = matplotlib.colors.Normalize(vmin=min(colorarray), 
                                            vmax=max(colorarray))
    cmap    = cm.Greys
    m       = cm.ScalarMappable(norm=norm, cmap=cmap)
    colors  = [m.to_rgba(c) for c in colorarray] 
    ax.bar(xvals,[1 for x in xvals], edgecolor="", color=colors)
    ax.tick_params(axis='y', which='both', left='off', right='off', 
                    labelleft='off') 
    ax.tick_params(axis='x', which='both', bottom='off', top='off', 
                    labelbottom='off')
    ax.set_xlim(xlimits)
    ax.set_ylim([0,1])
    ax.set_ylabel('Score', fontsize=10)

#==============================================================================
def fillplot(ax=None, xvals=None, yvals=None, xlimits=None):
    #This is the rank metric fill plot
    ax.fill_between(xvals,0,yvals,facecolor='grey',edgecolor="")
    ax.tick_params(axis='y', which='both', left='on', right='off', 
                    labelleft='on')
    ax.tick_params(axis='x', which='both', bottom='off', top='off', 
                    labelbottom='on')
    ylim = math.fabs(max([x for x in yvals if -500 < x < 500],key=abs))
    ax.set_ylim(ylimits)
    ax.yaxis.set_ticks([int(-ylim),0,int(ylim)])
    ax.set_xlim(xlimits)
    ax.set_xlabel('Relative Rank (n='+str(len(xvals))+')', 
                    fontsize=14)
    ax.set_ylabel('Rank Metric',fontsize=10)
    # try:
    #     ax2.axvline(len(updistancehist)+1,color='green',alpha=0.25)
    # except ValueError:
    #     pass
    # try:
    #     ax2.axvline(len(xvals) - len(downdistancehist), color='purple', 
    #                 alpha=0.25)
    # except ValueError:
    #     pass

#==============================================================================
def leftover_plotting_scripts():

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

#==============================================================================
def distance_heatmap_plot(figuredir=None, motif_file=None, q1_distances=None, 
                        q2_distances=None, q3_distances=None, 
                        q4_distances=None, bins=None, dpi=None):
    '''
    '''
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

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
def simulation_plot(ax=None, simulated=None, observed=None, title=None):
    '''
    '''
    maximum = max(simulated)
    minimum = min(simulated)
    ax.hist(simulated,bins=100)
    width = (maximum-minimum)/100.0
    rect = ax.bar(observed,ax.get_ylim()[1],color='red',width=width*2)[0]
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width()/2., 1.05*height, 
                'Observed', ha='center', va='bottom')

    ax.set_xlim([min(minimum,observed)-(width*40), \
                max(maximum,observed)+(width*40)])

    ax.set_ylim([0,(1.05*height)+5])
    ax.tick_params(axis='y', which='both', left='off', right='off', 
                    labelleft='on')

    ax.tick_params(axis='x', which='both', bottom='off', top='off', 
                    labelbottom='on')

    if title != None:
        plt.title(title,fontsize=14)
    ax.set_ylabel('Number of Simulations',fontsize=14)
    ax.set_xlabel('Enrichment',fontsize=14)
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
    import math
    import matplotlib.pyplot as plt

    import config

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
    import matplotlib.pyplot as plt

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
    import matplotlib.pyplot as plt
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
    import matplotlib.pyplot as plt
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