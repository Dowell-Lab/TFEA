#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''This module only contains scripts that plot things in matplotlib. 
    All such scripts are in this module. No exceptions.
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
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter
matplotlib.use('Agg')
matplotlib.rcParams['savefig.dpi'] = 100
import subprocess
import warnings
import pathlib
import ujson
from statistics import mean, median

import numpy as np
from scipy import stats

from TFEA import exceptions

## GC Decorator
import gc
import functools
def force_gc(func):
    @functools.wraps(func)
    def force_gc_decorator(*args, **kwargs):
        # Do something before
        value = func(*args, **kwargs)
        # Do something after
        gc.collect(2)
        return value
    return force_gc_decorator

#Functions
#==============================================================================
@force_gc
def plot_individual_graphs(use_config=True, distances=None, figuredir=None, 
                            fimo_motifs=None, largewindow=1500, score=None, 
                            pvals=None, fcs=None, 
                            cumscore=None, sim_auc=None, auc=None,
                            meta_profile_dict=None, label1=None, label2=None, 
                            motif=None, offset=None, plot_format=None):
    '''This function plots all TFEA related graphs for an individual motif
    '''
    if use_config:
        from TFEA import config
        pvals = config.vars['PVALS']
        fcs = config.vars['FCS']
        figuredir=config.vars['FIGUREDIR']
        largewindow=config.vars['LARGEWINDOW']
        fimo_motifs = config.vars['FIMO_MOTIFS']
        label1 = config.vars['LABEL1']
        label2 = config.vars['LABEL2']
        meta_profile_dict = config.vars['META_PROFILE']
        plot_format = config.vars['PLOT_FORMAT']

    #Create MEME logos
    if fimo_motifs:
        meme_logo(fimo_motifs, motif, figuredir, plot_format=plot_format)
    else:
        print("No MEME database inputted, logos will not be displayed.")

    #Filter distances into quartiles for plotting purposes
    q1 = int(round(np.percentile(np.arange(1, len(distances),1), 25)))
    q2 = int(round(np.percentile(np.arange(1, len(distances),1), 50)))
    q3 = int(round(np.percentile(np.arange(1, len(distances),1), 75)))
    q1_distances = [x for x in distances[:q1] if x != '.']
    q1_meta_retain = [i for i,x in enumerate(distances[:q1]) if x != '.']
    q2_distances = [x for x in distances[q1:q2] if x != '.']
    q2_meta_retain = [i for i,x in enumerate(distances[q1:q2]) if x != '.']
    q3_distances = [x for x in distances[q2:q3] if x != '.']
    q3_meta_retain = [i for i,x in enumerate(distances[q2:q3]) if x != '.']
    q4_distances = [x for x in distances[q3:] if x != '.']
    q4_meta_retain = [i for i,x in enumerate(distances[q3:]) if x != '.']

    
    #Get log pval to plot for rank metric
    logpval = list()
    for x,y in zip(pvals,fcs):
        try:
            if y < 1:
                logpval.append(math.log(x,10))
            else:
                logpval.append(-math.log(x,10))
        except:
            logpval.append(0.0)


    #Set up plotting space
    plt.figure(figsize=(15.5,12))
    outer_gs = gridspec.GridSpec(2, 1, height_ratios=[2,1])
    enrichment_gs = gridspec.GridSpecFromSubplotSpec(4, 1, 
                                                    subplot_spec=outer_gs[0], 
                                                    height_ratios=[4, 1, 4, 2], 
                                                    hspace=.1)
    meta_gs = gridspec.GridSpecFromSubplotSpec(3, 4, subplot_spec=outer_gs[1], 
                                                    hspace=0.1, wspace=0.3)

    #Create Enrichment Plot
    xvals = np.linspace(start=0, stop=1, num=len(cumscore))
    xlimits = [0, 1]
    lineplot_ax = plt.subplot(enrichment_gs[0])
    barplot_ax = plt.subplot(enrichment_gs[1])
    scatterplot_ax = plt.subplot(enrichment_gs[2])
    figure_title = motif.split('.bed')[0] + ' Enrichment Plot'
    lineplot(title=figure_title, ax=lineplot_ax, xvals=xvals, yvals=cumscore, 
                xlimits=xlimits)
    
    barplot(ax=barplot_ax, xvals=xvals, colorarray=score, xlimits=xlimits)

    scatter_data = [(i,x) for i,x in zip(xvals, distances) if x != '.']
    scatter_x = [i for i,_ in scatter_data]
    scatter_y = [x for _,x in scatter_data]
    if len(logpval) != 0:
        scatterplot(ax=scatterplot_ax, xvals=scatter_x, yvals=scatter_y, 
                    xlimits=xlimits, largewindow=largewindow)
        fillplot_ax = plt.subplot(enrichment_gs[3])
        fillplot(ax=fillplot_ax, xvals=xvals, yvals=logpval, xlimits=xlimits)
    else:
        scatterplot(ax=scatterplot_ax, xvals=scatter_x, yvals=scatter_y, 
                    xlimits=xlimits, largewindow=largewindow, xlabel=True)

    #Initiate meta plots
    if type(meta_profile_dict) == pathlib.PosixPath: #or type(meta_profile_dict) == dict and len(meta_profile_dict) != 0:
        
        q1posprofile1 = [ujson.loads((meta_profile_dict / f'q1posprofile1_{i}').read_text()) for i in q1_meta_retain]
        q1posprofile1 = [x for x in map(mean, zip(*q1posprofile1))]
        
        q1negprofile1 = [ujson.loads((meta_profile_dict / f'q1negprofile1_{i}').read_text()) for i in q1_meta_retain]
        q1negprofile1 = [x for x in map(mean, zip(*q1negprofile1))]
        
        q1posprofile2 = [ujson.loads((meta_profile_dict / f'q1posprofile2_{i}').read_text()) for i in q1_meta_retain]
        q1posprofile2 = [x for x in map(mean, zip(*q1posprofile2))]
        
        q1negprofile2 = [ujson.loads((meta_profile_dict / f'q1negprofile2_{i}').read_text()) for i in q1_meta_retain]
        q1negprofile2 = [x for x in map(mean, zip(*q1negprofile2))]
        
        q2posprofile1 = [ujson.loads((meta_profile_dict / f'q2posprofile1_{i}').read_text()) for i in q2_meta_retain]
        q2posprofile1 = [x for x in map(mean, zip(*q2posprofile1))]
        
        q2negprofile1 = [ujson.loads((meta_profile_dict / f'q2negprofile1_{i}').read_text()) for i in q2_meta_retain]
        q2negprofile1 = [x for x in map(mean, zip(*q2negprofile1))]
        
        q2posprofile2 = [ujson.loads((meta_profile_dict / f'q2posprofile2_{i}').read_text()) for i in q2_meta_retain]
        q2posprofile2 = [x for x in map(mean, zip(*q2posprofile2))]
        
        q2negprofile2 = [ujson.loads((meta_profile_dict / f'q2negprofile2_{i}').read_text()) for i in q2_meta_retain]
        q2negprofile2 = [x for x in map(mean, zip(*q2negprofile2))]
        
        q3posprofile1 = [ujson.loads((meta_profile_dict / f'q3posprofile1_{i}').read_text()) for i in q3_meta_retain]
        q3posprofile1 = [x for x in map(mean, zip(*q3posprofile1))]
        
        q3negprofile1 = [ujson.loads((meta_profile_dict / f'q3negprofile1_{i}').read_text()) for i in q3_meta_retain]
        q3negprofile1 = [x for x in map(mean, zip(*q3negprofile1))]
        
        q3posprofile2 = [ujson.loads((meta_profile_dict / f'q3posprofile2_{i}').read_text()) for i in q3_meta_retain]
        q3posprofile2 = [x for x in map(mean, zip(*q3posprofile2))]
        
        q3negprofile2 = [ujson.loads((meta_profile_dict / f'q3negprofile2_{i}').read_text()) for i in q3_meta_retain]
        q3negprofile2 = [x for x in map(mean, zip(*q3negprofile2))]
        
        q4posprofile1 = [ujson.loads((meta_profile_dict / f'q4posprofile1_{i}').read_text()) for i in q4_meta_retain]
        q4posprofile1 = [x for x in map(mean, zip(*q4posprofile1))]
        
        q4negprofile1 = [ujson.loads((meta_profile_dict / f'q4negprofile1_{i}').read_text()) for i in q4_meta_retain]
        q4negprofile1 = [x for x in map(mean, zip(*q4negprofile1))]
        
        q4posprofile2 = [ujson.loads((meta_profile_dict / f'q4posprofile2_{i}').read_text()) for i in q4_meta_retain]
        q4posprofile2 = [x for x in map(mean, zip(*q4posprofile2))]
        
        q4negprofile2 = [ujson.loads((meta_profile_dict / f'q4negprofile2_{i}').read_text()) for i in q4_meta_retain]
        q4negprofile2 = [x for x in map(mean, zip(*q4negprofile2))]
        # #UJSON
        # meta_profile_dict = ujson.loads(meta_profile_dict.read_text())
        
        #Pickle
        # with open(meta_profile_dict, 'rb') as f:
        #     meta_profile_dict = pickle.load(f)

        # #Pure Python
        # temp_dict = {}
        # key = ''
        # with open(meta_profile_dict) as F:
        #     for line in F:
        #         if '#' in line[0]:
        #             key = line.strip()[1:]
        #             temp_dict[key] = []
        #         else:
        #             temp_dict[key].append([float(x) for x in line.strip().split(',')])
        # meta_profile_dict = temp_dict

    #Create Meta Plot
    # if len(meta_profile_dict) != 0:
        ax3 = plt.subplot(meta_gs[:2, 0])
        ax4 = plt.subplot(meta_gs[:2, 1])
        ax5 = plt.subplot(meta_gs[:2, 2])
        ax6 = plt.subplot(meta_gs[:2, 3])
        xvals = range(-int(largewindow),int(largewindow))
        # q1posprofile1 = meta_profile_dict['q1posprofile1']
        # q1posprofile1 = [x for i,x in enumerate(q1posprofile1) if i in q1_meta_retain]
        # q1posprofile1 = [x for x in map(mean, zip(*q1posprofile1))]
        # q1negprofile1 = meta_profile_dict['q1negprofile1']
        # q1negprofile1 = [x for i,x in enumerate(q1negprofile1) if i in q1_meta_retain]
        # q1negprofile1 = [x for x in map(mean, zip(*q1negprofile1))]
        # q1posprofile2 = meta_profile_dict['q1posprofile2']
        # q1posprofile2 = [x for i,x in enumerate(q1posprofile2) if i in q1_meta_retain]
        # q1posprofile2 = [x for x in map(mean, zip(*q1posprofile2))]
        # q1negprofile2 = meta_profile_dict['q1negprofile2']
        # q1negprofile2 = [x for i,x in enumerate(q1negprofile2) if i in q1_meta_retain]
        # q1negprofile2 = [x for x in map(mean, zip(*q1negprofile2))]
        # q2posprofile1 = meta_profile_dict['q2posprofile1']
        # q2posprofile1 = [x for i,x in enumerate(q2posprofile1) if i in q2_meta_retain]
        # q2posprofile1 = [x for x in map(mean, zip(*q2posprofile1))]
        # q2negprofile1 = meta_profile_dict['q2negprofile1']
        # q2negprofile1 = [x for i,x in enumerate(q2negprofile1) if i in q2_meta_retain]
        # q2negprofile1 = [x for x in map(mean, zip(*q2negprofile1))]
        # q2posprofile2 = meta_profile_dict['q2posprofile2']
        # q2posprofile2 = [x for i,x in enumerate(q2posprofile2) if i in q2_meta_retain]
        # q2posprofile2 = [x for x in map(mean, zip(*q2posprofile2))]
        # q2negprofile2 = meta_profile_dict['q2negprofile2']
        # q2negprofile2 = [x for i,x in enumerate(q2negprofile2) if i in q2_meta_retain]
        # q2negprofile2 = [x for x in map(mean, zip(*q2negprofile2))]
        # q3posprofile1 = meta_profile_dict['q3posprofile1']
        # q3posprofile1 = [x for i,x in enumerate(q3posprofile1) if i in q3_meta_retain]
        # q3posprofile1 = [x for x in map(mean, zip(*q3posprofile1))]
        # q3negprofile1 = meta_profile_dict['q3negprofile1']
        # q3negprofile1 = [x for i,x in enumerate(q3negprofile1) if i in q3_meta_retain]
        # q3negprofile1 = [x for x in map(mean, zip(*q3negprofile1))]
        # q3posprofile2 = meta_profile_dict['q3posprofile2']
        # q3posprofile2 = [x for i,x in enumerate(q3posprofile2) if i in q3_meta_retain]
        # q3posprofile2 = [x for x in map(mean, zip(*q3posprofile2))]
        # q3negprofile2 = meta_profile_dict['q3negprofile2']
        # q3negprofile2 = [x for i,x in enumerate(q3negprofile2) if i in q3_meta_retain]
        # q3negprofile2 = [x for x in map(mean, zip(*q3negprofile2))]
        # q4posprofile1 = meta_profile_dict['q4posprofile1']
        # q4posprofile1 = [x for i,x in enumerate(q4posprofile1) if i in q4_meta_retain]
        # q4posprofile1 = [x for x in map(mean, zip(*q4posprofile1))]
        # q4negprofile1 = meta_profile_dict['q4negprofile1']
        # q4negprofile1 = [x for i,x in enumerate(q4negprofile1) if i in q4_meta_retain]
        # q4negprofile1 = [x for x in map(mean, zip(*q4negprofile1))]
        # q4posprofile2 = meta_profile_dict['q4posprofile2']
        # q4posprofile2 = [x for i,x in enumerate(q4posprofile2) if i in q4_meta_retain]
        # q4posprofile2 = [x for x in map(mean, zip(*q4posprofile2))]
        # q4negprofile2 = meta_profile_dict['q4negprofile2']
        # q4negprofile2 = [x for i,x in enumerate(q4negprofile2) if i in q4_meta_retain]
        # q4negprofile2 = [x for x in map(mean, zip(*q4negprofile2))]

        ylim = [min(q1negprofile1+q1negprofile2+q2negprofile1+q2negprofile2
                        +q3negprofile1+q3negprofile2+q4negprofile1+q4negprofile2),
                    max(q1posprofile1+q1posprofile2+q2posprofile1+q2posprofile2
                        +q3posprofile1+q3posprofile2+q4posprofile1+q4posprofile2)]
        maxy = abs(max(ylim, key=abs))
        ylim = [-maxy, maxy]

        metaplot(q1posprofile1, q1negprofile1, q1posprofile2, q1negprofile2, ax=ax3, 
                    xvals=xvals, label1=label1, label2=label2, ylim=ylim, 
                    title='Q1 (n=' + str(len(q1_meta_retain)) + ')', 
                    largewindow=largewindow)
        ax3.legend(loc='best', fontsize='small', frameon=False)
        metaplot(q2posprofile1, q2negprofile1, q2posprofile2, q2negprofile2, ax=ax4, 
                    xvals=xvals, label1=label1, label2=label2, ylim=ylim, 
                    title='Q2 (n=' + str(len(q2_meta_retain)) + ')', 
                    largewindow=largewindow)
        metaplot(q3posprofile1, q3negprofile1, q3posprofile2, q3negprofile2, ax=ax5, 
                    xvals=xvals, label1=label1, label2=label2, ylim=ylim, 
                    title='Q3 (n=' + str(len(q3_meta_retain)) + ')', 
                    largewindow=largewindow)
        metaplot(q4posprofile1, q4negprofile1, q4posprofile2, q4negprofile2, ax=ax6, 
                    xvals=xvals, label1=label1, label2=label2, ylim=ylim, 
                    title='Q4 (n=' + str(len(q4_meta_retain)) + ')', 
                    largewindow=largewindow)

    #Initiate heatmaps
    ax7 = plt.subplot(meta_gs[2, 0])
    ax8 = plt.subplot(meta_gs[2, 1])
    ax9 = plt.subplot(meta_gs[2, 2])
    ax10 = plt.subplot(meta_gs[2, 3])
    
    bins = 100
    xlim = [int(-largewindow), int(largewindow)]
    # if len(meta_profile_dict) != 0:
    if type(meta_profile_dict) == pathlib.PosixPath:
        heatmap(q1_distances, ax=ax7, xlim=xlim, bins=bins, largewindow=largewindow)
        heatmap(q2_distances, ax=ax8, xlim=xlim, bins=bins, largewindow=largewindow)
        heatmap(q3_distances, ax=ax9, xlim=xlim, bins=bins, largewindow=largewindow)
        heatmap(q4_distances, ax=ax10, xlim=xlim, bins=bins, largewindow=largewindow)
    else:
        heatmap(q1_distances, ax=ax7, xlim=xlim, bins=bins, largewindow=largewindow, 
                    title=f'Q1 (n={len(q1_distances)})')
        heatmap(q2_distances, ax=ax8, xlim=xlim, bins=bins, largewindow=largewindow, 
                    title=f'Q2 (n={len(q2_distances)})')
        heatmap(q3_distances, ax=ax9, xlim=xlim, bins=bins, largewindow=largewindow, 
                    title=f'Q3 (n={len(q3_distances)})')
        heatmap(q4_distances, ax=ax10, xlim=xlim, bins=bins, largewindow=largewindow, 
                    title=f'Q4 (n={len(q4_distances)})')

    plt.tight_layout()
    plt.savefig(os.path.join(figuredir, motif + f'_enrichment_plot.{plot_format}'), 
                format=plot_format)#, dpi=dpi, bbox_inches='tight')
    plt.close()

    #Simulation Plot
    F = plt.figure(figsize=(9,7.75))
    ax = F.add_subplot(111)
    bins=100
    maximum = max(sim_auc)
    minimum = min(sim_auc)

    ax.hist(sim_auc,bins=bins, linewidth=0)
    width = (maximum-minimum)/100.0
    ylim_max = ax.get_ylim()[1]

    if offset != 0:
        ax.bar(auc,ylim_max,color='red',width=width*2, linewidth=0, label='Non-Corrected E-Score')[0]
        # height = non_corrected_rect.get_height()
        # ax.text(non_corrected_rect.get_x() + non_corrected_rect.get_width()/2., 
        #         1.05*height, 'Non-Corected AUC', ha='center', va='bottom')

    ax.bar(auc-offset,ylim_max,color='green',width=width*2, linewidth=0, label='Corrected E-Score')[0]

    ax.legend(frameon=False)

    # height = corrected_rect.get_height()

    # ax.text(corrected_rect.get_x() + corrected_rect.get_width()/2., 
    #         1.05*height, 'Corrected AUC', ha='center', va='bottom')
    ax.set_xlim([min(minimum,auc, auc-offset)-(width*40), max(maximum,auc,auc-offset)+(width*40)])

    ax.set_ylim([0,(1.05*ylim_max)+5])
    ax.tick_params(axis='y', which='both', left=False, right=False, 
                    labelleft=True)

    ax.tick_params(axis='x', which='both', bottom=False, top=False, 
                    labelbottom=True)

    ax.set_title('Distribution of Simulated E-Scores', fontsize=14)
    ax.set_ylabel('Number of Simulations', fontsize=14)
    ax.set_xlabel('E-Score', fontsize=14)

    plt.tight_layout()
    F.savefig(os.path.join(figuredir, motif + f'_simulation_plot.{plot_format}'), 
                format=plot_format)#, dpi=dpi, bbox_inches='tight')
    plt.close(F)
    return

#==============================================================================
@force_gc
def plot_global_MA(results, p_cutoff=None, title=None, xlabel=None, 
                    ylabel=None, x_index=None, y_index=None, c_index=None,
                    p_index=None, savepath=None, ylimits=None, 
                    plot_format=None):
    '''Plot an MA plot. Note: x-values will be transformed to log10

    Parameters
    ----------
    results : list of lists
        contains calculated enrichment scores for all TFs of interest specified
        by the user
    '''
    F = plt.figure(figsize=(7,6))
    ax = plt.subplot(111)
    clean_results = [i for i in results if i[y_index] == i[y_index] and i[x_index] == i[x_index]]
    ylist = [i[y_index] for i in clean_results]
    # import sys
    # print(clean_results, file=sys.stderr)
    xlist = [math.log(i[x_index], 10) if i[x_index] != 0 else 0 for i in clean_results]
    if c_index is not None:
        clist = [i[c_index] for i in clean_results]
        clist = [y-c for y,c in zip(ylist,clist)]
        max_c = abs(max([x for x in clist if x == x], key=abs))
        # import sys
        # print("MA c-list:", clist, file=sys.stderr)
        # print("MA max_c:", max_c, file=sys.stderr)
        scatter = ax.scatter(xlist, ylist, edgecolor='', c=clist, s=50, cmap='viridis',
                                vmax=max_c, vmin=-max_c)
    else:
        scatter = ax.scatter(xlist, ylist, edgecolor='', color='navy', s=50)

    if p_index is not None:
        plist = [i[p_index] for i in clean_results]
        sigx = [x for x,p in zip(xlist,plist) if p<p_cutoff]
        sigy = [y for y,p in zip(ylist,plist) if p<p_cutoff]
        if c_index is not None:
            sigc = [scatter.to_rgba(c) for c,p in zip(clist,plist) if p<p_cutoff]
            ax.scatter(sigx, sigy, c=sigc, #marker='x', 
                        edgecolor='r',  linewidth=2, s=50)
            if p_cutoff < -3:
                legend = ax.scatter([1], [0], color='white', edgecolor='r',  s=50, 
                                    label=f'p < 1e{int(p_cutoff*np.log10(np.e))}')
            else:
                legend = ax.scatter([1], [0], color='white', edgecolor='r',  s=50, 
                                    label=f'p < {str("%.3g" % np.e**p_cutoff)}')
        else:
            ax.scatter(sigx, sigy, color='red', edgecolor='',  s=50)
            if p_cutoff < -3:
                legend = ax.scatter([1], [0], color='red', edgecolor='',  s=50, 
                                    label=f'p < 1e{int(p_cutoff*np.log10(np.e))}')
            else:
                legend = ax.scatter([1], [0], color='red', edgecolor='',  s=50, 
                                    label=f'p < {str("%.3g" % np.e**p_cutoff)}')
            
        ax.legend(loc='best', frameon=False)
        legend.remove()

    # ylist = [i[1] for i in results]
    # xlist = [math.log(i[2], 10) if i[2] != 0 else 0 for i in results]
    # plist = [i[-1] for i in results]

    # sigx = [x for x, p in zip(xlist, plist) if p < p_cutoff]
    # sigy = [y for y, p in zip(ylist, plist) if p < p_cutoff]

    # ax.scatter(xlist, ylist, color='navy', edgecolor='', s=50)
    # ax.scatter(sigx, sigy, color='red', edgecolor='', s=50)

    ax.set_title(title, fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_xlabel(xlabel, fontsize=14)
    ax.tick_params(axis='y', which='both', left=True, right=False, 
                    labelleft=True)
    ax.tick_params(axis='x', which='both', bottom=True, top=False, 
                    labelbottom=True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if ylimits is not None:
        ax.set_ylim(ylimits)

    plt.tight_layout()
    F.savefig(str(savepath), format=plot_format)#, dpi=dpi, bbox_inches='tight')
    plt.close()
    return

#==============================================================================
@force_gc
def plot_global_volcano(results, p_cutoff=None, title=None, xlabel=None, 
                        ylabel=None, savepath=None, plot_format=None):
    '''This function plots graphs that are displayed on the main results.html 
        filethat correspond to results relating to all analyzed TFs.

    Parameters
    ----------
    results : list of lists
        contains calculated enrichment scores for all TFs of interest specified
        by the user
    '''
    xlist = [i[1] for i in results]
    ylist = [-i[-1] if i[-1] != 0 else 0 for i in results]
    plist = [i[-1] for i in results]

    sigx = [x for x, p in zip(xlist, plist) if p < p_cutoff]
    sigy = [p for x, p in zip(ylist, plist) if p < p_cutoff]
    sigy = [-p if p != 0 else 0 for p in sigy]

    F = plt.figure(figsize=(7,6))
    ax = plt.subplot(111)
    ax.scatter(xlist, ylist, color='navy', edgecolor='', s=50)
    ax.scatter(sigx, sigy, color='red', edgecolor='', s=50)
    ax.set_title(title, fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_xlabel(xlabel, fontsize=14)
    ax.axhline(-p_cutoff, 10, linestyle='--', color='black')
    ax.tick_params(axis='y', which='both', left=True, right=False, 
                    labelleft=True)

    ax.tick_params(axis='x', which='both', bottom=True, top=False, 
                    labelbottom=True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    F.savefig(str(savepath), format=plot_format)#, dpi=dpi, bbox_inches='tight')
    plt.close()
    return

#==============================================================================
@force_gc
def plot_global_z_v(results, p_cutoff=None, title=None, xlabel=None, 
                        ylabel=None, savepath=None, dpi=100, x_index=None,
                        y_index=None, p_index=None, s_index=None, c_index=None, 
                        plot_format=None):
    ylist = [i[y_index] for i in results]
    xlist = [math.log(i[x_index], 10) if i[x_index] != 0 else 0 for i in results]
    plist = [i[p_index] for i in results]

    if s_index is not None:
        slist = [i[s_index] for i in results]
        print(slist)
        slist = [((s/max([sigs for sigs, p in zip(slist,plist) if p < p_cutoff]))**2)*100 if p < p_cutoff else 50 for s,p in zip(slist,plist) ]
        print(slist)
    else:
        slist = [50 for i in results]
    
    # print([tup for tup in zip(slist,xlist,ylist,[i[0] for i in results])])
    
    clist = np.zeros((len(results),4))
    clist = np.array([[1.0, 0, 0, 0] if p < p_cutoff else [0, 0, 0.4, 0] for p in plist])
    if c_index is not None:
        tempclist = [i[c_index] for i in results]
        max_val = max(tempclist)
        tempclist = [math.sqrt(c/max_val) for c in tempclist]
        clist[:,3] = tempclist
    else:
        clist[:,3] = 1.0


    F = plt.figure(figsize=(7,6))
    ax = plt.subplot(111)
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')
        ax.scatter(xlist, ylist, color=clist, edgecolor='', s=slist)
    ax.set_title(title, fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_xlabel(xlabel, fontsize=14)

    plt.tight_layout()
    F.savefig(str(savepath), format=plot_format)#, dpi=dpi, bbox_inches='tight')
    plt.close()
    return

#==============================================================================
@force_gc
def plot_global_gc(results, p_cutoff=None, title=None, xlabel=None, 
                        ylabel=None, savepath=None, dpi=100, x_index=None,
                        y_index=None, p_index=None, c_index=None, 
                        linear_regression=None, ylimits=None, plot_format=None):

    F = plt.figure(figsize=(7,6))
    ax = plt.subplot(111)
    #Remove Nans from results
    clean_results = [i for i in results if i[y_index] == i[y_index] and i[x_index] == i[x_index]]
    ylist = [i[y_index] for i in clean_results]
    xlist = [i[x_index] for i in clean_results]
    clist = [i[c_index] for i in clean_results]
    clist = [c-y for y,c in zip(ylist,clist)]
    max_c = abs(max([x for x in clist if x == x], key=abs))
    scatter = ax.scatter(xlist, ylist, edgecolor='', c=clist, s=50, cmap='viridis',
                            vmax=max_c, vmin=-max_c)
    cbar = plt.colorbar(scatter)
    cbar.set_label('E-Score Correction', rotation=270, labelpad=20)
    cbar.outline.set_visible(False)

    if p_index is not None:
        plist = [i[p_index] for i in clean_results]
        sigx = [x for x,p in zip(xlist,plist) if p<p_cutoff]
        sigy = [y for y,p in zip(ylist,plist) if p<p_cutoff]
        sigc = [scatter.to_rgba(c) for c,p in zip(clist,plist) if p<p_cutoff]
        ax.scatter(sigx, sigy, c=sigc, #marker='x', 
                        edgecolor='r', linewidth=2, s=50)
    
    if linear_regression is not None:
        slope, intercept, r_value, p_value, _ = linear_regression
        s = ("y = (" + str("%.2g" % slope) + ")x + " + str("%.2g" % intercept)
            + "\nR$^2$ = " + str("%.2g" % r_value**2)
            + "\np-val = " + str("%.2g" % p_value))
        ax.plot([0,1],[intercept, slope+intercept], color='r', alpha=0.5, label=s, 
                linewidth=5)

        ax.legend(loc='best', frameon=False)

    ax.axhline(0, linestyle='--', alpha=0.5, linewidth=2, c='k')
    ax.set_xlim([0,1])
    ax.set_title(title, fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_xlabel(xlabel, fontsize=14)
    ax.tick_params(axis='y', which='both', left=True, right=False, 
                    labelleft=True)

    ax.tick_params(axis='x', which='both', bottom=True, top=False, 
                    labelbottom=True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if ylimits is not None:
        ax.set_ylim(ylimits)

    plt.tight_layout()
    F.savefig(str(savepath), format=plot_format)#, dpi=dpi, bbox_inches='tight')
    plt.close()
    return

#==============================================================================
@force_gc
def meme_logo(motif_file, motif_ID, figuredir, plot_format=None):
    '''Runs meme2images that creates logo images
    '''
    meme2images_command = ['meme2images', '-rc', '-eps', '-motif', motif_ID, 
                            motif_file, figuredir]
    motif_ID = motif_ID.replace('.', '_')
    imagemagick_command = ['convert', figuredir / ('logo'+motif_ID+'.eps'), 
                            figuredir / (f'logo{motif_ID}.png')]
    imagemagick_rc_command = ['convert', figuredir / ('logo_rc'+motif_ID+'.eps'), 
                            figuredir / (f'logo_rc{motif_ID}.png')]
    try:
        subprocess.check_output(meme2images_command, stderr=subprocess.PIPE)
        subprocess.check_output(imagemagick_command, stderr=subprocess.PIPE)
        subprocess.check_output(imagemagick_rc_command, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode(), flush=True, file=sys.stdout)
    return

#==============================================================================
@force_gc
def metaplot(posprofile1, negprofile1, posprofile2, negprofile2, ax=None, 
                xvals=None, label1=None, label2=None, title=None, ylim=None, 
                largewindow=None):
    if len(posprofile1) != 0:
        ax.plot(xvals,posprofile1,color='#7570b3',label=label1)
    else:
        ax.plot(xvals, [0 for x in xvals], color='#7570b3',label=label1)
    if len(negprofile1) != 0:
        ax.plot(xvals,negprofile1,color='#7570b3')
    else:
        ax.plot(xvals,[0 for x in xvals],color='#7570b3')
    if len(posprofile2) != 0:
        ax.plot(xvals,posprofile2,color='#d76127',label=label2)
    else:
        ax.plot(xvals,[0 for x in xvals],color='#d76127',label=label2)
    if len(negprofile2) != 0:
        ax.plot(xvals,negprofile2,color='#d76127')
    else:
        ax.plot(xvals,[0 for x in xvals],color='#d76127')
    ax.set_title(title,fontsize=14)
    ax.tick_params(axis='y', which='both', left=False, right=False, 
                    labelleft=True)
    ax.tick_params(axis='x', which='both', bottom=False, top=False, 
                    labelbottom=False)
    ax.set_ylabel('Reads per Millions Mapped',fontsize=10)
    ax.set_ylim(ylim)
    ax.set_xlim([-int(largewindow), int(largewindow)])
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3g'))
    return

#==============================================================================
@force_gc
def heatmap(distances, ax=None, xlim=None, bins=None, title=None, 
            largewindow=None):
    counts, edges = np.histogram(distances, bins=bins)
    edges = (edges[1:]+edges[:-1])/2.0
    norm    = matplotlib.colors.Normalize(vmin=min(counts), 
                                            vmax=max(counts))
    cmap    = cm.YlOrRd
    m       = cm.ScalarMappable(norm=norm, cmap=cmap)
    colors  = [m.to_rgba(c) for c in counts] 
    
    ax.bar(edges,np.ones((len(edges),)), color=colors, 
                width=(edges[-1]-edges[0])/len(edges), edgecolor=colors)
    ax.set_ylim([0,1])
    ax.set_xlim(xlim)
    ax.tick_params(axis='y', which='both', left=False, right=False, 
                    labelleft=False) 
    ax.tick_params(axis='x', which='both', bottom=True, top=False, 
                    labelbottom=True)
    ax.set_xlabel('Motif Distance to Center (bp)')
    # ax.set_xticks(ax.get_xticks())
    # locs = [item.get_text() for item in ax.get_xticklabels()]
    # import sys
    # print(locs, file=sys.stderr)
    # locs = [str(float(x)/1000.0) for x in locs]
    # ax.set_xticklabels(locs)
    ax.set_xticks([-int(largewindow), 0, int(largewindow)])
    if title is not None:
        ax.set_title(title, fontsize=14)
    return

#==============================================================================
@force_gc
def lineplot(title=None, ax=None, xvals=None, yvals=None, xlimits=None):
    #This is the enrichment score plot (i.e. line plot)
    ax.plot(xvals,yvals,color='green')
    ax.plot([0, 1],[0, 1], '--', alpha=0.75)
    ax.set_title(title, fontsize=14)
    ax.set_ylabel('Enrichment Score (ES)', fontsize=10)
    ax.tick_params(axis='y', which='both', left=True, right=False, 
                    labelleft=True)
    ax.tick_params(axis='x', which='both', bottom=False, top=False, 
                    labelbottom=False)
    ax.set_ylim([0,1])
    ax.set_xlim(xlimits)
    return

#==============================================================================
@force_gc
def scatterplot(ax=None, xvals=None, yvals=None, xlimits=None, largewindow=None, 
                xlabel=False):
    #This is the barplot right below the enrichment score line plot
    ax.scatter(xvals, yvals, edgecolor="", color="black", 
                    s=10, alpha=0.25)
    ax.tick_params(axis='y', which='both', left=True, right=False, 
                    labelleft=True) 
    if xlabel:
        ax.tick_params(axis='x', which='both', bottom=False, top=False, 
                    labelbottom=True)
        ax.set_xlabel('Relative Rank (n='+str(len(xvals))+')', fontsize=14)
    else:
        ax.tick_params(axis='x', which='both', bottom=False, top=False, 
                        labelbottom=False)

    ax.set_yticks([-int(largewindow),0,int(largewindow)])
    ax.set_yticklabels([str(largewindow/-1000.0), '0', str(largewindow/1000.0)])
    ax.set_xlim(xlimits)
    ax.set_ylim([-int(largewindow),int(largewindow)])
    ax.set_ylabel('Distance (kb)', fontsize=10)
    return

#==============================================================================
@force_gc
def barplot(ax=None, xvals=None, colorarray=None, xlimits=None):
    norm    = matplotlib.colors.Normalize(vmin=min(colorarray), 
                                            vmax=max(colorarray))
    cmap    = cm.Greys
    m       = cm.ScalarMappable(norm=norm, cmap=cmap)
    colors  = [m.to_rgba(c) for c in colorarray] 
    ax.bar(xvals, [1 for x in xvals], width=1.0/len(xvals), edgecolor="", color=colors)
    ax.tick_params(axis='y', which='both', left=False, right=False, 
                    labelleft=False) 
    ax.tick_params(axis='x', which='both', bottom=False, top=False, 
                    labelbottom=False)
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
    ax.set_ylabel('Score', fontsize=10)
    return

#==============================================================================
@force_gc
def fillplot(ax=None, xvals=None, yvals=None, xlimits=None, ylimits=None):
    #This is the rank metric fill plot
    posvals = [(x,y) for x,y in zip(xvals,yvals) if y > 0]
    negvals = [(x,y) for x,y in zip(xvals,yvals) if y < 0]
    ax.fill_between([x for x,_ in posvals], 0, [y for _,y in posvals],facecolor='#d76127',edgecolor="")
    ax.fill_between([x for x,_ in negvals], 0, [y for _,y in negvals],facecolor='#7570b3',edgecolor="")
    ax.tick_params(axis='y', which='both', left=True, right=False, 
                    labelleft=True)
    ax.tick_params(axis='x', which='both', bottom=False, top=False, 
                    labelbottom=True)
    ylim = math.fabs(max([x for x in yvals if -500 < x < 500],key=abs))
    ax.yaxis.set_ticks([int(-ylim), 0, int(ylim)])
    ax.set_xlim(xlimits)
    ax.set_xlabel('Relative Rank (n='+str(len(xvals))+')', fontsize=14)
    ax.set_ylabel('Rank Metric', fontsize=10)
    # try:
    #     ax2.axvline(len(updistancehist)+1,color='green',alpha=0.25)
    # except ValueError:
    #     pass
    # try:
    #     ax2.axvline(len(xvals) - len(downdistancehist), color='purple', 
    #                 alpha=0.25)
    # except ValueError:
    #     pass
    return

#==============================================================================
@force_gc
def plot_deseq_MA(deseq_file=None, label1=None, label2=None, figuredir=None, 
                    dpi=100, basemean_cut=0, plot_format=None):
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

    c = np.linspace(0, 1, len(x))

    #Creates an MA-Plot of the region expression
    F = plt.figure(figsize=(7,6))
    ax = plt.subplot(111)
    plt.scatter(x=x,y=y,c=c,edgecolor='', cmap="RdYlGn")
    ax.set_title("DE-Seq MA-Plot",fontsize=14)
    ax.set_ylabel("Log2 Fold-Change ("+label2+"/"+label1+")",fontsize=14)
    ax.set_xlabel("Log10 Average Expression",fontsize=14)
    ax.tick_params(axis='y', which='both', left=True, right=False, 
                    labelleft=True)
    ax.tick_params(axis='x', which='both', bottom=True, top=False, 
                    labelbottom=True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if basemean_cut != 0:
        ax.axvline(math.log(basemean_cut, 10), linestyle='--', c='k', alpha=0.5)
    cbar = plt.colorbar()
    cbar.ax.invert_yaxis() 
    cbar.set_ticks([0, 0.25, 0.5, 0.75, 1], ['0', '0.25', '0.5', '0.75', '1'])
    cbar.set_label('Relative Rank (n=' + str(len(x)) + ')', rotation=270, 
                        labelpad=20)
    cbar.outline.set_visible(False)
    plt.tight_layout()
    F.savefig(os.path.join(figuredir, f'DESEQ_MA_Plot.{plot_format}'), format=plot_format)#, dpi=dpi)
                # bbox_inches='tight')
    plt.close()
    return

if __name__ == "__main__":
    results_file = '/Users/joru1876/Google_Drive/Colorado_University/Jonathan/TFEA_outputs/Allen2014/v5_outputs/20190620_DMSO_Nutlin_fimohits/results.txt'
    results = []
    with open(results_file) as F:
        for line in F:
            if '#' not in line[0]:
                linelist = [float(l) if 'H11MO' not in l else l for l in line.strip('\n').split('\t')]
                results.append(linelist)
    plot_global_MA(results, p_cutoff=0.001, title='MA Plot', 
                    c_index=4,
                    x_index=2,
                    y_index=1,
                    p_index=-1,
                    ylimits=[-0.5,0.5],
                    xlabel="GC-content", 
                    ylabel="E-Score", 
                    savepath='/Users/joru1876/Google_Drive/Colorado_University/Jonathan/TFEA_outputs/Allen2014/v5_outputs/20190620_DMSO_Nutlin_fimohits/newMA_plot.png', 
                    plot_format='png')

    plot_global_gc(results, p_cutoff=0.001, title='GC-Plot', xlabel='G-C', 
                        ylabel='y-axis', 
                        savepath='/Users/joru1876/Google_Drive/Colorado_University/Jonathan/TFEA_outputs/Allen2014/v5_outputs/20190620_DMSO_Nutlin_fimohits/newGC_plot.png', 
                        plot_format='png', 
                        x_index=3,
                        y_index=1, 
                        c_index=4,
                        p_index=-1,
                        ylimits=[-0.5,0.5])
