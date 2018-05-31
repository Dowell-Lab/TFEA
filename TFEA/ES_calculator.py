__author__ = 'Jonathan Rubin'

import matplotlib
matplotlib.use('Agg')
import sys
import math
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import gridspec
from scipy.stats import norm as normal
from scipy.stats import poisson
import time
from multiprocessing import Pool
import config
import motif_distance
import meta_eRNA
import create_html
import combine_bed
import meme
import main

def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])
    
    return newdir

def run(args):
    MOTIF_FILE,millions_mapped = args
    ranked_center_distance_file = motif_distance.run(config.RANKED_CENTER_FILE,config.MOTIF_HITS+MOTIF_FILE)
    if config.FIMO:
        ranked_fullregions_file = combine_bed.get_regions()
        ranked_fasta_file = combine_bed.getfasta(ranked_fullregions_file)
        background_file = combine_bed.get_bgfile(ranked_fasta_file)
        meme.fimo(background_file,config.SINGLEMOTIF.strip('.bed'),ranked_fasta_file)
        meme.meme2images(config.SINGLEMOTIF)

    #Initiate some variables
    ES = list()
    Eval = 0.0
    distances = list()
    ind = list()
    negatives = 0.0
    positives = 0.0
    distance_sum = 0.0
    scatterx = list()
    upsigscatterx = list()
    downsigscatterx = list()
    scattery = list()
    upsigscattery = list()
    downsigscattery = list()
    logpval = list()
    updistancehist = list()
    downdistancehist = list()
    middledistancehist = list()

    ##RFS addition
    ranks = list()
    updistancehist2 = list()
    downdistancehist2 = list()
    middledistancehist2 = list()
    upregions = list()
    downregions = list()
    middleregions = list()

    #First parse file containing motif distance and region rank. Also count total negatives to be used later
    with open(ranked_center_distance_file) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            distance = float(line[-1])
            pval = float(line[3])
            rank = int(line[5])
            fc = float(line[4])

            ##get distances based on quantile ranks
            ##added by RFS
            #5/22/18: Added negative distances, added lists for regions to plot meta_eRNA later
            ranks.append(rank) 

            if rank < np.percentile(ranks,25):
                updistancehist2.append(distance)
                upregions.append(line[:3])
            elif rank > np.percentile(ranks,75):
                downdistancehist2.append(distance)
                downregions.append(line[:3])
            else:
                middledistancehist2.append(distance)
                middleregions.append(line[:3])

            #Get absolute value of distance to be used by later plots
            distance = math.fabs(distance)

            #If a motif is within a SMALLWINDOW, then it is counted as a 'positive' hit
            if 0 <= distance <= config.SMALLWINDOW:
                positives += 1.0
                value = math.exp(-distance)
                distances.append(value)
                ind.append(rank)
                distance_sum += value
                scatterx.append(rank)
                scattery.append(distance)
                if fc > 1:
                    if pval < config.PVALCUTOFF:
                    	upsigscatterx.append(rank)
                    try:
                        logpval.append(-math.log(pval,10))
                    except ValueError:
                        logpval.append(500.0)
                    
                else:
                    if pval < config.PVALCUTOFF:
                    	downsigscatterx.append(rank)
                    try:
                        logpval.append(math.log(pval,10))
                    except ValueError:
                        logpval.append(-500.0)
                    
            #If a motif is not within a SMALLWINDOW, it is counted as a 'negative' hit
            elif distance > config.SMALLWINDOW:
                distances.append(-1)
                ind.append(rank)
                negatives += 1.0
                scatterx.append(rank)
                scattery.append(distance)
                if fc > 1:
                    if pval < config.PVALCUTOFF:
                        upsigscatterx.append(rank)
                    try:
                        logpval.append(-math.log(pval,10))
                    except ValueError:
                        logpval.append(500.0)
                else:
                    if pval < config.PVALCUTOFF:
                        downsigscatterx.append(rank)
                    try:
                        logpval.append(math.log(pval,10))
                    except ValueError:
                        logpval.append(-500.0)
                    
                    
    
    if len(distances) == 0 or positives < 10 or negatives < 10:
        return "no hits"        

    #actualES calculation:
    try:
        neg = -1.0/negatives
    except:
        neg = -1.0

    #1. Replace negatives values in distances with -1.0/negatives
    distances = [neg if x==-1 else x for x in distances]

    #2. Reorder list of motif distances based on rank
    distances = [x for _,x in sorted(zip(ind,distances))]

    #3. Go through distances and add appropriately to cumulative sum (ES) list
    for distance in distances:
        if distance != neg:
            Eval += distance/distance_sum
            ES.append(Eval)
        else:
            Eval += distance
            ES.append(Eval)

    #4. The enrichment score is the maximum deviation from 0
    actualES = max(ES,key=abs)

    #To get NES, first simulate 1000 permuations of region ranks
    simES = simulate(distances,distance_sum,neg)

    #NES is the actualES divided by the mean ES of all permutations with the same sign as actualES
    #p-value is caluclated with the theoretical normal distribution
    if actualES < 0:
        simESsubset = [x for x in simES if x < 0]
        mu = np.mean(simESsubset)
        NES = -(actualES/mu)
        sigma = np.std(simESsubset)
        p = normal.cdf(actualES,mu,sigma)
    else:
        simESsubset = [x for x in simES if x > 0]
        mu = np.mean(simESsubset)
        NES = actualES/mu
        sigma = np.std(simESsubset)
        p = 1-normal.cdf(actualES,mu,sigma)

    #Plot results for significant hits while list of simulated ES scores is in memory
    if 'HO_' in MOTIF_FILE:
        os.system("scp " + config.LOGOS + MOTIF_FILE.split('.bed')[0].split('HO_')[1] + "_direct.png " + config.FIGUREDIR)
        os.system("scp " + config.LOGOS + MOTIF_FILE.split('.bed')[0].split('HO_')[1] + "_revcomp.png " + config.FIGUREDIR)
    else:
        os.system("scp " + config.LOGOS + MOTIF_FILE.split('.bed')[0] + "_direct.png " + config.FIGUREDIR)
        os.system("scp " + config.LOGOS + MOTIF_FILE.split('.bed')[0] + "_revcomp.png " + config.FIGUREDIR)

    logpval = [x for _,x in sorted(zip(scatterx,logpval))]
    scattery = [x for _,x in sorted(zip(scatterx,scattery))]

    #Plots the enrichment plot which contains three subplots:
    #   1. Typical ES vs. region rank (GSEA-style)
    #   2. A scatter plot with distance to motif on y-axis and rank on x-axis
    #   3. A plot similar to GSEAs 'phenotype label' plot to show where the fc crosses 1. On the y-axis is log10(pval) (positive if
    #       fc > 1, else negative) on the x-axis is regions ranked
    F = plt.figure(figsize=(15.5,6))
    xvals = range(1,len(ES)+1)
    limits = [1,len(ES)]
    gs = gridspec.GridSpec(3, 1, height_ratios=[3, 1, 1])
    ax0 = plt.subplot(gs[0])
    ax0.plot(xvals,ES,color='green')
    ax0.axhline(0, color='black',linestyle='dashed')
    ax0.set_title('Enrichment Plot: '+ MOTIF_FILE.split('.bed')[0],fontsize=14)
    ax0.set_ylabel('Enrichment Score (ES)', fontsize=10)
    ax0.tick_params(axis='y', which='both', left='on', right='off', labelleft='on')
    ax0.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
    ylims = ax0.get_ylim()
    ymax = math.fabs(max(ylims,key=abs))
    ax0.set_ylim([-ymax,ymax])
    ax0.set_xlim(limits)
    ax1 = plt.subplot(gs[1])
    ax1.scatter(xvals,scattery,edgecolor="",color="black",s=10,alpha=0.25)
    ax1.axhline(config.SMALLWINDOW, color='red',alpha=0.25)
    if config.DRAWPVALCUTOFF != False:
        ax1.scatter(upsigscatterx,updistancehist2,edgecolor="",color="green",s=10,alpha=0.5)
        ax1.scatter(downsigscatterx,downdistancehist2,edgecolor="",color="blue",s=10,alpha=0.5)
    ax1.tick_params(axis='y', which='both', left='off', right='off', labelleft='on')
    ax1.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
    ax1.set_xlim(limits)
    ax1.set_ylim([0,config.LARGEWINDOW])
    plt.yticks([0,config.LARGEWINDOW],['0',str(float(config.LARGEWINDOW)/1000.0)])
    ax1.set_ylabel('Distance (kb)', fontsize=10)
    ax2 = plt.subplot(gs[2])
    ax2.fill_between(xvals,0,logpval,facecolor='grey',edgecolor="")
    ax2.tick_params(axis='y', which='both', left='on', right='off', labelleft='on')
    ax2.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
    ylim = math.fabs(max([x for x in logpval if -500 < x < 500],key=abs))
    ax2.set_ylim([-ylim,ylim])
    ax2.yaxis.set_ticks([int(-ylim),0,int(ylim)])
    ax2.set_xlim(limits)
    ax2.set_xlabel('Rank in Ordered Dataset', fontsize=14)
    ax2.set_ylabel('Rank Metric',fontsize=10)
    try:
        ax2.axvline(len(updistancehist2)+1,color='green',alpha=0.25)
    except ValueError:
        pass
    try:
        ax2.axvline(len(xvals) - len(downdistancehist2), color='purple', alpha=0.25)
    except ValueError:
        pass

    plt.savefig(config.FIGUREDIR + MOTIF_FILE.split('.bed')[0] + '_enrichment_plot.png',bbox_inches='tight')
    plt.cla()

    #Plots the distribution of simulated ESs and in a red bar plots the observed ES
    F = plt.figure(figsize=(7,6))
    ax2 = plt.subplot(111)
    maximum = max(simES)
    minimum = min(simES)
    ax2.hist(simES,bins=100)
    width = (maximum-minimum)/100.0
    rect = ax2.bar(actualES,ax2.get_ylim()[1],color='red',width=width*2)[0]
    height = rect.get_height()
    ax2.text(rect.get_x() + rect.get_width()/2., 1.05*height, 'Observed ES', ha='center', va='bottom')
    ax2.set_xlim([min(minimum,actualES)-(width*40),max(maximum,actualES)+(width*40)])
    ax2.set_ylim([0,(1.05*height)+5])
    ax2.tick_params(axis='y', which='both', left='off', right='off', labelleft='on')
    ax2.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
    plt.title('Distribution of Simulated Enrichment Scores',fontsize=14)
    ax2.set_ylabel('Number of Simulations',fontsize=14)
    ax2.set_xlabel('Enrichment Score (ES)',fontsize=14)
    plt.savefig(config.FIGUREDIR + MOTIF_FILE.split('.bed')[0] + '_simulation_plot.png',bbox_inches='tight')
    plt.cla()

    #5/29/18: Commented this out in favor of heatmaps with meta_eRNA plotss
    #Plots the distribution of motif distances with a red line at h
    # F = plt.figure(figsize=(6.5,6))
    # gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 1])
    # ax0 = plt.subplot(gs[0])
    # binwidth = LARGEWINDOW/100.0
    # ax0.hist(updistancehist2,bins=np.arange(0,int(LARGEWINDOW)+binwidth,binwidth),color='green')
    # ax0.set_title('Distribution of Motif Distance for: fc > 1',fontsize=14)
    # ax0.axvline(h,color='red',alpha=0.5)
    # ax0.tick_params(axis='y', which='both', left='off', right='off', labelleft='on')
    # ax0.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
    # ax0.set_xlim([0,H])
    # ax0.set_xlabel('Distance (bp)',fontsize=14)
    # ax0.set_ylabel('Hits',fontsize=14)
    # ax1 = plt.subplot(gs[2])
    # ax1.hist(downdistancehist2,bins=np.arange(0,int(LARGEWINDOW)+binwidth,binwidth),color='purple')
    # ax1.axvline(h,color='red',alpha=0.5)
    # ax1.set_title('Distribution of Motif Distance for: fc < 1',fontsize=14)
    # ax1.tick_params(axis='y', which='both', left='off', right='off', labelleft='on')
    # ax1.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
    # ax1.set_xlim([0,LARGEWINDOW])
    # ax1.set_ylabel('Hits',fontsize=14)
    # ax1.set_xlabel('Distance (bp)',fontsize=14)
    # ax2 = plt.subplot(gs[1])
    # ax2.hist(middledistancehist2,bins=np.arange(0,int(H)+binwidth,binwidth),color='blue')
    # ax2.set_title('Distribution of Motif Distance for: middle',fontsize=14)
    # ax2.axvline(h,color='red',alpha=0.5)
    # ax2.tick_params(axis='y', which='both', left='off', right='off', labelleft='on')
    # ax2.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
    # ax2.set_xlim([0,H])
    # ax2.set_xlabel('Distance (bp)',fontsize=14)
    # ax2.set_ylabel('Hits',fontsize=14)
    # plt.tight_layout()
    # plt.savefig(config.FIGUREDIR + MOTIF_FILE.split('.bed')[0] + '_distance_distribution.png',bbox_inches='tight')
    # plt.cla()

    #Plots a meta_eRNA
    #5/22/18: Major changes, now three meta_eRNA plots are created with associated heatmaps underneath
    posprofile_up1, negprofile_up1, posprofile_up2, negprofile_up2 = meta_eRNA.run3(upregions,millions_mapped)
    posprofile_middle1, negprofile_middle1, posprofile_middle2, negprofile_middle2 = meta_eRNA.run3(middleregions,millions_mapped)
    posprofile_down1, negprofile_down1, posprofile_down2, negprofile_down2 = meta_eRNA.run3(downregions,millions_mapped)
    F = plt.figure(figsize=(15.5,6))
    gs = gridspec.GridSpec(2, 3, height_ratios=[3, 1])
    xvals = range(-int(config.LARGEWINDOW),int(config.LARGEWINDOW))

    #Plot meta_eRNA for up regions
    ax0 = plt.subplot(gs[0,0])
    line1, = ax0.plot(xvals,posprofile_up1,color='blue',label=config.LABEL1)
    ax0.plot(xvals,negprofile_up1,color='blue')
    line2, = ax0.plot(xvals,posprofile_up2,color='red',label=config.LABEL2)
    ax0.plot(xvals,negprofile_up2,color='red')
    ax0.axvline(0,color='black',alpha=0.25)
    ax0.axhline(0,color='black',linestyle='dashed')
    ax0.tick_params(axis='y', which='both', left='off', right='off', labelleft='on')
    ax0.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
    ax0.legend(loc=1)
    ax0.set_title('Meta Plot of eRNA Signal (upper quartile)',fontsize=14,color='green')
    ax0.set_ylabel('Normalized Read Coverage',fontsize=14)

    #Plot heatmap underneath meta_eRNA for up regions
    ax1 = plt.subplot(gs[1,0])
    yvals,edges = np.histogram(updistancehist2,bins=len(xvals))
    edges = (edges[1:]+edges[:-1])/2. #bascially np.histogram gives you a list of start and stops of the bins, so we are just taking the center of the bins, so it matches the number of counts
    maximum = float(max(yvals))
    yvals = [float(y)/maximum for y in yvals]
    norm = matplotlib.colors.Normalize(vmin=min(yvals), vmax=max(yvals))
    cmap = cm.YlOrRd
    m = cm.ScalarMappable(norm=norm, cmap=cmap)
    colors = [m.to_rgba(c) for c in yvals] 
    ax1.bar(edges,[1 for i in range(len(xvals))], color=colors, width=(edges[-1]-edges[0])/len(edges), edgecolor=colors)
    ax1.set_xlim(xvals)
    ax1.tick_params(axis='y', which='both', left='off', right='off', labelleft='off')
    ax1.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
    ax1.set_xlabel('Distance to eRNA Origin (bp)')

    #Plot meta_eRNA for middle regions
    ax3 = plt.subplot(gs[0,1])
    line1, = ax3.plot(xvals,posprofile_middle1,color='blue',label=config.LABEL1)
    ax3.plot(xvals,negprofile_middle1,color='blue')
    line2, = ax3.plot(xvals,posprofile_middle2,color='red',label=config.LABEL2)
    ax3.plot(xvals,negprofile_middle2,color='red')
    ax3.axvline(0,color='black',alpha=0.25)
    ax3.axhline(0,color='black',linestyle='dashed')
    ax3.tick_params(axis='y', which='both', left='off', right='off', labelleft='on')
    ax3.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
    # ax3.legend(loc=1)
    ax3.set_title('Meta Plot of eRNA Signal (middle quartiles)',fontsize=14,color='blue')
    ax3.set_ylabel('Normalized Read Coverage',fontsize=14)

    #Plot heatmap underneath meta_eRNA for middle regions
    ax4 = plt.subplot(gs[1,1])
    yvals,edges = np.histogram(middledistancehist2,bins=len(xvals))
    edges = (edges[1:]+edges[:-1])/2. #bascially np.histogram gives you a list of start and stops of the bins, so we are just taking the center of the bins, so it matches the number of counts
    maximum = float(max(yvals))
    yvals = [float(y)/maximum for y in yvals]
    norm = matplotlib.colors.Normalize(vmin=min(yvals), vmax=max(yvals))
    cmap = cm.YlOrRd
    m = cm.ScalarMappable(norm=norm, cmap=cmap)
    colors = [m.to_rgba(c) for c in yvals] 
    ax4.bar(edges,[1 for i in range(len(xvals))], color=colors, width=(edges[-1]-edges[0])/len(edges), edgecolor=colors)
    ax4.set_xlim(xvals)
    ax4.tick_params(axis='y', which='both', left='off', right='off', labelleft='off')
    ax4.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
    ax4.set_xlabel('Distance to eRNA Origin (bp)')

    #Plot meta_eRNA for down regions
    ax5 = plt.subplot(gs[0,2])
    line1, = ax5.plot(xvals,posprofile_down1,color='blue',label=config.LABEL1)
    ax5.plot(xvals,negprofile_down1,color='blue')
    line2, = ax5.plot(xvals,posprofile_down2,color='red',label=config.LABEL2)
    ax5.plot(xvals,negprofile_down2,color='red')
    ax5.axvline(0,color='black',alpha=0.25)
    ax5.axhline(0,color='black',linestyle='dashed')
    ax5.tick_params(axis='y', which='both', left='off', right='off', labelleft='on')
    ax5.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
    # ax5.legend(loc=1)
    ax5.set_title('Meta Plot of eRNA Signal (lower quartiles)',fontsize=14,color='purple')
    ax5.set_ylabel('Normalized Read Coverage',fontsize=14)
    ax5.set_xlabel('Distance to eRNA Origin (bp)')

    #Plot heatmap underneath meta_eRNA for down regions
    ax6 = plt.subplot(gs[1,2])
    yvals,edges = np.histogram(downdistancehist2,bins=len(xvals))
    edges = (edges[1:]+edges[:-1])/2. #bascially np.histogram gives you a list of start and stops of the bins, so we are just taking the center of the bins, so it matches the number of counts
    maximum = float(max(yvals))
    yvals = [float(y)/maximum for y in yvals]
    norm = matplotlib.colors.Normalize(vmin=min(yvals), vmax=max(yvals))
    cmap = cm.YlOrRd
    m = cm.ScalarMappable(norm=norm, cmap=cmap)
    colors = [m.to_rgba(c) for c in yvals] 
    ax6.bar(edges,[1 for i in range(len(xvals))], color=colors, width=(edges[-1]-edges[0])/len(edges), edgecolor=colors)
    ax6.set_xlim(xvals)
    ax6.tick_params(axis='y', which='both', left='off', right='off', labelleft='off')
    ax6.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
    ax6.set_xlabel('Distance to eRNA Origin (bp)')


    plt.savefig(config.FIGUREDIR + MOTIF_FILE.split('.bed')[0] + '_meta_eRNA.png',bbox_inches='tight')
    plt.cla()

    #5/22/18: commented this out for now
    # os.system("rm " + ranked_center_distance_file)

    return [MOTIF_FILE.split('.bed')[0],actualES,NES,p,positives,negatives]

def simulate(distances,distance_sum,neg,N=1000):
    #Simulate 1000 permuations of region ranks
    simES = list()
    for i in range(N):
        Eval = 0.0
        maximum = 0.0
        minimum = 0.0
        #Here we actually shuffle the regions
        np.random.shuffle(distances)

        #Then we calculate an ES just like before
        for distance in distances:
            if distance != neg:
                Eval += distance/distance_sum
                if Eval > maximum:
                    maximum = Eval
            else:
                Eval += distance
                if Eval < minimum:
                    minimum = Eval
        simES.append(max(maximum,minimum,key=abs))

    return simES

def FDR(TFresults):
    #This function iterates through the results and calculates an FDR for each TF motif. Also creates a moustache plot.
    FDRlist = list()
    NESlist = list()
    sigx = list()
    sigy = list()
    pvals = list()
    ##pvalsNA = [1 if str(x) == 'nan' else x for x in pvals]
    totals = list()
    sigtotals = list()
    positives = list()

    TFresults = [x for x in TFresults if x != "no hits"]
    TFresults = sorted(TFresults, key=lambda x: x[3])
    for i in range(len(TFresults)):
        NES = TFresults[i][2]
        PVAL = TFresults[i][3]
        POS = TFresults[i][4]
        NEG = TFresults[i][5]
        pvals.append(PVAL)
        total = POS+NEG
        #Using classical FDR calculation ((pvalue*(# hypotheses tested))/rank of p-value)
        FDR = ((float(i)+1.0)/float(len(TFresults)))*config.FDRCUTOFF
        # FDR = (PVAL*len(TFresults))/float(i+1.0)
        TFresults[i].append(FDR)
        FDRlist.append(FDR)
        NESlist.append(NES)
        totals.append(math.log(float(total)))
        positives.append(math.log(float(POS)))

        # mu,sigma = TFresults[i][4]
        # if NES > 0:
        #     F = 1-norm.cdf(NES,mu,sigma)
        #     NESsubset = [x for x in NESlist if x > 0]
        #     N = len(NESsubset)
        #     Nmu = np.mean(NESsubset)
        #     Nsigma = np.std(NESsubset)
        #     D = 1-norm.cdf(NES,Nmu,Nsigma)
        #     q = (F*N)/D
        # else:
        #     F = norm.cdf(NES,mu,sigma)
        #     NESsubset = [x for x in NESlist if x < 0]
        #     N = len(NESsubset)
        #     Nmu = np.mean(NESsubset)
        #     Nsigma = np.std(NESsubset)
        #     D = norm.cdf(NES,Nmu,Nsigma)
        #     q = (F*N)/D

        if PVAL < FDR:
            sigx.append(NES)
            sigy.append(FDR)
            sigtotals.append(math.log(float(total)))

    ##pvalsNA = [1 if str(x) == 'nan' else x for x in pvals]
        
    create_html.createTFtext(TFresults,config.FIGUREDIR)

    #Creates a scatter plot of NES vs. number of motif hits within H
    F = plt.figure(figsize=(7,6))
    ax = plt.subplot(111)
    #5/22/18: Using positive values only instead of total
    ax.scatter(NESlist,positives,color='blue',edgecolor='')
    ax.scatter(sigx,sigtotals,color='red',edgecolor='')
    ax.set_title("TFEA NES MA-Plot",fontsize=14)
    ax.set_xlabel("Normalized Enrichment Score (NES)",fontsize=14)
    ax.set_ylabel("Number of Motif Hits within " + str(int(H)) + "bp (log10)",fontsize=14)
    ax.tick_params(axis='y', which='both', left='off', right='off', labelleft='on')
    ax.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
    plt.savefig(config.FIGUREDIR + 'TFEA_NES_MA_Plot.png',bbox_inches='tight')
    plt.cla()

    #Creates a moustache plot of the global FDRs vs. NESs
    F = plt.figure(figsize=(7,6))
    ax = plt.subplot(111)
    ax.scatter(NESlist,FDRlist,color='black',edgecolor='')
    ax.scatter(sigx,sigy,color='red',edgecolor='')
    ax.set_title("TFEA Moustache Plot",fontsize=14)
    ax.set_xlabel("Normalized Enrichment Score (NES)",fontsize=14)
    ax.set_ylabel("False Discovery Rate (FDR)",fontsize=14)
    limit = math.fabs(max(newNESlist,key=abs))
    ax.set_xlim([-limit,limit])
    ax.tick_params(axis='y', which='both', left='off', right='off', labelleft='on')
    ax.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
    plt.savefig(config.FIGUREDIR + 'TFEA_Results_Moustache_Plot.png',bbox_inches='tight')
    plt.cla()

    #Creates a histogram of p-values
    F = plt.figure(figsize=(7,6))
    ax = plt.subplot(111)
    binwidth = 1.0/100.0
    ax.hist(pvals,bins=np.arange(0,1.0+binwidth,binwidth),color='green')
    ax.set_title("TFEA P-value Histogram",fontsize=14)
    ax.set_xlabel("P-value",fontsize=14)
    ax.set_ylabel("Count",fontsize=14)
    ax.tick_params(axis='y', which='both', left='off', right='off', labelleft='on')
    ax.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
    plt.savefig(config.FIGUREDIR + 'TFEA_Pval_Histogram.png',bbox_inches='tight')
    plt.cla()

    return TFresults


