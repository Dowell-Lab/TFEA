__author__ = 'Jonathan Rubin'

import matplotlib
matplotlib.use('Agg')
import sys
import math
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.stats import norm
import time
from multiprocessing import Pool
# import HTSeq as hts
from config import *
import motif_distance
import meta_eRNA

def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])
    
    return newdir

# def pool(MOTIF_HITS,ranked_center_distance_file,figuredir,logos):
#     global ranked_center_distance_file
#     ranked_center_distance_file = ranked_center_distance_file
#     global figuredir
#     figuredir = figuredir
#     global logos
#     logos=logos
#     p = Pool(32)
#     TFresults = p.map(ES_calculator.run,os.listdir(config.MOTIF_HITS))
#     return TFresults

def run(args):
    MOTIF_FILE,ranked_center_distance_file,ranked_center_file,figuredir,logos = args
    ranked_center_distance_file = motif_distance.run(ranked_center_file,MOTIF_HITS+MOTIF_FILE)
    #Initiate some variables
    H = 1500.0
    h = 150.0
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

    #First parse file containing motif distance and region rank. Also count total negatives to be used later
    with open(ranked_center_distance_file) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            distance = float(line[-1])
            pval = float(line[3])
            rank = int(line[5])
            fc = float(line[4])
            if 0 <= distance <= h:
                positives += 1.0
                value = math.exp(-distance)
                distances.append(value)
                ind.append(rank)
                distance_sum += value
                scatterx.append(rank)
                scattery.append(distance)
                if fc > 1:
                    if pval < PVALCUTOFF:
                        updistancehist.append(distance)
                        upsigscatterx.append(rank)
                    try:
                        logpval.append(-math.log(pval,10))
                    except ValueError:
                        logpval.append(500.0)
                    
                else:
                    if pval < PVALCUTOFF:
                        downdistancehist.append(distance)
                        downsigscatterx.append(rank)
                    try:
                        logpval.append(math.log(pval,10))
                    except ValueError:
                        logpval.append(-500.0)
                    

            elif h < distance <= H:
                distances.append(-1)
                ind.append(rank)
                negatives += 1.0
                scatterx.append(rank)
                scattery.append(distance)
                if fc > 1:
                    if pval < PVALCUTOFF:
                        updistancehist.append(distance)
                        upsigscatterx.append(rank)
                    try:
                        logpval.append(-math.log(pval,10))
                    except ValueError:
                        logpval.append(500.0)
                    
                else:
                    if pval < PVALCUTOFF:
                        downdistancehist.append(distance)
                        downsigscatterx.append(rank)
                    try:
                        logpval.append(math.log(pval,10))
                    except ValueError:
                        logpval.append(-500.0)
                    
                    
    
    if len(distances) == 0:
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
    simES = simulate(H,distances,distance_sum,neg)

    #NES is the actualES divided by the mean ES of all permutations with the same sign as actualES
    #p-value is caluclated empirically (i.e. (# of simulated ES scores larger than actualES)/(rest of simulated ES scores))
    if actualES < 0:
        simESsubset = [x for x in simES if x < 0]
        mu = np.mean(simESsubset)
        NES = -(actualES/mu)
        sigma = np.std(simESsubset)
        p = norm.cdf(actualES,mu,sigma)
    else:
        simESsubset = [x for x in simES if x > 0]
        mu = np.mean(simESsubset)
        NES = actualES/mu
        sigma = np.std(simESsubset)
        p = 1-norm.cdf(actualES,mu,sigma)

    #Plot results for significant hits while list of simulated ES scores is in memory
    # if p < FDRCUTOFF or SINGLEMOTIF != False:
    #For human:
    # os.system("scp " + logos + MOTIF_FILE.split('.bed')[0].split('HO_')[1] + "_direct.png " + figuredir)
    # os.system("scp " + logos + MOTIF_FILE.split('.bed')[0].split('HO_')[1] + "_revcomp.png " + figuredir)

    #For mouse:
    os.system("scp " + logos + MOTIF_FILE.split('.bed')[0] + "_direct.png " + figuredir)
    os.system("scp " + logos + MOTIF_FILE.split('.bed')[0] + "_revcomp.png " + figuredir)

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
    ax1.axhline(h, color='red',alpha=0.25)
    if DRAWPVALCUTOFF != False:
        ax1.scatter(upsigscatterx,updistancehist,edgecolor="",color="green",s=10,alpha=0.5)
        ax1.scatter(downsigscatterx,downdistancehist,edgecolor="",color="blue",s=10,alpha=0.5)
        # ax1.axvline(PVALCUTOFF,linestyle='dashed')
        # ax1.text(PVALCUTOFF,H+H/10,str(PVALCUTOFF),ha='center',va='bottom')
    ax1.tick_params(axis='y', which='both', left='off', right='off', labelleft='on')
    ax1.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
    ax1.set_xlim(limits)
    ax1.set_ylim([0,H])
    # ax1.yaxis.set_ticks([0,H])
    plt.yticks([0,H],['0',str(float(H)/1000.0)])
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
        ax2.axvline(len(upsigscatterx)+1,color='green',alpha=0.25)
    except ValueError:
        pass
    try:
        ax2.axvline(len(xvals) - len(downsigscatterx),color='purple',alpha=0.25)
    except ValueError:
        pass

    plt.savefig(figuredir + MOTIF_FILE.split('.bed')[0] + '_enrichment_plot.png',bbox_inches='tight')
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
    plt.savefig(figuredir + MOTIF_FILE.split('.bed')[0] + '_simulation_plot.png',bbox_inches='tight')
    plt.cla()

    #Plots the distribution of motif distances with a red line at h
    F = plt.figure(figsize=(7,6))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])
    ax0 = plt.subplot(gs[0])
    binwidth = H/100.0
    ax0.hist(updistancehist,bins=np.arange(0,int(H)+binwidth,binwidth),color='green')
    ax0.set_title('Distribution of Motif Distance for: fc > 1, pval < ' + str(PVALCUTOFF),fontsize=14)
    ax0.axvline(h,color='red',alpha=0.5)
    ax0.tick_params(axis='y', which='both', left='off', right='off', labelleft='on')
    ax0.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
    ax0.set_xlim([0,H])
    ax0.set_xlabel('Distance (bp)',fontsize=14)
    ax0.set_ylabel('Hits',fontsize=14)
    ax1 = plt.subplot(gs[1])
    ax1.hist(downdistancehist,bins=np.arange(0,int(H)+binwidth,binwidth),color='purple')
    ax1.axvline(h,color='red',alpha=0.5)
    ax1.set_title('Distribution of Motif Distance for: fc < 1, pval < ' + str(PVALCUTOFF),fontsize=14)
    ax1.tick_params(axis='y', which='both', left='off', right='off', labelleft='on')
    ax1.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
    ax1.set_xlim([0,H])
    ax1.set_ylabel('Hits',fontsize=14)
    ax1.set_xlabel('Distance (bp)',fontsize=14)
    plt.tight_layout()
    plt.savefig(figuredir + MOTIF_FILE.split('.bed')[0] + '_distance_distribution.png',bbox_inches='tight')
    plt.cla()

    #Plots a meta_eRNA
    posprofile1, negprofile1, posprofile2, negprofile2 = meta_eRNA.run2(ranked_center_distance_file)
    F = plt.figure(figsize=(15.5,6))
    ax0 = plt.subplot(111)
    xvals = range(-int(H),int(H))
    line1, = ax0.plot(xvals,posprofile1,color='blue',label=LABEL1)
    ax0.plot(xvals,negprofile1,color='blue')
    line2, = ax0.plot(xvals,posprofile2,color='red',label=LABEL2)
    ax0.plot(xvals,negprofile2,color='red')
    ax0.legend(loc=1)
    ax0.axvline(0,color='black',alpha=0.25)
    ax0.axhline(0,color='black',linestyle='dashed')
    ax0.set_title('Meta Plot of eRNA Signal',fontsize=14)
    ax0.tick_params(axis='y', which='both', left='off', right='off', labelleft='on')
    ax0.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
    ax0.set_ylabel('Normalized Read Coverage',fontsize=14)
    ax0.set_xlabel('Distance to eRNA Origin (bp)')
    plt.savefig(figuredir + MOTIF_FILE.split('.bed')[0] + '_meta_eRNA.png',bbox_inches='tight')
    plt.cla()

    os.system("rm " + ranked_center_distance_file)

    # return [MOTIF_FILE.split('.bed')[0],actualES,NES,p,(simNESmu,simNESsigma)]
    return [MOTIF_FILE.split('.bed')[0],actualES,NES,p,positives,negatives]

def simulate(H,distances,distance_sum,neg,N=1000):
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

def FDR(TFresults,NESlist,figuredir):
    #This function iterates through the results and calculates an FDR for each TF motif. Also creates a moustache plot.
    FDRlist = list()
    newNESlist = list()
    sigx = list()
    sigy = list()
    pvals = list()
    totals = list()
    sigtotals = list()

    for i in range(len(TFresults)):
        NES = TFresults[i][2]
        PVAL = TFresults[i][3]
        POS = TFresults[i][4]
        NEG = TFresults[i][5]
        pvals.append(PVAL)
        total = POS+NEG
        #Using classical FDR calculation ((pvalue*(# hypotheses tested))/rank of p-value)
        FDR = ((float(i)+1.0)/float(len(TFresults)))*FDRCUTOFF
        # FDR = (PVAL*len(TFresults))/float(i+1.0)
        TFresults[i].append(FDR)
        FDRlist.append(FDR)
        newNESlist.append(NES)
        totals.append(math.log(float(total)))

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


    #Creates a scatter plot of NES vs. number of motif hits within H
    F = plt.figure(figsize=(7,6))
    ax = plt.subplot(111)
    ax.scatter(newNESlist,totals,color='blue',edgecolor='')
    ax.scatter(sigx,sigtotals,color='red',edgecolor='')
    ax.set_title("TFEA NES MA-Plot",fontsize=14)
    ax.set_xlabel("Normalized Enrichment Score (NES)",fontsize=14)
    ax.set_ylabel("Number of Motif Hits within " + str(int(H)) + "bp (log10)",fontsize=14)
    ax.tick_params(axis='y', which='both', left='off', right='off', labelleft='on')
    ax.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
    plt.savefig(figuredir + 'TFEA_NES_MA_Plot.png',bbox_inches='tight')
    plt.cla()

    #Creates a moustache plot of the global FDRs vs. NESs
    F = plt.figure(figsize=(7,6))
    ax = plt.subplot(111)
    ax.scatter(newNESlist,FDRlist,color='black',edgecolor='')
    ax.scatter(sigx,sigy,color='red',edgecolor='')
    ax.set_title("TFEA Moustache Plot",fontsize=14)
    ax.set_xlabel("Normalized Enrichment Score (NES)",fontsize=14)
    ax.set_ylabel("False Discovery Rate (FDR)",fontsize=14)
    limit = math.fabs(max(newNESlist,key=abs))
    ax.set_xlim([-limit,limit])
    ax.tick_params(axis='y', which='both', left='off', right='off', labelleft='on')
    ax.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
    plt.savefig(figuredir + 'TFEA_Results_Moustache_Plot.png',bbox_inches='tight')
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
    plt.savefig(figuredir + 'TFEA_Pval_Histogram.png',bbox_inches='tight')
    plt.cla()

    return TFresults


