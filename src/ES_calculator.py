__author__ = 'Jonathan Rubin'

import matplotlib
import sys
matplotlib.use('Agg')
import math
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from random import shuffle
from scipy.stats import norm
import time
from config import FDRCUTOFF

def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])
    
    return newdir

def run(MOTIF_FILE,ranked_center_distance_file,figuredir,filedir):
    #Initiate some variables
    H = 1500.0
    ES = list()
    Eval = 0.0
    distances = list()
    ind = list()
    negatives = 0.0
    distance_sum = 0.0

    #First parse file containing motif distance and region rank. Also count total negatives to be used later
    with open(ranked_center_distance_file) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            distance = float(line[-1])
            rank = int(line[4])
            if distance > H:
                distances.append(-1)
                ind.append(rank)
                negatives += 1.0
            else:
                distances.append(H-distance)
                ind.append(rank)
                distance_sum += H-distance

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
    a = time.time()
    simES = simulate(H,distances,distance_sum,neg)
    print "Simulation done in: ", time.time()-a, "s"

    #NES is the actualES divided by the mean ES of all permutations with the same sign as actualES
    #p-value is caluclated empirically (i.e. (# of simulated ES scores larger than actualES)/(rest of simulated ES scores))
    if actualES < 0:
        simESsubset = [x for x in simES if x < 0]
        mu = np.mean(simESsubset)
        sigma = np.std(simESsubset)
        NES = -(actualES/mu)
    else:
        simESsubset = [x for x in simES if x > 0]
        mu = np.mean(simESsubset)
        NES = actualES/mu

    # #Now calculate an NES for each simulated ES
    # simNES = list()
    # for singleES in simES:
    #     if singleES < 0:
    #         mu = np.mean(negativesubset)
    #         simNES.append(-(singleES/mu))
    #     else:
    #         mu = np.mean(positivesubset)
    #         simNES.append(singleES/mu)

    #This section calculates the theoretical p-value based on the mean and standard deviation of the 1000 simulations
    #The p-value is then the CDF where x = actualES. Test is two tailed, hence: min(p,1-p)
    mu = np.mean(simES)
    sigma = np.std(simES)
    p = norm.cdf(actualES,mu,sigma)
    p = min(p,1-p)

    #Plot results for significant hits while list of simulated ES scores is in memory
    if p < FDRCUTOFF:
        #Smooth hits to plot later
        hits = [1.0*(x/H) if x!=neg else 0 for x in distances]
        hitlength = len(hits)
        width = hitlength/1000
        # YlOrRd = plt.get_cmap('YlOrRd')
        # hist,_=np.histogram(hits,bins=width)
        # newhits = [0]*len(hits)
        # for i in range(0,len(hits),window):
        #     if sum(hits[i:i+window]) > 1:
        #         newhits[i+window/2] = 1

        newhits = list()
        for i in range(0,hitlength,width):
            newhits.append(float(sum(hits[i:i+width])))

        newhitmax = float(max(newhits))
        newhitlength = len(newhits)

        alphas = [(x/newhitmax) for x in newhits]
        rgba_colors = np.zeros((newhitlength,4))
        rgba_colors[:,0] = 0.0
        rgba_colors[:,1] = 0.0
        rgba_colors[:,2] = 0.0
        rgba_colors[:,3] = alphas

        F = plt.figure(figsize=(15,6))
        gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
        ax0 = plt.subplot(gs[0])
        ax0.plot(range(1,len(ES)+1),ES,color='green')
        ax1 = plt.subplot(gs[1])
        print range(1,hitlength+1,width)[:10], ([1]*newhitlength)[:10], rgba_colors[:10]
        ax1.bar(range(1,hitlength+1,width), [1]*newhitlength ,color=rgba_colors,edgecolor = "none")
        ax0.tick_params(
            axis='y',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            left='off',        # ticks along the bottom edge are off
            right='off',       # ticks along the top edge are off
            labelleft='on')
        ax0.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom='off',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            labelbottom='off')
        ax1.tick_params(
            axis='y',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            left='off',        # ticks along the bottom edge are off
            right='off',       # ticks along the top edge are off
            labelleft='off')
        ax1.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom='off',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            labelbottom='on')
        ax0.set_title('Enrichment Plot: '+ MOTIF_FILE,fontsize=14)
        ylims = ax0.get_ylim()
        ymax = math.fabs(max(ylims,key=abs))
        ax0.set_ylim([-ymax,ymax])
        ax0.set_xlim([1,hitlength])
        ax1.set_xlim([1,hitlength])
        ax1.set_xlabel('Rank in Ordered Dataset', fontsize=14)
        ax1.set_ylabel('Hits', fontsize=14)
        ax0.set_ylabel('Enrichment Score (ES)', fontsize=14)
        plt.savefig(figuredir + MOTIF_FILE + '_enrichment_plot.svg',bbox_inches='tight')
        plt.close()

        F = plt.figure()
        ax2 = plt.subplot(111)
        maximum = max(simES)
        minimum = min(simES)
        ax2.hist(simES,bins=100)
        width = (maximum-minimum)/100.0
        ax2.bar(actualES,ax2.get_ylim()[1]-1.0,color='red',width=width*2)
        ax2.set_xlim([min(minimum,actualES)-(width*4),max(maximum,actualES)+(width*4)])
        ax2.tick_params(
            axis='y',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            left='off',        # ticks along the bottom edge are off
            right='off',       # ticks along the top edge are off
            labelleft='on')
        ax2.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom='off',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            labelbottom='on')
        plt.title('Distribution of Simulated Enrichment Scores',fontsize=14)
        ax2.set_ylabel('Number of Simulations',fontsize=14)
        ax2.set_xlabel('Enrichment Score (ES)',fontsize=14)
        plt.savefig(figuredir + MOTIF_FILE + '_simulation_plot.svg')
        plt.close()

    return [MOTIF_FILE,actualES,NES,p]

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

def FDR(TFresults,figuredir):
    #This function iterates through the results and calculates an FDR for each TF motif. Also creates a moustache plot.
    FDRlist = list()
    newNESlist = list()
    LABELS = list()
    sigx = list()
    sigy = list()
    for i in range(len(TFresults)):
        NES = TFresults[i][2]
        PVAL = TFresults[i][3]

        #Using classical FDR calculation ((pvalue*(# hypotheses tested))/rank of p-value)
        FDR = (PVAL*len(TFresults))/float(i+1.0)
        TFresults[i].append(FDR)
        FDRlist.append(FDR)
        newNESlist.append(NES)
        if FDR < FDRCUTOFF:
            sigx.append(NES)
            sigy.append(FDR)

    F = plt.figure()
    plt.scatter(newNESlist,FDRlist,color='black',edgecolor='')
    plt.scatter(sigx,sigy,color='red',edgecolor='')
    plt.title("TFEA Results Moustache Plot",fontsize=14)
    plt.xlabel("Normalized Enrichment Score (NES)",fontsize=14)
    plt.ylabel("False Discovery Rate (FDR)",fontsize=14)
    limit = math.fabs(max(newNESlist,key=abs))
    plt.xlim([-limit,limit])
    plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        left='off',        # ticks along the bottom edge are off
        right='off',       # ticks along the top edge are off
        labelleft='on')
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='on')
    plt.savefig(figuredir + 'TFEA_Results_Moustache_Plot.svg')
    return TFresults


