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

def run(MOTIF_FILE,ranked_center_distance_file,figuredir,filedir,total_hits):
    #Initiate some variables
    H = 1500.0
    ES = list()
    Eval = 0.0
    distances = list()
    ind = list()
    total = 0.0
    negatives = 0.0
    distance_sum = 0.0

    #First parse file containing motif distance and region rank. Also count total negatives to be used later
    with open(ranked_center_distance_file) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            distance = float(line[-1])
            rank = int(line[4])
            total += 1
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
    simES = simulate(H,distances,distance_sum,total,neg)
    print "Simulation done in: ", time.time()-a, "s"

    #NES is the actualES divided by the mean ES of all permutations with the same sign as actualES
    #p-value is caluclated empirically (i.e. (# of simulated ES scores larger than actualES)/(rest of simulated ES scores))
    if actualES < 0:
        simESsubset = [x for x in simES if x < 0]
        mu = np.mean(simESsubset)
        NES = -(actualES/mu)
    else:
        simESsubset = [x for x in simES if x > 0]
        mu = np.mean(simESsubset)
        NES = actualES/mu

    #Now calculate an NES for each simulated ES
    simNES = list()
    for singleES in simES:
        if singleES < 0:
            simESsubset = [x for x in simES if x < 0]
            mu = np.mean(simESsubset)
            simNES.append(-(singleES/mu))
        else:
            simESsubset = [x for x in simES if x > 0]
            mu = np.mean(simESsubset)
            simNES.append(singleES/mu)

    #This section calculates the theoretical p-value based on the mean and standard deviation of the 1000 simulations
    #The p-value is then the CDF where x = actualES. Test is two tailed, hence: min(p,1-p)
    mu = np.mean(simES)
    sigma = np.std(simES)
    p = norm.cdf(actualES,mu,sigma)
    p = min(p,1-p)

    #Plot results for significant hits while list of simulated ES scores is in memory
    if p < FDRCUTOFF:
        #Smooth hits to plot later
        hits = [1 if x!=neg else 0 for x in distances]
        window = 100
        newhits = [0]*len(hits)
        for i in range(0,len(hits),window):
            if sum(hits[i:i+window]) > window/2:
                newhits[i+window/2] = 1

        F = plt.figure()
        gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
        ax0 = plt.subplot(gs[0])
        ax0.plot(range(1,len(ES)+1),ES,color='green')
        ax1 = plt.subplot(gs[1])
        ax1.bar(range(1,len(ES)+1), newhits ,color='black')
        # plt.plot(range(1,len(ES)+1),ES)
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
        ax1.set_xlabel('Rank in Ordered Dataset', fontsize=14)
        ax1.set_ylabel('Hits', fontsize=14)
        ax0.set_ylabel('Enrichment Score (ES)', fontsize=14)
        plt.savefig(figuredir + MOTIF_FILE + '_enrichment_plot.svg')
        plt.close()

        F = plt.figure()
        ax2 = plt.subplot(111)
        ax2.hist(simES,bins=100)
        width = (max(simES)-min(simES))/100.0
        ax2.bar(actualES,ax2.get_ylim()[1]-1.0,color='red',width=width)
        ax2.set_xlim([min(min(simES),actualES)+width,max(max(simES),actualES)+width])
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

def simulate(H,distances,distance_sum,total,neg,N=1000):
    #Simulate 1000 permuations of region ranks
    simES = list()
    for i in range(N):
        Eval = 0.0
        ES = list()
        #Here we actually shuffle the regions
        np.random.shuffle(distances)

        #Then we calculate an ES just like before
        for distance in distances:
            if distance != neg:
                Eval += distance/distance_sum
                ES.append(Eval)
            else:
                Eval += distance
                ES.append(Eval)
        simES.append(max(ES,key=abs))

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


