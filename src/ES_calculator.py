__author__ = 'Jonathan Rubin'

import matplotlib
import sys
matplotlib.use('Agg')
import math
import os
import numpy as np
import matplotlib.pyplot as plt
import pybedtools as py
from random import shuffle
from scipy.stats import norm
import time

def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])
    
    return newdir

def run_deseq(ranked_file,figuredir):
    #Home directory
    homedir = os.path.dirname(os.path.realpath(__file__))

    #Directory where all temp files will be stored
    filedir = parent_dir(homedir) + '/files/'

    ranked_files = [filedir + 'ranked_file_up.bed',filedir + 'ranked_file_down.bed']
    for ranked_file in ranked_files:
        H = 1500
        ES = list()
        vals = list()
        redvals = list()
        redind = list()
        ind = list()
        i = 0
        D = py.BedTool(filedir + 'SRR1105737-1_bidir_predictions.bed')
        N = py.BedTool(filedir + 'SRR1105739-1_bidir_predictions.bed')
        R = py.BedTool(ranked_file)
        R.intersect((N-D),c=True).saveas(ranked_file)
        pvalcut1 = list()
        pvalcut2 = list()
        with open(ranked_file) as F:
            for line in F:
                line = line.strip('\n').split('\t')
                val = float(line[4])
                pval = float(line[3])
                if pval > 0.05 and len(pvalcut1) == 0:
                    pvalcut1.append(i)
                    print "pval < 0.05"
                    print 'positives:',len(vals)
                    print 'total',i
                if pval > 0.1 and len(pvalcut2) == 0:
                    pvalcut2.append(i)
                    print "pval < 0.01"
                    print 'positives',len(vals)
                    print 'total',i
                if val > H:
                    pass
                else:
                    vals.append(math.fabs(val-H)/H)
                    ind.append(i)
                    if line[-2] != "0":
                        redvals.append(math.fabs(val-H)/H)
                        redind.append(i)
                i += 1

        print len(vals)
        F = plt.figure(figsize=(30,5))
        # cbar = plt.colorbar(colors)
        plt.scatter(ind,vals,edgecolor="")
        plt.scatter(redind,redvals,c='r',edgecolor="")
        plt.axvline(pvalcut1,linestyle='dashed')
        plt.axvline(pvalcut2,linestyle='dashed')
        plt.tick_params(
            axis='y',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            left='off',      # ticks along the bottom edge are off
            right='off',         # ticks along the top edge are off
            labelleft='off')
        plt.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom='off',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            labelbottom='off')
        plt.xlim((ind[0],ind[-1]))
        plt.savefig(figuredir + ranked_file.split('/')[-1] + 'test.png')

        plt.close()

def run(ranked_center_distance_file,figuredir,filedir,total_hits):
    H = 1500.0
    ES = list()
    Eval = 0.0
    distances = list()
    ind = list()
    total = 0.0
    negatives = 0.0
    distance_sum = 0.0
    with open(ranked_center_distance_file) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            distance = float(line[-1])
            rank = int(line[4])
            total += 1
            if distance > H:
                distances.append(-1)
                ind.append(r)
                negatives += 1.0
            else:
                distances.append(H-distance)
                ind.append(r)
                distance_sum += H-distance

    neg = -1.0/negatives
    distances = [neg if x==-1 else x for x in distances]
    distances = [distances[i] for i in ind]
    print distances[:10]
    for distance in distances:
        if distance != neg:
            Eval += distance/distance_sum
            ES.append(Eval)
        else:
            Eval += distance
            ES.append(Eval)
    
    # F = plt.figure(figsize=(30,5))
    # # cbar = plt.colorbar(colors)
    # plt.scatter(ind,vals,edgecolor="")
    # plt.tick_params(
    #     axis='y',          # changes apply to the x-axis
    #     which='both',      # both major and minor ticks are affected
    #     left='off',        # ticks along the bottom edge are off
    #     right='off',       # ticks along the top edge are off
    #     labelleft='off')
    # plt.tick_params(
    #     axis='x',          # changes apply to the x-axis
    #     which='both',      # both major and minor ticks are affected
    #     bottom='off',      # ticks along the bottom edge are off
    #     top='off',         # ticks along the top edge are off
    #     labelbottom='off')
    # plt.xlim((ind[0],ind[-1]))
    # plt.savefig(figuredir + ranked_file.split('/')[-1] + '.png')

    # plt.close()
    actualES = max(ES,key=abs)
    a = time.time()
    simES = simulate(H,ind,distances,distance_sum,total,negatives)
    print "Simulation done in: ", time.time()-a, "s"
    mu = np.mean(simES)
    sigma = np.std(simES)
    NES = actualES/mu
    p = norm.cdf(actualES,mu,sigma)
    p = min(p,1-p)
    return actualES,NES,p

def simulate(H,ind,distances,distance_sum,total,negatives,N=1000):
    simES = list()
    for i in range(1,N):
        b = "Simulations done: " + str(i)
        print "\r",
        sys.stdout.write("\033[K")
        print "\r",b,
        sys.stdout.flush()
        Eval = 0.0
        ES = list()
        neg = -1.0/negatives
        distances = [neg if x==-1 else x for x in distances]
        np.random.shuffle(ind)
        distances = [distances[i] for i in ind]
        for distance in distances:
            if distance != neg:
                Eval += distance/distance_sum
                ES.append(Eval)
            else:
                Eval += distance
                ES.append(Eval)
        simES.append(max(ES,key=abs))
    print "\n"
    return simES

