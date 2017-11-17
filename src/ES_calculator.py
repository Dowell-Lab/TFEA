__author__ = 'Jonathan Rubin'

import matplotlib
matplotlib.use('Agg')
import math
import os
import numpy as np
import matplotlib.pyplot as plt
import pybedtools as py

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
                    print len(vals)
                    print i
                if pval > 0.1 and len(pvalcut2) == 0:
                    pvalcut2.append(i)
                    print "pval < 0.01"
                    print len(vals)
                    print i
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
    H = 1500
    ES = list()
    Eval = 0
    vals = list()
    ind = list()
    total = 0
    distance_sum = 0
    with open(ranked_center_distance_file) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            val = float(line[-1])
            r = int(line[4])
            total += 1
            if val > H:
                pass
            else:
                vals.append(val)
                ind.append(r)
                distance_sum += H-val

    for i in range(total):
        if i in ind:
            Eval += float((H-vals[ind.index(i)]))/float(distance_sum)
            ES.append(Eval)
        else:
            Eval += -1.0/(total-len(ind))
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

    return max(ES,key=abs)

