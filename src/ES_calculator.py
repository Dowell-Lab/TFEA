__author__ = 'Jonathan Rubin'

import matplotlib
import sys
matplotlib.use('Agg')
import math
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
# import pybedtools as py
from random import shuffle
from scipy.stats import norm
import time
from config import FDRCUTOFF

def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])
    
    return newdir

# def run_deseq(ranked_file,figuredir):
#     #Home directory
#     homedir = os.path.dirname(os.path.realpath(__file__))

#     #Directory where all temp files will be stored
#     filedir = parent_dir(homedir) + '/files/'

#     ranked_files = [filedir + 'ranked_file_up.bed',filedir + 'ranked_file_down.bed']
#     for ranked_file in ranked_files:
#         H = 1500
#         ES = list()
#         vals = list()
#         redvals = list()
#         redind = list()
#         ind = list()
#         i = 0
#         D = py.BedTool(filedir + 'SRR1105737-1_bidir_predictions.bed')
#         N = py.BedTool(filedir + 'SRR1105739-1_bidir_predictions.bed')
#         R = py.BedTool(ranked_file)
#         R.intersect((N-D),c=True).saveas(ranked_file)
#         pvalcut1 = list()
#         pvalcut2 = list()
#         with open(ranked_file) as F:
#             for line in F:
#                 line = line.strip('\n').split('\t')
#                 val = float(line[4])
#                 pval = float(line[3])
#                 if pval > 0.05 and len(pvalcut1) == 0:
#                     pvalcut1.append(i)
#                     print "pval < 0.05"
#                     print 'positives:',len(vals)
#                     print 'total',i
#                 if pval > 0.1 and len(pvalcut2) == 0:
#                     pvalcut2.append(i)
#                     print "pval < 0.01"
#                     print 'positives',len(vals)
#                     print 'total',i
#                 if val > H:
#                     pass
#                 else:
#                     vals.append(math.fabs(val-H)/H)
#                     ind.append(i)
#                     if line[-2] != "0":
#                         redvals.append(math.fabs(val-H)/H)
#                         redind.append(i)
#                 i += 1

#         print len(vals)
#         F = plt.figure(figsize=(30,5))
#         # cbar = plt.colorbar(colors)
#         plt.scatter(ind,vals,edgecolor="")
#         plt.scatter(redind,redvals,c='r',edgecolor="")
#         plt.axvline(pvalcut1,linestyle='dashed')
#         plt.axvline(pvalcut2,linestyle='dashed')
#         plt.tick_params(
#             axis='y',          # changes apply to the x-axis
#             which='both',      # both major and minor ticks are affected
#             left='off',      # ticks along the bottom edge are off
#             right='off',         # ticks along the top edge are off
#             labelleft='off')
#         plt.tick_params(
#             axis='x',          # changes apply to the x-axis
#             which='both',      # both major and minor ticks are affected
#             bottom='off',      # ticks along the bottom edge are off
#             top='off',         # ticks along the top edge are off
#             labelbottom='off')
#         plt.xlim((ind[0],ind[-1]))
#         plt.savefig(figuredir + ranked_file.split('/')[-1] + 'test.png')

#         plt.close()

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
        # try:
        #     p = float(sum([x for x in simESsubset if x < actualES]))/float(sum([x for x in simESsubset if x > actualES]))
        # except ZeroDivisionError:
        #     print [x for x in simESsubset if x > actualES]
        #     p = 1.0
    else:
        simESsubset = [x for x in simES if x > 0]
        mu = np.mean(simESsubset)
        NES = actualES/mu
        # try:
        #     p = float(sum([x for x in simESsubset if x > actualES]))/float(sum([x for x in simESsubset if x < actualES]))
        # except ZeroDivisionError:
        #     print [x for x in simESsubset if x < actualES]
        #     p = 1.0

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
    # if p < FDRCUTOFF:
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


    # return [MOTIF_FILE,actualES,NES,p,simNES]
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
        # simNESlist = TFresults[i][4]
        PVAL = TFresults[i][3]


        # #Here we calculate a theoretical p-value for the actual value of this TF motif's NES against all other actual NESs (from the other TF motifs)
        # mu = np.mean(NESlist)
        # sigma = np.std(NESlist)
        # actualp = norm.cdf(NES,mu,sigma)
        # actualp = min(actualp,1-actualp)

        # #Here we calculate a theoretical p-value for the actual value of this TF motif's NES against the simulated NESs for this motif
        # simmu = np.mean(simNESlist)
        # simsigma = np.std(simNESlist)
        # simp = norm.cdf(NES,simmu,simsigma)
        # simp = min(simp,1-simp)

        # #The FDR is then the NES p-value against the distribution of simulated NESs divided by the NES p-value against the distribution of measured NESs
        # FDR1 = simp/actualp
        # TFresults[i].append(FDR1)


        #This is identical to the above calculation but uses EMPIRICAL p-value instead of theoretical
        # if NES < 0:
        #     q = ((sum([x for x in simNESlist if x<NES])/len(simNESlist))*sum([x for x in NESlist if x<0]))/sum([x for x in NESlist if x<NES])
        # else:
        #     q = ((sum([x for x in simNESlist if x>NES])/len(simNESlist))*sum([x for x in NESlist if x>0]))/sum([x for x in NESlist if x>NES])
        # TFresults[i].append(q)


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

def plot(TFresults):
    return

if __name__ == "__main__":
    # NESlist = np.random.randn(1000)
    # FDRlist = np.random.randn(1000)
    # sigx = NESlist[:10]
    # sigy = FDRlist[:10]
    # F = plt.figure()
    # plt.scatter(NESlist,FDRlist,color='black',edgecolor='')
    # plt.scatter(sigx,sigy,color='red',edgecolor='')
    # plt.title("TFEA Results Moustache Plot",fontsize=14)
    # plt.xlabel("Normalized Enrichment Score (NES)",fontsize=14)
    # plt.ylabel("False Discovery Rate (FDR)",fontsize=14)
    # plt.tick_params(
    #     axis='y',          # changes apply to the x-axis
    #     which='both',      # both major and minor ticks are affected
    #     left='off',        # ticks along the bottom edge are off
    #     right='off',       # ticks along the top edge are off
    #     labelleft='on')
    # plt.tick_params(
    #     axis='x',          # changes apply to the x-axis
    #     which='both',      # both major and minor ticks are affected
    #     bottom='off',      # ticks along the bottom edge are off
    #     top='off',         # ticks along the top edge are off
    #     labelbottom='on')
    # plt.savefig('../figures/TFEA_Results_Moustache_Plot.svg')
    # # plt.show()

    # ES = np.random.randn(1000)
    # F = plt.figure()
    # gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
    # ax0 = plt.subplot(gs[0])
    # ax0.plot(range(1,len(ES)+1),ES,color='green')
    # ax1 = plt.subplot(gs[1])
    # ax1.bar(range(1,len(ES)+1), [1 if x>2 else 0 for x in ES],color='black')
    # # plt.plot(range(1,len(ES)+1),ES)
    # ax0.tick_params(
    #     axis='y',          # changes apply to the x-axis
    #     which='both',      # both major and minor ticks are affected
    #     left='off',        # ticks along the bottom edge are off
    #     right='off',       # ticks along the top edge are off
    #     labelleft='on')
    # ax0.tick_params(
    #     axis='x',          # changes apply to the x-axis
    #     which='both',      # both major and minor ticks are affected
    #     bottom='off',      # ticks along the bottom edge are off
    #     top='off',         # ticks along the top edge are off
    #     labelbottom='off')
    # ax1.tick_params(
    #     axis='y',          # changes apply to the x-axis
    #     which='both',      # both major and minor ticks are affected
    #     left='off',        # ticks along the bottom edge are off
    #     right='off',       # ticks along the top edge are off
    #     labelleft='off')
    # ax1.tick_params(
    #     axis='x',          # changes apply to the x-axis
    #     which='both',      # both major and minor ticks are affected
    #     bottom='off',      # ticks along the bottom edge are off
    #     top='off',         # ticks along the top edge are off
    #     labelbottom='on')
    # ax0.set_title('Enrichment Plot: ',fontsize=14)
    # ax1.set_xlabel('Rank in Ordered Dataset', fontsize=14)
    # ax1.set_ylabel('Hits', fontsize=14)
    # ax0.set_ylabel('Enrichment Score (ES)', fontsize=14)
    # plt.savefig('../figures/HO_P53_HUMAN.H10MO.B.bed_enrichment_plot.svg')
    # # plt.show()

    MOTIF_FILE,ES,NES,PVAL,FDR = ['HO_P53_HUMAN.H10MO.B.bed',0.182143716966,6.22622338072,7.43595407471e-10,0.0]
    outfile = open('../figures/' + MOTIF_FILE + '.results.html','w')
    outfile.write("""<!DOCTYPE html>
    <html>
    <head>
    <title>"""+MOTIF_FILE+""" Results</title>
    <style>
    table {
        font-family: arial, sans-serif;
        border-collapse: collapse;
        width: 100%;
    }

    td, th {
        border: 1px solid #dddddd;
        text-align: left;
        padding: 8px;
    }

    tr:nth-child(even) {
        background-color: #dddddd;
    }
    </style>
    </head>
    <body>
    <h1>"""+MOTIF_FILE+""" Results</h1>
    <div>
        <div id="Positively Enriched" style="float: left;">
            <table> 
                <tr>
                    <th>TF Motif</th>
                    <th>ES</th> 
                    <th>NES</th>
                    <th>P-value</th>
                    <th>FDR</th>
                </tr>
                <tr>
                    <td>"""+MOTIF_FILE+"""</td>
                    <td>"""+str("%.3f" % ES)+"""</td>
                    <td>"""+str("%.3f" % NES)+"""</td>
                    <td>"""+str("%.4g" % PVAL)+"""</td>
                    <td>"""+str("%.4g" % FDR)+"""</td>
                </tr>
        </div>
    </div>
    <img src="../figures/""" + MOTIF_FILE + """_enrichment_plot.svg" alt="Enrichment Plot">
    </body>
    </html>""")
    outfile.close()


