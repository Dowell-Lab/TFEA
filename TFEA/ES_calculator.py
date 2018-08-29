__author__ = 'Jonathan Rubin'

global np
import matplotlib
matplotlib.use('Agg')
import sys
import math
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import gridspec
from scipy.stats import norm 
import time
from multiprocessing import Pool
import config
import motif_distance
import meta_eRNA
import create_html
import combine_bed
import meme
import main
import seaborn as sns


def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])
    
    return newdir

def permutations(distances, permutations=1000):
    '''generating permutations of the distances and calculating AUC.

    Parameters
    ----------
    distances : list or array
        normalized distances 
        
    permutations : int
        number of times to permute (default=1000)
        
    Returns
    -------
    es_permute : list 
        list of AUC calculated for permutations 
       
    '''

    es_permute = []
    triangle_area = 0.5*(len(distances))
    for i in range(permutations):
        random_distances = np.random.permutation(distances)
        cum_distances = np.cumsum(random_distances)
        es = np.trapz(cum_distances)
        auc = es - triangle_area
        es_permute.append(auc)
    return es_permute



def run(args):

    '''This function calculates the AUC for any TF based on the TF motif hits relative to the bidirectionals.
    The calculated AUC is used as a proxy for the enrichemnt of the TFs.

    Parameters                                                                                                                                                                           
    ----------                                                                                                                                                                           
    bedfile : tab delimeted                                                                                                                                                            
        tf bedfile with ranled distances                                                                                                                                                 
                                                                                                                                       
    Returns                                                                                                                                                                             
    -------                                                                                                                                                                             
    AUC : float                                                                                                                                                                   
        the AUC for a given tf

    p-value : float                                                                                                                                                                              imperical p-value calcutaed from 1000 simulations
                 
    Raises
    ------
    ValueError
        when tf does not have any hits

    '''

    MOTIF_FILE,millions_mapped = args
    if config.FIMO:
        ranked_fullregions_file = combine_bed.get_regions()
        ranked_fasta_file = combine_bed.getfasta(ranked_fullregions_file)
        background_file = combine_bed.get_bgfile(ranked_fasta_file)
        fimo_file = meme.fimo(background_file,MOTIF_FILE,ranked_fasta_file)
        meme.meme2images(MOTIF_FILE)
        ranked_center_distance_file = motif_distance.fimo_distance(fimo_file,MOTIF_FILE)
    else:
        ranked_center_distance_file = motif_distance.run(config.RANKED_CENTER_FILE,config.MOTIF_HITS+MOTIF_FILE)

    distances = []
    ranks = []
    pval = []
    pos = 0
    fc = []
    H = 1500.0
    h = 150.0

    with open(ranked_center_distance_file) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            distance = int(line[-1])
            rank = int(line[5])
            pval.append(float(line[3]))
            fc.append(float(line[4]))
            ranks.append(rank)
            distances.append(distance)
            if distance < H:
                pos+=1

    #sort distances based on the ranks from TF bed file
    #and calculate the absolute distance
    sorted_distances = [x for _,x in sorted(zip(ranks, distances))]
    distances_abs = [abs(x) for x in sorted_distances] 

    #filter any TFs/files without and hits
    if len(distances_abs) == 0:
        return "no hits"

    #Get -exp() of distance and get cumulative scores
    score = [math.exp(-x) for x in distances_abs] 
    total = float(sum(score))
    normalized_score = [x/total for x in score]
    cumscore = np.cumsum(normalized_score)

    #The AUC is the relative to the "random" line
    actualES = np.trapz(cumscore) - (0.5*len(cumscore))

    #Calculate random AUC
    simES = permutations(normalized_score)

    ##significance calculator                                                                                                                                                            
    mu = np.mean(simES)
    NES = actualES/abs(mu)
    sigma = np.std(simES)


    if actualES > 0:
        p = 1-norm.cdf(actualES,mu,sigma)
    else:
        p = norm.cdf(actualES,mu,sigma)   

    #Filter distances into quartiles for plotting purposes
    q1 = round(np.percentile(np.arange(1, len(distances_abs),1), 25))
    q3 = round(np.percentile(np.arange(1, len(distances_abs),1), 75))
    updistancehist = distances_abs[0:int(q1)]
    middledistancehist =  distances_abs[int(q1):int(q3)]
    downdistancehist = distances_abs[int(q3):len(distances_abs)]

    
    #Get log pval to plot for rank metric
    sorted_pval = [x for _,x in sorted(zip(ranks, pval))]
    sorted_fc = [x for _,x in sorted(zip(ranks, fc))]
    try:
        # logpval = [math.log(x,10) if y > 1 else -math.log(x,10) for x,y in zip(sorted_pval,sorted_fc)]
        logpval = [math.log(x,10) for x in sorted_pval]
    except ValueError:
        logpval = sorted_pval

    # #Plot results for significant hits while list of simulated ES scores is in memory                              
    if 'HO_' in MOTIF_FILE:
        os.system("scp " + config.LOGOS + MOTIF_FILE.split('.bed')[0].split('HO_')[1] + "_direct.png " + config.FIGUREDIR)
        os.system("scp " + config.LOGOS + MOTIF_FILE.split('.bed')[0].split('HO_')[1] + "_revcomp.png " + config.FIGUREDIR)
    else:
        os.system("scp " + config.LOGOS + MOTIF_FILE.split('.bed')[0] + "_direct.png " + config.FIGUREDIR)
        os.system("scp " + config.LOGOS + MOTIF_FILE.split('.bed')[0] + "_revcomp.png " + config.FIGUREDIR)


    #Begin plotting section
    F = plt.figure(figsize=(15.5,8))
    xvals = range(1,len(cumscore)+1)
    limits = [1,len(cumscore)]
    gs = gridspec.GridSpec(4, 1, height_ratios=[2, 2, 1, 1])

    #This is the enrichment score plot (i.e. line plot)
    ax0 = plt.subplot(gs[0])
    ax0.plot(xvals,cumscore,color='green')
    ax0.plot([0, len(cumscore)],[0, 1], '--', alpha=0.75)
    ax0.set_title('Enrichment Plot: ',fontsize=14)
    ax0.set_ylabel('Enrichment Score (ES)', fontsize=10)
    ax0.tick_params(axis='y', which='both', left='on', right='off', labelleft='on')
    ax0.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
    ylims = ax0.get_ylim()
    ymax = math.fabs(max(ylims,key=abs))
    ax0.set_ylim([0,ymax])
    ax0.set_xlim(limits)

    #This is the distance scatter plot right below the enrichment score plot
    ax1 = plt.subplot(gs[1])
    ax1.scatter(xvals,distances,edgecolor="",color="black",s=10,alpha=0.25)
    # ax1.axhline(config.SMALLWINDOW, color='red',alpha=0.25)
    # ax1.axhline(-config.SMALLWINDOW, color='red',alpha=0.25)
    ax1.tick_params(axis='y', which='both', left='off', right='off', labelleft='on')
    ax1.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
    ax1.set_xlim(limits)
    ax1.set_ylim([-int(config.LARGEWINDOW),int(config.LARGEWINDOW)])
    plt.yticks([-int(config.LARGEWINDOW),0,int(config.LARGEWINDOW)],[str(-int(config.LARGEWINDOW)/1000.0),'0',str(int(config.LARGEWINDOW)/1000.0)])
    ax1.set_ylabel('Distance (kb)', fontsize=10)

    #This is the rank metric plot
    ax2 = plt.subplot(gs[3])
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
        ax2.axvline(len(updistancehist)+1,color='green',alpha=0.25)
    except ValueError:
        pass
    try:
        ax2.axvline(len(xvals) - len(downdistancehist), color='purple', alpha=0.25)
    except ValueError:
        pass

    #This is the GC content plot
    ax3 = plt.subplot(gs[2])
    # ax3.set_ylim([-config.SMALLWINDOW,config.SMALLWINDOW])
    ax3.set_xlim(limits)
    GC_ARRAY = np.array(config.GC_ARRAY).transpose()
    sns.heatmap(GC_ARRAY, cbar=False, xticklabels='auto',yticklabels='auto') #, cbar_ax=F.add_axes([1, 1, .03, .4]))
    ax3.set_ylim([-int(config.LARGEWINDOW),int(config.LARGEWINDOW)])
    plt.yticks([-int(config.LARGEWINDOW),0,int(config.LARGEWINDOW)],[str(-int(config.LARGEWINDOW)/1000.0),'0',str(int(config.LARGEWINDOW)/1000.0)])
    # plt.imshow(GC_ARRAY, cmap='hot', interpolation='nearest')
    ax3.tick_params(axis='y', which='both', left='on', right='off', labelleft='on')
    ax3.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
    # ax3.yaxis.set_ticks([int(-config.SMALLWINDOW),0,int(config.SMALLWINDOW)])
    # ax3.set_xlabel('Rank in Ordered Dataset', fontsize=14)
    ax3.set_ylabel('Position (kb)',fontsize=10)


    plt.savefig(config.FIGUREDIR + MOTIF_FILE.split('.bed')[0] + '_enrichment_plot.png',bbox_inches='tight')
    plt.cla()

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


    #Plots the distribution of motif distances with a red line at h                                                                                                                                                                   
    F = plt.figure(figsize=(6.5,6))
    gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 1])
    ax0 = plt.subplot(gs[0])
    binwidth = H/100.0
    ax0.hist(updistancehist,bins=np.arange(0,int(H)+binwidth,binwidth),color='green')
    ax0.set_title('Distribution of Motif Distance for: fc > 1',fontsize=14)
    ax0.axvline(h,color='red',alpha=0.5)
    ax0.tick_params(axis='y', which='both', left='off', right='off', labelleft='on')
    ax0.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
    ax0.set_xlim([0,H])
    ##ax0.set_xlabel('Distance (bp)',fontsize=14)
    ax0.set_ylabel('Hits',fontsize=14)
    ax1 = plt.subplot(gs[2])
    ax1.hist(downdistancehist,bins=np.arange(0,int(H)+binwidth,binwidth),color='purple')
    ax1.axvline(h,color='red',alpha=0.5)
    ax1.set_title('Distribution of Motif Distance for: fc < 1',fontsize=14)
    ax1.tick_params(axis='y', which='both', left='off', right='off', labelleft='on')
    ax1.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
    ax1.set_xlim([0,H])
    ax1.set_ylabel('Hits',fontsize=14)
    ax1.set_xlabel('Distance (bp)',fontsize=14)
    ax2 = plt.subplot(gs[1])
    ax2.hist(middledistancehist,bins=np.arange(0,int(H)+binwidth,binwidth),color='blue')
    ax2.set_title('Distribution of Motif Distance for: middle',fontsize=14)
    ax2.axvline(h,color='red',alpha=0.5)
    ax2.tick_params(axis='y', which='both', left='off', right='off', labelleft='on')
    ax2.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
    ax2.set_xlim([0,H])
    ##ax2.set_xlabel('Distance (bp)',fontsize=14)
    ax2.set_ylabel('Hits',fontsize=14)
    plt.savefig(config.FIGUREDIR + MOTIF_FILE.split('.bed')[0] + '_distance_distribution.png',bbox_inches='tight')                                                                                                             
    plt.cla()

    return [MOTIF_FILE.split('.bed')[0],actualES,NES,p,pos]

def PADJ(TFresults):
    #This function iterates through the results and calculates a p-adj for each TF motif. Also creates a moustache plot.                                                                                                                         
    PADJlist = list()
    NESlist = list()
    ESlist = list()
    sigx = list()
    sigy = list()
    pvals = list()
    MAy = list()
    MAx = list()
    ##pvalsNA = [1 if str(x) == 'nan' else x for x in pvals]                                                                                                                                                                                    
    totals = list()
    sigtotals = list()
    positives = list()
    h = 150.0
    H = 1500.0

    TFresults = [x for x in TFresults if x != "no hits"]
    TFresults = sorted(TFresults, key=lambda x: x[3])
    for i in range(len(TFresults)):
        ES = TFresults[i][1]
        NES = TFresults[i][2]
        PVAL = TFresults[i][3]
        try:
            POS = math.log(TFresults[i][4],10)
        except:
            POS = 0.0
        # NEG = TFresults[i][5]
        pvals.append(PVAL)
        # total = POS+NEG
        #Using classical FDR calculation ((pvalue*(# hypotheses tested))/rank of p-value)                                                                                                                                                       
        # FDR = ((float(i)+1.0)/float(len(TFresults)))*config.FDRCUTOFF
        # FDR = (PVAL*len(TFresults))/float(i+1.0)                   
        #Using Bonferroni Correction
        PADJ = 1 if PVAL*len(TFresults) > 1 else PVAL*len(TFresults)
        TFresults[i].append(PADJ)
        PADJlist.append(PADJ)
        NESlist.append(ES)
        ESlist.append(ES)
        # totals.append(math.log(float(total)))
        positives.append(POS)

        if PADJ < config.PADJCUTOFF:
            sigx.append(ES)
            sigy.append(PADJ)
            MAy.append(ES)
            MAx.append(POS)
            #sigtotals.append(math.log(float(total)))

    #Creates a moustache plot of the global PADJs vs. ESs                                                                                                                                                                                       
    F = plt.figure(figsize=(7,6))
    ax = plt.subplot(111)
    ax.scatter(ESlist,PADJlist,color='black',edgecolor='')
    ax.scatter(sigx,sigy,color='red',edgecolor='')
    ax.set_title("TFEA Moustache Plot",fontsize=14)
    ax.set_xlabel("Enrichment Score (ES)",fontsize=14)
    ax.set_ylabel("Adjusted p-value (PADJ)",fontsize=14)
    xlimit = math.fabs(max(ESlist,key=abs))
    ylimit = math.fabs(max(PADJlist,key=abs))
    ax.set_xlim([-xlimit,xlimit])
    ax.set_ylim([0,ylimit])
    ax.tick_params(axis='y', which='both', left='off', right='off', labelleft='on')
    ax.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
    plt.savefig(config.FIGUREDIR + 'TFEA_Results_Moustache_Plot.png',bbox_inches='tight')
    plt.cla()

    #Creates a histogram of p-values                                                                                                                                                                                                            
    F = plt.figure(figsize=(7,6))
    ax = plt.subplot(111)
    binwidth = 1.0/100.0
    ax.hist(pvals,bins=np.arange(0,0.5+binwidth,binwidth),color='green')
    ax.set_title("TFEA P-value Histogram",fontsize=14)
    ax.set_xlabel("P-value",fontsize=14)
    ax.set_ylabel("Count",fontsize=14)
    ax.tick_params(axis='y', which='both', left='off', right='off', labelleft='on')
    ax.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
    plt.savefig(config.FIGUREDIR + 'TFEA_Pval_Histogram.png',bbox_inches='tight')
    plt.cla()

    #Creates an MA-plot with NES on Y-axis and positive hits on X-axis                                                                                                                                                                                                      
    F = plt.figure(figsize=(7,6))
    ax = plt.subplot(111)
    ax.scatter(positives,ESlist,color='black',edgecolor='')
    ax.scatter(MAx,MAy,color='red',edgecolor='')
    ax.set_title("TFEA MA-Plot",fontsize=14)
    ax.set_ylabel("Normalized Enrichment Score (NES) " + config.LABEL2 + "/" + config.LABEL1,fontsize=14)
    ax.set_xlabel("Hits Log10",fontsize=14)
    ax.tick_params(axis='y', which='both', left='off', right='off', labelleft='on')
    ax.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
    plt.savefig(config.FIGUREDIR + 'TFEA_NES_MA_Plot.png',bbox_inches='tight')
    plt.cla()


    return TFresults





