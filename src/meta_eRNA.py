__author__ = 'Jonathan Rubin'

import matplotlib
matplotlib.use('Agg')
import sys
import math
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import time
import HTSeq as hts
import numpy as np
from config import *

def run(ranked_center_distance_file):
    H=1500
    #This section of the code is new (2/5/18) and will draw a metagene plot with a motif heatmap underneath
    #As of 2/6/18 this section has been moved to a new file and is incomplete
    regions=list()
    with open(ranked_center_distance_file) as F:
            for line in F:
                line = line.strip('\n').split('\t')
                chrom,start,stop = line[:3]
                distance = float(line[-1])
                if 0 < distance < H:
                    regions.append(hts.GenomicInterval(chrom,int(start)-int(H),int(start)+int(H),'.'))
    #This section populates HTSeq variables that will read bed and bam files
    hts_bam1 = list()
    hts_bam2 = list()
    for bam in BAM1:
        hts_bam1.append(hts.BAM_Reader(bam))
    for bam in BAM2:
        hts_bam2.append(hts.BAM_Reader(bam))

    #This section creates two profile variables (one for each condition tested) that stores normalized coverage data
    profile1 = np.zeros(2*int(H), dtype='float64')
    profile2 = np.zeros(2*int(H), dtype='float64')
    for region in regions:
        window = hts.GenomicInterval(region.chrom, region.start, region.end, ".")
        for bam1 in hts_bam1:
            coverage1 = hts.GenomicArray( "auto", stranded=False, typecode="d" )
            for almnt in bam1:
                if almnt.aligned:
                    coverage1[almnt.iv] += 1
            wincvg1 = np.fromiter(coverage1[window], dtype='float64', count=2*H)
            profile1 += [x/sum(wincvg1) for x in wincvg1]
        for bam2 in hts_bam2:
            coverage2 = hts.GenomicArray( "auto", stranded=False, typecode="d" )
            for almnt in bam2:
                if almnt.aligned:
                    coverage2[almnt.iv] += 1
            wincvg2 = np.fromiter(coverage2[window], dtype='float64', count=2*H)
            profile2 += [x/sum(wincvg2) for x in wincvg2]

    F = plt.figure(figsize=(16.5,6))
    xvals = range(-H,H)
    gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1])
    ax0 = plt.subplot(gs[0])
    ax0.plot(xvals,profile1,color='green')
    ax1 = plt.subplot(gs[1])
    ax1.plot(xvals,profile2,color='red')
    plt.savefig('/scratch/Users/joru1876/TFEA_files/Allen2014/TFEA_output-0/meta_eRNA_test.png')

def run2(ranked_center_distance_file):

    regions=list()
    with open(ranked_center_distance_file) as F:
            for line in F:
                line = line.strip('\n').split('\t')
                chrom,start,stop = line[:3]
                distance = float(line[-1])
                if 0 <= distance <= H:
                    regions.append(hts.GenomicInterval(chrom,int(start)-int(H),int(start)+int(H),'.'))


    hts_bam1 = list()
    hts_bam2 = list()
    for bam in BAM1:
        hts_bam1.append(hts.BAM_Reader(bam))
    for bam in BAM2:
        hts_bam2.append(hts.BAM_Reader(bam))

    posprofile1 = np.zeros(2*int(H))   
    negprofile1 = np.zeros(2*int(H))
    posprofile2 = np.zeros(2*int(H))   
    negprofile2 = np.zeros(2*int(H))
    for window in regions:
        avgposprofile1 = np.zeros(2*int(H))
        avgnegprofile1 = np.zeros(2*int(H))
        for sortedbamfile in hts_bam1:
            tempposprofile = np.zeros(2*int(H))
            tempnegprofile = np.zeros(2*int(H))
            for almnt in sortedbamfile[ window ]:
                if almnt.iv.strand == '+':
                    start_in_window = almnt.iv.start - window.start
                    end_in_window   = almnt.iv.end - window.end + 2*int(H)
                    start_in_window = max( start_in_window, 0 )
                    end_in_window = min( end_in_window, 2*int(H) )
                    tempposprofile[ start_in_window : end_in_window ] += 1.0
                if almnt.iv.strand == '-':
                    start_in_window = almnt.iv.start - window.start
                    end_in_window   = almnt.iv.end - window.end + 2*int(H)
                    start_in_window = max( start_in_window, 0 )
                    end_in_window = min( end_in_window, 2*int(H) )
                    tempnegprofile[ start_in_window : end_in_window ] += -1.0
            pos_sum = np.sum(tempposprofile)
            neg_sum = np.sum(tempnegprofile)
            if pos_sum != 0:
                tempposprofile = [x/pos_sum for x in tempposprofile]
            if neg_sum != 0:
                tempnegprofile = [-(x/neg_sum) for x in tempnegprofile]
            avgposprofile1 = [x+y for x,y in zip(avgposprofile1,tempposprofile)]
            avgnegprofile1 = [x+y for x,y in zip(avgnegprofile1, tempnegprofile)]
        repnumber = len(hts_bam1)
        avgposprofile1 = [x/repnumber for x in avgposprofile1]
        avgnegprofile1 = [x/repnumber for x in avgnegprofile1]
        posprofile1 = [x+y for x,y in zip(posprofile1,avgposprofile1)]
        negprofile1 = [x+y for x,y in zip(negprofile1, avgnegprofile1)]

        avgposprofile2 = np.zeros(2*int(H))
        avgnegprofile2 = np.zeros(2*int(H))
        for sortedbamfile in hts_bam2:
            tempposprofile = np.zeros(2*int(H))
            tempnegprofile = np.zeros(2*int(H))
            for almnt in sortedbamfile[ window ]:
                if almnt.iv.strand == '+':
                    start_in_window = almnt.iv.start - window.start
                    end_in_window   = almnt.iv.end - window.end + 2*int(H)
                    start_in_window = max( start_in_window, 0 )
                    end_in_window = min( end_in_window, 2*int(H) )
                    tempposprofile[ start_in_window : end_in_window ] += 1.0
                if almnt.iv.strand == '-':
                    start_in_window = almnt.iv.start - window.start
                    end_in_window   = almnt.iv.end - window.end + 2*int(H)
                    start_in_window = max( start_in_window, 0 )
                    end_in_window = min( end_in_window, 2*int(H) )
                    tempnegprofile[ start_in_window : end_in_window ] += -1.0
            pos_sum = np.sum(tempposprofile)
            neg_sum = np.sum(tempnegprofile)
            if pos_sum != 0:
                tempposprofile = [x/pos_sum for x in tempposprofile]
            if neg_sum != 0:
                tempnegprofile = [-(x/neg_sum) for x in tempnegprofile]
            avgposprofile2 = [x+y for x,y in zip(avgposprofile2,tempposprofile)]
            avgnegprofile2 = [x+y for x,y in zip(avgnegprofile2, tempnegprofile)]
        repnumber = len(hts_bam2)
        avgposprofile2 = [x/repnumber for x in avgposprofile2]
        avgnegprofile2 = [x/repnumber for x in avgnegprofile2]
        posprofile2 = [x+y for x,y in zip(posprofile2,avgposprofile2)]
        negprofile2 = [x+y for x,y in zip(negprofile2, avgnegprofile2)]

    return posprofile1, negprofile1, posprofile2, negprofile2


if __name__ == "__main__":
    posprofile1, negprofile1, posprofile2, negprofile2 = run2('/scratch/Users/joru1876/TFEA_files/IRIS/TFEA_output-0/temp_files/ranked_file.center.sorted.distance.bed')
    print sum(posprofile1), sum(negprofile1), sum(posprofile2), sum(negprofile2)

    F = plt.figure(figsize=(15.5,6))
    ax0 = plt.subplot(111)
    xvals = range(-int(H),int(H))
    line1, = ax0.plot(xvals,posprofile1,color='blue',label=LABEL1)
    ax0.plot(xvals,negprofile1,color='blue')
    line2, = ax0.plot(xvals,posprofile2,color='red',label=LABEL2)
    ax0.plot(xvals,negprofile2,color='red')
    ax0.legend(loc=1)
    ax0.set_title('Meta Plot of eRNA Signal',fontsize=14)
    ax0.tick_params(axis='y', which='both', left='off', right='off', labelleft='on')
    ax0.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
    ax0.set_ylabel('Normalized Read Coverage',fontsize=14)
    ax0.set_xlabel('Distance to eRNA Origin (bp)')
    plt.savefig('/scratch/Users/joru1876/TFEA_files/IRIS/TFEA_output-0/meta_eRNA.png',bbox_inches='tight')
    plt.cla()