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

if __name__ == "__main__":
    run('/scratch/Users/joru1876/TFEA_files/Allen2014/TFEA_output-0/temp_files/ranked_file.center.sorted.distance.bed')