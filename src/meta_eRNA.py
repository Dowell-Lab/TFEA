__author__ = 'Jonathan Rubin'

import HTSeq as hts
from config import *

def run(ranked_center_distance_file):
    H=1500.0
    #This section of the code is new (2/5/18) and will draw a metagene plot with a motif heatmap underneath
    #As of 2/6/18 this section has been moved to a new file and is incomplete
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
    profile1 = np.zeros(2*int(H), dtype='i')
    profile2 = np.zeros(2*int(H), dtype='i')
    for region in regions:
        window = hts.GenomicInterval(region.chrom, region.start, region.end, ".")
        for bam1 in hts_bam1:
            coverage1 = hts.GenomicArray( "auto", stranded=False, typecode="i" )
            print "doing coverage..."
            for almnt in bam1:
                if almnt.aligned:
                    coverage1[almnt.iv] += 1
            print "done\nadding cov to window..."
            wincvg1 = np.fromiter(coverage1[window], dtype='i', count=2*H)
            print "done\nadding window to profile..."
            profile1 += [float(x)/float(sum(wincvg1)) for x in wincvg1]
            print profile1
        for coverage2 in hts_bam2:
            wincvg2 = np.fromiter(coverage2[window], dtype='i', count=2*H)
            profile2 += [x/sum(wincvg2) for x in wincvg2]

    print profile1,profile2
if __name__ == "__main__":