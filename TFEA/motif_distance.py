__author__ = 'Jonathan Rubin'

import os
import math
from config import FILEDIR,LARGEWINDOW

#Takes a single ranked_file and returns bed file with the centers for each region and a bed file with the distance to the closest motif for each centered region
def run(ranked_center_file,MOTIF_PATH):
    #Get closest motif hit to center base of each region
    # os.system("bedtools closest -d -t first -a " + ranked_center_file.split('.bed')[0] + ".sorted.bed -b " + MOTIF_PATH + " > " + ranked_center_file.split('.bed')[0] + ".sorted.distance.bed")
    os.system("bedtools closest -D ref -t first -a " + ranked_center_file.split('.bed')[0] + ".sorted.bed -b " + MOTIF_PATH + " > " + '/' + '/'.join(ranked_center_file.split('/')[:-1]) + '/' + MOTIF_PATH.split('/')[-1] + ".sorted.distance.bed")

    return '/' + '/'.join(ranked_center_file.split('/')[:-1]) + '/' + MOTIF_PATH.split('/')[-1] + ".sorted.distance.bed"

def fimo_distance(fimo_file,MOTIF_FILE):
    d = dict()
    with open(fimo_file) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            pval,fc,rank = line[0].split(',')
            start,stop = line[1:2]
            distance = ((int(start)+int(stop))/2)-LARGEWINDOW
            score = line[5]
            if rank not in d:
                d[rank] = [rank,pval,fc,score,distance]
            else:
                prev_distance = math.fabs(d[rank][-1])
                if prev_distance > math.fabs(distance):
                    d[rank] = [rank,pval,fc,score,distance]

    outname = FILEDIR+MOTIF_FILE+'.sorted.distance.bed'
    outfile = open(outname,'w')
    for i in range(len(d)):
        rank = str(i)
        outfile.write('\t'.join(d[rank])+'\n')
    outfile.close()

    return outname


