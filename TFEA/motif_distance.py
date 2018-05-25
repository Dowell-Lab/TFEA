__author__ = 'Jonathan Rubin'

import os

#Takes a single ranked_file and returns bed file with the centers for each region and a bed file with the distance to the closest motif for each centered region
def run(ranked_center_file,MOTIF_PATH):
    #Get closest motif hit to center base of each region
    # os.system("bedtools closest -d -t first -a " + ranked_center_file.split('.bed')[0] + ".sorted.bed -b " + MOTIF_PATH + " > " + ranked_center_file.split('.bed')[0] + ".sorted.distance.bed")
    os.system("bedtools closest -D ref -t first -a " + ranked_center_file.split('.bed')[0] + ".sorted.bed -b " + MOTIF_PATH + " > " + '/' + '/'.join(ranked_center_file.split('/')[:-1]) + '/' + MOTIF_PATH.split('/')[-1] + ".sorted.distance.bed")

    return '/' + '/'.join(ranked_center_file.split('/')[:-1]) + '/' + MOTIF_PATH.split('/')[-1] + ".sorted.distance.bed"