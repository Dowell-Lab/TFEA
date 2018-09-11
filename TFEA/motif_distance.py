__author__ = 'Jonathan Rubin'

import os
import math
import numpy as np
import combine_bed, config

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
            distance = ((int(start)+int(stop))/2)-config.LARGEWINDOW
            score = line[5]
            if rank not in d:
                d[rank] = [rank,pval,fc,score,distance]
            else:
                prev_distance = math.fabs(d[rank][-1])
                if prev_distance > math.fabs(distance):
                    d[rank] = [rank,pval,fc,score,distance]

    outname = config.FILEDIR+MOTIF_FILE+'.sorted.distance.bed'
    outfile = open(outname,'w')
    for i in range(len(d)):
        rank = str(i)
        outfile.write('\t'.join(d[rank])+'\n')
    outfile.close()

    return outname

#Function to calculate GC content for all regions - creates panel underneath enrichment score figure in html output
def get_gc(ranked_file,window=int(config.LARGEWINDOW),bins=1000):
    '''This function calculates gc content over all eRNAs. It uses the LARGEWINDOW variable within
    the config file instead of the whole region. It performs a running average of window size = total_regions/1000

    Parameters                                                                                                                                                                           
    ----------                                                                                                                                                                           
    ranked_bed_file : tab delimited bedfile with regions ranked by some metric of 
        differential transcription
    window : an int. The size of the window (in bp) that you wish to calculate gc content for
    bins : an int. The number of equal sized bins to compute for the gc array                                                                                                                                           
                                                                                                                                       
    Returns                                                                                                                                                                             
    -------                                                                                                                                                                             
    AUC : float                                                                                                                                                                   
        the AUC for a given tf

    p-value : float                                                                                                                                         
                 
    Raises
    ------
    None

    '''

    #First, create a bed file with the correct coordinates centered on the given regions with the specified window size on either side
    outfile = open(config.FILEDIR+"ranked_file.windowed.bed",'w')
    with open(config.FILEDIR + "ranked_file.bed") as F:
        for line in F:
            line = line.strip('\n').split('\t')
            chrom,start,stop = line[:3]
            center = (int(start)+int(stop))/2
            newstart = center - window
            newstop = center + window
            outfile.write(chrom + '\t' + str(newstart) + '\t' + str(newstop) + '\t' + '\t'.join(line[3:]) + '\n')
    outfile.close()

    #Convert the bed file created above into a fasta file using meme (contained within the combine_bed file for some reason)
    ranked_file_windowed_fasta = combine_bed.getfasta(config.FILEDIR+"ranked_file.windowed.bed")

    #Create a gc_array which simply contains all sequences in fasta file collapsed into 1.0 for G/C and 0.0 for A/T
    gc_array = []
    with open(ranked_file_windowed_fasta) as F:
        for line in F:
            if '>' not in line:
                line = line.strip('\n')
                gc_array.append(convert_sequence_to_array(line))

    #The length of each bin is equal to the total positions over the number of desired bins
    binwidth = len(gc_array)/bins

    #Collapse the gc_array into an array containing the correct number of bins (specified by user)
    ##First, step through the gc_array with binwidth step size (i)
    for i in range(0,len(gc_array),binwidth):
        ##Initialize a position_average list that will store the mean value for each position (along the window)
        position_average = []
        ##Now we step through the total window size (window*2) position by position
        for k in range(window*2):
            ##new_array stores for each position in the window, a binwidth amount of data points to be averaged
            new_array = []
            ##Now, if we are not at the end of the gc_array, we will step through for each position, a binwidth amount of values and store them in new_array
            if i+binwidth < len(gc_array):
                for j in range(i,i+binwidth):
                    new_array.append(gc_array[j][k])
            else:
                for j in range(i,len(gc_array)):
                    new_array.append(gc_array[j][k])
            ##Finally, we simply append the average of the new_array into the position_average list
            position_average.append(np.mean(new_array))
        ##And this poition_average list is appended to the actual GC_ARRAY within the config file for later use
        config.GC_ARRAY.append(position_average)





def convert_sequence_to_array(sequence):

    '''This function converts a DNA sequence (ACGT alphabet) to an array, collapsing GCs
    to 1's and ATs to 0's

    Parameters                                                                                                                                                                           
    ----------                                                                                                                                                                           
    sequence : a string containing ACGT character                                                                                                                                               
                                                                                                                                       
    Returns                                                                                                                                                                             
    -------                                                                                                                                                                             
    array : list
        of floats corresponding to 1 or 0 depending on whether the 
        input sequence was GC or AT at a specific site       

    Raises
    ------
    ValueError
        when a character is not ACTG                                                                                                                                                  

    '''
    array = []
    for character in sequence:
        character = character.upper()
        if character == 'G' or character == 'C':
            array.append(1.0)
        elif character == 'A' or character == 'T':
            array.append(0.0)
        else:
            print "Warning: Character not recognized"
            array.append(0.0)

    return array 


if __name__ == "__main__":
    #=====================TEST=============================
    # sequence = 'actgACTG'
    # array = convert_sequence_to_array(sequence)
    # print "Sequence to Array passed test: ",array == [0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0]


    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np


    gc_array = [[0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0] for i in range(1009)]
    gc_array.append([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    bins = 100
    binwidth = len(gc_array)/bins
    window = len(gc_array[0])/2
    new_GC_array = []

    #Collapse the gc_array into an array containing the correct number of bins (specified by user)
    ##First, step through the gc_array with binwidth step size (i)
    for i in range(0,len(gc_array),binwidth):
        ##Initialize a position_average list that will store the mean value for each position (along the window)
        position_average = []
        ##Now we step through the total window size (window*2) position by position
        for k in range(window*2):
            ##new_array stores for each position in the window, a binwidth amount of data points to be averaged
            new_array = []
            ##Now, if we are not at the end of the gc_array, we will step through for each position, a binwidth amount of values and store them in new_array
            if i+binwidth < len(gc_array):
                for j in range(i,i+binwidth):
                    new_array.append(gc_array[j][k])
            else:
                for j in range(i,len(gc_array)):
                    new_array.append(gc_array[j][k])
            ##Finally, we simply append the average of the new_array into the position_average list
            position_average.append(np.mean(new_array))
        ##And this poition_average list is appended to the actual GC_ARRAY within the config file for later use
        new_GC_array.append(position_average)

    new_GC_array = np.array(new_GC_array)
    sns.heatmap(new_GC_array.transpose())
    # plt.imshow(new_GC_array.transpose(), cmap='hot', interpolation='nearest')
    plt.show()


