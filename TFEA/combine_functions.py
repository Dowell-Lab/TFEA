#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''This file contains a list of functions associated with combining regions
    of interest (replicates/conditions) in a user-specified way
'''

#==============================================================================
__author__ = 'Jonathan D. Rubin and Rutendo F. Sigauke'
__credits__ = ['Jonathan D. Rubin', 'Rutendo F. Sigauke', 'Jacob T. Stanley',
                'Robin D. Dowell']
__maintainer__ = 'Jonathan D. Rubin'
__email__ = 'Jonathan.Rubin@colorado.edu'

#Imports
#==============================================================================
import os
import subprocess

#Functions
#==============================================================================
def merge_bed(beds=None, tempdir=None):
    '''Concatenates, sorts, and merges (bedtools) a list of bed files. Outputs 
        into the tempdir directory created by TFEA

    Parameters
    ----------
    beds : list or array
        full paths to bed files (strings)
        
    tempdir : string
        full path to tempdir directory in output directory (created by TFEA)

    Returns
    -------
    combined_input_merged_bed : string 
        full path to a bed file containing the merged regions inputted by the 
        user 
    '''
    combined_input_merged_bed = os.path.join(tempdir, 
                                                "combined_input.merge.bed")

    merge_bed_command = "cat " + " ".join(beds) 
                        + " | bedtools sort -i stdin | bedtools merge -i stdin > " 
                        + combined_input_merged_bed

    subprocess.call(merge_bed_command, shell=True)

    return combined_input_merged_bed

#==============================================================================
def tfit_clean_merge(beds=None, tempdir=None, size_cut=500):
    '''Takes in a list of bed files outputted by Tfit and for each region in
        each bed file splits the region based on its genomic size and the 
        size_cut variable. The 'large regions' get merged via bedtools, the 
        'small regions' get instersected via bedtools, then expanded to be
        size_cut in length, then merged with the 'large regions'

    Parameters
    ----------
    beds : list or array
        full paths to bed files (strings)
        
    tempdir : string
        full path to tempdir directory in output directory (created by TFEA)

    Returns
    -------
    combined_input_merged_bed : string 
        full path to a bed file containing the merged regions inputted by the 
        user 
    '''
    #First define a large_regions file and a small_regions list that will
    #store small regions files
    large_regions = os.path.join(tempdir, 'large_regions.bed')
    small_regions_list = list()
    large_regions_file = open(large_regions,'w')

    #Loop through the bed files and accordingly either append regions to 
    # the large regions file or to individual small regions files (append these
    # to the list)
    for bedfile in beds:
        small_regions = os.path.join(tempdir, 
            bedfile.split('/')[-1].split('.bed')[0] + '.small_regions.bed')
        small_regions_list.append(small_regions)
        small_regions_file = open(small_regions,'w')
        with open(bedfile) as F:
            for line in F:
                if '#' not in line:
                    chrom,start,stop = line.strip('\n').split('\t')[:3]
                    start = int(start)
                    stop = int(stop)
                    if stop-start >= size_cut:
                        large_regions_file.write('\t'.join([chrom,str(start),
                                                            str(stop)]) + '\n')
                    else:
                        small_regions_file.write('\t'.join([chrom,str(start),
                                                            str(stop)]) + '\n')
        small_regions_file.close()

    large_regions_file.close()

    #Perform bedtools intersect on the small regions
    command = ("bedtools intersect -a " + small_regions_list[0] + " -b " 
                + small_regions_list[1])
    for bedfile in small_regions_list[2:]:
        command = command + " | bedtools intersect -a stdin -b " + bedfile
    command = (command + " > " 
                + os.path.join(tempdir, "small_regions.intersect.bed"))
    # print command
    subprocess.call(command, shell=True)

    #Expand the intersected regions so they are size_cut in length
    small_regions_expanded = os.path.join(tempdir, 
                                        "small_regions.intersect.expanded.bed")
    small_regions_expanded_outfile = open(small_regions_expanded,'w')
    with open(os.path.join(tempdir, "small_regions.intersect.bed")) as F:
        for line in F:
            chrom,start,stop = line.strip('\n').split('\t')[:3]
            start = int(start)
            stop = int(stop)
            center = (start+stop)/2
            new_start = center - size_cut/2
            new_stop = center + size_cut/2
            small_regions_expanded_outfile.write('\t'.join([chrom, 
                                                            str(new_start), 
                                                            str(new_stop)])
                                                                + '\n')
    small_regions_expanded_outfile.close()

    #Define the output file
    combined_input_merged_bed = os.path.join(tempdir, 
                                                "combined_input.merge.bed")

    #Concatenate the large and small regions, sort, and merge
    command = "cat " + large_regions + " " + small_regions_expanded 
                + " | bedtools sort -i stdin | bedtools merge -i stdin > "
                +  combined_input_merged_bed
    subprocess.call(command, shell=True)

    return combined_input_merged_bed

#==============================================================================
def intersect_merge_bed(bed1=None, bed2=None, tempdir=None):
    '''Takes in two lists of bed files, each containing replicates for one 
        condition. Intersects replicates using bedtools then merges intersected
        regions.

    Parameters
    ----------
    bed1 : list or array
        full paths to bed files for condition 1 (strings)
    
    bed2 : list or array
        full paths to bed files for condition 2 (strings)
        
    tempdir : string
        full path to tempdir directory in output directory (created by TFEA)

    Returns
    -------
    combined_input_merged_bed : string 
        full path to a bed file containing the merged regions inputted by the 
        user 
    '''
    #Define the output file
    combined_input_merged_bed = os.path.join(tempdir, 
                                                "combined_input.merge.bed")

    if len(bed1) > 1:
        #Build command to perform bedtools intersect on condition1 beds
        intersect1 = ("bedtools intersect -a " + bed1[0] + " -b " + bed1[1])
        for bedfile in bed1[2:]:
            intersect1 = (intersect1 + " | bedtools intersect -a stdin -b " 
                        + bedfile)
    else:
        intersect1 = "cat " + bed1[0]

    if len(bed2) > 1:
        #Build command to perform bedtools intersect on condition2 beds
        intersect2 = ("bedtools intersect -a " + bed2[0] + " -b " + bed2[1])
        for bedfile in bed2[2:]:
            intersect2 = (intersect2 + " | bedtools intersect -a stdin -b " 
                        + bedfile)
    else:
        intersect2 = "cat " + bed2[0]

    #Build full command which pipes both intersect commands into cat, then 
    # sorts and merges this resulting bed file
    command = ("cat <(" + intersect1 + ") <(" + intersect2 
                + ") | bedtools sort -i stdin | bedtools merge -i stdin > " 
                + combined_input_merged_bed)
    
    #Need to use subprocess here because this command is bash not sh
    subprocess.call(['bash', '-c', command])

    return combined_input_merged_bed