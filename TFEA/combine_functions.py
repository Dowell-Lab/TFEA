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

    merge_bed_command = ("cat " + " ".join(beds) 
                        + " | bedtools sort -i stdin | bedtools merge -i stdin > " 
                        + combined_input_merged_bed)

    subprocess.call(merge_bed_command, shell=True)

    return combined_input_merged_bed

#==============================================================================
def intersect_all(beds=None, tempdir=None):
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
    combined_input_intersect_bed = os.path.join(tempdir, 
                                                "combined_input.intersect.bed")
    
    if len(beds) > 2:
        intersect_bed_command = ("bedtools intersect -a " + beds[0] 
                            + " -b " + beds[1])
        for bed in beds[2:]:
            intersect_bed_command = (intersect_bed_command 
                                    + " | bedtools intersect -a stdin -b "
                                    + bed)
    else:
        intersect_bed_command = ("bedtools intersect -a " + beds[0] 
                            + " -b " + beds[1])

    intersect_bed_command = (intersect_bed_command + " | bedtools sort -i stdin "
                            + " | bedtools merge -i stdin > " 
                            + combined_input_intersect_bed)

    # intersect_bed_command = ("bedtools intersect -a " + beds[0] 
    #                         + " -b " + " ".join(beds[1:]) 
    #                         + " | bedtools sort -i stdin | bedtools merge -i stdin > " 
    #                         + combined_input_intersect_bed)

    subprocess.call(intersect_bed_command, shell=True)

    return combined_input_intersect_bed

#==============================================================================
def tfit_remove_small(beds=None, tempdir=None, size_cut=200):
    '''
    '''
    combined_input_merged_bed = os.path.join(tempdir, 
                                                "combined_input.merge.bed")
    large_regions = os.path.join(tempdir, 'large_regions.bed')
    large_regions_file = open(large_regions,'w')

    for bedfile in beds:
        with open(bedfile) as F:
            for line in F:
                if '#' not in line:
                    chrom,start,stop = line.strip('\n').split('\t')[:3]
                    start = int(start)
                    stop = int(stop)
                    if stop-start >= size_cut:
                        large_regions_file.write('\t'.join([chrom,str(start),
                                                            str(stop)]) + '\n')

    large_regions_file.close()

    merge_bed_command = ("bedtools sort -i " + large_regions 
                        + " | bedtools merge -i stdin > " 
                        + combined_input_merged_bed)

    subprocess.call(merge_bed_command, shell=True)

    return combined_input_merged_bed

#==============================================================================
def tfit_clean(beds=None, tempdir=None, size_cut=500):
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
    command = ("cat " + large_regions + " " + small_regions_expanded 
                + " | bedtools sort -i stdin | bedtools merge -i stdin > "
                +  combined_input_merged_bed)
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

#==============================================================================
def combine_md(bed1=None, bed2=None, tempdir=None, method=None, size_cut=200, 
                largewindow=None, centerbed=False):
    md_bedfile1 = os.path.join(tempdir, 'md_bedfile1.bed')
    md_bedfile2 = os.path.join(tempdir, 'md_bedfile2.bed')
    if method == 'intersect/merge' or method == 'intersect all':
        if len(bed1) > 1:
            bed1_command = ("bedtools intersect -a "
                                + bed1[0] 
                                + " -b " + ' '.join(bed1[1:])
                                + " > " + md_bedfile1)
            bed2_command = ("bedtools intersect -a "
                                + bed2[0] 
                                + " -b " + ' '.join(bed2[1:])
                                + " > " + md_bedfile2)
        else:
            bed1_command = "cat " + bed1[0] + " > " + md_bedfile1
            bed2_command = "cat " + bed2[0] + " > " + md_bedfile2

    if method == 'merge all':
        bed1_command = ("cat " + ' '.join(bed1) 
            + " | bedtools sort -i stdin | bedtools merge -i stdin > " 
            + md_bedfile1)
        bed2_command = ("cat " + ' '.join(bed2) 
            + " | bedtools sort -i stdin | bedtools merge -i stdin > " 
            + md_bedfile2)

    if method == 'tfit remove small':
        for bedfile in bed1+bed2:
            outfile = os.path.join(tempdir, os.path.basename(bedfile))
            with open(outfile, 'w') as ofile:
                with open(bedfile) as F:
                    for line in F:
                        if '#' not in line:
                            chrom,start,stop = line.strip('\n').split('\t')[:3]
                            start = int(start)
                            stop = int(stop)
                            if stop-start >= size_cut:
                                ofile.write('\t'.join([chrom,str(start),
                                                        str(stop)]) + '\n')
        
        bed1_command = ("cat " + ' '.join([os.path.join(tempdir, os.path.basename(bedfile)) for bedfile in bed1])
                        + " | bedtools sort -i stdin | bedtools merge -i stdin > " 
                        + md_bedfile1)

        bed2_command = ("cat " + ' '.join([os.path.join(tempdir, os.path.basename(bedfile)) for bedfile in bed2])
                        + " | bedtools sort -i stdin | bedtools merge -i stdin > " 
                        + md_bedfile2)
    
    subprocess.call(bed1_command, shell=True)

    subprocess.call(bed2_command, shell=True)
    md_centerfile1 = os.path.join(tempdir, 'md_bedfile1.center.bed') 
    md_centerfile2 = os.path.join(tempdir, 'md_bedfile2.center.bed') 
    with open(md_bedfile1) as F:
        with open(md_centerfile1, 'w') as outfile:
            rank = 1
            for line in F:
                if '#' not in line[0]:
                    linelist = line.strip('\n').split('\t')
                    chrom, start, stop = linelist[:3]
                    if len(linelist) > 3:
                        rest = linelist[3:] + [str(rank)]
                    else:
                        rest = [str(rank)]
                    center = (int(start) + int(stop)) / 2
                    if centerbed:
                        outfile.write('\t'.join([chrom, str(int(center)), str(int(center+1))]) 
                                    + '\t' + ','.join(rest) + '\n')
                    else:
                        outfile.write('\t'.join([chrom, str(int(center-largewindow)), str(int(center+largewindow))]) 
                                    + '\t' + ','.join(rest) + '\n')
                    rank += 1
    md_sortedcenterfile1 = os.path.join(tempdir, 'md_bedfile1.center.sorted.bed')
    bedtools_sort = "bedtools sort -i " + md_centerfile1 + " > " + md_sortedcenterfile1
    subprocess.check_output(bedtools_sort, shell=True)
    
    with open(md_bedfile2) as F:
        with open(md_centerfile2, 'w') as outfile:
            rank = 1
            for line in F:
                if '#' not in line[0]:
                    linelist = line.strip('\n').split('\t')
                    chrom, start, stop = linelist[:3]
                    if len(linelist) > 3:
                        rest = linelist[3:] + [str(rank)]
                    else:
                        rest = [str(rank)]
                    center = (int(start) + int(stop)) / 2
                    if centerbed:
                        outfile.write('\t'.join([chrom, str(int(center)), str(int(center+1))]) 
                                    + '\t' + ','.join(rest) + '\n')
                    else:
                        outfile.write('\t'.join([chrom, str(int(center-largewindow)), str(int(center+largewindow))]) 
                                    + '\t' + ','.join(rest) + '\n')                    
                    rank += 1
    md_sortedcenterfile2 = os.path.join(tempdir, 'md_bedfile2.center.sorted.bed')
    bedtools_sort = "bedtools sort -i " + md_centerfile2 + " > " + md_sortedcenterfile2
    subprocess.check_output(bedtools_sort, shell=True)
    return md_sortedcenterfile1, md_sortedcenterfile2