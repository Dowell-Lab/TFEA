#!/usr/bin/env python3
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
import sys
import time
import datetime
import subprocess
from pathlib import Path

import pybedtools
from pybedtools import BedTool, featurefuncs

from TFEA import config
from TFEA import multiprocess
from TFEA import exceptions

#Main Script
#==============================================================================
def main(use_config=True, bed1=None, bed2=None, method=None, tempdir=None, 
        md=None, largewindow=None, scanner=None, debug=False, label1=None, 
        label2=None, jobid=None):
    '''This is the main script of the combine function that is called within
        TFEA. Default arguments are assigned to variables within config.vars.

    Parameters
    ----------
    use_config : boolean
        Whether to use a config module to assign variables.
    bed1 : list
        A list of strings specifying full paths to bed files corresponding to
        a single condition (replicates)
    bed2 : list
        A list of strings specifying full paths to bed files corresponding to
        a single condition (replicates)
    method : str
        Method for combining input bed files into a single bed file
    tempdir : str
        Full path to a directory where files will be saved
    md : boolean
        Whether md-score bed files are generated
    largewindow : int
        Half-length of window size to use when generating md-score related
        bed files
    scanner : str
        Scanner method to use in SCANNER module. Only needed if md also
        specified. If equal to 'genome hits', md bed files generated will be 
        only contain one base and be centered at the middle of the region

    Returns
    -------
    None - Assigns varaibles within config if use_config set to True

    Raises
    ------
    FileEmptyError
        If any resulting file is empty
    '''
    start_time = time.time()
    if use_config:
        bed1 = config.vars['BED1']
        bed2 = config.vars['BED2']
        method = config.vars['COMBINE']
        tempdir = config.vars['TEMPDIR']
        md = config.vars['MD']
        md_bedfile1 = config.vars['MD_BEDFILE1']
        md_bedfile2 = config.vars['MD_BEDFILE2']
        largewindow = config.vars['LARGEWINDOW']
        scanner = config.vars['SCANNER']
        label1 = config.vars['LABEL1']
        label2 = config.vars['LABEL2']
        debug = config.vars['DEBUG']
        jobid = config.vars['JOBID']


    print("Combining Regions...", end=' ', flush=True, file=sys.stderr)

    if md_bedfile1 and md_bedfile2:
        centered_md_bedfile1 = tempdir / 'md_bedfile1.centered.bed'
        centered_md_bedfile2 = tempdir / 'md_bedfile2.centered.bed'
        md = md and (not md_bedfile1 or not md_bedfile2) #Boolean to determine whether to generate MD bed files
        md_pybedtool1 = BedTool(str(md_bedfile1))
        md_pybedtool1.each(center_feature).each(extend_feature, size=largewindow).remove_invalid().saveas(centered_md_bedfile1)
        md_pybedtool2 = BedTool(str(md_bedfile2))
        md_pybedtool2.each(center_feature).each(extend_feature, size=largewindow).remove_invalid().saveas(centered_md_bedfile2)
        if use_config:
            config.vars['MD_BEDFILE1'] = centered_md_bedfile1 
            config.vars['MD_BEDFILE2'] = centered_md_bedfile2
        

    #Use MuMerge to merge bed files
    if method == 'mumerge':
        mumerge_input = tempdir / 'mumerge_input.txt'
        combined_file = tempdir / 'combined_file.mumerge'
        #Write MuMerge input file
        # with open(mumerge_input, 'w') as F:
        #     F.write("#file\tsampid\tgroup\n")
        #     for i,bedpath in enumerate(bed1, 1):
        #         F.write(f'{bedpath}\t{label1}{i}\t{label1}\n')
        #     for i,bedpath in enumerate(bed2, 1):
        #         F.write(f'{bedpath}\t{label2}{i}\t{label2}\n')
        
        #MuMerge Command - output to combined_file.mumerge.bed
        combined_file = mumerge(mumerge_input, combined_file, bed1=bed1, 
                                bed2=bed2, label1=label1, label2=label2)
        clean_combined_file = tempdir / 'combined_file.mumerge.clean.bed'
        combined_pybedtool = BedTool(str(combined_file))
        combined_pybedtool.remove_invalid().saveas(clean_combined_file)
        combined_file = clean_combined_file
        # combined_file = Path(str(combined_file) + '_MUMERGE.bed')

        #Perform simple merge same as merge all for md bed files
        if md:
            md_bedfile1 = tempdir / "md_bedfile1.mumerge"
            md_mumerge_input1 = tempdir / "md_mumerge_input1.txt"
            md_bedfile1 = mumerge(md_mumerge_input1, md_bedfile1, bed1=bed1, 
                                    label1=label1, label2=label2)
            md_pybedtool1 = BedTool(str(md_bedfile1))
            md_bedfile1 = tempdir / "md_bedfile1.mumerge.final.bed"
            md_pybedtool1.each(center_feature).each(extend_feature, size=largewindow).remove_invalid().saveas(md_bedfile1)
            md_bedfile2 = tempdir / "md_bedfile2.mumerge"
            md_mumerge_input2 = tempdir / "md_mumerge_input2.txt"
            md_bedfile2 = mumerge(md_mumerge_input2, md_bedfile2, bed2=bed2, 
                                    label1=label1, label2=label2)
            md_pybedtool2 = BedTool(str(md_bedfile2))
            md_bedfile2 = tempdir / "md_bedfile2.mumerge.final.bed"
            md_pybedtool2.each(center_feature).each(extend_feature, size=largewindow).remove_invalid().saveas(md_bedfile2)

            # md_merged_bed1 = merge_bed(beds=bed1).each(featurefuncs.extend_fields, 4)
            # md_merged_bed2 = merge_bed(beds=bed2).each(featurefuncs.extend_fields, 4)
            # md_merged_bed1.each(center_feature).each(extend_feature, size=largewindow).remove_invalid().saveas(md_bedfile1)
            # md_merged_bed2.each(center_feature).each(extend_feature, size=largewindow).remove_invalid().saveas(md_bedfile2)
        

    #Merge all bed regions, for MD merge condition replicates
    elif method == 'mergeall':
        combined_file = tempdir / "combined_file.mergeall.bed"
        merged_bed = merge_bed(beds=bed1+bed2)
        # merged_bed.each(center_feature).each(extend_feature, size=largewindow).saveas(combined_file)
        merged_bed.remove_invalid().saveas(combined_file)
        if md:
            md_bedfile1 = tempdir / "md_bedfile1.merge.bed"
            md_bedfile2 = tempdir / "md_bedfile2.merge.bed"
            # md_merged_bed1 = merge_bed(beds=bed1).each(featurefuncs.extend_fields, 4).each(featurefuncs.rename, '1')
            # md_merged_bed2 = merge_bed(beds=bed2).each(featurefuncs.extend_fields, 4).each(featurefuncs.rename, '1')
            md_merged_bed1 = merge_bed(beds=bed1).each(featurefuncs.extend_fields, 4)
            md_merged_bed2 = merge_bed(beds=bed2).each(featurefuncs.extend_fields, 4)
            md_merged_bed1.each(center_feature).each(extend_feature, size=largewindow).remove_invalid().saveas(md_bedfile1)
            # md_merged_bed1.saveas(md_bedfile1)
            md_merged_bed2.each(center_feature).each(extend_feature, size=largewindow).remove_invalid().saveas(md_bedfile2)
            # md_merged_bed2.saveas(md_bedfile2)

    elif method == 'tfitclean':
        # combined_file = tfit_clean(beds=bed1+bed2, tempdir=tempdir)
        combined_file = tempdir / "combined_file.tfitclean.bed"
        size_cut = 200
        cleaned_bed = clean_bed(beds=bed1+bed2, size_cut=size_cut)
        # cleaned_bed.each(center_feature).each(extend_feature, size=largewindow).saveas(combined_file)
        cleaned_bed.remove_invalid().saveas(combined_file)
        if md:
            md_bedfile1 = tempdir / "md_bedfile1.clean.bed"
            md_bedfile2 = tempdir / "md_bedfile2.clean.bed"
            md_cleaned_bed1 = clean_bed(beds=bed1)
            md_cleaned_bed2 = clean_bed(beds=bed2)
            # md_cleaned_bed1.each(center_feature).each(extend_feature, size=largewindow).saveas(md_bedfile1)
            md_cleaned_bed1.saveas(combined_file)
            # md_cleaned_bed2.each(center_feature).each(extend_feature, size=largewindow).saveas(md_bedfile2)
            md_cleaned_bed2.saveas(combined_file)


    #Intersect all bed regions, for MD intersect condition replicates
    elif method == 'intersectall':
        combined_file = tempdir / 'combined_file.intersectall.bed'
        intersected_bed = intersect_bed(beds=bed1+bed2)
        # intersected_bed.each(center_feature).each(extend_feature, size=largewindow).saveas(combined_file)
        intersected_bed.remove_invalid().saveas(combined_file)
        if md:
            md_bedfile1 = tempdir / "md_bedfile1.intersect.bed"
            md_bedfile2 = tempdir / "md_bedfile2.intersect.bed"
            md_intersected_bed1 = intersect_bed(beds=bed1)
            md_intersected_bed2 = intersect_bed(beds=bed2)
            md_intersected_bed1.each(center_feature).each(extend_feature, size=largewindow).remove_invalid().saveas(md_bedfile1)
            # md_intersected_bed1.saveas(combined_file)
            md_intersected_bed2.each(center_feature).each(extend_feature, size=largewindow).remove_invalid().saveas(md_bedfile2)
            # md_intersected_bed2.saveas(combined_file)

    #Merge all regions, filter small regions. For MD perform this for each condition
    elif method == 'tfitremovesmall':
        # combined_file = tfit_remove_small(beds=bed1+bed2, tempdir=tempdir)
        size_cut = 200
        combined_file = tempdir / "combined_file.mergeallnosmall.bed"
        merged_bed = merge_bed(beds=bed1+bed2)
        # merged_bed.filter(lambda b: b.stop - b.start > size_cut).each(center_feature).each(extend_feature, size=largewindow).saveas(combined_file)
        merged_bed.filter(lambda b: b.stop - b.start > size_cut).saveas(combined_file)
        if md:
            md_bedfile1 = tempdir / "md_bedfile1.merge.bed"
            md_bedfile2 = tempdir / "md_bedfile2.merge.bed"
            md_merged_bed1 = merge_bed(beds=bed1)
            md_merged_bed2 = merge_bed(beds=bed2)
            # md_merged_bed1.filter(lambda b: b.stop - b.start > size_cut).each(center_feature).each(extend_feature, size=largewindow).saveas(md_bedfile1)
            md_merged_bed1.filter(lambda b: b.stop - b.start > size_cut).saveas(combined_file)
            # md_merged_bed2.filter(lambda b: b.stop - b.start > size_cut).each(center_feature).each(extend_feature, size=largewindow).saveas(md_bedfile2)
            md_merged_bed2.filter(lambda b: b.stop - b.start > size_cut).saveas(combined_file)

    #Intersect replicates, merge conditions. For MD intersect condition replicates
    elif method == 'intersect/merge':
        # combined_file = intersect_merge_bed(bed1=bed1, bed2=bed2, tempdir=tempdir)
        combined_file = tempdir / 'combined_file.intermerge.bed'
        intersected_bed1 = intersect_bed(beds=bed1)
        intersected_bed2 = intersect_bed(beds=bed2)
        merged_bed = intersected_bed1.cat(intersected_bed2).merge().sort()
        # merged_bed.each(center_feature).each(extend_feature, size=largewindow).saveas(combined_file)
        merged_bed.remove_invalid().saveas(combined_file)
        if md:
            md_bedfile1 = tempdir / "md_bedfile1.intersect.bed"
            md_bedfile2 = tempdir / "md_bedfile2.intersect.bed"
            md_intersected_bed1 = intersect_bed(beds=bed1).each(featurefuncs.extend_fields, 4).each(featurefuncs.rename, '1')
            md_intersected_bed2 = intersect_bed(beds=bed2).each(featurefuncs.extend_fields, 4).each(featurefuncs.rename, '1')
            md_intersected_bed1.each(center_feature).each(extend_feature, size=largewindow).remove_invalid().saveas(md_bedfile1)
            # md_intersected_bed1.saveas(md_bedfile1)
            md_intersected_bed2.each(center_feature).each(extend_feature, size=largewindow).remove_invalid().saveas(md_bedfile2)
            # md_intersected_bed2.saveas(md_bedfile2)
    
    else:
        raise exceptions.InputError("Error: COMBINE option not recognized.")

    #Check to make sure no files are empty
    if os.stat(combined_file).st_size == 0:
        raise exceptions.FileEmptyError("Error in COMBINE module. Resulting bed file is empty.")

    if md:
        if os.stat(md_bedfile1).st_size == 0 or os.stat(md_bedfile2).st_size == 0:
            raise exceptions.FileEmptyError("Error in COMBINE module. Resulting md bed file is empty.")
        if use_config:
            #Assign MD_BEDFILE variables in config
            config.vars['MD_BEDFILE1'] = md_bedfile1 
            config.vars['MD_BEDFILE2'] = md_bedfile2

    #Assign COMBINED_FILE variable in config
    if use_config:
        config.vars['COMBINED_FILE'] = combined_file
    
    #Record time, print
    total_time = time.time() - start_time
    if use_config:
        config.vars['COMBINEtime'] = total_time
    print("done in: " + str(datetime.timedelta(seconds=int(total_time))), 
            ". Processing", len(combined_file.read_text().split('\n')), "regions", file=sys.stderr)

    if debug:
        multiprocess.current_mem_usage(jobid)

#Functions
#==============================================================================
def merge_bed(beds=None):
    '''Concatenates, sorts, and merges (bedtools) a list of bed files. Outputs 
        into the tempdir directory created by TFEA

    Parameters
    ----------
    beds : list or array
        full paths to bed files (python Path objects from pathlib)

    Returns
    -------
    merged_bed : BedTool object 
        resulting merged bed object 
    '''
    parent_bed = BedTool(str(beds[0]))
    for bed in beds[1:]:
        parent_bed = parent_bed.cat(str(bed))
    merged_bed = parent_bed.sort().merge().sort()

    return merged_bed

#==============================================================================
def intersect_bed(beds=None):
    '''Intersects, sorts, and merges a list of bed files

    Parameters
    ----------
    beds : list or array
        full paths to bed files (python Path objects from pathlib)

    Returns
    -------
    intersected_bed : BedTool object 
        resulting intersected bed object 
    '''
    parent_bed = BedTool(beds[0].as_posix())
    for bed in beds[1:]:
        parent_bed = parent_bed.intersect(BedTool(bed.as_posix()))
    intersected_bed = parent_bed.sort()

    return intersected_bed

#==============================================================================
def clean_bed(beds=None, size_cut=None):
    '''This function separates a list of beds into small and large regions based
        on a size_cut, intersects small regions, merges large regions, and then
        merges the small and large regions.

    Parameters
    ----------
    beds : list or array
        full paths to bed files (python Path objects from pathlib)

    size_cut : int
        cutoff value to separate large and small regions

    Returns
    -------
    clean_bed : BedTool object 
        resulting clean bed object 
    '''
    small_regions = list()
    large_regions = list()
    for bed in beds:
        bed = BedTool(bed)
        small_regions.append(bed.filter(lambda b: b.stop - b.start < size_cut))
        large_regions.append(bed.filter(lambda b: b.stop - b.start > size_cut))
    small_bed = intersect_bed(beds=small_regions)
    large_bed = merge_bed(beds=large_regions)
    clean_bed = large_bed.cat(small_bed).merge().sort()

    return clean_bed

#==============================================================================
def mumerge(input_file, output_basename, bed1=[], bed2=[], label1=None, 
            label2=None,
            mumerge_path=Path(__file__).absolute().parent / 'mumerge.py'):
    '''This function runs MuMerge, a script written by Jacob T. Stanley that 
        merges a list of bed files in a probabilistic way.
        
    Parameters
    ----------
    input_file: path to .txt file
        A .txt file formatted according to MuMerge specifications. From doc:
            Input file containing bedfiles, sample ID's, and replicate groupings. Input
            file (indicated by the '-i' flag) should be of the following (tab delimited)
            format:

            #file   sampid  group
            /full/file/path/filename1.bed   sampid1 A
            /full/file/path/filename2.bed   sampid2 B
            ...
            
            Header line indicated by '#' character must be included and fields must
            follow the same order as non-header lines. The order of subsequent lines does
            matter. 'group' identifiers should group files that are technical/biological
            replicates. Different experimental conditions should recieve different 'group'
            identifiers. The 'group' identifier can be of type 'int' or 'str'. If 'sampid'
            is not specified, then default sample ID's will be used.
            
    output_basename: Path to output file without file extension
        From doc:
            Output file basename (full path, sans extension).
            WARNING: will overwrite any existing file)'''
    with open(input_file, 'w') as F:
            F.write("#file\tsampid\tgroup\n")
            for i,bedpath in enumerate(bed1, 1):
                F.write(f'{bedpath}\t{label1}{i}\t{label1}\n')
            for i,bedpath in enumerate(bed2, 1):
                F.write(f'{bedpath}\t{label2}{i}\t{label2}\n')
        
    mumerge_command = ['python3', mumerge_path, '-i', input_file, '-o', output_basename]
    try:
        subprocess.check_output(mumerge_command, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        raise exceptions.SubprocessError(e.stderr.decode())
    combined_file = Path(str(output_basename) + '_MUMERGE.bed')

    return combined_file

#Pybedtools each() functions
#==============================================================================
def center_feature(feature):
    '''A function for use with pybedtools.BedTool.each() that centers each
        feature

    Parameters
    ----------
    feature : Interval object
        an Interval object from pybedtools library
    
    Returns
    -------
    feature : Interval object
        the modified feature
    '''
    center = int((feature.start + feature.stop)/2)
    feature.start = center
    feature.stop = center + 1

    return feature

#==============================================================================
def extend_feature(feature, size=0):
    '''A function for use with pybedtools.BedTool.each() that extends each
        feature by size
    
    Parameters
    ----------
    feature : Interval object
        an Interval object from pybedtools library

    size : int
        size to extend feature on each size

    Returns
    -------
    feature : Interval object
        the modified feature
    '''
    feature.start = feature.start - size
    if feature.start < 0:
        feature.start = 0
    feature.stop = feature.stop + size
    return feature

    