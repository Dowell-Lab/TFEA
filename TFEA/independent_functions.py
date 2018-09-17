#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''This file contains a list of independent functions that do not call other 
    functions
'''
#==============================================================================
__author__ = 'Jonathan D. Rubin and Rutendo Sigauke'
__credits__ = ['Jonathan D. Rubin', 'Rutendo Sigauke', 'Jacob Stanley',
                'Robin Dowell']
__maintainer__ = 'Jonathan D. Rubin'
__email__ = 'Jonathan.Rubin@colorado.edu'
#==============================================================================
import matplotlib
matplotlib.use('Agg')
import os
import math
import datetime
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.cm as cm
from matplotlib import gridspec
from scipy.stats import norm 
import numpy as np
import config
#==============================================================================
#Functions
#==============================================================================
def make_out_directories(dirs=False,config=None):
    '''Creates output directories in a user-specified location where all TFEA 
        outputs will go.

    Parameters
    ----------
    dirs : boolean
        determines whether output folders will be created or not (default: 
        False)
        
    config : dict
        a configparser object that contains variables within the config file

    Returns
    -------
    output : string 
        full path to the parent output directory

    tempdir : string
        full path to the temporary directory where files are stored

    figuredir : string
        full path to the directory containing figures and plots
        
    e_and_o : string
        full path to the directory that stores stdout and stderr files
    '''
    #Output directory
    output = config['DATA']['OUTPUT'].strip("'")
    label1 = config['DATA']['LABEL1'].strip("'")
    label2 = config['DATA']['LABEL2'].strip("'")
    outfoldername = 'TFEA_'+label1+'-'+label2+'_'
    if dirs:
        if not os.path.isdir(output + outfoldername + '0/'):
            output = output + outfoldername + '0/'
            os.makedirs(output)
        else:
            outputfolders = list()
            for folder in os.listdir(output):
                if outfoldername in folder:
                    outputfolders.append(int(folder.split('_')[-1]))
            output = output + outfoldername + str(max(outputfolders)+1) + '/'
            os.makedirs(output)
    else:
        outputfolders = list()
        for folder in os.listdir(output):
            if outfoldername in folder:
                outputfolders.append(int(folder.split('_')[-1]))
        output = output + outfoldername + str(max(outputfolders)) + '/' 


    #Temporary files will go in this directory
    tempdir = output + 'temp_files/'
    if dirs:
        if not os.path.isdir(tempdir):
            os.makedirs(tempdir)

    #Error and out files will go in this directory
    e_and_o = output + 'e_and_o/'
    if dirs:
        if not os.path.isdir(e_and_o):
            os.makedirs(e_and_o)


    #Directory where plots used in html file will be stored.
    figuredir = output + 'plots/'
    if dirs:
        if not os.path.isdir(figuredir):
            os.makedirs(figuredir)

    return output,tempdir,figuredir,e_and_o
#==============================================================================

#==============================================================================
def parse_config(srcdirectory=str(), config=str(), output=str(), tempdir=str(),
                figuredir=str()):
    '''Creates the config.py file which is used in many aspects of TFEA. Within
        this config.py file, it writes all variables provided in the config
        parameter and also writes output, tempdir, and figuredir full paths

    Parameters
    ----------
    srcdirectory : string
        full path to TFEA source directory
        
    config : dict
        a configparser object that contains variables within the config file

    output : string 
        full path to the parent output directory

    tempdir : string
        full path to the temporary directory where files are stored

    figuredir : string
        full path to the directory containing figures and plots

    Returns
    -------
    None
    '''
    with open(os.path.join(srcdirectory,'config.py'),'w') as outfile:
        for key in config:
            for item in config[key]:
                outfile.write(item.upper()+'='+config[key][item]+'\n')

        outfile.write('OUTPUTDIR="'+output+'"\n')
        outfile.write('TEMPDIR="'+tempdir+'"\n')
        outfile.write('FIGUREDIR="'+figuredir+'"\n')

        #Path to count file. Can be changed if using your own count file.
        #Generated in count_reads function
        count_file = tempdir + "count_file.header.bed"
        outfile.write('COUNT_FILE="'+count_file+'"\n')
        
        #Path to ranked center file. 
        ranked_center_file = tempdir + "ranked_file.center.bed"
        outfile.write('RANKED_CENTER_FILE="'+ranked_center_file+'"\n')
#==============================================================================

#==============================================================================
def verify_config_file():
    '''Verifies that all necessary variables are present within the inputted 
        config file and that they are the correct variable types.
    
    Parameters
    ----------
    None

    Returns
    -------
    None

    Raises
    ------
    TypeError
        When variable type does not match required type

    NameError
        When variable is not found within config
    '''
    try:
        if type(config.COMBINE) != bool:
            raise TypeError('COMBINE variable must be a boolean. This switch \
                            determines whether TFEA merges bed files within \
                            BED input.')
    except NameError:
        raise NameError('COMBINE variable not found in config.ini file. This \
                            switch determines whether TFEA merges bed files \
                            within BED input.')

    try:
        if type(config.COUNT) != bool:
            raise TypeError('COUNT variable must be a boolean. This switch \
                                determines whether TFEA performs read \
                                counting over BED regions.')
    except NameError:
        raise NameError('COUNT variable not found in config.ini file. This \
                            switch determines whether TFEA performs read \
                            counting over BED regions.')

    try:
        if type(config.DESEQ) != bool:
            raise TypeError('DESEQ variable must be a boolean. This switch \
                                determines whether TFEA performs DE-Seq \
                                analysis on counted BED regions.')
    except NameError:
        raise NameError('DESEQ variable not found in config.ini file. This \
                            switch determines whether TFEA performs DE-Seq \
                            analysis on counted BED regions.')

    try:
        if type(config.CALCULATE) != bool:
            raise TypeError('CALCULATE variable must be a boolean. This \
                                switch determines whether TFEA performs its \
                                standard enrichment score calculation and \
                                plotting.')
    except NameError:
        raise NameError('CALCULATE variable not found in config.ini file. \
                            This switch determines whether TFEA performs its \
                            standard enrichment score calculation and \
                            plotting.')

    try:
        if type(config.POOL) != bool:
            raise TypeError('POOL variable must be a boolean. This switch \
                                determines whether TFEA runs the analysis \
                                in parallel using the multiprocessing library \
                                in python.')
    except NameError:
        raise NameError('POOL variable not found in config.ini file. This \
                            switch determines whether TFEA runs the analysis \
                            in parallel using the multiprocessing library in \
                            python.')

    try:
        if type(config.SINGLEMOTIF) != bool and type(config.SINGLEMOTIF) != str:
            raise TypeError('SINGLEMOTIF variable must be a boolean or \
                                string. This switch determines whether TFEA \
                                performs its analysis on a single motif or \
                                all. If not False, set to a string matching a \
                                motif name.')
    except NameError:
        raise NameError('SINGLEMOTIF variable not found in config.ini file. \
                            This switch determines whether TFEA performs its \
                            analysis on a single motif or all. If not False, \
                            set to a string matching a motif name.')

    try:
        if type(config.FIMO) != bool:
            raise TypeError('FIMO variable must be a boolean. This switch \
                                determines whether TFEA uses FIMO to get \
                                motif hits or whether a database of motif hit \
                                calls (bed format) is used.')
    except NameError:
        raise NameError('FIMO variable not found in config.ini file. This \
                            switch determines whether TFEA uses FIMO to get \
                            motif hits or whether a database of motif hit \
                            calls (bed format) is used.')

    try:
        if type(config.TEMP) != bool:
            raise TypeError('TEMP variable must be a boolean. This switch \
                                determines whether TFEA saves large temporary \
                                files. If True, temporary files will be \
                                stored in the temp_files directory. Warning: \
                                There will be many large files.')
    except NameError:
        raise NameError('TEMP variable not found in config.ini file. This \
                            switch determines whether TFEA saves large \
                            temporary files. If True, temporary files will be \
                            stored in the temp_files directory. Warning: \
                            There will be many large files.')

    try:
        if type(config.PLOT) != bool:
            raise TypeError('PLOT variable must be a boolean. This switch \
                                determines whether TFEA outputs plots for \
                                all motifs provided or just significant ones. \
                                Warning: Setting this to True will slow down \
                                TFEA and create large output folders.')
    except NameError:
        raise NameError('PLOT variable not found in config.ini file. This \
                            switch determines whether TFEA outputs plots for \
                            all motifs provided or just significant ones. \
                            Warning: Setting this to True will slow down TFEA \
                            and create large output folders.')

    try:
        if type(config.OUTPUT) != str:
            raise TypeError('OUTPUT variable must be a string. Determines \
                                where TFEA stores output.')
    except NameError:
        raise NameError('OUTPUT variable not found in config.ini file. \
                            Determines where TFEA stores output.')

    try:
        if type(config.BEDS) != list:
            raise TypeError('BED variable must be a list. Input a list of \
                                regions (bed-format) to perform analysis \
                                over. If merging not desired, simply input a \
                                single BED file within a python list.')
    except NameError:
        raise NameError('BED variable not found in config.ini file. Input a \
                            list of regions (bed-format) to perform analysis \
                            over. If merging not desired, simply input a \
                            single BED file within a python list.')

    try:
        if type(config.BAM1) != list:
            raise TypeError('BAM1 variable must be a list. Input a list of \
                                BAM files to obtain read depth and coverage \
                                over regions of interest. One or multiple \
                                bams can be specified but they must be within \
                                a python list.')
    except NameError:
        raise NameError('BAM1 variable not found in config.ini file. Input a \
                            list of BAM files to obtain read depth and \
                            coverage over regions of interest. One or \
                            multiple bams can be specified but they must be \
                            within a python list.')

    try:
        if type(config.LABEL1) != str:
            raise TypeError('LABEL1 variable must be a string. Define a \
                                treatment or condition to label BAM1 files. \
                                Used in plotting and output folder naming.')
    except NameError:
        raise NameError('LABEL1 variable not found in config.ini file. Define \
                            a treatment or condition to label BAM1 files. \
                            Used in plotting and output folder naming.')

    try:
        if type(config.BAM2) != list:
            raise TypeError('BAM2 variable must be a list. Input a list of \
                                BAM files to obtain read depth and coverage \
                                over regions of interest. One or multiple \
                                bams can be specified but they must be within \
                                a python list.')
    except NameError:
        raise NameError('BAM2 variable not found in config.ini file. Input a \
                            list of BAM files to obtain read depth and \
                            coverage over regions of interest. One or \
                            multiple bams can be specified but they must be \
                            within a python list.')

    try:
        if type(config.LABEL2) != str:
            raise TypeError('LABEL2 variable must be a string. Define a \
                                treatment or condition to label BAM1 files. \
                                Used in plotting and output folder naming.')
    except NameError:
        raise NameError('LABEL2 variable not found in config.ini file. Define \
                            a treatment or condition to label BAM1 files. \
                            Used in plotting and output folder naming.')

    try:
        if type(config.PADJCUTOFF) != float:
            raise TypeError('PADJCUTOFF variable must be a float. Provide a \
                                p-adjusted cutoff value to determine which \
                                TFs are plotted and called as significant.')
    except NameError:
        raise NameError('PADJCUTOFF variable not found in config.ini file. \
                                Provide a p-adjusted cutoff value to \
                                determine which TFs are plotted and called as \
                                significant.')

    try:
        if type(config.LARGEWINDOW) != float:
            raise TypeError('LARGEWINDOW variable must be a float. Provide a \
                                larger window to be used by TFEA for plotting \
                                and GC-content calculation.')
    except NameError:
        raise NameError('LARGEWINDOW variable not found in config.ini file. \
                                Provide a larger window to be used by TFEA \
                                for plotting and GC-content calculation.')

    try:
        if type(config.SMALLWINDOW) != float:
            raise TypeError('SMALLWINDOW variable must be a float. Provide a \
                                smaller window to be used by TFEA for \
                                plotting.')
    except NameError:
        raise NameError('SMALLWINDOW variable not found in config.ini file. \
                            Provide a smaller window to be used by TFEA for \
                            plotting.')

    try:
        if type(config.MOTIF_HITS) != str:
            raise TypeError('MOTIF_HITS variable must be a string. Provide \
                                the full path to a folder containing bed \
                                files with motif hits across the genome. \
                                These can be generated using TFEAs compile \
                                module.')
    except NameError:
        raise NameError('MOTIF_HITS variable not found in config.ini file. \
                            Provide the full path to a folder containing bed \
                            files with motif hits across the genome. These \
                            can be generated using TFEAs compile module.')

    try:
        if type(config.GENOMEFASTA) != str:
            raise TypeError('GENOMEFASTA variable must be a string. Provide \
                                the full path to a fasta file containing a \
                                genome of interest to perform motif scannning \
                                over.')
    except NameError:
        raise NameError('GENOMEFASTA variable not found in config.ini file. \
                            Provide the full path to a fasta file containing \
                            a genome of interest to perform motif scannning \
                            over.')

    try:
        if type(config.MOTIFDATABASE) != str:
            raise TypeError('MOTIFDATABASE variable must be a string. Provide \
                                the full path to a database of motifs in meme \
                                format.')
    except NameError:
        raise NameError('MOTIFDATABASE variable not found in config.ini file. \
                            Provide the full path to a database of motifs in \
                            meme format.')

    try:
        if type(config.LOGOS) != str:
            raise TypeError('LOGOS variable must be a string. Provide the \
                                full path to a directory containing meme \
                                formatted motif logos whose name correspond \
                                motifs within MOTIFDATABASE.')
    except NameError:
        raise NameError('LOGOS variable not found in config.ini file. Provide \
                                the full path to a directory containing meme \
                                formatted motif logos whose name correspond \
                                motifs within MOTIFDATABASE.')

    print "Config file verified, all inputs present and correct type."
#==============================================================================

#==============================================================================
def verify_config_object(config=object()):
    '''Verifies the names and values of variables within a configParser object
        to make sure that all necessary variables are present and they are the
        correct type
    
    Parameters
    ----------
    config : configParser object
        a configParser object storing the variables contained within the 
        user-inputted config.ini file
    
    Returns
    -------
    config_dict : dictionary
        a dictionary with keys corresponding to variable names and values to
        user-inputted values

    Raises
    ------
    TypeError
        When variable type does not match required type

    NameError
        When variable is not found within config
    '''
    #Initiallize the output dictionary
    config_dict = dict()

    print config

    #get a list of (name, value) pairs for each option in the given section
    config_items = list()
    config_items = [config_items + config.items(section) for section in config]

    print config_items

    #pull out only names, give them all uppercase so the user input can be 
    #case-insensitive
    names = [name.upper() for (name,value) in config_items]

    #pull out corresponding value pair for each variable name
    values = [value for (name,value) in config_items]

    #Try to find the index of a given variable
    try:
        #If variable name found, verify it's type
        combine = values[names.index('COMBINE')]
        if type(combine) != bool:
            #If it's type is incorrect raise a TypeError
            raise TypeError('COMBINE variable must be a boolean. This switch \
                            determines whether TFEA merges bed files within \
                            BED input.')
        else:
            config_dict['COMBINE'] = combine
    #Catch if the variable name is not within the name list
    except ValueError:
        #Raise a NameError
        raise NameError('COMBINE variable not found in config.ini file. This \
                            switch determines whether TFEA merges bed files \
                            within BED input.')

    try:
        count = values[names.index('COUNT')]
        if type(count) != bool:
            raise TypeError('COUNT variable must be a boolean. This switch \
                                determines whether TFEA performs read \
                                counting over BED regions.')
        else:
            config_dict['COUNT'] = count
    except ValueError:
        raise NameError('COUNT variable not found in config.ini file. This \
                            switch determines whether TFEA performs read \
                            counting over BED regions.')

    try:
        deseq = values[names.index('DESEQ')]
        if type(deseq) != bool:
            raise TypeError('DESEQ variable must be a boolean. This switch \
                                determines whether TFEA performs DE-Seq \
                                analysis on counted BED regions.')
        else:
            config_dict['DESEQ'] = deseq
    except ValueError:
        raise NameError('DESEQ variable not found in config.ini file. This \
                            switch determines whether TFEA performs DE-Seq \
                            analysis on counted BED regions.')

    try:
        calculate = values[names.index('CALCULATE')]
        if type(calculate) != bool:
            raise TypeError('CALCULATE variable must be a boolean. This \
                                switch determines whether TFEA performs its \
                                standard enrichment score calculation and \
                                plotting.')
        else:
            config_dict['CALCULATE']
    except ValueError:
        raise NameError('CALCULATE variable not found in config.ini file. \
                            This switch determines whether TFEA performs its \
                            standard enrichment score calculation and \
                            plotting.')

    try:
        pool = values[names.index('POOL')]
        if type(pool) != bool:
            raise TypeError('POOL variable must be a boolean. This switch \
                                determines whether TFEA runs the analysis \
                                in parallel using the multiprocessing library \
                                in python.')
        else:
            config_dict['POOL'] = pool
    except ValueError:
        raise NameError('POOL variable not found in config.ini file. This \
                            switch determines whether TFEA runs the analysis \
                            in parallel using the multiprocessing library in \
                            python.')

    try:
        singlemotif = values[names.index('SINGLEMOTIF')]
        if type(singlemotif) != bool and type(singlemotif) != str:
            raise TypeError('SINGLEMOTIF variable must be a boolean or \
                                string. This switch determines whether TFEA \
                                performs its analysis on a single motif or \
                                all. If not False, set to a string matching a \
                                motif name.')
        else:
            config_dict['SINGLEMOTIF'] = singlemotif
    except ValueError:
        raise NameError('SINGLEMOTIF variable not found in config.ini file. \
                            This switch determines whether TFEA performs its \
                            analysis on a single motif or all. If not False, \
                            set to a string matching a motif name.')

    try:
        fimo = values[names.index('FIMO')]
        if type(fimo) != bool:
            raise TypeError('FIMO variable must be a boolean. This switch \
                                determines whether TFEA uses FIMO to get \
                                motif hits or whether a database of motif hit \
                                calls (bed format) is used.')
        else:
            config_dict['FIMO'] = fimo
    except ValueError:
        raise NameError('FIMO variable not found in config.ini file. This \
                            switch determines whether TFEA uses FIMO to get \
                            motif hits or whether a database of motif hit \
                            calls (bed format) is used.')

    try:
        temp = values[names.index('TEMP')]
        if type(temp) != bool:
            raise TypeError('TEMP variable must be a boolean. This switch \
                                determines whether TFEA saves large temporary \
                                files. If True, temporary files will be \
                                stored in the temp_files directory. Warning: \
                                There will be many large files.')
        else:
            config_dict['TEMP'] = temp
    except ValueError:
        raise NameError('TEMP variable not found in config.ini file. This \
                            switch determines whether TFEA saves large \
                            temporary files. If True, temporary files will be \
                            stored in the temp_files directory. Warning: \
                            There will be many large files.')

    try:
        plot = values[names.index('PLOT')]
        if type(plot) != bool:
            raise TypeError('PLOT variable must be a boolean. This switch \
                                determines whether TFEA outputs plots for \
                                all motifs provided or just significant ones. \
                                Warning: Setting this to True will slow down \
                                TFEA and create large output folders.')
        else:
            config_dict['PLOT'] = plot
    except ValueError:
        raise NameError('PLOT variable not found in config.ini file. This \
                            switch determines whether TFEA outputs plots for \
                            all motifs provided or just significant ones. \
                            Warning: Setting this to True will slow down TFEA \
                            and create large output folders.')

    try:
        output = values[names.index('OUTPUT')] 
        if type(output) != str:
            raise TypeError('OUTPUT variable must be a string. Determines \
                                where TFEA stores output.')
        else:
            config_dict['OUTPUT'] = output
    except ValueError:
        raise NameError('OUTPUT variable not found in config.ini file. \
                            Determines where TFEA stores output.')

    try:
        beds = values[names.index('BEDS')]
        if type(beds) != list:
            raise TypeError('BED variable must be a list. Input a list of \
                                regions (bed-format) to perform analysis \
                                over. If merging not desired, simply input a \
                                single BED file within a python list.')
        else:
            config_dict['BEDS'] = beds
    except ValueError:
        raise NameError('BED variable not found in config.ini file. Input a \
                            list of regions (bed-format) to perform analysis \
                            over. If merging not desired, simply input a \
                            single BED file within a python list.')

    try:
        bam1 = values[names.index('BAM1')]
        if type(bam1) != list:
            raise TypeError('BAM1 variable must be a list. Input a list of \
                                BAM files to obtain read depth and coverage \
                                over regions of interest. One or multiple \
                                bams can be specified but they must be within \
                                a python list.')
        else:
            config_dict['BAM1'] = bam1
    except ValueError:
        raise NameError('BAM1 variable not found in config.ini file. Input a \
                            list of BAM files to obtain read depth and \
                            coverage over regions of interest. One or \
                            multiple bams can be specified but they must be \
                            within a python list.')

    try:
        label1 = values[names.index('LABEL1')] 
        if type(label1) != str:
            raise TypeError('LABEL1 variable must be a string. Define a \
                                treatment or condition to label BAM1 files. \
                                Used in plotting and output folder naming.')
        else:
            config_dict['LABEL1'] = label1
    except ValueError:
        raise NameError('LABEL1 variable not found in config.ini file. Define \
                            a treatment or condition to label BAM1 files. \
                            Used in plotting and output folder naming.')

    try:
        bam2 = values[names.index('BAM2')]
        if type(bam2) != list:
            raise TypeError('BAM2 variable must be a list. Input a list of \
                                BAM files to obtain read depth and coverage \
                                over regions of interest. One or multiple \
                                bams can be specified but they must be within \
                                a python list.')
        else:
            config_dict['BAM2'] = bam2
    except ValueError:
        raise NameError('BAM2 variable not found in config.ini file. Input a \
                            list of BAM files to obtain read depth and \
                            coverage over regions of interest. One or \
                            multiple bams can be specified but they must be \
                            within a python list.')

    try:
        label2 = values[names.index('LABEL2')] 
        if type(label2) != str:
            raise TypeError('LABEL2 variable must be a string. Define a \
                                treatment or condition to label BAM1 files. \
                                Used in plotting and output folder naming.')
        else:
            config_dict['LABEL2'] = label2
    except ValueError:
        raise NameError('LABEL2 variable not found in config.ini file. Define \
                            a treatment or condition to label BAM1 files. \
                            Used in plotting and output folder naming.')

    try:
        padj_cutoff = values[names.index('PADJCUTOFF')] 
        if type(values[names.index('PADJCUTOFF')]) != float:
            raise TypeError('PADJCUTOFF variable must be a float. Provide a \
                                p-adjusted cutoff value to determine which \
                                TFs are plotted and called as significant.')
        else:
            config_dict['PADJCUTOFF'] = padj_cutoff
    except ValueError:
        raise NameError('PADJCUTOFF variable not found in config.ini file. \
                                Provide a p-adjusted cutoff value to \
                                determine which TFs are plotted and called as \
                                significant.')

    try:
        largewindow = values[names.index('LARGEWINDOW')] 
        if type(values[names.index('LARGEWINDOW')]) != float:
            raise TypeError('LARGEWINDOW variable must be a float. Provide a \
                                larger window to be used by TFEA for plotting \
                                and GC-content calculation.')
        else:
            config_dict['LARGEWINDOW'] = largewindow
    except ValueError:
        raise NameError('LARGEWINDOW variable not found in config.ini file. \
                                Provide a larger window to be used by TFEA \
                                for plotting and GC-content calculation.')

    try:
        smallwindow = values[names.index('SMALLWINDOW')] 
        if type(values[names.index('SMALLWINDOW')]) != float:
            raise TypeError('SMALLWINDOW variable must be a float. Provide a \
                                smaller window to be used by TFEA for \
                                plotting.')
        else:
            config_dict['SMALLWINDOW'] = smallwindow
    except ValueError:
        raise NameError('SMALLWINDOW variable not found in config.ini file. \
                            Provide a smaller window to be used by TFEA for \
                            plotting.')

    try:
        motif_hits = values[names.index('MOTIF_HITS')] 
        if type(values[names.index('MOTIF_HITS')]) != str:
            raise TypeError('MOTIF_HITS variable must be a string. Provide \
                                the full path to a folder containing bed \
                                files with motif hits across the genome. \
                                These can be generated using TFEAs compile \
                                module.')
        else:
            config_dict['MOTIF_HITS'] = motif_hits
    except ValueError:
        raise NameError('MOTIF_HITS variable not found in config.ini file. \
                            Provide the full path to a folder containing bed \
                            files with motif hits across the genome. These \
                            can be generated using TFEAs compile module.')

    try:
        genomefasta = values[names.index('GENOMEFASTA')] 
        if type(values[names.index('GENOMEFASTA')]) != str:
            raise TypeError('GENOMEFASTA variable must be a string. Provide \
                                the full path to a fasta file containing a \
                                genome of interest to perform motif scannning \
                                over.')
        else:
            config_dict['GENOMEFASTA'] = genomefasta
    except ValueError:
        raise NameError('GENOMEFASTA variable not found in config.ini file. \
                            Provide the full path to a fasta file containing \
                            a genome of interest to perform motif scannning \
                            over.')

    try:
        motifdatabase = values[names.index('MOTIFDATABASE')] 
        if type(values[names.index('MOTIFDATABASE')]) != str:
            raise TypeError('MOTIFDATABASE variable must be a string. Provide \
                                the full path to a database of motifs in meme \
                                format.')
        else:
            config_dict['MOTIFDATABASE'] = motifdatabase
    except ValueError:
        raise NameError('MOTIFDATABASE variable not found in config.ini file. \
                            Provide the full path to a database of motifs in \
                            meme format.')

    try:
        logos = values[names.index('LOGOS')] 
        if type(values[names.index('LOGOS')]) != str:
            raise TypeError('LOGOS variable must be a string. Provide the \
                                full path to a directory containing meme \
                                formatted motif logos whose name correspond \
                                motifs within MOTIFDATABASE.')
        else:
            config_dict['LOGOS'] = logos
    except ValueError:
        raise NameError('LOGOS variable not found in config.ini file. Provide \
                                the full path to a directory containing meme \
                                formatted motif logos whose name correspond \
                                motifs within MOTIFDATABASE.')

    print "Config file verified, all inputs present and correct type."
    return config_dict
#==============================================================================

#==============================================================================
def sbatch_submit(srcdirectory=str(),configpath=str(),script=str(),email=str(),
                    config=None):
    '''Submits an sbatch job using the configuration provided to the script 
        variable that runs TFEA with options specified within a give config 
        file

    Parameters
    ----------
    srcdirectory : string
        full path to TFEA source directory
        
    configpath : string
        full path to the config file

    script : string
        full path to the sbatch script containing configuration options

    email : string
        user e-mail to receive sbatch job notifications

    config : dict
        a configparser object

    Returns
    -------
    None
    '''
    #First make output directories, only save path to the e_and_o folder
    _, _, _, e_and_o = make_out_directories(dirs=True, config=config)

    #Submit the sbatch job
    os.system("sbatch --error=" + e_and_o + "%x.err --output=" + e_and_o 
                + "%x.out --mail-user="+email+" --export=src="+srcdirectory
                + ",config=" +configpath+ " " + script)
#==============================================================================

#==============================================================================
def merge_bed(beds=list(),tempdir=str()):
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
    os.system("cat " + " ".join(beds) + " > " + tempdir + "combined_input.bed")

    os.system("sort -k1,1 -k2,2n " + tempdir + "combined_input.bed > " 
                + tempdir + "combined_input.sorted.bed")

    os.system("bedtools merge -i " + tempdir + "combined_input.sorted.bed > " 
                + tempdir + "combined_input.merge.bed")

    combined_input_merged_bed = tempdir + "combined_input.merge.bed"

    return combined_input_merged_bed
#==============================================================================

#==============================================================================
def getfasta(bedfile=str(), genomefasta=str(), tempdir=str()):
    '''Converts a bed file to a fasta file using bedtools. Outputs into the 
        tempdir directory created by TFEA.

    Parameters
    ----------
    bedfile : string
        full path to a bed file

    genomefasta : string
        full path to a fasta file for the genome of interest
        
    tempdir : string
        full path to temp directory in output directory (created by TFEA)

    Returns
    -------
    ranked_file_fasta : string 
        full path to a fasta file containing the inputted bed file regions in 
        fasta format 
    '''
    os.system("bedtools getfasta -name -fi "+genomefasta+" -bed "+bedfile
                + " -fo " + tempdir + "ranked_file.fullregions.fa")

    ranked_file_fasta = tempdir + "ranked_file.fullregions.fa"

    return ranked_file_fasta
#==============================================================================

#==============================================================================
def get_bgfile(fastafile=str(), tempdir=str()):
    '''Obtains a zero order markov background model (used in FIMO) from a fasta
        file. Outputs into the tempdir directory created by TFEA.

    Parameters
    ----------
    bedfile : string
        full path to a bed file

    genomefasta : string
        full path to a fasta file for the genome of interest
        
    tempdir : string
        full path to temp directory in output directory (created by TFEA)

    Returns
    -------
    markov_background : string 
        full path to a meme-formatted zero order markov background file 
        (txt file)
    '''
    os.system("fasta-get-markov " + fastafile + " " + tempdir
                + "markov_background.txt")

    markov_background = tempdir + "markov_background.txt"

    return markov_background
#==============================================================================

#==============================================================================
def get_regions(tempdir=str(), ranked_center_file=str(), largewindow=float()):
    '''Takes in a bed file that contains regions centered on the user-inputted 
        bed files and outputs a 'full regions' bed file which simply adds a 
        user-defined window to either side of these centered regions.

    Parameters
    ----------
    tempdir : string
        full path to temp directory in output directory (created by TFEA)

    ranked_center_file : string
        full path to a bed file with centered coordinates
        
    largewindow : float
        half the desired window size

    Returns
    -------
    ranked_full_regions : string 
        full path to a bed file with full regions to be compared with TFEA
    '''
    outfile = open(tempdir+'ranked_file.fullregions.bed','w')
    with open(ranked_center_file) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            chrom,start,stop = line[:3]
            start = str(int(start)-int(largewindow))
            stop = str(int(stop)+int(largewindow))
            pval,fc,rank = line[3:]
            name = ','.join([rank,pval,fc])
            outfile.write('\t'.join([chrom,start,stop,name]) + '\n')

    ranked_full_regions = tempdir + 'ranked_file.fullregions.bed'
    return ranked_full_regions
#==============================================================================

#==============================================================================
def count_reads(bedfile=str(), bam1=list(), bam2=list(), tempdir=str(), 
                label1=str(), label2=str()):
    '''Counts reads across regions in a given bed file using bam files inputted
        by a user

    Parameters
    ----------
    bedfile : string
        full path to a bed file containing full regions of interest which will 
        be counted using bedtools multicov

    bam1 : list or array
        a list of full paths to bam files pertaining to a single condition 
        (i.e. replicates of a single treatment)

    bam2 : list or array
        a list of full paths to bam files pertaining to a single condition 
        (i.e. replicates of a single treatment)

    tempdir : string
        full path to temp directory in output directory (created by TFEA)

    label1 : string
        the name of the treatment or condition corresponding to bam1 list

    label2 : string
        the name of the treatment or condition corresponding to bam2 list

    Returns
    -------
    None
    '''
    #This os.system call runs bedtools multicov to count reads in all specified
    #BAMs for given regions in BED
    os.system("bedtools multicov -bams " + " ".join(bam1) + " " 
                + " ".join(bam2) + " -bed " + bedfile + " > " + tempdir 
                + "count_file.bed")

    #This section adds a header to the count_file and reformats it to remove 
    #excess information and add a column with the region for later use
    outfile = open(tempdir + "count_file.header.bed",'w')
    outfile.write("chrom\tstart\tstop\tregion\t" 
                    + '\t'.join([label1]*len(bam1)) + "\t" 
                    + '\t'.join([label2]*len(bam2)) + "\n")

    with open(tempdir + "count_file.bed") as F:
        for line in F:
            line = line.strip('\n').split('\t')
            chrom,start,stop = line[:3]
            counts = line[-(len(bam1)+len(bam2)):]
            outfile.write('\t'.join([chrom,start,stop]) + "\t" 
                            + chrom + ":" + start + "-" + stop + "\t"
                            + '\t'.join(counts) + "\n")
#==============================================================================

#==============================================================================
def write_deseq_script(bam1=list(), bam2=list(), tempdir=str(), 
                        count_file=str(), label1=str(), label2=str()):
    '''Writes an R script within the tempdir directory in TFEA output to run 
        either DE-Seq or DE-Seq2 depending on the number of user-inputted 
        replicates.

    Parameters
    ----------
    bedfile : string
        full path to a bed file containing full regions of interest which will
        be counted using bedtools multicov

    bam1 : list or array
        a list of full paths to bam files pertaining to a single condition 
        (i.e. replicates of a single treatment)

    bam2 : list or array
        a list of full paths to bam files pertaining to a single condition 
        (i.e. replicates of a single treatment)

    tempdir : string
        full path to temp directory in output directory (created by TFEA)

    label1 : string
        the name of the treatment or condition corresponding to bam1 list

    label2 : string
        the name of the treatment or condition corresponding to bam2 list

    Returns
    -------
    None
    '''
    #If more than 1 replicate, use DE-Seq2
    if (len(bam1) > 1 and len(bam2) > 1):
        outfile = open(tempdir + 'DESeq.R','w')
        outfile.write('sink("'+tempdir+'DESeq.Rout")\n')
        outfile.write('library("DESeq2")\n')
        outfile.write('data <- read.delim("'+count_file+'", sep="\t", \
                        header=TRUE)\n')
        outfile.write('countsTable <- subset(data, select=c('
                +', '.join([str(i) for i in range(5,5+len(bam1)+len(bam2))])
                +'))\n')

        outfile.write('rownames(countsTable) <- data$region\n')
        outfile.write('conds <- as.data.frame(c(' 
                        + ', '.join(['"'+label1+'"']*len(bam1)) 
                        + ', ' 
                        + ', '.join(['"'+label2+'"']*len(bam2)) 
                        + '))\n')

        outfile.write('colnames(conds) <- c("treatment")\n')
        outfile.write('ddsFullCountTable <- DESeqDataSetFromMatrix(\
                                                    countData = countsTable, \
                                                    colData = conds, \
                                                    design = ~ treatment)\n')

        outfile.write('dds <- DESeq(ddsFullCountTable)\n')
        outfile.write('res1 <- results(dds,alpha = 0.05, \
                                        contrast=c("treatment",\
                                                        "'+label2+'",\
                                                        "'+label1+'"))\
                                                        \n')

        outfile.write('resShrink <- lfcShrink(dds, res = res1, \
                                                contrast = c("treatment",\
                                                "'+label2+'",\
                                                "'+label1+'"))\n')

        outfile.write('resShrink$fc <- 2^(resShrink$log2FoldChange)\n')
        outfile.write('res <- resShrink[c(1:3,7,4:6)]\n')
        outfile.write('write.table(res, file = "'
                        +tempdir+'DESeq.res.txt", \
                        append = FALSE, sep= "\t" )\n')
        outfile.write('sink()')
    else:
        outfile = open(tempdir + 'DESeq.R','w')
        outfile.write('sink("'+tempdir+'DESeq.Rout")\n')
        outfile.write('library("DESeq")\n')
        outfile.write('data <- read.delim("'+count_file+'", sep="\t", \
                        header=TRUE)\n')

        outfile.write('countsTable <- subset(data, select=c('
            +', '.join([str(i) for i in range(5,5+len(bam1)+len(bam2))])
            +'))\n')

        outfile.write('rownames(countsTable) <- data$region\n')
        outfile.write('conds <- c(' + ', '.join(['"'+label1+'"']*len(bam1)) 
                        + ', ' 
                        + ', '.join(['"'+label2+'"']*len(bam2)) 
                        + ')\n')

        outfile.write('cds <- newCountDataSet( countsTable, conds )\n')
        outfile.write('cds <- estimateSizeFactors( cds )\n')
        outfile.write('sizeFactors(cds)\n')                                                               
        outfile.write('cds <- estimateDispersions( cds ,method="blind", \
                        sharingMode="fit-only")\n')

        outfile.write('res <- nbinomTest( cds, "'+label1+'", "'+label2+'" )\n')
        outfile.write('rownames(res) <- res$id\n')                      
        outfile.write('write.table(res, file = "'+tempdir+'DESeq.res.txt", \
                        append = FALSE, sep= "\t" )\n')

        outfile.write('sink()')
    outfile.close()
#==============================================================================

#==============================================================================
def plot_deseq_MA(deseq_file=str(), label1=str(), label2=str(), 
                    figuredir=str()):
    '''Plots the DE-Seq MA-plot using the full regions of interest and saves it
    to the figuredir directory created in TFEA output folder

    Parameters
    ----------
    deseqfile : string
        full path to the deseq file (specifically .res.txt)

    label1 : string
        the name of the treatment or condition corresponding to bam1 list

    label2 : string
        the name of the treatment or condition corresponding to bam2 list

    figuredir : string
        full path to figure directory in output directory (created by TFEA)

    Returns
    -------
    None
    '''
    up_x = list()
    up_y = list()
    up_p = list()
    dn_x = list()
    dn_y = list()
    dn_p = list()
    with open(deseq_file,'r') as F:
        header = F.readline().strip('\n').split('\t')
        basemean_index = header.index('"baseMean"')
        log2fc_index = header.index('"log2FoldChange"')
        for line in F:
            line = line.strip('\n').split('\t')
            try:
                log2fc = float(line[log2fc_index+1])
                basemean = math.log(float(line[basemean_index+1]),10)
                pval = float(line[-2])
                if log2fc > 0:
                    up_x.append(basemean)
                    up_y.append(log2fc)
                    up_p.append(pval)
                else:
                    dn_x.append(basemean)
                    dn_y.append(log2fc)
                    dn_p.append(pval)
            except:
                pass

    x = [x for _,x in sorted(zip(up_p,up_x))] \
        + [x for _,x in sorted(zip(dn_p,dn_x),reverse=True)]

    y = [y for _,y in sorted(zip(up_p,up_y))] \
        + [y for _,y in sorted(zip(dn_p,dn_y),reverse=True)]

    c = plt.cm.RdYlGn(np.linspace(0, 1, len(x)))
    #Creates an MA-Plot of the region expression
    F = plt.figure(figsize=(7,6))
    ax = plt.subplot(111)
    plt.scatter(x=x,y=y,color=c,edgecolor='')
    ax.set_title("DE-Seq MA-Plot",fontsize=14)
    ax.set_ylabel("Log2 Fold-Change ("+label2+"/"+label1+")",fontsize=14)
    ax.set_xlabel("Log10 Average Expression",fontsize=14)
    ax.tick_params(axis='y', which='both', left='off', right='off', 
                    labelleft='on')
    ax.tick_params(axis='x', which='both', bottom='off', top='off', 
                    labelbottom='on')
    plt.savefig(figuredir + 'DESEQ_MA_Plot.png',bbox_inches='tight')
#==============================================================================

#==============================================================================
def permutations(distances=list(), permutations=1000):
    '''Generates permutations of the distances and calculates AUC for each 
        permutation.

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
#==============================================================================

#==============================================================================
def padj_bonferroni(TFresults=list()):
    '''This function iterates through TFEA results, removes TFs that returned 
        "no hits" and calculates a p-adj using the Bonferroni Correction for 
        each TF motif appending it to the given TFresults array

    Parameters
    ----------
    TFresults : list of lists
        contains calculated enrichment scores for all TFs of interest specified
        by the user
        
    Returns
    -------
    TFresults : list of lists
        same as input with an additional p-adjusted value appended to each TF
    '''
    TFresults = [x for x in TFresults if x != "no hits"]
    for i in range(len(TFresults)):
        PVAL = TFresults[i][-1]
        #Using Bonferroni Correction
        PADJ = 1 if PVAL*len(TFresults) > 1 else PVAL*len(TFresults)
        TFresults[i].append(PADJ)

    return TFresults
#==============================================================================

#==============================================================================
def motif_distance_bedtools_closest(ranked_center_file=str(), 
                                    motif_path=str()):
    '''Calculates nearest motif hit from a bed file. TFEA provides this 
        function with a bed file containing the center of the inputted regions.

    Parameters
    ----------
    TFresults : list of lists
        contains calculated enrichment scores for all TFs of interest specified
        by the user
        
    Returns
    -------
    motif_distance_bed_sorted : string
        full path to where the sorted motif distance file was outputted
    '''
    os.system("bedtools closest -D ref -t first -a " 
                + ranked_center_file.split('.bed')[0] + ".sorted.bed -b " 
                + motif_path + " > " + '/' 
                + '/'.join(ranked_center_file.split('/')[:-1]) + '/' 
                + motif_path.split('/')[-1] + ".sorted.distance.bed")

    motif_distance_bed_sorted = '/' \
                                + '/'.join(ranked_center_file.split('/')[:-1])\
                                + '/' + motif_path.split('/')[-1]\
                                + ".sorted.distance.bed"

    return motif_distance_bed_sorted
#==============================================================================

#==============================================================================
def fimo(tempdir=str(), motifdatabase=str(), bgfile=str(), motif=str(), 
            fastafile=str()):
    '''This function runs fimo on a given fastafile for a single motif in a 
        provided motif database. The output is cut and sorted to convert into 
        a sorted bed file

    Parameters
    ----------
    tempdir : string
        full path to temp directory in output directory (created by TFEA)

    motifdatabase : string
        full path to a motif database file in meme format

    bgfile : string
        full path to a markov background model

    motif : string
        the name of a motif that matches a motif within motifdatabase

    fastafile : string
        full path to a fasta file that fimo will perform motif scanning on
        
    Returns
    -------
    fimo_out : string
        full path to where fimo output which is stored within the tempdir 
        directory.
    '''
    os.system("fimo --text --bgfile "+bgfile+" --motif "+motif+" " +
              motifdatabase+" "+fastafile+" > "+tempdir+motif+".txt")
    os.system("cut -d$'\\t' -f 3- "+tempdir+motif +
              ".txt > " + tempdir+motif+".bed")

    fimo_out = tempdir+motif+".bed"

    return fimo_out
#==============================================================================

#==============================================================================
def meme2images(motifdatabase=str(), figuredir=str(), motif=str()):
    '''This function creates meme logos for use in the output html

    Parameters
    ----------
    motifdatabase : string
        full path to a motif database file in meme format

    motif : string
        the name of a motif that matches a motif within motifdatabase

    figuredir : string
        full path to figure directory in output directory (created by TFEA)
        
    Returns
    -------
    None
    '''
    os.system("meme2images -png -rc --motif " + motif + " " + motifdatabase
              + " " + figuredir)
#==============================================================================

#==============================================================================
def fasta_markov(tempdir=str(), fastafile=str(), order='0'):
    '''This function runs meme's fasta-get-markov function that generates a 
        background markov file (for use with fimo) from a fasta file.
    Parameters
    ----------
    tempdir : string
        full path to temp directory in output directory (created by TFEA)

    fastafile : string
        full path to fasta file that will be used to generate the markov
        background model file

    order : string
        an integer formatted as a string where a user may specify what order
        markov model they would like (default='0')
        
    Returns
    -------
    None
    '''
    os.system("fasta-get-markov -m " + order + " " + fastafile +
              " > " + tempdir + "Markov_Background.txt")
#==============================================================================

#==============================================================================
def fimo_parse(largewindow=float(), tempdir=str(), fimo_file=str(), 
                motif_file=str()):
    '''Parses a fimo output file and writes into a new file that is formatted
        in a way that can be parsed within existing TFEA functions
    Parameters
    ----------
    largewindow : float
        the size of the larger window to perform TFEA. Specified by user in
        config file

    tempdir : string
        full path to temp directory in output directory (created by TFEA)

    fimo_file : string
        full path to the fimo output file to be parsed by this function

    motif_file : string
        the name of the motif being parsed, this function will create a file
        using this motif_file string
        
    Returns
    -------
    outname : string
        the full path to the file to be used by other TFEA functions
    '''
    d = dict()
    with open(fimo_file) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            pval, fc, rank = line[0].split(',')
            start, stop = line[1:2]
            distance = ((int(start)+int(stop))/2)-largewindow
            score = line[5]
            if rank not in d:
                d[rank] = [rank, pval, fc, score, distance]
            else:
                prev_distance = math.fabs(d[rank][-1])
                if prev_distance > math.fabs(distance):
                    d[rank] = [rank, pval, fc, score, distance]

    outname = tempdir+motif_file+'.sorted.distance.bed'
    outfile = open(outname, 'w')
    for i in range(len(d)):
        rank = str(i)
        outfile.write('\t'.join(d[rank])+'\n')
    outfile.close()

    return outname
#==============================================================================

#==============================================================================
def convert_sequence_to_array(sequence=str()):
    '''This function converts a DNA sequence (ACGT alphabet) to an array, 
    collapsing GCs to 1's and ATs to 0's

    Parameters
    ----------
    sequence : string
        a string containing ACGT character
        
    Returns
    -------
    array : list
        a list of floats corresponding to 1.0 or 0.0 depending on whether the 
        input sequence was GC or AT at a specific site

    Raises
    ------
    Warning : str
        when a character is not ACTG. Value is given 0.0
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
#==============================================================================

#==============================================================================
def rank_deseqfile(deseq_file=str(), tempdir=str()):
    '''This function parses a DE-seq output file and creates a new file with 
        the center of each region ranked by p-value 
    
    Parameters
    ----------
    deseq_file : string
        full path to a DE-Seq output file
    
    tempdir : string
        full path to the tempdir directory in the output directory (created by 
        TFEA)
        
    Returns
    -------
    ranked_center_file : string
        full path to a bed file that contains the center of regions of interest
        ranked via DE-Seq p-value
    '''
    up = list()
    down = list()
    with open(deseq_file) as F:
        header = F.readline().strip('\n').split('\t')
        fc_index = [i for i in range(len(header)) if header[i]=='"fc"' or header[i]=='"foldChange"'][0]
        for line in F:
            line = line.strip('\n').split('\t')
            if line[fc_index] != 'NA':
                try:
                    pval = format(float(line[-2]),'.12f')
                except ValueError:
                    pval = format(1.0,'.12f')
                region = line[0].split(':')
                chrom = region[0]
                coordinates = region[1].split('-')
                start = coordinates[0]
                stop = coordinates[1]
                chrom = chrom.strip('"')
                stop = stop.strip('"')
                fc = float(line[fc_index+1])
                if fc < 1:
                    down.append((chrom,start,stop,pval,str(fc)))
                else:
                    up.append((chrom,start,stop,pval,str(fc)))

    #Save ranked regions in a bed file (pvalue included)
    outfile = open(tempdir + "ranked_file.bed",'w')
    r=1
    for region in sorted(up, key=lambda x: x[3]):
        outfile.write('\t'.join(region) + '\t' + str(r) + '\n')
        r += 1
    for region in sorted(down, key=lambda x: x[3], reverse=True):
        outfile.write('\t'.join(region) + '\t' + str(r) + '\n')
        r += 1
    outfile.close()

    #Get center base for each region
    outfile = open(tempdir+"ranked_file.center.bed",'w')
    with open(tempdir + "ranked_file.bed") as F:
        for line in F:
            line = line.strip('\n').split('\t')
            chrom,start,stop = line[:3]
            center = (int(start)+int(stop))/2
            outfile.write(chrom + '\t' + str(center) + '\t' + str(center+1) 
                            + '\t' + '\t'.join(line[3:]) + '\n')
    outfile.close()


    os.system("sort -k1,1 -k2,2n " + tempdir+"ranked_file.center.bed" + " > " 
                + tempdir + "ranked_file.center.sorted.bed")

    ranked_center_file = tempdir + "ranked_file.center.sorted.bed"

    return ranked_center_file
#==============================================================================

#==============================================================================
def samtools_flagstat(args):
    '''Performs samtools flagstat on a bam file. Then parses the samtools 
        flagstat output and returns a the millions mapped reads for the given
        bam file.

    Parameters
    ----------
    args : tuple
        contains arguments for this function:
        bam : string
            full paths to a bam files
        tempdir : string
            full path to the tempdir directory within the output directory 
            (created by TFEA)
        
    Returns
    -------
    millions_mapped : float
        millions mapped reads for the given bam file
    '''
    bam, tempdir = args
    filename = bam.split('/')[-1]
    os.system("samtools flagstat " + bam + " > " + tempdir + filename 
                + ".flagstat")
    with open(tempdir+filename+".flagstat") as F:
        lines = F.readlines()
        millions_mapped = float(lines[4].strip('\n').split(' ')[0])/1000000.0

        return millions_mapped
#==============================================================================

#==============================================================================
def enrichment_plot(largewindow=float(),
                    smallwindow=float(), figuredir=str(),
                    cumscore=list(), sorted_distances=list(), logpval=list(), 
                    updistancehist=list(), downdistancehist=list(), 
                    gc_array=list(), motif_file=''):
    '''This function plots the TFEA enrichment plot.

    Parameters
    ----------
    largewindow : float
        a user specified larger window used for plotting purposes and to do
        some calculations regarding the user-provided regions of interest

    smallwindow : float
        a user specified smaller window used for plotting purposes and to do
        some calculations regarding the user-provided regions of interest

    figuredir : string
        the full path to the figuredir within the ouptut directory containing
        all figures and plots

    cumscore : list
        the cumulative score of the running sum as we walk through the ranked
        regions

    sorted_distances : list or array
        a sorted (based on rank) list of motif distances. Negative corresponds
        to upstream of region

    logpval : list or array
        a way to visualize the ranking of the regions of interest. It is
        the log10 of the p-value with the sign (positive or negative) based on
        whether the fold change of the region is over 1 or less than 1.

    updistancehist : list or array
        the first quartile of ranked regions. These are presumably higher in
        condition1

    downdistancehist : list or array
        the fourth quartile of ranked regions. These are presumably higher in
        condition2

    gc_array : list
        an array of gc richness of regions of interest. It is recommended that
        this array be no larger than 1000 bins.

    motif_file : string
        the name of the motif thats associated with all the input data. Used 
        for figure naming purposes.

    Returns
    -------
    None
    '''
    #Begin plotting section
    plt.figure(figsize=(15.5,8))
    xvals = range(1,len(cumscore)+1)
    limits = [1,len(cumscore)]
    gs = gridspec.GridSpec(4, 1, height_ratios=[2, 2, 1, 1])

    #This is the enrichment score plot (i.e. line plot)
    ax0 = plt.subplot(gs[0])
    ax0.plot(xvals,cumscore,color='green')
    ax0.plot([0, len(cumscore)],[0, 1], '--', alpha=0.75)
    ax0.set_title('Enrichment Plot: ',fontsize=14)
    ax0.set_ylabel('Enrichment Score (ES)', fontsize=10)
    ax0.tick_params(axis='y', which='both', left='on', right='off', 
                    labelleft='on')
    ax0.tick_params(axis='x', which='both', bottom='off', top='off', 
                    labelbottom='off')
    ylims = ax0.get_ylim()
    ymax = math.fabs(max(ylims,key=abs))
    ax0.set_ylim([0,ymax])
    ax0.set_xlim(limits)

    #This is the distance scatter plot right below the enrichment score 
    #plot
    ax1 = plt.subplot(gs[1])
    ax1.scatter(xvals,sorted_distances,edgecolor="",color="black",s=10,
                alpha=0.25)
    ax1.tick_params(axis='y', which='both', left='off', right='off', 
                    labelleft='on')
    ax1.tick_params(axis='x', which='both', bottom='off', top='off', 
                    labelbottom='off')
    ax1.set_xlim(limits)
    ax1.set_ylim([-int(largewindow),int(largewindow)])
    plt.yticks([-int(largewindow),0,int(largewindow)],
                [str(-int(largewindow)/1000.0),'0',\
                str(int(largewindow)/1000.0)])
    ax1.set_ylabel('Distance (kb)', fontsize=10)

    #This is the rank metric plot
    ax2 = plt.subplot(gs[3])
    ax2.fill_between(xvals,0,logpval,facecolor='grey',edgecolor="")
    ax2.tick_params(axis='y', which='both', left='on', right='off', 
                    labelleft='on')
    ax2.tick_params(axis='x', which='both', bottom='off', top='off', 
                    labelbottom='on')
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
        ax2.axvline(len(xvals) - len(downdistancehist), color='purple', 
                    alpha=0.25)
    except ValueError:
        pass

    #This is the GC content plot
    ax3 = plt.subplot(gs[2])
    ax3.set_xlim(limits)
    sns.heatmap(gc_array, cbar=False, xticklabels='auto',
                yticklabels='auto')

    plt.yticks([0,int(largewindow),int(largewindow*2)],
                [str(-int(largewindow)/1000.0),'0',\
                str(int(largewindow)/1000.0)])

    ax3.tick_params(axis='y', which='both', left='on', right='off', 
                    labelleft='on')

    ax3.tick_params(axis='x', which='both', bottom='off', top='off', 
                    labelbottom='off')

    ax3.set_ylabel('GC content per kb',fontsize=10)

    plt.savefig(figuredir + motif_file.split('.bed')[0] 
                + '_enrichment_plot.png',bbox_inches='tight')

    plt.cla()
#==============================================================================

#==============================================================================
def simulation_plot(figuredir=str(), simES=list(), actualES=float(),
                        motif_file=str()):
    '''This function plots the simulated 'enrichment' scores against the
        observed 'enrichment' score

    Parameters
    ----------
    figuredir : string
        the full path to the figuredir within the ouptut directory containing
        all figures and plots

    simES : list
        a list of simulated 'enrichment' scores calculated by randomizing
        region rank

    actualES : float
        the observed 'enchrichment' score. Can be calculated in many different
        ways..

    motif_file : string
        the name of the motif thats associated with all the input data. Used 
        for figure naming purposes.

    Returns
    -------
    None
    '''
    F = plt.figure(figsize=(7,6))
    ax2 = plt.subplot(111)
    maximum = max(simES)
    minimum = min(simES)
    ax2.hist(simES,bins=100)
    width = (maximum-minimum)/100.0
    rect = ax2.bar(actualES,ax2.get_ylim()[1],color='red',width=width*2)[0]
    height = rect.get_height()
    ax2.text(rect.get_x() + rect.get_width()/2., 1.05*height, 
                'Observed ES', ha='center', va='bottom')

    ax2.set_xlim([min(minimum,actualES)-(width*40), \
                max(maximum,actualES)+(width*40)])

    ax2.set_ylim([0,(1.05*height)+5])
    ax2.tick_params(axis='y', which='both', left='off', right='off', 
                    labelleft='on')

    ax2.tick_params(axis='x', which='both', bottom='off', top='off', 
                    labelbottom='on')

    plt.title('Distribution of Simulated Enrichment Scores',fontsize=14)
    ax2.set_ylabel('Number of Simulations',fontsize=14)
    ax2.set_xlabel('Enrichment Score (ES)',fontsize=14)
    plt.savefig(figuredir + motif_file.split('.bed')[0] 
                + '_simulation_plot.png',bbox_inches='tight')

    plt.cla()
#==============================================================================

#==============================================================================
def distance_distribution_plot(largewindow=float(),
                                smallwindow=float(), 
                                figuredir=str(), 
                                updistancehist=list(), 
                                middledistancehist=list(),
                                downdistancehist=list(), motif_file=str()):
    '''This function plots histograms of motif distances to region centers 
        based on ranking quartiles

    Parameters
    ----------
    largewindow : float
        a user specified larger window used for plotting purposes and to do
        some calculations regarding the user-provided regions of interest

    smallwindow : float
        a user specified smaller window used for plotting purposes and to do
        some calculations regarding the user-provided regions of interest

    figuredir : string
        the full path to the figuredir within the ouptut directory containing
        all figures and plots

    updistancehist : list or array
        the first quartile of ranked regions. These are presumably higher in
        condition1

    middledistancehist : list or array
        the second and third quartile of ranked regions. These are presumably 
        unchanging regions

    downdistancehist : list or array
        the fourth quartile of ranked regions. These are presumably higher in
        condition2

    motif_file : string
        the name of the motif thats associated with all the input data. Used 
        for figure naming purposes.
    '''
    plt.figure(figsize=(6.5,6))
    gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 1])
    ax0 = plt.subplot(gs[0])
    binwidth = largewindow/100.0
    ax0.hist(updistancehist,
                bins=np.arange(0,int(largewindow)+binwidth,binwidth),
                color='green')
    ax0.set_title('Distribution of Motif Distance for: fc > 1',fontsize=14)
    ax0.axvline(smallwindow,color='red',alpha=0.5)
    ax0.tick_params(axis='y', which='both', left='off', right='off', 
                    labelleft='on')
    ax0.tick_params(axis='x', which='both', bottom='off', top='off', 
                    labelbottom='off')
    ax0.set_xlim([0,largewindow])
    ax0.set_ylabel('Hits',fontsize=14)
    ax1 = plt.subplot(gs[2])
    ax1.hist(downdistancehist,
                bins=np.arange(0,int(largewindow)+binwidth,binwidth),
                color='purple')

    ax1.axvline(smallwindow,color='red',alpha=0.5)
    ax1.set_title('Distribution of Motif Distance for: fc < 1',fontsize=14)
    ax1.tick_params(axis='y', which='both', left='off', right='off', 
                    labelleft='on')

    ax1.tick_params(axis='x', which='both', bottom='off', top='off', 
                    labelbottom='on')

    ax1.set_xlim([0,largewindow])
    ax1.set_ylabel('Hits',fontsize=14)
    ax1.set_xlabel('Distance (bp)',fontsize=14)
    ax2 = plt.subplot(gs[1])
    ax2.hist(middledistancehist,
                bins=np.arange(0,int(largewindow)+binwidth,binwidth),
                color='blue')

    ax2.set_title('Distribution of Motif Distance for: middle',fontsize=14)
    ax2.axvline(smallwindow,color='red',alpha=0.5)
    ax2.tick_params(axis='y', which='both', left='off', right='off', 
                    labelleft='on')
    ax2.tick_params(axis='x', which='both', bottom='off', top='off', 
                    labelbottom='off')
    ax2.set_xlim([0,largewindow])
    ax2.set_ylabel('Hits',fontsize=14)
    plt.savefig(figuredir + motif_file.split('.bed')[0] 
                + '_distance_distribution.png',bbox_inches='tight')
    
    plt.cla()
#==============================================================================

#==============================================================================
def moustache_plot(figuredir=str(),ESlist=list(),PADJlist=list(),
                    sigx=list(),sigy=list()):

    '''This function plots a moustache plot for all motifs. In the x-axis, 
        all observed 'enrichment' scores are plotted against the adjusted
        pvalue of each motif

    Parameters
    ----------
    figuredir : string
        the full path to the figuredir within the ouptut directory containing
        all figures and plots

    ESlist : list or array
        a list of 'enrichment' scores to be plotted on the x-axis

    PADJlist : list or array
        a list of p-adjusted values to be plotted on the y-axis

    sigx : list or array
        a list of significant motifs to be colored red

    sigy : list or array
        a list of significant motifs to be colored red

    Returns
    -------
    None
    '''
    plt.figure(figsize=(7,6))
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
    ax.tick_params(axis='y', which='both', left='off', right='off', 
                    labelleft='on')

    ax.tick_params(axis='x', which='both', bottom='off', top='off', 
                    labelbottom='on')

    plt.savefig(figuredir + 'TFEA_Results_Moustache_Plot.png',
                    bbox_inches='tight')
    plt.cla()
#==============================================================================

#==============================================================================
def pval_histogram_plot(figuredir=str(), PVALlist=list()):
    '''This function plots a histogram of p-values for all motifs

    Parameters
    ----------
    figuredir : string
        the full path to the figuredir within the ouptut directory containing
        all figures and plots

    PVALlist : list or array
        a list of p-values corresponding to the observed 'enrichment' score
        compared to the distribution of simulated 'enrichment' scores.
    
    Returns
    -------
    None
    '''
    plt.figure(figsize=(7,6))
    ax = plt.subplot(111)
    binwidth = 1.0/100.0
    ax.hist(PVALlist,bins=np.arange(0,0.5+binwidth,binwidth),color='green')
    ax.set_title("TFEA P-value Histogram",fontsize=14)
    ax.set_xlabel("P-value",fontsize=14)
    ax.set_ylabel("Count",fontsize=14)
    ax.tick_params(axis='y', which='both', left='off', right='off', 
                    labelleft='on')

    ax.tick_params(axis='x', which='both', bottom='off', top='off', 
                    labelbottom='on')

    plt.savefig(figuredir + 'TFEA_Pval_Histogram.png',
                    bbox_inches='tight')
    plt.cla()
#==============================================================================

#==============================================================================
def MA_plot(figuredir=str(), label1=str(), label2=str(), POSlist=list(), 
                ESlist=list(), MAx=list(), MAy=list()):
    '''This function plots an 'MA' plot with the 'enrichment' score on the 
        y-axis and the number of hits within the largewindow in the x-axis

    Parameters
    ----------
    figuredir : string
        the full path to the figuredir within the ouptut directory containing
        all figures and plots

    label1 : string
        a label that corresponds to condition1

    label2 : string
        a label that corresponds to condition2

    POSlist : list or array
        a list of 'positive' hits for each motif defined as being within a 
        largewindow

    ESlist : list or array
        a list of 'enrichment' scores for each motif

    MAx : list or array
        a list of x-values corresponding to significant motifs to be colored 
        red

    MAx : list or array
        a list of y-values corresponding to significant motifs to be colored 
        red

    Returns
    -------
    None
    '''
    plt.figure(figsize=(7,6))
    ax = plt.subplot(111)
    ax.scatter(POSlist,ESlist,color='black',edgecolor='')
    ax.scatter(MAx,MAy,color='red',edgecolor='')
    ax.set_title("TFEA MA-Plot",fontsize=14)
    ax.set_ylabel("Normalized Enrichment Score (NES) " + label2 + "/" 
                    + label1, fontsize=14)

    ax.set_xlabel("Hits Log10",fontsize=14)
    ax.tick_params(axis='y', which='both', left='off', right='off', 
                    labelleft='on')

    ax.tick_params(axis='x', which='both', bottom='off', top='off', 
                    labelbottom='on')

    plt.savefig(figuredir + 'TFEA_NES_MA_Plot.png',bbox_inches='tight')
    plt.cla()
#==============================================================================

#==============================================================================
def create_text_output(outputdir=str(), TFresults=list()):
    '''Creates a .txt output of results from TFEA

    Parameters
    ----------
    outputdir : string
        the full path to TFEA output directory (created by TFEA)

    TFresults : list or array
        a list of lists contining 'enrichment' scores, normalized 'enrichment' 
        scores, p-value, p-adj, and number of hits for each individual motif

    Returns
    -------
    None
    '''
    TFresults = sorted(TFresults, key=lambda x: x[3])
    outfile = open(outputdir + 'results.txt', 'w')
    outfile.write('TF-Motif\tES\tNES\tP-value\tPADJ\tHITS\n')
    for val in TFresults:
        outfile.write('\t'.join([str(val[i]) for i in range(len(val))]) + '\n')
    outfile.close()
#==============================================================================

#==============================================================================
def create_html_output(TFresults=list(), config_dict=dict(), 
                        COMBINEtime=float(), COUNTtime=float(), 
                        DESEQtime=float(), CALCULATEtime=float()):
    '''Creates the main html output and also individual html outputs for each
        motif
    
    Parameters
    ----------
    outputdir : string
        the full path to the output directory created by TFEA

    beds : list or array
        a list of full paths to bed files to be considered as regions of 
        interest

    label1 : string
        an informative label describing sample corresponding to condition1

    label2 : string
        an informative label describing sample corresponding to condition2

    bam1 : list or array
        a list of full paths to bam files corresponding to condition1

    bam2 : list or array
        a list of full paths to bam files corresponding to condition2

    singlemotif : boolean or string
        either False if all motifs should be considered in TFEA or the name of
        a specific motif to be analyzed

    motif_hits : string
        the full path to a directory containing motif hits across the genome

    output : string
        the full path to a user-specified output directory. TFEA will create
        a new folder within this directory - this is called outputdir

    padj_cutoff : float
        the cutoff value for determining significance

    plot : boolean
        a switch that controls whether all motifs are plotted or just 
        significant ones defined by the p-adj cutoff

    combine : boolean
        a switch that determines whether bed files within the beds variable
        get combined and merged using bedtools

    count : boolean
        a switch that controls whether reads are counted over the regions of
        interest

    deseq : boolean
        a switch that controls whether DE-Seq is performed on the inputted
        regions that have been counted over

    calculate : boolean
        a switch that determines whether the TFEA calculation is performed

    TFresults : list or array
        a list of lists contining 'enrichment' scores, normalized 'enrichment' 
        scores, p-value, p-adj, and number of hits for each individual motif

    COMBINEtime : float
        the time it took to combine and merge the bed files using bedtools

    COUNTtime : float
        the time it took to count reads over regions of interest

    DESEQtime : float
        the time it took to perform DE-Seq using the counts file

    CALCULATEtime : float
        the time it took to perform TFEA

    Returns
    -------
    None
    '''
    outputdir = config_dict['OUTPUTDIR'] 
    beds = config_dict['BEDS']
    label1 = config_dict['LABEL1']
    label2 = config_dict['LABEL2']
    bam1 = config_dict['BAM1']
    bam2 = config_dict['BAM2']
    singlemotif = config_dict['SINGLEMOTIF']
    motif_hits = config_dict['MOTIF_HITS']
    output = config_dict['OUTPUT']
    padj_cutoff = config_dict['PADJCUTOFF']
    plot = config_dict['PLOT']
    combine = config_dict['COMBINE']
    count = config_dict['COUNT']
    deseq = config_dict['DESEQ']
    calculate = config_dict['CALCULATE']
    
    #Creates results.txt which is a tab-delimited text file with the results    
    TFresults = sorted(TFresults, key=lambda x: x[5])

    #summary.html contains all user-defined variables, and also information 
    #about module used
    outfile = open(outputdir+'summary.html','w')
    outfile.write("""<!DOCTYPE html>
            <html>
            <head>
            <title>Variables Used</title>
            </head>
            <body>
                <h1>Variables Used</h1>
                <p>BEDS = """+str(beds)+"""</p>
                <p>LABEL1 = """+label1+"""</p>
                <p>LABEL2 = """+label2+"""</p>
                <p>BAM1 = """+str(bam1)+"""</p>
                <p>BAM2 = """+str(bam2)+"""</p>
                <p>SINGLEMOTIF = """+str(singlemotif)+"""</p>
                <p>MOTIF_HITS = """+str(motif_hits)+"""</p>
                <p>OUTPUT = """+output+"""
            </body>""")

    #For each TF motif with an PADJ value less than a cutoff, an html file is 
    #created to be used in results.html
    positivelist = [x[0] for x in TFresults if x[2] > 0 and (plot or x[-1] < padj_cutoff)]
    negativelist = [x[0] for x in TFresults if x[2] < 0 and (plot or x[-1] < padj_cutoff)]
    for i in range(len(TFresults)):
        MOTIF_FILE,ES,NES,PVAL,POS,PADJ = TFresults[i] 
        
        if NES > 0:
            try:
                NEXT_MOTIF = positivelist[positivelist.index(MOTIF_FILE)+1]
            except IndexError:
                NEXT_MOTIF = positivelist[0]
            try:
                PREV_MOTIF = positivelist[positivelist.index(MOTIF_FILE)-1]
            except IndexError:
                PREV_MOTIF = positivelist[len(positivelist)]
        else:
            try:
                NEXT_MOTIF = negativelist[negativelist.index(MOTIF_FILE)+1]
            except IndexError:
                NEXT_MOTIF = negativelist[0]
            try:
                PREV_MOTIF = negativelist[negativelist.index(MOTIF_FILE)-1]
            except IndexError:
                PREV_MOTIF = negativelist[len(negativelist)]
        if plot or PADJ < padj_cutoff:
            outfile = open(outputdir + 'plots/' + MOTIF_FILE 
                            + '.results.html','w')
            outfile.write("""<!DOCTYPE html>
    <html>
    <head>
    <title>"""+MOTIF_FILE+""" Results</title>
    <style>
        table {
            font-family: arial, sans-serif;
            border-collapse: collapse;
            width: 100%;
        }

        .row {
        display: flex; /* equal height of the children */
        width: 100%;
        padding-bottom: 50px
        }

        img {
            max-width: 100%;
            max-height: 100%;
        }

        td, th {
            border: 1px solid #dddddd;
            text-align: left;
            padding: 8px;
        }

        tr:nth-child(even) {
            background-color: #dddddd;
        }
    </style>
    </head>
    <body style="width:1300px; margin:0 auto;">
        <div>
            <div style="float:left">
                <a href="./"""+PREV_MOTIF+""".results.html">PREV</a>
            </div>
            <div style="float:right">
                <a href="./"""+NEXT_MOTIF+""".results.html">NEXT</a>
            </div>
            <div style="text-align:center">
                <a href="../results.html">ALL</a>
        </div>
        <div class="row">
        </div>
            <h1>"""+MOTIF_FILE+""" Results</h1>
        <div>
            <div style="float: middle; width: 1300px; padding-bottom:25px; \
                padding-top:25px">
                <table> 
                    <tr>
                        <th>TF Motif</th>
                        <th>ES</th> 
                        <th>NES</th>
                        <th>P-value</th>
                        <th>PADJ</th>
                        <th>HITS</th>
                    </tr>
                    <tr>
                        <td>"""+MOTIF_FILE+"""</td>
                        <td>"""+str("%.2f" % ES)+"""</td>
                        <td>"""+str("%.2f" % NES)+"""</td>
                        <td>"""+str("%.4g" % PVAL)+"""</td>
                        <td>"""+str("%.4g" % PADJ)+"""</td>
                        <td>"""+str(POS)+"""</td>
                    </tr>
                </table>
            </div>
        </div>
        <div>
            <div style="float: left; width 1250px; padding-bottom:50px; \
                padding-top:50px">
                <img src="./"""+MOTIF_FILE+"""_enrichment_plot.png" \
                    alt="Enrichment Plot">
            </div>
        </div>
        <!--<div class="row">
            <div style="float: left; width 1250px; padding-bottom:50px; \
                padding-top:50px">
                <img src="./"""+MOTIF_FILE+"""_meta_eRNA.png" alt="Meta Plot">
            </div>
        </div>-->
        <div class="row">
            <div style="float: right; width: 600px">
                <p>Forward:</p>
                <img src="./"""+MOTIF_FILE.split('HO_')[-1]+"""_direct.png" \
                    alt="Forward Logo">
                <p></p>
                <p>Reverse:</p>
                <img src="./"""+MOTIF_FILE.split('HO_')[-1]+"""_revcomp.png" \
                    alt="Reverse Logo">
            </div>
            <div style="float:left; width: 600px">
                <img src="./"""+MOTIF_FILE+"""_simulation_plot.png" \
                    alt="Simulation Plot">
            </div>
        </div>

    </body>
    </html>""")
            outfile.close()
            PREV_MOTIF = MOTIF_FILE

    outfile = open(outputdir+'results.html','w')
    outfile.write("""<!DOCTYPE html>
    <html>
    <head>
    <title>TFEA Results</title>
    <style>
        table {
            font-family: arial, sans-serif;
            border-collapse: collapse;
            width: 100%;
        }

        .row {
        display: flex; /* equal height of the children */
        width: 100%;
        padding-bottom: 50px
        }

        img {
            max-width: 100%;
            max-height: 100%;
        }

        td, th {
            border: 1px solid #dddddd;
            text-align: left;
            padding: 8px;
        }

        tr:nth-child(even) {
            background-color: #dddddd;
        }
    </style>
    </head>
    <body style="width:1300px; margin:0 auto;">

        <h1>TFEA Results """ +label1+ """ vs. """ +label2
            +"""</h1>
        <div class="row">
            <div style="float: left; width: 45%">
                <img src="./plots/TFEA_NES_MA_Plot.png" alt="NES MA-Plot">
            </div>
            <div style="float: left; width: 45%">
                <img src="./plots/DESEQ_MA_Plot.png" alt="DE-Seq MA-Plot">
            </div>
        </div>
        <div class="row">
            <div style="float: left; width:45%">
                <img src="./plots/TFEA_Pval_Histogram.png" alt="TFEA P-Value \
                    Histogram">
            </div>
            <div id="Summary of Variables Used" style="float: right; \
                width: 45%">
                <p><a href="./Summary.html">Full Summary of Variables Used\
                </a></p>
                <p><b>PADJ < """ + str(padj_cutoff) + """</b></p>
                <table>
                    <tr>
                        <th>Module</th>
                        <th>Switch</th>
                        <th>Time (hh:mm:ss)</th>
                    </tr>
                    <tr>
                        <td>COMBINE</td>
                        <td>"""+str(combine)+"""</td>
                        <td>"""+str(datetime.timedelta(
                                    seconds=int(COMBINEtime)))
                        +"""</td>
                    </tr>
                    <tr>
                        <td>COUNT</td>
                        <td>"""+str(count)+"""</td>
                        <td>"""+str(datetime.timedelta(seconds=int(COUNTtime)))
                        +"""</td>
                    </tr>
                    <tr>
                        <td>DESEQ</td>
                        <td>"""+str(deseq)+"""</td>
                        <td>"""+str(datetime.timedelta(seconds=int(DESEQtime)))
                        +"""</td>
                    </tr>
                    <tr>
                        <td>CALCULATE</td>
                        <td>"""+str(calculate)+"""</td>
                        <td>"""+str(datetime.timedelta(
                                    seconds=int(CALCULATEtime)))
                        +"""</td>
                    </tr>
                    <tr>
                        <td><b>TOTAL</b></td>
                        <td> </td>
                        <td>"""+str(datetime.timedelta(
                            seconds=int(COMBINEtime)
                                    +int(COUNTtime)
                                    +int(DESEQtime)
                                    +int(CALCULATEtime)))
                            +"""</td>
                    </tr>
                </table>   
            </div>
        </div>
        <div>
            <div id="Positive Enrichment Score" style="float: left; width:45%">
                <h1>Positive Enrichment Score</h1>
                <table> 
                    <tr>
                        <th>TF Motif</th>
                        <th>ES</th> 
                        <th>NES</th>
                        <th>P-value</th>
                        <th>PADJ</th>
                        <th>HITS</th>
                    </tr>
                """)

    for MOTIF_FILE,ES,NES,PVAL,POS,PADJ in TFresults:
        if NES > 0:
            if PADJ < padj_cutoff:
                outfile.write("""
            <tr style="color: red;">
                <td><a href="./plots/"""+MOTIF_FILE+""".results.html">"""
                    +MOTIF_FILE+"""</td>
                <td>"""+str("%.2f" % ES)+"""</td>
                <td>"""+str("%.2f" % NES)+"""</td>
                <td>"""+str("%.3g" % PVAL)+"""</td>
                <td>"""+str("%.3g" % PADJ)+"""</td>
                <td>"""+str(POS)+"""</td>
            </tr>
                    """)
            elif plot:
                outfile.write("""
            <tr>
                <td><a href="./plots/"""+MOTIF_FILE+""".results.html">"""
                    +MOTIF_FILE+"""</td>
                <td>"""+str("%.2f" % ES)+"""</td>
                <td>"""+str("%.2f" % NES)+"""</td>
                <td>"""+str("%.3g" % PVAL)+"""</td>
                <td>"""+str("%.3g" % PADJ)+"""</td>
                <td>"""+str(POS)+"""</td>
            </tr>
                    """)

            else:
                outfile.write("""
            <tr>
                <td>"""+MOTIF_FILE+"""</td>
                <td>"""+str("%.2f" % ES)+"""</td>
                <td>"""+str("%.2f" % NES)+"""</td>
                <td>"""+str("%.3g" % PVAL)+"""</td>
                <td>"""+str("%.3g" % PADJ)+"""</td>
                <td>"""+str(POS)+"""</td>
            </tr>
                    """)


    outfile.write("""            
        </table>
    </div>

    <div id="Negative Enrichment Score" style="float: right; width: 45%">
        <h1>Negative Enrichment Score</h1>
        <table> 
            <tr>
                <th>TF Motif</th>
                <th>ES</th> 
                <th>NES</th>
                <th>P-value</th>
                <th>PADJ</th>
                <th>HITS</th>
            </tr>
                """)

    for MOTIF_FILE,ES,NES,PVAL,POS,PADJ in TFresults:
        if NES < 0:
            if PADJ < padj_cutoff:
                outfile.write("""
            <tr style="color: red;">
                <td><a href="./plots/"""+MOTIF_FILE+""".results.html">"""
                    +MOTIF_FILE+"""</td>
                <td>"""+str("%.2f" % ES)+"""</td>
                <td>"""+str("%.2f" % NES)+"""</td>
                <td>"""+str("%.3g" % PVAL)+"""</td>
                <td>"""+str("%.3g" % PADJ)+"""</td>
                <td>"""+str(POS)+"""</td>
            </tr>
                    """)
            elif plot:
                outfile.write("""
            <tr>
                <td><a href="./plots/"""+MOTIF_FILE+""".results.html">"""
                    +MOTIF_FILE+"""</td>
                <td>"""+str("%.2f" % ES)+"""</td>
                <td>"""+str("%.2f" % NES)+"""</td>
                <td>"""+str("%.3g" % PVAL)+"""</td>
                <td>"""+str("%.3g" % PADJ)+"""</td>
                <td>"""+str(POS)+"""</td>
            </tr>
                    """)
            else:
                outfile.write("""
            <tr>
                <td>"""+MOTIF_FILE+"""</td>
                <td>"""+str("%.2f" % ES)+"""</td>
                <td>"""+str("%.2f" % NES)+"""</td>
                <td>"""+str("%.3g" % PVAL)+"""</td>
                <td>"""+str("%.3g" % PADJ)+"""</td>
                <td>"""+str(POS)+"""</td>
            </tr>
                    """)

    outfile.write("""        
            </table>
        </div>
        </div>

    </body>
    </html>""")

    outfile.close()
#==============================================================================

#==============================================================================
def get_TFresults_from_txt(results_file=str()):
    '''This function parses a results_file into a TFresults array

    Parameters
    ----------
    results_file : string
        file containing TFEA results
        
    Returns
    -------
    TFresults : array
        results for all motifs within the results_file
                 
    Raises
    ------
    None
    '''
    TFresults = list()
    with open(results_file) as F:
        F.readline()
        for line in F:
            line = line.strip('\n').split('\t')
            TFresults.append(line)
    return TFresults
#==============================================================================

#==============================================================================
def create_single_motif_html(outputdir=str(), results=list()):
    MOTIF_FILE,ES,NES,PVAL,POS = results
    outfile = open(os.path.join(outputdir, MOTIF_FILE + '.results.html'),'w')
    outfile.write("""<!DOCTYPE html>
        <html>
        <head>
        <title>"""+MOTIF_FILE+""" Results</title>
        <style>
        table {
            font-family: arial, sans-serif;
            border-collapse: collapse;
            width: 100%;
        }

        td, th {
            border: 1px solid #dddddd;
            text-align: left;
            padding: 8px;
        }

        tr:nth-child(even) {
            background-color: #dddddd;
        }
        </style>
        </head>
        <body style="width:1300px; margin:0 auto;">
            <h1>"""+MOTIF_FILE+""" Results</h1>
        <div>
            <div style="float: middle; width: 1300px; overflow:scroll; \
                padding-bottom:25px; padding-top:25px">
                <table> 
                    <tr>
                        <th>TF Motif</th>
                        <th>ES</th> 
                        <th>NES</th>
                        <th>P-value</th>
                        <th>Total Hits</th>
                       
                    </tr>
                    <tr>
                        <td>"""+MOTIF_FILE+"""</td>
                        <td>"""+str("%.2f" % ES)+"""</td>
                        <td>"""+str("%.2f" % NES)+"""</td>
                        <td>"""+str("%.4g" % PVAL)+"""</td>
                        <td>"""+str(POS)+"""</td>
                       
                    </tr>
                </table>
            </div>
        </div>
        <div>
            <div style="float: left; width 1250px; padding-bottom:50px; \
                padding-top:50px">
                <img src="./plots/"""+MOTIF_FILE+"""_enrichment_plot.png" \
                    alt="Enrichment Plot">
            </div>
        </div>
        <!--<div>
            <div style="float: left; width 1250px; padding-bottom:50px; \
                padding-top:50px">
                <img src="./plots/"""+MOTIF_FILE+"""_meta_eRNA.png" \
                    alt="Meta Plot">
            </div>
        </div>-->
        <div>
            <div style="float: right; width: 600px; overflow: scroll">
                <p>Forward:</p>
                <img src="./plots/"""+MOTIF_FILE+"""_direct.png" \
                    alt="Forward Logo">
                <p></p>
                <p>Reverse:</p>
                <img src="./plots/"""+MOTIF_FILE+"""_revcomp.png" \
                    alt="Reverse Logo">
            </div>
            <div style="float:left; width: 600px overflow:scroll">
                <img src="./plots/"""+MOTIF_FILE+"""_simulation_plot.png" \
                    alt="Simulation Plot">
            </div>
        </div>
        
        </body>
        </html>""")
    outfile.close()
#==============================================================================

#==============================================================================