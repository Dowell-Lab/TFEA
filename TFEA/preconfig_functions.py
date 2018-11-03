#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''This file contains a list of independent functions that do not call other 
    functions
'''
#==============================================================================
__author__ = 'Jonathan D. Rubin and Rutendo F. Sigauke'
__credits__ = ['Jonathan D. Rubin', 'Rutendo F. Sigauke', 'Jacob T. Stanley',
                'Robin D. Dowell']
__maintainer__ = 'Jonathan D. Rubin'
__email__ = 'Jonathan.Rubin@colorado.edu'
#==============================================================================
import os
import shutil
#==============================================================================
#Functions used before 'import config' statement
#==============================================================================
def make_out_directories(dirs=False, config_object=None):
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
    output = config_object['DATA']['OUTPUT'].strip("'")
    name = config_object['DATA']['NAME'].strip("'")
    outfoldername = name
    output = os.path.join(output, outfoldername)
    if dirs:
        if not os.path.isdir(output):
            os.makedirs(output)
        else:
            shutil.rmtree(output)
            os.makedirs(output)

    #Temporary files will go in this directory
    tempdir = os.path.join(output, 'temp_files')
    if dirs:
        if not os.path.isdir(tempdir):
            os.makedirs(tempdir)

    #Error and out files will go in this directory
    e_and_o = os.path.join(output, 'e_and_o')
    if dirs:
        if not os.path.isdir(e_and_o):
            os.makedirs(e_and_o)


    #Directory where plots used in html file will be stored.
    figuredir = os.path.join(output, 'plots')
    if dirs:
        if not os.path.isdir(figuredir):
            os.makedirs(figuredir)

    return output,tempdir,figuredir,e_and_o
#==============================================================================

#==============================================================================
def parse_config(srcdirectory=None, config_object=None, output=None, 
                tempdir=None, figuredir=None):
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
    outfile =  open(os.path.join(srcdirectory,'config.py'),'w')
    for key in config_object:
        for item in config_object[key]:
            outfile.write(item.upper()+'='+config_object[key][item]+'\n')

    outfile.write('OUTPUTDIR="'+output+'"\n')
    outfile.write('TEMPDIR="'+tempdir+'"\n')
    outfile.write('FIGUREDIR="'+figuredir+'"\n')
    
    outfile.close()
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
    import config
    try:
        if type(config.COMBINE) != bool:
            raise TypeError('COMBINE variable must be a boolean.')
    except NameError:
        raise NameError('COMBINE variable not found in config.ini file.')

    try:
        if type(config.COUNT) != bool:
            raise TypeError('COUNT variable must be a boolean.')
    except NameError:
        raise NameError('COUNT variable not found in config.ini file.')

    try:
        if type(config.DESEQ) != bool:
            raise TypeError('DESEQ variable must be a boolean.')
    except NameError:
        raise NameError('DESEQ variable not found in config.ini file.')

    try:
        if type(config.CALCULATE) != bool:
            raise TypeError('CALCULATE variable must be a boolean.')
    except NameError:
        raise NameError('CALCULATE variable not found in config.ini file.')

    try:
        if type(config.HOMER) != bool:
            raise TypeError('HOMER variable must be a boolean.')
    except NameError:
        raise NameError('HOMER variable not found in config.ini file.')

    try:
        if type(config.SINGLEMOTIF) != bool and type(config.SINGLEMOTIF) != str:
            raise TypeError('SINGLEMOTIF variable must be a boolean or string.')
    except NameError:
        raise NameError('SINGLEMOTIF variable not found in config.ini file.')
    try:
        if type(config.FIMO) != bool:
            raise TypeError('FIMO variable must be a boolean.')
    except NameError:
        raise NameError('FIMO variable not found in config.ini file.')

    try:
        if type(config.TEMP) != bool:
            raise TypeError('TEMP variable must be a boolean.')
    except NameError:
        raise NameError('TEMP variable not found in config.ini file.')

    try:
        if type(config.PLOT) != bool:
            raise TypeError('PLOT variable must be a boolean.')
    except NameError:
        raise NameError('PLOT variable not found in config.ini file.')

    try:
        if type(config.OUTPUT) != str:
            raise TypeError('OUTPUT variable must be a string.')
    except NameError:
        raise NameError('OUTPUT variable not found in config.ini file.')

    try:
        if type(config.BED1) != list:
            raise TypeError('BED1 variable must be a list.')
    except NameError:
        raise NameError('BED1 variable not found in config.ini file.')
    
    try:
        if type(config.BED2) != list:
            raise TypeError('BED2 variable must be a list.')
    except NameError:
        raise NameError('BED2 variable not found in config.ini file.')

    try:
        if type(config.BAM1) != list:
            raise TypeError('BAM1 variable must be a list.')
    except NameError:
        raise NameError('BAM1 variable not found in config.ini file.')
    try:
        if type(config.LABEL1) != str:
            raise TypeError('LABEL1 variable must be a string.')
    except NameError:
        raise NameError('LABEL1 variable not found in config.ini file.')

    try:
        if type(config.BAM2) != list:
            raise TypeError('BAM2 variable must be a list.')
    except NameError:
        raise NameError('BAM2 variable not found in config.ini file.')

    try:
        if type(config.LABEL2) != str:
            raise TypeError('LABEL2 variable must be a string.')
    except NameError:
        raise NameError('LABEL2 variable not found in config.ini file.')

    try:
        if type(config.PADJCUTOFF) != float:
            raise TypeError('PADJCUTOFF variable must be a float.')
    except NameError:
        raise NameError('PADJCUTOFF variable not found in config.ini file.')

    try:
        if type(config.LARGEWINDOW) != float:
            raise TypeError('LARGEWINDOW variable must be a float.')
    except NameError:
        raise NameError('LARGEWINDOW variable not found in config.ini file.')

    try:
        if type(config.SMALLWINDOW) != float:
            raise TypeError('SMALLWINDOW variable must be a float.')
    except NameError:
        raise NameError('SMALLWINDOW variable not found in config.ini file.')

    if config.FIMO != True:
        try:
            if type(config.MOTIF_GENOMEWIDE_HITS) != str:
                raise TypeError('MOTIF_GENOMEWIDE_HITS variable must be a string.')
        except NameError:
            raise NameError('MOTIF_GENOMEWIDE_HITS variable not found in config.ini file.')

    if config.HOMER != True:
        try:
            if type(config.HOMER_MOTIF_FILE) != str:
                raise TypeError('HOMER_MOTIF_FILE variable must be a string.')
        except NameError:
            raise NameError('HOMER_MOTIF_FILE variable not found in config.ini file.')

    try:
        if type(config.GENOMEFASTA) != str:
            raise TypeError('GENOMEFASTA variable must be a string.')
    except NameError:
        raise NameError('GENOMEFASTA variable not found in config.ini file.')

    try:
        if type(config.MOTIFDATABASE) != str:
            raise TypeError('MOTIFDATABASE variable must be a string.')
    except NameError:
        raise NameError('MOTIFDATABASE variable not found in config.ini file.')

    try:
        if type(config.LOGOS) != str:
            raise TypeError('LOGOS variable must be a string.')
    except NameError:
        raise NameError('LOGOS variable not found in config.ini file.')
    
    try:
        if type(config.DPI) != float and config.DPI != None:
            raise TypeError('DPI variable must be a float or None.')
    except NameError:
        raise NameError('DPI variable not found in config.ini file.')

    print "Config file verified, all inputs present and correct type."
#==============================================================================

#==============================================================================
