#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''This file is a graveyard for older functions that are no longer in use
'''
#==============================================================================
__author__ = 'Jonathan D. Rubin and Rutendo F. Sigauke'
__credits__ = ['Jonathan D. Rubin', 'Rutendo F. Sigauke', 'Jacob T. Stanley',
                'Robin D. Dowell']
__maintainer__ = 'Jonathan D. Rubin'
__email__ = 'Jonathan.Rubin@colorado.edu'
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
            raise TypeError(('COMBINE variable must be a boolean. This switch '
                            'determines whether TFEA merges bed files within '
                            'BED input.'))
    except NameError:
        raise NameError(('COMBINE variable not found in config.ini file. This '
                            'switch determines whether TFEA merges bed files '
                            'within BED input.'))

    try:
        if type(config.COUNT) != bool:
            raise TypeError(('COUNT variable must be a boolean. This switch '
                                'determines whether TFEA performs read '
                                'counting over BED regions.'))
    except NameError:
        raise NameError(('COUNT variable not found in config.ini file. This '
                            'switch determines whether TFEA performs read '
                            'counting over BED regions.'))

    try:
        if type(config.DESEQ) != bool:
            raise TypeError(('DESEQ variable must be a boolean. This switch '
                                'determines whether TFEA performs DE-Seq '
                                'analysis on counted BED regions.'))
    except NameError:
        raise NameError(('DESEQ variable not found in config.ini file. This '
                            'switch determines whether TFEA performs DE-Seq '
                            'analysis on counted BED regions.'))

    try:
        if type(config.CALCULATE) != bool:
            raise TypeError(('CALCULATE variable must be a boolean. This '
                                'switch determines whether TFEA performs its '
                                'standard enrichment score calculation and '
                                'plotting.'))
    except NameError:
        raise NameError(('CALCULATE variable not found in config.ini file. '
                            'This switch determines whether TFEA performs its '
                            'standard enrichment score calculation and '
                            'plotting.'))

    try:
        if type(config.POOL) != bool:
            raise TypeError(('POOL variable must be a boolean. This switch '
                                'determines whether TFEA runs the analysis '
                                'in parallel using the multiprocessing library '
                                'in python.'))
    except NameError:
        raise NameError(('POOL variable not found in config.ini file. This '
                            'switch determines whether TFEA runs the analysis '
                            'in parallel using the multiprocessing library in '
                            'python.'))

    try:
        if type(config.SINGLEMOTIF) != bool and type(config.SINGLEMOTIF) != str:
            raise TypeError(('SINGLEMOTIF variable must be a boolean or '
                                'string. This switch determines whether TFEA '
                                'performs its analysis on a single motif or '
                                'all. If not False, set to a string matching '
                                'a motif name.'))
    except NameError:
        raise NameError(('SINGLEMOTIF variable not found in config.ini file. '
                            'This switch determines whether TFEA performs its '
                            'analysis on a single motif or all. If not False, '
                            'set to a string matching a motif name.'))

    try:
        if type(config.FIMO) != bool:
            raise TypeError(('FIMO variable must be a boolean. This switch '
                                'determines whether TFEA uses FIMO to get '
                                'motif hits or whether a database of motif hit '
                                'calls (bed format) is used.'))
    except NameError:
        raise NameError(('FIMO variable not found in config.ini file. This '
                            'switch determines whether TFEA uses FIMO to get '
                            'motif hits or whether a database of motif hit '
                            'calls (bed format) is used.'))

    try:
        if type(config.TEMP) != bool:
            raise TypeError(('TEMP variable must be a boolean. This switch '
                                'determines whether TFEA saves large temporary '
                                'files. If True, temporary files will be '
                                'stored in the temp_files directory. Warning: '
                                'There will be many large files.'))
    except NameError:
        raise NameError(('TEMP variable not found in config.ini file. This '
                            'switch determines whether TFEA saves large '
                            'temporary files. If True, temporary files will be '
                            'stored in the temp_files directory. Warning: '
                            'There will be many large files.'))

    try:
        if type(config.PLOT) != bool:
            raise TypeError(('PLOT variable must be a boolean. This switch '
                                'determines whether TFEA outputs plots for '
                                'all motifs provided or just significant ones. '
                                'Warning: Setting this to True will slow down '
                                'TFEA and create large output folders.'))
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
        if type(config.BED1) != list:
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
        if type(config.BED2) != list:
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
def verify_config_object(config=None):
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
    config_items = [config.items(section) for section in config]
    config_items = [item for sublist in config_items for item in sublist]

    print config_items

    #pull out only names, give them all uppercase so the user input can be 
    #case-insensitive
    names = [name.upper().encode("utf-8") for (name,value) in config_items]

    #pull out corresponding value pair for each variable name
    values = [value.encode("utf-8") for (name,value) in config_items]

    #Try to find the index of a given variable
    try:
        #If variable name found, verify it's type
        combine = values[names.index('COMBINE')]
        if combine != 'True' or combine != 'False':
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
def pvalue_global_auc(TFresults=None,auc_index=None):
    '''Calculates a p-value based on a distribution of youden_ranks

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
    auc_list = [x[auc_index] for x in TFresults]
    ##significance calculator                                                                                                                                                            
    mu = np.mean(auc_list)
    sigma = np.std(auc_list)
    for i in range(len(TFresults)):
        actual_auc = TFresults[i][auc_index]
        p = min(norm.cdf(actual_auc,mu,sigma), 
                1-norm.cdf(actual_auc,mu,sigma))
        TFresults[i].append(p)

    return TFresults
#==============================================================================

#==============================================================================
def distance_distribution_plot(largewindow=None, smallwindow=None, 
                                figuredir=None, updistancehist=None, 
                                middledistancehist=None, downdistancehist=None,
                                motif_file=None):
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
    plt.savefig(os.path.join(figuredir, motif_file.split('.bed')[0] 
                + '_distance_distribution.png'),bbox_inches='tight')
    
    plt.cla()
#==============================================================================
