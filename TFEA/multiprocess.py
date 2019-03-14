#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''This module contains scripts to perform parallel processing using python's
    built-in multiprocessing module
'''

#==============================================================================
__author__ = 'Jonathan D. Rubin'
__credits__ = ['Jonathan D. Rubin', 'Rutendo F. Sigauke', 'Jacob T. Stanley',
                'Robin D. Dowell']
__maintainer__ = 'Jonathan D. Rubin'
__email__ = 'Jonathan.Rubin@colorado.edu'

#Imports
#==============================================================================
from os import getpid
import multiprocessing as mp
import subprocess
import sys
import config
from pathlib import Path
from contextlib import closing
#Main Script
#==============================================================================
def main(function=None, args=None, kwargs=None, debug=False):
    '''This is the main script of the multiprocessing module. It is written
        to simplify code in other modules. It performs parallel processing
        using python's built-in multiprocessing module on a function given
        keyword arguments shared between threads and arguments given to each
        thread

    Parameters
    ----------
    function : function object
        A python function object to process some arguments
    args : list
        A list of arguments to be input into given function
    kwargs : dict
        A dictionary of keyword arguments to input into given function
    debug : boolean
        Whether to print memory usage information of running processes

    Returns
    -------
    results : list
        A list of the result of processing given args and kwargs with given 
        function
    '''
    # #Get max number of cpus, max out at 64
    cpus = mp.cpu_count()
    if cpus > 64:
        cpus = 64
    # p = mp.Pool(processes=cpus)
    # with closing(mp.Pool(processes=cpus)) as p:
    with mp.get_context("spawn").Pool() as p:
    # p = mp.Pool()
    # results = list()
        results = p.starmap(helper, [(function, arg, kwargs, debug) for arg in args])
        p.close()
        p.join()
        
    # for i, arg in enumerate(args):
    #     results.append(p.apply_async(function, (arg,), kwargs).get())
    #     if debug:
    #         print("Process number: " + str(i), file=sys.stderr)
    #         current_mem_usage()
    # p.close()
    # p.join()
    # p.terminate()

    return results

def helper(function, arg, kwargs, debug):
    '''This function serves to unpack keyword arguments provided to this module
        and feed it in appropriately to the desired function

    Parameters
    ----------
    function : function object
        A python function object to process some arguments
    arg : object
        A positional argument to pass to function
    kwargs : dict
        keyword arguments to pass to function

    Returns
    -------
    result : object
        The result of running the function with given argument and keyword
        arguments.
    '''
    if debug:
        print("Process", getpid(), "running", function.__name__, file=sys.stderr)
        current_mem_usage()
    result = function(arg, **kwargs)

    return result

#Functions
#==============================================================================
def current_mem_usage():
    '''Prints current memory usage. Adapted from scripts by Chris Slocum and 
        Fred Cirera
    '''
    suffix = 'B'
    # num = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    job = (Path(config.TEMPDIR) / 'jobid.txt').read_text().strip('\n')
    num = subprocess.check_output(["sstat", "--format=maxRSS", "-j", job+'.batch']).decode('UTF-8').split('\n')[-2].strip()
    unit = num[-1]
    if unit == 'K':
        byte = float(num.strip('K'))*1024.0
    elif unit == 'M':
        byte = float(num.strip('M'))*1024.0*1024.0
    elif unit == 'G':
        byte = float(num.strip('G'))*1024.0*1024.0*1024.0
    for unit in ['','K','M','G','T','P','E','Z']:
        if abs(byte) < 1024.0:
            print("Memory Usage: %3.1f%s%s" % (byte, unit, suffix), file=sys.stderr)
            return
        byte /= 1024.0
    print("Memory Usage: %.1f%s%s" % (byte, 'Y', suffix), file=sys.stderr)
    return