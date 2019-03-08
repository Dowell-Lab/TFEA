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
import resource
import multiprocessing as mp
import sys
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
    #Get max number of cpus, max out at 64
    cpus = mp.cpu_count()
    if cpus > 64:
        cpus = 64
    
    p = mp.Pool(processes=cpus)
    results = list()
    for i, arg in enumerate(args):
        results.append(p.apply_async(function, (arg,), kwargs).get())
        if debug:
            print("Process number: " + str(i), file=sys.stderr)
            current_mem_usage()
    p.close()
    p.join()

    return results

#Functions
#==============================================================================
def current_mem_usage():
    '''Prints current memory usage. Adapted from scripts by Chris Slocum and 
        Fred Cirera
    '''
    suffix = 'B'
    num = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.
    for unit in ['','K','M','G','T','P','E','Z']:
        if abs(num) < 1024.0:
            print("Memory Usage: %3.1f%s%s" % (num, unit, suffix), file=sys.stderr)
            return
        num /= 1024.0
    print("Memory Usage: %.1f%s%s" % (num, 'Y', suffix), file=sys.stderr)
    return