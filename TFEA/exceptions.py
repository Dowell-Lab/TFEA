#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''This file contains exceptions used within TFEA
'''

#==============================================================================
__author__ = 'Jonathan D. Rubin and Rutendo F. Sigauke'
__credits__ = ['Jonathan D. Rubin', 'Rutendo F. Sigauke', 'Jacob T. Stanley',
                'Robin D. Dowell']
__maintainer__ = 'Jonathan D. Rubin'
__email__ = 'Jonathan.Rubin@colorado.edu'

#Exceptions
#==============================================================================
class Error(Exception):
   """Base class for other exceptions"""
   pass

#==============================================================================
class FileEmptyError(Error):
   """Raised when resulting file is empty"""
   pass

#==============================================================================
class InputError(Error):
   """Raised when the input value is incorrect type or missing"""
   pass

#==============================================================================
class OutputError(Error):
   """Raised when the output of a module or function is incorrect"""
   pass

#==============================================================================
class SubprocessError(Error):
   """Raised when there is an error with a subprocess call"""
   pass