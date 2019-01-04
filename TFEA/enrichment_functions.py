#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''This file contains a list of functions associated with outputting a ranked
    bed file based on some metric.
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
import numpy as np

#Functions
#==============================================================================
def area_under_curve(distances, smallwindow=None, output_type=None, 
                        permutations=None):
    #sort distances based on the ranks from TF bed file
    #and calculate the absolute distance
    motif = distances[0]
    distances = distances[1:]
    distances_abs = [abs(x) for x in distances] 

    #Get -exp() of distance and get cumulative scores
    #Filter distances into quartiles to get middle distribution
    q1 = round(np.percentile(np.arange(1, len(distances_abs),1), 25))
    q3 = round(np.percentile(np.arange(1, len(distances_abs),1), 75))
    middledistancehist =  [x for x in distances_abs[int(q1):int(q3)] if x <= smallwindow]
    average_distance = float(sum(middledistancehist))/float(len(middledistancehist))
    
    score = [math.exp(-float(x)/average_distance) if x <= smallwindow else 0.0 for x in distances_abs]
    total = float(sum(score))
    

    #Filter any TFs/files without any hits
    if total == 0:
        return "no hits"

    binwidth = 1.0/float(len(distances_abs))
    normalized_score = [(float(x)/total)*binwidth for x in score]
    cumscore = np.cumsum(normalized_score)

    trend = np.arange(0,1,1.0/float(len(ranks)))
    trend = [x*binwidth for x in trend]

    #The AUC is the relative to the "random" line
    actualES = np.trapz(cumscore) - np.trapz(trend)

    #Calculate random AUC
    simES = permute_auc(distances=normalized_score, trend=trend, 
                        permutations=permutations)

    ##significance calculator                                                                                                                                                            
    mu = np.mean(simES)
    NES = actualES/abs(mu)
    sigma = np.std(simES)

    p = min(norm.cdf(actualES,mu,sigma), 1-norm.cdf(actualES,mu,sigma))

    if math.isnan(p):
        p = 1.0

    if output_type == 'html':
        import plotting_functions
        plotting_score = [(float(x)/total) for x in score]
        plotting_cumscore = np.cumsum(plotting_score)

#==============================================================================
def permute_auc(distances=None, trend=None, permutations=None):
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
    triangle_area = np.trapz(trend)
    for i in range(permutations):
        random_distances = np.random.permutation(distances)
        cum_distances = np.cumsum(random_distances)
        es = np.trapz(cum_distances)
        auc = es - triangle_area
        es_permute.append(auc)

    return es_permute