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
def area_under_curve(distances, output_type=None, permutations=None):
    #sort distances based on the ranks from TF bed file
    #and calculate the absolute distance
    motif = distances[0]
    distances = distances[1:]
    distances_abs = [abs(x) for x in distances] 

    #Get -exp() of distance and get cumulative scores
    #Filter distances into quartiles to get middle distribution
    q1 = round(np.percentile(np.arange(1, len(distances_abs),1), 25))
    q3 = round(np.percentile(np.arange(1, len(distances_abs),1), 75))
    middledistancehist =  [x for x in distances_abs[int(q1):int(q3)] if x != '.']
    average_distance = float(sum(middledistancehist))/float(len(middledistancehist))
    
    score = [math.exp(-float(x)/average_distance) if x != '.' else 0.0 for x in distances_abs]
    hits = len(score)
    total = float(sum(score))
    

    #Filter any TFs/files without any hits
    if total == 0:
        return [motif, 0.0, 0, 1.0]

    binwidth = 1.0/float(len(distances_abs))
    normalized_score = [(float(x)/total)*binwidth for x in score]
    cumscore = np.cumsum(normalized_score)

    trend = np.arange(0,1,1.0/float(len(ranks)))
    trend = [x*binwidth for x in trend]

    #The AUC is the relative to the "random" line
    auc = np.trapz(cumscore) - np.trapz(trend)

    #Calculate random AUC
    sim_auc = permute_auc(distances=normalized_score, trend=trend, 
                        permutations=permutations)

    #Calculate p-value                                                                                                                                                          
    mu = np.mean(sim_auc)
    sigma = np.std(sim_auc)
    p = min(norm.cdf(auc,mu,sigma), 1-norm.cdf(auc,mu,sigma))
    if math.isnan(p):
        p = 1.0

    if output_type == 'html':
        import output_functions
        plotting_score = [(float(x)/total) for x in score]
        plotting_cumscore = np.cumsum(plotting_score)
        output_functions.plot_individual_graphs(plot=None, padj_cutoff=None,
                            figuredir=None, logos=None, 
                            largewindow=None, score=None, 
                            smallwindow=None,
                            distances_abs=None, sorted_distances=None,
                            ranks=None, pvals=None, fc=None, 
                            cumscore=None, motif_file=None, p=None,
                            simES=None, actualES=None, gc_array=None,
                            meta_profile_dict=None, label1=None, label2=None)
    
    return [motif, auc, hits, p]

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

#==============================================================================
def padj_bonferroni(results):
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
    for i in range(len(results)):
        PVAL = results[i][-1]
        #Using Bonferroni Correction
        PADJ = 1 if PVAL*len(results) > 1 else PVAL*len(results)
        results[i].append(PADJ)

    return results