#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''This file contains functions related to performing enrichment analysis
'''
#==============================================================================
__author__ = 'Jonathan D. Rubin and Rutendo F. Sigauke'
__credits__ = ['Jonathan D. Rubin', 'Rutendo F. Sigauke', 'Jacob T. Stanley',
                'Robin D. Dowell']
__maintainer__ = 'Jonathan D. Rubin'
__email__ = 'Jonathan.Rubin@colorado.edu'
__version__ = '4.0'
#==============================================================================
#Imports
#==============================================================================
import math
import traceback
import numpy as np
from scipy.stats import norm
#==============================================================================
#Functions
#==============================================================================
def AUC_enrichment(ranked_list=None, cutoff=None):
    '''This function calculates the AUC for any TF based on the TF motif hits 
        relative to the bidirectionals. The calculated AUC is used as a proxy 
        for the enrichemnt of the TFs.

    Parameters
    ----------
    ranked_list : list or array
        a ranked list of ints or floats corresponding to motif distances

    cutoff : int or float
        a cutoff value after which, the cumulative sum to add should be 0

    Returns
    -------
    AUC : float
        the area under the curve (AUC) for a given tf

    p-value : float
        theoretical p-value calculated by comparing the observed value  to the
        distribution of all simulations (default:1000)
                 
    Raises
    ------
    ValueError
        when tf does not have any hits
    '''
    try:
        #Get -exp() of distance and get cumulative scores
        #Filter distances into quartiles to get middle distribution
        q1 = round(np.percentile(np.arange(1, len(ranked_list),1), 25))
        q3 = round(np.percentile(np.arange(1, len(ranked_list),1), 75))
        middledistancehist =  [x for x in ranked_list[int(q1):int(q3)]]
        average_distance = float(sum(middledistancehist))/float(len(middledistancehist))
        score = [math.exp(-float(x)/average_distance) if x <= cutoff else 0.0 for x in ranked_list]
        print(np.cumsum(score))
        total = float(sum(score))
        

        #Filter any TFs/files without any hits
        if total == 0:
            return "no hits"

        binwidth = 1.0/float(len(ranked_list))
        normalized_score = [(float(x)/total)*binwidth for x in score]
        cumscore = np.cumsum(normalized_score)

        #np.arange returns a closed, open interval (ie. would not include
        #1). Therefore the list is shifted by 1 and a 1.0 is added to the end
        trend = list(np.arange(0,1,1.0/float(len(score)))[1:]) + [1.0]
        trend = [x*binwidth for x in trend]

        #The AUC is the relative to the "random" line
        auc = np.trapz(cumscore) - np.trapz(trend)

        sim_auc = [x for x in AUC_permute(ranked_list=ranked_list, trend=trend)]

        ##significance calculator                                                                                                                                                            
        mu = np.mean(sim_auc)
        sigma = np.std(sim_auc)

        p = min(norm.cdf(auc,mu,sigma), 1-norm.cdf(auc,mu,sigma))

        if math.isnan(p):
            p = 1.0

    except Exception as e:
        # This prints the type, value, and stack trace of the
        # current exception being handled.
        print(traceback.print_exc())
        raise e

    return auc, p
#==============================================================================

#==============================================================================
def AUC_permute(ranked_list=None, trend=None, permutations=1000):
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
    triangle_area = np.trapz(trend)
    for i in range(permutations):
        random_list = np.random.permutation(ranked_list)
        cum_values = np.cumsum(random_list)
        es = np.trapz(cum_values)
        auc = es - triangle_area

    yield auc
#==============================================================================
#Testing
#==============================================================================
if __name__ == "__main__":
    ranked_list = [1,1,1,1,1]
    auc, p = AUC_enrichment(ranked_list=ranked_list, cutoff=5)
    print(auc, p)