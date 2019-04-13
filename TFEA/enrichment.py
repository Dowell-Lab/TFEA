#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''This module calculates enrichment statistics for TF motifs across inputted
    regions of interest. This module takes as input the output of the SCANNER
    module. It outputs a list of lists with enrichment statistics including
    a p-value.
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
import math
import time
import datetime
import traceback
import numpy as np
import math
from scipy.stats import norm
from scipy.stats import anderson_ksamp

import config
import multiprocess

#Main Script
#==============================================================================
def main(motif_distances=config.vars.MOTIF_DISTANCES, 
        md_distances1=config.vars.MD_DISTANCES1, md_distances2=config.vars.MD_DISTANCES2, 
        mdd_distances1=config.vars.MDD_DISTANCES1, mdd_distances2=config.vars.MDD_DISTANCES2, 
        enrichment=config.vars.ENRICHMENT, output_type=config.vars.OUTPUT_TYPE, 
        permutations=config.vars.PERMUTATIONS, debug=config.vars.DEBUG, 
        largewindow=config.vars.LARGEWINDOW, smallwindow=config.vars.SMALLWINDOW, 
        md=config.vars.MD):
    '''This is the main script of the ENRICHMENT module. It takes as input
        a list of distances outputted from the SCANNER module and calculates
        an enrichment score, a p-value, and in some instances an adjusted 
        p-value for each motif.
    
    Parameters
    ----------
    motif_distances : list of lists
        A list containing a list for each motif scanned. For each motif, the 
        list begins with the motif name as a string and is followed by int
        values corresponding to the motif distance for each region (ranked). A
        '.' value means the motif was not within the given region
    md_distances1 : list of lists
        A list containing a list for each motif scanned. For each motif, the 
        list begins with the motif name as a string and is followed by int
        values corresponding to the motif distance for each region (ranked). A
        '.' value means the motif was not within the given region
    md_distances2 : list of lists
        A list containing a list for each motif scanned. For each motif, the 
        list begins with the motif name as a string and is followed by int
        values corresponding to the motif distance for each region (ranked). A
        '.' value means the motif was not within the given region
    enrichment : str
        The type of enrichment analysis to perform
    output_type : str
        Determines what some functions will output. At this point, this is mostly
        intended for debug purposes.
    permutations : int
        Number of random shuffling permutations to perform to calculate a 
        p-value
    debug : boolean
        Whether to print debug statements specifically within the multiprocess
        module
    largewindow : int
        A distance cutoff value used within auc_bgcorrect
    smallwindow : int
        A distance cutoff value used within the md score analysis
    
    Returns
    -------
    results : list of lists
        A list of lists corresponding to enrichment statistics for each motif
    md_results : list of lists
        A list of lists corresponding to md-score statistics for each motif
    '''
    ENRICHMENTtime = time.time()
    print("Calculating enrichment...", end=' ', flush=True, file=sys.stderr)
    if enrichment == 'auc':
        auc_keywords = dict(output_type=output_type, 
                            permutations=permutations)
        results = multiprocess.main(function=area_under_curve, 
                                    args=motif_distances, kwargs=auc_keywords,
                                    debug=debug)

        padj_bonferroni(results)
    elif enrichment == 'anderson-darling':
        p = Pool(processes=cpus, maxtasksperchild=10)
        results = list()
        for distances in motif_distances:
            results.append(p.apply_async(anderson_darling, 
                            (distances,)).get())
        p.close()
        p.join()

    elif enrichment == 'auc_bgcorrect':
        auc_bgcorrect_keywords = dict(output_type=output_type, 
                                        permutations=permutations,
                                        largewindow=largewindow)
        results = multiprocess.main(function=area_under_curve_bgcorrect, 
                                    args=motif_distances, 
                                    kwargs=auc_bgcorrect_keywords,
                                    debug=debug)

        padj_bonferroni(results)

    if md:
        md_results = md(md_distances1=md_distances1, md_distances2=md_distances2)
        config.vars.MD_RESULTS = md_results
    if mdd:
        mdd_results =  md(md_distances1=mdd_distances1, md_distances2=mdd_distances2)
        config.vars.MDD_RESULTS = mdd_results


    config.vars.RESULTS = results


    ENRICHMENTtime = time.time()-ENRICHMENTtime
    print("done in: " + str(datetime.timedelta(seconds=int(SCANNERtime))), file=sys.stderr)

    if config.vars.DEBUG:
        multiprocess.current_mem_usage()

#==============================================================================
def md(md_distances1=None, md_distances2=None, output_type=config.vars.OUTPUT_TYPE, 
        debug=config.vars.DEBUG, smallwindow=config.vars.SMALLWINDOW):
    md_keywords = dict(smallwindow=smallwindow)
    md_results = multiprocess.main(function=md_score, 
                    args=zip(sorted(md_distances1),sorted(md_distances2)), 
                    kwargs=md_keywords,
                    debug=debug)

    md_results = md_score_p(md_results)

    return md_results

#Functions
#==============================================================================
def area_under_curve(distances, output_type=None, permutations=None):
    '''Calculates an enrichment score using the area under the curve. This
        method is not as sensitive to artifacts as other methods. It works well
        as an asymmetry detector and will be good at picking up cases where
        most of the motif localization changes happen at the most differentially
        transcribed regions.
    '''
    try:
        #sort distances based on the ranks from TF bed file
        #and calculate the absolute distance
        motif = distances[0]
        distances = distances[1:]
        distances_abs = [abs(x)  if x != '.' else x for x in distances]

        hits = len([x for x in distances_abs if x != '.'])

        #Filter any TFs/files without any hits
        if hits == 0:
            return [motif, 0.0, 0, 1.0]

        #Get -exp() of distance and get cumulative scores
        #Filter distances into quartiles to get middle distribution
        q1 = int(round(len(distances)*.25))
        q3 = int(round(len(distances)*.75))
        middledistancehist =  [x for x in distances_abs[int(q1):int(q3)] if x != '.']
        try:
            average_distance = float(sum(middledistancehist))/float(len(middledistancehist))
        except ZeroDivisionError:
            return [motif, 0.0, 0, 1.0]
        
        score = [math.exp(-float(x)/average_distance) if x != '.' else 0.0 for x in distances_abs]
        total = sum(score)
        if output_type == 'score':
            return score

        binwidth = 1.0/float(len(distances_abs))
        normalized_score = [(float(x)/total)*binwidth for x in score]
        cumscore = np.cumsum(normalized_score)
        trend = np.append(np.arange(0,1,1.0/float(len(cumscore)))[1:], 1.0)
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
            from TFEA import plotting_functions
            plotting_score = [(float(x)/total) for x in score]
            plotting_cumscore = np.cumsum(plotting_score)
            plotting_functions.plot_individual_graphs(plot=config.vars.PLOTALL, 
                                padj_cutoff=config.vars.PADJCUTOFF,
                                figuredir=config.vars.FIGUREDIR, 
                                logos=config.vars.MOTIF_LOGOS, 
                                largewindow=config.vars.LARGEWINDOW, 
                                score=plotting_score, 
                                smallwindow=None,
                                distances_abs=None, sorted_distances=None,
                                ranks=None, pvals=None, fc=None, 
                                cumscore=None, motif_file=None, p=None,
                                simES=None, actualES=None, gc_array=None,
                                meta_profile_dict=None, label1=None, label2=None)
        
        if output_type == 'lineplot_output':
            plotting_score = [0] + [(float(x)/total) for x in score]
            plotting_cumscore = np.cumsum(plotting_score)
            return plotting_cumscore
        
        if output_type == 'simulation_output':
            return auc, sim_auc
    except Exception as e:
        # This prints the type, value, and stack trace of the
        # current exception being handled.
        print(traceback.print_exc())
        raise e
    return [motif, auc, hits, p]

#==============================================================================
def area_under_curve_bgcorrect(distances, output_type=None, permutations=None, 
                                largewindow=None):
    #sort distances based on the ranks from TF bed file
    #and calculate the absolute distance
    motif = distances[0]
    distances = distances[1:]
    distances_abs = [abs(x)  if x != '.' else x for x in distances]

    hits = len([x for x in distances_abs if x != '.'])

    #Filter any TFs/files without any hits
    if hits == 0:
        return [motif, 0.0, 0, 1.0]

    #Get -exp() of distance and get cumulative scores
    #Filter distances into quartiles to get middle distribution
    q1 = int(round(len(distances)*.25))
    q3 = int(round(len(distances)*.75))
    middledistancehist = [x for x in distances_abs[int(q1):int(q3)] if x != '.']
    expectation = np.histogram(middledistancehist, bins=largewindow)
    expectation_sum = np.sum(expectation)
    expectation = [x/expectation_sum for x in expectation]
    
    score = [1/expectation[x] if x != '.' else 0.0 for x in distances_abs]
    total = sum(score)
    if output_type == 'score':
        return score

    binwidth = 1.0/float(len(distances_abs))
    normalized_score = [(float(x)/total)*binwidth for x in score]
    cumscore = np.cumsum(normalized_score)
    trend = np.append(np.arange(0,1,1.0/float(len(cumscore)))[1:], 1.0)
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
        from TFEA import plotting_functions
        plotting_score = [(float(x)/total) for x in score]
        plotting_cumscore = np.cumsum(plotting_score)
        plotting_functions.plot_individual_graphs(plot=None, padj_cutoff=None,
                            figuredir=None, logos=None, 
                            largewindow=None, score=None, 
                            smallwindow=None,
                            distances_abs=None, sorted_distances=None,
                            ranks=None, pvals=None, fc=None, 
                            cumscore=None, motif_file=None, p=None,
                            simES=None, actualES=None, gc_array=None,
                            meta_profile_dict=None, label1=None, label2=None)
    
    if output_type == 'lineplot_output':
        plotting_score = [(float(x)/total) for x in score]
        plotting_cumscore = np.cumsum(plotting_score)
        return plotting_cumscore
    
    if output_type == 'simulation_output':
        return auc, sim_auc

    return [motif, auc, hits, p]

#==============================================================================
def max_GSEA(distances, output_type=None, permutations=None, cutoff=None):
    '''Takes the max of the enrichment curve. This method is good at detecting
        significant motif instances in a row. It is susceptible to instances
        where a motif appears generally throughout all regions.
    '''
    motif = distances[0]
    distances = [x if x != '.' else cutoff+1 for x in distances[1:]]
    hits = [1 if abs(x) < cutoff else -1 for x in distances]
    pos_total = float(sum([1 for x in hits if x != -1]))
    neg_total = float(sum([1 for x in hits if x == -1]))

    Eval = 0
    ES = list()
    for distance in hits:
        if distance == -1:
            Eval += -1.0/neg_total
            ES.append(Eval)
        else:
            Eval += 1.0/pos_total
            ES.append(Eval)
    
    actualES = max(ES, key=abs)

    #To get NES, first simulate 1000 permuations of region ranks
    simES = permute_max_GSEA(hits=hits, permutations=1000)

    #NES is the actualES divided by the mean ES of all permutations with the 
    # same sign as actualES p-value is caluclated empirically 
    # (i.e. (# of simulated ES scores larger than actualES)/(rest of simulated ES scores))
    if actualES < 0:
        simESsubset = [x for x in simES if x < 0]
        mu = np.mean(simESsubset)
        NES = -(actualES/mu)
        sigma = np.std(simESsubset)
        p = norm.cdf(actualES, mu, sigma)
    else:
        simESsubset = [x for x in simES if x > 0]
        mu = np.mean(simESsubset)
        NES = actualES/mu
        sigma = np.std(simESsubset)
        p = 1-norm.cdf(actualES, mu, sigma)

    if math.isnan(p):
        p = 1.0

    if output_type == 'html':
        from TFEA import plotting_functions
        plotting_functions.plot_individual_graphs(plot=None, padj_cutoff=None,
                            figuredir=None, logos=None, 
                            largewindow=None, score=None, 
                            smallwindow=None,
                            distances_abs=None, sorted_distances=None,
                            ranks=None, pvals=None, fc=None, 
                            cumscore=None, motif_file=None, p=None,
                            simES=None, actualES=None, gc_array=None,
                            meta_profile_dict=None, label1=None, label2=None)
    
    if output_type == 'lineplot_output':
        return ES
    
    if output_type == 'simulation_output':
        return actualES, simES

    return [motif, actualES, NES, pos_total, p]

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
def permute_max_GSEA(hits=None, permutations=1000):
    import numpy as np
    #Simulate 1000 permuations of region ranks
    pos_total = float(sum([1 for x in hits if x != -1]))
    neg_total = float(sum([1 for x in hits if x == -1]))
    simES = list()
    for i in range(permutations):
        np.random.shuffle(hits)
        Eval = 0
        ES = list()
        for distance in hits:
            if distance == -1:
                Eval += -1.0/neg_total
                ES.append(Eval)
            else:
                Eval += 1.0/pos_total
                ES.append(Eval)
        simES.append(max(ES, key=abs))

    return simES

#==============================================================================
def padj_bonferroni(results, pvalindex=-1):
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
        PVAL = results[i][pvalindex]
        #Using Bonferroni Correction
        PADJ = 1 if PVAL*len(results) > 1 else PVAL*len(results)
        results[i].append(PADJ)

    return results

#==============================================================================
def anderson_darling(distances):
    #sort distances based on the ranks from TF bed file
    #and calculate the absolute distance
    motif = distances[0]
    distances = distances[1:]
    hits = sum([1 if x != '.' else 0 for x in distances])

    #Get -exp() of distance and get cumulative scores
    #Filter distances into quartiles to get middle distribution
    q1 = int(round(len(distances)*.25))
    q3 = int(round(len(distances)*.75))
    expected_distribution =  [x for x in distances[int(q1):int(q3)] if x != '.']
    observed_distribution = [x for x in distances if x != '.']

    stat, _, sig = anderson_ksamp([observed_distribution, expected_distribution])

    return [motif, stat, hits, sig]

#==============================================================================
def md_score(distances, smallwindow=None):
    '''
    '''
    distances1, distances2 = distances
    motif = distances1[0]
    distances1 = [abs(float(x)) for x in distances1[1:] if x != '.']
    distances2 = [abs(float(x)) for x in distances2[1:] if x != '.']
    d1_total = float(len(distances1))
    d2_total = float(len(distances2))
    # print(motif)
    # print("distances1: ", distances1[:100])
    # print("md1sum: ", sum([1.0 if d <= smallwindow else 0.0 for d in distances1]))
    # print("md1tot: ", d1_total)
    # print("distances2: ", distances2[:100])
    # print("md2sum: ", sum([1.0 if d <= smallwindow else 0.0 for d in distances2]))
    # print("md2tot: ", d2_total)
    try:
        md1 = sum([1.0 if d <= smallwindow else 0.0 for d in distances1])/d1_total
        md2 = sum([1.0 if d <= smallwindow else 0.0 for d in distances2])/d2_total
    except ZeroDivisionError:
        return [motif, 0, 0, 1, 1]

    return [motif, md1, md2, d1_total, d2_total]

#==============================================================================
def md_score_p(results):
    '''
    '''
    mean = float(sum([x[2]-x[1] for x in results]))/float(len(results))
    for i in range(len(results)):
        motif = results[i][0]
        md1=float(results[i][1])
        md2=float(results[i][2])
        N1=float(results[i][3])
        N2=float(results[i][4])
        if (N1+N2)/2.0 != 0:
            p=((md1*N1)+(md2*N2))/(N1+N2)
            SE=(p*(1-p))*((1/N1)+(1/N2))
            try:
                z = (md2-md1-mean)/math.sqrt(SE)
            except ZeroDivisionError:
                z=0
            cdf=norm.cdf(z,0,1)
            p=min(cdf,1-cdf)*2
            md = md2-md1-mean
            total = N1 + N2
        
        results[i] = [motif, md, total, p]

    return results

#TESTS
#==============================================================================
if __name__ == "__main__":
    import numpy as np
    from TFEA import enrichment_functions

    distances = [x for x in range(900, 1000)] + ['.' for x in range(40)]
    np.random.shuffle(distances)
    distances = [1 for x in range(20)] + distances
    distances = ['test'] + distances
    # results = enrichment_functions.max_GSEA(distances, output_type=None, 
    #                                         permutations=1000, cutoff=150)
    # print(results)

    # results = enrichment_functions.area_under_curve(distances, 
    #                                             output_type='score', 
    #                                             permutations=1000)
    

    distances_reverse = ['reverse'] + [x for x in reversed(distances[1:])]

    # results_rev = enrichment_functions.area_under_curve(distances, 
                                                # output_type='score', 
                                                # permutations=1000)

    # results = anderson_darling(distances)
    # results_rev = anderson_darling(distances_reverse)
    # print(results, results_rev)

    results = md_score([distances, distances], smallwindow=150)
    # results = ['motif', 0.1, 0.2, 1500, 1500]
    # results = md_score_p([results])
    print(results)


    # # results_rev = [x for x in reversed(results_rev)]
    # for i in range(len(results)):
    #     if results[i] != results_rev[-(i+1)]:
    #         print(results[i], results_rev[-(i+1)])

    # print(np.sum(results), np.sum(results_rev), np.sum([x for x in reversed(results)]))

    # print(results)
    # print(results_rev)
    # print(results == [x for x in reversed(results_rev)])
    # print(set(results) - set([x for x in reversed(results_rev)]))
    # print(set([x for x in reversed(results_rev)]) - set(results))

    # print(lineplot_output)