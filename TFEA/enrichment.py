#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''This module calculates enrichment statistics for TF motifs across inputted
    regions of interest. This module takes as input the output of the SCANNER
    module. It outputs a list of lists with enrichment statistics including
    a p-value.
'''

#==============================================================================
__author__ = ['Jonathan D. Rubin', 'Rutendo F. Sigauke', 'Hope A. Townsend']
__credits__ = ['Jonathan D. Rubin', 'Rutendo F. Sigauke', 'Jacob T. Stanley',
               'Hope A. Townsend', 'Robin D. Dowell']
__maintainer__ = 'Jonathan D. Rubin'
__email__ = 'Jonathan.Rubin@colorado.edu'

#Imports
#==============================================================================
import os
import sys
import math
import time
import datetime
import traceback
import numpy as np
import math
import pandas
import pathlib
import ujson
import shutil
from scipy import stats
from scipy.interpolate import splrep, splev
import re
from multiprocessing import Manager

from TFEA import config
from TFEA import multiprocess
from TFEA import plot
from TFEA import exceptions

#Main Script
#==============================================================================
def main(use_config=True, motif_distances=None, md_distances1=None, 
            md_distances2=None, mdd_distances1=None, mdd_distances2=None, 
            enrichment=None, output_type=None, permutations=None, debug=None, 
            largewindow=None, smallwindow=None, md=None, mdd=None, cpus=None, 
            jobid=None, pvals=None, fcs=None, p_cutoff=None, figuredir=None, 
            plotall=False, fimo_motifs=None, meta_profile_dict=None, 
            label1=None, label2=None, dpi=None, motif_fpkm={}, bootstrap=False,
            gc=None, plot_format=None):
    '''This is the main script of the ENRICHMENT module. It takes as input
        a list of distances outputted from the SCANNER module and calculates
        an enrichment score, a p-value, and in some instances an adjusted 
        p-value for each motif.
    
    Parameters
    ----------
    use_config : boolean
        Whether to use a config module to assign variables.
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
    start_time = time.time()
    if use_config:
        motif_distances = config.vars['MOTIF_DISTANCES']
        md_distances1 = config.vars['MD_DISTANCES1']
        md_distances2 = config.vars['MD_DISTANCES2']
        mdd_distances1 = config.vars['MDD_DISTANCES1']
        mdd_distances2 = config.vars['MDD_DISTANCES2']
        enrichment = config.vars['ENRICHMENT']
        permutations = config.vars['PERMUTATIONS']
        debug = config.vars['DEBUG']
        largewindow = config.vars['LARGEWINDOW']
        smallwindow = config.vars['SMALLWINDOW']
        ranked_file = config.vars['RANKED_FILE']
        pvals = config.vars['PVALS']
        fcs = config.vars['FCS']
        md = config.vars['MD']
        mdd = config.vars['MDD']
        cpus = config.vars['CPUS']
        jobid = config.vars['JOBID']
        p_cutoff = np.log(config.vars['PADJCUTOFF'])
        figuredir = config.vars['FIGUREDIR']
        plotall = config.vars['PLOTALL']
        fimo_motifs = config.vars['FIMO_MOTIFS']
        meta_profile_dict = config.vars['META_PROFILE']
        label1 = config.vars['LABEL1']
        label2 = config.vars['LABEL2']
        output_type = config.vars['OUTPUT_TYPE']
        bootstrap = config.vars['BOOTSTRAP']
        gc = config.vars['GC']
        plot_format = config.vars['PLOT_FORMAT']
        try:
            motif_fpkm = config.vars['MOTIF_FPKM']
        except:
            motif_fpkm = {}

    print("Calculating enrichment...", flush=True, file=sys.stderr)

    results = None
    md_results = None
    mdd_results = None

    if enrichment == 'auc':
        gc_correct = {}
        linear_regression = None
        # If no gc correction, still calculate the regression to graph
        if gc:
            print('\tCorrecting GC:', file=sys.stderr)
            auc_keywords = dict(fimo_motifs=fimo_motifs)
            motif_gc_auc = multiprocess.main(function=get_auc_gc, 
                                        args=motif_distances, kwargs=auc_keywords,
                                        debug=debug, jobid=jobid, cpus=cpus)

            #Calculate linear regression based on AUC and GC content of motifs
            varx = np.array([i[2] for i in motif_gc_auc])
            vary = np.array([i[1] for i in motif_gc_auc])
            mask = ~np.isnan(varx) & ~np.isnan(vary)
            linear_regression = [x for x in stats.linregress(varx[mask], vary[mask])]
            slope, intercept, _, _, _ = linear_regression
            for key, _, gc in motif_gc_auc:
                offset = slope*gc + intercept
                gc_correct[key] = offset
        else:
            print('\tNOT correcting for GC:', file=sys.stderr)
            auc_keywords = dict(fimo_motifs=fimo_motifs)
            motif_gc_auc = multiprocess.main(function=get_auc_gc, 
                                        args=motif_distances, kwargs=auc_keywords,
                                        debug=debug, jobid=jobid, cpus=cpus)

            #Calculate linear regression based on AUC and GC content of motifs
            varx = np.array([i[2] for i in motif_gc_auc])
            vary = np.array([i[1] for i in motif_gc_auc])
            mask = ~np.isnan(varx) & ~np.isnan(vary)
            linear_regression = [x for x in stats.linregress(varx[mask], vary[mask])]
            for key, _, gc in motif_gc_auc:
                gc_correct[key] = 0


        print('\tCalculating E-Score:', file=sys.stderr)
        # manager = Manager()
        # meta_profile_dict = manager.dict(meta_profile_dict)
        auc_keywords = dict(permutations=permutations, use_config=use_config, 
                        output_type=output_type, pvals=pvals, plotall=plotall, 
                        p_cutoff=p_cutoff, figuredir=figuredir, 
                        largewindow=largewindow, fimo_motifs=fimo_motifs, 
                        meta_profile_dict=meta_profile_dict, label1=label1, 
                        label2=label2, fcs=fcs, motif_fpkm=motif_fpkm, 
                        tests=len(motif_distances), bootstrap=bootstrap, 
                        gc_correct=gc_correct, plot_format=plot_format, 
                           ranked_file=ranked_file)
        results = multiprocess.main(function=auc_simulate_and_plot, 
                                    args=motif_distances, kwargs=auc_keywords,
                                    debug=debug, jobid=jobid, cpus=cpus)
                                    
        if gc:
            plot.plot_global_gc(results, p_cutoff=p_cutoff, 
                                    title='TFEA GC-Plot', 
                                    xlabel='Motif GC-content',
                                    ylabel='Non-corrected E-Score', 
                                    savepath=figuredir / ('TFEA_GC.' + plot_format), 
                                    linear_regression=linear_regression,
                                    plot_format=plot_format, 
                                    x_index=4,
                                    y_index=1, 
                                    c_index=2,
                                    p_index=-1,
                                    ylimits=[-1,1])
        else:
            plot.plot_global_gc(results, p_cutoff=p_cutoff, 
                                    title='TFEA GC-Plot (BUT correction not applied)', 
                                    xlabel='Motif GC-content',
                                    ylabel='Non-corrected E-Score', 
                                    savepath=figuredir / ('TFEA_GC.' + plot_format), 
                                    linear_regression=linear_regression,
                                    plot_format=plot_format, 
                                    x_index=4,
                                    y_index=1, 
                                    c_index=2,
                                    p_index=-1,
                                    ylimits=[-1,1])

        
        # results = list()
        # for motif_distance in motif_distances:
        #     results.append(area_under_curve(motif_distance, **auc_keywords))
        # padj_bonferroni(results)
    # elif enrichment == 'anderson-darling':
    #     results = multiprocess.main(function=anderson_darling, 
    #                                     args=motif_distances, debug=debug, 
    #                                     jobid=jobid, cpus=cpus)

    # elif enrichment == 'auc_bgcorrect':
    #     print('\tTFEA:', file=sys.stderr)
    #     auc_bgcorrect_keywords = dict(permutations=permutations)
    #     results = multiprocess.main(function=area_under_curve_bgcorrect, 
    #                                 args=motif_distances, 
    #                                 kwargs=auc_bgcorrect_keywords,
    #                                 debug=debug, jobid=jobid, cpus=cpus)

    #     padj_bonferroni(results)
    else:
        raise exceptions.InputError("Enrichment option not recognized or supported.")

    if md:
        print('\tMD:', file=sys.stderr)
        md_results = calculate_md(md_distances1=md_distances1, 
                                    md_distances2=md_distances2, 
                                    smallwindow=smallwindow, 
                                    jobid=jobid, cpus=cpus, debug=debug)
        if use_config:
            config.vars['MD_RESULTS'] = md_results
    if mdd:
        print('\tMDD:', file=sys.stderr)
        mdd_results =  calculate_md(md_distances1=mdd_distances1, 
                                        md_distances2=mdd_distances2, 
                                        smallwindow=smallwindow,
                                        jobid=jobid, cpus=cpus, debug=debug)
        if use_config:
            config.vars['MDD_RESULTS'] = mdd_results

    if use_config:
        config.vars['RESULTS'] = results

    total_time = time.time() - start_time
    if use_config:
        config.vars['ENRICHMENTtime'] = total_time

    #Remove large meta profile file
    # meta_profile_file.unlink()
    if type(meta_profile_dict) == pathlib.PosixPath:
        shutil.rmtree(meta_profile_dict, ignore_errors=True)

    print("done in: " + str(datetime.timedelta(seconds=int(total_time))), file=sys.stderr)

    if debug:
        multiprocess.current_mem_usage(jobid)

    return results, md_results, mdd_results

#==============================================================================
def calculate_md(md_distances1=None, md_distances2=None, smallwindow=None, 
                    jobid=None, cpus=None, debug=None):
    md_keywords = dict(smallwindow=smallwindow)
    md_results = multiprocess.main(function=md_score, 
                    args=list(zip(sorted(md_distances1),sorted(md_distances2))), 
                    kwargs=md_keywords,
                    debug=debug, jobid=jobid, cpus=cpus)

    md_results = md_score_p(md_results)

    return md_results

#Functions
#==============================================================================
def get_auc_gc(distances, fimo_motifs=None):
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
        if fimo_motifs is not None:
            gc = get_gc(motif=motif, motif_database=fimo_motifs)
        nan = float('Nan')
        distances = distances[1:]
        distances_abs = [abs(x)  if x != '.' else x for x in distances]

        hits = len([x for x in distances_abs if x != '.'])

        #Filter any TFs/files without any hits
        if hits == 0:
            return [motif, nan, gc]

        #Get -exp() of distance and get cumulative scores
        #Filter distances into quartiles to get middle distribution
        q1 = int(round(len(distances)*.25))
        q3 = int(round(len(distances)*.75))
        middledistancehist =  [x for x in distances_abs[int(q1):int(q3)] if x != '.']
        try:
            average_distance = float(sum(middledistancehist))/float(len(middledistancehist))
        except ZeroDivisionError:
            return [motif, nan, gc]
        
        score = [math.exp(-float(x)/average_distance) if x != '.' else 0.0 for x in distances_abs]
        total = sum(score)

        binwidth = 1.0/float(len(distances_abs))
        normalized_score = [(float(x)/total)*binwidth for x in score]
        cumscore = np.cumsum(normalized_score)
        trend = np.append(np.arange(0,1,1.0/float(len(cumscore) - 1)), 1.0)
        trend = [x*binwidth for x in trend]

        #The AUC is the relative to the "random" line
        auc = (np.trapz(cumscore) - np.trapz(trend))*2

    except Exception as e:
        # This prints the type, value, and stack trace of the
        # current exception being handled.
        print(traceback.print_exc())
        raise e
    return [motif, auc, gc]

#==============================================================================
def get_null_data(ranked_file):
    """
    Get the start and the end of the region to consider as the null hypothesis taccording to 
    p-values being above 0.9.
    Input: Ranked File WITH a column called "fc,p-value,rank"
    Output: (pval_09[0], pval_09[-1]): integer positions between which to consider for null data
    """
    # read in the ranked file
    ranked = pandas.read_csv(ranked_file, sep="\t")
    # get the fc and p-values
    split = ranked["fc,p-value,rank"].str.split(",", n=2, expand=True)
    ranked['fc'] = split[0].astype(float)
    ranked['pval'] = split[1].astype(float)
    # get where the pvalues are above 0.9
    pval_09 = np.where(ranked['pval'] > 0.9)[0]
    print(f'\nRegions considered for middle (Num:{len(pval_09)}, Percentage:{round(100*len(pval_09)/ranked.shape[0],2)}, Start:{pval_09[0]}, End:{pval_09[-1]}', 
          flush=True, file=sys.stderr)
    return((pval_09[0], pval_09[-1]))

#==============================================================================
def get_FC1_perc(fcs):
    """
    Get the percentile at which the FC is about 1 (meaning the expected middle region with the most things that 
    are not changing. We do this by finding the first position where the fc is <=1. It will first check 
    if there are any negative values
    Input: fcs (list or numpy array of fold changes (can be LFC or FC))
    Output: index at which the FC is about 1 (things aren't changing at all)
    """
    # get if LFC instead of FC
    fcs = np.array(fcs)
    below0_indices = np.where(fcs <= 0)[0]
    if (len(below0_indices) > 0):
        print("There are negative FCs indicating that this is likely log fold changes. If this is not the case, then there is a problem with the input. Downstream analysis will be performed assuming LFCs.", 
             file=sys.stderr)
        fc1_index = below0_indices[0]
    else:
        fc1_index = np.where(fcs <= 1)[0][0]
    print(f"The position at which FC is 1 or <1 is {fc1_index}", flush=True, file=sys.stderr) 
    return((fc1_index / (len(fcs) - 1)) )

#==============================================================================
def auc_simulate_and_plot(distances, use_config=True, output_type=None, 
                        permutations=None, pvals=None,
                        plotall=None, p_cutoff=None, figuredir=None, 
                        largewindow=None, fimo_motifs=None, 
                        meta_profile_dict=None, label1=None, label2=None, 
                        dpi=None, fcs=None, tests=None, motif_fpkm=None, 
                        bootstrap=False, gc_correct=None, plot_format=None, 
                         ranked_file=None):
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
        nan = float('Nan')
        gc = nan
        if fimo_motifs:
            gc = get_gc(motif=motif, motif_database=fimo_motifs)
        else:
            gc = nan
        try:
            fpkm = motif_fpkm[motif]
        except KeyError:
            fpkm = nan
        distances = distances[1:]
        distances_abs = [abs(x)  if x != '.' else x for x in distances]
        hits = len([x for x in distances_abs if x != '.'])
        #Filter any TFs/files without any hits
        if hits == 0:
            return [motif, 0, 0, hits, gc, fpkm, 0, 0, nan, nan, nan, nan, nan]

        #Get -exp() of distance and get cumulative scores
        # Get the distribution for null (usually just middle)
        # (q1, q3) = get_null_data(ranked_file)
        # print(f"The Q1 and Q3 if using center where pval > 0.8 are {q1} and {q3}", flush=True, file=sys.stdout) 
        # TODO: put so if null_window pushes to limits outside the numbers then fix it
        q1 = int(round(len(distances)*.25))
        q3 = int(round(len(distances)*.75))
        print(f"The Q1 and Q3 if using middle are {q1} and {q3}", flush=True, file=sys.stdout) 
        middledistancehist =  [x for x in distances_abs[int(q1):int(q3)] if x != '.']
        #In the case where there are no hits in the middle two quartiles, then
        #don't perform computation
        if len(middledistancehist) == 0:
            return [motif, 0, 0, hits, gc, fpkm, 0, 0, nan, nan, nan, nan, nan]
        try:
            average_distance = np.mean(middledistancehist)
        except ZeroDivisionError:
            return [motif, 0, 0, hits, gc, fpkm, 0, 0, nan, nan, nan, nan, nan]
        
        score = [math.exp(-float(x)/average_distance) if x != '.' else 0.0 for x in distances_abs]
        total = np.sum(score)

        binwidth = 1.0/float(len(distances_abs))
        normalized_score = np.multiply(np.divide(score, total), binwidth)
        #[(float(x)/total)*binwidth for x in score]
        cumscore = np.cumsum(normalized_score)
        trend = np.append(np.arange(0,1,1.0/float(len(cumscore) - 1)), 1.0)
        trend = np.multiply(trend, binwidth)
        # trend = [x*binwidth for x in trend]

        #The AUC is the relative to the "random" line
        auc = (np.trapz(cumscore) - np.trapz(trend))*2
        offset = 0
        if motif in gc_correct:
            offset = gc_correct[motif]
        corrected_auc = auc - offset

        #Calculate random AUC
        if bootstrap:
            sim_auc = permute_auc_bootstrap(original_distances=distances, trend=trend, 
                                permutations=permutations, bootstrap=bootstrap)
        else:
            sim_auc = permute_auc(distances=normalized_score, trend=trend, 
                                permutations=permutations)

        #Calculate p-value
        mu = np.mean(sim_auc)
        sigma = np.std(sim_auc)
        p = min(stats.norm.logcdf(auc,mu,sigma), stats.norm.logsf(auc,mu,sigma))
        if math.isnan(p):
            p = 0
        p = p+np.log(tests) if p+np.log(tests) <= 0 else 0
        p = p*math.log(np.e, 10)

        corrected_p = min(stats.norm.logcdf(corrected_auc,mu,sigma), stats.norm.logsf(corrected_auc,mu,sigma))
        if math.isnan(corrected_p):
            corrected_p = 0
        corrected_p = corrected_p+np.log(tests) if corrected_p+np.log(tests) <= 0 else 0
        corrected_p = corrected_p*math.log(np.e, 10)

        # GET THE LEADING EDGE
        if plotall or (p < p_cutoff) or (corrected_p < p_cutoff):
            print("Getting Leading Edge for", motif)
            le_mb, le_se, le_mb_stdev, le_se_stdev, frac_back = get_lead_edges(cumscores=cumscore, Enr_score=auc, num_motif_regions=hits)
        else:
            le_mb = nan
            le_mb_stdev = nan
            le_se = nan
            le_se_stdev = nan
            frac_back = nan

        if plotall or (output_type=='html' and p < p_cutoff) or (output_type=='html' and corrected_p < p_cutoff):
            from TFEA import plot
            plotting_score = np.divide(score, total)
            # [(float(x)/total) for x in score]
            plotting_cumscore = np.cumsum(plotting_score)
            plot.plot_individual_graphs(motif=motif, distances=distances, 
                                        figuredir=figuredir, 
                                        fimo_motifs=fimo_motifs, 
                                        largewindow=largewindow, 
                                        score=plotting_score, 
                                        use_config=use_config, 
                                        pvals=pvals, fcs=fcs, 
                                        cumscore=plotting_cumscore, 
                                        sim_auc=sim_auc, auc=auc,
                                        meta_profile_dict=meta_profile_dict, 
                                        label1=label1, label2=label2, 
                                        offset=offset, plot_format=plot_format, 
                                        le_mb=le_mb, le_se=le_se)
        
        # once plotted fix the le_mb and le_se to fit the opposite direction if a "repressor"
        if auc < 0:
            le_mb = len(distances) - le_mb
            le_se = len(distances) - le_se

    except Exception as e:
        # This prints the type, value, and stack trace of the
        # current exception being handled.
        print(traceback.print_exc())
        raise e
    ## SAVE THE TEMPORARAY FILES FOR LEADING EDGE ##
    list_of_tuples = list(zip(distances, distances_abs))
    df = pandas.DataFrame(list_of_tuples, columns = ['distances', 'distances_abs'])
    df["score"]=score
    df["normalized_score"]=normalized_score
    df["cumscore"]=cumscore
    tempdir = config.vars['TEMPDIR']
    df.to_csv(str(tempdir)+"/withAUC"+str(motif)+"__dis_andmore.csv")
    print("MOTIF:" + str(motif) + " AUC:" + str(auc) + " " + str(corrected_auc) + 
    " Events:" + str(hits) + " Pval: " + str(p) + " " + str(corrected_p) + 
    "\n\tLEMB: " + str(le_mb) + " Stdev: " + str(le_mb_stdev) + "\tLESE: " + str(le_se) + "Stdev: " + str(le_se_stdev) + "\tFrac: " + str(frac_back), file=sys.stderr)
    return [motif, auc, corrected_auc, hits, gc, fpkm, p, corrected_p, le_mb, le_mb_stdev, le_se, le_se_stdev, frac_back]

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
    es_permute = np.zeros(permutations)
    triangle_area = np.trapz(trend)
    for i in range(permutations):
        random_distances = np.random.permutation(distances)
        cum_distances = np.cumsum(random_distances)
        es = np.trapz(cum_distances)
        auc = es - triangle_area
        es_permute[i] = auc*2

    return es_permute

#==============================================================================
def permute_auc_bootstrap(original_distances=None, trend=None, permutations=None, bootstrap=False):
    '''Generates permutations of the original_ and calculates AUC for each 
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
    hits = [x for x in original_distances if x != '.']
    if bootstrap and bootstrap <= len(original_distances):
        if bootstrap < len(hits):
            #Determine where there are hits
            hit_indexes = [i for i,x in enumerate(original_distances) if x != '.']

            #Subsample down hits to bootstrap number
            new_hit_indexes = np.random.choice(hit_indexes, bootstrap, replace=False)

            #Only keep hit if it is within subsampled indexes, otherwise replace it with '.'
            new_distances = [x if i in new_hit_indexes else '.' for i,x in enumerate(original_distances)]
        else:
            #Get all numerical hit values
            hits = [x for x in original_distances if x != '.']

            #Generate a list of hit values of size bootstrap
            new_hits = iter(np.random.choice(hits, bootstrap, replace=True))

            #Determine where the hits will be in the new distances array
            new_hit_indexes = np.random.choice(len(original_distances), bootstrap, replace=False)

            #Place the hits in their appropriate place
            new_distances = [next(new_hits) if i in new_hit_indexes else '.' for i in range(len(original_distances))]

    #Get -exp() of distance and get cumulative scores
    #Filter distances into quartiles to get middle distribution
    distances_abs = [abs(x) if x != '.' else x for x in new_distances]
    q1 = int(round(len(new_distances)*.25))
    q3 = int(round(len(new_distances)*.75))
    middledistancehist =  [x for x in distances_abs[int(q1):int(q3)] if x != '.']
    average_distance = float(sum(middledistancehist))/float(len(middledistancehist))
    
    score = [math.exp(-float(x)/average_distance) if x != '.' else 0.0 for x in distances_abs]
    total = sum(score)

    binwidth = 1.0/float(len(distances_abs))
    normalized_score = [(float(x)/total)*binwidth for x in score]

    es_permute = []
    triangle_area = np.trapz(trend)

    normalized_distances = normalized_score
    for i in range(permutations):
        random_distances = np.random.permutation(normalized_distances)
        cum_distances = np.cumsum(random_distances)
        es = np.trapz(cum_distances)
        auc = es - triangle_area
        es_permute.append(auc*2)

    return es_permute

#==============================================================================
def max_GSEA(distances, use_config=False, output_type=None, cutoff=None,
                        permutations=None, pvals=None,
                        plotall=None, p_cutoff=None, figuredir=None, 
                        largewindow=None, fimo_motifs=None, 
                        meta_profile_dict=None, label1=None, label2=None, 
                        dpi=None, fcs=None):
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
        p = stats.norm.cdf(actualES, mu, sigma)
    else:
        simESsubset = [x for x in simES if x > 0]
        mu = np.mean(simESsubset)
        NES = actualES/mu
        sigma = np.std(simESsubset)
        p = 1-stats.norm.cdf(actualES, mu, sigma)

    if math.isnan(p):
        p = 1.0

    # if output_type == 'html':
    #     from TFEA import plot
    #     plotting_score = [(float(x)/total) for x in score]
    #     plotting_cumscore = np.cumsum(plotting_score)
    #     plot.plot_individual_graphs(distances=distances, figuredir=figuredir, 
    #                     fimo_motifs=fimo_motifs, largewindow=largewindow, 
    #                     score=plotting_score, use_config=use_config,
    #                     dpi=dpi, pvals=pvals, fcs=fcs, 
    #                     cumscore=plotting_cumscore, sim_auc=sim_auc, auc=auc,
    #                     meta_profile_dict=meta_profile_dict, label1=label1, 
    #                     label2=label2)
    
    if output_type == 'lineplot_output':
        return ES
    
    if output_type == 'simulation_output':
        return actualES, simES

    return [motif, actualES, NES, pos_total, p]

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

    stat, _, sig = stats.anderson_ksamp([observed_distribution, expected_distribution])

    return [motif, stat, hits, sig]

#==============================================================================
def md_score(distances, smallwindow=None):
    '''Calculate md score
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
    '''Calculate a p-value for md-score results
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
            cdf=stats.norm.cdf(z,0,1)
            p=min(cdf,1-cdf)*2
            md = md2-md1-mean
            total = (N1 + N2)/2.0
        
        results[i] = [motif, md, total, p]

    return results

#==============================================================================
def get_gc(motif=None, motif_database=None, alphabet=['A','C','G','T']):
    '''
    Obtain a pssm model from a meme formatted database file. Warning: If there are multiple
    motif matches, this function will return the last match in the database.
    
    Parameters
    ----------
    motif_database : string
        full path to a MEME formatted (.meme) file
    motif : string
        the name of a motif that exactly matches a motif in motif_database
        
    Returns
    -------
    PSSM : list or array
        a list of lists where each corresponding to position then alphabet probability
    '''
    motif_hit = False
    PSSM = list()
    with open(motif_database,'r') as F:
        for line in F:
            if 'MOTIF' in line:
                if motif in line:
                    motif_hit = True
                else:
                    motif_hit = False
            elif motif_hit and 'URL' not in line and 'letter-probability' not in line and line != '\n':
                acgt_probabilities = [float(x) for x in line.strip('\n').split()]
                total_prob = sum(acgt_probabilities)
                acgt_probabilities = [x/total_prob for x in acgt_probabilities] #Convert to probabilities
                PSSM.append(acgt_probabilities)

    gc = 0
    for base in PSSM:
        gc += base[alphabet.index('C')]
        gc += base[alphabet.index('G')]
    
    gc = gc/float(len(PSSM))
    
    return gc
#==============================================================================
###################### LEADING EDGE FUNCTIONS ########################
#==============================================================================
###########################
###    Main Function    ###
###########################
def get_lead_edges(cumscores, Enr_score, num_motif_regions, increment_value=0.5):
    """
    Gets the median Leading Edge across multiple splines using both a Matchback and
    Stalled Enrichment leading edge method.

    Parameters:
    * cumscores: cumulative enrichment scores from enrichment curve
    * Enr_score: final AUC of TF
    * num_motif_regions: number of regions with motif hits
    * increment_value: positive float of 0.5 or 1 (default 0.5): number of base splines to increment
         by (e.g. 1 means go from 9e-11 to 8e-11, .5 means go from 9e-11 to 8.5e-11)

    Output:
    * le_mb: (int) Final median LE according to the place where slope matches null background expectation
    * le_se: (int) Final median LE according to stalled enrichment across splines
    * le_mb_stdev: (float, rounded to 3 decimals) Standard deviation of le_mb across splines (If only one spline used, NaN)
    * le_se_stdev: (float, rounded to 3 decimals) Standard deviation of le_se across splines (If only one spline used, NaN)
    * final_frac_back: (float, rounded to 2 decimals) Median fraction of regions with enrichment curve
            slope above background (meant to help pinpoint FPs)
    """
    num_regions = len(cumscores)
    # print("Num regions:", num_regions, "Num Motif Regions:", num_motif_regions, 
    #  "\tPer", round(num_motif_regions/num_regions, 2)*100, "%")
    # Get initial spline and incrementation
    num_der_limit = get_num_der_limit(num_regions)
    spline_use, spline_power = get_start_spline(num_motif_regions)
    # ITERATE THROUGH THE SPLINES
    # Lists ot hold things
    spline_list = []
    #der_list = []; second_der_list = []
    LE_mb_list = []; LE_se_list = []; LE_se_type_list = []
    frac_back_list = []
    # Default starting values
    num_der = 0
    ranks = np.arange(0, len(cumscores), 1)
    background_slope=cumscores[-1]/num_regions
    quint_point = int(num_regions/5)
    fourquint_point = num_regions - quint_point
    # keep track of number of times tried 
    trials = 0
    # If there are more than 14 max/min in second derivative, stop to avoid noise
    while (num_der < num_der_limit) and (trials < 20):
        # get the smoothing information for the enrichment curve
        tck = splrep(ranks, cumscores, k=5, s=spline_use*10**spline_power)
        # evaluate the spline's first, second, and third derivative
        y_spline_prime = splev(ranks, tck, der=1)
        y_spline_2ndder = splev(ranks, tck, der=2)
        y_spline_3rdder = splev(ranks, tck, der=3)
        # determine if include this for LE prediction based on # curves in 2nd der near LE
        if Enr_score > 0:
            num_der = count_min_max(y_spline_2ndder[0:quint_point])
        else:
            num_der = count_min_max(y_spline_2ndder[fourquint_point:len(y_spline_2ndder)])
        # If reasonable amount of noise
        #print("NUM DER", num_der, "SPLINE", spline_use*10**spline_power)
        trials = trials + 1
        if num_der < num_der_limit:
            ## GET THE MIN POSITION WHERE SLOPE MATCHES BACKGROUND
            stand_list = y_spline_prime - background_slope
            matchback, frac_back = get_matchback(stand_list, Enr_score)
            ## GET THE 2nd DER BASED LE (Stalled_enrichment)
            curve_pos, peak_type = get_LE_stalled_enrichment(num_regions, y_spline_2ndder, y_spline_3rdder, Enr_score)
            # If the 2nd derivative is at a decent position
            if ((curve_pos < quint_point) & (Enr_score > 0)) or ((curve_pos > fourquint_point) & (Enr_score < 0)):
                # print("SPLINE Adding", spline_use, spline_power, Enr_score)
                # # get graph related lists
                # spline_list.append(spline_use*10**spline_power)
                # der_list.append(y_spline_prime)
                # second_der_list.append(y_spline_2ndder)
                # Save LEs
                LE_mb_list.append(matchback)
                LE_se_list.append(curve_pos)
                LE_se_type_list.append(peak_type)
                # Save fraction below background
                frac_back_list.append(frac_back)
        elif len(LE_mb_list) == 0:
            og_spline = spline_use*10**spline_power
            # if couldn't get any , move back the spline and change the num_der
            spline_use = spline_use + 3
            if spline_use > 9:
                spline_use = spline_use - 9 # so 10 becomes 1 and 11 2
                spline_power = spline_power + 1
            print("HAD TO START OVER with new spline of ", spline_use*10**spline_power, " from ", og_spline, num_regions)
            num_der=0
        ## ASSESSING THE NEXT SPLINE
        spline_use, spline_power = get_next_spline_use(spline_use, increment_value, spline_power)
    # If nothing then just say NA
    if len(LE_mb_list) == 0:
        print("REACHED GREATER THAN 20 TRIALS WITH NO LUCK")
        le_mb = np.nan
        le_se = np.nan
        le_mb_stdev = np.nan
        le_se_stdev = np.nan
        final_frac_back = np.nan
    elif len(LE_mb_list) == 1:
        le_mb_stdev = np.nan
        le_se_stdev = np.nan
        le_mb = int(np.median(LE_mb_list))
        le_se = int(np.median(LE_se_list))
        final_frac_back = round(np.median(frac_back_list), 2)
    else:
        le_mb = int(np.median(LE_mb_list))
        le_se = int(np.median(LE_se_list))
        le_mb_stdev = round(np.std(LE_mb_list),2)
        le_se_stdev = round(np.std(LE_se_list),2)
        final_frac_back = round(np.median(frac_back_list), 2)
    
    
    return le_mb, le_se, le_mb_stdev, le_se_stdev, final_frac_back

########################
##  GETTING SPLINES  ###
########################
def get_next_spline_use(spline_use, increment_value, spline_power):
    """
    Gets the next spline smoothing parameter to use given the current
    spline, power, and increment value. Can handle shifting exponents
    (e.g. going from 1e-10 to 9e-11 with increment of 1).

    Parameters:
    * spline_use: current base value of spline (not including exponent)
    * increment_value: the value by which we change the spline_use.
    * spline_power: The current exponent used 

    Outputs:
    * spline_use: new base value of spline to use
    * spline_power: new exponent of spline to use
    """
    spline_use = spline_use - increment_value
    # if got to 0, make it the next increment
    if spline_use < 1:
        spline_use = 10 - increment_value
        spline_power = spline_power - 1
    # if changing powers
    if spline_use < 0:
        spline_use = 10 + spline_use # so if 5e-10 to 9e-11, 5-6=-1 10+-1=9 & power-1
        spline_power = spline_power - 1
    return spline_use, spline_power

# print("Should be 9 and -11", get_next_spline_use(1, 1, -10))
# print("Should be 9.5 and -11", get_next_spline_use(1, .5, -10))

def get_num_der_limit(num_regions):
    """
    Gets the limit to the number of curves that should be within the first part
    of the enrichment curve's second derivative (according to 
    if AUC is + or -)
    """
    if num_regions > 60000:
        num_der_limit = 6
    else:
        num_der_limit = 4
    return num_der_limit

def get_start_spline(num_motif_regions):
    """
    Gets the starting spline smoothing parameter based on the number
    of regions with a motif in them.
    Parameters
    * num_motif_regions: (int) number of regions with a motif in them
    Outputs
    * start_spline: (int) base of exponent for smoothing parameter (e.g. 1 for 1e-X)
    * spline_power: (int) exponent for smoothing parameter (e.g. -11 for Xe-11)
    """
    if num_motif_regions > 30000:
        start_spline = 3
        spline_power = -12
    elif num_motif_regions > 20000:
        start_spline = 4
        spline_power = -12
    elif num_motif_regions > 10000:
        start_spline = 5
        spline_power = -12
    elif num_motif_regions > 8000:
        start_spline = 2
        spline_power = -11
    elif num_motif_regions > 7000:
        start_spline = 3
        spline_power = -11
    elif num_motif_regions > 5000:
        start_spline = 4
        spline_power = -11
    elif num_motif_regions > 4000:
        start_spline = 5
        spline_power = -11
    elif num_motif_regions > 3000:
        start_spline = 7
        spline_power = -11
    elif num_motif_regions > 2500:
        start_spline = 8
        spline_power = -11
    elif num_motif_regions > 2000:
        start_spline = 9
        spline_power = -11
    elif num_motif_regions > 1500:
        start_spline = 4
        spline_power = -11
    elif num_motif_regions > 1000:
        start_spline = 1
        spline_power = -10
    else:
        start_spline = 2
        spline_power = -10
    # elif num_motif_regions > 600: # help
    #     start_spline = 3
    #     spline_power = -10
    # elif num_motif_regions > 500: # help
    #     start_spline = 6
    #     spline_power = -10
    # elif num_motif_regions > 200:
    #     start_spline = 7
    #     spline_power = -10
    # elif num_motif_regions > 100: # help
    #     start_spline = 8
    #     spline_power = -10
    return start_spline, spline_power

########################################
###    Supporter Functions For LE    ###
########################################
def find_curve_point(lst, rev=False):
    """
    Get the earliest minima and maxima positions. If none then returns None
    Used for get_LE_stalled_enrichment to find where the first (in the case
    of an activator) or last (in case of a repressor) local minimum and
    maximum are.

    Parameters:
    lst: list of numbers (usually a derivative) 
    rev: If True then reverses the list 

    Outputs:
    {"min":min_index, "max":max_index}:
        * min_index is the first index at which values change from negative to positive
            (In the case of a derivative -- local minimum)
        * max_index is the first index at which values change from positive to negative
            (In the case of a derivative -- local max)
    """
    # Convert the list to a NumPy array
    if rev:
        lst = lst[::-1]
    arr = np.array(lst)
    min_index=None
    max_index=None
    # Ensure the array has at least 2 elements
    if len(arr) < 2:
        return {"min":None, "max":None}
    # MAX: Check for transitions from positive to negative (up then down)
    for i in range(len(arr) - 1):
        if i==0:
            i=1
        if arr[i] > 0 and arr[i + 1] < 0:
            max_index =  i   # Return the index right before the positive number
            break
        elif arr[i] > 0 and arr[i + 1] == 0 and arr[i + 2] < 0:
            max_index = i + 1 # return the 0 positon
            break
    # MIN: Check for transitions from negative to positive (down then up)
    for i in range(len(arr) - 1):
        if arr[i] < 0 and arr[i + 1] > 0:
            min_index = i   # Return the index right before the positive number
            break
        elif arr[i] < 0 and arr[i + 1] == 0 and arr[i + 2] > 0:
            min_index = i + 1 # return the 0 positon
            break
    if rev:
        if min_index is not None:
            min_index = len(arr) - min_index - 1
        if max_index is not None:
            max_index = len(arr) - max_index - 1
    return {"min":min_index, "max":max_index}

# # Example usage
# my_list = [-3, -2, 0, 4, -1, 2]
# print("Should be 2, 3", find_curve_point(my_list)) 
# print("Should be 4, 2", find_curve_point(my_list, rev=True)) 
# my_list = [-3, -2, 0, 4, -1, 2, -3, -4]
# print("Should be 5, 6", find_curve_point(my_list, rev=True)) 

# my_list = [-1, -0.000002, 2, 3, 0, -1]
# print("Should be 1, 4", find_curve_point(my_list))
# print("Should be 4, 2", find_curve_point(my_list, rev=True)) 

def get_curve_elbow(num_regions, y_spline_secder, Enr_score=0):
    """
    Gets the first (activator) or last (repressor) index where the max distance 
    from the background line is.
    If an activator (Enr_score > 0) then looks in the first "half" (len/2.5).
    If a repressor (Enr_score < 0) then looks in the second "half." If a 
    Used for get_LE_stalled_enrichment.
    *Note*: get_dist_point_line is the main time bottleneck.

    Parameters:
    * num_regions: Number of ranked regions being assessed
    * y_spline_secder: Spline of second derivative of the enrichment curve
    * Enr_score: AUC of TF to determine if activator or repressor
    """
    half_point = int(num_regions/2.5)
    tolerance = 1e-15
    # Get set up for point-line distance computations
    first_point = np.array([1, y_spline_secder[0]])
    end_point = np.array([num_regions, y_spline_secder[num_regions - 1]])
    # get the distances of the second derivative from the null line --> find elbow
    mfg_dists = [
        get_dist_point_line(np.array([k+1, y_spline_secder[k]]), first_point, end_point)
        for k in range(num_regions)
    ]
    if Enr_score > 0:
        # Get the indices that have the maximum distances
        max_dist_indices = np.where(np.isclose(mfg_dists[0:half_point], max(mfg_dists[0:half_point]), atol=tolerance))[0]
        return(np.min(max_dist_indices))
    else:
        # "half point" is the opposite direction
        half_point = num_regions - half_point
        # Get the indices that have the maximum distances
        max_dist_indices = np.where(np.isclose(mfg_dists[half_point:len(mfg_dists)], 
                                                  max(mfg_dists[half_point:len(mfg_dists)]), atol=tolerance))[0] 
        #only consider indices 
        return(np.max(max_dist_indices)+half_point)


# Function to calculate the distance from a point to a line
def get_dist_point_line(point, line_coord1, line_coord2):
    """
    Get the distance from a point to a line with two coordinates: line_coord1 and line_coord2
    """
    line_vec = line_coord2 - line_coord1
    point_vec = point - line_coord1
    line_len = np.linalg.norm(line_vec)
    line_unitvec = line_vec / line_len
    proj_length = np.dot(point_vec, line_unitvec)
    proj_point = line_coord1 + proj_length * line_unitvec
    distance = np.linalg.norm(point - proj_point)
    return distance

# print("Should be 0.1403663161257098", get_dist_point_line(np.array([2, 0.3]), np.array([1, 0]), np.array([20, 3])))
# print("Should be 3.3593115504912983", get_dist_point_line(np.array([3, 4]), np.array([1, 0]), np.array([20, 5])))
# print("Should be 4.377284747609873", get_dist_point_line(np.array([3, -4]), np.array([1, 0]), np.array([20, 5])))


def count_min_max(arr):
    """
    Count the number of curves (min/max) based on the 
    number of sign changes in the derivative (arr)
    """
    signs = np.sign(arr)  # Get the sign of each element (-1, 0, or 1)
    sign_changes = np.diff(signs) != 0  # Check where the sign changes
    return int(np.count_nonzero(sign_changes))  # Count the number of changes

# # Example usage
# arr = np.array([-3, -2, 0, 2, -1, 5, -6, 7])
# print("Should be:", count_min_max(arr))
# arr = np.array([7, -1, 2, 4, 1, -2, 3, 9, 10, -3, 0, 4, 2, -4, -5, 2, 3, -6])
# print("Number of min or max:", count_min_max(arr))

#################################
###     Individual Methods    ###
#################################
def get_LE_stalled_enrichment(num_regions, y_spline_2ndder, y_spline_3rdder, Enr_score=0):
    """
    This function returns the median LE across the spline methods where the individual
       LE for each spline is considered accordingly:
            1. "MAX": If the 2nd derivative starts out negative, it will return the first max after which the 2nd derivative
                either stabilizes or becomes more negative.
            2. "MIN": If the 2nd derivative starts at a positive value, meaning the enrichment scores are increasing, 
               it will return the first minimum after which the the 2nd derivative starts stabilizing.
            3. "ELBOW": If there are neither maxima/minima OR the above case leads to a LE call 
                after the halfway point (hence indicating a noise is captured), the elbow of the
                2nd derivative is used (Details in get_curve_elbow)
    
    Inputs:
    * num_regions: number of regions under study
    * y_spline_2ndder: spline based y values for the 2nd derivative
    * y_spline_3rdder: spline based y values for the 3rd derivative
    * Enr_score: The enrichment score (to determine if positive or negative)

    Outputs:
    * curve_pos: The final LE according to the considerations above
    * peak_type: The methodology used to get the final LE according to the considerations above
    """
    half_point = int(num_regions/2.5)
    if Enr_score > 0:
        # Get the first min and max of the 2nd derivative from 3rd derivative
        curve_indices_dict = find_curve_point(y_spline_3rdder)
        # If the 2nd derivative starts out positive (use position 2 for noise), get the min
        if y_spline_2ndder[2] > 0:
            # If no clear minimum or the minimum is > half point, use elbow
            if curve_indices_dict["min"] is None or curve_indices_dict["min"] > half_point:
                curve_pos = get_curve_elbow(num_regions, y_spline_2ndder, Enr_score)
                peak_type = "ELBOW"
            else:
                curve_pos = curve_indices_dict["min"]
                peak_type = "MIN"
        else:
            # If starts out negative, get the first max or elbow
            curve_pos = get_curve_elbow(num_regions, y_spline_2ndder, Enr_score)
            if curve_indices_dict["max"] is None or curve_indices_dict["max"] > curve_pos:
                peak_type = "ELBOW"
            else:
                curve_pos = curve_indices_dict["max"]
                peak_type = "MAX"
    else:
        # "half point" is in opposite direction
        half_point = num_regions - half_point
        # Get the "first" min and max of the 2nd derivative from 3rd derivative (reversed)
        curve_indices_dict = find_curve_point(y_spline_3rdder, rev=True)
        # If the 2nd derivative "starts negative" then just get the maximum
        if y_spline_2ndder[-2] < 0:
            # If no clear minimum or the minimum is > half point (so <), use elbow
            if curve_indices_dict["min"] is None or curve_indices_dict["min"] < half_point:
                curve_pos = get_curve_elbow(num_regions, y_spline_2ndder)
                peak_type = "ELBOW"
            else:
                curve_pos = curve_indices_dict["min"]
                peak_type = "MIN"
        else:
            # If starts out negative, get the first max or elbow
            curve_pos = get_curve_elbow(num_regions, y_spline_2ndder)
            # print("CURVE INDICES after negative", curve_indices_dict)
            # print("CURVE POS", curve_pos)
            if curve_indices_dict["max"] is None or curve_indices_dict["max"] < curve_pos:
                peak_type = "ELBOW"
            else:
                curve_pos = curve_indices_dict["max"]
                peak_type = "MAX"
        
    # RETURN THE POSITION AND TYPE
    return(curve_pos, peak_type)


def get_matchback(stand_list, Enr_score=0):
    """
    Gets the first instance where the cumulative enrichment change lowers to that
    expected from background with t.
     Enr_score > 0: ACT: want min index where starts higher then goes below background (so from positive to negative) 
     Enr_score < 0: REP: want max index where starts lower than goes above background (so from neg to positive)

     Outputs:
     * matchback:
     * frac_background: Fraction of tREs at which slope is > than that of background
    """
    # change to a numpy array
    arr = np.array(stand_list)
    frac_background = np.sum(arr > 0)/len(stand_list)
    if Enr_score > 0:
        # Thinking of activator - where goes from higher to lower
        transitions = np.where((arr[:-1] > 0) & (arr[1:] <= 0))[0]
        if transitions.size == 0:
            matchback = 0
        else:
            matchback = np.min(transitions) + 1
    else:
        # Thinking of activator - where goes from lower to higher
        transitions = np.where((arr[:-1] < 0) & (arr[1:] >= 0))[0]
        if transitions.size == 0:
            matchback = 0
        else:
            matchback = np.max(transitions) + 1
    return matchback, frac_background
        

# stand_list = [1, 3, 5, 0, -2, -5, 2, 4, 5, -2]
# print("Should get 3 and .6", get_matchback(stand_list, Enr_score=.1))
# print("Should get 6 and .3", get_matchback(stand_list, Enr_score=-.1))

# stand_list = [-2, 3, 4, 5, -3, -4, 1, 2, 3, 4]
# print("Should get 4 and .7", get_matchback(stand_list, Enr_score=.1))
# print("Should get 6 and .3", get_matchback(stand_list, Enr_score=-.1))
# # ACT: want min index where starts higher then goes below background (so from positive to negative) 
# # REP: want max index where starts lower than goes above background (so from neg to positive)




