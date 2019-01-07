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
def txt_output(results=None, outputdir=None, enrichment=None):
    with open(os.path.join(outputdir, 'results.txt'), 'w') as outfile:
        outfile.write('TF Motif\tEnrichment(' + enrichment + ')\tMotif Hits\tp-value\tp-adj\n')
        for line in sorted(results key=lambda x: x[-2]): #Sort by p-value
            outfile.write('\t'.join(line)+'\n')

#==============================================================================
def html_output(results=None, module_list=None):
    '''Creates the main html output and also individual html outputs for each
        motif
    
    Parameters
    ----------
    outputdir : string
        the full path to the output directory created by TFEA

    beds : list or array
        a list of full paths to bed files to be considered as regions of 
        interest

    label1 : string
        an informative label describing sample corresponding to condition1

    label2 : string
        an informative label describing sample corresponding to condition2

    bam1 : list or array
        a list of full paths to bam files corresponding to condition1

    bam2 : list or array
        a list of full paths to bam files corresponding to condition2

    singlemotif : boolean or string
        either False if all motifs should be considered in TFEA or the name of
        a specific motif to be analyzed

    motif_hits : string
        the full path to a directory containing motif hits across the genome

    output : string
        the full path to a user-specified output directory. TFEA will create
        a new folder within this directory - this is called outputdir

    padj_cutoff : float
        the cutoff value for determining significance

    plot : boolean
        a switch that controls whether all motifs are plotted or just 
        significant ones defined by the p-adj cutoff

    combine : boolean
        a switch that determines whether bed files within the beds variable
        get combined and merged using bedtools

    count : boolean
        a switch that controls whether reads are counted over the regions of
        interest

    deseq : boolean
        a switch that controls whether DE-Seq is performed on the inputted
        regions that have been counted over

    calculate : boolean
        a switch that determines whether the TFEA calculation is performed

    TFresults : list or array
        a list of lists contining 'enrichment' scores, normalized 'enrichment' 
        scores, p-value, p-adj, and number of hits for each individual motif

    COMBINEtime : float
        the time it took to combine and merge the bed files using bedtools

    COUNTtime : float
        the time it took to count reads over regions of interest

    DESEQtime : float
        the time it took to perform DE-Seq using the counts file

    CALCULATEtime : float
        the time it took to perform TFEA

    Returns
    -------
    None
    '''
    #Using a config file
    import config
    outputdir = config.OUTPUTDIR
    label1 = config.LABEL1
    label2 = config.LABEL2
    padj_cutoff = config.PADJCUTOFF
    plot = config.PLOTALL

    outfile = open(os.path.join(outputdir, 'results.html'),'w')
    outfile.write("""<!DOCTYPE html>
    <html>
    <head>
    <title>TFEA Results """ + label1 + """ vs. """ + label2 +"""</title>
    <style>
        table {
            font-family: arial, sans-serif;
            border-collapse: collapse;
            width: 100%;
        }

        .row {
        display: flex; /* equal height of the children */
        width: 100%;
        padding-bottom: 50px
        }

        img {
            max-width: 100%;
            max-height: 100%;
        }

        td, th {
            border: 1px solid #dddddd;
            text-align: left;
            padding: 8px;
        }

        tr:nth-child(even) {
            background-color: #dddddd;
        }
    </style>
    </head>
    <body style="width:1300px; margin:0 auto;">

        <h1>TFEA Results """ + label1 + """ vs. """ + label2 + """</h1>
        <div class="row">
            <div style="float: left; width: 45%">
                <img src="./plots/TFEA_NES_MA_Plot.png" alt="NES MA-Plot">
            </div>
            <div style="float: left; width: 45%">
                <img src="./plots/DESEQ_MA_Plot.png" alt="DE-Seq MA-Plot">
            </div>
        </div>
        <div class="row">
            <div style="float: left; width:45%">
                <img src="./plots/TFEA_Pval_Histogram.png" alt="TFEA P-Value \
                    Histogram">
            </div>
            <div id="Summary of Variables Used" style="float: right; \
                width: 45%">
                <p><a href="./Summary.html">Full Summary of Variables Used\
                </a></p>
                <p><b>PADJ < """ + str(padj_cutoff) + """</b></p>
                <table>
                    <tr>
                        <th>Module</th>
                        <th>Value</th>
                        <th>Time (hh:mm:ss)</th>
                    </tr>""")
    for (module, value, time) in module_list:
        outfile.write("""<tr>
                        <td>"""+module+"""</td>
                        <td>"""+str(value)+"""</td>
                        <td>"""+str(datetime.timedelta(
                                    seconds=int(time)))
                        +"""</td>
                    </tr>""")
    outfile.write("""
                </table>   
            </div>
        </div>
        <div>
            <div id="Positive Enrichment Score" style="float: left; width:45%">
                <h1>Positive Enrichment Score</h1>
                <table> 
                    <tr>
                        <th>TF Motif</th>
                        <th>ENRICHMENT("""+config.ENRICHMENT+""")</th>
                        <th>P-value</th>
                        <th>PADJ</th>
                        <th>HITS</th>
                    </tr>
                """)

    for [motif, enrichment, hits, p-value, p-adj] in results:
        if enrichment > 0:
            if p-adj < padj_cutoff:
                outfile.write("""
            <tr style="color: red;">
                <td><a href="./plots/"""+motif+""".results.html">"""
                    +motif+"""</td>
                <td>"""+str("%.4f" % enrichment)+"""</td>
                <td>"""+str("%.3g" % p-value)+"""</td>
                <td>"""+str("%.3g" % p-adj)+"""</td>
                <td>"""+str(hits)+"""</td>
            </tr>
                    """)
            elif plot:
                outfile.write("""
            <tr>
                <td><a href="./plots/"""+motif+""".results.html">"""
                    +motif+"""</td>
                <td>"""+str("%.4f" % enrichment)+"""</td>
                <td>"""+str("%.3g" % p-value)+"""</td>
                <td>"""+str("%.3g" % p-adj)+"""</td>
                <td>"""+str(hits)+"""</td>
            </tr>
                    """)

            else:
                outfile.write("""
            <tr>
                <td>"""+motif+"""</td>
                <td>"""+str("%.4f" % enrichment)+"""</td>
                <td>"""+str("%.3g" % p-value)+"""</td>
                <td>"""+str("%.3g" % p-adj)+"""</td>
                <td>"""+str(hits)+"""</td>
            </tr>
                    """)


    outfile.write("""            
        </table>
    </div>

    <div id="Negative Enrichment Score" style="float: right; width: 45%">
        <h1>Negative Enrichment Score</h1>
        <table> 
            <tr>
                <th>TF Motif</th>
                <th>ENRICHMENT("""+config.ENRICHMENT+""")</th>
                <th>P-value</th>
                <th>PADJ</th>
                <th>HITS</th>
            </tr>
                """)

    for [motif, enrichment, hits, p-value, p-adj] in results:
        if enrichment < 0:
            if p-adj < padj_cutoff:
                outfile.write("""
            <tr style="color: red;">
                <td><a href="./plots/"""+motif+""".results.html">"""
                    +motif+"""</td>
                <td>"""+str("%.4f" % enrichment)+"""</td>
                <td>"""+str("%.3g" % p-value)+"""</td>
                <td>"""+str("%.3g" % p-adj)+"""</td>
                <td>"""+str(hits)+"""</td>
            </tr>
                    """)
            elif plot:
                outfile.write("""
            <tr>
                <td><a href="./plots/"""+motif+""".results.html">"""
                    +motif+"""</td>
                <td>"""+str("%.4f" % enrichment)+"""</td>
                <td>"""+str("%.3g" % p-value)+"""</td>
                <td>"""+str("%.3g" % p-adj)+"""</td>
                <td>"""+str(hits)+"""</td>
            </tr>
                    """)
            else:
                outfile.write("""
            <tr>
                <td>"""+motif+"""</td>
                <td>"""+str("%.4f" % enrichment)+"""</td>
                <td>"""+str("%.3g" % p-value)+"""</td>
                <td>"""+str("%.3g" % p-adj)+"""</td>
                <td>"""+str(hits)+"""</td>
            </tr>
                    """)

    outfile.write("""        
            </table>
        </div>
        </div>

    </body>
    </html>""")

    outfile.close()

#==============================================================================
def summary_html_output(config_object=None, outputdir=None):
    with open(os.path.join(outputdir,'summary.html'),'w') as outfile:
        outfile.write("""<!DOCTYPE html>
                <html>
                <head>
                <title>Variables Used</title>
                </head>
                <body>
                    <h1>Variables Used</h1>""")
        for key in config_object:
            for item in config_object[key]:
                outfile.write('<p>' + item.upper() + '=' 
                                + config_object[key][item] + '</p>\n')
        outfile.write('</body>')

#==============================================================================
def create_motif_result_html(results=None):
    '''Creates the main html output and also individual html outputs for each
        motif
    
    Parameters
    ----------
    outputdir : string
        the full path to the output directory created by TFEA

    beds : list or array
        a list of full paths to bed files to be considered as regions of 
        interest

    label1 : string
        an informative label describing sample corresponding to condition1

    label2 : string
        an informative label describing sample corresponding to condition2

    bam1 : list or array
        a list of full paths to bam files corresponding to condition1

    bam2 : list or array
        a list of full paths to bam files corresponding to condition2

    singlemotif : boolean or string
        either False if all motifs should be considered in TFEA or the name of
        a specific motif to be analyzed

    motif_hits : string
        the full path to a directory containing motif hits across the genome

    output : string
        the full path to a user-specified output directory. TFEA will create
        a new folder within this directory - this is called outputdir

    padj_cutoff : float
        the cutoff value for determining significance

    plot : boolean
        a switch that controls whether all motifs are plotted or just 
        significant ones defined by the p-adj cutoff

    combine : boolean
        a switch that determines whether bed files within the beds variable
        get combined and merged using bedtools

    count : boolean
        a switch that controls whether reads are counted over the regions of
        interest

    deseq : boolean
        a switch that controls whether DE-Seq is performed on the inputted
        regions that have been counted over

    calculate : boolean
        a switch that determines whether the TFEA calculation is performed

    TFresults : list or array
        a list of lists contining 'enrichment' scores, normalized 'enrichment' 
        scores, p-value, p-adj, and number of hits for each individual motif

    COMBINEtime : float
        the time it took to combine and merge the bed files using bedtools

    COUNTtime : float
        the time it took to count reads over regions of interest

    DESEQtime : float
        the time it took to perform DE-Seq using the counts file

    CALCULATEtime : float
        the time it took to perform TFEA

    Returns
    -------
    None
    '''
    #Using a config file
    import config
    outputdir = config.OUTPUTDIR
    padj_cutoff = config.PADJCUTOFF
    singlemotif = config.SINGLEMOTIF
    plot = config.PLOTALL

    #For each TF motif with an PADJ value less than a cutoff, an html file is 
    #created to be used in results.html
    positivelist = [x[0] for x in TFresults 
                    if x[2] > 0 and (plot or x[-1] < padj_cutoff)]
    negativelist = [x[0] for x in TFresults 
                    if x[2] < 0 and (plot or x[-1] < padj_cutoff)]

    for i in range(len(TFresults)):
        motif, enrichment, hits, p-value, p-adj = results[i] 
        if plot or p-adj < padj_cutoff or singlemotif:
            if enrichment > 0:
                try:
                    NEXT_MOTIF = positivelist[positivelist.index(motif)+1]
                except IndexError:
                    NEXT_MOTIF = positivelist[0]
                try:
                    PREV_MOTIF = positivelist[positivelist.index(motif)-1]
                except IndexError:
                    PREV_MOTIF = positivelist[len(positivelist)]
            else:
                try:
                    NEXT_MOTIF = negativelist[negativelist.index(motif)+1]
                except IndexError:
                    NEXT_MOTIF = negativelist[0]
                try:
                    PREV_MOTIF = negativelist[negativelist.index(motif)-1]
                except IndexError:
                    PREV_MOTIF = negativelist[len(negativelist)]
            direct_logo = motif.strip('HO_') + "_direct.png"
            reverse_logo = motif.strip('HO_') + "_revcomp.png"
            outfile = open(os.path.join(outputdir, 'plots', motif 
                            + '.results.html'),'w')
            outfile.write("""<!DOCTYPE html>
    <html>
    <head>
    <title>"""+motif+""" Results</title>
    <style>
        table {
            font-family: arial, sans-serif;
            border-collapse: collapse;
            width: 100%;
        }

        .row {
        display: flex; /* equal height of the children */
        width: 100%;
        padding-bottom: 50px
        }

        img {
            max-width: 100%;
            max-height: 100%;
        }

        td, th {
            border: 1px solid #dddddd;
            text-align: left;
            padding: 8px;
        }

        tr:nth-child(even) {
            background-color: #dddddd;
        }
    </style>
    </head>
    <body style="width:1300px; margin:0 auto;">
        <div>
            <div style="float:left">
                <a href="./"""+PREV_MOTIF+""".results.html">PREV</a>
            </div>
            <div style="float:right">
                <a href="./"""+NEXT_MOTIF+""".results.html">NEXT</a>
            </div>
            <div style="text-align:center">
                <a href="../results.html">ALL</a>
        </div>
        <div class="row">
        </div>
            <h1>"""+motif+""" Results</h1>
        <div>
            <div style="float: middle; width: 1300px; padding-bottom:25px; \
                padding-top:25px">
                <table> 
                    <tr>
                        <th>TF Motif</th>
                        <th>AUC</th> 
                        <th>P-value</th>
                        <th>PADJ</th>
                        <th>HITS</th>
                    </tr>
                    <tr>
                        <td>"""+motif+"""</td>
                        <td>"""+str("%.2f" % enrichment)+"""</td>
                        <td>"""+str("%.4g" % p-value)+"""</td>
                        <td>"""+str("%.4g" % p-adj)+"""</td>
                        <td>"""+str(POS)+"""</td>
                    </tr>
                </table>
            </div>
        </div>
        <div>
            <div style="float: left; width 1250px; padding-bottom:50px; \
                padding-top:50px">
                <img src="./"""+motif+"""_enrichment_plot.png" \
                    alt="Enrichment Plot">
            </div>
        </div>
        <div class="row">
            <div style="float: right; width: 600px">
                <p>Forward:</p>
                <img src="./"""+direct_logo+"""" \
                    alt="Forward Logo">
                <p></p>
                <p>Reverse:</p>
                <img src="./"""+reverse_logo+"""" \
                    alt="Reverse Logo">
            </div>
            <div style="float:left; width: 600px">
                <img src="./"""+motif+"""_simulation_plot.png" \
                    alt="Simulation Plot">
            </div>
        </div>

    </body>
    </html>""")
            outfile.close()
            PREV_MOTIF = motif

#==============================================================================
def plot_individual_graphs(plot=None, padj_cutoff=None,
                            figuredir=None, logos=None, 
                            largewindow=None, score=None, 
                            smallwindow=None,
                            distances_abs=None, sorted_distances=None,
                            ranks=None, pvals=None, fc=None, 
                            cumscore=None, motif_file=None, p=None,
                            simES=None, actualES=None, gc_array=None,
                            meta_profile_dict=None, label1=None, label2=None):
    '''This function plots all TFEA related graphs for an individual motif

    Parameters
    ----------
    plot : boolean
        switch to determine whether all motifs are plotted or just significant
        ones
    
    padj_cutoff : float
        the cutoff value that determines signficance. This is technically
        evaluated against the non-adjusted p-value since p-adj has not been
        calculated yet at this point
    
    figuredir : string
        the full path to the figuredir within the ouptut directory containing
        all figures and plots

    logos : string
        full path to a directory containing meme logos. These are copied to the
        output directory to be displayed in the html results.

    largewindow : float
        a user specified larger window used for plotting purposes and to do
        some calculations regarding the user-provided regions of interest

    smallwindow : float
        a user specified smaller window used for plotting purposes and to do
        some calculations regarding the user-provided regions of interest

    distances_abs : list or array
        a sorted (based on rank) list of absolute motif to region center 
        distances

    sorted_distances : list or array
        a sorted (based on rank) list of motif distances. Negative corresponds
        to upstream of region

    ranks : list or array
        a list of ranks that correspond to each region of interest. The ranking
        comes from DE-Seq p-value

    pvals : list or array
        a list of DE-Seq p-values for all regions

    p : float
        the p-value of the given motif based on simulations

    simES : list
        a list of simulated 'enrichment' scores calculated by randomizing
        region rank

    actualES : float
        the observed 'enchrichment' score. Can be calculated in many different
        ways..

    gc_array : list
        an array of gc richness of regions of interest. It is recommended that
        this array be no larger than 1000 bins.

    Returns
    -------
    None
    '''
    #Only plot things if user selects to plot all or if the pvalue is less than
    #the cutoff
    if plot or p < padj_cutoff:
        if config.HOMER:
            pass
        else:
            #Get motif logos from logo directory
            os.system("scp '" + logos 
                        + motif_file.strip('.sorted.distance.bed').strip('HO_') 
                        + "_direct.png' " + figuredir)

            os.system("scp '" + logos 
                        + motif_file.strip('.sorted.distance.bed').strip('HO_') 
                        + "_revcomp.png' " + figuredir)
        

        #Filter distances into quartiles for plotting purposes
        q1 = int(round(np.percentile(np.arange(1, len(sorted_distances),1), 25)))
        q2 = int(round(np.percentile(np.arange(1, len(sorted_distances),1), 50)))
        q3 = int(round(np.percentile(np.arange(1, len(sorted_distances),1), 75)))
        q1_distances = [x for x in sorted_distances[:q1] if x <= largewindow]
        q2_distances = [x for x in sorted_distances[q1:q2] if x <= largewindow]
        q3_distances = [x for x in sorted_distances[q2:q3] if x <= largewindow]
        q4_distances = [x for x in sorted_distances[q3:] if x <= largewindow]

        
        #Get log pval to plot for rank metric
        sorted_pval = [x for _,x in sorted(zip(ranks, pvals))]
        sorted_fc = [x for _,x in sorted(zip(ranks, fc))]
        logpval = list()
        for x,y in zip(sorted_pval,sorted_fc):
            try:
                if y < 1:
                    logpval.append(math.log(x,10))
                else:
                    logpval.append(-math.log(x,10))
            except:
                logpval.append(0.0)

        #Plot the enrichment plot
        independent_functions.enrichment_plot(largewindow=largewindow,
                                            smallwindow=smallwindow,
                                            figuredir=figuredir,
                                            cumscore=cumscore, 
                                            sorted_distances=sorted_distances, 
                                            logpval=logpval,  
                                            gc_array=gc_array, score=score, 
                                            motif_file=motif_file, 
                                            q1_distances=q1_distances, 
                                            q2_distances=q2_distances, 
                                            q3_distances=q3_distances, 
                                            q4_distances=q4_distances,
                                            meta_profile_dict=meta_profile_dict,
                                            dpi=None, label1=label1, 
                                            label2=label2)

        #Plot the simulation plot
        independent_functions.simulation_plot(figuredir=figuredir, simES=simES,
                                            actualES=actualES,
                                            motif_file=motif_file)

        # #Plot the distance distribution histograms plot
        # independent_functions.meta_plot(figuredir=figuredir, 
        #                                     meta_profile_dict=meta_profile_dict,
        #                                     motif_file=motif_file, 
        #                                     q1_distances=q1_distances, 
        #                                     q2_distances=q2_distances, 
        #                                     q3_distances=q3_distances, 
        #                                     q4_distances=q4_distances, 
        #                                     bins=100, 
        #                                     dpi=None)

#==============================================================================
def plot_global_graphs(padj_cutoff=None, label1=None, label2=None, 
                        figuredir=None, TFresults=None):
    '''This function plots graphs that are displayed on the main results.html 
        filethat correspond to results relating to all analyzed TFs.

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

    ESlist = [i[1] for i in TFresults]
    NESlist = [i[2] for i in TFresults]
    PVALlist = [i[3] for i in TFresults]
    POSlist = [i[4] for i in TFresults]
    PADJlist = [i[5] for i in TFresults]

    POSlist = [math.log(x,10) if x > 0 else 0.0 for x in POSlist]

    sigx = [x for x, p in zip(ESlist, PADJlist) if p < padj_cutoff]
    sigy = [p for x, p in zip(ESlist, PADJlist) if p < padj_cutoff]

    MAy = sigx
    MAx = [x for x, p in zip(POSlist, PADJlist) if p < padj_cutoff]

    #Creates a moustache plot of the global PADJs vs. ESs                                                                                                                                                                                       
    independent_functions.moustache_plot(figuredir=figuredir, ESlist=ESlist,
                                            PADJlist=PVALlist, sigx=sigx, 
                                            sigy=sigy)

    #Creates a histogram of p-values                                                                                                                                                                                                            
    independent_functions.pval_histogram_plot(figuredir=figuredir, 
                                                PVALlist=PVALlist)

    #Creates an MA-plot with NES on Y-axis and positive hits on X-axis                                                                                                                                                                                                      
    independent_functions.MA_plot(POSlist=POSlist, ESlist=ESlist, MAx=MAx, 
                                    MAy=MAy, figuredir=figuredir, label1=label1,
                                    label2=label2)