__author__ = 'Jonathan Rubin'

from config import BATCH,COMBINE,COUNT,RANK,DISTANCE,CALCULATE,FDRCUTOFF,BEDS,BAM1,BAM2,SINGLEMOTIF,DATABASE,GENOME,MEMEDB,MOTIF_HITS,DESEQFILE,CELLTYPE,LABEL1,LABEL2

def run(TFresults,outputdir,COMBINEtime,COUNTtime,RANKtime,DISTANCEtime,CALCULATEtime):
    #summary.html contains all user-defined variables, and also information about module used
    outfile = open(outputdir+'summary.html','w')
    outfile.write("""<!DOCTYPE html>
            <html>
            <head>
            <title>Variables Used</title>
            </head>
            <body>
                <h1>Variables Used</h1>
                <p>BEDS = """+str(BEDS)+"""</p>
                <p>BAM1 = """+str(BAM1)+"""</p>
                <p>BAM2 = """+str(BAM2)+"""</p>
                <p>SINGLEMOTIF = """+str(SINGLEMOTIF)+"""</p>
                <p>DATABASE = """+str(DATABASE)+"""</p>
                <p>GENOME = """+str(GENOME)+"""</p>
                <p>MEMEDB = """+str(MEMEDB)+"""</p>
                <p>MOTIF_HITS = """+str(MOTIF_HITS)+"""</p>
                <p>DESEQFILE = """+str(DESEQFILE)+"""</p>
                <p>CELLTYPE = """+str(CELLTYPE)+"""</p>
            </body>""")

    #For each TF motif with an FDR value less than a cutoff, an html file is created to be used in results.html
    for MOTIF_FILE,ES,NES,PVAL,FDR in TFresults:
        if FDR < FDRCUTOFF:
            outfile = open(outputdir + 'plots/' + MOTIF_FILE + '.results.html','w')
            outfile.write("""<!DOCTYPE html>
            <html>
            <head>
            <title>"""+MOTIF_FILE+""" Results</title>
            <style>
            table {
                font-family: arial, sans-serif;
                border-collapse: collapse;
                width: 100%;
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
            <body>
            <h1>"""+MOTIF_FILE+""" Results</h1>
            <div>
                <div id="Positively Enriched" style="float: left; width: 1000px; overflow:scroll">
                    <table> 
                        <tr>
                            <th>TF Motif</th>
                            <th>ES</th> 
                            <th>NES</th>
                            <th>P-value</th>
                            <th>FDR</th>
                        </tr>
                        <tr>
                            <td>"""+MOTIF_FILE+"""</td>
                            <td>"""+str("%.3f" % ES)+"""</td>
                            <td>"""+str("%.3f" % NES)+"""</td>
                            <td>"""+str("%.4g" % PVAL)+"""</td>
                            <td>"""+str("%.4g" % FDR)+"""</td>
                        </tr>
            <img src="./"""+MOTIF_FILE+"""_enrichment_plot.svg" alt="Enrichment Plot">
            <img src="./"""+MOTIF_FILE+"""_simulation_plot.svg" alt="Simulation Plot">
            </body>
            </html>""")
            outfile.close()

    outfile = open(outputdir+'results.html','w')
    outfile.write("""<!DOCTYPE html>
    <html>
    <head>
    <title>TFEA Results</title>
    <style>
    table {
        font-family: arial, sans-serif;
        border-collapse: collapse;
        width: 100%;
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
    <body style="width: 1250px; overflow:scroll">

    <h1>TFEA Results """ +LABEL1+ """ vs. """ +LABEL2+ """</h1>
    <div>
        <div style="float: left">
            <img src="./plots/TFEA_Results_Moustache_Plot.svg" alt="Moustache Plot (FDR vs. NES)">
        </div>
        <div id="Summary of Variables Used" style="float: right; width: 400px">
            <p><a href="./Summary.html">Full Summary of Variables Used</a></p>
            <p>FDR cutoff = """ + str(FDRCUTOFF) + """</p>
            <p>Genome = """ + GENOME + """</p>
            <p>Cell Type = """ + CELLTYPE + """</p>
            <table>
                <tr>
                    <th>Module</th>
                    <th>Switch</th>
                    <th>Time (s)</th>
                </tr>
                <tr>
                    <td>BATCH</td>
                    <td>"""+str(BATCH)+"""</td>
                    <td>N/A</td>
                </tr>
                <tr>
                    <td>COMBINE</td>
                    <td>"""+str(COMBINE)+"""</td>
                    <td>"""+str(COMBINEtime)+"""</td>
                </tr>
                <tr>
                    <td>COUNT</td>
                    <td>"""+str(COUNT)+"""</td>
                    <td>"""+str(COUNTtime)+"""</td>
                </tr>
                <tr>
                    <td>RANK</td>
                    <td>"""+str(RANK)+"""</td>
                    <td>"""+str(RANKtime)+"""</td>
                </tr>
                <tr>
                    <td>DISTANCE</td>
                    <td>"""+str(DISTANCE)+"""</td>
                    <td>"""+str(DISTANCEtime)+"""</td>
                </tr>
                <tr>
                    <td>CALCULATE</td>
                    <td>"""+str(CALCULATE)+"""</td>
                    <td>"""+str(CALCULATEtime)+"""</td>
                </tr>
            </table>
                
        </div>
    </div>
    <div>
        <div id="Positively Enriched" style="float: left; width: 600px; overflow:scroll">
            <h1>Positively Enriched</h1>
            <table> 
                <tr>
                    <th>TF Motif</th>
                    <th>ES</th> 
                    <th>NES</th>
                    <th>P-value</th>
                    <th>FDR</th>
                </tr>
                """)

    for MOTIF_FILE,ES,NES,PVAL,FDR in TFresults:
        if NES > 0:
            if FDR < FDRCUTOFF:
                outfile.write("""                <tr>
                        <td><a href="./plots/"""+MOTIF_FILE+""".results.html">"""+MOTIF_FILE+"""</td>
                        <td>"""+str("%.3f" % ES)+"""</td>
                        <td>"""+str("%.3f" % NES)+"""</td>
                        <td>"""+str("%.4g" % PVAL)+"""</td>
                        <td>"""+str("%.4g" % FDR)+"""</td>
                    </tr>
                    """)
            else:
                outfile.write("""                <tr>
                        <td>"""+MOTIF_FILE+"""</td>
                        <td>"""+str("%.3f" % ES)+"""</td>
                        <td>"""+str("%.3f" % NES)+"""</td>
                        <td>"""+str("%.4g" % PVAL)+"""</td>
                        <td>"""+str("%.4g" % FDR)+"""</td>
                    </tr>
                    """)


    outfile.write("""            </table>
        </div>

        <div id="Negatively Enriched" style="float: right; width: 600px; overflow:scroll">
            <h1>Negatively Enriched</h1>
            <table> 
                <tr>
                    <th>TF Motif</th>
                    <th>ES</th> 
                    <th>NES</th>
                    <th>P-value</th>
                    <th>FDR</th>
                </tr>
                """)

    for MOTIF_FILE,ES,NES,PVAL,FDR in TFresults:
        if NES < 0:
            if FDR < FDRCUTOFF:
                outfile.write("""                <tr>
                        <td><a href="./plots/"""+MOTIF_FILE+""".results.html">"""+MOTIF_FILE+"""</td>
                        <td>"""+str("%.3f" % ES)+"""</td>
                        <td>"""+str("%.3f" % NES)+"""</td>
                        <td>"""+str("%.4g" % PVAL)+"""</td>
                        <td>"""+str("%.4g" % FDR)+"""</td>
                    </tr>
                    """)
            else:
                outfile.write("""                <tr>
                        <td>"""+MOTIF_FILE+"""</td>
                        <td>"""+str("%.3f" % ES)+"""</td>
                        <td>"""+str("%.3f" % NES)+"""</td>
                        <td>"""+str("%.4g" % PVAL)+"""</td>
                        <td>"""+str("%.4g" % FDR)+"""</td>
                    </tr>
                    """)

    outfile.write("""        </table>
        </div>
    </div>

    </body>
    </html>""")

    outfile.close()

if __name__ == "__main__":
    TFresults = [['HO_P53_HUMAN.H10MO.B.bed',0.182143716966,6.22622338072,7.43595407471e-10,0.0]]
    filedir = './'
    figuredir = './'
    run(filedir,TFresults,figuredir)