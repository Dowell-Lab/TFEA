__autor__ = 'Jonathan Rubin'

import os
import math
from scipy.stats import norm
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import pandas as pd


def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])

    return newdir

def run(MDS1,MDS2,savedir):
    name1 = MDS1.split('/')[-1].split('_')[0]
    name2 = MDS2.split('/')[-1].split('_')[0]
    d = dict()
    with open(MDS1) as F:
        for line in F:
            if '#Binned' in line:
                break
            if '#' not in line[0]:
                line = line.strip().split()
                d[line[0]] = [line[1].split(','),line[2].split(',')]

    with open(MDS2) as F:
        for line in F:
            if '#Binned' in line:
                break
            if '#' not in line[0]:
                line = line.strip().split()
                d[line[0]].append(line[1].split(','))
                d[line[0]].append(line[2].split(','))

    namelist = ['NON','TSS','COMB']
    

    for i in range(len(namelist)):
        name = namelist[i]
        X = list()
        Y = list()
        X2 = list()
        Y2 = list()
        X3 = list()
        Y3 = list()
        diff = list()
        siglist = list()
        siglist2 = list()
        genelist = list()
        tf_motifs = list()

        ps = list()
        zs = list()
        for key in d:
            mdj=float(d[key][1][i])
            mdk=float(d[key][3][i])
            Nj=float(d[key][0][i])
            Nk=float(d[key][2][i])
            diff.append(float(mdj)-float(mdk))
            if 'PPARA' in key:
                    print key,mdj-mdk,math.log((Nj+Nk)/2.0,10)
        mean = sum(diff)/len(diff)
        print mean
        for key in d:
            mdj=float(d[key][1][i])
            mdk=float(d[key][3][i])
            Nj=float(d[key][0][i])
            Nk=float(d[key][2][i])
            # if key.split('.')[0].split('_')[1] == 'SRF':
            #     X2.append(mdj-mdk-mean)
            #     Y2.append(math.log((Nj+Nk)/2.0,10))
            if (Nj+Nk)/2.0 > 10:
                p=((mdj*Nj)+(mdk*Nk))/(Nj+Nk)
                SE=(p*(1-p))*((1/Nj)+(1/Nk))
                Y.append(mdj-mdk-mean)
                X.append(math.log((Nj+Nk)/2.0,10))
                tf_motifs.append(key)
                genelist.append(key.split('.')[0].split('_')[1])
                try:
                    z = (mdj-mdk-mean)/math.sqrt(SE)
                except ZeroDivisionError:
                    z=0
                zs.append(z)
                cdf=norm.cdf(z,0,1)
                p=min(cdf,1-cdf)*2
                ps.append(p)
                if p < 0.01:
                    if mdj-mdk-mean > 0:
                        X2.append(math.log((Nj+Nk)/2.0,10))
                        Y2.append(mdj-mdk-mean)
                        siglist.append(key.split('.')[0].split('_')[1])
                        ##print 'up', key.split('.')[0].split('_')[0], math.log((Nj+Nk)/2.0,10), mdj-mdk-mean
                        ##print 'up', key.split('.')[0].split('_')[0], "%.4f" % p
                    else:
                        X3.append(math.log((Nj+Nk)/2.0,10))
                        Y3.append(mdj-mdk-mean)
                        siglist2.append(key.split('.')[0].split('_')[1])
                        # print 'down', key.split('.')[0].split('_')[0], math.log((Nj+Nk)/2.0,10), mdj-mdk-mean
                        ##print 'down', key.split('.')[0].split('_')[0], "%.4f" % p
                if 'PPARA' in key:
                    print key,mdj-mdk,math.log((Nj+Nk)/2.0,10)
        if name == 'NON':
            MDSsum = zip(genelist,diff ,X, Y, ps, zs,tf_motifs)
            MDSsummary = pd.DataFrame(MDSsum,columns=['TF','DeltaMDS','log_meanHits_X','Y','pval','ztest','TF_motif'])
            ##MDSsummary.to_csv(savedir+name1+'_'+name2 +'_MDSout_04232018.tsv',sep='\t')
            MDSsummary.to_csv(savedir+name1+'_'+name2 +'_MDSout_04302018.tsv',sep='\t')

            #print MDSsummary                                                                                                                                                     

if __name__ == "__main__":
    #Home directory
    homedir = os.path.dirname(os.path.realpath(__file__))

    #File directory
    ##filedir = parent_dir(homedir) + '/files/'
    ##MDS1 = parent_dir(homedir) + '/MDS_files/J12_MDS.tsv'
    ##MDSdir = '/scratch/Shares/dowell/md_score_paper/mdscore_tsv_files/human/recent/'
    MDSdir = '/scratch/Shares/dowell/md_score_paper/mdscore_tsv_files/human/archive/150/'

    #MDS1 = MDSdir + 'SRR1105737_MDS.tsv'
    #MDS1 = MDSdir + 'SRR1105736_MDS.tsv'
    #MDS1 = MDSdir + 'SRR1015583_MDS.tsv'
    MDS1 = MDSdir + 'SRR653421_MDS.tsv'
    # MDS1 = MDSdir + 'SRR1015584_MDS.tsv'
    ##MDS2 = parent_dir(homedir) + '/MDS_files/J32_MDS.tsv'
    #MDS2 = MDSdir + 'SRR1105739_MDS.csv'
    #MDS2 = MDSdir + 'SRR1105739_MDS.tsv'
    #MDS2 = MDSdir + 'SRR1105738_MDS.tsv'
    #MDS2 = MDSdir + 'SRR1015587_MDS.tsv'
    MDS2 = MDSdir + 'SRR653425_MDS.tsv'
    # MDS2 = MDSdir + 'SRR1015588_MDS.tsv'
    ##savedir = parent_dir(homedir) + '/figures/'
    savedir = '/scratch/Users/rusi2317/projects/rotation/output/MDS_compare/'
    run(MDS2,MDS1,savedir)
