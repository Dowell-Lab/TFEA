__autor__ = 'Rutendo F. Sigauke'

import pandas as pd 
import scipy
from scipy import stats



def sortTFEA_MDS(tfea_res, mds_res, outfile):
    """
    take a tsv files for tfea and mds score, sort by the pvalues and 
    use the rank of the both to calculate the spearman correlation.

    """

    tfs = list()
    nes = list()
    pval = list()
    #name = 'Allen_2014'
    #name = 'Luo_2014'
    name = 'Hah_2013'
    exp = 'DMSO_nutlin'

    with open(tfea_res) as tfea:
        for line in tfea:
            line = line.strip('\n')
            tfs.append(line.split('\t')[0])
            nes.append(line.split('\t')[2])
            pval.append(line.split('\t')[3])
        TFEAsum = zip(tfs, nes, pval)
        TFEAsummary = pd.DataFrame(TFEAsum,columns=['TF_motif','NES','pvals'])
        TFEAsummary2 = TFEAsummary.drop(TFEAsummary.index[:1])
        TFEAsummary2['TF'] = TFEAsummary2['TF_motif'].str.split('_', 1).str[1].str.split('_',1).str[0]
        TFEAsummary2['pvals'] = TFEAsummary2['pvals'].astype('float')
        
        ##split the NES into positive and negative values
        ##now sort the two df by the pvalues
        positive = TFEAsummary2[(TFEAsummary2['NES'].astype(float) >= 0.0)]
        negative = TFEAsummary2[TFEAsummary2['NES'].astype(float) < 0.0]
        positive_sort = positive.sort_values(by=['pvals'])
        negative_sort = negative.sort_values(by=['pvals'],ascending=False)
        
        dframes = [positive_sort, negative_sort]
        
        sorted_pval_tfea = pd.concat(dframes)
        #sorted_pval_tfea.insert(0, 'Rank1', range(1, 1 + len(sorted_pval_tfea)))
        sorted_pval_tfea.insert(0, 'RankTFEA', range(1, 1 + len(sorted_pval_tfea)))
        
        
    tf = list()
    Y = list()
    X = list()
    pvals = list()
    tf_motif =list()
    deltaMDS = list()
    
    with open(mds_res) as mds:
        for line in mds:
            line = line.strip('\n')
            tf.append(line.split('\t')[1])
            deltaMDS.append(line.split('\t')[2])
            X.append(line.split('\t')[3])
            Y.append(line.split('\t')[4])
            pvals.append(line.split('\t')[5])
            tf_motif.append(line.split('\t')[7])

        MDSsum = zip(tf, deltaMDS,X, Y, pvals, tf_motif)
        MDSsummary = pd.DataFrame(MDSsum,columns=['TF','deltaMDS','log_meanHits_X','Y','pvals', 'TF_motif'])
        MDSsummary2 = MDSsummary.drop(MDSsummary.index[:1])
        MDSsummary2['pvals'] =MDSsummary2['pvals'].astype('float')
        
        ##split the NES into positive and negative values
        ##now sort the two df by the pvalues
        positive = MDSsummary2[(MDSsummary2['deltaMDS'].astype(float) >= 0.0)]
        negative = MDSsummary2[MDSsummary2['deltaMDS'].astype(float) < 0.0]
        positive_sort = positive.sort_values(by=['pvals'])
        negative_sort = negative.sort_values(by=['pvals'],ascending=False)
        
        dframes = [positive_sort, negative_sort]

        sorted_pval_mds = pd.concat(dframes)
        #sorted_pval_mds.insert(0, 'Rank2', range(1, 1 + len(sorted_pval_mds)))
        sorted_pval_mds.insert(0, 'RankMDS', range(1, 1 + len(sorted_pval_mds)))

        df = pd.merge(sorted_pval_mds, sorted_pval_tfea, on='TF_motif',how='inner')
        print df

        #df.to_csv(outfile + name +'MDS_TFEA_Ranks.tsv', sep='\t')
        #df.to_csv(outfile + name +'MDS_TFEA_Ranks_12.tsv', sep='\t')
        #df.to_csv(outfile + name +'MDS_TFEA_Ranks_15.tsv', sep='\t')
        #df.to_csv(outfile + name +'MDS_TFEA_Ranks_16.tsv', sep='\t')
        # df.to_csv(outfile + name +'MDS_TFEA_Ranks_3.tsv', sep='\t')
        #df.to_csv(outfile + name +'MDS_TFEA_Ranks_4.tsv', sep='\t')
        #df.to_csv(outfile + name +'MDS_TFEA_Ranks_20.tsv', sep='\t')
        #df.to_csv(outfile + name +'MDS_TFEA_Ranks_20.tsv', sep='\t')
        #df.to_csv(outfile + name +'MDS_TFEA_Ranks_26.tsv', sep='\t')
        #df.to_csv(outfile + name +'MDS_TFEA_Ranks_25.tsv', sep='\t')
        ##df.to_csv(outfile + name +'MDS_TFEA_Ranks_5.tsv', sep='\t')

        ##spearman_corr = scipy.stats.spearmanr(df["Rank1"].values, df["Rank2"].values)
        spearman_corr = scipy.stats.spearmanr(df["RankTFEA"].values, df["RankMDS"].values)
        ##spearman_corr = scipy.stats.pearsonr(df["Rank1"].values, df["Rank2"].values) 
        
        print spearman_corr


if __name__ == "__main__":

    #tfea_res = '/scratch/Users/rusi2317/projects/rotation/output/TFEA/Allen2014/TFEA_output-11/results.txt'
    #tfea_res = '/scratch/Users/rusi2317/projects/rotation/output/TFEA/Allen2014/TFEA_output-12/results.txt' 
    #tfea_res = '/scratch/Users/rusi2317/projects/rotation/output/TFEA/Allen2014/TFEA_output-15/results.txt'
    #tfea_res = '/scratch/Users/rusi2317/projects/rotation/output/TFEA/Luo2014/TFEA_output-3/results.txt'
    #tfea_res = '/scratch/Users/rusi2317/projects/rotation/output/TFEA/Luo2014/TFEA_output-4/results.txt'
    ###tfea_res = '/scratch/Users/rusi2317/projects/rotation/output/TFEA/Luo2014/TFEA_output-20/results.txt'
    #tfea_res = '/scratch/Users/rusi2317/projects/rotation/output/TFEA/Allen2014/TFEA_output-20/results.txt'
    ##tfea_res = '/scratch/Users/rusi2317/projects/rotation/output/TFEA/Allen2014/TFEA_output-26/results.txt'
    tfea_res = '/scratch/Users/rusi2317/projects/rotation/output/TFEA/Allen2014/TFEA_output-33/results.txt'
    #tfea_res = '/scratch/Users/rusi2317/projects/rotation/output/TFEA/Luo2014/TFEA_output-25/results.txt'
    ##tfea_res = '/scratch/Users/rusi2317/projects/rotation/output/TFEA/Hah2013/TFEA_output-5/results.txt'
    ##tfea_res = '/scratch/Users/rusi2317/projects/rotation/output/TFEA/Allen2014/TFEA_output-16/results.txt'
    ##mds_res = '/scratch/Users/rusi2317/projects/rotation/output/MDS_compare/SRR1015587SRR1015583_MDSout.tsv'
    #mds_res = '/scratch/Users/rusi2317/projects/rotation/output/MDS_compare/SRR1015587SRR1015583_MDSout2.tsv'
    ##mds_res = '/scratch/Users/rusi2317/projects/rotation/output/MDS_compare/SRR1015587_SRR1015583_MDSout_04202018.tsv'
    mds_res = '/scratch/Users/rusi2317/projects/rotation/output/MDS_compare/SRR1105739_SRR1105737_MDSout_04232018.tsv'
    ##mds_res = '/scratch/Users/rusi2317/projects/rotation/output/MDS_compare/SRR653425_SRR653421_MDSout_04302018.tsv'
    outfile ='/scratch/Users/rusi2317/projects/rotation/output/spearman_ranks/'
    sortTFEA_MDS(tfea_res, mds_res, outfile)


