__autor__ = 'Rutendo F. Sigauke'

import pandas as pd

def MDScalc(MDS1):

    with open(MDS1) as F:
        
        nm1 = MDS1.split('/')[-1].split('_')[0]
        
        tf = list()
        motifs = list()
        hsum = list()
        Hsum = list()

        for line in F:
            if 'HO_' in line:
                tf.append(line.split(',')[0])
                motifs.append(line.split(',')[1:3000])
        for motif in motifs:
            Hsum.append(sum(int(i) for i in motif))
            hsum.append(sum(int(i) for i in motif[149:299]))
        
        mds = [float(x)/float(y) if x > 0 and y > 0 else 0 for x,y in zip(hsum, Hsum)] 

        MDSzip = zip(tf, hsum, Hsum, mds)
        MDS = pd.DataFrame(MDSzip,columns=['TF_motif','h_hits','H_hits','MDS'])
        MDS.to_csv(savedir + nm1 + '_MDSout.tsv',sep='\t')

if __name__ == "__main__":
    MDSdir = '/scratch/Shares/dowell/md_score_paper/mdscore_tsv_files/human/recent/'
    ##MDS1 = MDSdir + 'SRR1105737_MDS.csv'#rep 1
    MDS1 = MDSdir + 'SRR1105739_MDS.csv' #rep2
    savedir = '/scratch/Users/rusi2317/projects/rotation/output/MDS/'
    MDScalc(MDS1)


