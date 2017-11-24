__author__ = 'Jonathan Rubin'

import pybedtools as py

def run(DMSO,Nutlin,deseq):
    regions = list()
    with open(deseq) as F:
        F.readline()
        for line in F:
            line = line.strip('\n').split('\t')
            regions.append((line[1],format(float(line[-2]),'.12f')))

    d = dict()
    with open(DMSO) as F:
        i = 0
        for line in F:
            if '#' not in line[0]:
                if '>' in line[0]:
                    interval = list()
                    line = line.strip('\n').split('|')
                    chrom = line[1].split(':')[0]
                    start,stop = line[1].split(':')[1].split('-')
                    if chrom not in d:
                        d[chrom] = list()
                    interval.append(chrom)
                    interval.append(int(start))
                    interval.append(int(stop))
                    i = 0
                if i == 2:
                    line = line.strip('\n').split('\t')
                    mu = float(line[1])
                    interval.append(mu)
                    d[interval[0]].append(interval[1:])
                i += 1

    print d




if __name__ == "__main__":
    DMSO = '/scratch/Shares/dowell/md_score_paper/tfit_bed_files/human/recent/SRR1105737-1_K_models_MLE.tsv'
    Nutlin = '/scratch/Shares/dowell/md_score_paper/tfit_bed_files/human/recent/SRR1105739-1_K_models_MLE.tsv'
    deseq = '/Users/joru1876/scratch_backup/TFEA_grant_figure/Tfit_SRR3739_bidir_predictions.sorted.merge.count.bed.id.bed.DMSONutlinnascent.res.txt'
    run(DMSO,Nutlin,deseq)