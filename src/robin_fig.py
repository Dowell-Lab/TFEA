__author__ = 'Jonathan Rubin'

import math

def get_intersect_mu(start1,stop1,int_list):
    for region in int_list:
        start2,stop2,mu1 = region
        if start1 < start2 < stop1 or start1 < stop2 < stop1 or (start2 < start1 and stop2 > stop1):
            return start2,stop2,mu1

    return False

def closest_distance(mu,center_list):
    results = list()
    for val in center_list:
        results.append(math.fabs(val-mu))

    return min(results)
        
def run(DMSO,Nutlin,deseq,P53):
    regions = list()
    with open(deseq) as F:
        F.readline()
        for line in F:
            line = line.strip('\n').split('\t')
            pval = format(float(line[-2]),'.12f')
            chrom,start,stop = line[1].split(',')
            if pval < 0.05:
                regions.append((chrom,int(start),int(stop),pval))

    sorted_regions = sorted(regions,key=lambda x: x[3])
    print regions

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

    d2 = dict()
    with open(Nutlin) as F:
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

    d3 = dict()
    with open(P53) as F:
        for line in F:
            if '#' not in line[0]:
                line = line.strip('\n').split('\t')
                chrom,start,stop = line[:3]
                center = (int(start)+int(stop))/2
                if chrom not in d3:
                    d3[chrom] = list()
                d3[chrom].append(center)

    print "Parsed all files"

    total = 0
    hits = 0
    for region in sorted_regions:
        total += 1
        chrom,start,stop,pval = region
        DMSO = d[chrom]
        Nutlin = d2[chrom]
        P53 = d3[chrom]

        value1 = get_intersect_mu(start,stop,DMSO)
        if value1 != False:
            start2,stop2,mu1 = value1
            value2 = get_intersect_mu(start2,stop2,Nutlin)
            if value2 != False:
                mu2 = value2[-1]
                truemu = (mu1+mu2)/2.0
                if closest_distance(truemu,P53) < 1500:
                    hits += 1
            else:
                if closest_distance(mu1,P53) < 1500:
                    hits += 1
        else:
            value2 = get_intersect_mu(start,stop,Nutlin)
            if value2 != False:
                if closest_distance(value2[-1],P53) < 1500:
                    hits += 1
            else:
                pass



    print "Total regions p < 0.05: ", total
    print "Total P53 hits H < 1500: ", hits


if __name__ == "__main__":
    DMSO = '/scratch/Shares/dowell/md_score_paper/tfit_bed_files/human/recent/SRR1105737-1_K_models_MLE.tsv'
    Nutlin = '/scratch/Shares/dowell/md_score_paper/tfit_bed_files/human/recent/SRR1105739-1_K_models_MLE.tsv'
    P53 = '/scratch/Shares/dowell/md_score_paper/PSSM_hits_genome_wide/pval-6/HO_P53_HUMAN.H10MO.B.bed'
    deseq = '/Users/joru1876/scratch_backup/TFEA_grant_figure/Tfit_SRR3739_bidir_predictions.sorted.merge.count.bed.id.bed.DMSONutlinnascent.res.txt'
    run(DMSO,Nutlin,deseq,P53)