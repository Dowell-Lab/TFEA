__author__ = 'Jonathan Rubin'

import os

def runfimo(ranked_file,filedir,MEMEDB,DATABASE,SINGLEMOTIF):
    #This os.system call uses bedtools to convert the ranked_file.bed into fasta format (ranked_file.fasta)
    # command = "bedtools getfasta -fi " + GENOME + " -bed " + filedir + "ranked_file.bed -fo " + filedir + "ranked_file.fasta"
    # print command
    # subprocess.call(command,shell=True)
    # exit_code = os.system("bedtools getfasta -fi " + GENOME + " -bed " + filedir + "ranked_file.bed -fo " + filedir + "ranked_file.fasta")
    # print exit_code
    # print os.environ

    if os.path.isdir(filedir + "fimo_out/"):
        os.system("rm -r " + filedir + "fimo_out/")
    if SINGLEMOTIF == False:
        os.system("fimo -o " + filedir + "fimo_out/ " + MEMEDB + DATABASE + " " + ranked_file)
    else:
        os.system("fimo -o " + filedir + "fimo_out/ " + "--motif " + SINGLEMOTIF + " " + MEMEDB + DATABASE + " " + ranked_file)

#Takes a single ranked_file and returns bed file with the centers for each region and a bed file with the distance to the closest motif for each centered region
def run(ranked_file,filedir,MOTIF_PATH):
    #Get center base for each region
    outfile = open(filedir+"ranked_file.center.bed",'w')
    with open(ranked_file) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            chrom,start,stop = line[:3]
            center = (int(start)+int(stop))/2
            outfile.write(chrom + '\t' + str(center) + '\t' + str(center+1) + '\t' + '\t'.join(line[3:]) + '\n')
    outfile.close()

    #Get closest motif hit to center base of each region
    os.system("sort -k1,1 -k2,2n " + filedir + "ranked_file.center.bed > " + filedir + "ranked_file.center.sorted.bed")
    os.system("bedtools closest -d -t first -a " + filedir + "ranked_file.center.sorted.bed -b " + MOTIF_PATH + " > " + filedir + "ranked_file.center.sorted.distance.bed")
