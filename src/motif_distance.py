__author__ = 'Jonathan Rubin'

import os
# from pybedtools import BedTool

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
            chrom,start,stop = line.strip('\n').split('\t')
            center = (int(start)+int(stop))/2
            outfile.write(chrom + '\t' + str(center) + '\t' + str(center+1) + '\n')

    #Get closest motif hit to center base of each region
    command = "bedtools closest -d -a " + filedir + "ranked_file.center.bed -b " + MOTIF_PATH + " > " + filedir + "ranked_file.center.distance.bed"
    exit_code = os.system(command)
    print exit_code
    print command

def run_pybedtools(ranked_file,filedir,MOTIF_HITS,SINGLEMOTIF):
    #Get center base for each region
    outfile = open(filedir+"ranked_file.center.bed",'w')
    with open(ranked_file) as F:
        for line in F:
            chrom,start,stop = line.strip('\n').split('\t')
            center = (int(start)+int(stop))/2
            outfile.write(chrom + '\t' + str(center) + '\t' + str(center+1) + '\n')

    #Get closest motif hit to center base of each region

    # if SINGLEMOTIF == False:
    #     for MOTIF_FILE in os.listdir(MOTIF_HITS):
    #         a = BedTool(filedir + "ranked_file.center.bed")
    #         b = BedTool(MOTIF_HITS + MOTIF_FILE)
    #         c = a.closest(b,d=True)
    #         c.saveas(filedir + "ranked_file.center.distance.bed")
    # else:
    #     for MOTIF_FILE in os.listdir(MOTIF_HITS):
    #         if MOTIF_FILE == SINGLEMOTIF:
    #             a = BedTool(filedir + "ranked_file.center.bed")
    #             b = BedTool(MOTIF_HITS + MOTIF_FILE)
    #             c = a.closest(b,d=True)
    #             c.saveas(filedir + "ranked_file.center.distance.bed")
    #             # command = "bedtools closest -d -a " + filedir + "ranked_file.center.bed -b " + MOTIF_HITS + MOTIF_FILE + " > " + filedir + "ranked_file.center.distance.bed"
    #             # exit_code = os.system(command)
    #             # print exit_code
    #             # print command