### DESCRIPTION ###
# This script outputs text files for regions within the leading edges (Plateaued E and Match Background) that 
#    also have a motif within a desired distance (None means all regions within leading edge are maintained)
#    Output is stored as [TF]_motif[motif_dist]_[MB|PE]_LE.bed at the results_dir.
#     header = #chrom, start, stop, name, distance(of motif to µ)

# Run with Rscript Get_LE.R

library(data.table)
#### EDIT HERE ####
TF="SRF.H12CORE.0.PSM.A"
# directory where the results from TFEA are held (and output will be saved)
results_dir="~/Desktop/TF_Enh_Linking/Resp_Env/Comb_UPM_WSP_ADP/results/TFEA/WSP_PRO/NonTSS_WSP30vsVeh_par_rat_Wd_GeneHK_noback"
# ranked file only needed if used as input to TFEA (otherwise TFEA outputs it)
ranked_file = "~/Desktop/TF_Enh_Linking/Resp_Env/Comb_UPM_WSP_ADP/data/TFEA_input/NonTSS_WSP30_vs_Veh_par_rat_Wd_GeneHK.bed"
ranked_file = ""
# Maximum distance for which a motif can be found from mu (midpoint) to be included
# put "None" if you want ALL within the leading edge
motif_dist=500 # "None"

### DO NOT EDIT BELOW ####

read_in_results <- function(file_name) {
    res = as.data.frame(fread(file_name, sep="\t"))
    colnames(res) <- c("TF", "Escore", "CorrEscore", "Events", "GC", "FPKM", "Padj", "CorrPadj", 
                      "MB_LE", "MB_LE_stdev", "PE_LE", "PE_LE_stdev", "Frac_abov_null")
    print(head(res))

    res$CorrPadj <- as.numeric(res$CorrPadj)
    res$Padj <- as.numeric(res$Padj)

    return(res)
    }

get_LE_calls <- function(resdir, ranked_file_use="", TF_use="", motif_dist_use="None") {
    if (TF_use == "") {
        stop("A TF name must be supplied")
    }
    # read in the results to get TF results
    res_calls = read_in_results(paste0(resdir, "/results.txt"))
    res_calls = res_calls[res_calls$TF == TF_use,]
    # read in the distance df (distance of motifs)
    dist_df = fread(paste0(resdir, "/temp_files/withAUC", TF_use, "__dis_andmore.csv"), sep=",")
    # read in the ranked file for the names
    if (ranked_file_use == "") {
        # Get the ranked names from the ranked_file.fa
        #system(sprintf(
        #'grep ">" %s/tempfiles/ranked_file.fa > %s/tempfiles/ranked_file.txt',
       # resdir, resdir
      #))
      outfile <- path.expand(file.path(resdir, "temp_files", "ranked_file.txt"))
      infile  <- path.expand(file.path(resdir, "temp_files", "ranked_file.fa"))
      print(infile)
      system2(
          "grep",
          args = c("'>'", file.path(infile)),
          stdout = file.path(outfile)
        )
        ranked_file_use = paste0(resdir, "/temp_files/ranked_file.txt")
        rank_df_orig = fread(ranked_file_use, sep=":", header=FALSE)
        rank_df_orig$V1 <- gsub(">", "", rank_df_orig$V1)
        rank_df = data.frame(chrom=rank_df_orig$V1,
                            start=sapply(strsplit(rank_df_orig$V2, "-"), "[", 1),
                            stop=sapply(strsplit(rank_df_orig$V2, "-"), "[", 2))

    } else {
      rank_df = fread(ranked_file_use, sep="\t")
    }
    
    # get a name for each
    dist_df$V1 <- seq(1, nrow(dist_df))
    dist_df$chrom = rank_df[,1]
    dist_df$start = rank_df[,2]
    dist_df$stop = rank_df[,3]
    dist_df$name = paste0(dist_df$chrom, ":", dist_df$start, "-", dist_df$stop)
    # only consider those within the leading edge (consider Enrichment score)
    if ((res_calls$Padj[1] > 0.05) | (res_calls$CorrPadj[1] > 0.05)) {
        message("WARNING: the selected TF does not have a significant enrichment, and therefore the LE is likely not accurate")
    }
    # Get leading edges
    PE_LE = res_calls$PE_LE[1]
    MB_LE = res_calls$MB_LE[1]
    # get enrichment score
    ES = res_calls$Escore[1]
    # Only consider calls within the LEs
    if (ES > 0) {
        PE_dist_df = dist_df[1:PE_LE,]
        MB_dist_df = dist_df[1:MB_LE,]
    } else {
        PE_LE = nrow(dist_df) - PE_LE
        MB_LE = nrow(dist_df) - MB_LE
        PE_dist_df = dist_df[PE_LE:nrow(dist_df),]
        MB_dist_df = dist_df[MB_LE:nrow(dist_df),]
    }
    # only consider those with motif within the desired distance
    if (motif_dist_use != "None") {
        # remove cases where no motif
        PE_dist_df = PE_dist_df[PE_dist_df$distances_abs != ".",]
        MB_dist_df = MB_dist_df[MB_dist_df$distances_abs != ".",]
        PE_dist_df$distances_abs <- as.numeric(PE_dist_df$distances_abs)
        # only get cases where motif within certain distance
        PE_dist_df = PE_dist_df[PE_dist_df$distances_abs <= as.numeric(motif_dist_use),]
        MB_dist_df = MB_dist_df[MB_dist_df$distances_abs <= as.numeric(motif_dist_use),]

    }
    # make bed files of regions AND Motif distance
    PE_bed = PE_dist_df[, c("chrom", "start", "stop", "name", "distances")]
    MB_bed = MB_dist_df[, c("chrom", "start", "stop", "name", "distances")]
    # save
    #[TF]_motif[motif_dist]_[MB|PE]_LE.bed
    PE_output_name = paste0(resdir, "/", TF_use, "_motif", as.character(motif_dist), "_PE_LE.txt")
    write.table(PE_bed, file=PE_output_name, row.names=FALSE, sep="\t", quote=FALSE)
    MB_output_name = paste0(resdir, "/", TF_use, "_motif", as.character(motif_dist), "_MB_LE.txt")
    write.table(MB_bed, file=MB_output_name, row.names=FALSE, sep="\t", quote=FALSE)
    message("Wrote files to: ", PE_output_name, " and ", MB_output_name)
    }

get_LE_calls(resdir=results_dir, ranked_file_use=ranked_file, TF_use=TF, motif_dist_use=motif_dist)