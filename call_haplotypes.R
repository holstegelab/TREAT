#!/usr/bin/R

#########################################
# Script to analyze TRF output for both #
# single-reads and assembly and come up #
# with a way of suggesting which is the #
# most likely haplotyping structure.    #
#########################################

# Libraries
    list.of.packages <- c("plyr", "data.table", "argparse", "stringr", "parallel", "spgs")
    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
    if(length(new.packages)) install.packages(new.packages)
    library(plyr)
    library(data.table)
    library(argparse)
    library(parallel)
    library(stringr)
    library(spgs)

# Functions
    # Function to make permutations of letters given a word, used to check motifs
    permutMotif = function(motif){
        # manage main motif
        all_perm = c(motif); motif_spl = unlist(strsplit(motif, ""))
        for (i in 1:length(motif_spl)){ 
            if (i == 1){
                tmp = paste(c(motif_spl[i:length(motif_spl)]), collapse="")
                all_perm = c(all_perm, tmp)
            } else {
                tmp = paste(c(motif_spl[i:length(motif_spl)], motif_spl[1:(i-1)]), collapse="")
                all_perm = c(all_perm, tmp)
            }
        }
        # then reverse complement
        motif_spl = unlist(strsplit(spgs::reverseComplement(motif, "dna", case = "upper"), ""))
        for (i in 1:length(motif_spl)){ 
            if (i == 1){
                tmp = paste(c(motif_spl[i:length(motif_spl)]), collapse="")
                all_perm = c(all_perm, tmp)
            } else {
                tmp = paste(c(motif_spl[i:length(motif_spl)], motif_spl[1:(i-1)]), collapse="")
                all_perm = c(all_perm, tmp)
            }
        }
        all_perm = all_perm[!duplicated(all_perm)]
        return(all_perm)  
    }

    # Function to perform guided haplotyping using phasing information using the sizes
    PhasingBased_haplotyping_size <- function(reads_df, sample_name, thr_mad){
        # split cn with haplotypes from those without
        phased = reads_df[!is.na(reads_df$HAPLOTYPE),]; nonPhased = reads_df[is.na(reads_df$HAPLOTYPE),]
        # if we're dealing with reference (only 1 haplotype), treat it differently -- the final df is all_res
        if (sample_name == 'reference'){ 
            reads_df$HAPLOTYPE = 1; reads_df$type = 'Assigned'; all_res = reads_df
        } else {
            # first check if there are non-phased reads
            if (nrow(nonPhased) >0){
                # if so, check if all reads are not phased
                if (nrow(phased) == 0){
                    # if so, we run the kmeans-based haplotyping
                    res = KMeansBased_haplotyping(reads = nonPhased$LENGTH_SEQUENCE, thr = 2, min_support = 2, thr_mad_orig = thr_mad, type = 'single_sample')
                    # then polish with the polisher without phasing
                    res_polished = polishHaplo_noPhasing(res, thr_mad)
                    # extract polished reads
                    reads_h1 = unlist(strsplit(as.character(res_polished$reads_h1), ',')); reads_h2 = unlist(strsplit(as.character(res_polished$reads_h2), ','))
                    reads_df$type = NA
                    reads_df$HAPLOTYPE = as.character(reads_df$HAPLOTYPE)
                    # make the final df with reads and haplotypes (depending on how many alleles were found)
                    if (reads_h2 == 'homozygous' && unique(!is.na(reads_h2))){
                        reads_df$type[which(reads_df$LENGTH_SEQUENCE %in% reads_h1)] = 'Assigned'
                        reads_df$HAPLOTYPE[which(reads_df$LENGTH_SEQUENCE %in% reads_h1)] = 1
                    } else if (unique(!is.na(reads_h2))){
                        reads_df$type[which(reads_df$LENGTH_SEQUENCE %in% c(reads_h1, reads_h2))] = 'Assigned'
                        reads_df$HAPLOTYPE[which(reads_df$LENGTH_SEQUENCE %in% reads_h1)] = 1; reads_df$HAPLOTYPE[which(reads_df$LENGTH_SEQUENCE %in% reads_h2)] = 2
                    }
                    # the final df is all_res
                    all_res = reads_df
                # there is phasing information available, but only haplo 1 was found
                } else if (length(unique(phased$HAPLOTYPE)) == 1) {
                    reads_df$type = NA; reads_df$type[!is.na(reads_df$HAPLOTYPE)] <- 'Phased'
                    # gather all reads within thr_mad from the haplotype
                    median_cn = median(phased$LENGTH_SEQUENCE)
                    thr_mad = max(1.5, ceiling(median(median_cn)*thr_mad))
                    # idea: look at reads without phase and try to assign them a phase based on thr_mad 
                    reads_df$type[which(abs(reads_df$LENGTH_SEQUENCE - median_cn) <= thr_mad)] <- 'Assigned'
                    reads_df$HAPLOTYPE[which(abs(reads_df$LENGTH_SEQUENCE - median_cn) <= thr_mad)] <- 1
                    # the remaining reads (those with mad > thr_mad) are haplotype 2
                    reads_df$HAPLOTYPE[is.na(reads_df$HAPLOTYPE)] = 2
                    reads_df$type[is.na(reads_df$HAPLOTYPE)] = 'Assigned'
                    # the final df is all_res
                    all_res = reads_df
                # if there is phasing information for both alleles
                } else {
                    phased$type = 'Phased'
                    # find centers (mean of copy number in the group of same haplotype)
                    h1 = mean(phased$LENGTH_SEQUENCE[which(phased$HAPLOTYPE == 1)])
                    h2 = mean(phased$LENGTH_SEQUENCE[which(phased$HAPLOTYPE == 2)])
                    # idea: assign reads without phase based on distance from those with a phase --> lowest distance gets the assignment
                    for (i in 1:nrow(nonPhased)){
                        tmp = nonPhased$LENGTH_SEQUENCE[i]; closest_haplo = ifelse(abs(tmp - h1) < abs(tmp - h2), 1, 2)
                        nonPhased$type[i] = 'Assigned'; nonPhased$HAPLOTYPE[i] = closest_haplo
                    }
                    # the final df is all_res
                    all_res = rbind(phased, nonPhased)
                }
            # if phasing information is available for all reads, we don't need to do anything
            } else {
                phased$type = 'Phased'; all_res = phased
            }
        }
        all_res$sample = sample_name
        return(all_res)
    }

    # Polish haplotypes after phasing is done
    polishHaplo_afterPhasing_size = function(phased_data, thr_mad){
        # split h1 from h2, reads and excluded
        h1 = phased_data[which(phased_data$HAPLOTYPE == 1),]; h2 = phased_data[which(phased_data$HAPLOTYPE == 2),]
        excluded = phased_data[is.na(phased_data$HAPLOTYPE),]; if (nrow(excluded) >=1) { excluded$polished_reads = 'exclude' }
        # assign haplotypes before polishing
        if (nrow(h1) >=1){ h1$haplo_value = median(h1$LENGTH_SEQUENCE) }
        if (nrow(h2) >=1){ h2$haplo_value = median(h2$LENGTH_SEQUENCE) }
        # using the thr_mad, split reads to keep from those to exclude
        h1_thr_mad = max(1.5, ceiling(median(h1$LENGTH_SEQUENCE)*thr_mad))
        h2_thr_mad = max(1.5, ceiling(median(h2$LENGTH_SEQUENCE)*thr_mad))
        h1$polished_reads = ifelse(abs(h1$LENGTH_SEQUENCE - median(h1$LENGTH_SEQUENCE)) <= h1_thr_mad, 'keep', 'exclude')
        h2$polished_reads = ifelse(abs(h2$LENGTH_SEQUENCE - median(h2$LENGTH_SEQUENCE)) <= h2_thr_mad, 'keep', 'exclude')
        # assign haplotypes after polishing
        if (nrow(h1) >= 1){ h1$polished_haplo_values = median(h1$LENGTH_SEQUENCE[which(h1$polished_reads == 'keep')]) }
        if (nrow(h2) >= 1){ h2$polished_haplo_values = median(h2$LENGTH_SEQUENCE[which(h2$polished_reads == 'keep')]) }
        # return combined object
        res = rbind.fill(h1, h2)
        # clean haplotypes -- sometimes there's haplotype 2 but not 1 --> make it 1
        if (2 %in% res$HAPLOTYPE && (!(1 %in% res$HAPLOTYPE))){ res$HAPLOTYPE = 1 }
        res = plyr::rbind.fill(res, excluded)
        return(res)
    }

    # Function to fit K-means clustering to the data to find haplotypes when no phasing information is available
    KMeansBased_haplotyping = function(reads, thr, min_support, thr_mad_orig, type){
        # define ploidy depending on the analysis type
        ploidy = ifelse(type %in% c("single_tissue", "single_sample"), 2, 3)
        # perform clustering and polishing until condition is met
        isOK = FALSE; isNA = FALSE; deleted = c()
        # main loop, condition --> isOK == TRUE
        while (isOK == FALSE){
            # first, let's check if we have enough reads
            if (length(reads) >= thr){
                # second, check whether the reads are all the same (or within acceptable distance to be called homozygous) -> homozygous
                median_value_cn = median(na.omit(reads)); reads_dist_median = abs(reads - median_value_cn)
                # check if distance is <= 1.5 or thr_mad % from the median size of the reads
                thr_mad = max(1.5, ceiling(median_value_cn * thr_mad_orig))
                if (max(reads_dist_median) <= thr_mad){
                    # if it is homozygous call, the centers of the k-means are the median values (repeted)
                    centers_kmeans = c(median_value_cn, median_value_cn)
                    # build df of reads and clusters
                    tmp_df = data.frame(reads = reads, cluster = 1); tmp_df2 = tmp_df
                    # update variables controlling the loop
                    isNA = FALSE; isOK = TRUE; mod = FALSE
                } else {
                    # check if there are exactly 2 haplotypes otherwise kmeans with k=3 will throw an error (in case k=3 is selected), otherwise this has no effect
                    if (length(unique(reads)) == 2){ ploidy = 2 }
                    # additional check here: in case we have 2 reads with non-homozygous genotype, doing k-means with k=2 will throw and error. In these cases, manually force the two alleles
                    if (length(reads) == 2){
                        centers_kmeans = reads; tmp_df = data.frame(reads = reads, cluster = c(1, 2))
                        tb_frq = as.data.frame(table(tmp_df$cluster)); tb_frq$Prop = tb_frq$Freq/sum(tb_frq$Freq); tb_frq  = tb_frq[order(tb_frq$Freq),]
                    } else {
                        # run k-means using as k the ploidy value
                        kmeans_res <- kmeans(x = reads, centers = ploidy, iter.max = 100)
                        # take centers of the k-means
                        centers_kmeans <- kmeans_res$centers[order(kmeans_res$centers)]
                        # check whether there is support for all clusters
                        tmp_df = data.frame(reads = reads, cluster = kmeans_res$cluster)
                        tb_frq = as.data.frame(table(tmp_df$cluster)); tb_frq$Prop = tb_frq$Freq/sum(tb_frq$Freq); tb_frq  = tb_frq[order(tb_frq$Freq),]
                    }
                    # if there are less than [min_support] reads per haplotype, need to remove the reads and re-do clustering
                    mod = TRUE
                    # flag haplotypes with less than [min_support] reads
                    tb_frq$support = ifelse(tb_frq$Freq < min_support, "NO", "YES")
                    # check haplotypes without enough support and in case delete the reads (add them to delete vector)
                    if ("NO" %in% tb_frq$support){
                        clusters_to_delete = tb_frq$Var1[which(tb_frq$support == "NO")]
                        deleted = c(deleted, tmp_df$reads[which(tmp_df$cluster %in% clusters_to_delete)])
                        # finally update read df for the next iteration of the loop (without excluded reads)
                        tmp_df2 = tmp_df[which(tmp_df$cluster != clusters_to_delete),]
                    } else {
                        # mod object controls for whether to re-do clustering or not
                        mod = FALSE
                    }
                    # see whether there are also outlier reads in the clusters
                    if (mod == FALSE){
                        # if we don't need to re-do the clustering, we can proceed and clean it
                        # cleaning is done based on median-absolute-deviation -- threshold is thr_mad, now in percentage depending on read size
                        tmp_df$distance_median_cluster = NA
                        # loop on every cluster and calculate the distance between reads and median value of the reads
                        for (cl in unique(tmp_df$cluster)){
                            tmp_df$distance_median_cluster[which(tmp_df$cluster == cl)] = abs(median(na.omit(tmp_df$reads[which(tmp_df$cluster == cl)])) - tmp_df$reads[which(tmp_df$cluster == cl)])
                        }
                        # check whether the median distance is greater than thr_mad, otherwise flag the read as to be removed
                        tmp_df$excluded = ifelse(tmp_df$distance_median_cluster > thr_mad, "YES", "NO")
                        # if has to be removed, add it to deleted object
                        deleted = c(deleted, tmp_df$reads[which(tmp_df$excluded == "YES")])
                        # finally exclude these outlier reads
                        tmp_df2 = tmp_df[which(tmp_df$excluded == "NO"), c("reads", "cluster")]
                    }
                    # finally check whether any modification was made to the set of reads
                    # in case of any modification, need to re-do the clustering again
                    if (nrow(tmp_df2) == nrow(tmp_df)){
                        mod = FALSE
                    } else {
                        tmp_df = tmp_df2; mod = TRUE
                    }
                }
                # update reads and modify the variables that check the re-execution of the loop
                reads = tmp_df$reads
                # always check whether we have enough reads to do stuff
                if (length(reads) >= thr & mod == TRUE){
                    isOK = FALSE; isNA = FALSE
                } else if (length(reads) >= thr & mod == FALSE){
                    isOK = TRUE; isNA = FALSE
                } else {
                    isOK = TRUE; isNA = TRUE
                }
            } else {
                isOK = TRUE; isNA = TRUE
            }
        }
        # after the loop is ended, we put things together in a df (final_df)
        if (isNA == TRUE){
            if (length(deleted) == 0){ excl = NA } else {excl = paste0(deleted, collapse = ",") }
            final_df = data.frame("reads" = paste0(reads, collapse = ","), "kmeans_centers" = NA, "reads_h1" = NA, "reads_h2" = NA, "excluded" = excl)
        } else {
            if (length(deleted) == 0){ excl = NA } else {excl = paste0(deleted, collapse = ",") }
            if (length(unique(tmp_df2$cluster)) >1){
                final_df = data.frame("reads" = paste0(reads, collapse = ","), "kmeans_centers" = paste0(centers_kmeans, collapse = ","), "reads_h1" = paste(tmp_df2$reads[which(tmp_df2$cluster == 1)], collapse = ","), "reads_h2" = paste(tmp_df2$reads[which(tmp_df2$cluster == 2)], collapse = ","), "excluded" = excl) 
            } else {
                final_df = data.frame("reads" = paste0(reads, collapse = ","), "kmeans_centers" = paste0(centers_kmeans, collapse = ","), "reads_h1" = paste(tmp_df2$reads[which(tmp_df2$cluster == 1)], collapse = ","), "reads_h2" = 'homozygous', "excluded" = excl) 
            }
        }
        return(final_df)
    }

    # Function to polish clustering of TR when no phasing information is available
    polishHaplo_noPhasing = function(data, thr_mad){
        # create columns to fille: haplotypes_polished, freq_h1 and fre1_h2
        data$haplotypes_polished = NA; data$freq_h1 = NA; data$freq_h2 = NA; 
        # convert factorized values of trf to characters
        data$reads_h1 = as.character(data$reads_h1); data$reads_h2 = as.character(data$reads_h2); data$excluded = as.character(data$excluded)
        # loop on data
        for (i in 1:nrow(data)){
            # skip if there's no kmeans-center
            if (!is.na(data$kmeans_centers[i])){
                # extract all values of reads of h1
                h1 = as.numeric(unlist(strsplit(data$reads_h1[i], ",")))
                # check if it is homozygous call
                if (data$reads_h2[i] == "homozygous"){
                    # extract also excluded reads
                    excl_reads = data$excluded[i]
                    # check if we have reads excluded or not
                    if (!is.na(excl_reads)){
                        # in case there are excluded reads, take them
                        excl_reads = as.numeric(unlist(strsplit(excl_reads, ",")))
                        new_excl = c()
                        # check if distance is <= 1.5 or thr_mad % from the median size of the reads
                        thr_mad = max(1.5, ceiling(median(excl_reads)*thr_mad))
                        # loop through excluded reads and see if they can be assigned back to some haplotype
                        for (x in 1:length(excl_reads)){
                            if (abs(median(h1) - excl_reads[x]) <= thr_mad){
                                # if the read is good, add it to h1 and re-estimate h1
                                h1 = paste(c(h1, excl_reads[x]), collapse = ",")
                                h1 = as.numeric(unlist(strsplit(h1, ",")))
                            } else {
                                new_excl = c(new_excl, excl_reads[x])
                            }
                        }
                        # at the end of the loop, re-define the excluded reads
                        new_excl = ifelse(length(new_excl) == 0, NA, paste(new_excl, collapse = ","))
                    } else {
                        new_excl = excl_reads
                    }
                    # at the end, assign the new values of h1, haplotypes, and the frequencies
                    data$reads_h1[i] = paste(h1, collapse = ","); data$excluded[i] = new_excl; 
                    data$haplotypes_polished[i] = paste0(median(h1), ",", median(h1))
                    data$freq_h1[i] = 1; data$freq_h2[i] = NA
                } else {
                    # here below in case it is NOT an homozygous call
                    h2 = as.numeric(unlist(strsplit(data$reads_h2[i], ",")))
                    # check if distance is <= 1.5 or thr_mad % from the median size of the reads
                    thr_mad = max(1.5, ceiling(median(c(h1, h2)) * thr_mad))
                    # if this is true, it is the same haplotype
                    if (abs(median(h1) - median(h2)) <= thr_mad){
                        # combine the two haplotypes
                        h1 = c(h1, h2); excl_reads = data$excluded[i]
                        # then check if excluded reads can be saved
                        if (!is.na(excl_reads)){
                            excl_reads = as.numeric(unlist(strsplit(excl_reads, ",")))
                            new_excl = c()
                            # loop through excluded reads and see if they can be assigned back to some haplotype
                            for (x in 1:length(excl_reads)){
                                if (abs(median(h1) - excl_reads[x]) <= thr_mad){
                                    # if we can save reads, add them and re-estimate reads again including the excluded
                                    h1 = paste(c(h1, excl_reads[x]), collapse = ",")
                                    h1 = as.numeric(unlist(strsplit(h1, ",")))
                                } else {
                                    new_excl = c(new_excl, excl_reads[x])
                                }
                            }
                            # at the end of the loop, re-define the excluded reads 
                            new_excl = ifelse(length(new_excl) == 0, NA, paste(new_excl, collapse = ","))
                        } else {
                            new_excl = excl_reads
                        }
                        # at the end, assign the new values of h1, haplotypes, and the frequencies
                        data$reads_h1[i] = paste(h1, collapse = ","); data$reads_h2[i] = "homozygous"; data$escluded[i] = new_excl
                        data$haplotypes_polished[i] = paste0(median(h1), ",", median(h1)); data$freq_h1[i] = 1; data$freq_h2[i] = NA
                    } else {
                        # otherwise, it is two real different haplotypes, and we check for excluded reads
                        excl_reads = data$excluded[i]
                        # then check if excluded reads can be saved
                        if (!is.na(excl_reads)){
                            excl_reads = as.numeric(unlist(strsplit(excl_reads, ",")))
                            new_excl = c()
                            # loop through excluded reads and see if they can be assigned back to some haplotype
                            for (x in 1:length(excl_reads)){
                                # check if distance is <= 1.5 or thr_mad % from the median size of the reads (h1)
                                if (abs(median(h1) - excl_reads[x]) <= thr_mad){
                                    # if we can save reads, add them and re-estimate reads again including the excluded
                                    h1 = paste(c(h1, excl_reads[x]), collapse = ",")
                                    h1 = as.numeric(unlist(strsplit(h1, ",")))
                                # check if distance is <= 1.5 or thr_mad % from the median size of the reads (h2)
                                } else if (abs(median(h2) - excl_reads[x]) <= thr_mad){
                                    # if we can save reads, add them and re-estimate reads again including the excluded
                                    h2 = paste(c(h2, excl_reads[x]), collapse = ",")
                                    h2 = as.numeric(unlist(strsplit(h2, ",")))
                                } else {
                                new_excl = c(new_excl, excl_reads[x])
                                }
                            }
                            # at the end of the loop, re-define the excluded reads
                            new_excl = ifelse(length(new_excl) == 0, NA, paste(new_excl, collapse = ","))
                        } else {
                            new_excl = excl_reads
                        }
                        # at the end, assign the new values of h1, h2, the haplotypes, and the frequencies
                        data$excluded[i] = new_excl; data$reads_h1[i] = paste(h1, collapse = ","); data$reads_h2[i] = paste(h2, collapse = ",")
                        data$haplotypes_polished[i] = paste0(median(h1), ",", median(h2))
                        data$freq_h1[i] = length(h1) / (length(h1) + length(h2)); data$freq_h2[i] = length(h2) / (length(h1) + length(h2))
                    }
                }
            }
        }
        return(data)
    }

    # Function to call haplotypes and combine information
    callHaplo_size <- function(data){
        # find all unique samples
        all_samples = unique(data$sample); all_regions = unique(data$REGION); all_haplo = list()
        # main loop across samples
        for (s in all_samples){
            for (r in all_regions){
                # gather all data together, excluded reads and good reads
                tmp = data[which(data$sample == s & data$REGION == r),]
                excl = tmp[which(tmp$polished_reads == 'exclude'),]; tmp = tmp[which(tmp$polished_reads != 'exclude'),]
                # first we managed to identify the haplotypes
                if (nrow(tmp) >0){
                    if (!(NA %in% tmp$HAPLOTYPE) & nrow(tmp) >0){
                        # if we have both haplotypes assigned, get the relative info: copies, motif, coef. of variation, length, info
                        if (1 %in% tmp$HAPLOTYPE & 2 %in% tmp$HAPLOTYPE){
                            h1_info = median(tmp$CONSENSUS_COPY_NUMBER[which(tmp$HAPLOTYPE == 1)])
                            h1_motif = paste(unique(tmp$CONSENSUS_MOTIF[which(tmp$HAPLOTYPE == 1)]), collapse = ',')
                            h1_var = sd(tmp$CONSENSUS_COPY_NUMBER[which(tmp$HAPLOTYPE == 1)]) / mean(tmp$CONSENSUS_COPY_NUMBER[which(tmp$HAPLOTYPE == 1)])
                            h1_len = paste(unique(tmp$LENGTH_SEQUENCE[which(tmp$HAPLOTYPE == 1)]), collapse = ',')
                            h1_motif_info = unique(tmp$CONSENSUS_MOTIF_COMPLEX[which(tmp$HAPLOTYPE == 1)])
                            h1_add_motifs = unique(tmp$ADDITIONAL_MOTIFS[which(tmp$HAPLOTYPE == 1)])
                            h2_info = median(tmp$CONSENSUS_COPY_NUMBER[which(tmp$HAPLOTYPE == 2)])
                            h2_motif = paste(unique(tmp$CONSENSUS_MOTIF[which(tmp$HAPLOTYPE == 2)]), collapse = ',')
                            h2_var = sd(tmp$CONSENSUS_COPY_NUMBER[which(tmp$HAPLOTYPE == 2)]) / mean(tmp$CONSENSUS_COPY_NUMBER[which(tmp$HAPLOTYPE == 2)])
                            h2_len = paste(unique(tmp$LENGTH_SEQUENCE[which(tmp$HAPLOTYPE == 2)]), collapse = ',')
                            h2_motif_info = unique(tmp$CONSENSUS_MOTIF_COMPLEX[which(tmp$HAPLOTYPE == 2)])
                            h2_add_motifs = unique(tmp$ADDITIONAL_MOTIFS[which(tmp$HAPLOTYPE == 2)])
                            data_type = unique(tmp$DATA_TYPE)
                        # otherwise if only h1 is present, only gather h1 info
                        } else if (1 %in% tmp$HAPLOTYPE){
                            h1_info = median(tmp$CONSENSUS_COPY_NUMBER[which(tmp$HAPLOTYPE == 1)])
                            h1_motif = paste(unique(tmp$CONSENSUS_MOTIF[which(tmp$HAPLOTYPE == 1)]), collapse = ',')
                            h1_var = sd(tmp$CONSENSUS_COPY_NUMBER[which(tmp$HAPLOTYPE == 1)]) / mean(tmp$CONSENSUS_COPY_NUMBER[which(tmp$HAPLOTYPE == 1)])
                            h1_len = paste(unique(tmp$LENGTH_SEQUENCE[which(tmp$HAPLOTYPE == 1)]), collapse = ',')
                            h1_motif_info = unique(tmp$CONSENSUS_MOTIF_COMPLEX[which(tmp$HAPLOTYPE == 1)])
                            h1_add_motifs = unique(tmp$ADDITIONAL_MOTIFS[which(tmp$HAPLOTYPE == 1)])
                            h2_info = NA; h2_motif = NA; h2_var = NA; h2_len = NA; h2_motif_info = NA; h2_add_motifs = NA
                            data_type = unique(tmp$DATA_TYPE)
                        # otherwise gather only h2 info
                        } else {
                            h2_info = median(tmp$CONSENSUS_COPY_NUMBER[which(tmp$HAPLOTYPE == 2)])
                            h2_motif = paste(unique(tmp$CONSENSUS_MOTIF[which(tmp$HAPLOTYPE == 2)]), collapse = ',')
                            h2_var = sd(tmp$CONSENSUS_COPY_NUMBER[which(tmp$HAPLOTYPE == 2)]) / mean(tmp$CONSENSUS_COPY_NUMBER[which(tmp$HAPLOTYPE == 2)])
                            h2_len = paste(unique(tmp$LENGTH_SEQUENCE[which(tmp$HAPLOTYPE == 2)]), collapse = ',')
                            h2_motif_info = unique(tmp$CONSENSUS_MOTIF_COMPLEX[which(tmp$HAPLOTYPE == 2)])
                            h2_add_motifs = unique(tmp$ADDITIONAL_MOTIFS[which(tmp$HAPLOTYPE == 2)])
                            h1_info = NA; h1_motif = NA; h1_var = NA; h1_len = NA; h1_motif_info = NA; h1_add_motifs = NA
                            data_type = unique(tmp$DATA_TYPE)
                        }
                        # then check excluded reads
                        if (nrow(excl) >0){ 
                            excl_info = paste(paste0(excl$CONSENSUS_MOTIF, '(', excl$CONSENSUS_COPY_NUMBER, ')'), collapse = '_')
                            excl_len = paste(excl$LENGTH_SEQUENCE, collapse = ',')
                        } else { 
                            excl_info = NA; excl_len = NA
                        }
                        # finally save a df with all info regarding TR
                        tmp_df = data.frame(SAMPLE = unique(tmp$sample), REGION = unique(tmp$REGION), H1 = h1_info, H2 = h2_info, 
                            H1_MOTIF = h1_motif, H2_MOTIF = h2_motif, H1_LENGTH = h1_len, H2_LENGTH = h2_len, COEF_VAR_H1 = h1_var, COEF_VAR_H2 = h2_var,
                            H1_MOTIF_TYPE = h1_motif_info, H2_MOTIF_TYPE = h2_motif_info, H1_ADDITIONAL_MOTIF = h1_add_motifs, H2_ADDITIONAL_MOTIF = h2_add_motifs,
                            EXCLUDED = excl_info, EXCLUDED_LEN = excl_len, DATA_TYPE = data_type)
                        # and add to the growing list
                        all_haplo[[(length(all_haplo) + 1)]] = tmp_df
                    # if haplotypes were not identified (i.e they are NAs)
                    } else {
                        excl_info = paste(paste0(excl$CONSENSUS_MOTIF, '(', excl$CONSENSUS_COPY_NUMBER, ')'), collapse = '_')
                        excl_len = paste(excl$LENGTH_SEQUENCE, collapse = ',')
                        # create the same df with NA values and add to growing list
                        tmp_df = data.frame(SAMPLE = unique(excl$sample), REGION = unique(excl$REGION), H1 = NA, H2 = NA, 
                            H1_MOTIF = NA, H2_MOTIF = NA, H1_LENGTH = NA, H2_LENGTH = NA, COEF_VAR_H1 = NA, COEF_VAR_H2 = NA,
                            H1_MOTIF_TYPE = NA, H2_MOTIF_TYPE = NA, H1_ADDITIONAL_MOTIF = NA, H2_ADDITIONAL_MOTIF = NA,
                            EXCLUDED = excl_info, EXCLUDED_LEN = excl_len, DATA_TYPE = data_type)
                        all_haplo[[(length(all_haplo) + 1)]] = tmp_df
                    }
                }
            }
        }
        # convert list to df at the end
        all_haplo = rbindlist(all_haplo)
        return(all_haplo)
    }

    # Function to call haplotypes when dealing with assembly-based TRF
    assemblyBased_size = function(data, sample_name, region, thr_mad){
        # assign variables to fill
        data$type = NA
        data$sample = sample_name
        # if we're dealing with assemblies, in theory it is simpler the haplotyping
        # we need to check whether we have enough data (>0 contigs)
        if (nrow(data) >0){
            # flag cases when there are multiple data (>2 contigs)
            if (nrow(data) >2){
                # in this case run the normal kmean function
                res = KMeansBased_haplotyping(reads = data$LENGTH_SEQUENCE, thr = 2, min_support = 1, thr_mad_orig = thr_mad, type = 'single_sample')
                # then polish with the polisher without phasing
                res_polished = polishHaplo_noPhasing(res, thr_mad)
                # extract polished reads
                reads_h1 = unlist(strsplit(as.character(res_polished$reads_h1), ',')); reads_h2 = unlist(strsplit(as.character(res_polished$reads_h2), ','))
                # make the final df with reads and haplotypes (depending on how many alleles were found)
                if (reads_h2 == 'homozygous' && !is.na(reads_h2)){
                    data$type[which(data$LENGTH_SEQUENCE %in% reads_h1)] = 'Assigned'; data$HAPLOTYPE[which(data$LENGTH_SEQUENCE %in% reads_h1)] = 1
                }
            } else {
                # then we need to check if it is an homozygous call (1 contig)
                if (nrow(data) == 1){
                    # homozygous call
                    data$HAPLOTYPE = 1
                    data$type = 'assembly'
                } else {
                    # heterozygous call
                    data$HAPLOTYPE = c(1, 2)
                    data$type = 'assembly'
                }
            }
        }
        return(data)
    }

    # Function to check consensus motif given a list of motifs -- more conscious approach
    consensusMotif_conscious = function(motifs){
        # first thing is to check whether we only have 1 single motif --> easy case
        if (length(na.omit(unique(motifs$UNIFORM_MOTIF))) == 1){
            # before assigning, let's check for NA, and if there are too many, exclude the call
            tb = data.frame(table(motifs$UNIFORM_MOTIF, exclude = F)); tb$Perc = tb$Freq / nrow(motifs)
            if (NA %in% tb$Var1 && tb$Perc[is.na(tb$Var1)] > 0.8){
                motifs$CONSENSUS_MOTIF = NA
                motifs$CONSENSUS_COPY_NUMBER = NA
                motifs$CONSENSUS_MOTIF_COMPLEX = NA
                motifs$ADDITIONAL_MOTIFS = NA
            } else {
                motifs$CONSENSUS_MOTIF = unique(na.omit(motifs$UNIFORM_MOTIF))
                motifs$CONSENSUS_COPY_NUMBER = motifs$COPIES_TRF
                motifs$CONSENSUS_MOTIF_COMPLEX = 'simple'
                motifs$ADDITIONAL_MOTIFS = NA
            }
        } else if (length(na.omit(unique(motifs$UNIFORM_MOTIF))) == 0){
            motifs$CONSENSUS_MOTIF = NA
            motifs$CONSENSUS_COPY_NUMBER = NA
            motifs$CONSENSUS_MOTIF_COMPLEX = NA
            motifs$ADDITIONAL_MOTIFS = NA
        } else {
            # calculate frequencies of each motif
            motifs_frq = data.frame(table(as.character(motifs$UNIFORM_MOTIF))); motifs_frq$Perc = motifs_frq$Freq / nrow(motifs); motifs_frq$motif_length = nchar(as.character(motifs_frq$Var1))
            # keep motifs with highest frequency
            motifs_frq = motifs_frq[order(-motifs_frq$Perc, -motifs_frq$motif_length),]
            consensus_motif = as.character(motifs_frq$Var1[1])
            # raise a flag when the frequency of most common motif is <0.8
            if (motifs_frq$Perc[1] < 0.8){ 
                # in this case: either the motif is really variable across reads, or there is a clear pattern, such as each read has two motifs that suggests these should be combined
                # then it depends on the number of duplicated reads, if this approximate 50%, we should further look, otherwise not
                dups = unique(motifs$READ_NAME[duplicated(motifs$READ_NAME)])
                if (length(dups)/nrow(motifs) >= 0.4){
                    # if there are many duplicates, there are 2 options: either there are multiple motifs that combine together (nested repeats), or there are two alternative motifs that can be estimated
                    # to distinguish the two cases, we can look at the percentage of coverage
                    tmp = motifs; tmp$coverage = (tmp$END_TRF - tmp$START_TRF) / tmp$LENGTH_SEQUENCE
                    # calculate also difference wrt 1
                    tmp$difference = abs(1 - tmp$coverage)
                    # if the coverage is very high, that means we are dealing with alternative motifs for the same TR: take the motif with highest coverage
                    if (TRUE %in% (tmp$coverage >0.90)){
                        #consensus_motif = unique(tmp$UNIFORM_MOTIF[which(tmp$coverage == max(tmp$coverage))])[1]
                        consensus_motif = unique(tmp$UNIFORM_MOTIF[which(tmp$coverage == min(tmp$coverage))])[1]
                        motifs$CONSENSUS_MOTIF = consensus_motif
                        motifs$CONSENSUS_COPY_NUMBER = NA; motifs$CONSENSUS_COPY_NUMBER[which(motifs$UNIFORM_MOTIF == consensus_motif)] = motifs$COPIES_TRF[which(motifs$UNIFORM_MOTIF == consensus_motif)]
                        motifs$CONSENSUS_COPY_NUMBER[is.na(motifs$CONSENSUS_COPY_NUMBER)] = median(motifs$CONSENSUS_COPY_NUMBER[which(motifs$UNIFORM_MOTIF == consensus_motif)])
                        motifs$CONSENSUS_MOTIF_COMPLEX = 'alternative'; motifs$ADDITIONAL_MOTIFS = paste(unique(motifs$UNIFORM_MOTIF[which(motifs$UNIFORM_MOTIF != consensus_motif)]), collapse = ',')
                        motifs = motifs[!duplicated(motifs$READ_NAME),]
                    } else {
                        # what to do in these cases? there are multiple TRF matches, with similar percentage not above 0.9 coverage. This suggests they should be merged
                        # check how many matches are left: if it is only 1, then we're done, otherwise probably we need to merge the motifs
                        if (nrow(tmp) == 1){
                            consensus_motif = unique(tmp$UNIFORM_MOTIF)
                            motifs$CONSENSUS_MOTIF = consensus_motif
                            motifs$CONSENSUS_COPY_NUMBER = NA; motifs$CONSENSUS_COPY_NUMBER[which(motifs$UNIFORM_MOTIF == consensus_motif)] = motifs$COPIES_TRF[which(motifs$UNIFORM_MOTIF == consensus_motif)]
                            motifs$CONSENSUS_COPY_NUMBER[is.na(motifs$CONSENSUS_COPY_NUMBER)] = median(motifs$CONSENSUS_COPY_NUMBER[which(motifs$UNIFORM_MOTIF == consensus_motif)])
                            motifs$CONSENSUS_MOTIF_COMPLEX = 'simple'; motifs$ADDITIONAL_MOTIFS = NA
                            motifs = motifs[!duplicated(motifs$READ_NAME),]                            
                        } else {
                            # these cases can be multiple motifs that need to be merged together, or a missing call
                            # check whether if summing the motifs coverage we reach a high value
                            if (sum(tmp$coverage) > 0.8){
                                # in this case it suggests that we should merge the motifs together. loop on each unique read
                                for (read in unique(tmp$READ_NAME)){
                                    tmp_duplicated_reads = tmp[which(tmp$READ_NAME == read),]
                                    # measure the overlap between the duplicates
                                    overlaps = list(); for (i in 1:nrow(tmp_duplicated_reads)){ overlaps[[(length(overlaps) + 1)]] = seq(tmp_duplicated_reads$START_TRF[i], tmp_duplicated_reads$END_TRF[i]) }
                                    # then calculate joined coverage and overlaps
                                    join_coverage = unlist(overlaps); overlaps = join_coverage[duplicated(join_coverage)]
                                    # calculate cumulative coverage (20 is the estimated padding)
                                    cum_coverage = (length(join_coverage) - length(overlaps)) / tmp_duplicated_reads$LENGTH_SEQUENCE[1]
                                    # check if there was a gain in doing this in terms of cumulative coverage compared to the original coverages
                                    if (max(cum_coverage, tmp_duplicated_reads$PERC_SEQ_COVERED_TR) == cum_coverage){
                                        # if there is a gain, probably we're dealing with a true variable motif
                                        merged_motif = paste(tmp_duplicated_reads$UNIFORM_MOTIF, collapse = '/'); motifs$CONSENSUS_MOTIF[which(motifs$READ_NAME == read)] = merged_motif
                                        motifs$CONSENSUS_COPY_NUMBER[which(motifs$READ_NAME == read)] = tmp_duplicated_reads$LENGTH_SEQUENCE[1]/(nchar(merged_motif)-1); motifs$CONSENSUS_MOTIF_COMPLEX = 'nested_motif'; motifs$ADDITIONAL_MOTIFS = paste(tmp$UNIFORM_MOTIF, collapse = ',')
                                    } else {
                                        # if there is no gain, probably we're dealing with a probable alternative motif
                                        motifs$CONSENSUS_MOTIF = NA
                                        motifs$CONSENSUS_COPY_NUMBER = NA; motifs$CONSENSUS_MOTIF_COMPLEX = 'motif_low_coverage'; motifs$ADDITIONAL_MOTIFS = paste(tmp$UNIFORM_MOTIF, collapse = ',')
                                    }
                                }
                            } else {
                                motifs$CONSENSUS_MOTIF = NA
                                motifs$CONSENSUS_COPY_NUMBER = NA; motifs$CONSENSUS_MOTIF_COMPLEX = 'motif_low_coverage'; motifs$ADDITIONAL_MOTIFS = paste(tmp$UNIFORM_MOTIF, collapse = ',')
                            }
                        }
                    }
                } else {
                    # otherwise remove the wrong motif and keep the consensus motif from reads with duplicated motifs
                    wrong_dup = motifs[which(motifs$UNIFORM_MOTIF != consensus_motif),]
                    for (i in 1:nrow(wrong_dup)){
                        good_dup = motifs[which(motifs$READ_NAME == wrong_dup$READ_NAME[i] & motifs$UNIFORM_MOTIF == consensus_motif),]
                        if (nrow(good_dup) >0){ motifs = motifs[-which(motifs$READ_NAME == wrong_dup$READ_NAME[i] & motifs$UNIFORM_MOTIF == wrong_dup$UNIFORM_MOTIF[i])] }
                    }
                    # the remaining reads: those with different motifs and no duplicate with consensus motif, will be forced to be with the consensus motif
                    # at the end we need to re-estimate copy-number when the motif is different and the read with the wrong motif has no duplicated with the consensus motif
                    motifs$CONSENSUS_MOTIF = consensus_motif
                    motifs$CONSENSUS_COPY_NUMBER = NA; motifs$CONSENSUS_COPY_NUMBER[which(motifs$UNIFORM_MOTIF == consensus_motif)] = motifs$COPIES_TRF[which(motifs$UNIFORM_MOTIF == consensus_motif)]
                    motifs$CONSENSUS_COPY_NUMBER[is.na(motifs$CONSENSUS_COPY_NUMBER)] = median(motifs$CONSENSUS_COPY_NUMBER[which(motifs$UNIFORM_MOTIF == consensus_motif)])
                    motifs$CONSENSUS_MOTIF_COMPLEX = 'variable'; motifs$ADDITIONAL_MOTIFS = paste(unique(motifs$UNIFORM_MOTIF[which(motifs$UNIFORM_MOTIF != consensus_motif)]), collapse = ',')
                }
            } else {
                # most times TRF finds two motifs for the same read, one of which is the consensus and the other is not
                # if the consensun motif is largely the most frequent, just delete the other motif for the duplicated read
                wrong_dup = motifs[which(motifs$UNIFORM_MOTIF != consensus_motif),]
                for (i in 1:nrow(wrong_dup)){
                    good_dup = motifs[which(motifs$READ_NAME == wrong_dup$READ_NAME[i] & motifs$UNIFORM_MOTIF == consensus_motif),]
                    if (nrow(good_dup) >0){ motifs = motifs[-which(motifs$READ_NAME == wrong_dup$READ_NAME[i] & motifs$UNIFORM_MOTIF == wrong_dup$UNIFORM_MOTIF[i])] }
                }
                # at the end we need to re-estimate copy-number when the motif is different and the read with the wrong motif has no duplicated with the consensus motif
                motifs$CONSENSUS_MOTIF = consensus_motif
                motifs$CONSENSUS_COPY_NUMBER = NA; motifs$CONSENSUS_COPY_NUMBER[which(motifs$UNIFORM_MOTIF == consensus_motif)] = motifs$COPIES_TRF[which(motifs$UNIFORM_MOTIF == consensus_motif)]
                motifs$CONSENSUS_COPY_NUMBER[is.na(motifs$CONSENSUS_COPY_NUMBER)] = median(motifs$CONSENSUS_COPY_NUMBER[which(motifs$UNIFORM_MOTIF == consensus_motif)])
                motifs$CONSENSUS_MOTIF_COMPLEX = 'simple'; motifs$ADDITIONAL_MOTIFS = NA
            }
        }
        return(motifs)
    }

    # Function to compare reads-spanning and assembly-based approached
    compareReadsSpanning_Asm = function(all_haplo, deviation){
        all_haplo$SAME_NAME = str_split_fixed(all_haplo$SAMPLE, '__', 2)[, 1]
        # exclude reference
        reference = all_haplo[which(all_haplo$SAME_NAME == 'reference'),]; all_haplo = all_haplo[which(all_haplo$SAME_NAME != 'reference'),]
        # loop on each sample
        all_samples = unique(all_haplo$SAME_NAME); all_regions = unique(all_haplo$REGION); comparison = list()
        for (s in all_samples){
            for (r in all_regions){
                tmp = all_haplo[which(all_haplo$SAME_NAME == s & all_haplo$REGION == r),]
                # split assembly and reads-spanning
                asm = tmp[which(tmp$DATA_TYPE == 'assembly'),]; spa = tmp[which(tmp$DATA_TYPE == 'reads-spanning'),]
                # for the assembly, keep the haplotype-aware contigs if there are multiple entries
                if (nrow(asm) >1){ asm = asm[grep('Assembly_haps', asm$SAMPLE),] }
                # first check if we have asm and reads-spanning
                if (nrow(asm) >0 & nrow(spa)){
                    # pairwise comparison of haplotypes. first extract data for asm and spa for h1 and h2
                    asm_h1 = asm[, c('H1', 'H1_MOTIF', 'H1_LENGTH')]; spa_h1 = spa[, c('H1', 'H1_MOTIF', 'H1_LENGTH')]
                    asm_h2 = asm[, c('H2', 'H2_MOTIF', 'H2_LENGTH')]; spa_h2 = spa[, c('H2', 'H2_MOTIF', 'H2_LENGTH')]
                    # then check for NAs
                    if (NA %in% c(asm_h1$H1, asm_h2$H2, spa_h1$H1, spa_h2$H2)){
                        # find the NA(s)
                        na_index = which(is.na(c(asm_h1$H1, asm_h2$H2, spa_h1$H1, spa_h2$H2)))
                        if (length(na_index) == 4){
                            colnames(asm) = paste0(colnames(asm), '_ASM'); asm$COMPARISON = 'missing_all_alleles'
                            colnames(spa) = paste0(colnames(spa), '_READS_SPANNING'); spa$H1_DIFF = NA; spa$H2_DIFF = NA
                        } else if (length(na_index) == 3){
                            good_one = which(!seq(1, 4) %in% na_index)
                            colnames(asm) = paste0(colnames(asm), '_ASM')
                            colnames(spa) = paste0(colnames(spa), '_READS_SPANNING'); spa$H1_DIFF = diff_tmp; spa$H2_DIFF = NA
                            if (good_one %in% c(1, 2)){ asm$COMPARISON = 'only_1_allele_from_assembly' } else { asm$COMPARISON = 'only_1_allele_from_reads_spanning' }
                        } else if (length(na_index) == 2){
                            if (na_index[1] %in% c(1, 2) & na_index[2] %in% c(1, 2)){
                                # in this case we are missing both asm alleles
                                colnames(asm) = paste0(colnames(asm), '_ASM'); asm$COMPARISON = 'missing_both_assembly_alleles'
                                colnames(spa) = paste0(colnames(spa), '_READS_SPANNING'); spa$H1_DIFF = NA; spa$H2_DIFF = NA
                            } else if (na_index[1] %in% c(3, 4) & na_index[2] %in% c(3, 4)){
                                # in this case we are missing both reads-spanning alleles
                                colnames(spa) = paste0(colnames(spa), '_READS_SPANNING')
                                colnames(asm) = paste0(colnames(asm), '_ASM'); asm$COMPARISON = 'missing_both_reads_spanning_alleles'
                                spa$H1_DIFF = NA; spa$H2_DIFF = NA                             
                            } else {
                                # in this case there's one NA in both asm and reads-spanning -- we need to check the other allele
                                # manage asm
                                if (1 %in% na_index){ nonmissing_asm = asm_h2 } else if (2 %in% na_index){ nonmissing_asm = asm_h1 }
                                # manage spanning
                                if (3 %in% na_index){ nonmissing_spa = spa_h2 } else if (4 %in% na_index){ nonmissing_spa = spa_h1 }
                                # check nonmissing allele
                                diff_tmp = abs(as.numeric(nonmissing_asm[, 1]) - as.numeric(nonmissing_spa[, 1]))
                                if (diff_tmp <= max(1.5, as.numeric(nonmissing_asm[, 1])*deviation)){ check = 'both_miss_1_allele_other_ok' } else { check = 'both_miss_1_allele_other_different' }
                                colnames(asm) = paste0(colnames(asm), '_ASM'); asm$COMPARISON = check
                                colnames(spa) = paste0(colnames(spa), '_READS_SPANNING'); spa$H1_DIFF = diff_tmp; spa$H2_DIFF = NA
                            }
                        } else {
                            if (na_index %in% c(1, 2)){
                                if (na_index == 1){ nonmissing = asm_h2 } else { nonmissing = asm_h1 }
                                # calculate distance between alleles based on copy number
                                diff_alleles = abs(c(as.numeric(nonmissing[, 1]) - spa_h1$H1, as.numeric(nonmissing[, 1]) - spa_h2$H2))
                                if (length(unique(diff_alleles)) == 1){
                                    colnames(asm) = paste0(colnames(asm), '_ASM'); asm$COMPARISON = 'assembly_miss_1_but_homozygous'
                                    colnames(spa) = paste0(colnames(spa), '_READS_SPANNING'); spa$H1_DIFF = diff_alleles[1]; spa$H2_DIFF = NA
                                } else {
                                    closest = which(diff_alleles == min(diff_alleles)); if (closest == 1){ closest_all = spa_h1 } else { closest_all = spa_h2 }
                                    if (diff_alleles[closest] <= max(1.5, as.numeric(closest_all[,1])*deviation)){
                                        colnames(asm) = paste0(colnames(asm), '_ASM'); asm$COMPARISON = 'assembly_miss_1_allele_other_ok'
                                        colnames(spa) = paste0(colnames(spa), '_READS_SPANNING'); spa$H1_DIFF = diff_alleles[closest]; spa$H2_DIFF = NA
                                    } else {
                                        colnames(asm) = paste0(colnames(asm), '_ASM'); asm$COMPARISON = 'assembly_miss_1_allele_other_different'
                                        colnames(spa) = paste0(colnames(spa), '_READS_SPANNING'); spa$H1_DIFF = diff_alleles[closest]; spa$H2_DIFF = NA
                                    }
                                }
                            } else {
                                if (na_index == 3){ nonmissing = spa_h2 } else { nonmissing = spa_h1 }
                                # calculate distance between alleles based on copy number
                                diff_alleles = abs(c(as.numeric(nonmissing[, 1]) - asm_h1$H1, as.numeric(nonmissing[, 1]) - asm_h2$H2))
                                if (length(unique(diff_alleles)) == 1){
                                        colnames(asm) = paste0(colnames(asm), '_ASM'); asm$COMPARISON = 'assembly_miss_1_but_homozygous'
                                        colnames(spa) = paste0(colnames(spa), '_READS_SPANNING'); spa$H1_DIFF = diff_alleles[1]; spa$H2_DIFF = NA
                                } else {
                                    closest = which(diff_alleles == min(diff_alleles)); if (closest == 1){ closest_all = asm_h1 } else { closest_all = asm_h2 }
                                    if (diff_alleles[closest] <= max(1.5, as.numeric(closest_all[,1])*deviation)){
                                        colnames(asm) = paste0(colnames(asm), '_ASM'); asm$COMPARISON = 'reads_spanning_miss_1_allele_other_ok'
                                        colnames(spa) = paste0(colnames(spa), '_READS_SPANNING'); spa$H1_DIFF = diff_alleles[closest]; spa$H2_DIFF = NA
                                    } else {
                                        colnames(asm) = paste0(colnames(asm), '_ASM'); asm$COMPARISON = 'reads_spanning_miss_1_allele_other_different'
                                        colnames(spa) = paste0(colnames(spa), '_READS_SPANNING'); spa$H1_DIFF = diff_alleles[closest]; spa$H2_DIFF = NA
                                    }
                                }
                            }
                        }
                    } else {
                        # calculate distance between alleles: h1
                        diff_h1 = abs(c(asm_h1$H1 - spa_h1$H1, asm_h1$H1 - spa_h2$H2))
                        closest_h1 = which(diff_h1 == min(diff_h1))[1]; if (closest_h1 == 1){ closest_all = spa_h1 } else { closest_all = spa_h2 }
                        if (diff_h1[closest_h1] <= max(1.5, as.numeric(closest_all[,1])*deviation)){ check_h1 = 'same' } else { check_h1 = 'different' }
                        # then h2, exclude the h1 allele
                        if (colnames(closest_all)[1] == 'H2'){ diff_h2 = abs(asm_h2$H2 - spa_h1$H1) } else { diff_h2 = abs(asm_h2$H2 - spa_h2$H2) }
                        closest_h2 = abs(closest_h1 - 1); if (closest_h2 == 1){ closest_all = spa_h1 } else { closest_all = spa_h2 }
                        if (diff_h2 <= max(1.5, as.numeric(closest_all[,1])*deviation)){ check_h2 = 'same' } else { check_h2 = 'different' }
                        if (check_h1 == 'same' & check_h2 == 'same'){ check = 'same_alleles' } else { check = 'alleles_do_not_overlap' }
                        colnames(asm) = paste0(colnames(asm), '_ASM'); asm$COMPARISON = check; colnames(spa) = paste0(colnames(spa), '_READS_SPANNING'); spa$H1_DIFF = diff_h1[closest_h1]; spa$H2_DIFF = diff_h2
                    }
                } else if (nrow(asm) >0){
                    colnames(asm) = paste0(colnames(asm), '_ASM'); asm$COMPARISON = 'missing_reads_spanning'                   
                    colnames(spa) = paste0(colnames(spa), '_READS_SPANNING'); spa$H1_DIFF = NA; spa$H2_DIFF = NA
                } else if (nrow(spa) >0){
                    colnames(spa) = paste0(colnames(spa), '_READS_SPANNING'); spa$COMPARISON = 'missing_assembly'
                    colnames(asm) = paste0(colnames(asm), '_ASM'); spa$H1_DIFF = NA; spa$H2_DIFF = NA
                } 
                tmp_combined_res = cbind(asm, spa)
                if (nrow(tmp_combined_res) >0){ comparison[[(length(comparison) + 1)]] = tmp_combined_res }
            }
        }
        comparison = rbindlist(comparison, use.names = T)
    }

    # Function to make final decision between assembly and reads-spanning
    makeDecision = function(all_haplo_annotated){
        all_haplo_annotated$DECISION = NA        
        # main loop to make decisions
        for (i in 1:nrow(all_haplo_annotated)){
            if (all_haplo_annotated$COMPARISON[i] %in% c('assembly_miss_1_but_homozygous', 'same_alleles')){
                all_haplo_annotated$DECISION[i] = 'both_fine'
            } else if (all_haplo_annotated$COMPARISON[i] %in% c('reads_spanning_miss_1_allele_other_ok', 'missing_reads_spanning', 'missing_both_reads_spanning_alleles', 'only_1_allele_from_assembly', 'reads_spanning_miss_1_allele_other_different')){
                all_haplo_annotated$DECISION[i] = 'better_assembly'
            } else if (all_haplo_annotated$COMPARISON[i] %in% c('assembly_miss_1_allele_other_ok', 'missing_assembly', 'only_1_allele_from_reads_spanning', 'assembly_miss_1_allele_other_different')){
                all_haplo_annotated$DECISION[i] = 'better_reads_spanning'
            } else if (all_haplo_annotated$COMPARISON[i] %in% c('alleles_do_not_overlap', 'both_miss_1_allele_other_different')){
                all_haplo_annotated$DECISION[i] = 'different_alleles'
            } else if (all_haplo_annotated$COMPARISON[i] %in% c('both_miss_1_allele_other_ok')){
                all_haplo_annotated$DECISION[i] = 'missing_1_allele'
            } else if (all_haplo_annotated$COMPARISON[i] %in% c('missing_all_alleles')){
                all_haplo_annotated$DECISION[i] = 'missing_all_alleles'
            } else {
                cat(paste0('**** !!! problem for ', i, '\n'))
            }
        }
        return(all_haplo_annotated)
    }

    # Function to merge motifs based on multiple processing
    mergeMotif_mp = function(motif){
        all_perm = permutMotif(as.character(motif))
        all_perm = all_perm[order(all_perm)]
        return(all_perm[1])
    }

    # Function to do haplotyping based on multiple processing
    haplotyping_mp = function(s, reads_span, all_regions, type, thr_mad){
        print(paste0('** processing sample --> ', s))
        # initialize dataframe for results
        tmp_res = data.frame()
        for (r in all_regions){
            # get data of the sample and the region of interest -- depending on type (reads_spanning or asm)
            tmp_data = reads_span[which(reads_span$SAMPLE_NAME == s & reads_span$REGION == r),]
            if (nrow(tmp_data) >0){
                reads_df = tmp_data[, c('COPIES_TRF', 'HAPLOTYPE', 'UNIFORM_MOTIF', 'REGION', 'PASSES', 'READ_QUALITY', 'LENGTH_SEQUENCE', 'READ_NAME', 'START_TRF', 'END_TRF', 'TRF_PERC_MATCH', 'TRF_PERC_INDEL')]
                if (type == 'reads_spanning'){
                    phased_data = PhasingBased_haplotyping_size(reads_df, sample_name = s, thr_mad)
                } else {
                    phased_data = assemblyBased_size(reads_df, sample_name = s, region = r, thr_mad)
                }
                polished_data = polishHaplo_afterPhasing_size(phased_data, thr_mad)
                tmp_res = rbind.fill(tmp_res, polished_data)
            }
        }
        return(tmp_res)
    }

    # Function to generate consensus motif using the majority rule consensus
    motif_generalization = function(h1){
        # idea is to calculate a consensus motif, then iterate on the reads and generate the consensus representation
        # to generate the consensus motif, list all possible motifs
        all_motifs = as.character(unique(h1$UNIFORM_MOTIF))
        motifs_to_use = c()
        # in general there can be two scenarios: 
        # 1. different trf matches cover different part of the sequence --> motifs should be combined
        # 2. different trf matches cover the same sequence --> motifs are variable and we should report 1 motif only
        # before going to the consensus motif generation, we should identify the scenario we are in
        # to do that, first, calculate the fraction of sequence covered by each read
        h1$COVERAGE_TR = (h1$END_TRF - h1$START_TRF + 1) / h1$polished_haplo_values
        # next, take the motif that covers the most of the TR
        best_motif = as.character(h1$UNIFORM_MOTIF[order(-h1$COVERAGE_TR)][1])
        motifs_to_use = best_motif
        # define a sequence of indexes representing the coverage of the best motif based on TRF_START and TRF_END
        best_motif_seq = seq(h1$START_TRF[order(-h1$COVERAGE_TR)][1], h1$END_TRF[order(-h1$COVERAGE_TR)][1])
        # then see if adding other motifs improves the sequence coverage
        alt_motifs = h1[which(h1$UNIFORM_MOTIF != best_motif),]
        if (nrow(alt_motifs) >=1){
            for (i in 1:nrow(alt_motifs)){
                # define a sequence of indexes representing the coverage of the best motif based on TRF_START and TRF_END
                alt_motif_seq = seq(alt_motifs$START_TRF[i], alt_motifs$END_TRF[i])
                # join the two list of indexes and remove duplicates
                combined_seq = c(best_motif_seq, alt_motif_seq)
                combined_seq = combined_seq[!duplicated(combined_seq)]
                # calculate coverage of this new combined sequence
                combined_seq_coverage = length(combined_seq) / h1$polished_haplo_values
                # check if the coverage of the combined sequence is higher than the single sequence
                if (combined_seq_coverage[i] > (max(h1$COVERAGE_TR) + 0.05)){
                    # if so, this is a motif we want to use
                    motifs_to_use = c(motifs_to_use, as.character(alt_motifs$UNIFORM_MOTIF[i]))
                }
            }
        }
        motifs_to_use = unique(motifs_to_use)
        # then go ahead with merging the motifs only if we have more than 1 motif
        if (length(motifs_to_use) >1){
            # for generating the consensus motif, we will use the "majority rule consensus"
            # first thing is to convert the characters into a list of n elements (as many as the sequence length)
            all_motifs_chars = lapply(X = motifs_to_use, FUN = function(x){strsplit(x, "")[[1]]})
            # find the maximum length of all motifs in the sequence
            max_len_motif <- max(sapply(all_motifs_chars, length))
            # make sure to have all motifs with the same size by adding "x"
            all_motifs_chars <- lapply(all_motifs_chars, function(y) { n_x <- max_len_motif - length(y); c(y, rep("x", n_x)) })
            # initialize the consensus motif variable
            consensus_motif = c()
            # iterate on the indexes of the motif
            for (pos in seq(1, max_len_motif)){
                # get the values of each motif at a given position
                pos_values <- unlist(lapply(all_motifs_chars, "[[", pos))
                # check if the values are all the same
                pos_value_consensus = ifelse(test = length(unique(pos_values)) == 1, yes = unique(pos_values), no = 'x')
                # then add to the consensus motif
                consensus_motif = c(consensus_motif, pos_value_consensus)
            }
            # finally, combine motif
            consensus_motif = paste(consensus_motif, collapse = "")
            # also combine the different alternative motifs
            motifs_to_use = motifs_to_use[order(motifs_to_use)]
            alternative_motifs = paste(motifs_to_use, collapse = ',')
            
            # In addition to the consensus motif, we can also find the specific motif
            # depending on the number of motifs that we find, we can use a clustering approach
            # for example, if there were 2 motifs found, that means that there are 2 defined TRF matches that (should) start at the same position
            h1 = h1[order(h1$START_TRF),]
            # fit k-means using the data sorted by start position of TRF, and the number of motifs is the k-value
            fit = kmeans(x = h1[, c('START_TRF', 'END_TRF')], centers = length(motifs_to_use), algorithm = 'Lloyd')
            specific_representation = c()
            # iterate over the clusters to determine the specific representation
            for (i in unique(fit$cluster)){
                # isolate matches on the cluster
                tmp_reads = h1[which(fit$cluster == i),]
                # there can be that there are still multiple motifs, but we want to select the most abundant one
                most_frequent <- names(which.max(table(tmp_reads$UNIFORM_MOTIF)))
                # the number of copies is the median value
                copy_n_median = median(tmp_reads$COPIES_TRF)
                # compose the specific motif
                specific_representation = paste0(specific_representation, '(', most_frequent, ')', copy_n_median)
            }
        } else {
            # in this case, we only have 1 motif, so we are OK
            consensus_motif = motifs_to_use
            alternative_motifs = NA
            specific_representation = motifs_to_use
        }
        # add combined motif to data and estimate number of copies
        h1$CONSENSUS_MOTIF = consensus_motif
        h1$CONSENSUS_MOTIF_EST_COPIES = h1$polished_haplo_values / nchar(consensus_motif)
        h1$CONSENSUS_MOTIF_ALT_MOTIFS = alternative_motifs
        h1$SPECIFIC_MOTIF = specific_representation
        # return the same object we used as input with additional columns
        return(h1)
    }
    
    # Function to generate consensus motif for each sample and region based on multiple processing
    generateConsens_mp = function(s, all_regions, all_res, motif_res_reference){
        motif_res = list()
        for (r in all_regions){
            # take all reads
            tmp = all_res[which(all_res$sample == s & all_res$REGION == r),]
            # check number of haplotypes
            n_haplo = unique(tmp$HAPLOTYPE)
            # treat haplotypes independently from each other
            if (1 %in% tmp$HAPLOTYPE & 2 %in% tmp$HAPLOTYPE){
                h1 = tmp[which(tmp$HAPLOTYPE == 1),]
                h2 = tmp[which(tmp$HAPLOTYPE == 2),]
                #h1_consensus_motif = consensusMotif_conscious(h1); h2_consensus_motif = consensusMotif_conscious(h2)
                # calculate consensus motif using the majority rule consensus
                h1_consensus_motif = motif_generalization(h1)
                h2_consensus_motif = motif_generalization(h2)
                # calculate copy number estimation wrt reference genome -- uniform for associations altough not exactly correct
                reference_motif = unique(motif_res_reference$CONSENSUS_MOTIF[which(motif_res_reference$REGION == r)])
                h1_consensus_motif$CONSENSUS_MOTIF_REF = reference_motif; h2_consensus_motif$CONSENSUS_MOTIF_REF = reference_motif
                h1_consensus_motif$CONSENSUS_MOTIF_REF_EST_COPIES = h1_consensus_motif$polished_haplo_values / nchar(reference_motif); h2_consensus_motif$CONSENSUS_MOTIF_REF_EST_COPIES = h2_consensus_motif$polished_haplo_values / nchar(reference_motif)
                tmp_res = rbind(h1_consensus_motif, h2_consensus_motif)
            } else if (1 %in% tmp$HAPLOTYPE){
                h1 = tmp[which(tmp$HAPLOTYPE == 1),]
                #tmp_res = consensusMotif_conscious(h1)
                # calculate consensus motif using the majority rule consensus
                tmp_res = motif_generalization(h1)
                # calculate copy number estimation wrt reference genome -- uniform for associations altough not exactly correct
                if (s == 'reference'){
                    tmp_res$CONSENSUS_MOTIF_REF = NA
                    tmp_res$CONSENSUS_MOTIF_REF_EST_COPIES = NA
                } else {
                    reference_motif = unique(motif_res_reference$CONSENSUS_MOTIF[which(motif_res_reference$REGION == r)])
                    tmp_res$CONSENSUS_MOTIF_REF = reference_motif
                    tmp_res$CONSENSUS_MOTIF_REF_EST_COPIES = tmp_res$polished_haplo_values / nchar(reference_motif)
                }
            } else if (2 %in% tmp$HAPLOTYPE){
                h2 = tmp[which(tmp$HAPLOTYPE == 2),]
                #tmp_res = consensusMotif_conscious(h2)
                # calculate consensus motif using the majority rule consensus
                tmp_res = motif_generalization(h2)
                # calculate copy number estimation wrt reference genome -- uniform for associations altough not exactly correct
                if (s == 'reference'){
                    tmp_res$CONSENSUS_MOTIF_REF = NA
                    tmp_res$CONSENSUS_MOTIF_REF_EST_COPIES = NA
                } else {
                    reference_motif = unique(motif_res_reference$CONSENSUS_MOTIF[which(motif_res_reference$REGION == r)])
                    tmp_res$CONSENSUS_MOTIF_REF = reference_motif
                    tmp_res$CONSENSUS_MOTIF_REF_EST_COPIES = tmp_res$polished_haplo_values / nchar(reference_motif)
                }
            } else {
                # in this case the clustering did not work, can be because the number of reads/contigs is too low
                tmp_res = tmp; tmp_res$COVERAGE_TR = NA; tmp_res$CONSENSUS_MOTIF = NA; tmp_res$CONSENSUS_MOTIF_EST_COPIES = NA; tmp_res$CONSENSUS_MOTIF_ALT_MOTIFS = NA
                tmp_res$CONSENSUS_MOTIF_REF = NA; tmp_res$CONSENSUS_MOTIF_REF_EST_COPIES = NA; tmp_res$SPECIFIC_MOTIF = NA
            }
            # add to results list
            motif_res[[(length(motif_res) + 1)]] = tmp_res
        }
        motif_res = rbindlist(motif_res, use.names = T)
        print(paste0('** done with ', s))
        return(motif_res)
    }
    
    # Function to call haplotypes and combine information
    callHaplo_mp <- function(s, all_regions, data){
        all_haplo = list()
        for (r in all_regions){
            # gather all data together, excluded reads and good reads
            tmp = data[which(data$sample == s & data$REGION == r),]
            excl = tmp[is.na(tmp$HAPLOTYPE),]; tmp = tmp[!is.na(tmp$HAPLOTYPE),]
            # also take unique reads
            tmp = tmp[!duplicated(tmp$READ_NAME),]
            # first we managed to identify the haplotypes
            if (nrow(tmp) >0){
                if (!(NA %in% tmp$HAPLOTYPE) & nrow(tmp) >0){
                    # if we have both haplotypes assigned, get the relative info: copies, motif, coef. of variation, length, info
                    if (1 %in% tmp$HAPLOTYPE & 2 %in% tmp$HAPLOTYPE){
                        # number of copies should be wrt consensus motif and reference motif
                        # copy number of consensus motif
                        h1_info = median(tmp$CONSENSUS_MOTIF_EST_COPIES[which(tmp$HAPLOTYPE == 1)])
                        h2_info = median(tmp$CONSENSUS_MOTIF_EST_COPIES[which(tmp$HAPLOTYPE == 2)])
                        # consensus motif
                        h1_motif = paste(unique(tmp$CONSENSUS_MOTIF[which(tmp$HAPLOTYPE == 1)]), collapse = ',')
                        h2_motif = paste(unique(tmp$CONSENSUS_MOTIF[which(tmp$HAPLOTYPE == 2)]), collapse = ',')
                        # size of haplotype
                        h1_len = paste(unique(tmp$polished_haplo_values[which(tmp$HAPLOTYPE == 1)]), collapse = ',')
                        h2_len = paste(unique(tmp$polished_haplo_values[which(tmp$HAPLOTYPE == 2)]), collapse = ',')
                        # specific motif of sample
                        h1_motif_specific = unique(tmp$SPECIFIC_MOTIF[which(tmp$HAPLOTYPE == 1)])
                        h2_motif_specific = unique(tmp$SPECIFIC_MOTIF[which(tmp$HAPLOTYPE == 2)])
                        # additional motifs        
                        h1_add_motifs = unique(tmp$CONSENSUS_MOTIF_ALT_MOTIFS[which(tmp$HAPLOTYPE == 1)])
                        h2_add_motifs = unique(tmp$CONSENSUS_MOTIF_ALT_MOTIFS[which(tmp$HAPLOTYPE == 2)])
                        # motif of reference 
                        ref_motif = unique(tmp$CONSENSUS_MOTIF_REF)
                        # copy number of reference motif
                        h1_ref_motif_cn = unique(tmp$CONSENSUS_MOTIF_REF_EST_COPIES[which(tmp$HAPLOTYPE == 1)])
                        h2_ref_motif_cn = unique(tmp$CONSENSUS_MOTIF_REF_EST_COPIES[which(tmp$HAPLOTYPE == 2)])
                        data_type = unique(tmp$DATA_TYPE)
                    # otherwise if only h1 is present, only gather h1 info
                    } else if (1 %in% tmp$HAPLOTYPE){
                        # copy number of consensus motif
                        h1_info = median(tmp$CONSENSUS_MOTIF_EST_COPIES[which(tmp$HAPLOTYPE == 1)])
                        # consensus motif
                        h1_motif = paste(unique(tmp$CONSENSUS_MOTIF[which(tmp$HAPLOTYPE == 1)]), collapse = ',')
                        # size of haplotype
                        h1_len = paste(unique(tmp$polished_haplo_values[which(tmp$HAPLOTYPE == 1)]), collapse = ',')
                        # specific motif of sample
                        h1_motif_specific = unique(tmp$SPECIFIC_MOTIF[which(tmp$HAPLOTYPE == 1)])
                        # additional motifs        
                        h1_add_motifs = unique(tmp$CONSENSUS_MOTIF_ALT_MOTIFS[which(tmp$HAPLOTYPE == 1)])
                        # motif of reference 
                        ref_motif = unique(tmp$CONSENSUS_MOTIF_REF)
                        # copy number of reference motif
                        h1_ref_motif_cn = unique(tmp$CONSENSUS_MOTIF_REF_EST_COPIES)
                        # fill in information for h1
                        h2_info = NA; h2_motif = NA; h2_len = NA; h2_motif_specific = NA; h2_add_motifs = NA; h2_ref_motif = NA; h2_ref_motif_cn = NA
                        data_type = unique(tmp$DATA_TYPE)
                    # otherwise gather only h2 info
                    } else {
                        # copy number of consensus motif
                        h2_info = median(tmp$CONSENSUS_MOTIF_EST_COPIES[which(tmp$HAPLOTYPE == 2)])
                        # consensus motif
                        h2_motif = paste(unique(tmp$CONSENSUS_MOTIF[which(tmp$HAPLOTYPE == 2)]), collapse = ',')
                        # size of haplotype
                        h2_len = paste(unique(tmp$polished_haplo_values[which(tmp$HAPLOTYPE == 2)]), collapse = ',')
                        # specific motif of sample
                        h2_motif_specific = unique(tmp$SPECIFIC_MOTIF[which(tmp$HAPLOTYPE == 2)])
                        # additional motifs
                        h2_add_motifs = unique(tmp$CONSENSUS_MOTIF_ALT_MOTIFS[which(tmp$HAPLOTYPE == 2)])
                        # motif of reference 
                        ref_motif = unique(tmp$CONSENSUS_MOTIF_REF)
                        # copy number of reference motif
                        h2_ref_motif_cn = unique(tmp$CONSENSUS_MOTIF_REF_EST_COPIES)
                        # fill in information for h1
                        h1_info = NA; h1_motif = NA; h1_len = NA; h1_motif_specific = NA; h1_add_motifs = NA; h1_ref_motif = NA; h1_ref_motif_cn = NA
                        data_type = unique(tmp$DATA_TYPE)
                    }
                    # then check excluded reads
                    if (nrow(excl) >0){ 
                        excl_info = paste(paste0(excl$CONSENSUS_MOTIF, '(', excl$CONSENSUS_MOTIF_EST_COPIES, ')'), collapse = '_')
                        excl_len = paste(excl$LENGTH_SEQUENCE, collapse = ',')
                    } else { 
                        excl_info = NA; excl_len = NA
                    }
                    # finally save a df with all info regarding TR
                    tmp_df = data.frame(SAMPLE = unique(tmp$sample), REGION = unique(tmp$REGION), H1_CONSENSUS = h1_info, H2_CONSENSUS = h2_info, 
                        H1_CONSENSUS_MOTIF = h1_motif, H2_CONSENSUS_MOTIF = h2_motif, H1_HAPLO_SIZE = h1_len, H2_HAPLO_SIZE = h2_len,
                        H1_SPECIFIC_MOTIF = h1_motif_specific, H2_SPECIFIC_MOTIF = h2_motif_specific, H1_ADDITIONAL_MOTIF = h1_add_motifs, H2_ADDITIONAL_MOTIF = h2_add_motifs,
                        REFERENCE_MOTIF = ref_motif, H1_REFERENCE_MOTIF_CN = h1_ref_motif_cn, H2_REFERENCE_MOTIF_CN = h2_ref_motif_cn,
                        EXCLUDED = excl_info, EXCLUDED_LEN = excl_len, DATA_TYPE = data_type)
                    # and add to the growing list
                    all_haplo[[(length(all_haplo) + 1)]] = tmp_df
                # if haplotypes were not identified (i.e they are NAs)
                } else {
                    excl_info = paste(paste0(excl$CONSENSUS_MOTIF, '(', excl$CONSENSUS_COPY_NUMBER, ')'), collapse = '_')
                    excl_len = paste(excl$LENGTH_SEQUENCE, collapse = ',')
                    # create the same df with NA values and add to growing list
                    tmp_df = data.frame(SAMPLE = unique(excl$sample), REGION = unique(excl$REGION), H1_CONSENSUS = NA, H2_CONSENSUS = NA, 
                        H1_CONSENSUS_MOTIF = NA, H2_CONSENSUS_MOTIF = NA, H1_HAPLO_SIZE = NA, H2_HAPLO_SIZE = NA,
                        H1_SPECIFIC_MOTIF = NA, H2_SPECIFIC_MOTIF = NA, H1_ADDITIONAL_MOTIF = NA, H2_ADDITIONAL_MOTIF = NA,
                        REFERENCE_MOTIF = ref_motif, H1_REFERENCE_MOTIF_CN = NA, H2_REFERENCE_MOTIF_CN = NA,
                        EXCLUDED = excl_info, EXCLUDED_LEN = excl_len, DATA_TYPE = data_type)
                    all_haplo[[(length(all_haplo) + 1)]] = tmp_df
                }
            }
        }
        all_haplo = rbindlist(all_haplo)
        print(paste0('** done with ', s))
        return(all_haplo)
    }
    
    # Function to compare reads-spanning and assembly-based approached based on multiple processing
    compareReadsSpanning_Asm_mp = function(s, all_regions, all_haplo, deviation){
        comparison = list()
        for (r in all_regions){
            tmp = all_haplo[which(all_haplo$SAME_NAME == s & all_haplo$REGION == r),]
            # split assembly and reads-spanning
            asm = tmp[which(tmp$DATA_TYPE == 'assembly'),]; spa = tmp[which(tmp$DATA_TYPE == 'reads-spanning'),]
            # for the assembly, keep the haplotype-aware contigs if there are multiple entries
            if (nrow(asm) >1){ asm = asm[grep('Assembly_haps', asm$SAMPLE),] }
            # first check if we have asm and reads-spanning
            if (nrow(asm) >0 & nrow(spa)){
                # pairwise comparison of haplotypes. first extract data for asm and spa for h1 and h2
                asm_h1 = asm[, c('H1', 'H1_MOTIF', 'H1_LENGTH')]; spa_h1 = spa[, c('H1', 'H1_MOTIF', 'H1_LENGTH')]
                asm_h2 = asm[, c('H2', 'H2_MOTIF', 'H2_LENGTH')]; spa_h2 = spa[, c('H2', 'H2_MOTIF', 'H2_LENGTH')]
                # then check for NAs
                if (NA %in% c(asm_h1$H1, asm_h2$H2, spa_h1$H1, spa_h2$H2)){
                    # find the NA(s)
                    na_index = which(is.na(c(asm_h1$H1, asm_h2$H2, spa_h1$H1, spa_h2$H2)))
                    if (length(na_index) == 4){
                        colnames(asm) = paste0(colnames(asm), '_ASM'); asm$COMPARISON = 'missing_all_alleles'
                        colnames(spa) = paste0(colnames(spa), '_READS_SPANNING'); spa$H1_DIFF = NA; spa$H2_DIFF = NA
                    } else if (length(na_index) == 3){
                        good_one = which(!seq(1, 4) %in% na_index)
                        colnames(asm) = paste0(colnames(asm), '_ASM')
                        colnames(spa) = paste0(colnames(spa), '_READS_SPANNING'); spa$H1_DIFF = diff_tmp; spa$H2_DIFF = NA
                        if (good_one %in% c(1, 2)){ asm$COMPARISON = 'only_1_allele_from_assembly' } else { asm$COMPARISON = 'only_1_allele_from_reads_spanning' }
                    } else if (length(na_index) == 2){
                        if (na_index[1] %in% c(1, 2) & na_index[2] %in% c(1, 2)){
                            # in this case we are missing both asm alleles
                            colnames(asm) = paste0(colnames(asm), '_ASM'); asm$COMPARISON = 'missing_both_assembly_alleles'
                            colnames(spa) = paste0(colnames(spa), '_READS_SPANNING'); spa$H1_DIFF = NA; spa$H2_DIFF = NA
                        } else if (na_index[1] %in% c(3, 4) & na_index[2] %in% c(3, 4)){
                            # in this case we are missing both reads-spanning alleles
                            colnames(spa) = paste0(colnames(spa), '_READS_SPANNING')
                            colnames(asm) = paste0(colnames(asm), '_ASM'); asm$COMPARISON = 'missing_both_reads_spanning_alleles'
                            spa$H1_DIFF = NA; spa$H2_DIFF = NA                             
                        } else {
                            # in this case there's one NA in both asm and reads-spanning -- we need to check the other allele
                            # manage asm
                            if (1 %in% na_index){ nonmissing_asm = asm_h2 } else if (2 %in% na_index){ nonmissing_asm = asm_h1 }
                            # manage spanning
                            if (3 %in% na_index){ nonmissing_spa = spa_h2 } else if (4 %in% na_index){ nonmissing_spa = spa_h1 }
                            # check nonmissing allele
                            diff_tmp = abs(as.numeric(nonmissing_asm[, 1]) - as.numeric(nonmissing_spa[, 1]))
                            if (diff_tmp <= max(1.5, as.numeric(nonmissing_asm[, 1])*deviation)){ check = 'both_miss_1_allele_other_ok' } else { check = 'both_miss_1_allele_other_different' }
                            colnames(asm) = paste0(colnames(asm), '_ASM'); asm$COMPARISON = check
                            colnames(spa) = paste0(colnames(spa), '_READS_SPANNING'); spa$H1_DIFF = diff_tmp; spa$H2_DIFF = NA
                        }
                    } else {
                        if (na_index %in% c(1, 2)){
                            if (na_index == 1){ nonmissing = asm_h2 } else { nonmissing = asm_h1 }
                            # calculate distance between alleles based on copy number
                            diff_alleles = abs(c(as.numeric(nonmissing[, 1]) - spa_h1$H1, as.numeric(nonmissing[, 1]) - spa_h2$H2))
                            if (length(unique(diff_alleles)) == 1){
                                colnames(asm) = paste0(colnames(asm), '_ASM'); asm$COMPARISON = 'assembly_miss_1_but_homozygous'
                                colnames(spa) = paste0(colnames(spa), '_READS_SPANNING'); spa$H1_DIFF = diff_alleles[1]; spa$H2_DIFF = NA
                            } else {
                                closest = which(diff_alleles == min(diff_alleles)); if (closest == 1){ closest_all = spa_h1 } else { closest_all = spa_h2 }
                                if (diff_alleles[closest] <= max(1.5, as.numeric(closest_all[,1])*deviation)){
                                    colnames(asm) = paste0(colnames(asm), '_ASM'); asm$COMPARISON = 'assembly_miss_1_allele_other_ok'
                                    colnames(spa) = paste0(colnames(spa), '_READS_SPANNING'); spa$H1_DIFF = diff_alleles[closest]; spa$H2_DIFF = NA
                                } else {
                                    colnames(asm) = paste0(colnames(asm), '_ASM'); asm$COMPARISON = 'assembly_miss_1_allele_other_different'
                                    colnames(spa) = paste0(colnames(spa), '_READS_SPANNING'); spa$H1_DIFF = diff_alleles[closest]; spa$H2_DIFF = NA
                                }
                            }
                        } else {
                            if (na_index == 3){ nonmissing = spa_h2 } else { nonmissing = spa_h1 }
                            # calculate distance between alleles based on copy number
                            diff_alleles = abs(c(as.numeric(nonmissing[, 1]) - asm_h1$H1, as.numeric(nonmissing[, 1]) - asm_h2$H2))
                            if (length(unique(diff_alleles)) == 1){
                                    colnames(asm) = paste0(colnames(asm), '_ASM'); asm$COMPARISON = 'assembly_miss_1_but_homozygous'
                                    colnames(spa) = paste0(colnames(spa), '_READS_SPANNING'); spa$H1_DIFF = diff_alleles[1]; spa$H2_DIFF = NA
                            } else {
                                closest = which(diff_alleles == min(diff_alleles)); if (closest == 1){ closest_all = asm_h1 } else { closest_all = asm_h2 }
                                if (diff_alleles[closest] <= max(1.5, as.numeric(closest_all[,1])*deviation)){
                                    colnames(asm) = paste0(colnames(asm), '_ASM'); asm$COMPARISON = 'reads_spanning_miss_1_allele_other_ok'
                                    colnames(spa) = paste0(colnames(spa), '_READS_SPANNING'); spa$H1_DIFF = diff_alleles[closest]; spa$H2_DIFF = NA
                                } else {
                                    colnames(asm) = paste0(colnames(asm), '_ASM'); asm$COMPARISON = 'reads_spanning_miss_1_allele_other_different'
                                    colnames(spa) = paste0(colnames(spa), '_READS_SPANNING'); spa$H1_DIFF = diff_alleles[closest]; spa$H2_DIFF = NA
                                }
                            }
                        }
                    }
                } else {
                    # calculate distance between alleles: h1
                    diff_h1 = abs(c(asm_h1$H1 - spa_h1$H1, asm_h1$H1 - spa_h2$H2))
                    closest_h1 = which(diff_h1 == min(diff_h1))[1]; if (closest_h1 == 1){ closest_all = spa_h1 } else { closest_all = spa_h2 }
                    if (diff_h1[closest_h1] <= max(1.5, as.numeric(closest_all[,1])*deviation)){ check_h1 = 'same' } else { check_h1 = 'different' }
                    # then h2, exclude the h1 allele
                    if (colnames(closest_all)[1] == 'H2'){ diff_h2 = abs(asm_h2$H2 - spa_h1$H1) } else { diff_h2 = abs(asm_h2$H2 - spa_h2$H2) }
                    closest_h2 = abs(closest_h1 - 1); if (closest_h2 == 1){ closest_all = spa_h1 } else { closest_all = spa_h2 }
                    if (diff_h2 <= max(1.5, as.numeric(closest_all[,1])*deviation)){ check_h2 = 'same' } else { check_h2 = 'different' }
                    if (check_h1 == 'same' & check_h2 == 'same'){ check = 'same_alleles' } else { check = 'alleles_do_not_overlap' }
                    colnames(asm) = paste0(colnames(asm), '_ASM'); asm$COMPARISON = check; colnames(spa) = paste0(colnames(spa), '_READS_SPANNING'); spa$H1_DIFF = diff_h1[closest_h1]; spa$H2_DIFF = diff_h2
                }
            } else if (nrow(asm) >0){
                colnames(asm) = paste0(colnames(asm), '_ASM'); asm$COMPARISON = 'missing_reads_spanning'                   
                colnames(spa) = paste0(colnames(spa), '_READS_SPANNING'); spa$H1_DIFF = NA; spa$H2_DIFF = NA
            } else if (nrow(spa) >0){
                colnames(spa) = paste0(colnames(spa), '_READS_SPANNING'); spa$COMPARISON = 'missing_assembly'
                colnames(asm) = paste0(colnames(asm), '_ASM'); spa$H1_DIFF = NA; spa$H2_DIFF = NA
            } 
            tmp_combined_res = cbind(asm, spa)
            if (nrow(tmp_combined_res) >0){ comparison[[(length(comparison) + 1)]] = tmp_combined_res }
        }
        comparison = rbindlist(comparison, use.names = T)
    }

# Manage arguments
    parser <- ArgumentParser()
    # add arguments: --reads_spannning is the trf output of the reads-spanning analysis
    parser$add_argument("--reads_spanning", default = 'None', help = "Output file(s) of TRF analysis on reads-spanning. Multiple files should be comman-separated.")
    # add arguments: --asm is the trf output of the assembly-based analysis
    parser$add_argument("--asm", default = 'None', help = "Output file(s) of TRF analysis on assembly. Multiple files should be comman-separated.")
    # add arguments: --phase
    parser$add_argument("--phase", default = 'None', help = "Output file(s) of PHASING analysis. Multiple files should be comman-separated.")
    # add arguments: --out
    parser$add_argument("--out", default = 'None', help = "Output directory where output will be placed.")
    # add arguments: --cpu
    parser$add_argument("--cpu", default = 2, help = "Number of CPU to use for parallel computation.")
    # add arguments: --deviation
    parser$add_argument("--deviation", default = 0.10, help = "Median absolute deviation to assign reads to the same allele.")
    # read arguments
    args <- parser$parse_args()
    inp_trf = args$reads_spanning; inp_trf = unlist(strsplit(inp_trf, ','))
    inp_asm = args$asm; inp_asm = unlist(strsplit(inp_asm, ','))
    inp_pha = args$phase; inp_pha = unlist(strsplit(inp_pha, ','))
    out_dir = args$out
    n_cpu = as.numeric(args$cpu)
    thr_mad = as.numeric(args$deviation)

    # check inputs and print summary
    if ((inp_trf[1] == 'None' & inp_asm[1] == 'None') | out_dir == 'None'){ RUN = FALSE } else { RUN = TRUE }
    if (RUN == FALSE){
        stop("Input error: Missing input file(s) or output directory detected!")
    } else {
        # decide analysis to perform based on the provided input(s) files
        if (inp_trf[1] == 'None' & inp_asm[1] != 'None'){
            anal_type = 'assembly'
        } else if (inp_trf[1] != 'None' & inp_asm == 'None'){
            anal_type = 'reads-spanning'
        } else {
            anal_type = 'reads-spanning + assembly + comparison'
        }
        #cat('\n** Haplotyping Tandem Repeats\n\n')
        cat(paste0('**** input TRF (reads-spanning) selected -> ', inp_trf, '\n'))
        cat(paste0('**** input TRF (assembly-based) selected -> ', inp_asm, '\n'))
        cat(paste0('**** input PHASING selected -> ', inp_pha, '\n'))
        cat(paste0('**** analysis type based on input(s) --> ', anal_type, '\n'))
        cat(paste0('**** output directory -> ', out_dir, '\n'))
    }

# Main script
    # 1. Read and combine all input TRF files
    cat('**** Analysis started\n')
    cat('****** Reading TRF input(s)\n')
    if (anal_type == 'reads-spanning'){
        # read reads-spanning TRF
        all_trf = list()
        for (i in 1:length(inp_trf)){
            # load trf datasets
            data = fread(inp_trf[i], h=T, stringsAsFactors = F)
            data$BATCH_TRF = i
            data$DATA_TYPE = 'reads-spanning'
            all_trf[[(length(all_trf) + 1)]] = data
        }
        all_trf = rbindlist(all_trf)
    } else if (anal_type == 'assembly'){
        # read assembly-based TRF        
        all_trf = list()
        for (i in 1:length(inp_asm)){
            # load trf datasets
            data = fread(inp_asm[i], h=T, stringsAsFactors = F)
            data$BATCH_TRF = i
            data$DATA_TYPE = 'assembly'
            all_trf[[(length(all_trf) + 1)]] = data
        }
        all_trf = rbindlist(all_trf)
    } else {
        # read reads-spanning TRF
        all_trf = list()
        for (i in 1:length(inp_trf)){
            # load trf datasets
            data = fread(inp_trf[i], h=T, stringsAsFactors = F)
            data$BATCH_TRF = i
            data$DATA_TYPE = 'reads-spanning'
            all_trf[[(length(all_trf) + 1)]] = data
        }
        # read assembly-based TRF
        for (i in 1:length(inp_asm)){
            # load trf datasets
            data = fread(inp_asm[i], h=T, stringsAsFactors = F)
            data$BATCH_TRF = i
            data$DATA_TYPE = 'assembly'
            all_trf[[(length(all_trf) + 1)]] = data
        }
        all_trf = rbindlist(all_trf)
    }

    # 2. Read and combine all input PHASING files
    if (inp_pha != 'None'){
        cat('****** Reading PHASING input\n')
        all_pha = list()
        for (i in 1:length(inp_pha)){
            # load trf datasets
            data = fread(inp_pha[i], h=T, stringsAsFactors = F)
            data$BATCH_PHASING = i
            all_pha[[(length(all_pha) + 1)]] = data
        }
        all_pha = rbindlist(all_pha)
    } else {
        cat('** PHASING file not provided. This could lead to loss of accuracy.\n')
        all_pha = data.frame(SAMPLE = as.character(), READ_ID = as.character(), HAPLOTYPE = as.character())
    }

    # 3. Merge trf and phases -- exclude reference before -- this in case reads-spanning are available
    reference = all_trf[which(all_trf$SAMPLE_NAME == 'reference'),]
    all_trf = all_trf[which(all_trf$SAMPLE_NAME != 'reference'),]
    trf_pha = merge(all_trf, all_pha, by.x = 'READ_NAME', by.y = 'READ_ID', all.x = T)
    trf_pha = rbind.fill(trf_pha, reference)
  
    # 4. Let's adjust the motifs -- generate a consensus -- essentially merging the same motifs together -- this for both analyses
    cat('****** Merging similar motifs together\n')
    all_motifs = data.frame(motif = unique(trf_pha$TRF_MOTIF), stringsAsFactors = F)
    all_motifs = all_motifs[!is.na(all_motifs$motif),]
    # manage the reference
    trf_pha$READ_NAME[is.na(trf_pha$READ_NAME)] <- 'reference'
    # alternative way using multiprocessing
    main_motifs = unlist(mclapply(all_motifs, mergeMotif_mp, mc.cores = n_cpu))
    all_motifs = data.frame(motif = all_motifs, UNIFORM_MOTIF = main_motifs)
    trf_pha = merge(trf_pha, all_motifs, by.x = 'TRF_MOTIF', by.y = 'motif')

    # to make sure there are no conflicts, change the read name
    trf_pha$UNIQUE_NAME = paste(trf_pha$READ_NAME, trf_pha$SAMPLE_NAME, trf_pha$REGION, sep = "___")

    # 5. add the TR size by substracting 20 (10*2 padding) to LENGTH_SEQUENCE
    trf_pha$TR_LENGHT = trf_pha$LENGTH_SEQUENCE

    # 6. good to exclude duplicated reads otherwise results would be biased towards sequences with more complex motifs (where TRF finds multiple matches)
    #trf_pha = trf_pha[which(trf_pha$SEQUENCE_WITH_WINDOW != ''),]
    trf_pha_nodup = trf_pha[!duplicated(trf_pha$UNIQUE_NAME),]; dups = trf_pha[duplicated(trf_pha$UNIQUE_NAME),]

    # 7. we do genotyping on the sizes
    cat('****** Genotyping TRs\n')
    # temporarily disable warnings
    defaultW <- getOption("warn"); options(warn = -1)
    # split reads-spanning and assembly-based
    reads_span = trf_pha_nodup[which(trf_pha_nodup$DATA_TYPE == 'reads-spanning'),]; asm = trf_pha_nodup[which(trf_pha_nodup$DATA_TYPE == 'assembly'),]
    # haplotyping is done with multiple processors (1 for each sample)
    if (anal_type %in% c('reads-spanning', 'reads-spanning + assembly + comparison')){
        cat(paste0('******** Processing reads-spanning data\n'))
        all_samples = unique(reads_span$SAMPLE_NAME); all_regions = unique(reads_span$REGION)
        res_reads_spanning = rbindlist(mclapply(all_samples, haplotyping_mp, reads_span = reads_span, all_regions = all_regions, type = 'reads_spanning', thr_mad = thr_mad, mc.cores = n_cpu), fill = T)
        res_reads_spanning$DATA_TYPE = 'reads-spanning'
    }
    if (anal_type %in% c('assembly', 'reads-spanning + assembly + comparison')){
        cat(paste0('******** Processing assembly-based data\n'))
        all_samples = unique(asm$SAMPLE_NAME); all_regions = unique(asm$REGION)
        res_asm = rbindlist(mclapply(all_samples, haplotyping_mp, reads_span = asm, all_regions = all_regions, type = 'asm', thr_mad = thr_mad, mc.cores = n_cpu), fill = T)
        res_asm$DATA_TYPE = 'assembly'
    }
    # Merge all results together
    if (anal_type == 'reads-spanning + assembly + comparison'){ all_res = rbind(res_reads_spanning, res_asm) } else if (anal_type == 'reads-spanning'){ all_res = res_reads_spanning } else if (anal_type == 'assembly'){ all_res = res_asm }
    options(warn = defaultW)
    # re-assign the unique identifier after haplotype calling
    all_res$UNIQUE_NAME = paste(all_res$READ_NAME, all_res$sample, all_res$REGION, sep = "___")

    # 8. now let's bring the duplicates in again and assign the correct haplotype based on the other duplicate
    # this is also now done using multiple processors
    dups_name = unique(dups$UNIQUE_NAME)
    tmp_df_dups = data.frame(COPIES_TRF = dups$COPIES_TRF, UNIFORM_MOTIF = dups$UNIFORM_MOTIF, REGION = dups$REGION, PASSES = dups$PASSES, READ_QUALITY = dups$READ_QUALITY,
        LENGTH_SEQUENCE = dups$LENGTH_SEQUENCE, READ_NAME = dups$READ_NAME, START_TRF = dups$START_TRF, END_TRF = dups$END_TRF, TRF_PERC_MATCH = dups$TRF_PERC_MATCH, TRF_PERC_INDEL = dups$TRF_PERC_INDEL,
        DATA_TYPE = dups$DATA_TYPE, UNIQUE_NAME = dups$UNIQUE_NAME)
    # add information about type, sample, haplo_value, polished_reads and polished_haplo_values
    tmp_info = all_res[, c('UNIQUE_NAME', 'type', 'sample', 'haplo_value', 'polished_reads', 'polished_haplo_values', 'HAPLOTYPE')]
    tmp_df_dups = merge(tmp_df_dups, tmp_info, by = 'UNIQUE_NAME')
    # finally combine the duplicate information with the unique information
    all_res_combined = rbind(all_res, tmp_df_dups)

    # 9. now we should look at the motifs -- implemented parallel computing
    cat('****** Generating consensus motifs\n')
    all_samples = unique(all_res$sample); all_regions = unique(all_res$REGION); motif_res = list()
    # first run on the reference genome
    motif_res_reference = generateConsens_mp(s = 'reference', all_regions, all_res = all_res_combined, motif_res_reference = NA)
    motif_res = rbindlist(mclapply(all_samples[which(all_samples != 'reference')], generateConsens_mp, all_regions = all_regions, all_res = all_res_combined, motif_res_reference = motif_res_reference, mc.cores = n_cpu), use.names=TRUE)

    # 10. finally call haplotypes -- implemented parallel computing
    cat('****** Haplotype calling\n')
    all_samples = unique(motif_res$sample); all_regions = unique(motif_res$REGION)
    all_haplo = rbindlist(mclapply(all_samples, callHaplo_mp, all_regions = all_regions, data = motif_res, mc.cores = n_cpu))

    # 11. if analysis type was different from the combined, we are done here: save output tables
    cat('****** Generating tables\n')
    out_dir = str_replace_all(out_dir, '/$', ''); if (!dir.exists(out_dir)){ system(paste0('mkdir ', out_dir)) }
    write.table(motif_res, paste0(out_dir, '/haplotyping_single_reads.txt'), quote=F, row.names = F, sep = '\t')
    write.table(all_haplo, paste0(out_dir, '/haplotyping_single_samples.txt'), quote=F, row.names = F, sep = '\t')

    # 12. if both reads-spanning and assembly were submitted, we should compare them and make a unified call
    # this will assume the names are the same. Will split with '.' and take the name before that.
    if (anal_type == 'reads-spanning + assembly + comparison'){
        cat('****** Comparing reads-spanning with assembly\n')
        # add a unique identifier to match samples reads-spanning and assembly results
        all_haplo$SAME_NAME = str_split_fixed(all_haplo$SAMPLE, '__', 2)[, 1]
        # exclude reference
        reference = all_haplo[which(all_haplo$SAME_NAME == 'reference'),]; all_haplo = all_haplo[which(all_haplo$SAME_NAME != 'reference'),]
        # run multiprocessing for each sample independently
        all_samples = unique(all_haplo$SAME_NAME); all_regions = unique(all_haplo$REGION)
        all_haplo_annotated = rbindlist(mclapply(all_samples, compareReadsSpanning_Asm_mp, all_regions = all_regions, all_haplo = all_haplo, deviation = thr_mad, mc.cores = n_cpu))
        # then make a final decision
        decisions = makeDecision(all_haplo_annotated)
        # finally also save this dataset
        write.table(decisions, paste0(out_dir, '/haplotyping_reads_spanning_VS_assembly.txt'), quote=F, row.names = F, sep = '\t')
    }
    cat('**** Analysis done!!\n')
