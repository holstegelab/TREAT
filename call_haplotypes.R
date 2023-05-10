#!/usr/bin/R

#########################################
# Script to analyze TRF output for both #
# single-reads and assembly and come up #
# with a way of suggesting which is the #
# most likely haplotyping structure.    #
#########################################

# Libraries
    list.of.packages <- c("plyr", "data.table", "argparse", "stringr", "parallel", "spgs", "seqinr")
    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
    if(length(new.packages)) install.packages(new.packages)
    library(plyr)
    library(data.table)
    library(argparse)
    library(parallel)
    library(stringr)
    library(spgs)

# Functions
    # Function to make permutations of letters given a word, used to check motifs -- in use
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

    # Function to perform guided haplotyping using phasing information using the sizes -- in use
    PhasingBased_haplotyping_size <- function(reads_df, sample_name, thr_mad, region){
        # exclude reads with NA at length
        reads_df = reads_df[!is.na(reads_df$LENGTH_SEQUENCE),]
        # if we are dealing with a Y chromosome entry, it is already phased and it's 1 haplotype
        chrom = str_replace_all(str_split_fixed(region, ':', 2)[, 1], 'chr', '')
        if (chrom == 'Y' || chrom == 'y'){ reads_df$HAPLOTYPE = as.character(reads_df$HAPLOTYPE); reads_df$HAPLOTYPE = 1 }
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

    # Polish haplotypes after phasing is done -- in use
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

    # Function to fit K-means clustering to the data to find haplotypes when no phasing information is available -- in use
    KMeansBased_haplotyping = function(reads, thr, min_support, thr_mad_orig, type, region){
        # define ploidy depending on the analysis type
        ploidy = ifelse(type %in% c("single_tissue", "single_sample"), 2, 3)
        # if we are dealing with a Y chromosome entry, it is already phased and it's 1 haplotype
        chrom = str_replace_all(str_split_fixed(region, ':', 2)[, 1], 'chr', '')
        if (chrom == 'Y' || chrom == 'y'){ ploidy = 1 }
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

    # Function to polish clustering of TR when no phasing information is available -- in use
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

    # Function to call haplotypes when dealing with assembly-based TRF -- in use
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

    # Function to merge motifs based on multiple processing -- in use
    mergeMotif_mp = function(motif){
        all_perm = permutMotif(as.character(motif))
        all_perm = all_perm[order(all_perm)]
        return(all_perm[1])
    }

    # Function to do haplotyping based on multiple processing -- in use
    haplotyping_mp = function(s, reads_span, all_regions, type, thr_mad){
        # initialize dataframe for results
        tmp_res = data.frame()
        for (r in all_regions){
            # get data of the sample and the region of interest -- depending on type (reads_spanning or asm)
            tmp_data = reads_span[which(reads_span$SAMPLE_NAME == s & reads_span$REGION == r),]
            # good idea to exclude NAs here
            tmp_data = tmp_data[!is.na(tmp_data$LENGTH_SEQUENCE),]
            if (nrow(tmp_data) >0){
                reads_df = tmp_data[, c('COPIES_TRF', 'HAPLOTYPE', 'UNIFORM_MOTIF', 'REGION', 'PASSES', 'READ_QUALITY', 'LENGTH_SEQUENCE', 'READ_NAME', 'START_TRF', 'END_TRF', 'TRF_PERC_MATCH', 'TRF_PERC_INDEL')]
                if (type == 'reads_spanning'){
                    phased_data = PhasingBased_haplotyping_size(reads_df, sample_name = s, thr_mad, region = r)
                } else {
                    phased_data = assemblyBased_size(data = reads_df, sample_name = s, region = r, thr_mad)
                }
                polished_data = polishHaplo_afterPhasing_size(phased_data, thr_mad)
                tmp_res = rbind.fill(tmp_res, polished_data)
            }
        }
        return(tmp_res)
    }

    # Function to align motifs and take only one -- in use
    aln_motifs = function(h1){
        # the problem is that TRF finds different motifs that represent be the same motif
        # example: AAAAG and AAAAAGAAAAG --> AAAG(2)
        # take all motifs
        all_motifs = as.character(unique(h1$UNIFORM_MOTIF))
        all_motifs = all_motifs[!is.na(all_motifs)]
        # do this only if there are more than 1 motif -- otherwise there's nothing to merge
        if (length(all_motifs) >1){
            # so we take the smaller motif
            smaller_motif = all_motifs[order(nchar(all_motifs))][1]
            other_motifs = all_motifs[order(nchar(all_motifs))][2:length(all_motifs)]
            df = data.frame(motif = smaller_motif, copies = 1, tomerge = 'yes', motif_to_merge = smaller_motif)
            # then we iterate over the other motifs
            for (m in other_motifs){
                # we count the occurrences of the small motif in the other motif
                copies_smaller_motif = str_count(m, smaller_motif)
                # also check reverse complement
                copies_smaller_motif_reverse = str_count(m, toupper(reverseComplement(x = smaller_motif)))
                copies_smaller_motif = max(copies_smaller_motif, copies_smaller_motif_reverse)
                # then generate the corresponding sequence of the small motif repetead as many times as it was found in the longer motif
                smaller_motif_seq = paste(rep(smaller_motif, copies_smaller_motif), collapse = '')
                # now we need to check how similar the motifs are -- allow 10% deviation max
                deviation = 1 - (nchar(smaller_motif_seq) / nchar(m))
                tomerge = ifelse(deviation < 0.10, 'yes', 'no')
                # create a dataframe
                df = rbind(df, data.frame(motif = m, copies = copies_smaller_motif, tomerge = tomerge, motif_to_merge = smaller_motif))
            }
            # finally we need to change the motifs to be merged in the data
            for (i in 1:nrow(df)){
                if (df$copies[i] != 1){
                    # update motif and copy number if it needs to be updated
                    if (df$tomerge[i] == 'yes'){
                        h1 = h1[which(h1$UNIFORM_MOTIF != df$motif[i])]
                    }
                }
            }
        }
        return(h1)
    }

    # Function to generate consensus motif using the majority rule consensus -- in use
    motif_generalization = function(h1){
        # idea is to calculate a consensus motif, then iterate on the reads and generate the consensus representation
        # the first step is to recognized whether different motifs are actually a repetition of a smaller motif
        h1 = aln_motifs(h1)
        motifs_to_use = c()
        # in general there can be two scenarios: 
        # 1. different trf matches cover different part of the sequence --> motifs should be combined
        # 2. different trf matches cover the same sequence --> motifs are variable and we should report 1 motif only
        # before going to the consensus motif generation, we should identify the scenario we are in
        # to do that, first, calculate the fraction of sequence covered by each read
        h1$COVERAGE_TR = (h1$END_TRF - h1$START_TRF + 1) / h1$LENGTH_SEQUENCE
        # sort results by the highest coverage value
        h1 = h1[order(-h1$COVERAGE_TR),]
        # next, take the motif that covers the most of the TR
        best_motif = as.character(h1$UNIFORM_MOTIF[order(-h1$COVERAGE_TR)][1])
        best_motif = best_motif[!is.na(best_motif)]
        # sometimes, there are no motifs, so we should skip the whole thing
        if (length(best_motif) >0){
            motifs_to_use = best_motif
            # define a sequence of indexes representing the coverage of the best motif based on TRF_START and TRF_END
            best_motif_seq = seq(h1$START_TRF[order(-h1$COVERAGE_TR)][1], h1$END_TRF[order(-h1$COVERAGE_TR)][1])
            # then see if adding other motifs improves the sequence coverage
            alt_motifs = h1[which(h1$UNIFORM_MOTIF != best_motif),]
            if (nrow(alt_motifs) >1){
                for (i in 1:nrow(alt_motifs)){
                    # define a sequence of indexes representing the coverage of the best motif based on TRF_START and TRF_END
                    alt_motif_seq = seq(alt_motifs$START_TRF[i], alt_motifs$END_TRF[i])
                    # join the two list of indexes and remove duplicates
                    combined_seq = c(best_motif_seq, alt_motif_seq)
                    combined_seq = combined_seq[!duplicated(combined_seq)]
                    # calculate coverage of this new combined sequence
                    combined_seq_coverage = length(combined_seq) / max(h1$LENGTH_SEQUENCE[1], alt_motifs$LENGTH_SEQUENCE[i])
                    # check if the coverage of the combined sequence is higher than the single sequence
                    if (combined_seq_coverage > (max(na.omit(h1$COVERAGE_TR)) + 0.10)){
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
                h1_noNA = h1[!is.na(h1$START_TRF),]
                fit = kmeans(x = h1_noNA[, c('START_TRF', 'END_TRF', 'TRF_PERC_MATCH')], centers = length(motifs_to_use), algorithm = 'Lloyd')
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
        } else {
            h1$CONSENSUS_MOTIF = NA
            h1$CONSENSUS_MOTIF_EST_COPIES = NA
            h1$CONSENSUS_MOTIF_ALT_MOTIFS = NA
            h1$SPECIFIC_MOTIF = NA

        }
        # return the same object we used as input with additional columns
        return(h1)
    }
    
    # Function to generate consensus motif for each sample and region based on multiple processing -- in use
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
        return(motif_res)
    }
    
    # Function to call haplotypes and combine information -- in use
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
        return(all_haplo)
    }
    
    # Function to compare reads-spanning and assembly-based approached based on multiple processing -- in use
    comparison_faster = function(s, all_regions, all_haplo, deviation){
        # initialize container
        comparison = list()
        # main loop over regions
        for (r in all_regions){
            # extract data of the sample in the region of interest
            tmp = all_haplo[which(all_haplo$SAME_NAME == s & all_haplo$REGION == r),]
            # split assembly and reads-spanning
            asm = tmp[which(tmp$DATA_TYPE == 'assembly'),]; spa = tmp[which(tmp$DATA_TYPE == 'reads-spanning'),]
            # for the assembly, keep the haplotype-aware contigs if there are multiple entries
            if (nrow(asm) >1){ asm = asm[grep('Assembly_haps', asm$SAMPLE),] }
            # first check if we have asm and reads-spanning
            if (nrow(asm) >0 & nrow(spa) >0){
                # adjust homozygous genotypes
                # asm
                asm$H2_CONSENSUS[is.na(asm$H2_CONSENSUS)] = asm$H1_CONSENSUS[is.na(asm$H2_CONSENSUS)]
                asm$H2_CONSENSUS_MOTIF[is.na(asm$H2_CONSENSUS_MOTIF)] = asm$H1_CONSENSUS_MOTIF[is.na(asm$H2_CONSENSUS_MOTIF)]
                asm$H2_HAPLO_SIZE[is.na(asm$H2_HAPLO_SIZE)] = asm$H1_HAPLO_SIZE[is.na(asm$H2_HAPLO_SIZE)]
                if (NA %in% asm$H1_CONSENSUS){
                    asm$H1_CONSENSUS[is.na(asm$H1_CONSENSUS)] = asm$H2_CONSENSUS[is.na(asm$H1_CONSENSUS)]
                    asm$H1_CONSENSUS_MOTIF[is.na(asm$H1_CONSENSUS_MOTIF)] = asm$H2_CONSENSUS_MOTIF[is.na(asm$H1_CONSENSUS_MOTIF)]
                    asm$H1_HAPLO_SIZE[is.na(asm$H1_HAPLO_SIZE)] = asm$H2_HAPLO_SIZE[is.na(asm$H1_HAPLO_SIZE)]
                }
                # spanning
                spa$H2_CONSENSUS[is.na(spa$H2_CONSENSUS)] = spa$H1_CONSENSUS[is.na(spa$H2_CONSENSUS)]
                spa$H2_CONSENSUS_MOTIF[is.na(spa$H2_CONSENSUS_MOTIF)] = spa$H1_CONSENSUS_MOTIF[is.na(spa$H2_CONSENSUS_MOTIF)]
                spa$H2_HAPLO_SIZE[is.na(spa$H2_HAPLO_SIZE)] = spa$H1_HAPLO_SIZE[is.na(spa$H2_HAPLO_SIZE)]
                if (NA %in% spa$H1_CONSENSUS){
                    spa$H1_CONSENSUS[is.na(spa$H1_CONSENSUS)] = spa$H2_CONSENSUS[is.na(spa$H1_CONSENSUS)]
                    spa$H1_CONSENSUS_MOTIF[is.na(spa$H1_CONSENSUS_MOTIF)] = spa$H2_CONSENSUS_MOTIF[is.na(spa$H1_CONSENSUS_MOTIF)]
                    spa$H1_HAPLO_SIZE[is.na(spa$H1_HAPLO_SIZE)] = spa$H2_HAPLO_SIZE[is.na(spa$H1_HAPLO_SIZE)]
                }
                # check if there are NAs (but there shouldn't be at this point)
                if (NA %in% c(asm$H1_HAPLO_SIZE, asm$H2_HAPLO_SIZE, spa$H1_HAPLO_SIZE, spa$H2_HAPLO_SIZE)){
                    print(paste0('help!! ', s, '-', r))
                } else {
                    # separate data of asm and spa, for both alleles
                    asm_h1 = asm[, c('H1_CONSENSUS', 'H1_CONSENSUS_MOTIF', 'H1_HAPLO_SIZE')]; spa_h1 = spa[, c('H1_CONSENSUS', 'H1_CONSENSUS_MOTIF', 'H1_HAPLO_SIZE')]
                    asm_h2 = asm[, c('H2_CONSENSUS', 'H2_CONSENSUS_MOTIF', 'H2_HAPLO_SIZE')]; spa_h2 = spa[, c('H2_CONSENSUS', 'H2_CONSENSUS_MOTIF', 'H2_HAPLO_SIZE')]
                    # make sure variables are numeric
                    spa_h1$H1_HAPLO_SIZE = as.numeric(spa_h1$H1_HAPLO_SIZE); spa_h2$H2_HAPLO_SIZE = as.numeric(spa_h2$H2_HAPLO_SIZE)
                    asm_h1$H1_HAPLO_SIZE = as.numeric(asm_h1$H1_HAPLO_SIZE); asm_h2$H2_HAPLO_SIZE = as.numeric(asm_h2$H2_HAPLO_SIZE)
                    # calculate distance between alleles: h1
                    diff_h1 = abs(c(asm_h1$H1_HAPLO_SIZE - spa_h1$H1_HAPLO_SIZE, asm_h1$H1_HAPLO_SIZE - spa_h2$H2_HAPLO_SIZE))
                    closest_h1 = which(diff_h1 == min(diff_h1))[1]; if (closest_h1 == 1){ closest_all = spa_h1 } else { closest_all = spa_h2 }
                    if (diff_h1[closest_h1] <= max(1.5, as.numeric(closest_all[, 3])*deviation)){ check_h1 = 'same' } else { check_h1 = 'different' }
                    # then h2, exclude the h1 allele
                    if (colnames(closest_all)[1] == 'H2'){ diff_h2 = abs(asm_h2$H2_HAPLO_SIZE - spa_h1$H1_HAPLO_SIZE) } else { diff_h2 = abs(asm_h2$H2_HAPLO_SIZE - spa_h2$H2_HAPLO_SIZE) }
                    closest_h2 = abs(closest_h1 - 1); if (closest_h2 == 1){ closest_all = spa_h1 } else { closest_all = spa_h2 }
                    if (diff_h2 <= max(1.5, as.numeric(closest_all[, 3])*deviation)){ check_h2 = 'same' } else { check_h2 = 'different' }
                    # summary of comparison
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
        comparison = rbindlist(comparison, use.names=TRUE)
        return(comparison)
    }

    # Function to make final decision between assembly and reads-spanning -- in use
    makeDecision_faster = function(all_haplo_annotated, n_cpu){
        # fix some NA in the columns of interes
        all_haplo_annotated[all_haplo_annotated == 'NA'] = NA
        all_haplo_annotated$H2_REFERENCE_MOTIF_CN_ASM[is.na(all_haplo_annotated$H2_REFERENCE_MOTIF_CN_ASM)] = all_haplo_annotated$H1_REFERENCE_MOTIF_CN_ASM[is.na(all_haplo_annotated$H2_REFERENCE_MOTIF_CN_ASM)]
        all_haplo_annotated$H2_REFERENCE_MOTIF_CN_READS_SPANNING[is.na(all_haplo_annotated$H2_REFERENCE_MOTIF_CN_READS_SPANNING)] = all_haplo_annotated$H1_REFERENCE_MOTIF_CN_READS_SPANNING[is.na(all_haplo_annotated$H2_REFERENCE_MOTIF_CN_READS_SPANNING)]
        # first, we take out all regions (the majority) where reads-spanning and assembly do agree
        regions_agree = all_haplo_annotated[which(all_haplo_annotated$COMPARISON == 'same_alleles'),]
        # set decision, alleles, motifs and copies
        regions_agree$DECISION = 'OK_BOTH'
        regions_agree$H1_HAPLO_FINAL = regions_agree$H1_HAPLO_SIZE_READS_SPANNING; regions_agree$H2_HAPLO_FINAL = regions_agree$H2_HAPLO_SIZE_READS_SPANNING
        regions_agree$H1_MOTIF_FINAL = regions_agree$H1_CONSENSUS_MOTIF_ASM; regions_agree$H2_MOTIF_FINAL = regions_agree$H2_CONSENSUS_MOTIF_ASM
        regions_agree$H1_COPIE_FINAL = regions_agree$H1_CONSENSUS_ASM; regions_agree$H2_COPIE_FINAL = regions_agree$H2_CONSENSUS_ASM
        regions_agree$H1_REFCOPY_FINAL = regions_agree$H1_REFERENCE_MOTIF_CN_ASM; regions_agree$H2_REFCOPY_FINAL = regions_agree$H2_REFERENCE_MOTIF_CN_ASM
        # sometimes, there are NAs (this happens when no motif was found either in assembly or the reads-spanning). Because we take values from assembly, this can be a problem. So here we correct that.
        regions_agree$H2_CONSENSUS_MOTIF_READS_SPANNING[is.na(regions_agree$H2_CONSENSUS_MOTIF_READS_SPANNING)] = regions_agree$H1_CONSENSUS_MOTIF_READS_SPANNING[is.na(regions_agree$H2_CONSENSUS_MOTIF_READS_SPANNING)]
        regions_agree$H1_MOTIF_FINAL[is.na(regions_agree$H1_MOTIF_FINAL)] = regions_agree$H1_CONSENSUS_MOTIF_READS_SPANNING[is.na(regions_agree$H1_MOTIF_FINAL)]
        regions_agree$H2_MOTIF_FINAL[is.na(regions_agree$H2_MOTIF_FINAL)] = regions_agree$H2_CONSENSUS_MOTIF_READS_SPANNING[is.na(regions_agree$H2_MOTIF_FINAL)]
        regions_agree$H1_COPIE_FINAL[is.na(regions_agree$H1_COPIE_FINAL)] = regions_agree$H1_CONSENSUS_READS_SPANNING[is.na(regions_agree$H1_COPIE_FINAL)]
        regions_agree$H2_COPIE_FINAL[is.na(regions_agree$H2_COPIE_FINAL)] = regions_agree$H2_CONSENSUS_READS_SPANNING[is.na(regions_agree$H2_COPIE_FINAL)]
        regions_agree$H1_REFCOPY_FINAL[is.na(regions_agree$H1_REFCOPY_FINAL)] = regions_agree$H1_REFERENCE_MOTIF_CN_READS_SPANNING[is.na(regions_agree$H1_REFCOPY_FINAL)]
        regions_agree$H2_REFCOPY_FINAL[is.na(regions_agree$H2_REFCOPY_FINAL)] = regions_agree$H2_REFERENCE_MOTIF_CN_READS_SPANNING[is.na(regions_agree$H2_REFCOPY_FINAL)]
        # then we look at the missings: split missing assembly from missing reads_spanning
        regions_missing_asm = all_haplo_annotated[which(all_haplo_annotated$COMPARISON == 'missing_assembly'),]
        if (nrow(regions_missing_asm) >0){ 
            regions_missing_asm$DECISION = 'PRIORITY_READS_SPANNING'
            regions_missing_asm$H1_HAPLO_FINAL = regions_missing_asm$H1_HAPLO_SIZE_READS_SPANNING; regions_missing_asm$H2_HAPLO_FINAL = regions_missing_asm$H2_HAPLO_SIZE_READS_SPANNING
            regions_missing_asm$H1_MOTIF_FINAL = regions_missing_asm$H1_CONSENSUS_MOTIF_READS_SPANNING; regions_missing_asm$H2_MOTIF_FINAL = regions_missing_asm$H2_CONSENSUS_MOTIF_READS_SPANNING
            regions_missing_asm$H1_COPIE_FINAL = regions_missing_asm$H1_CONSENSUS_READS_SPANNING; regions_missing_asm$H2_COPIE_FINAL = regions_missing_asm$H2_CONSENSUS_READS_SPANNING
            regions_missing_asm$H1_REFCOPY_FINAL = regions_missing_asm$H1_REFERENCE_MOTIF_CN_READS_SPANNING; regions_missing_asm$H2_REFCOPY_FINAL = regions_missing_asm$H2_REFERENCE_MOTIF_CN_READS_SPANNING
            # fix NA
            regions_missing_asm$H2_HAPLO_FINAL[is.na(regions_missing_asm$H2_HAPLO_FINAL)] = regions_missing_asm$H1_HAPLO_FINAL[is.na(regions_missing_asm$H2_HAPLO_FINAL)]
            regions_missing_asm$H2_MOTIF_FINAL[is.na(regions_missing_asm$H2_MOTIF_FINAL)] = regions_missing_asm$H1_MOTIF_FINAL[is.na(regions_missing_asm$H2_MOTIF_FINAL)]
            regions_missing_asm$H2_COPIE_FINAL[is.na(regions_missing_asm$H2_COPIE_FINAL)] = regions_missing_asm$H1_COPIE_FINAL[is.na(regions_missing_asm$H2_COPIE_FINAL)]
            regions_missing_asm$H2_REFCOPY_FINAL[is.na(regions_missing_asm$H2_REFCOPY_FINAL)] = regions_missing_asm$H1_REFCOPY_FINAL[is.na(regions_missing_asm$H2_REFCOPY_FINAL)]
        }
        regions_missing_spa = all_haplo_annotated[which(all_haplo_annotated$COMPARISON == 'missing_reads_spanning'),]
        if (nrow(regions_missing_spa) >0){ 
            regions_missing_spa$DECISION = 'PRIORITY_ASSEMBLY' 
            regions_missing_spa$H1_HAPLO_FINAL = regions_missing_spa$H1_HAPLO_SIZE_ASM; regions_missing_spa$H2_HAPLO_FINAL = regions_missing_spa$H2_HAPLO_SIZE_ASM
            regions_missing_spa$H1_MOTIF_FINAL = regions_missing_spa$H1_CONSENSUS_MOTIF_ASM; regions_missing_spa$H2_MOTIF_FINAL = regions_missing_spa$H2_CONSENSUS_MOTIF_ASM
            regions_missing_spa$H1_COPIE_FINAL = regions_missing_spa$H1_CONSENSUS_ASM; regions_missing_spa$H2_COPIE_FINAL = regions_missing_spa$H2_CONSENSUS_ASM
            regions_missing_spa$H1_REFCOPY_FINAL = regions_missing_spa$H1_REFERENCE_MOTIF_CN_ASM; regions_missing_spa$H2_REFCOPY_FINAL = regions_missing_spa$H2_REFERENCE_MOTIF_CN_ASM
            # fix NA
            regions_missing_spa$H2_HAPLO_FINAL[is.na(regions_missing_spa$H2_HAPLO_FINAL)] = regions_missing_spa$H1_HAPLO_FINAL[is.na(regions_missing_spa$H2_HAPLO_FINAL)]
            regions_missing_spa$H2_MOTIF_FINAL[is.na(regions_missing_spa$H2_MOTIF_FINAL)] = regions_missing_spa$H1_MOTIF_FINAL[is.na(regions_missing_spa$H2_MOTIF_FINAL)]
            regions_missing_spa$H2_COPIE_FINAL[is.na(regions_missing_spa$H2_COPIE_FINAL)] = regions_missing_spa$H1_COPIE_FINAL[is.na(regions_missing_spa$H2_COPIE_FINAL)]
            regions_missing_spa$H2_REFCOPY_FINAL[is.na(regions_missing_spa$H2_REFCOPY_FINAL)] = regions_missing_spa$H1_REFCOPY_FINAL[is.na(regions_missing_spa$H2_REFCOPY_FINAL)]
        }
        # then we look at the inconsistent regions. In general, prioritize more variability (2 alleles instead of 1)
        inconsistent = all_haplo_annotated[which(all_haplo_annotated$COMPARISON == 'alleles_do_not_overlap'),]
        inconsistent$DECISION = NA
        # loop over the inconsistent regions
        inconsistent_solved = rbindlist(mclapply(1:nrow(inconsistent), solve_inconsistency, inconsistent = inconsistent, mc.cores=n_cpu))
        # fix NA
        inconsistent_solved$H2_HAPLO_FINAL[is.na(inconsistent_solved$H2_HAPLO_FINAL)] = inconsistent_solved$H1_HAPLO_FINAL[is.na(inconsistent_solved$H2_HAPLO_FINAL)]
        inconsistent_solved$H2_MOTIF_FINAL[is.na(inconsistent_solved$H2_MOTIF_FINAL)] = inconsistent_solved$H1_MOTIF_FINAL[is.na(inconsistent_solved$H2_MOTIF_FINAL)]
        inconsistent_solved$H2_COPIE_FINAL[is.na(inconsistent_solved$H2_COPIE_FINAL)] = inconsistent_solved$H1_COPIE_FINAL[is.na(inconsistent_solved$H2_COPIE_FINAL)]
        inconsistent_solved$H2_REFCOPY_FINAL[is.na(inconsistent_solved$H2_REFCOPY_FINAL)] = inconsistent_solved$H1_REFCOPY_FINAL[is.na(inconsistent_solved$H2_REFCOPY_FINAL)]
        # combine data again
        haplo_final = rbind(regions_agree, regions_missing_asm, regions_missing_spa, inconsistent_solved, use.names=TRUE, fill=TRUE)
        return(haplo_final)
    }

    # Function to solve inconsistency of calls from assembly and reads spanning -- in use
    solve_inconsistency = function(i, inconsistent){
        # select entry of interest
        tmp = inconsistent[i, ]
        # extract alleles of asm and spa
        asm_alleles = as.numeric(c(tmp$H1_HAPLO_SIZE_ASM, tmp$H2_HAPLO_SIZE_ASM))
        spa_alleles = as.numeric(c(tmp$H1_HAPLO_SIZE_READS_SPANNING, tmp$H2_HAPLO_SIZE_READS_SPANNING))
        # an easy check if by looking at unique alleles (if homozygous, the 2 alleles will be the same, otherwise different). we can then prioritize when 2 alleles are present.
        n_unique_asm = length(unique(asm_alleles)); n_unique_spa = length(unique(spa_alleles))
        if (n_unique_asm > n_unique_spa){
            # more alleles in assembly than in reads-spanning
            tmp$DECISION = 'PRIORITY_ASSEMBLY'
            tmp$H1_HAPLO_FINAL = tmp$H1_HAPLO_SIZE_ASM; tmp$H2_HAPLO_FINAL = tmp$H2_HAPLO_SIZE_ASM
            tmp$H1_MOTIF_FINAL = tmp$H1_CONSENSUS_MOTIF_ASM; tmp$H2_MOTIF_FINAL = tmp$H2_CONSENSUS_MOTIF_ASM
            tmp$H1_COPIE_FINAL = tmp$H1_CONSENSUS_ASM; tmp$H2_COPIE_FINAL = tmp$H2_CONSENSUS_ASM
            tmp$H1_REFCOPY_FINAL = tmp$H1_REFERENCE_MOTIF_CN_ASM; tmp$H2_REFCOPY_FINAL = tmp$H2_REFERENCE_MOTIF_CN_ASM
        } else if (n_unique_asm < n_unique_spa){
            # more alleles in reads-spanning than in assembly
            tmp$DECISION = 'PRIORITY_READS_SPANNING'
            tmp$H1_HAPLO_FINAL = tmp$H1_HAPLO_SIZE_READS_SPANNING; tmp$H2_HAPLO_FINAL = tmp$H2_HAPLO_SIZE_READS_SPANNING
            tmp$H1_MOTIF_FINAL = tmp$H1_CONSENSUS_MOTIF_READS_SPANNING; tmp$H2_MOTIF_FINAL = tmp$H2_CONSENSUS_MOTIF_READS_SPANNING
            tmp$H1_COPIE_FINAL = tmp$H1_CONSENSUS_READS_SPANNING; tmp$H2_COPIE_FINAL = tmp$H2_CONSENSUS_READS_SPANNING
            tmp$H1_REFCOPY_FINAL = tmp$H1_REFERENCE_MOTIF_CN_READS_SPANNING; tmp$H2_REFCOPY_FINAL = tmp$H2_REFERENCE_MOTIF_CN_READS_SPANNING
        } else {
            # same number of alleles, but different
            tmp$DECISION = 'UNCERTAIN'
            tmp$H1_HAPLO_FINAL = NA; tmp$H2_HAPLO_FINAL = NA
            tmp$H1_MOTIF_FINAL = NA; tmp$H2_MOTIF_FINAL = NA
            tmp$H1_COPIE_FINAL = NA; tmp$H2_COPIE_FINAL = NA
            tmp$H1_REFCOPY_FINAL = NA; tmp$H2_REFCOPY_FINAL = NA
        }
        return(tmp)
    }

    # Function to create VCF file -- in use
    createVCFfile = function(data, reference, type, all_trf_nodup_coverage){
        # Create a data frame with the variant information -- this will be different for certain regions (also whether based on reads-spanning only or assembly only) and uncertain regions
        if (type == 'comparison'){
            # However, some info are common to all regions
            data$REGION_COMBINED = ifelse(is.na(data$REGION_ASM), data$REGION_READS_SPANNING, data$REGION_ASM)
            data$SAMPLE_COMBINED = ifelse(is.na(data$SAMPLE_ASM), data$SAMPLE_READS_SPANNING, data$SAMPLE_ASM)
            # adjust the sample names
            data$SAMPLE_COMBINED = str_replace_all(data$SAMPLE_COMBINED, '__Assembly_haps_aln.bam', '')
            data$SAMPLE_COMBINED = str_replace_all(data$SAMPLE_COMBINED, '__Assembly_p_ctg_aln.bam', '')
            data$SAMPLE_COMBINED = str_replace_all(data$SAMPLE_COMBINED, '__rawReads.bam', '')
        } else {
            data$REGION_COMBINED = data$REGION
            data$SAMPLE_COMBINED = data$SAMPLE
            # adjust the sample names
            data$SAMPLE_COMBINED = str_replace_all(data$SAMPLE_COMBINED, '__Assembly_haps_aln.bam', '')
            data$SAMPLE_COMBINED = str_replace_all(data$SAMPLE_COMBINED, '__Assembly_p_ctg_aln.bam', '')
            data$SAMPLE_COMBINED = str_replace_all(data$SAMPLE_COMBINED, '__rawReads.bam', '')
        }
        # select unique regions
        if ('reads-spanning' %in% reference$DATA_TYPE){ reference = reference[which(reference$DATA_TYPE == 'reads-spanning'),] }
        # Chromosome -- from reference
        chroms = str_split_fixed(reference$REGION, ':', 2)[, 1]
        # Position
        pos = as.numeric(str_split_fixed(str_split_fixed(reference$REGION, ':', 2)[, 2], '-', 2)[, 1])
        # Id is the REGION_COMBINED_COLUMN
        ids = reference$REGION
        # REF is the size in the reference genome
        ref_size = reference$H1_HAPLO_SIZE
        # ALT we set to . as there are many alleles possible, and we report them in the SAMPLE field
        alt_size = rep('.', length(ref_size))
        # QUAL we also set to .
        qual = rep('.', length(ref_size))
        # FILTER we set to PASS, then each sample will have a personal entry of quality
        filter = rep('PASS', length(ref_size))
        # INFO we set put the information from the reference genome (consensus motif and number of copies)
        info = paste(reference$H1_CONSENSUS_MOTIF, reference$H1_CONSENSUS, sep = ';')
        # FORMAT is important: for each sample, we put the filter (PASS_BOTH, PASS_ASM, PASS_RSP, UNCERTAIN), the sizes (H1_HAPLO|H2_HAPLO), MOTIF (H1_MOTIF|H2_MOTIF), COPIES (H1_CN|H2_CN)
        format = rep('QC;GT;MOTIF;CN;CN_REF', length(ref_size))
        # SAMPLE is the actual sample information, we use a function to extract that for each sample
        all_samples = unique(data$SAMPLE_COMBINED)
        samples_info = lapply(all_samples, extractInfoVCF, data = data, ids = ids)
        # FINALLY COMPILE THE DATAFRAME FOR VCF file
        df <- data.frame(CHROM = chroms, POS = pos, ID = ids, REF = ref_size, ALT = alt_size, QUAL = qual, FILTER = filter, INFO = info, FORMAT = format)
        # add the samples to the df
        for (i in seq_along(samples_info)) { df <- cbind(df, samples_info[[i]]); colnames(df)[ncol(df)] <- all_samples[i] }
        # write header
        header <- c("##fileformat=VCFv4.2", "##INFO=<ID=REFERENCE_INFO,Number=2,Type=String,Description=\"Motif observed in the reference genome (GRCh38), and relative number of motif repetitions.\">", "##FORMAT=<ID=QC,Number=1,Type=String,Description=\"Quality summary of TREAT genotyping. PASS_BOTH: genotype agreed between reads-spanning and assembly. PASS_RSP: priority to genotype from reads-spanning. PASS_ASM: priority to genotype from assembly.\">", "##FORMAT=<ID=GT,Number=2,Type=String,Description=\"Phased size of the tandem repeats. H1_size | H2_size\">", "##FORMAT=<ID=MOTIF,Number=2,Type=String,Description=\"Phased consensus motif found in the sample. H1_motif | H2_motif\">", "##FORMAT=<ID=CN,Number=2,Type=String,Description=\"Phased number of copies of the motif found in the sample. H1_copies | H2_copies\">", "##FORMAT=<ID=CN_REF,Number=2,Type=String,Description=\"Phased estimation of the reference motif as found in the sample. H1_motif_ref | H2_motif_ref\">", paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", paste(all_samples, collapse = "\t"), sep = "\t"))
        # Combine header and data frame
        vcf <- c(header, apply(df, 1, paste, collapse = "\t"))
        return(vcf)
    }

    # Function to extract sample-specific information -- in use
    extractInfoVCF = function(s, data, ids){
        # extract data for the sample
        tmp = data[which(data$SAMPLE_COMBINED == s),]
        ids_df = data.frame(REGION_COMBINED = ids)
        # match the order of the reference genome
        tmp = merge(ids_df, tmp, by = 'REGION_COMBINED', all.x = T)
        tmp_match <- tmp[order(match(tmp$REGION_COMBINED, ids)), ]

        # extract QC
        tmp_match$FILTER = ifelse(tmp_match$DECISION == 'OK_BOTH', 'PASS_BOTH', NA)
        tmp_match$FILTER[which(tmp_match$DECISION == 'PRIORITY_READS_SPANNING')] = 'PASS_RSP'
        tmp_match$FILTER[which(tmp_match$DECISION == 'PRIORITY_ASSEMBLY')] = 'PASS_ASM'
        tmp_match$FILTER[is.na(tmp_match$FILTER)] = 'UNCERTAIN'
        # extract GT
        tmp_match$GT = paste(tmp_match$H1_HAPLO_FINAL, tmp_match$H2_HAPLO_FINAL, sep = "|")
        # extract MOTIF
        tmp_match$MOTIF = paste(tmp_match$H1_MOTIF_FINAL, tmp_match$H2_MOTIF_FINAL, sep = "|")
        # extract COPIES
        tmp_match$CN = paste(tmp_match$H1_COPIE_FINAL, tmp_match$H2_COPIE_FINAL, sep = "|")
        # extract COPIES OF REFERENCE MOTIF
        tmp_match$CN_REF = paste(tmp_match$H1_REFCOPY_FINAL, tmp_match$H2_REFCOPY_FINAL, sep = "|")
        # combine all info
        combined_info = paste(tmp_match$FILTER, tmp_match$GT, tmp_match$MOTIF, tmp_match$CN, tmp_match$CN_REF, sep=";")
        return(combined_info)
    }

    # -----------------------------------------------------------------------------------------------------------------------------------------------------------------#

    # Function to make final decision between assembly and reads-spanning -- NOT in use, but maybe useful for debugging in the near future
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
    
    # Function to compare reads-spanning and assembly-based approached based on multiple processing -- NOT in use, but maybe useful for debugging in the near future
    compareReadsSpanning_Asm_mp = function(s, all_regions, all_haplo, deviation){
        comparison = list()
        for (r in all_regions){
            tmp = all_haplo[which(all_haplo$SAME_NAME == s & all_haplo$REGION == r),]
            # split assembly and reads-spanning
            asm = tmp[which(tmp$DATA_TYPE == 'assembly'),]; spa = tmp[which(tmp$DATA_TYPE == 'reads-spanning'),]
            # for the assembly, keep the haplotype-aware contigs if there are multiple entries
            if (nrow(asm) >1){ asm = asm[grep('Assembly_haps', asm$SAMPLE),] }
            # first check if we have asm and reads-spanning
            if (nrow(asm) >0 & nrow(spa) >0){
                # pairwise comparison of haplotypes. first extract data for asm and spa for h1 and h2
                asm_h1 = asm[, c('H1_CONSENSUS', 'H1_CONSENSUS_MOTIF', 'H1_HAPLO_SIZE')]; spa_h1 = spa[, c('H1_CONSENSUS', 'H1_CONSENSUS_MOTIF', 'H1_HAPLO_SIZE')]
                asm_h2 = asm[, c('H2_CONSENSUS', 'H2_CONSENSUS_MOTIF', 'H2_HAPLO_SIZE')]; spa_h2 = spa[, c('H2_CONSENSUS', 'H2_CONSENSUS_MOTIF', 'H2_HAPLO_SIZE')]
                # then check for NAs
                if (NA %in% c(asm_h1$H1_HAPLO_SIZE, asm_h2$H2_HAPLO_SIZE, spa_h1$H1_HAPLO_SIZE, spa_h2$H2_HAPLO_SIZE)){
                    # find the NA(s)
                    na_index = which(is.na(c(asm_h1$H1_HAPLO_SIZE, asm_h2$H2_HAPLO_SIZE, spa_h1$H1_HAPLO_SIZE, spa_h2$H2_HAPLO_SIZE)))
                    if (length(na_index) == 4){
                        colnames(asm) = paste0(colnames(asm), '_ASM'); asm$COMPARISON = 'missing_all_alleles'
                        colnames(spa) = paste0(colnames(spa), '_READS_SPANNING'); spa$H1_DIFF = NA; spa$H2_DIFF = NA
                    } else if (length(na_index) == 3){
                        good_one = which(!seq(1, 4) %in% na_index)
                        colnames(asm) = paste0(colnames(asm), '_ASM')
                        colnames(spa) = paste0(colnames(spa), '_READS_SPANNING')
                        spa$H1_DIFF = NA; spa$H2_DIFF = NA
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
                            diff_tmp = abs(as.numeric(nonmissing_asm[, 3]) - as.numeric(nonmissing_spa[, 3]))
                            if (diff_tmp <= max(1.5, as.numeric(nonmissing_asm[, 1])*deviation)){ check = 'both_miss_1_allele_other_ok' } else { check = 'both_miss_1_allele_other_different' }
                            colnames(asm) = paste0(colnames(asm), '_ASM'); asm$COMPARISON = check
                            colnames(spa) = paste0(colnames(spa), '_READS_SPANNING'); spa$H1_DIFF = diff_tmp; spa$H2_DIFF = NA
                        }
                    } else {
                        if (na_index %in% c(1, 2)){
                            if (na_index == 1){ nonmissing = asm_h2 } else { nonmissing = asm_h1 }
                            # calculate distance between alleles based on copy number
                            diff_alleles = abs(c(as.numeric(nonmissing[, 3]) - as.numeric(spa_h1$H1_HAPLO_SIZE), as.numeric(nonmissing[, 3]) - as.numeric(spa_h2$H2_HAPLO_SIZE)))
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
                    # make sure variables are numeric
                    spa_h1$H1_HAPLO_SIZE = as.numeric(spa_h1$H1_HAPLO_SIZE); spa_h2$H2_HAPLO_SIZE = as.numeric(spa_h2$H2_HAPLO_SIZE)
                    asm_h1$H1_HAPLO_SIZE = as.numeric(asm_h1$H1_HAPLO_SIZE); asm_h2$H2_HAPLO_SIZE = as.numeric(asm_h2$H2_HAPLO_SIZE)
                    # calculate distance between alleles: h1
                    diff_h1 = abs(c(asm_h1$H1_HAPLO_SIZE - spa_h1$H1_HAPLO_SIZE, asm_h1$H1_HAPLO_SIZE - spa_h2$H2_HAPLO_SIZE))
                    closest_h1 = which(diff_h1 == min(diff_h1))[1]; if (closest_h1 == 1){ closest_all = spa_h1 } else { closest_all = spa_h2 }
                    if (diff_h1[closest_h1] <= max(1.5, as.numeric(closest_all[, 3])*deviation)){ check_h1 = 'same' } else { check_h1 = 'different' }
                    # then h2, exclude the h1 allele
                    if (colnames(closest_all)[1] == 'H2'){ diff_h2 = abs(asm_h2$H2_HAPLO_SIZE - spa_h1$H1_HAPLO_SIZE) } else { diff_h2 = abs(asm_h2$H2_HAPLO_SIZE - spa_h2$H2_HAPLO_SIZE) }
                    closest_h2 = abs(closest_h1 - 1); if (closest_h2 == 1){ closest_all = spa_h1 } else { closest_all = spa_h2 }
                    if (diff_h2 <= max(1.5, as.numeric(closest_all[, 3])*deviation)){ check_h2 = 'same' } else { check_h2 = 'different' }
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
    #inp_pha = args$phase; inp_pha = unlist(strsplit(inp_pha, ','))
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
        #cat(paste0('**** input PHASING selected -> ', inp_pha, '\n'))
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
            all_trf = rbindlist(all_trf, use.names=TRUE)
        }
  
    # 2. Let's adjust the motifs -- generate a consensus -- essentially merging the same motifs together at the single read level -- this for both analyses
        cat('****** Merging similar motifs together\n')
        all_motifs = data.frame(motif = unique(all_trf$TRF_MOTIF), stringsAsFactors = F)
        all_motifs = all_motifs[!is.na(all_motifs$motif),]
        # alternative way using multiprocessing
        main_motifs = unlist(mclapply(all_motifs, mergeMotif_mp, mc.cores = n_cpu))
        all_motifs = data.frame(motif = all_motifs, UNIFORM_MOTIF = main_motifs)
        all_trf = merge(all_trf, all_motifs, by.x = 'TRF_MOTIF', by.y = 'motif', all.x = T)
        # to make sure there are no conflicts, change the read name
        all_trf$UNIQUE_NAME = paste(all_trf$READ_NAME, all_trf$SAMPLE_NAME, all_trf$REGION, sep = "___")

    # 3. add the TR size
        # change two colnames
        colnames(all_trf)[which(colnames(all_trf) == 'HAPLOTAG')] = 'HAPLOTYPE'
        colnames(all_trf)[which(colnames(all_trf) == 'LEN_SEQUENCE_FOR_TRF')] = 'LENGTH_SEQUENCE'
        all_trf$TR_LENGHT = all_trf$LENGTH_SEQUENCE

    # 4. good to exclude duplicated reads otherwise results would be biased towards sequences with more complex motifs (where TRF finds multiple matches)
        all_trf_nodup = all_trf[!duplicated(all_trf$UNIQUE_NAME),]; dups = all_trf[duplicated(all_trf$UNIQUE_NAME),]
        # also extract the number of reads per sample, per region, per data type
        all_trf_nodup_coverage = data.frame(table(all_trf_nodup$SAMPLE_NAME, all_trf_nodup$REGION, all_trf_nodup$DATA_TYPE)); colnames(all_trf_nodup_coverage) = c('SAMPLE_NAME', 'REGION', 'DATA_TYPE', 'COVERAGE')

    # 5. we do genotyping on the sizes
        cat('****** Genotyping TRs\n')
        # temporarily disable warnings
        defaultW <- getOption("warn"); options(warn = -1)
        # split reads-spanning and assembly-based
        reads_span = all_trf_nodup[which(all_trf_nodup$DATA_TYPE == 'reads-spanning'),]; asm = all_trf_nodup[which(all_trf_nodup$DATA_TYPE == 'assembly'),]
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

    # 6. now let's bring the duplicates in again and assign the correct haplotype based on the other duplicate
        # this is also now done using multiple processors
        dups_name = unique(dups$UNIQUE_NAME)
        tmp_df_dups = data.frame(COPIES_TRF = dups$COPIES_TRF, UNIFORM_MOTIF = dups$UNIFORM_MOTIF, REGION = dups$REGION, PASSES = dups$PASSES, READ_QUALITY = dups$READ_QUALITY, LENGTH_SEQUENCE = dups$LENGTH_SEQUENCE, READ_NAME = dups$READ_NAME, START_TRF = dups$START_TRF, END_TRF = dups$END_TRF, TRF_PERC_MATCH = dups$TRF_PERC_MATCH, TRF_PERC_INDEL = dups$TRF_PERC_INDEL, DATA_TYPE = dups$DATA_TYPE, UNIQUE_NAME = dups$UNIQUE_NAME)
        # add information about type, sample, haplo_value, polished_reads and polished_haplo_values
        tmp_info = all_res[, c('UNIQUE_NAME', 'type', 'sample', 'haplo_value', 'polished_reads', 'polished_haplo_values', 'HAPLOTYPE')]
        tmp_df_dups = merge(tmp_df_dups, tmp_info, by = 'UNIQUE_NAME')
        # finally combine the duplicate information with the unique information
        all_res_combined = rbind(all_res, tmp_df_dups)
        #saved_all = all_res_combined

    # 7. now we should look at the motifs -- implemented parallel computing
        cat('****** Generating consensus motifs\n')
        all_samples = unique(all_res$sample); all_regions = unique(all_res$REGION); motif_res = list()
        # first run on the reference genome
        motif_res_reference = generateConsens_mp(s = 'reference', all_regions, all_res = all_res_combined, motif_res_reference = NA)
        motif_res = rbindlist(mclapply(all_samples[which(all_samples != 'reference')], generateConsens_mp, all_regions = all_regions, all_res = all_res_combined, motif_res_reference = motif_res_reference, mc.cores = n_cpu), use.names=TRUE)
        motif_res = rbind(motif_res, motif_res_reference)
        # add actual sequence to motif_res
        raw_seqs = all_trf[, c('READ_NAME', 'SEQUENCE_WITH_PADDINGS')]; raw_seqs = raw_seqs[!duplicated(raw_seqs$READ_NAME),]
        motif_res = merge(motif_res, raw_seqs, by = 'READ_NAME')

    # 8. finally call haplotypes -- implemented parallel computing
        cat('****** Haplotype calling\n')
        all_samples = unique(motif_res$sample); all_regions = unique(motif_res$REGION)
        # haplotype calling should be done in reads-spanning and assembly separately, otherwise if the haplotypes are different there will be problems
        if (anal_type == 'reads-spanning + assembly + comparison'){
            # reads spanning
            all_haplo_rs = rbindlist(mclapply(all_samples, callHaplo_mp, all_regions = all_regions, data = motif_res[which(motif_res$DATA_TYPE == 'reads-spanning'),], mc.cores = n_cpu))
            # assembly
            all_haplo_asm = rbindlist(mclapply(all_samples, callHaplo_mp, all_regions = all_regions, data = motif_res[which(motif_res$DATA_TYPE == 'assembly'),], mc.cores = n_cpu))
            # combine them
            all_haplo = rbind(all_haplo_rs, all_haplo_asm)
        } else {
            all_haplo = rbindlist(mclapply(all_samples, callHaplo_mp, all_regions = all_regions, data = motif_res, mc.cores = n_cpu))
        }
        # add a unique identifier to match samples reads-spanning and assembly results
        all_haplo$SAME_NAME = str_split_fixed(all_haplo$SAMPLE, '__', 2)[, 1]
        # exclude reference
        reference = all_haplo[which(all_haplo$SAME_NAME == 'reference'),]; all_haplo = all_haplo[which(all_haplo$SAME_NAME != 'reference'),]

    # 9. if both reads-spanning and assembly were submitted, we should compare them and make a unified call
        if (anal_type == 'reads-spanning + assembly + comparison'){
            cat('****** Comparing reads-spanning with assembly\n')
            # run multiprocessing for each sample independently
            all_samples = unique(all_haplo$SAME_NAME); all_regions = unique(all_haplo$REGION)
            all_haplo_annotated = rbindlist(mclapply(all_samples, comparison_faster, all_regions = all_regions, all_haplo = all_haplo, deviation = thr_mad, mc.cores = n_cpu))
            # then make a final decision
            all_haplo_final = makeDecision_faster(all_haplo_annotated, n_cpu)
            cat('****** Generating tables\n')
            # finally we need to make VCFs for compatibility
            vcf = createVCFfile(all_haplo_final, reference, type='comparison', all_trf_nodup_coverage)
            # finally also save this dataset
            out_dir = str_replace_all(out_dir, '/$', ''); if (!dir.exists(out_dir)){ system(paste0('mkdir ', out_dir)) }
            write.table(vcf, paste0(out_dir, '/samples_genotypes.vcf'), quote=F, row.names = F, col.names = F, sep = '\t')
            # also output the single-reads results
            write.table(motif_res, paste0(out_dir, '/haplotyping_single_reads.txt'), quote=F, row.names = F, sep = '\t')
        } else {
            # otherwise, do not do the comparison, generate the tables and exit
            cat('****** Generating tables\n')
            out_dir = str_replace_all(out_dir, '/$', ''); if (!dir.exists(out_dir)){ system(paste0('mkdir ', out_dir)) }
            # adjust for making the vcf
            all_haplo$DECISION = ifelse(anal_type == 'reads_spanning', 'PRIORITY_READS_SPANNING', 'PRIORITY_ASSEMBLY')
            # fix NA
            all_haplo[all_haplo == 'NA'] = NA
            all_haplo$H2_HAPLO_SIZE[is.na(all_haplo$H2_HAPLO_SIZE)] = all_haplo$H1_HAPLO_SIZE[is.na(all_haplo$H2_HAPLO_SIZE)]
            all_haplo$H2_SPECIFIC_MOTIF[is.na(all_haplo$H2_SPECIFIC_MOTIF)] = all_haplo$H1_SPECIFIC_MOTIF[is.na(all_haplo$H2_SPECIFIC_MOTIF)]
            all_haplo$H2_CONSENSUS[is.na(all_haplo$H2_CONSENSUS)] = all_haplo$H1_CONSENSUS[is.na(all_haplo$H2_CONSENSUS)]
            all_haplo$H2_CONSENSUS_MOTIF[is.na(all_haplo$H2_CONSENSUS_MOTIF)] = all_haplo$H1_CONSENSUS_MOTIF[is.na(all_haplo$H2_CONSENSUS_MOTIF)]
            all_haplo$H2_REFERENCE_MOTIF_CN[is.na(all_haplo$H2_REFERENCE_MOTIF_CN)] = all_haplo$H1_REFERENCE_MOTIF_CN[is.na(all_haplo$H2_REFERENCE_MOTIF_CN)]
            # set the right columns
            colnames(all_haplo) = c('SAMPLE', 'REGION', 'H1_COPIE_FINAL', 'H2_COPIE_FINAL', 'H1_MOTIF_FINAL', 'H2_MOTIF_FINAL', 'H1_HAPLO_FINAL', 'H2_HAPLO_FINAL', 'H1_SPECIFIC_MOTIF', 'H2_SPECIFIC_MOTIF', 'H1_ADDITIONAL_MOTIF', 'H2_ADDITIONAL_MOTIF', 'REFERENCE_MOTIF', 'H1_REFCOPY_FINAL', 'H2_REFCOPY_FINAL', 'EXCLUDED', 'EXCLUDED_LEN', 'DATA_TYPE', 'SAME_NAME', 'DECISION')
            vcf = createVCFfile(all_haplo, reference, type='single', all_trf_nodup_coverage)
            # write outputs
            write.table(vcf, paste0(out_dir, '/samples_genotypes.vcf'), quote=F, row.names = F, col.names = F, sep = '\t')
            write.table(motif_res, paste0(out_dir, '/haplotyping_single_reads.txt'), quote=F, row.names = F, sep = '\t')
        }
    cat('**** Analysis done!!\n')
