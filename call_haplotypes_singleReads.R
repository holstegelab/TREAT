#!/usr/bin/R

# Libraries
    library(plyr)
    library(data.table)
    library(argparse)

# Functions
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
                    # run k-means using as k the ploidy value
                    kmeans_res <- kmeans(x = reads, centers = ploidy, iter.max = 100)
                    # take centers of the k-means
                    centers_kmeans <- kmeans_res$centers[order(kmeans_res$centers)]
                    # check whether there is support for all clusters
                    tmp_df = data.frame(reads = reads, cluster = kmeans_res$cluster)
                    tb_frq = as.data.frame(table(tmp_df$cluster)); tb_frq$Prop = tb_frq$Freq/sum(tb_frq$Freq); tb_frq  = tb_frq[order(tb_frq$Freq),]
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

    # Polish haplotypes after phasing is done
    polishHaplo_afterPhasing = function(phased_data, thr_mad){
        # split h1 from h2, reads and excluded
        h1 = phased_data[which(phased_data$HAPLOTYPE == 1),]; h2 = phased_data[which(phased_data$HAPLOTYPE == 2),]
        excluded = phased_data[is.na(phased_data$HAPLOTYPE),]; if (nrow(excluded) >=1) { excluded$polished_reads = 'exclude' }
        # assign haplotypes before polishing
        if (nrow(h1) >=1){ h1$haplo_value = median(h1$COPIES_TRF) }
        if (nrow(h2) >=1){ h2$haplo_value = median(h2$COPIES_TRF) }
        # using the thr_mad, split reads to keep from those to exclude
        h1_thr_mad = max(1.5, ceiling(median(h1$COPIES_TRF)*thr_mad))
        h2_thr_mad = max(1.5, ceiling(median(h2$COPIES_TRF)*thr_mad))
        h1$polished_reads = ifelse(abs(h1$COPIES_TRF - median(h1$COPIES_TRF)) <= h1_thr_mad, 'keep', 'exclude')
        h2$polished_reads = ifelse(abs(h2$COPIES_TRF - median(h2$COPIES_TRF)) <= h2_thr_mad, 'keep', 'exclude')
        # assign haplotypes after polishing
        if (nrow(h1) >= 1){ h1$polished_haplo_values = median(h1$COPIES_TRF[which(h1$polished_reads == 'keep')]) }
        if (nrow(h2) >= 1){ h2$polished_haplo_values = median(h2$COPIES_TRF[which(h2$polished_reads == 'keep')]) }
        # return combined object
        res = rbind.fill(h1, h2)
        # clean haplotypes -- sometimes there's haplotype 2 but not 1 --> make it 1
        if (2 %in% res$HAPLOTYPE && (!(1 %in% res$HAPLOTYPE))){ res$HAPLOTYPE = 1 }
        res = plyr::rbind.fill(res, excluded)
        return(res)
    }

    # Function to perform guided haplotyping using phasing information
    PhasingBased_haplotyping <- function(reads_df, sample_name, thr_mad){
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
                    res = KMeansBased_haplotyping(reads = nonPhased$COPIES_TRF, thr = 2, min_support = 1, thr_mad_orig = 0.015, type = 'single_sample')
                    # then polish with the polisher without phasing
                    res_polished = polishHaplo_noPhasing(res, 0.015)
                    # extract polished reads
                    reads_h1 = unlist(strsplit(as.character(res_polished$reads_h1), ',')); reads_h2 = unlist(strsplit(as.character(res_polished$reads_h2), ','))
                    reads_df$type = NA
                    # make the final df with reads and haplotypes (depending on how many alleles were found)
                    if (reads_h2 == 'homozygous' && !is.na(reads_h2)){
                        reads_df$type[which(reads_df$COPIES_TRF %in% reads_h1)] = 'Assigned'; reads_df$HAPLOTYPE[which(reads_df$COPIES_TRF %in% reads_h1)] = 1
                    } else if (!is.na(reads_h2)){
                        reads_df$type[which(reads_df$COPIES_TRF %in% c(reads_h1, reads_h2))] = 'Assigned'
                        reads_df$HAPLOTYPE[which(reads_df$COPIES_TRF %in% reads_h1)] = 1; reads_df$HAPLOTYPE[which(reads_df$COPIES_TRF %in% reads_h2)] = 2
                    }
                    # the final df is all_res
                    all_res = reads_df
                # there is phasing information available, but only haplo 1 was found
                } else if (length(unique(phased$HAPLOTYPE)) == 1) {
                    reads_df$type = NA; reads_df$type[!is.na(reads_df$HAPLOTYPE)] <- 'Phased'
                    # gather all reads within thr_mad from the haplotype
                    median_cn = median(phased$COPIES_TRF)
                    thr_mad = max(1.5, ceiling(median(median_cn)*thr_mad))
                    # idea: look at reads without phase and try to assign them a phase based on thr_mad 
                    reads_df$type[which(abs(reads_df$COPIES_TRF - median_cn) <= thr_mad)] <- 'Assigned'
                    reads_df$HAPLOTYPE[which(abs(reads_df$COPIES_TRF - median_cn) <= thr_mad)] <- 1
                    # the remaining reads (those with mad > thr_mad) are haplotype 2
                    reads_df$HAPLOTYPE[is.na(reads_df$HAPLOTYPE)] = 2
                    reads_df$type[is.na(reads_df$HAPLOTYPE)] = 'Assigned'
                    # the final df is all_res
                    all_res = reads_df
                # if there is phasing information for both alleles
                } else {
                    phased$type = 'Phased'
                    # find centers (mean of copy number in the group of same haplotype)
                    h1 = mean(phased$COPIES_TRF[which(phased$HAPLOTYPE == 1)])
                    h2 = mean(phased$COPIES_TRF[which(phased$HAPLOTYPE == 2)])
                    # idea: assign reads without phase based on distance from those with a phase --> lowest distance gets the assignment
                    for (i in 1:nrow(nonPhased)){
                        tmp = nonPhased$COPIES_TRF[i]; closest_haplo = ifelse(abs(tmp - h1) < abs(tmp - h2), 1, 2)
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

    # Function to clean TRF matches per read based on percentage of coverage and percentage of TRF match
    cleanMotifsTRF <- function(entire_set){
        # split exact matches and non-exact matches based on motif match type
        exact_match = entire_set[which(entire_set$MATCH_TYPE %in% c('perm_motif', 'exact_motif')),]
        non_exact_matches = entire_set[which(entire_set$MATCH_TYPE %in% c('different_length', 'different_motif')),]
        # calculate percentage of sequence covered by TRF, for exact and non-exact matches
        exact_match$PERC_SEQ_COVERED_TR = (nchar(exact_match$PADDING_AFTER) + nchar(exact_match$PADDING_BEFORE) + nchar(exact_match$SEQUENCE_TRF)) / exact_match$LENGTH_SEQUENCE
        non_exact_matches$PERC_SEQ_COVERED_TR = (nchar(non_exact_matches$PADDING_AFTER) + nchar(non_exact_matches$PADDING_BEFORE) + nchar(non_exact_matches$SEQUENCE_TRF)) / non_exact_matches$LENGTH_SEQUENCE
            #ggplot(exact_match, aes(x = PERC_SEQ_COVERED_TR)) + geom_density()
            #plot(exact_match$PERC_SEQ_COVERED_TR)
        # split well covered samples from non-well covered (EXACT MATCH)
        well_covered = exact_match[which(exact_match$PERC_SEQ_COVERED_TR >= 0.95),]
        missing_cove = exact_match[which(exact_match$PERC_SEQ_COVERED_TR < 0.95),]
        
        # look at DUPLICATED READS poorly covered reads
        # these are reads with (1) multiple trf matches exact-motif 
            #print(paste0(length(well_covered$READ_NAME[duplicated(well_covered$READ_NAME)]), ' duplicates in well-covered reads'))
            #print(paste0(length(missing_cove$READ_NAME[duplicated(missing_cove$READ_NAME)]), ' duplicates in poorly-covered reads'))
        poorly_covered_duplicates = missing_cove[which(missing_cove$READ_NAME %in% missing_cove$READ_NAME[duplicated(missing_cove$READ_NAME)]),]
        results_poor_duplicates = list()
        # loop on duplicated reads to see if they are genuinely measuring 2 motifs or not
        for (d in unique(poorly_covered_duplicates$READ_NAME)){
            tmp = poorly_covered_duplicates[which(poorly_covered_duplicates$READ_NAME == d), ]; overlaps = list()
            # loop on each of the duplicates and measure the overlap between them
            for (i in 1:nrow(tmp)){ overlaps[[(length(overlaps) + 1)]] = seq(tmp$START_TRF[i], tmp$END_TRF[i]) }
            # then calculate joined coverage and overlaps
            join_coverage = unlist(overlaps); overlaps = join_coverage[duplicated(join_coverage)]
            # calculate cumulative coverage (50 is the estimated padding)
            cum_coverage = (length(join_coverage) + 50 - length(overlaps)) / tmp$LENGTH_SEQUENCE[1]
            # check if there was a gain in doing this in terms of cumulative coverage compared to the original coverages
            if (max(cum_coverage, tmp$PERC_SEQ_COVERED_TR) == cum_coverage){
                # if there is a gain, probably we're dealing with a true variable motif
                tmp$EXITUS = 'true_variable_motif'; tmp$CUM_PERC_MATCH = cum_coverage
            } else {
                # if there is no gain, probably we're dealing with a probable alternative motif
                tmp$EXITUS = 'probable_alternative_motif'; tmp$CUM_PERC_MATCH = NA
            }
            results_poor_duplicates[[(length(results_poor_duplicates) + 1)]] = tmp
        }
        # at the end, combine all duplicates-checks together
        results_poor_duplicates = data.table::rbindlist(results_poor_duplicates, fill=T)
        # merge poorly-covered (duplicated) with well-covered
        well_covered$EXITUS = 'true_exact_motif'; well_covered$CUM_PERC_MATCH = well_covered$PERC_SEQ_COVERED_TR
        well_covered = rbind(well_covered, results_poor_duplicates)  
        # exclude duplicates from the poorly covered reads
        missing_cove = missing_cove[which(!(missing_cove$READ_NAME %in% well_covered$READ_NAME)),]
        print(paste0(length(missing_cove$READ_NAME[duplicated(missing_cove$READ_NAME)]), ' duplicates in poorly-covered reads now'))
        
        # (2) can be exact-matches with missing coverage + non-exact matches 
            #print(paste0(length(which(unique(missing_cove$READ_NAME) %in% non_exact_matches$READ_NAME)), ' reads poorly covered have non-exact matches'))
        results_poor_nonexact = list()
        # loop on poorly covered with non-exact matches
        for (d in unique(missing_cove$READ_NAME)){
            tmp = plyr::rbind.fill(missing_cove[which(missing_cove$READ_NAME == d), ], non_exact_matches[which(non_exact_matches$READ_NAME == d),])
            overlaps = list()
            # loop on each of the matches and measure the overlap between them
            for (i in 1:nrow(tmp)){ overlaps[[(length(overlaps) + 1)]] = seq(tmp$START_TRF[i], tmp$END_TRF[i]) }
            # calculate joined coverage and overlaps
            join_coverage = unlist(overlaps); overlaps = join_coverage[duplicated(join_coverage)]
            # calculate cumulative coverage (50 is the estimated padding)
            cum_coverage = (length(join_coverage) + 50 - length(overlaps)) / tmp$LENGTH_SEQUENCE[1]
            # check if there was a gain in doing this
            if (max(cum_coverage, na.omit(tmp$PERC_SEQ_COVERED_TR)) == cum_coverage){
                # if there is a gain, probably we're dealing with a true variable motif
                tmp$EXITUS = 'true_variable_motif'; tmp$CUM_PERC_MATCH = cum_coverage
            } else {
                # if there is no gain, probably we're dealing with a probable alternative motif
                tmp$EXITUS = 'probable_alternative_motif'; tmp$CUM_PERC_MATCH = NA
            }
            results_poor_nonexact[[(length(results_poor_nonexact) + 1)]] = tmp  
        }
        # at the end, combine all duplicates-checks together
        results_poor_nonexact = data.table::rbindlist(results_poor_nonexact, fill=T)
        # merge results poor covered + non-exact matches with exact results
        well_covered = rbind(well_covered, results_poor_nonexact)
        non_exact_matches = non_exact_matches[which(!(non_exact_matches$READ_NAME %in% missing_cove$READ_NAME)),]
        # exclude these non-exact matches from poorly covered reads
        missing_cove = missing_cove[which(!(missing_cove$READ_NAME %in% well_covered$READ_NAME)),]

        # (3) or can be non-exact matches + exact-matches
        results_nonexact = list()
        # loop on non-exact matches
        for (d in unique(non_exact_matches$READ_NAME)){
            # given a read, take exact-matches and non-exact matches for that read
            tmp_exact = well_covered[which(well_covered$READ_NAME == d),]
            tmp_non_exact = non_exact_matches[which(non_exact_matches$READ_NAME == d),]
            tmp_both = rbind.fill(tmp_exact, tmp_non_exact)
            # how to keep only genuine trf match? look at coverage of sequence, and if both coverage is 1 --> keep the one with highest PC_MATCH
            if (nrow(tmp_exact) == 1 & nrow(tmp_non_exact) == 1){
                if (tmp_exact$PERC_SEQ_COVERED_TR == 1 & tmp_non_exact$PERC_SEQ_COVERED_TR == 1){
                    # 1 exact match and 1 non-exact match for the same read
                    if (tmp_exact$PC_MATCH_TRF > tmp_non_exact$PC_MATCH_TRF){
                        # exact match cover more TR than non exact match --> likely wrong motif
                        tmp_non_exact$EXITUS = 'likely_wrong_motif'; tmp_non_exact$CUM_PERC_MATCH = tmp_non_exact$PERC_SEQ_COVERED_TR
                    } else {
                        # non exact match has higher PC_MATCH than exact match
                        tmp_non_exact$EXITUS = 'likely_wrong_motif'; tmp_non_exact$CUM_PERC_MATCH = tmp_non_exact$PERC_SEQ_COVERED_TR
                    }
                } else {
                    # exact match of non-exact match do not cover the whole TR
                    if (tmp_exact$PERC_SEQ_COVERED_TR == 1){
                        # exact match does --> non-exact match is likely a wrong motif
                        tmp_non_exact$EXITUS = 'likely_wrong_motif'; tmp_non_exact$CUM_PERC_MATCH = tmp_non_exact$PERC_SEQ_COVERED_TR
                    } else {
                        # non-exact match covers but exact matcch doesn't --> never happened thus far
                        cat(paste0('**** something is weird with percentage of coverage of non-exact of ', d, '\n'))
                    }
                }
            } else if (nrow(tmp_exact) >1){
                # if there are multiple exact matches, probably the trf match is correct
                tmp_non_exact$EXITUS = 'true_alternative_motif'; tmp_non_exact$CUM_PERC_MATCH = tmp_non_exact$PERC_SEQ_COVERED_TR
            } else if (nrow(tmp_exact) == 0 & nrow(tmp_non_exact) == 1){
                # no exact matches --> likely true alternative motif
                tmp_non_exact$EXITUS = 'true_alternative_motif'
                tmp_non_exact$CUM_PERC_MATCH = tmp_non_exact$PERC_SEQ_COVERED_TR
            } else if (nrow(tmp_non_exact) > 1){
                # no exact matches but multiple non-exact matches --> check overlap and gain in keeping both matches
                overlaps = list(); for (i in 1:nrow(tmp_non_exact)){ overlaps[[(length(overlaps) + 1)]] = seq(tmp_non_exact$START_TRF[i], tmp_non_exact$END_TRF[i]) }
                # calculate joined coverage and overlaps
                join_coverage = unlist(overlaps); overlaps = join_coverage[duplicated(join_coverage)]
                # calculate cumulative coverage (50 is the estimated padding)
                cum_coverage = (length(join_coverage) + 50 - length(overlaps)) / tmp_non_exact$LENGTH_SEQUENCE[1]
                # check if there was a gain in doing this
                if (max(cum_coverage, na.omit(tmp_non_exact$PERC_SEQ_COVERED_TR)) == cum_coverage){
                    tmp_non_exact$EXITUS = 'true_variable_motif'; tmp_non_exact$CUM_PERC_MATCH = cum_coverage
                } else {
                    tmp_non_exact = tmp_non_exact[order(-tmp_non_exact$PERC_SEQ_COVERED_TR),]
                    tmp_non_exact$EXITUS[1] = 'true_alternative_motif'
                    tmp_non_exact$CUM_PERC_MATCH[1] = cum_coverage
                    tmp_non_exact$EXITUS[2:nrow(tmp_non_exact)] = rep('likely_wrong_motif', nrow(tmp_non_exact)-1)
                    tmp_non_exact$CUM_PERC_MATCH[2:nrow(tmp_non_exact)] = rep(NA, nrow(tmp_non_exact)-1)
                }
            } else { 
                cat(paste0('**** somthing is weird for ', d, '\n')) 
            }
            results_nonexact[[(length(results_nonexact) + 1)]] = tmp_non_exact
        }
        results_nonexact = rbindlist(results_nonexact, fill=T)
  
        # split non-exact results based on what we did: exclude likely_wrong_motif
        non_exact_tokeep = results_nonexact[which(results_nonexact$EXITUS != 'likely_wrong_motif'),]
        well_covered = rbind(well_covered, non_exact_tokeep)
        non_exact_toexclude = results_nonexact[which(results_nonexact$EXITUS == 'likely_wrong_motif'),]
  
        # make list of results
        res = list(well_covered, non_exact_toexclude)
        return(res)
    }

    # Function to call haplotypes and combine information
    callHaplo <- function(data){
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
                            h1_info = unique(tmp$polished_haplo_values[which(tmp$HAPLOTYPE == 1)])
                            h1_motif = paste(unique(tmp$UNIFORM_MOTIF[which(tmp$HAPLOTYPE == 1)]), collapse = ',')
                            h1_var = sd(tmp$COPIES_TRF[which(tmp$HAPLOTYPE == 1)]) / mean(tmp$COPIES_TRF[which(tmp$HAPLOTYPE == 1)])
                            h1_len = paste(unique(tmp$LENGTH_SEQUENCE[which(tmp$HAPLOTYPE == 1)]), collapse = ',')
                            h2_info = unique(tmp$polished_haplo_values[which(tmp$HAPLOTYPE == 2)])
                            h2_motif = paste(unique(tmp$UNIFORM_MOTIF[which(tmp$HAPLOTYPE == 2)]), collapse = ',')
                            h2_var = sd(tmp$COPIES_TRF[which(tmp$HAPLOTYPE == 2)]) / mean(tmp$COPIES_TRF[which(tmp$HAPLOTYPE == 2)])
                            h2_len = paste(unique(tmp$LENGTH_SEQUENCE[which(tmp$HAPLOTYPE == 2)]), collapse = ',')
                        # otherwise if only h1 is present, only gather h1 info
                        } else if (1 %in% tmp$HAPLOTYPE){
                            h1_info = unique(tmp$polished_haplo_values[which(tmp$HAPLOTYPE == 1)])
                            h1_motif = paste(unique(tmp$UNIFORM_MOTIF[which(tmp$HAPLOTYPE == 1)]), collapse = ',')
                            h1_var = sd(tmp$COPIES_TRF[which(tmp$HAPLOTYPE == 1)]) / mean(tmp$COPIES_TRF[which(tmp$HAPLOTYPE == 1)])
                            h1_len = paste(unique(tmp$LENGTH_SEQUENCE[which(tmp$HAPLOTYPE == 1)]), collapse = ',')
                            h2_info = NA; h2_motif = NA; h2_var = NA; h2_len = NA
                        # otherwise gather only h2 info
                        } else {
                            h2_info = unique(tmp$polished_haplo_values[which(tmp$HAPLOTYPE == 2)])
                            h2_motif = paste(unique(tmp$UNIFORM_MOTIF[which(tmp$HAPLOTYPE == 2)]), collapse = ',')
                            h2_var = sd(tmp$COPIES_TRF[which(tmp$HAPLOTYPE == 2)]) / mean(tmp$COPIES_TRF[which(tmp$HAPLOTYPE == 2)])
                            h2_len = paste(unique(tmp$LENGTH_SEQUENCE[which(tmp$HAPLOTYPE == 2)]), collapse = ',')
                            h1_info = NA; h1_motif = NA; h1_var = NA; h1_len = NA
                        }
                        # then check excluded reads
                        if (nrow(excl) >0){ 
                            excl_info = paste(paste0(excl$UNIFORM_MOTIF, '(', excl$COPIES_TRF, ')'), collapse = '_')
                            excl_len = paste(excl$LENGTH_SEQUENCE, collapse = ',')
                        } else { 
                            excl_info = NA; excl_len = NA
                        }
                        # finally save a df with all info regarding TR
                        tmp_df = data.frame(SAMPLE = unique(tmp$sample), REGION = unique(tmp$REGION), H1 = h1_info, H2 = h2_info, 
                            H1_MOTIF = h1_motif, H2_MOTIF = h2_motif, H1_LENGTH = h1_len, H2_LENGTH = h2_len, COEF_VAR_H1 = h1_var, COEF_VAR_H2 = h2_var,
                            EXCLUDED = excl_info, EXCLUDED_LEN = excl_len)
                        # and add to the growing list
                        all_haplo[[(length(all_haplo) + 1)]] = tmp_df
                    # if haplotypes were not identified (i.e they are NAs)
                    } else {
                        excl_info = paste(paste0(excl$UNIFORM_MOTIF, '(', excl$COPIES_TRF, ')'), collapse = '_')
                        excl_len = paste(excl$LENGTH_SEQUENCE, collapse = ',')
                        # create the same df with NA values and add to growing list
                        tmp_df = data.frame(SAMPLE = unique(excl$sample), REGION = unique(excl$REGION), H1 = NA, H2 = NA, 
                            H1_MOTIF = NA, H2_MOTIF = NA, H1_LENGTH = NA, H2_LENGTH = NA, COEF_VAR_H1 = NA, COEF_VAR_H2 = NA,
                            EXCLUDED = excl_info, EXCLUDED_LEN = excl_len)
                        all_haplo[[(length(all_haplo) + 1)]] = tmp_df
                    }
                }
            }
        }
        # convert list to df at the end
        all_haplo = rbindlist(all_haplo)
        return(all_haplo)
    }

    # Function to call merge TRF calls when multiple matches exist for the same read
    mergeMotifs = function(data){
        # isolate reads with alternative motif
        var_motifs = data[which(data$EXITUS == 'true_variable_motif'),]
        exact_motif = data[which(data$EXITUS != 'true_variable_motif'),]
        # normally, these reads should all have a duplicate
        single_reads = unique(var_motifs$READ_NAME[duplicated(var_motifs$READ_NAME)])
        combined_res = list()
        # loop on the duplicated reads
        for (read in single_reads){
            # take both duplicates
            tmp = var_motifs[which(var_motifs$READ_NAME == read),]; tmp = tmp[order(tmp$START_TRF),]
            # split (for simplicity) simple cases (2 dups only) or more complex cases (>2 dups)
            # below more than 2 duplicates (that means, TRF found 3 matches in the same read)
            if (nrow(tmp) >2){
                # ind is the index to loop across duplicates
                ind = 1; multi_results = list()
                # loop across all duplicates
                while (ind <= (nrow(tmp)-1)){
                    # check if there are results (that is, if it is not the first iteration)
                    if (length(multi_results) >0){
                        # estimate start and end of line x-1 and line x
                        orig_start = as.numeric(multi_results[[ind-1]]$start); orig_end = as.numeric(multi_results[[ind-1]]$end); next_start = tmp$START_TRF[ind+1]; next_end = tmp$END_TRF[ind+1]
                        # motif length of line x
                        mot_len = tmp$LENGTH_MOTIF_TRF[ind]
                        # motif and copies of line x-1
                        unif_motif_orig = multi_results[[ind-1]]$combined_mot; copies_orig = as.numeric(multi_results[[ind-1]]$copies)
                        # motif and copies of line x+1
                        unif_motif_next = tmp$UNIFORM_MOTIF[ind+1]; copies_next = tmp$COPIES_TRF[ind+1]
                        # get sequence of line x
                        seq_orig = tmp$SEQUENCE_WITH_WINDOW[ind]
                    # otherwise at the first iteration
                    } else {
                        # estimate start and end of line x and x+1
                        orig_start = tmp$START_TRF[ind]; orig_end = tmp$END_TRF[ind]; next_start = tmp$START_TRF[ind+1]; next_end = tmp$END_TRF[ind+1]
                        # motif length of line x 
                        mot_len = tmp$LENGTH_MOTIF_TRF[ind]
                        # motif and copies of line x
                        unif_motif_orig = tmp$UNIFORM_MOTIF[ind]; copies_orig = tmp$COPIES_TRF[ind]
                        # motif and copies of line x+1
                        unif_motif_next = tmp$UNIFORM_MOTIF[ind+1]; copies_next = tmp$COPIES_TRF[ind+1]
                        # sequence of line x
                        seq_orig = tmp$SEQUENCE_WITH_WINDOW[ind]
                    }
                    # check if there's overlap between line x and line x+1
                    if (next_start < orig_end){
                        # in this case we need to remove some copies from first match 
                        copies_overlap = abs(next_start - orig_end)/mot_len
                        # if it is not the first iteratioon
                        if (length(multi_results) >0){
                            # combine number of copies
                            copies_info = c(as.character(unif_motif_orig), paste0(unif_motif_next, '(', copies_next, ')'))
                        # if it is the first iteration
                        } else {
                            # combine number of copies
                            copies_info = c(paste0(unif_motif_orig, '(', copies_orig - copies_overlap, ')'), paste0(unif_motif_next, '(', copies_next, ')'))
                        }
                        # save data: new copies, new start and end, new sequence information
                        new_copies = copies_orig - copies_overlap + copies_next
                        new_start = orig_start; new_end = next_end; motif_length = 'mixed'
                        sequence = seq_orig
                        new_padd_before = substr(sequence, 1, new_start-1)
                        new_sequence = substr(sequence, new_start, new_end)
                        new_padd_after = substr(sequence, new_end+1, nchar(sequence))
                    # if the TRF matches do not overlap
                    } else {
                        # save new copy information
                        if (length(multi_results) >0){
                            copies_info = c(as.character(unif_motif_orig), paste0(unif_motif_next, '(', copies_next, ')'))
                        } else {
                            copies_info = c(paste0(unif_motif_orig, '(', copies_orig, ')'), paste0(unif_motif_next, '(', copies_next, ')'))
                        }
                        # finally save data
                        new_copies = copies_orig + copies_next
                        new_start = next_start; new_end = next_end; motif_length = 'mixed'
                        sequence = seq_orig 
                        new_padd_before = substr(sequence, 1, orig_start-1)
                        new_sequence = substr(sequence, orig_start, next_end)
                        new_padd_after = substr(sequence, next_end+1, nchar(sequence))
                    }
                    # finally combine data together in a df and add to a growing list
                    combined_motifs = paste(copies_info, collapse = '/')
                    tmp_df = data.frame(combined_mot = combined_motifs, copies = new_copies, start = new_start, end = new_end, length_motif = motif_length, seq = new_sequence, padd_before = new_padd_before, padd_after = new_padd_after)
                    multi_results[[(length(multi_results) + 1)]] = tmp_df
                    ind = ind + 1
                }
                # combine results at the end of the loop. generate a new line with all combined information
                multi_results = rbindlist(multi_results); new_res = tmp[1, ]
                new_res$START_TRF = paste(multi_results$start, collapse = ','); new_res$END_TRF = paste(multi_results$end, collapse = ',')
                new_res$LENGTH_MOTIF_TRF = motif_length;
                new_res$COPIES_TRF = max(multi_results$copies); new_res$PC_MATCH_TRF = NA; new_res$PC_INDEL_TRF = NA; new_res$MOTIF_TRF = multi_results$combined_mot[nrow(multi_results)]
                new_res$MATCH_TYPE = 'mixed'; new_res$UNIFORM_MOTIF = multi_results$combined_mot[nrow(multi_results)]; new_res$PERC_SEQ_COVERED_TR = new_res$CUM_PERC_MATCH
                new_res$PADDING_BEFORE = multi_results$padd_before[nrow(multi_results)]; new_res$SEQUENCE_TRF = multi_results$seq[nrow(multi_results)]; new_res$PADDING_AFTER = multi_results$padd_after[nrow(multi_results)]
                combined_res[[(length(combined_res) + 1)]] = new_res
            # below is when there are only 2 duplicates (majority of cases)
            } else {
                # get the start and end of both matches
                orig_start = tmp$START_TRF[1]; orig_end = tmp$END_TRF[1]; next_start = tmp$START_TRF[2]; next_end = tmp$END_TRF[2]
                # check if they overlap
                if (next_start < orig_end){
                    # in this case we need to remove some copies from first match, then combine the matches
                    copies_overlap = abs(next_start - orig_end)/tmp$LENGTH_MOTIF_TRF[1]
                    copies_info = c(paste0(tmp$UNIFORM_MOTIF[1], '(', tmp$COPIES_TRF[1] - copies_overlap, ')'), paste0(tmp$UNIFORM_MOTIF[2], '(', tmp$COPIES_TRF[2], ')'))
                    new_copies = tmp$COPIES_TRF[1] - copies_overlap + tmp$COPIES_TRF[2]
                    new_start = orig_start; new_end = next_end; motif_length = 'mixed'
                    sequence = tmp$SEQUENCE_WITH_WINDOW[1]; 
                    new_padd_before = substr(sequence, 1, new_start-1)
                    new_sequence = substr(sequence, new_start, new_end)
                    new_padd_after = substr(sequence, new_end+1, nchar(sequence))
                # if they don't overlap
                } else {
                    copies_info = c(paste0(tmp$UNIFORM_MOTIF[1], '(', tmp$COPIES_TRF[1], ')'), paste0(tmp$UNIFORM_MOTIF[2], '(', tmp$COPIES_TRF[2], ')'))
                    new_copies = tmp$COPIES_TRF[1] + tmp$COPIES_TRF[2]
                    new_start = paste(orig_start, next_start, sep = ","); new_end = paste(orig_end, next_end, sep = ","); motif_length = 'mixed'
                    sequence = tmp$SEQUENCE_WITH_WINDOW[1]; 
                    new_padd_before = substr(sequence, 1, tmp$START_TRF[1]-1)
                    new_sequence = substr(sequence, tmp$START_TRF[1], tmp$END_TRF[2])
                    new_padd_after = substr(sequence, tmp$END_TRF[2]+1, nchar(sequence))
                }
                # combine matches at the end
                combined_motifs = paste(copies_info, collapse = '/')
                new_res = tmp[1, ]; new_res$START_TRF = new_start; new_res$END_TRF = new_end; new_res$LENGTH_MOTIF_TRF = motif_length;
                new_res$COPIES_TRF = new_copies; new_res$PC_MATCH_TRF = NA; new_res$PC_INDEL_TRF = NA; new_res$MOTIF_TRF = combined_motifs
                new_res$MATCH_TYPE = 'mixed'; new_res$UNIFORM_MOTIF = combined_motifs; new_res$PERC_SEQ_COVERED_TR = new_res$CUM_PERC_MATCH
                new_res$PADDING_BEFORE = new_padd_before; new_res$SEQUENCE_TRF = new_sequence; new_res$PADDING_AFTER = new_padd_after
                combined_res[[(length(combined_res) + 1)]] = new_res
            }
        }
        # at the very end, combine results in a df and later with other exact matches
        combined_res = rbindlist(combined_res)
        combined_res = rbind(combined_res, exact_motif)
        return(combined_res)
    }

    # Function to call haplotypes when dealing with assembly-based TRF
    assemblyBased = function(data, sample_name, region){
        # assign variables to fill
        data$type = NA
        data$sample = sample_name
        # if we're dealing with assemblies, in theory it is simpler the haplotyping
        # we need to check whether we have enough data (>0 contigs)
        if (nrow(data) >0){
            # flag cases when there are multiple data (>2 contigs)
            if (nrow(data) >2){
                # in this case run the normal kmean function
                res = KMeansBased_haplotyping(reads = data$COPIES_TRF, thr = 2, min_support = 1, thr_mad_orig = 0.015, type = 'single_sample')
                # then polish with the polisher without phasing
                res_polished = polishHaplo_noPhasing(res, 0.015)
                # extract polished reads
                reads_h1 = unlist(strsplit(as.character(res_polished$reads_h1), ',')); reads_h2 = unlist(strsplit(as.character(res_polished$reads_h2), ','))
                # make the final df with reads and haplotypes (depending on how many alleles were found)
                if (reads_h2 == 'homozygous' && !is.na(reads_h2)){
                    data$type[which(data$COPIES_TRF %in% reads_h1)] = 'Assigned'; data$HAPLOTYPE[which(data$COPIES_TRF %in% reads_h1)] = 1
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

# Manage arguments
    parser <- ArgumentParser()
    # add arguments: --trf is the trf output
    parser$add_argument("--trf", default = 'None', help = "Output file(s) of TRF analysis. Multiple files should be comman-separated.")
    # add arguments: --phase
    parser$add_argument("--phase", default = 'None', help = "Output file(s) of PHASING analysis. Multiple files should be comman-separated.")
    # add arguments: --trf_type
    parser$add_argument("--trf_type", default = 'None', help = "Whether TRF results are from single_reads or assembly")
    # add arguments: --out
    parser$add_argument("--out", default = 'None', help = "Output directory where output will be placed.")
    # read arguments
    args <- parser$parse_args()
    inp_trf = args$trf; inp_trf = unlist(strsplit(inp_trf, ','))
    inp_pha = args$phase; inp_pha = unlist(strsplit(inp_pha, ','))
    inp_typ = args$trf_type
    out_dir = args$out
    # check inputs and print summary
    if (inp_trf[1] == 'None' | out_dir == 'None' | inp_typ == 'None'){ RUN = FALSE } else { RUN = TRUE }
    if (RUN == FALSE){
        stop("Input error: No TRF input, input type or output directory detected!")
    } else {
        cat('** Downstream analysis of TRF\n\n')
        cat(paste0('**** input TRF selected -> ', inp_trf, '\n'))
        cat(paste0('**** input TRF type -> ', inp_typ, '\n'))
        cat(paste0('**** input PHASING selected -> ', inp_pha, '\n'))
        cat(paste0('**** output directory -> ', out_dir, '\n\n'))
    }

# Main script
    # Read and combine all input TRF files
    cat('** Analysis started\n')
    cat('** Reading TRF input\n')
    all_trf = list()
    for (i in 1:length(inp_trf)){
        # load trf datasets
        data = fread(inp_trf[i], h=T, stringsAsFactors = F)
        data$BATCH_TRF = i
        all_trf[[(length(all_trf) + 1)]] = data
    }
    all_trf = rbindlist(all_trf)
    
    # Read and combine all input PHASING files
    if (inp_pha != 'None'){
        cat('** Reading PHASING input\n')
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

    # Merge trf and phases -- exclude reference before
    reference = all_trf[which(all_trf$SAMPLE_NAME == 'reference'),]
    all_trf = all_trf[which(all_trf$SAMPLE_NAME != 'reference'),]
    trf_pha = merge(all_trf, all_pha, by.x = 'READ_NAME', by.y = 'READ_ID', all.x = T)
    trf_pha = rbind.fill(trf_pha, reference)
  
    # Let's adjust the motifs
    cat('** Generating consensus motifs\n')
    all_motifs = data.frame(motif = unique(trf_pha$MOTIF_TRF), stringsAsFactors = F); all_motifs$main_motif = NA
    all_motifs = all_motifs[!is.na(all_motifs$motif),]
    trf_pha$UNIFORM_MOTIF = NA
    # manage the reference
    trf_pha$READ_NAME[is.na(trf_pha$READ_NAME)] <- 'reference'
    for (i in 1:nrow(all_motifs)){
        if (is.na(all_motifs$main_motif[i])){
            all_perm = permutMotif(as.character(all_motifs$motif[i]))
            all_motifs$main_motif[i] = all_motifs$motif[i]
            all_motifs$main_motif[which(all_motifs$motif %in% all_perm)] = all_motifs$motif[i]
            trf_pha$UNIFORM_MOTIF[which(trf_pha$MOTIF_TRF %in% all_perm)] = all_motifs$motif[i]
        }
    }
    # to make sure there are no conflicts, change the read name
    trf_pha$READ_NAME = paste(trf_pha$READ_NAME, trf_pha$SAMPLE_NAME, sep = "__")

    # Manage when 2 reads report results
    cat('** Cleaning TRF matches\n')
    res_cleaning = cleanMotifsTRF(trf_pha)
    tokeep = res_cleaning[[1]]; toexclude = res_cleaning[[2]]
      
    # now genotype
    cat('** Genotyping TRs\n')
    all_samples = unique(tokeep$SAMPLE_NAME); all_regions = unique(tokeep$REGION); all_res = list()
    for (s in all_samples){
        print(paste0('** processing sample --> ', s))
        for (r in all_regions){
            # get data of the sample and the region of interest
            tmp_data = tokeep[which(tokeep$SAMPLE_NAME == s & tokeep$REGION == r),]
            if (nrow(tmp_data) >0){
                if (inp_typ == 'single_reads'){
                    if ('true_variable_motif' %in% tmp_data$EXITUS){
                        # here I need to merge the possible alternative motifs
                        tmp_data = mergeMotifs(tmp_data)
                    }
                    reads_df = tmp_data[, c('COPIES_TRF', 'HAPLOTYPE', 'UNIFORM_MOTIF', 'REGION', 'PASSES', 'READ_QUALITY', 'LENGTH_SEQUENCE')]
                    phased_data = PhasingBased_haplotyping(reads_df, sample_name = s, thr_mad = 0.015)
                    polished_data = polishHaplo_afterPhasing(phased_data, 0.015)
                    all_res[[(length(all_res) + 1)]] = polished_data
                } else {
                    if ('true_variable_motif' %in% tmp_data$EXITUS){
                        # merge possible alternative motifs
                        tmp_data = mergeMotifs(tmp_data)
                    }
                    reads_df = tmp_data[, c('COPIES_TRF', 'HAPLOTYPE', 'UNIFORM_MOTIF', 'REGION', 'PASSES', 'READ_QUALITY', 'LENGTH_SEQUENCE')]
                    phased_data = assemblyBased(reads_df, sample_name = s, region = r)
                    polished_data = polishHaplo_afterPhasing(phased_data, 0.015)
                    all_res[[(length(all_res) + 1)]] = polished_data
                }
            }
        }
    }
    # Merge all results together
    all_res = data.table::rbindlist(all_res, fill=T)
  
    # Call haplotypes and calculate coefficient of variation
    cat('** Merging haplotypes\n')
    all_haplo = callHaplo(data = all_res)
  
    # Saving data (rdata object + plain text file)
    trf_matches_kept = tokeep; trf_matches_excluded = toexclude
    trf_matches_genotypes_per_read = all_res
    trf_matches_genotypes_per_sample = all_haplo
    all_results = list(trf_matches_kept, trf_matches_excluded, trf_matches_genotypes_per_read, trf_matches_genotypes_per_sample)
    # rdata object with all results
    save(all_results, file = paste0(out_dir, '/haplotyping.RData'))
    # plain text for genotypes per read and per sample
    write.table(trf_matches_genotypes_per_read, paste0(out_dir, '/Haplotypes_per_read.txt'), quote=F, row.names=F, sep = '\t')
    write.table(trf_matches_genotypes_per_sample, paste0(out_dir, '/Haplotypes_per_sample.txt'), quote=F, row.names=F, sep = '\t')