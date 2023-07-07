#!/usr/bin/env python3

# libraries
print('**** Loading libraries')
import pandas as pd
import sys
from Bio.Seq import reverse_complement
import multiprocessing
from functools import partial
import numpy as np
import math
import os
import scipy.stats as stats
import statistics
from sklearn.cluster import KMeans
import numpy as np
import shutil
import warnings
import gzip

# FUNCTIONS
# function to make permutations
def permutMotif(motif):
    # standard motif
    motif_cons = [motif[i::] + motif[0:i] for i in range(len(motif))]
    # reverse complement
    rc_motif = reverse_complement(motif)
    rc_motif_cons = [rc_motif[i::] + rc_motif[0:i] for i in range(len(rc_motif))]
    # combine
    comb_motifs = motif_cons + rc_motif_cons
    # take the first in ascending order
    sel_motif = sorted(comb_motifs)[0]
    return sel_motif

# function to guide haplotyping
def haplotyping(r, s, thr_mad, data_nodup, type, dup_df, reference_motif_dic, intervals, min_support):
    if r in intervals:
        print('****** done %s%% of the regions' %(intervals.index(r)*5+5), end = '\r')
    # data of interest
    sbs = data_nodup[data_nodup['REGION'] == r].copy()
    # exclude nas
    sbs = sbs.dropna(subset=['LEN_SEQUENCE_FOR_TRF'])
    # check if there are rows
    if sbs.shape[0] >0:
        if type == 'assembly':
            # identify haplotypes
            phased_sbs = assemblyBased_size(sbs, r)
            # polish haplotypes
            pol_sbs = polishHaplo_asm(phased_sbs, r)
            # add duplicates
            all_sbs = addDups(pol_sbs, dup_df, s, r)
            # finally look at the motif
            final_sbs_h1 = sampleMotifs(r, all_sbs, reference_motif_dic, 1)
            final_sbs_h2 = sampleMotifs(r, all_sbs, reference_motif_dic, 2)
            final_sbs = pd.concat([final_sbs_h1, final_sbs_h2], axis=0)
            # prepare for file writing
            tmp_vcf, tmp_seq = prepareOutputs(final_sbs, reference_motif_dic, r, type, 'None')
            return tmp_vcf, tmp_seq
        elif type == 'reads':
            # check minimum support: minimum support is for alleles --> 2*min_support is the total minimum coverage required for autosomal regions. For sex-regions, we will use min_support directly
            # find chromosome to adapt coverage
            chrom = r.split(':')[0]
            minimum_coverage = min_support if chrom in ['chrY', 'Y'] else min_support*2
            # check if there is minimum support
            if sbs.shape[0] >= minimum_coverage:
                # identify haplotypes
                warnings.filterwarnings("ignore", category=RuntimeWarning)
                pol_sbs = readBased_size(sbs, r, chrom)
                warnings.resetwarnings()
                # add duplicates
                all_sbs = addDups(pol_sbs, dup_df, s, r, type)
                # finally look at the motif
                final_sbs_h1, depth_h1 = sampleMotifs(r, all_sbs, reference_motif_dic, 0, type)
                final_sbs_h2, depth_h2 = sampleMotifs(r, all_sbs, reference_motif_dic, 1, type)
                final_sbs = pd.concat([final_sbs_h1, final_sbs_h2], axis=0)
                # prepare for file writing
                tmp_vcf, tmp_seq = prepareOutputs(final_sbs, reference_motif_dic, r, type, [depth_h1, depth_h2])
                return tmp_vcf, tmp_seq
            else:
                tmp_vcf = [chrom, r.split(':')[1].split('-')[0], r, reference_motif_dic[r][1], '.', '.', 'LOW_COVERAGE', '%s;%s' %(reference_motif_dic[r][0], reference_motif_dic[r][2]), 'QC;GT;MOTIF;CN;CN_REF;DP', 'QC_ISSUE;NA|NA;NA|NA;NA|NA;NA|NA;%s|0' %(sbs.shape[0])]
                tmp_seq = []
                return tmp_vcf, tmp_seq
    else:
        chrom = r.split(':')[0]
        tmp_vcf = [chrom, r.split(':')[1].split('-')[0], r, reference_motif_dic[r][1], '.', '.', 'LOW_COVERAGE', '%s;%s' %(reference_motif_dic[r][0], reference_motif_dic[r][2]), 'QC;GT;MOTIF;CN;CN_REF;DP', 'QC_ISSUE;NA|NA;NA|NA;NA|NA;NA|NA;%s|0' %(sbs.shape[0])]
        tmp_seq = []
        return tmp_vcf, tmp_seq

# function to call haplotypes for reads
def readBased_size(sbs, r, chrom):
    # get phasing information
    haplotag_list = [x for x in list(sbs['HAPLOTAG']) if not isinstance(x, float) or not math.isnan(x)]
    if len(haplotag_list) == 0:
        haplo_center, haplo_id, deleted = kmeans_haplotyping(sbs, min_support, thr_mad, chrom, r)
        # check if there were excluded reads
        if len(deleted) == 0:
            #sbs['HAPLO_CONFINT'] = haplo_confint
            sbs['HAPLOTAG'] = haplo_id
            sbs['POLISHED_HAPLO'] = haplo_center
        else:
            deleted_rows = sbs[sbs['LEN_SEQUENCE_FOR_TRF'].isin(deleted)].copy()
            deleted_rows['HAPLOTAG'] = 'NA'
            deleted_rows['POLISHED_HAPLO'] = 'NA'
            kept_rows = sbs[~sbs['LEN_SEQUENCE_FOR_TRF'].isin(deleted)].copy()
            kept_rows['HAPLOTAG'] = haplo_id
            kept_rows['POLISHED_HAPLO'] = haplo_center
            sbs = pd.concat([kept_rows, deleted_rows], axis=0)
        return sbs
    else:
        print('!! haplotags found at %s' %(r))

# function to check deviations within each haplotype
def checkDeviation(read_lengths):
    # calculate median and thr_mad % value of the median
    median_value = statistics.median(read_lengths)
    median_value_boundary = median_value*thr_mad
    # and distances from median
    median_distances = [abs(x - median_value) for x in read_lengths]
    return median_value_boundary, median_distances

# function to find median and confint for homozygous calls
def findHaplo_homozygous(read_lengths, deleted, thr_mad):
    # find haplotypes
    centers_kmeans = [statistics.median(read_lengths) for x in range(len(read_lengths))]
    # if there are deleted reads, then try to recover them
    if len(deleted) >0:
        median_value = list(set(centers_kmeans))[0]
        keep = [x for x in deleted if abs(x - median_value) < median_value*thr_mad]
        # if there are reads to keep, add them and recalculate
        if len(keep) >0:
            read_lengths = read_lengths + keep
            deleted = [x for x in deleted if x not in keep]
            centers_kmeans = [statistics.median(read_lengths) for x in range(len(read_lengths))]
    # find confidence intervals
    centers_confint = ['_'.join([str(element) for element in list(stats.t.interval(0.95, len(read_lengths)-1, statistics.mean(read_lengths), statistics.stdev(read_lengths)/len(read_lengths)**0.5))]) for x in range(len(read_lengths))]
    # find haplotype assignment
    haplo_list = [1 for x in range(len(read_lengths))]
    return centers_kmeans, centers_confint, haplo_list, deleted

# function to find median and confint for heterozygous calls
def findHaplo_hetero(read_lengths, haplo_list, centers_kmeans):
    # find haplotypes
    centers = [centers_kmeans[x] for x in haplo_list]
    # confidence intervals
    reads_h1 = [read_lengths[i] for i, x in enumerate(haplo_list) if x == 0]
    reads_h2 = [read_lengths[i] for i, x in enumerate(haplo_list) if x == 1]
    #confint_h1 =  list(stats.t.interval(0.95, len(reads_h1)-1, statistics.mean(reads_h1), statistics.stdev(reads_h1)/len(reads_h1)**0.5))
    #confint_h2 =  list(stats.t.interval(0.95, len(reads_h2)-1, statistics.mean(reads_h2), statistics.stdev(reads_h2)/len(reads_h2)**0.5))
    return centers, haplo_list

# function to do kmeans and give centers and haplotype lists
def kmeans(read_lengths, ploidy):
    # do kmeans with the ploidy as the number of clusters
    my_array = np.array(read_lengths).reshape(-1, 1)
    # perform k-means clustering with 2 clusters
    kmeans = KMeans(n_clusters=ploidy).fit(my_array)
    centers_kmeans = [center for sublist in kmeans.cluster_centers_.tolist() for center in sublist]
    haplo_list = kmeans.labels_.tolist()
    return centers_kmeans, haplo_list

# function to check whether there's support for haplotypes
def checkSupport(read_lengths, haplo_list, min_support):
    haplo_check = []
    for haplo in list(set(haplo_list)):
        indexes = [i for i in range(len(haplo_list)) if haplo_list[i] == haplo]
        if len(indexes) >= min_support:
            haplo_check.append('ok')
        else:
            haplo_check.append('delete')
    # check if both haplotypes were ok
    if 'delete' in haplo_check:
        index_delete = haplo_check.index('delete')
        # get index of the relative reads
        index_delete_reads = [i for i in range(len(haplo_list)) if haplo_list[i] == index_delete]
        # update reads
        deleted = [x for i, x in enumerate(read_lengths) if i in index_delete_reads]
        read_lengths = [x for i, x in enumerate(read_lengths) if i not in index_delete_reads]
    else:
        deleted = []
    loop = True if len(deleted) >0 else False                
    return read_lengths, deleted, loop

# function to fit kmeans for finding haplotypes
def kmeans_haplotyping(sbs, min_support, thr_mad, chrom, r):
    # extract read lengths
    read_lengths = [x for x in list(sbs['LEN_SEQUENCE_FOR_TRF']) if not isinstance(x, float) or not math.isnan(x)]
    # define ploidy: for autosomal chromosomes and X it's 2, otherwise 1
    ploidy = 1 if chrom in ['chrY', 'Y'] else 2
    # define variables to control the loop
    loop = True
    deleted_all = []
    # main loop
    while loop == True:
        # check if we have minimum support in terms of coverage
        if len(read_lengths) >= min_support:
            # check if we have a homozygous call (reads are similar in size)
            median_value_boundary, median_distances = checkDeviation(read_lengths)
            # decide whether it's homozygous or heterozygous
            call = 'homo' if (max(median_distances) <= median_value_boundary or ploidy == 1) else 'hetero'
            if call == 'homo':
                # homozygous, we're done
                centers_kmeans, centers_confint, haplo_list, deleted = findHaplo_homozygous(read_lengths, deleted_all, thr_mad)
                deleted_all = deleted
                loop = False
            else:
                # heterozygous, do kmeans
                centers_kmeans, haplo_list = kmeans(read_lengths, ploidy)
            # now we need to check the support for the haplotypes (if not homozygous)
            if loop == True:
                read_lengths, deleted, loop = checkSupport(read_lengths, haplo_list, min_support)
                deleted_all.extend(deleted)
                # if no loop is required, find haplotype stats
                if loop == False:
                    centers_kmeans, haplo_list = findHaplo_hetero(read_lengths, haplo_list, centers_kmeans)
                    centers_confint = []
        else:
            centers_kmeans = []; centers_confint = []; haplo_list = []; deleted_all = []
            loop = False
    return centers_kmeans, haplo_list, deleted_all

# function to make data for vcf writing
def prepareOutputs(final_sbs, reference_motif_dic, r, type, depths):
    # prepare data for VCF
    chrom, start, end = [r.split(':')[0]] + r.split(':')[-1].split('-')
    ref_motif, ref_len, ref_copies = reference_motif_dic[r]
    info_field = '%s;%s' %(ref_motif, ref_copies)
    format_field = 'QC;GT;MOTIF;CN;CN_REF;DP'
    sam_gt = '|'.join(list(final_sbs['POLISHED_HAPLO'].map(str))) if final_sbs.shape[0] >1 else list(final_sbs["POLISHED_HAPLO"].map(str).apply(lambda x: x + "|" + x))[0]
    sam_mot = '|'.join(list(final_sbs['CONSENSUS_MOTIF'].map(str))) if final_sbs.shape[0] >1 else list(final_sbs["CONSENSUS_MOTIF"].map(str).apply(lambda x: x + "|" + x))[0]
    sam_cop = '|'.join(list(final_sbs['CONSENSUS_MOTIF_COPIES'].map(str))) if final_sbs.shape[0] >1 else list(final_sbs["CONSENSUS_MOTIF_COPIES"].map(str).apply(lambda x: x + "|" + x))[0]
    sam_cop_ref = '|'.join(list(final_sbs['REFERENCE_MOTIF_COPIES'].map(str))) if final_sbs.shape[0] >1 else list(final_sbs["REFERENCE_MOTIF_COPIES"].map(str).apply(lambda x: x + "|" + x))[0]
    if type == 'reads':
        sam_depth = '|'.join([str(x) for x in depths])
    else:
        sam_depth = '|'.join(list(final_sbs["READ_NAME"].str.split("_").str[1])) if final_sbs.shape[0] >1 else list(final_sbs["READ_NAME"].str.split("_").str[1].apply(lambda x: x + "|" + x))[0]
    sample_field = '%s;%s;%s;%s;%s;%s' %('PASS_ASM', sam_gt, sam_mot, sam_cop, sam_cop_ref, sam_depth)
    res_vcf = [chrom, start, r, ref_len, '.', '.', 'PASS', info_field, format_field, sample_field]
    # then prepare the sequence output
    if type == 'reads':
        final_sbs['DEPTH'] = [x for x in depths if x != 0]
        final_sbs['type'] = type
    else:
        final_sbs['DEPTH'] = final_sbs["READ_NAME"].str.split("_").str[1]
    final_sbs['MOTIF_REF'] = ref_motif
    keep_cols = ['READ_NAME', 'HAPLOTAG', 'REGION', 'PASSES', 'READ_QUALITY', 'LEN_SEQUENCE_FOR_TRF', 'START_TRF', 'END_TRF', 'type', 'SAMPLE_NAME', 'POLISHED_HAPLO', 'DEPTH', 'CONSENSUS_MOTIF', 'CONSENSUS_MOTIF_COPIES', 'MOTIF_REF', 'REFERENCE_MOTIF_COPIES', 'SEQUENCE_WITH_PADDINGS', 'SEQUENCE_FOR_TRF']
    # subset dataframe to keep only selected columns
    df_subset = final_sbs.loc[:, keep_cols]
    res_seq = df_subset.values.tolist()
    return res_vcf, res_seq

# function to call haplotypes for assembly
def assemblyBased_size(sbs, r):
    n_contigs = sbs.shape[0]
    # check number of contigs
    if n_contigs == 1:
        # homozygous call
        sbs['HAPLOTAG'] = 1; sbs['type'] = 'assembly'
    elif n_contigs == 2:
        # heterozygous call
        sbs['HAPLOTAG'] = [1, 2]; sbs['type'] = 'assembly'
    elif sbs.shape[0] >2:
        print('\n!! More than 2 haps for --> %s' %(r))
    return sbs

# function to polish haplotypes
def polishHaplo_asm(phased_sbs, r):
    # find haplotypes
    n_haplo = len(list(phased_sbs['HAPLOTAG'].dropna().unique()))
    n_contigs = phased_sbs.shape[0]
    if n_haplo == 1:
        pol_sbs = phased_sbs
        pol_sbs['POLISHED_HAPLO'] = pol_sbs['LEN_SEQUENCE_FOR_TRF']
    elif n_haplo >1 and n_contigs == 2:
        pol_sbs = phased_sbs
        pol_sbs['POLISHED_HAPLO'] = pol_sbs['LEN_SEQUENCE_FOR_TRF']
    elif n_haplo >1 and n_contigs >2:
        print('\n!! More than 1 haplo and more than 2 contigs for %s' %(r))
    return pol_sbs

# function to add duplicates back
def addDups(pol_sbs, dup_df, s, r, type):
    pol_sbs['type'] = type
    # subset of the duplicates of that sample and that region
    sbs_dups = dup_df[(dup_df['SAMPLE_NAME'] == s) & (dup_df['REGION'] == r)].copy()
    if sbs_dups.shape[0] >0:
        n_haplo = len([x for x in list(pol_sbs['HAPLOTAG'].dropna().unique()) if x != 'NA'])
        if n_haplo == 1:
            sbs_dups['HAPLOTAG'] = [x for x in list(pol_sbs['HAPLOTAG'].dropna().unique()) if x != 'NA'][0]
            sbs_dups['type'] = type
            sbs_dups['POLISHED_HAPLO'] = [x for x in list(pol_sbs['POLISHED_HAPLO'].dropna().unique()) if x != 'NA'][0]
            return pd.concat([pol_sbs, sbs_dups], axis=0)
        else:
            h1_size = pol_sbs.loc[pol_sbs['HAPLOTAG'] == 0, 'POLISHED_HAPLO'].unique()
            h2_size = pol_sbs.loc[pol_sbs['HAPLOTAG'] == 1, 'POLISHED_HAPLO'].unique()
            target = list(sbs_dups['LEN_SEQUENCE_FOR_TRF'])
            haplo = []; size = []
            for x in target:
                tmp_haplo, tmp_size = assignHaplotag_asm(h1_size, h2_size, x)
                haplo.append(tmp_haplo); size.append(tmp_size)
            sbs_dups['type'] = type
            sbs_dups['POLISHED_HAPLO'] = [float(x) for x in size]
            sbs_dups['HAPLOTAG'] = [float(x) for x in haplo]
            combined = pd.concat([pol_sbs, sbs_dups], axis=0)
            combined = combined.drop_duplicates()
            return combined
    else:
        return pol_sbs

# function to assign haplotag based on closest size
def assignHaplotag_asm(h1_size, h2_size, target):
    dist_h1 = abs(h1_size - target)
    dist_h2 = abs(h2_size - target)
    if dist_h1 < dist_h2:
        return 1, h1_size
    else:
        return 2, h2_size

# function to look at reference motifs
def referenceMotifs(r, ref, intervals):
    if r in intervals:
        print('****** done %s%% of the regions' %(intervals.index(r)*5+5), end = '\r')
    # subset of reference data
    sbs = ref[ref['REGION'] == r].copy()
    # if there's only 1 motif, we are done
    if sbs.shape[0] == 1:
        sbs['CONSENSUS_MOTIF'] = sbs['UNIFORM_MOTIF']
        sbs['CONSENSUS_MOTIF_COPIES'] = sbs['COPIES_TRF']
    elif sbs.shape[0] >1:
        # first align motifs
        sbs = sbs[sbs['HAPLOTAG'] == 1].copy()
        sbs = motif_generalization(sbs, r)
    # in the end, take only what we need to bring along
    tmp_dic = {row["REGION"]: [row["CONSENSUS_MOTIF"], row["POLISHED_HAPLO"], row['CONSENSUS_MOTIF_COPIES']] for _, row in sbs.iterrows()}
    return tmp_dic

# function to generate consensus motif using majority rule
def motif_generalization(haplo_data, r):
    # calculate fraction of sequence covered
    haplo_data['COVERAGE_TR'] = (haplo_data['END_TRF'] - haplo_data['START_TRF'] + 1) / haplo_data['LEN_SEQUENCE_FOR_TRF']
    haplo_data = haplo_data.sort_values('COVERAGE_TR', ascending=False)
    # if >95% of the sequence is covered, stop
    high_coverage = haplo_data[haplo_data['COVERAGE_TR'] >0.95].copy()
    if high_coverage.shape[0] == 1:
        best_motif = list(high_coverage['UNIFORM_MOTIF'])[0]
        copies = list(high_coverage['COPIES_TRF'])[0]
        start = list(high_coverage['START_TRF'])[0]
        end = list(high_coverage['END_TRF'])[0]
    elif high_coverage.shape[0] >1:
        high_coverage['MOTIF_SCORE'] = high_coverage['COVERAGE_TR'] * high_coverage['TRF_SCORE']
        high_coverage = high_coverage.sort_values('MOTIF_SCORE', ascending=False)
        best_motif = high_coverage['UNIFORM_MOTIF'].iloc[0]
        copies = high_coverage['COPIES_TRF'].iloc[0]
        start = high_coverage['START_TRF'].iloc[0]
        end = high_coverage['END_TRF'].iloc[0]
    elif haplo_data['motif'].isna().all():
        best_motif, copies, indexes, start, end = 'NA', 'NA', 'NA', 'NA', 'NA'
    else:
        best_motif, copies, indexes = combineMotifs(haplo_data, r)
        start, end = indexes[0], indexes[-1]
    # then combine with haplotype data
    haplo_data = haplo_data.iloc[[0]]
    haplo_data['CONSENSUS_MOTIF'] = best_motif
    haplo_data['CONSENSUS_MOTIF_COPIES'] = copies
    haplo_data['START_TRF'] = start
    haplo_data['END_TRF'] = end
    return haplo_data

# function to combine motifs
def combineMotifs(haplo_data, r):
    # take the most representative motif
    haplo_data['MOTIF_SCORE'] = haplo_data['COVERAGE_TR'] * haplo_data['TRF_SCORE']
    haplo_data = haplo_data.sort_values('MOTIF_SCORE', ascending=False)
    best_motif = haplo_data['UNIFORM_MOTIF'].iloc[0]
    best_motif_copies = haplo_data['COPIES_TRF'].iloc[0]
    best_motif_index = range(int(haplo_data['START_TRF'].iloc[0]), int(haplo_data['END_TRF'].iloc[0]))
    # loop on the other motifs
    for i in range(1, haplo_data.shape[0]):
        try:
            # take corresponding motif and range
            tmp_motif = haplo_data['UNIFORM_MOTIF'].iloc[i]
            tmp_motif_copies = haplo_data['COPIES_TRF'].iloc[i]
            tmp_motif_index = range(int(haplo_data['START_TRF'].iloc[i]), int(haplo_data['END_TRF'].iloc[i]))
            # find the length of the intersection
            intersection_length = max(0, min(best_motif_index.stop, tmp_motif_index.stop) - max(best_motif_index.start, tmp_motif_index.start))
            # calculate the percentage of overlap
            fraction_overlap = intersection_length / len(tmp_motif_index)
            # if the fraction is below 90%, then we should combine the motifs
            if fraction_overlap < 0.90:
                # update best_motif and index: this is a combination of the two motifs
                if best_motif == tmp_motif:
                    # here the simple case: motif is the same
                    best_motif_index = range(min(best_motif_index.start, tmp_motif_index.start), max(best_motif_index.stop, tmp_motif_index.stop))
                    best_motif_copies = len(list(best_motif_index))/len(best_motif)
                else:
                    # here the more tricky case: motifs are different
                    not_in_best = [x for x in tmp_motif_index if x not in best_motif_index]
                    # find list of consecutive numbers
                    consecutive_seq = is_consecutive(not_in_best)
                    # iterate over each subsequence: if a motif can be fit, add it, otherwise skip
                    for subs in consecutive_seq:
                        if len(subs) > len(tmp_motif):
                            sub_copies = len(subs)/len(tmp_motif)
                            best_motif = '%s+%s' %(best_motif, tmp_motif)
                            best_motif_copies = '%s+%s' %(best_motif_copies, sub_copies)
                            best_motif_index = range(min(best_motif_index.start, subs[0]), max(best_motif_index.stop, subs[-1]))
            else:
                best_motif_index = best_motif_index
        except:
            pass
    return best_motif, best_motif_copies, best_motif_index

# function to find sequences of consecutive integers
def is_consecutive(numbers):
    current_seq = []
    sequences = []
    # iterate over numbers and find consecutive sequences
    for i, n in enumerate(numbers):
        if i == 0 or n != numbers[i-1]+1:
            # start new sequence
            current_seq = [n]
            sequences.append(current_seq)
        else:
            # continue current sequence
            current_seq.append(n)
    return sequences

# function to look at the motif of the samples
def sampleMotifs(r, all_sbs, reference_motif_dic, haplo, type):
    # subset of data
    haplo_data = all_sbs[all_sbs['HAPLOTAG'] == haplo].copy()
    haplo_depth = len(set(list(haplo_data['UNIQUE_NAME'])))
    # if there's only 1 motif, we are done
    if haplo_data.shape[0] == 1:
        haplo_data['CONSENSUS_MOTIF'] = haplo_data['UNIFORM_MOTIF']
        haplo_data['CONSENSUS_MOTIF_COPIES'] = haplo_data['COPIES_TRF']
    elif haplo_data.shape[0] >1:
        haplo_data = motif_generalization(haplo_data, r)
    # finally wrt reference motif
    try:
        ref_motif = reference_motif_dic[r][0]
        if haplo_data['CONSENSUS_MOTIF'].iloc[0] == ref_motif:
            haplo_data['REFERENCE_MOTIF_COPIES'] = haplo_data['CONSENSUS_MOTIF_COPIES']
        elif isinstance(ref_motif, (int, str)):
            haplo_data['REFERENCE_MOTIF_COPIES'] = haplo_data['POLISHED_HAPLO'] / len(ref_motif)
        else:
            haplo_data['REFERENCE_MOTIF_COPIES'] = 'nan'
    except:
        haplo_data['REFERENCE_MOTIF_COPIES'] = 'nan'
    if type == 'reads':
        return haplo_data, haplo_depth
    else:
        return haplo_data

# function to prepare intervals
def prepareIntervals(all_regions):
    intervals = []
    for percentile in np.arange(0.05, 1, 0.05).tolist():
        index = int(len(all_regions) * percentile)
        intervals.append(all_regions[index])
    return intervals

# function to write header of vcf file
def writeVCFheader(vcf_file, samples):
    # open file
    outf = open(vcf_file, 'w')
    # write header
    outf.write('##fileformat=VCFv4.2\n##INFO=<ID=REFERENCE_INFO,Number=2,Type=String,Description="Motif observed in the reference genome (GRCh38), and relative number of motif repetitions."\n##FORMAT=<ID=QC,Number=1,Type=String,Description="Quality summary of TREAT genotyping. PASS_BOTH: genotype agreed between reads-spanning and assembly. PASS_RSP: genotype from reads-spanning. PASS_ASM: genotype from assembly."\n##FORMAT=<ID=GT,Number=2,Type=String,Description="Phased size of the tandem repeats. H1_size | H2_size"\n##FORMAT=<ID=MOTIF,Number=2,Type=String,Description="Phased consensus motif found in the sample. H1_motif | H2_motif"\n##FORMAT=<ID=CN,Number=2,Type=String,Description="Phased number of copies of the motif found in the sample. H1_copies | H2_copies"\n##FORMAT=<ID=CN_REF,Number=2,Type=String,Description="Phased estimation of the reference motif as found in the sample. H1_motif_ref | H2_motif_ref"\n##FORMAT=<ID=DP,Number=1,Type=String,Description="Phased depth found in the sample. H1_depth | H2_depth"\n')
    outf.close()
    return    

# function to write outputs
def writeOutputs(df_vcf, df_seq, seq_file, vcf_file):
    # write header of vcf file
    writeVCFheader(vcf_file, all_samples)
    with open(vcf_file, mode='a') as file:
        df_vcf.to_csv(file, header=True, index=False, sep='\t')
    # then compress it
    os.system('gzip %s' %(vcf_file))
    os.system('rm %s' %(vcf_file))
    # write sequence file
    df_seq.to_csv(seq_file, header=True, index=False, sep='\t', compression = 'gzip')
    return

# 1. arguments
inpf = sys.argv[1]
outd = sys.argv[2]
n_cpu = int(sys.argv[3])
thr_mad = float(sys.argv[4])
type = sys.argv[5]
min_support = int(sys.argv[6])

# 2. read data
print('** Read data')
data = pd.read_csv(inpf, sep='\t', compression='gzip')

# 3. adjust the motif -- merge the same motifs
print('** Adjust motifs')
all_motifs = data['TRF_MOTIF'].dropna().unique()
main_motifs = [permutMotif(motif) for motif in all_motifs]
motifs_df = pd.DataFrame({'motif' : all_motifs, 'UNIFORM_MOTIF' : main_motifs})
data = pd.merge(data, motifs_df, left_on='TRF_MOTIF', right_on='motif', how='left')
data['UNIQUE_NAME'] = data.apply(lambda row: str(row['READ_NAME']) + '___' + str(row['SAMPLE_NAME']) + '___' + str(row['REGION']) + '___' + str(row['LEN_SEQUENCE_FOR_TRF']), axis = 1)

# 4. extract the reference and work on the motifs
print('** Reference motifs                                     ')
ref = data[data['SAMPLE_NAME'] == 'reference'].copy()
all_regions = list(ref['REGION'].dropna().unique())
intervals = prepareIntervals(all_regions)
ref['HAPLOTAG'] = 1; ref['POLISHED_HAPLO'] = ref['LEN_SEQUENCE_FOR_TRF']
all_regions = list(ref['REGION'].dropna().unique())
pool = multiprocessing.Pool(processes=n_cpu)
motif_fun = partial(referenceMotifs, ref = ref, intervals = intervals)
motif_res = pool.map(motif_fun, all_regions)
# combine dictionaries
reference_motif_dic = {k: v for d in motif_res for k, v in d.items()}

# 5. add unique name and split duplicated from unique reads
data = data[data['SAMPLE_NAME'] != 'reference']
data_nodup = data.drop_duplicates(subset = 'UNIQUE_NAME')
dup_df = data[data.duplicated(subset = 'UNIQUE_NAME', keep=False)]

# 6. genotyping on the sizes
print('** Genotyping                                         ')
all_samples = data_nodup['SAMPLE_NAME'].dropna().unique()
all_regions = list(data_nodup['REGION'].dropna().unique())
intervals = prepareIntervals(all_regions)
sample_res = []
for s in all_samples:
    print('**** %s                      ' %(s))
    sbs = data_nodup[(data_nodup['SAMPLE_NAME'] == s)]
    sbs_dups = dup_df[(dup_df['SAMPLE_NAME'] == s)]
    #for r in all_regions:
    #    haplo_results = haplotyping(r, s, thr_mad, sbs, type, sbs_dups, reference_motif_dic, intervals, min_support)
    #    sample_res.append(haplo_results)
    pool = multiprocessing.Pool(processes=n_cpu)
    haplo_fun = partial(haplotyping, s = s, thr_mad = thr_mad, data_nodup = sbs, type = type, dup_df = sbs_dups, reference_motif_dic = reference_motif_dic, intervals = intervals, min_support = 3)
    haplo_results = pool.map(haplo_fun, all_regions)
    sample_res.append(haplo_results)

# 7. compose the vcf and seq dataframes
df_vcf = pd.DataFrame([x[0] for x in sample_res[0]], columns=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', all_samples[0]])
df_seq = pd.DataFrame([x for y in sample_res[0] for x in y[1]], columns=['READ_NAME', 'HAPLOTAG', 'REGION', 'PASSES', 'READ_QUALITY', 'LEN_SEQUENCE_FOR_TRF', 'START_TRF', 'END_TRF', 'type', 'SAMPLE_NAME', 'POLISHED_HAPLO', 'DEPTH', 'CONSENSUS_MOTIF', 'CONSENSUS_MOTIF_COPIES', 'MOTIF_REF', 'REFERENCE_MOTIF_COPIES', 'SEQUENCE_WITH_PADDINGS', 'SEQUENCE_FOR_TRF'])
if len(sample_res) >1:
    for i in range(1, len(sample_res)):
        tmp = pd.DataFrame([[x[0][2], x[0][-1]] for x in sample_res[i]], columns=['ID', all_samples[i]])
        # merge with main df
        df_vcf = pd.merge(df_vcf, tmp, on='ID', how='outer')
        # then compose the seq dataframe
        tmp = pd.DataFrame([x for y in sample_res[i] for x in y[1]], columns=['READ_NAME', 'HAPLOTAG', 'REGION', 'PASSES', 'READ_QUALITY', 'LEN_SEQUENCE_FOR_TRF', 'START_TRF', 'END_TRF', 'type', 'SAMPLE_NAME', 'POLISHED_HAPLO', 'DEPTH', 'CONSENSUS_MOTIF', 'CONSENSUS_MOTIF_COPIES', 'MOTIF_REF', 'REFERENCE_MOTIF_COPIES', 'SEQUENCE_WITH_PADDINGS', 'SEQUENCE_FOR_TRF'])
        # add to main df
        df_seq = pd.concat([df_seq, tmp], axis=0)

# 8. write outputs: vcf file and raw sequences
print('** Producing outputs: VCF file and table with sequences                        ')
seq_file = '%s/sample.seq.txt.gz' %(outd)
vcf_file = '%s/sample.vcf' %(outd)
writeOutputs(df_vcf, df_seq, seq_file, vcf_file)

