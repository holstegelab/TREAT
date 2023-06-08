import pandas as pd
import sys
from Bio.Seq import reverse_complement
import multiprocessing
from functools import partial
import numpy as np
import math

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
def haplotyping(r, s, thr_mad, data_nodup, type, dup_df, reference_motif_dic, intervals):
    if r in intervals:
        print('****** done %s%% of the regions' %(intervals.index(r)*5+5), end = '\r')
    # data of interest
    sbs = data_nodup[data_nodup['REGION'] == r]
    # exclude nas
    sbs = sbs.dropna(subset=['LEN_SEQUENCE_FOR_TRF'])
    # check if there are rows
    if sbs.shape[0] >0:
        if type == 'asm':
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
            tmp_vcf, tmp_seq = prepareOutputs(final_sbs, reference_motif_dic, r)
            return tmp_vcf, tmp_seq
        elif type == 'reads_spanning':
            return None

# function to make data for vcf writing
def prepareOutputs(final_sbs, reference_motif_dic, r):
    # prepare data for VCF
    chrom, start, end = [r.split(':')[0]] + r.split(':')[-1].split('-')
    ref_motif, ref_len, ref_copies = reference_motif_dic[r]
    info_field = '%s;%s' %(ref_motif, ref_copies)
    format_field = 'QC;GT;MOTIF;CN;CN_REF;DP'
    sam_gt = '|'.join(list(final_sbs['POLISHED_HAPLO'].map(str))) if final_sbs.shape[0] >1 else list(final_sbs["POLISHED_HAPLO"].map(str).apply(lambda x: x + "|" + x))[0]
    sam_mot = '|'.join(list(final_sbs['CONSENSUS_MOTIF'].map(str))) if final_sbs.shape[0] >1 else list(final_sbs["CONSENSUS_MOTIF"].map(str).apply(lambda x: x + "|" + x))[0]
    sam_cop = '|'.join(list(final_sbs['CONSENSUS_MOTIF_COPIES'].map(str))) if final_sbs.shape[0] >1 else list(final_sbs["CONSENSUS_MOTIF_COPIES"].map(str).apply(lambda x: x + "|" + x))[0]
    sam_cop_ref = '|'.join(list(final_sbs['REFERENCE_MOTIF_COPIES'].map(str))) if final_sbs.shape[0] >1 else list(final_sbs["REFERENCE_MOTIF_COPIES"].map(str).apply(lambda x: x + "|" + x))[0]
    sam_depth = '|'.join(list(final_sbs["READ_NAME"].str.split("_").str[1])) if final_sbs.shape[0] >1 else list(final_sbs["READ_NAME"].str.split("_").str[1].apply(lambda x: x + "|" + x))[0]
    sample_field = '%s;%s;%s;%s;%s;%s' %('PASS_ASM', sam_gt, sam_mot, sam_cop, sam_cop_ref, sam_depth)
    res_vcf = [chrom, start, r, ref_len, '.', '.', 'PASS', info_field, format_field, sample_field]
    # then prepare the sequence output
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
def addDups(pol_sbs, dup_df, s, r):
    # subset of the duplicates of that sample and that region
    sbs_dups = dup_df[(dup_df['SAMPLE_NAME'] == s) & (dup_df['REGION'] == r)].copy()
    if sbs_dups.shape[0] >0:
        n_haplo = len(list(pol_sbs['HAPLOTAG'].dropna().unique()))
        if n_haplo == 1:
            sbs_dups['HAPLOTAG'] = list(pol_sbs['HAPLOTAG'].dropna().unique())[0]
            sbs_dups['type'] = 'assembly'
            sbs_dups['POLISHED_HAPLO'] = list(pol_sbs['POLISHED_HAPLO'].dropna().unique())[0]
            return pd.concat([pol_sbs, sbs_dups], axis=0)
        else:
            h1_size = pol_sbs.loc[pol_sbs['HAPLOTAG'] == 1, 'POLISHED_HAPLO'].unique()
            h2_size = pol_sbs.loc[pol_sbs['HAPLOTAG'] == 2, 'POLISHED_HAPLO'].unique()
            target = list(sbs_dups['LEN_SEQUENCE_FOR_TRF'])
            haplo = []; size = []
            for x in target:
                tmp_haplo, tmp_size = assignHaplotag_asm(h1_size, h2_size, x)
                haplo.append(tmp_haplo); size.append(tmp_size)
            sbs_dups['type'] = 'assembly'
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
        best_motif = high_coverage['UNIFORM_MOTIF'].iloc[[0]]
        copies = high_coverage['COPIES_TRF'].iloc[[0]]
        start = high_coverage['START_TRF'].iloc[[0]]
        end = high_coverage['END_TRF'].iloc[[0]]
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
def sampleMotifs(r, all_sbs, reference_motif_dic, haplo):
    # subset of data
    haplo_data = all_sbs[all_sbs['HAPLOTAG'] == haplo].copy()
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
    # write sequence file
    df_seq.to_csv(seq_file, header=True, index=False, sep='\t')
    return

# 1. arguments
inpf = sys.argv[1]
outd = sys.argv[2]
n_cpu = int(sys.argv[3])
thr_mad = float(sys.argv[4])

# 2. read data
print('** Read data')
data = pd.read_csv(inpf, sep='\t')

# 3. adjust the motif -- merge the same motifs
print('** Adjust motifs')
all_motifs = data['TRF_MOTIF'].dropna().unique()
main_motifs = [permutMotif(motif) for motif in all_motifs]
motifs_df = pd.DataFrame({'motif' : all_motifs, 'UNIFORM_MOTIF' : main_motifs})
data = pd.merge(data, motifs_df, left_on='TRF_MOTIF', right_on='motif', how='left')
data['UNIQUE_NAME'] = data.apply(lambda row: row['READ_NAME'] + '___' + row['SAMPLE_NAME'] + '___' + row['REGION'] + '___' + str(row['LEN_SEQUENCE_FOR_TRF']), axis = 1)

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

# 5. add TR size
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
    print('**** %s' %(s))
    sbs = data_nodup[(data_nodup['SAMPLE_NAME'] == s)]
    sbs_dups = dup_df[(dup_df['SAMPLE_NAME'] == s)]
    pool = multiprocessing.Pool(processes=n_cpu)
    haplo_fun = partial(haplotyping, s = s, thr_mad = thr_mad, data_nodup = sbs, type = 'asm', dup_df = sbs_dups, reference_motif_dic = reference_motif_dic, intervals = intervals)
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
seq_file = '%s/sample.seq.txt' %(outd)
vcf_file = '%s/sample.vcf' %(outd)
writeOutputs(df_vcf, df_seq, seq_file, vcf_file)

