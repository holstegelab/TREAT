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
        print('****** done %s%% of the regions' %(intervals.index(r)*5+5), end='\r')
    # data of interest
    sbs = data_nodup[data_nodup['REGION'] == r]
    # exclude nas
    sbs = sbs.dropna(subset=['LEN_SEQUENCE_FOR_TRF'])
    # check if there are rows
    if sbs.shape[0] >0:
        if type == 'asm':
            # identify haplotypes
            phased_sbs = assemblyBased_size(sbs)
            # polish haplotypes
            pol_sbs = polishHaplo_asm(phased_sbs)
            # add duplicates
            all_sbs = addDups(pol_sbs, dup_df, s, r)
            # finally look at the motif
            final_sbs_h1 = sampleMotifs(r, all_sbs, reference_motif_dic, 1)
            final_sbs_h2 = sampleMotifs(r, all_sbs, reference_motif_dic, 2)
            final_sbs = pd.concat([final_sbs_h1, final_sbs_h2], axis=0)
            return final_sbs
        elif type == 'reads_spanning':
            return None

# function to call haplotypes for assembly
def assemblyBased_size(sbs):
    n_contigs = sbs.shape[0]
    # check number of contigs
    if n_contigs == 1:
        # homozygous call
        sbs['HAPLOTAG'] = 1; sbs['type'] = 'assembly'
    elif n_contigs == 2:
        # heterozygous call
        sbs['HAPLOTAG'] = [1, 2]; sbs['type'] = 'assembly'
    elif sbs.shape[0] >2:
        print('more than 2 haps')
    return sbs

# function to polish haplotypes
def polishHaplo_asm(phased_sbs):
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
        print('more than 1 haplo and more than 2 contigs')
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
        print('****** done %s%% of the regions' %(intervals.index(r)*5+5), end='\r')
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
    return sbs

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
                copies = len(list(upd_motif_index))/len(best_motif)
            else:
                # here the more tricky case: motifs are different
                not_in_best = [x for x in tmp_motif_index if x not in best_motif_index]
                # check if these are consecutive numbers
                if is_consecutive(not_in_best) == True:
                    not_in_best_copies = len(not_in_best)/len(tmp_motif)
                    best_motif = '%s+%s' %(best_motif, tmp_motif)
                    copies = '%s+%s' %(best_motif_copies, not_in_best_copies)
                    best_motif_index = range(min(best_motif_index.start, tmp_motif_index.start), max(best_motif_index.stop, tmp_motif_index.stop))
                else:
                    print(r)
                    copies = best_motif_copies
                    best_motif_index = best_motif_index
        else:
            copies = best_motif_copies; best_motif_index = best_motif_index
    return best_motif, copies, best_motif_index

# function to check whether a sequence is consecutive
def is_consecutive(nums):
    # check if the list is empty or has a single element
    if len(nums) <= 1:
        return True    
    # check if the difference between the maximum and minimum element
    # is equal to the length of the list minus 1
    return max(nums) - min(nums) == len(nums) - 1 and set(nums) == set(range(min(nums), max(nums)+1))

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

# function to write vcf file
def writeVcf(sam_info, ref_info, outf, r, samples):
    # prepare data
    chrom, start, end = [r.split(':')[0]] + r.split(':')[-1].split('-')
    ref_motif, ref_len, ref_copies = ref_info
    info_field = '%s;%s' %(ref_motif, ref_copies)
    format_field = 'QC;GT;MOTIF;CN;CN_REF;DP'
    outf.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %(chrom, start, r, ref_len, '.', '.', 'PASS', info_field, format_field))
    for s in samples:
        outf.write('\t')
        s_info = sam_info[sam_info['SAMPLE_NAME'] == s]
        sam_gt = '|'.join(list(sam_info['POLISHED_HAPLO'].map(str))) if sam_info.shape[0] >1 else list(sam_info["POLISHED_HAPLO"].map(str).apply(lambda x: x + "|" + x))[0]
        sam_mot = '|'.join(list(sam_info['CONSENSUS_MOTIF'].map(str))) if sam_info.shape[0] >1 else list(sam_info["CONSENSUS_MOTIF"].map(str).apply(lambda x: x + "|" + x))[0]
        sam_cop = '|'.join(list(sam_info['CONSENSUS_MOTIF_COPIES'].map(str))) if sam_info.shape[0] >1 else list(sam_info["CONSENSUS_MOTIF_COPIES"].map(str).apply(lambda x: x + "|" + x))[0]
        sam_cop_ref = '|'.join(list(sam_info['REFERENCE_MOTIF_COPIES'].map(str))) if sam_info.shape[0] >1 else list(sam_info["REFERENCE_MOTIF_COPIES"].map(str).apply(lambda x: x + "|" + x))[0]
        sam_depth = '|'.join(list(sam_info["READ_NAME"].str.split("_").str[1])) if sam_info.shape[0] >1 else list(sam_info["READ_NAME"].str.split("_").str[1].apply(lambda x: x + "|" + x))[0]
        sample_field = '%s;%s;%s;%s;%s;%s' %('PASS_ASM', sam_gt, sam_mot, sam_cop, sam_cop_ref, sam_depth)
        outf.write('%s' %(sample_field))
    outf.write('\n')
    return 

# function to write sequence file
def writeSeq(sam_info, outf, ref_info):
    # add depth and reference motif
    sam_info['DEPTH'] = sam_info["READ_NAME"].str.split("_").str[1]
    sam_info['MOTIF_REF'] = ref_info[0]
    keep_cols = ['READ_NAME', 'HAPLOTAG', 'REGION', 'PASSES', 'READ_QUALITY', 'LEN_SEQUENCE_FOR_TRF', 'START_TRF', 'END_TRF', 'type', 'SAMPLE_NAME', 'POLISHED_HAPLO', 'DEPTH', 'CONSENSUS_MOTIF', 'CONSENSUS_MOTIF_COPIES', 'MOTIF_REF', 'REFERENCE_MOTIF_COPIES', 'SEQUENCE_WITH_PADDINGS', 'SEQUENCE_FOR_TRF']
    # subset dataframe to keep only selected columns
    df_subset = sam_info.loc[:, keep_cols]
    df_subset.to_csv(outf, index=False, header=False, sep="\t")
    return

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
    outf.write('##fileformat=VCFv4.2\n##INFO=<ID=REFERENCE_INFO,Number=2,Type=String,Description="Motif observed in the reference genome (GRCh38), and relative number of motif repetitions."\n##FORMAT=<ID=QC,Number=1,Type=String,Description="Quality summary of TREAT genotyping. PASS_BOTH: genotype agreed between reads-spanning and assembly. PASS_RSP: genotype from reads-spanning. PASS_ASM: genotype from assembly."\n##FORMAT=<ID=GT,Number=2,Type=String,Description="Phased size of the tandem repeats. H1_size | H2_size"\n##FORMAT=<ID=MOTIF,Number=2,Type=String,Description="Phased consensus motif found in the sample. H1_motif | H2_motif"\n##FORMAT=<ID=CN,Number=2,Type=String,Description="Phased number of copies of the motif found in the sample. H1_copies | H2_copies"\n##FORMAT=<ID=CN_REF,Number=2,Type=String,Description="Phased estimation of the reference motif as found in the sample. H1_motif_ref | H2_motif_ref"\n##FORMAT=<ID=DP,Number=1,Type=String,Description="Phased depth found in the sample. H1_depth | H2_depth"\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' %('\t'.join(samples)))
    outf.close()
    return    

# function to write header of the sequence file
def writeSeqheader(seq_file, samples):
    # open file
    outf = open(seq_file, 'w')
    # write header
    outf.write('READ_NAME\tHAPLOTYPE\tREGION\tPASSES\tREAD_QUALITY\tLENGTH_SEQUENCE\tSTART_TRF\tEND_TRF\tTYPE\tSAMPLE\tPOLISHED_HAPLO\tDEPTH\tCONSENSUS_MOTIF\tCONSENSUS_MOTIF_EST_COPIES\tCONSENSUS_MOTIF_REF\tCONSENSUS_MOTIF_REF_EST_COPIES\tSEQUENCE_WITH_PADDINGS\tSEQUENCE_FOR_TRF\n')
    outf.close()
    return

# function to write outputs
def writeOutputs(sample_res, seq_file, vcf_file, reference_motif_dic):
    # extract regions
    all_regions = list(sample_res['REGION'].dropna().unique())
    all_samples = list(sample_res['SAMPLE_NAME'].dropna().unique())
    # write header of vcf file
    writeVCFheader(vcf_file, all_samples)
    outvcf = open(vcf_file, 'a')
    # write header of sequence file
    writeSeqheader(seq_file, all_samples)
    outseq = open(seq_file, 'a')
    # loop across samples and regions
    for r in all_regions:
        # get reference info and sample info
        ref_info = reference_motif_dic[r]
        sam_info = sample_res[sample_res['REGION'] == r].copy()
        # write to vcf
        writeVcf(sam_info, ref_info, outvcf, r, all_samples)
        # write to sequence file
        writeSeq(sam_info, outseq, ref_info)
    outvcf.close()
    outseq.close()
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
motif_ref = pd.concat(motif_res, axis=0)
# probably good to create a dictionary as well, will be faster to access it
reference_motif_dic = {row["REGION"]: [row["CONSENSUS_MOTIF"], row["POLISHED_HAPLO"], row['CONSENSUS_MOTIF_COPIES']] for _, row in motif_ref.iterrows()}

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
    haplo_results_combined = pd.concat(haplo_results, axis=0)
    sample_res.append(haplo_results_combined)
sample_res = pd.concat(sample_res, axis=0)

# 7. write outputs: vcf file and raw sequences
print('** Producing outputs: VCF file and table with sequences                        ')
seq_file = '%s/sample.seq.txt' %(outd)
vcf_file = '%s/sample.vcf' %(outd)
writeOutputs(sample_res, seq_file, vcf_file, reference_motif_dic)

