#######################################################################
#######################################################################
# reads spanning is likely the way to go. The reason is that assembly still doesn't help with the haplotypes estimation -- this is now with all the regions -- updated 2023-10-31

# Libraries
import re
import time
import pandas as pd
import math
import pysam
import numpy as np
from scipy.stats import binom_test
import os
from collections import Counter

# Functions
# function to read all alleles from database
def readDatabase(database_path):
    fname = open(database_path).readlines()
    i = 0
    # set variables
    alleles = []
    alleles_ref = {}
    while i<len(fname):
        # rstrip and split line
        line = fname[i].rstrip().split()
        # check the size: if size 1 it can either be the name of a region, or a sequence. check if there are numbers to distinguish the two
        if len(line) <=1:
            pass
        # if first element is Marker, then we can pass
        elif line[0] == 'Marker':
            pass
        # otherwise this line contains alleles. check whether the size is 2 (allele_name and allele) or 3 (allele_name, allele and sequence)
        elif line[0] == 'REF':
            alleles_ref[fname[i-2].rstrip()] = line[-1]
        else:
            allele_name, allele = line[0:2]
            if allele_name == 'DYS461+DYS460':
                allele_name = 'DYS461-460'
            elif allele_name == 'Y-GATA-H4':
                allele_name = 'GATAH4'
            alleles.append([allele_name, allele])
        i += 1
    # convert to dataframe
    alleles_df = pd.DataFrame(alleles, columns=['allele_name', 'allele'])
    # exclude the reference alleles (likely we have them as well)
    alleles_df = alleles_df[alleles_df['allele_name'] != 'REF'].copy()
    return alleles_df, alleles_ref

# Function that given a motif, extend it such as AG[2] becomes AGAG -- examples: chrY:20471984-20472026 -- chr4:154587714-154587823
def extendMotif(input_str):
    # extract sequence without modification
    input_str = input_str.split('_')[0]
    # extract sequence
    i = 0; motif = ''; num_rep = ''; output_str = ''
    while i < len(input_str):
        if input_str[i] == '[':
            pass
        elif input_str[i] == ']':
            output_str = output_str + motif * int(num_rep)
            num_rep = ''
            motif = ''
        elif input_str[i] in [str(x) for x in range(10)]:
            num_rep = num_rep + input_str[i]
        else:
            motif = motif + input_str[i]
        i += 1
    return output_str

# function to match on single reads
def manageRs(sample_info_rs, pattern, allele_name):
    res = []
    # extract all haplotypes by treat
    haplo = list(set(list(sample_info_rs['HAPLOTAG'])))
    # check haplotypes
    for h in haplo:
        # subset of sequences for that haplotype
        tmp = sample_info_rs[sample_info_rs['HAPLOTAG'] == h].copy()
        haplo_value = list(set(list(tmp['POLISHED_HAPLO'])))[0]
        # extract sequences
        seqs = list(tmp['SEQUENCE_WITH_PADDINGS'])
        # define counter for the sequences that match the allele
        counter = 0
        # iterate over the sequences
        for seq in seqs:
            start_index = seq.find(pattern)
            # check if motif was found
            if start_index != -1:
                counter += 1
        if counter >0:
            res.append([allele_name, len(pattern), counter, start_index, h, len(seqs), haplo_value])
    return res

# function to prioritize motif
def prioritizeMot(df, type):
    newdf = []
    # haplotypes
    haplo = list(set(list(df['Haplotype_treat'])))
    for h in haplo:
        tmp = df[df['Haplotype_treat'] == h].copy()
        if type == 'asm':
            tmp['Priority'] = tmp['Length_allele_database'].apply(lambda x: 'Priority' if x == tmp['Length_allele_database'].max() else 'No')
        elif type == 'rs':
            weight_size = 0.7
            weight_cov = 0.3
            tmp['Score'] = tmp['Length_allele_database']*weight_size + tmp['Coverage']*weight_cov
            # sort by score
            tmp = tmp.sort_values(by='Score')
            #tmp['Priority'] = tmp['Score'].apply(lambda x: 'Priority' if x == tmp['Score'].max() else 'No')
            tmp['Priority'] = np.where(tmp['Score'] == tmp['Score'].max(), 'Priority', np.where(tmp['Score'] == tmp['Score'].nlargest(2).iloc[-1], 'Second_Priority', 'No'))
        newdf.append(tmp)
    combined = pd.concat(newdf)
    return combined

# Define function to apply to depth column for minimum support per match from database
def check_value(x, min_support):
    if x >= min_support:
        return 'keep'
    else:
        return 'low_coverage'

# Define function to apply to depth column for minimum support per allele
def check_imbalance(x):
    if 1 in x and 2 in x:
        # count how many reads with alleles matched are from haplotype 1 and 2
        h1_count = x.count(1)
        h2_count = x.count(2)
        # calculate imbalance of alleles
        allele_imbalance = min(h1_count, h2_count) / max(h1_count, h2_count)
        # make binomial test to see whthere there's a significant imbalance
        num_successes = h1_count
        num_trials = h1_count + h2_count
        null_prob = 0.5
        p_value = binom_test(num_successes, num_trials, null_prob)
        if p_value <0.05:
            return 'significant_allele_imbalance'
        else:
            return 'alleles_balanced'
    else:
        return 'alleles_balanced'

# Define function to apply to check validation
def check_validation(x, validation):
    if x in validation:
        return 'match'
    else:
        return 'mismatch'

# Function to do some filtering on the motif matched: there should be a minimum support of 50% of all reads in the haplotype for a given motif
def minSupportMotif(df_rs, min_support):
    qced_data = []
    # iterate over rows of haplotype reads
    for index, row in df_rs.iterrows():
        motif_coverage = row['Coverage']
        hap_coverage = row['Depth']
        # check if motif is seen in at least 50% of reads
        if motif_coverage/hap_coverage >= min_support:
            qced_data.append('keep')
        elif motif_coverage/hap_coverage > 1:
            print('weird, check!!')
            break
        else:
            qced_data.append('remove')
    # combined results
    df_rs['QC_motif_minSupport'] = qced_data
    return df_rs

# maybe first do a round of QC based on total coverage and allelic imbalance, then go ahead with the remaining regions
def general_QC_reads(tmp, coverage, region, gender):
    # autosomal or sex chromosoms
    region_type = 'sex' if 'chrY' in region else 'autosomal'
    if region_type == 'sex':
        pass
    else:
        region_type = 'sex' if 'chrX' in region else 'autosomal'
    if region_type == 'autosomal' and tmp.shape[0] < coverage:
        return 'low_coverage: %s' %(region)
    else:
        # list haplotypes
        haplo_list = list(set(list(tmp['HAPLOTAG'])))
        all_haplo = list(tmp['HAPLOTAG'])
        # check if the chromosome is sexual and adapt accordingly
        if region_type == 'sex' and 'chrX' in region and gender == 'male' and len(list(set(all_haplo))) >1:
            all_haplo = [0 for x in all_haplo]          # this case is a male, chrX --> there should be only 1 haplotype (sometimes there're more just because of noise)
            haplo_list = [0]
        elif region_type == 'sex' and 'chrY' in region and gender == 'male' and len(list(set(all_haplo))) >1:
            all_haplo = [0 for x in all_haplo]          # this case is a male, chrY --> there should be only 1 haplotype (sometimes there're more just because of noise)
            haplo_list = [0]
        check_haplo = []
        # check coverage per haplotype
        for haplo in haplo_list:
            count_haplo = all_haplo.count(haplo)
            if count_haplo >= (coverage/2):
                check_haplo.append('ok') 
            else:
                check_haplo.append('no') 
        # in case coverage is ok, cehck allelic imbalance
        if 'no' not in check_haplo:
            imb = check_imbalance(haplo_list)
            return imb
        else:
            return 'low_coverage: %s' %(region)

# function to do sex check
def sex_check(sample, data_rs_info):
    # subset to sample of interest
    sb = data_rs_info[data_rs_info['SAMPLE_NAME'] == sample]
    # check if there are any chrY reads
    sb_regions = list(set(list(sb['REGION'])))
    gender = 'male' if True in list(sb['REGION'].str.contains('chrY')) else 'female'
    return gender

# function to clean database to exclude the SNP modifications (this will be checked later)
def clean_database(database_df):
    # list regions
    tmp_regions = list(set(list(database_df['allele_name'])))
    # create dataframe for output
    cleaned_db = []
    modification_db = []
    # iterate over regions
    for r in tmp_regions:
        sb = database_df[database_df['allele_name'] == r].copy()
        # convert to list
        sb_list = sb.values.tolist()
        # iterate over list
        for x in sb_list:
            tmp_allele = x[-1].split('_')
            if len(tmp_allele) >2:
                tmp_allele_nomod = '_'.join(tmp_allele[0:2])
                tmp_allele_nomod2 = tmp_allele_nomod.replace('[]', '')
                cleaned_db.append([x[0], tmp_allele_nomod2])
                modification_db.append(x)
            else:
                tmp_allele_cleaned = '_'.join(tmp_allele[0:2]).replace('[]', '')
                cleaned_db.append([x[0], tmp_allele_cleaned]) 
    # convert lists to dataframes
    cleaned_df = pd.DataFrame(cleaned_db, columns = ['allele_name', 'allele'])
    modification_df = pd.DataFrame(modification_db, columns = ['allele_name', 'allele'])
    # remove duplicated alleles
    cleaned_df['locus'] = cleaned_df['allele_name'] + '_' + cleaned_df['allele']
    modification_df['locus'] = modification_df['allele_name'] + '_' + modification_df['allele']
    cleaned_df_nodup = cleaned_df.drop_duplicates(subset='locus').copy()
    modification_df_nodup = modification_df.drop_duplicates(subset='locus').copy()
    # delete locus column
    cleaned_df_nodup.drop('locus', axis=1, inplace=True)
    modification_df_nodup.drop('locus', axis=1, inplace=True)
    return cleaned_df_nodup, modification_df_nodup

# function to check whether multiple alleles with the same size are prioritized -- maybe it's a snp
def checkSameSize(priority_df_rs):
    tocheck = []
    all_samples = list(set(list(priority_df_rs['Sample'])))
    all_reg = list(set(list(priority_df_rs['Region'])))
    for s in all_samples:
        for r in all_reg:
            tmp = priority_df_rs[(priority_df_rs['Sample'] == s) & (priority_df_rs['Region'] == r)]
            haplotypes = list(set(list(tmp['Haplotype_treat'])))
            for h in haplotypes:
                tmptmp = tmp[tmp['Haplotype_treat'] == h]
                tmptmp_mot = [x.split('_')[1] for x in list(tmptmp['Allele_name'])]
                tmptmp_mot_ext = [extendMotif(x) for x in tmptmp_mot]
                tmptmp_mot_ext_len = [len(x) for x in tmptmp_mot_ext]
                prior_index = list(tmptmp['Priority']).index('Priority')
                if len(list(set(tmptmp_mot_ext_len))) < len(tmptmp_mot_ext_len):
                    # check indexes -- has to be duplicate with the prioritized motif
                    duplicate_indices = [[i, j] for i, x in enumerate(tmptmp_mot_ext_len) for j, y in enumerate(tmptmp_mot_ext_len) if i < j and x == y]
                    for x in duplicate_indices:
                        if prior_index in x:
                            tmpdup = [x.split('_')[1] for x in list(tmptmp['Allele_name'])]
                            if tmpdup[x[0]] != tmpdup[x[1]]:
                                tocheck.append([s, r])
    return tocheck

# MAIN
# define path
PATH = 'path/to/TREAT/output/folder'
# read database of alleles
database_df_all, alleles_ref = readDatabase('/path/to/database_allele.txt')
database_df, database_modifications = clean_database(database_df_all)
database_df['allele_name'] = database_df['allele_name'].replace('GATAH4', 'Y-GATA-H4')
# read data with sequences
data_rs = pd.read_csv('%s/sample.raw.txt.gz' %(PATH), sep='\t')

# read data with motif structures and genomic regions
data_motif_genomicRegions = pd.read_csv('/path/to/forensic_hg38.extended.bed', sep='\t', header=None)
# add locus column and merge the two data
data_motif_genomicRegions.columns = ['chrom', 'start', 'end', 'id']
data_motif_genomicRegions['locus'] = data_motif_genomicRegions.apply(lambda row: f"{row['chrom']}:{row['start']}-{row['end']}", axis=1)
# merge bed file with data
data_rs_info = pd.merge(data_rs, data_motif_genomicRegions, left_on='REGION', right_on='locus')
data_rs_info['id'] = data_rs_info['id'].replace('GATAH4', 'Y-GATA-H4')
data_rs_info.shape
data_rs_info.head()

# parameters for qc
min_support = 0.30
coverage = [6, 8, 10, 12]
failed_dic = {}

# iterate over coverage thresholds
for cov in coverage:
    print('*** Coverage threshold is: %s' %(cov))
    # extract all regions
    all_regions = list(set(list(data_rs_info['locus'])))
    all_samples = list(set(list(data_rs_info['SAMPLE_NAME'])))
    # qc on data based on coverage and allelic imbalance
    data_rs_info_qc = []
    data_rs_info_failed_qc = []
    for sample in all_samples:
        # check if sample is male or female
        gender = sex_check(sample, data_rs_info)
        print('** %s' %(sample))
        for region in all_regions:
            # extract data
            tmp = data_rs_info[(data_rs_info['SAMPLE_NAME'] == sample) & (data_rs_info['REGION'] == region)]
            # exclude duplicates based on read name
            tmp = tmp.drop_duplicates(subset='READ_NAME').copy()
            # exclude any read without the haplotag (deviating too much)
            tmp = tmp.dropna(subset=['HAPLOTAG'])
            # run qc
            qc_result = general_QC_reads(tmp, cov, region, gender)
            # check qc results
            if qc_result == 'alleles_balanced':
                data_rs_info_qc.append(tmp)
            else:
                data_rs_info_failed_qc.append(tmp)
    # combine results
    data_rs_info_qc = pd.concat(data_rs_info_qc, axis=0)
    data_rs_info_qc.shape
    data_rs_info_failed_qc = pd.concat(data_rs_info_failed_qc, axis=0)
    data_rs_info_failed_qc.shape
    # check regions that failed qc
    data_rs_info_failed_qc_list = data_rs_info_failed_qc.values.tolist()
    failed_dic[cov] = {}
    for x in data_rs_info_failed_qc_list:
        failed_sample, failed_region = x[0:2]
        if failed_sample in failed_dic[cov].keys():
            if failed_region in failed_dic[cov][failed_sample]:
                pass
            else:
                failed_dic[cov][failed_sample].append(failed_region)
        else:
            failed_dic[cov][failed_sample] = []
            failed_dic[cov][failed_sample].append(failed_region)
    # get the list of samples and regions after qc
    all_regions = list(set(list(data_rs_info_qc['locus'])))
    all_samples = list(set(list(data_rs_info_qc['SAMPLE_NAME'])))
    # container for results
    all_res_rs = []
    priority_res_rs = []
    missing_regions_in_database = []
    # iterate over the samples and regions
    for sample in all_samples:
        print('** working on sample %s' %(sample))
        for region_interest in all_regions:
            # get the name of the region
            region_id = data_rs_info_qc[data_rs_info_qc['locus'] == region_interest].copy()
            region_id = list(set(list(region_id['id'])))[0].split(';')[1].split('=')[1]
            # get the motifs from the database
            if region_id not in list(database_df['allele_name']):
                missing_regions_in_database.append(region_id)
            database_alleles = database_df[database_df['allele_name'] == region_id]
            allele_list = list(database_alleles['allele'])
            # get corresponding data for the sample
            sample_info_rs = data_rs_info_qc[(data_rs_info_qc['SAMPLE_NAME'] == sample) & (data_rs_info_qc['locus'] == region_interest)].copy()
            # save the number of haplotypes here
            haplo_before_qc = list(set(list(sample_info_rs['HAPLOTAG'])))
            # define container for results
            results_rs = []
            if len(allele_list) >0 and sample_info_rs.shape[0] >0:
                # extract sequences
                sample_seqs_rs = list(sample_info_rs['SEQUENCE_WITH_PADDINGS'])
                #tmp_h1 = sample_info_rs[sample_info_rs['HAPLOTYPE'] == 1]
                #tmp_h2 = sample_info_rs[sample_info_rs['HAPLOTYPE'] == 2]
                #tmp_h1_seq = list(tmp_h1['SEQUENCE_WITH_PADDINGS'])
                #tmp_h2_seq = list(tmp_h2['SEQUENCE_WITH_PADDINGS'])
                # then iterate over the alleles and check occurrence in the sequence
                for allele in allele_list:
                    # restructure the allele
                    try:
                        allele_seq = extendMotif('_'.join(allele.split('_')[1::]))
                        res_rs = manageRs(sample_info_rs, allele_seq, allele)
                        results_rs.extend(res_rs)
                    except:
                        pass
            # convert to df
            if len(results_rs) >0:
                df_rs = pd.DataFrame(results_rs, columns=['Allele_name', 'Length_allele_database', 'Coverage', 'Index_start', 'Haplotype_treat', 'Depth', 'Haplotype_size'])
                df_rs['Region'] = region_interest
                df_rs['Region_ID'] = region_id
                df_rs['Type'] = 'reads_spanning'
                df_rs['Sample'] = sample
                if df_rs.shape[0] >0:
                    # then check if there's enough support for the allele matches from the database
                    df_rs_qc = minSupportMotif(df_rs, min_support)
                    # when we remove motifs because of lack of support, it can be that one haplotype is not represented anymore. check that and in case, exclude from further analysis
                    df_rs_qc = df_rs_qc[df_rs_qc['QC_motif_minSupport'] == 'keep']
                    haplo_after_qc = list(set(list(df_rs_qc['Haplotype_treat'])))
                    check_dropout = [x in haplo_after_qc for x in haplo_before_qc]
                    if False not in check_dropout:
                        #df_rs['Support'] = df_rs['Coverage'].apply(lambda x: check_value(x, min_support))
                        #df_rs_qc = df_rs[df_rs['Support'] == 'keep']
                        if df_rs_qc.shape[0] >0:
                            df_rs_mot = prioritizeMot(df_rs_qc, 'rs')
                            df_priority_rs = df_rs_mot[df_rs_mot['Priority'].isin(['Priority', 'Second_Priority'])]
                            priority_res_rs.append(df_priority_rs)
                    else:
                        print('%s at %s' %(sample, region_interest))
                all_res_rs.append(df_rs)
    missing_regions_in_database = list(set(missing_regions_in_database))
    # join dataframes
    priority_df_rs = pd.concat(priority_res_rs, axis=0)
    all_df_rs = pd.concat(all_res_rs, axis=0)
    # check if sometimes the size of the matching allele is the same but there's a snp maybe
    regions_to_check = checkSameSize(priority_df_rs)
    # check these regions afterwards -- turns out only 1 case need a modification (PH6 -- chr1...)
    # read validation data
    validation_data = pd.read_csv('/project/holstegelab/Share/nicco/workspaces/projects/automatic_pipeline/Peter_forensic/new_data/validation_data.txt', sep="\t")
    # compare
    results_comparison_rs = []
    for sample in all_samples:
        # extract samples and relative regions
        tmp_sample_rs = priority_df_rs[priority_df_rs['Sample'] == sample]
        tmp_regions = list(tmp_sample_rs['Region'].dropna().unique())
        for region_interest in tmp_regions:
            # get the region id
            region_id = tmp_sample_rs[tmp_sample_rs['Region'] == region_interest]
            region_id = list(set(list(region_id['Region_ID'])))[0]
            # get validation data for the sample and the region
            tmp_validation = validation_data[(validation_data['Markername'] == region_id)]
            tmp_validation_alleles = [x for x in list(tmp_validation['Allelename_%s' %(sample.split('.')[0])]) if isinstance(x, str) or not math.isnan(x)]
            # added on 30-11-23 to exclude modifications from the matching -- modifications will be checked later
            tmp_validation_alleles = ['_'.join(x.split('_')[0:2]) for x in tmp_validation_alleles]
            # get prioritized alleles
            tmp_prior_rs = tmp_sample_rs[tmp_sample_rs['Region'] == region_interest].copy()
            tmp_prior_rs = tmp_prior_rs[tmp_prior_rs['Priority'] == 'Priority'].copy()
            tmp_prior_alleles_rs = list(tmp_prior_rs['Allele_name'])
            if len(tmp_validation_alleles) >0:
                if len(tmp_prior_alleles_rs) >0:
                    validated = []
                    tmp_prior_rs['Validation'] = tmp_prior_rs['Allele_name'].apply(lambda x: check_validation(x, tmp_validation_alleles))
                    # then check rs
                    match = 0
                    for x in tmp_validation_alleles:
                        if x in tmp_prior_alleles_rs:
                            match += 1
                    match_pc = match/len(tmp_validation_alleles)
                    tmp_prior_rs['Match_fraction'] = match_pc
                    results_comparison_rs.append(tmp_prior_rs)
    results_comparison_rs_df = pd.concat(results_comparison_rs, axis=0)
    # check modifications like snps here
    #results_comparison_rs_df_modification = checkAlleleModification(results_comparison_rs_df, modification_db)
    # write output
    outname = '/project/holstegelab/Share/nicco/workspaces/projects/automatic_pipeline/Peter_forensic/new_data/treat_otter/new_data_reads/20231211_haplodev_05/comparison_reads_validation_haplo0.05_cov%s.txt' %(cov)
    results_comparison_rs_df.to_csv(outname, sep='\t', index=False, decimal=',')
    outname2 = '/project/holstegelab/Share/nicco/workspaces/projects/automatic_pipeline/Peter_forensic/new_data/treat_otter/new_data_reads/20231211_haplodev_05/cov%s_failed_qc.txt' %(cov)
    data_rs_info_failed_qc.to_csv(outname2, sep='\t', index=False, decimal=',')

# AmelogeninX we should look specifically
amelox = data_rs[data_rs['REGION'] == 'chrX:11296893-11296955']
ameloy = data_rs[data_rs['REGION'] == 'chrY:6869852-6869958']
for cov in coverage:
    print('\n*** Coverage threshold is: %s' %(cov))
    # extract all regions
    all_samples_x = list(set(list(amelox['SAMPLE_NAME'])))
    all_samples_y = list(set(list(ameloy['SAMPLE_NAME'])))
    # qc on data based on coverage and allelic imbalance
    for sample in list(set(all_samples_x + all_samples_y)):
        # extract data
        tmp_x = amelox[(amelox['SAMPLE_NAME'] == sample)]
        tmp_y = ameloy[(ameloy['SAMPLE_NAME'] == sample)]
        # exclude duplicates based on read name
        tmp_x = tmp_x.drop_duplicates(subset='READ_NAME').copy()
        tmp_y = tmp_y.drop_duplicates(subset='READ_NAME').copy()
        # exclude any read without the haplotag (deviating too much)
        tmp_x = tmp_x.dropna(subset=['HAPLOTAG'])
        tmp_y = tmp_y.dropna(subset=['HAPLOTAG'])
        # check if it's homozygous or heterozygous
        tmp_geno = 'homo' if tmp_y.shape[0] == 0 else 'hetero'
        cov_thr = cov/2 if tmp_geno == 'hetero' else cov
        # check amelo x first -- check if there are enough reads and if the deletion exists
        gender_match_x = []
        if tmp_x.shape[0] >= cov_thr:
            #print('** N=%s reads for %s' %(tmp_x.shape[0], sample))
            #print('** N=%s haplotypes for %s' %(len(set(tmp_x['HAPLOTAG'])), sample))
            # get sequences for haplotypes
            seqs = list(tmp_x['SEQUENCE_WITH_PADDINGS'])
            # before it was AGATGTT
            for x in seqs:
                if 'CACTTT' in x:
                    gender_match_x.append('female')
                else:
                    gender_match_x.append('male')
        # then ameloY
        gender_match_y = []
        if tmp_y.shape[0] >= cov_thr:
            #print('** N=%s reads for %s' %(tmp_y.shape[0], sample))
            #print('** N=%s haplotypes for %s' %(len(set(tmp_y['HAPLOTAG'])), sample))
            # get sequences for haplotypes
            seqs = list(tmp_y['SEQUENCE_WITH_PADDINGS'])
            # before it was AGATGTT
            for x in seqs:
                if 'CACTTT' in x:
                    gender_match_y.append('male')
                else:
                    gender_match_y.append('female')
        # call by majority
        counter_x = Counter(gender_match_x)
        gender_call_x = counter_x.most_common(1)[0][0] if len(counter_x) >0 else 'None'
        counter_y = Counter(gender_match_y)
        gender_call_y = counter_y.most_common(1)[0][0] if len(counter_y) >0 else 'None'
        # compare ameloX and ameloY
        if gender_call_x == 'None' or gender_call_y == 'None':
            if gender_call_x == 'None' and gender_call_y == 'None':
                print('%s has no AmeloX and AmeloY' %(sample))
            elif gender_call_x == 'None':
                print('%s has no AmeloX -- [%s]' %(sample, len(gender_match_y)))
            else:
                print('%s has no AmeloY -- [%s]' %(sample, len(gender_match_x)))
        else:
            if gender_call_x == gender_call_y:
                print('%s is %s -- [%s - %s]' %(sample, gender_call_x, len(gender_match_x), len(gender_match_y)))
            else:
                print('%s is discordant%s' %(sample, gender_call_x))

# manual check of some alleles where we see a mismatch
# PH6 - DYS570
sb = data_rs[(data_rs['SAMPLE_NAME'] == 'PH6.ontarget.sup.hg38') & (data_rs['REGION'] == 'chrY:6993188-6993260')]
sb_seq = list(sb['SEQUENCE_WITH_PADDINGS'])
tmp1 = validation_data[validation_data['Markername'] == 'DYS570']
validated_motif = extendMotif(list(tmp1['Allelename_PH6'])[0].split('_')[1])
match_count = 0
for s in sb_seq:
    if validated_motif in s:
        match_count += 1
print('%s out of %s sequences match the right allele' %(match_count, len(sb_seq)))

# PH6 - FGA
sb = data_rs[(data_rs['SAMPLE_NAME'] == 'PH6.ontarget.sup.hg38') & (data_rs['REGION'] == 'chr4:154587714-154587823')]
sb_seq = list(sb['SEQUENCE_WITH_PADDINGS'])
tmp1 = validation_data[validation_data['Markername'] == 'FGA']
validated_motifs = [extendMotif(x.split('_')[1]) for x in list(tmp1['Allelename_PH6'])]
for motif in validated_motifs:
    match_count = 0
    for s in sb_seq:
        if motif in s:
            match_count += 1
    print('%s out of %s sequences match the right allele (%s)' %(match_count, len(sb_seq), motif))

# PH7 - DYS481
sb = data_rs[(data_rs['SAMPLE_NAME'] == 'PH7.ontarget.sup.hg38') & (data_rs['REGION'] == 'chrY:8558313-8558408')]
sb_seq = list(sb['SEQUENCE_WITH_PADDINGS'])
tmp1 = validation_data[validation_data['Markername'] == 'DYS481']
validated_motifs = [extendMotif(x.split('_')[1]) for x in list(tmp1['Allelename_PH7'])]
for motif in validated_motifs:
    match_count = 0
    for s in sb_seq:
        if motif in s:
            match_count += 1
    print('%s out of %s sequences match the right allele (%s)' %(match_count, len(sb_seq), motif))

# PH7 - DXS10135
sb = data_rs[(data_rs['SAMPLE_NAME'] == 'PH7.ontarget.sup.hg38') & (data_rs['REGION'] == 'chrX:9338302-9338520')]
sb_seq = list(sb['SEQUENCE_WITH_PADDINGS'])
tmp1 = validation_data[validation_data['Markername'] == 'DXS10135']
validated_motifs = [extendMotif(x.split('_')[1]) for x in list(tmp1['Allelename_PH7'])]
for motif in validated_motifs:
    match_count = 0
    for s in sb_seq:
        if motif in s:
            match_count += 1
    print('%s out of %s sequences match the right allele (%s)' %(match_count, len(sb_seq), motif))

# PH7 - D13S317
sb = data_rs[(data_rs['SAMPLE_NAME'] == 'PH7.ontarget.sup.hg38') & (data_rs['REGION'] == 'chr13:82148023-82148100')]
sb['LEN_SEQUENCE_WITH_PADDINGS']
sb_seq = list(sb['SEQUENCE_WITH_PADDINGS'])
tmp1 = validation_data[validation_data['Markername'] == 'DYS481']
validated_motifs = [extendMotif(x.split('_')[1]) for x in list(tmp1['Allelename_PH7'])]
for motif in validated_motifs:
    match_count = 0
    for s in sb_seq:
        if motif in s:
            match_count += 1
    print('%s out of %s sequences match the right allele (%s)' %(match_count, len(sb_seq), motif))

# PH11 - DYS570
sb = data_rs[(data_rs['SAMPLE_NAME'] == 'PH11.ontarget.sup.hg38') & (data_rs['REGION'] == 'chrY:6993188-6993260')]
sb_seq = list(sb['SEQUENCE_WITH_PADDINGS'])
tmp1 = validation_data[validation_data['Markername'] == 'DYS570']
validated_motifs = [extendMotif(x.split('_')[1]) for x in list(tmp1['Allelename_PH11'])]
for motif in validated_motifs:
    match_count = 0
    for s in sb_seq:
        if motif in s:
            match_count += 1
    print('%s out of %s sequences match the right allele (%s)' %(match_count, len(sb_seq), motif))

# PH11 - DYS612
sb = data_rs[(data_rs['SAMPLE_NAME'] == 'PH11.ontarget.sup.hg38') & (data_rs['REGION'] == 'chrY:13640628-13640935')]
sb_seq = list(sb['SEQUENCE_WITH_PADDINGS'])
tmp1 = validation_data[validation_data['Markername'] == 'DYS612']
validated_motifs = [extendMotif(x.split('_')[1]) for x in list(tmp1['Allelename_PH11'])]
for motif in validated_motifs:
    match_count = 0
    for s in sb_seq:
        if motif in s:
            match_count += 1
    print('%s out of %s sequences match the right allele (%s)' %(match_count, len(sb_seq), motif))

# PH11 - DYS481
sb = data_rs[(data_rs['SAMPLE_NAME'] == 'PH11.ontarget.sup.hg38') & (data_rs['REGION'] == 'chrY:8558313-8558408')]
sb_seq = list(sb['SEQUENCE_WITH_PADDINGS'])
tmp1 = validation_data[validation_data['Markername'] == 'DYS481']
validated_motifs = [extendMotif(x.split('_')[1]) for x in list(tmp1['Allelename_PH11'])]
for motif in validated_motifs:
    match_count = 0
    for s in sb_seq:
        if motif in s:
            match_count += 1
    print('%s out of %s sequences match the right allele (%s)' %(match_count, len(sb_seq), motif))

# PH11 - DYS576
sb = data_rs[(data_rs['SAMPLE_NAME'] == 'PH11.ontarget.sup.hg38') & (data_rs['REGION'] == 'chrY:7185316-7185388')]
sb_seq = list(sb['SEQUENCE_WITH_PADDINGS'])
tmp1 = validation_data[validation_data['Markername'] == 'DYS576']
validated_motifs = [extendMotif(x.split('_')[1]) for x in list(tmp1['Allelename_PH11'])]
for motif in validated_motifs:
    match_count = 0
    for s in sb_seq:
        if motif in s:
            match_count += 1
    print('%s out of %s sequences match the right allele (%s)' %(match_count, len(sb_seq), motif))

# PH11 - FGA
sb = data_rs[(data_rs['SAMPLE_NAME'] == 'PH11.ontarget.sup.hg38') & (data_rs['REGION'] == 'chr4:154587714-154587823')]
sb_seq = list(sb['SEQUENCE_WITH_PADDINGS'])
tmp1 = validation_data[validation_data['Markername'] == 'FGA']
validated_motifs = [extendMotif(x.split('_')[1]) for x in list(tmp1['Allelename_PH11'])]
for motif in validated_motifs:
    match_count = 0
    for s in sb_seq:
        if motif in s:
            match_count += 1
    print('%s out of %s sequences match the right allele (%s)' %(match_count, len(sb_seq), motif))

# PH11 - D2S1338
sb = data_rs[(data_rs['SAMPLE_NAME'] == 'PH11.ontarget.sup.hg38') & (data_rs['REGION'] == 'chr2:218014852-218014954')]
sb_seq = list(sb['SEQUENCE_WITH_PADDINGS'])
tmp1 = validation_data[validation_data['Markername'] == 'D2S1338']
validated_motifs = [extendMotif(x.split('_')[1]) for x in list(tmp1['Allelename_PH11'])]
for motif in validated_motifs:
    match_count = 0
    for s in sb_seq:
        if motif in s:
            match_count += 1
    print('%s out of %s sequences match the right allele (%s)' %(match_count, len(sb_seq), motif))

# PH1 - DYS570
sb = data_rs[(data_rs['SAMPLE_NAME'] == 'PH1.ontarget.sup.hg38') & (data_rs['REGION'] == 'chrY:6993188-6993260')]
sb_seq = list(sb['SEQUENCE_WITH_PADDINGS'])
tmp1 = validation_data[validation_data['Markername'] == 'DYS570']
validated_motifs = [extendMotif(x.split('_')[1]) for x in list(tmp1['Allelename_PH1'])]
for motif in validated_motifs:
    match_count = 0
    for s in sb_seq:
        if motif in s:
            match_count += 1
    print('%s out of %s sequences match the right allele (%s)' %(match_count, len(sb_seq), motif))

# PH1 - DYS576
sb = data_rs[(data_rs['SAMPLE_NAME'] == 'PH1.ontarget.sup.hg38') & (data_rs['REGION'] == 'chrY:7185316-7185388')]
sb_seq = list(sb['SEQUENCE_WITH_PADDINGS'])
tmp1 = validation_data[validation_data['Markername'] == 'DYS576']
validated_motifs = [extendMotif(x.split('_')[1]) for x in list(tmp1['Allelename_PH1'])]
for motif in validated_motifs:
    match_count = 0
    for s in sb_seq:
        if motif in s:
            match_count += 1
    print('%s out of %s sequences match the right allele (%s)' %(match_count, len(sb_seq), motif))

# PH1 - FGA
sb = data_rs[(data_rs['SAMPLE_NAME'] == 'PH1.ontarget.sup.hg38') & (data_rs['REGION'] == 'chr4:154587714-154587823')]
sb_seq = list(sb['SEQUENCE_WITH_PADDINGS'])
tmp1 = validation_data[validation_data['Markername'] == 'FGA']
validated_motifs = [extendMotif(x.split('_')[1]) for x in list(tmp1['Allelename_PH1'])]
for motif in validated_motifs:
    match_count = 0
    for s in sb_seq:
        if motif in s:
            match_count += 1
    print('%s out of %s sequences match the right allele (%s)' %(match_count, len(sb_seq), motif))

# PH1 - PENTAD
sb = data_rs[(data_rs['SAMPLE_NAME'] == 'PH1.ontarget.sup.hg38') & (data_rs['REGION'] == 'chr21:43636190-43636272')]
sb_seq = list(sb['SEQUENCE_WITH_PADDINGS'])
tmp1 = validation_data[validation_data['Markername'] == 'PentaD']
validated_motifs = [extendMotif(x.split('_')[1]) for x in list(tmp1['Allelename_PH1'])]
for motif in validated_motifs:
    match_count = 0
    for s in sb_seq:
        if motif in s:
            match_count += 1
    print('%s out of %s sequences match the right allele (%s)' %(match_count, len(sb_seq), motif))

# PH3 - D2S1338
sb = data_rs[(data_rs['SAMPLE_NAME'] == 'PH3.ontarget.sup.hg38') & (data_rs['REGION'] == 'chr2:218014852-218014954')]
sb_seq = list(sb['SEQUENCE_WITH_PADDINGS'])
tmp1 = validation_data[validation_data['Markername'] == 'D2S1338']
validated_motifs = [extendMotif(x.split('_')[1]) for x in list(tmp1['Allelename_PH3'])]
for motif in validated_motifs:
    match_count = 0
    for s in sb_seq:
        try:
            if motif in s:
                match_count += 1
        except:
            pass
    print('%s out of %s sequences match the right allele (%s)' %(match_count, len(sb_seq), motif))

# PH15 - DYS570
sb = data_rs[(data_rs['SAMPLE_NAME'] == 'PH15.ontarget.sup.hg38') & (data_rs['REGION'] == 'chrY:6993188-6993260')]
sb_seq = list(sb['SEQUENCE_WITH_PADDINGS'])
tmp1 = validation_data[validation_data['Markername'] == 'DYS570']
validated_motifs = [extendMotif(x.split('_')[1]) for x in list(tmp1['Allelename_PH15'])]
for motif in validated_motifs:
    match_count = 0
    for s in sb_seq:
        try:
            if motif in s:
                match_count += 1
        except:
            pass
    print('%s out of %s sequences match the right allele (%s)' %(match_count, len(sb_seq), motif))

# PH15 - DYS612
sb = data_rs[(data_rs['SAMPLE_NAME'] == 'PH15.ontarget.sup.hg38') & (data_rs['REGION'] == 'chrY:13640628-13640935')]
sb_seq = list(sb['SEQUENCE_WITH_PADDINGS'])
tmp1 = validation_data[validation_data['Markername'] == 'DYS612']
validated_motifs = [extendMotif(x.split('_')[1]) for x in list(tmp1['Allelename_PH15'])]
for motif in validated_motifs:
    match_count = 0
    for s in sb_seq:
        try:
            if motif in s:
                match_count += 1
        except:
            pass
    print('%s out of %s sequences match the right allele (%s)' %(match_count, len(sb_seq), motif))

# PH15 - FGA
sb = data_rs[(data_rs['SAMPLE_NAME'] == 'PH15.ontarget.sup.hg38') & (data_rs['REGION'] == 'chr4:154587714-154587823')]
sb_seq = list(sb['SEQUENCE_WITH_PADDINGS'])
tmp1 = validation_data[validation_data['Markername'] == 'FGA']
validated_motifs = [extendMotif(x.split('_')[1]) for x in list(tmp1['Allelename_PH15'])]
for motif in validated_motifs:
    match_count = 0
    for s in sb_seq:
        try:
            if motif in s:
                match_count += 1
        except:
            pass
    print('%s out of %s sequences match the right allele (%s)' %(match_count, len(sb_seq), motif))

# PH15 - D18S51
sb = data_rs[(data_rs['SAMPLE_NAME'] == 'PH15.ontarget.sup.hg38') & (data_rs['REGION'] == 'chr18:63281665-63281773')]
sb_seq = list(sb['SEQUENCE_WITH_PADDINGS'])
tmp1 = validation_data[validation_data['Markername'] == 'D18S51']
validated_motifs = [extendMotif(x.split('_')[1]) for x in list(tmp1['Allelename_PH15'])]
for motif in validated_motifs:
    match_count = 0
    for s in sb_seq:
        try:
            if motif in s:
                match_count += 1
        except:
            pass
    print('%s out of %s sequences match the right allele (%s)' %(match_count, len(sb_seq), motif))

# PH15 - D3S1358
sb = data_rs[(data_rs['SAMPLE_NAME'] == 'PH15.ontarget.sup.hg38') & (data_rs['REGION'] == 'chr3:45540735-45540803')]
sb_seq = list(sb['SEQUENCE_WITH_PADDINGS'])
tmp1 = validation_data[validation_data['Markername'] == 'D3S1358']
validated_motifs = [extendMotif(x.split('_')[1]) for x in list(tmp1['Allelename_PH15'])]
for motif in validated_motifs:
    match_count = 0
    for s in sb_seq:
        try:
            if motif in s:
                match_count += 1
        except:
            pass
    print('%s out of %s sequences match the right allele (%s)' %(match_count, len(sb_seq), motif))

# now check why some regions are excluded
    # TPOX
        failed_dic[6]['PH6.ontarget.sup.hg38']
        data_motif_genomicRegions.loc[data_motif_genomicRegions['id'].str.contains('TPOX', case=False), 'locus']
        tmp = data_rs_info[(data_rs_info['SAMPLE_NAME'] == 'PH6.ontarget.sup.hg38') & (data_rs_info['REGION'] == 'chr2:1489649-1489684')]
        tmp = tmp.drop_duplicates(subset='READ_NAME').copy()
        tmp = tmp.dropna(subset=['HAPLOTAG'])
        general_QC_reads(tmp, cov, region, gender)
        # excluded because of QC -- low coverage
    # DXS10074
        data_motif_genomicRegions.loc[data_motif_genomicRegions['id'].str.contains('DXS10074', case=False), 'locus']
        tmp = data_rs_info[(data_rs_info['SAMPLE_NAME'] == 'PH6.ontarget.sup.hg38') & (data_rs_info['REGION'] == 'chrX:67757300-67757464')]
        tmp = tmp.drop_duplicates(subset='READ_NAME').copy()
        tmp = tmp.dropna(subset=['HAPLOTAG'])
        general_QC_reads(tmp, cov, region, gender)
        # excluded because of QC -- dropout after allele matching
    # DXS10103
        data_motif_genomicRegions.loc[data_motif_genomicRegions['id'].str.contains('DXS10103', case=False), 'locus']
        tmp = data_rs_info[(data_rs_info['SAMPLE_NAME'] == 'PH1.ontarget.sup.hg38') & (data_rs_info['REGION'] == 'chrX:134284908-134285040')]
        tmp = tmp.drop_duplicates(subset='READ_NAME').copy()
        tmp = tmp.dropna(subset=['HAPLOTAG'])
        general_QC_reads(tmp, cov, region, gender)
        # excluded because of QC -- low coverage
    # DYS385a-b
        data_motif_genomicRegions.loc[data_motif_genomicRegions['id'].str.contains('DYS385a-b', case=False), 'locus']
        tmp = data_rs_info[(data_rs_info['SAMPLE_NAME'] == 'PH6.ontarget.sup.hg38') & (data_rs_info['REGION'] == 'chrY:18639698-18639898')]
        tmp = tmp.drop_duplicates(subset='READ_NAME').copy()
        tmp = tmp.dropna(subset=['HAPLOTAG'])
        general_QC_reads(tmp, cov, region, gender)
        # excluded because of QC -- no match in database -- also is heterozygous
    # DYS389II
        data_motif_genomicRegions.loc[data_motif_genomicRegions['id'].str.contains('DYS389II', case=False), 'locus']
        tmp = data_rs_info[(data_rs_info['SAMPLE_NAME'] == 'PH6.ontarget.sup.hg38') & (data_rs_info['REGION'] == 'chrY:12500387-12500633')]
        tmp = tmp.drop_duplicates(subset='READ_NAME').copy()
        tmp = tmp.dropna(subset=['HAPLOTAG'])
        general_QC_reads(tmp, cov, region, gender)
        # excluded because of QC -- no match in database
    # DYS438
        data_motif_genomicRegions.loc[data_motif_genomicRegions['id'].str.contains('DYS438', case=False), 'locus']
        tmp = data_rs_info[(data_rs_info['SAMPLE_NAME'] == 'PH11.ontarget.sup.hg38') & (data_rs_info['REGION'] == 'chrY:12825887-12825949')]
        tmp = tmp.drop_duplicates(subset='READ_NAME').copy()
        tmp = tmp.dropna(subset=['HAPLOTAG'])
        # excluded because of QC -- low coverage
    # DYS461-460
        data_motif_genomicRegions.loc[data_motif_genomicRegions['id'].str.contains('DYS461-460', case=False), 'locus']
        tmp = data_rs_info[(data_rs_info['SAMPLE_NAME'] == 'PH6.ontarget.sup.hg38') & (data_rs_info['REGION'] == 'chrY:18888802-18889046')]
        tmp = tmp.drop_duplicates(subset='READ_NAME').copy()
        tmp = tmp.dropna(subset=['HAPLOTAG'])
        general_QC_reads(tmp, cov, 'chrY:18888802-18889046', gender)
        # excluded because of QC -- low coverage
    # Y-GATA-H4
        data_motif_genomicRegions.loc[data_motif_genomicRegions['id'].str.contains('Y-GATA-H4', case=False), 'locus']
        tmp = data_rs_info[(data_rs_info['SAMPLE_NAME'] == 'PH6.ontarget.sup.hg38') & (data_rs_info['REGION'] == 'chrY:16631210-16631912')]
        tmp = tmp.drop_duplicates(subset='READ_NAME').copy()
        tmp = tmp.dropna(subset=['HAPLOTAG'])
        general_QC_reads(tmp, cov, 'chrY:16631210-16631912', gender)
        # excluded because of QC -- low coverage

# Finally check the heterozygous chrY 
    # DYF387S1
        data_motif_genomicRegions.loc[data_motif_genomicRegions['id'].str.contains('DYF387S1', case=False), 'locus']
        # PH6
            tmp = data_rs_info[(data_rs_info['SAMPLE_NAME'] == 'PH6.ontarget.sup.hg38') & (data_rs_info['REGION'] == 'chrY:23785347-23785514')]
            tmp = tmp.drop_duplicates(subset='READ_NAME').copy()
            reads = list(tmp['SEQUENCE_WITH_PADDINGS'])
            # CE35_A[5]AAGA[2]AAGG[1]TAGG[1]AAGG[3]AAGA[2]AAGG[1]AAGA[2]AAGG[10]AAGA[13]A[3]T[1]A[8] -- CE40_A[5]AAGA[2]AAGG[1]TAGG[1]AAGG[3]AAGA[2]AAGG[1]AAGA[2]AAGG[10]AAGA[18]A[3]T[1]A[8]
            a1 = extendMotif('A[5]AAGA[2]AAGG[1]TAGG[1]AAGG[3]AAGA[2]AAGG[1]AAGA[2]AAGG[10]AAGA[13]A[3]T[1]A[8]')
            a2 = extendMotif('A[5]AAGA[2]AAGG[1]TAGG[1]AAGG[3]AAGA[2]AAGG[1]AAGA[2]AAGG[10]AAGA[18]A[3]T[1]A[8]')

        general_QC_reads(tmp, cov, 'chrY:16631210-16631912', gender)

# Check specifically SE33 -- it's not in the database, so we have to do with validation only
    # define samples to look at
    import os
    import re
    sample_list = ['PH1', 'PH3', 'PH6', 'PH7', 'PH11', 'PH15']
    validation = {'PH1' : ['CE23.2_CTTT[12]TTCT[12]CTTT[3]CT[2]TTCT[2]TTTT[1]CTTT[2]TTCT[1]TCCT[3]TTCT[1]CT[5]CTTT[3]CTAA[1]CT[3]TTGT[1]CT[2]TTCT[3]'], 'PH3' : ['CE28.2_CTTT[19]TTCT[10]CTTT[3]CT[2]TTCT[2]TTTT[1]CTTT[2]TTCT[1]TCCT[3]TTCT[1]CT[5]CTTT[3]CTAA[1]CT[3]TTGT[1]CT[2]TTCT[3]', 'CE27.2_CTTT[16]TTCT[12]CTTT[3]CT[2]TTCT[2]TTTT[1]CTTT[2]TTCT[1]TCCT[3]TTCT[1]CT[5]CTTT[3]CTAA[1]CT[3]TTGT[1]CT[2]TTCT[3]'], 'PH6' : ['CE16_CTTT[16]CT[2]TTCT[3]CTTT[3]TTCT[2]TT[1]TTCT[1]TCCT[3]TTCT[1]CT[5]CTTT[3]CTAA[1]CT[3]TTGT[1]CT[2]TTCT[3]_-4C>T', 'CE18_CTTT[18]CT[2]TTCT[3]CTTT[3]TTCT[2]TT[1]TTCT[1]TCCT[3]TTCT[1]CT[5]CTTT[3]CTAA[1]CT[3]TTGT[1]CT[2]TTCT[3]_-4C>T'], 'PH7' : ['CE22_CTTT[22]CT[2]TTCT[3]CTTT[3]TTCT[2]TT[1]TTCT[1]TCCT[3]TTCT[1]CT[5]CTTT[3]CTAA[1]CT[3]TTGT[1]CT[2]TTCT[3]_-4C>T', 'CE29.2_CTTT[20]TTCT[10]CTTT[3]CT[2]TTCT[2]TTTT[1]CTTT[2]TTCT[1]TCCT[3]TTCT[1]CT[5]CTTT[3]CTAA[1]CT[3]TTGT[1]CT[2]TTCT[3]'], 'PH11' : ['CE19_CTTT[19]CT[2]TTCT[3]CTTT[3]TTCT[2]TT[1]TTCT[1]TCCT[3]TTCT[1]CT[5]CTTT[3]CTAA[1]CT[3]TTGT[1]CT[2]TTCT[3]_-4C>T'], 'PH15' : ['CE26.2_CTTT[18]TTCT[9]CTTT[3]CT[2]TTCT[2]TTTT[1]CTTT[2]TTCT[1]TCCT[3]TTCT[1]CT[5]CTTT[3]CTAA[1]CT[3]TTGT[1]CT[2]TTCT[3]', 'CE27.2_CTTT[16]TTCT[12]CTTT[3]CT[2]TTCT[2]TTTT[1]CTTT[2]TTCT[1]TCCT[3]TTCT[1]CT[5]CTTT[3]CTAA[1]CT[3]TTGT[1]CT[2]TTCT[3]']}
    # define a function that given 2 sequences, find the maximum overlap between them
    def longest_common_substring(a, b):
        # Initialize a matrix to store the lengths of common substrings
        matrix = [[0] * (len(b) + 1) for _ in range(len(a) + 1)]
        # Variables to keep track of the longest substring and its end position
        max_length = 0
        end_position = 0
        # Fill in the matrix and find the longest common substring
        for i in range(1, len(a) + 1):
            for j in range(1, len(b) + 1):
                if a[i - 1] == b[j - 1]:
                    matrix[i][j] = matrix[i - 1][j - 1] + 1
                    if matrix[i][j] > max_length:
                        max_length = matrix[i][j]
                        end_position = i
        # Extract the longest common substring
        longest_substring = a[end_position - max_length:end_position]
        return longest_substring
    # dictionary for the sequences
    seqs_per_sample = {}
    validation_results = {}
    # iterate over samples
    for s in sample_list:
        print(s)
        # run treat to extract sequences -- bed file is larger to maximize probability of covering the sequence
        samtools_cmd = 'samtools view ../../../%s.ontarget.sup.hg38.bam chr6:88277144-88277400 | cut -f1,10' %(s)
        tmp_seqs = [x.rstrip() for x in os.popen(samtools_cmd)]
        # save the sequences
        seqs_per_sample[s] = tmp_seqs
        # extend the validated motif
        tmp_validation = validation[s]
        tmp_validated_alleles = [extendMotif(x.split('_')[1]) for x in tmp_validation]
        # then iterate over the sequences and try to match the validated motif
        tmp_validation_results = []
        for rname_seq in tmp_seqs:
            rname, seq = rname_seq.split('\t')
            for val in tmp_validated_alleles:
                if val in seq:
                    tmp_validation_results.append([rname, val, seq, 'full_match', 1, None, None, None, None])
                else:
                    tmp_common = longest_common_substring(seq, val)
                    # get indexes on the sequence
                    start_seq, end_seq = seq.find(tmp_common), seq.find(tmp_common) + len(tmp_common)
                    # get indexes on the validated allele
                    start_val, end_val = val.find(tmp_common), val.find(tmp_common) + len(tmp_common)
                    # then calculate the fraction of the validated motif that is covered
                    fraction_covered = len(tmp_common) / len(val)
                    # save
                    tmp_validation_results.append([rname, val, seq, 'partial_match', fraction_covered, start_seq, end_seq, start_val, end_val])
        validation_results[s] = tmp_validation_results
    # save table
    fout = open('se33_table_output.txt', 'w')
    fout.write('SAMPLE\tREAD_NAME\tVALIDATION\tSEQUENCE\tMATCH\tFRACTION_MATCH\tINDEX_START_SEQ\tINDEX_END_SEQ\tINDEX_START_VAL\tINDEX_END_VAL\n')
    for s in validation_results.keys():
        for x in validation_results[s]:
            fout.write('%s\t%s\n' %(s, '\t'.join([str(k) for k in x])))
    fout.close()
