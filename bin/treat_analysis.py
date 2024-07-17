#!/usr/bin/env python

# TREAT analysis
# outlier-based analysis
# case-control analysis

# Libraries
import argparse
import sys
import os
import gzip
import numpy as np
from scipy.stats import chi2
from statsmodels.stats.multitest import fdrcorrection
from statsmodels.api import Logit
import pandas as pd
import warnings

###########################################################
# Define the parser
parser = argparse.ArgumentParser(description='TREAT: Tandem REpeat Haplotyping Toolkit ~ Analysis framework for outlier-based and case-control analyses')
###########################################################

###########################################################
# Define the arguments for analysis
# required arguments
parser.add_argument("-a", "--analysis", required = True, choices=["outlier", "case-control"], help = "Type of analysis to perform: [outlier / case-control]")
parser.add_argument("-v", "--vcf", default = 'None', help = "VCF file output of TREAT. Multiple files should be comma-separated.", required = True)
# optional arguments
# add arguments: --file is the case-control labels (mandatory if case-control is selected)
parser.add_argument("-l", "--labels", required=False, help = "Labels for case-control association. Only valid if case-control analysis is selected. This should be a tab-separated file with 2 columns and no header, reporting sample name and a binary label.", default='None')
# add arguments: --out is the output directory
parser.add_argument("-o", "--outDir", default = './', help = "Output directory where output will be placed. Default is the current directory.", required = False)
# add arguments: --outname is the name of the output file
parser.add_argument("-n", "--outName", default = 'treat_outlier_analysis.txt', help = "Name of the plot. Default name is treat_analysis_output.txt", required = False)
# add arguments: --region is the name of the region to analyze
parser.add_argument("-r", "--region", default = 'all', help = "Name of the region to analyze. By default, all regions will be analyzed.", required = False)
# add arguments: --madThr is the value to call outliers
parser.add_argument("-t", "--madThr", default = 3, help = "Median Absolute deviation value used to call outliers.", required = False)
# add arguments: --cpu is the number of cpus to use
parser.add_argument("-c", "--cpu", default = 1, help = "Number of CPUs to use (Default: 1)", required = False)
###########################################################

###########################################################
# Parse the arguments
args = parser.parse_args()
###########################################################

# Functions
# Function to check output directory
def checkOutDir(out_dir):
    if os.path.exists(out_dir):
        print('** Output directory exists: will output files there.')
        return out_dir
    elif os.path.exists(os.path.dirname(out_dir)):
        print('** Output directory does not exist. Will create and output there.')
        os.system('mkdir %s' %(out_dir))
        return out_dir
    else:
        print('!! Output directory does not exists and cannot be created. Please check. Exiting.')
        sys.exit(1)

# Function to check input vcf
def checkVCF(inp_vcf):
    if os.path.exists(inp_vcf):
        print('** Input VCF file exists.')
        return inp_vcf
    else:
        print('!! Input VCF does not exist. Please check. Exiting.')
        sys.exit(1)

# Function to check input regions
def checkRegion(region):
    if region == 'all':
        print('** All regions will be analyzed.')
    else:
        region = region.split(',')
        print('** %s regions will be analyzed.' %(len(region)))
    return region

# Function to check labels
def checkLabels(labels):
    if os.path.exists(labels):
        print('** Label file exists.')
        labels_dic = {}
        with open(labels) as finp:
            for line in finp:
                line = line.rstrip().split()
                sample, cc_status = line
                if cc_status in labels_dic.keys():
                    labels_dic[cc_status].append(sample)
                else:
                    labels_dic[cc_status] = [sample]
        if len(labels_dic.keys()) !=2:
            print('!! More than 2 unique labels in the case-control labels. Make sure there is no header and file is correctly formatted')
            sys.exit(1)
        else:
            label1_name = labels_dic.keys()
            print('** Found %s samples with label %s' %(len(labels_dic[list(labels_dic.keys())[0]]), list(labels_dic.keys())[0]))
            print('** Found %s samples with label %s' %(len(labels_dic[list(labels_dic.keys())[1]]), list(labels_dic.keys())[1]))
            return labels_dic
    else:
        print('!! Input VCF does not exist. Please check. Exiting.')
        sys.exit(1)

# Function for outlier analysis
def outlier_analysis(inp_vcf, region, madThr, cpu, out_dir, out_name):
    # initialize output file
    fout = open('%s/%s' %(out_dir, out_name), 'w')
    fout.write('REGION\tALLELE\tMEAN_ALL\tOUTLIER_SAMPLE\tOUTLIER_SIZE\tOUTLIER_DIST\tOUTLIER_P\n')
    # list for sample names
    sample_names = []
    # iterate over regions
    with gzip.open(inp_vcf) as finp:
        for line in finp:
            if line.startswith(b'##'):
                pass
            elif line.startswith(b'#CHROM'):
                line = line.rstrip().split()
                sample_names = [x.decode('utf-8') for x in line[9::]] + ['GRCh38']
            else:
                line = line.rstrip().split()
                chrom, pos, id, ref, format = [x.decode('utf-8') for x in [line[0], line[1], line[2], line[3], line[8]]]
                # check input file is from treat
                if format != 'QC;GT;GT_LEN;MOTIF;CN;CN_REF;DP':
                    print('!! VCF was likely not produced with TREAT. Please check.')
                    sys.exit(1)
                else:
                    # check if 1 region need to be processed or all of them
                    if (isinstance(region, list) == True and id in region) or (isinstance(region, list) == False):
                        # extract genotypes per sample
                        allele_sizes = [x.decode('utf-8').split(';')[2] for x in line[9::]] + [len(ref)]
                        # extract list of shorter alleles
                        short_alleles = [giveAllele(x, 'short') for x in allele_sizes]
                        # extract list of longer alleles
                        long_alleles = [giveAllele(x, 'long') for x in allele_sizes]
                        # extract list of sum of both alleles
                        join_alleles = [giveAllele(x, 'sum') for x in allele_sizes]
                        # find outliers for shorter alleles
                        res_short = scoreOutliers(short_alleles, madThr)
                        # find outliers for longer alleles
                        res_long = scoreOutliers(long_alleles, madThr)
                        # find outliers for sum of both alleles
                        res_join = scoreOutliers(join_alleles, madThr)
                        # write to file
                        writeToFile(res_short, fout, 'SHORT', id, sample_names)
                        writeToFile(res_long, fout, 'LONG', id, sample_names)
                        writeToFile(res_join, fout, 'JOIN', id, sample_names)
                    else:
                        pass
    fout.close()
    return '** Analysis completed. Cheers!'

# Function to write outliers to file
def writeToFile(res, fout, type, id, sample_names):
    # check if there are outliers, otherwise skip
    if len(res[1]) >0:
        for i in range(len(res[1])):
            fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(id, type, res[0], sample_names[res[1][i]], res[-1][i], res[-2][i], res[2][i]))
    else:
        pass

# Function for modified version of zscore
def scoreOutliers(alleles, madThr):
    try:
        # compute mean and covariance
        #mean_alleles = np.nanmean(alleles)
        #variance = np.nanvar(alleles, ddof=1)
        median_alleles = np.nanmedian(alleles)
        mad = np.nanmedian(np.abs(alleles - np.nanmedian(alleles))) + median_alleles*0.10
        # calculate the Mahalanobis distance
        #mahalanobis_dist = np.abs(alleles - mean_alleles) / np.sqrt(variance)
        mahalanobis_dist = np.abs(alleles - median_alleles) / (mad*madThr)
        # remove NaN values from mahalanobis_dist for further calculations
        valid_idx = ~np.isnan(mahalanobis_dist)
        valid_mahalanobis_dist = mahalanobis_dist[valid_idx]
        # compute p-values using chi-squared distribution
        p_values = 1 - chi2.cdf(mahalanobis_dist**2, df=1)
        # correct p-values using FDR method
        _, p_values_corrected_valid = fdrcorrection(p_values[valid_idx])
        # make sure to have all pvalues in the right place for every sample
        p_values_corrected = np.full(p_values.shape, np.nan)
        p_values_corrected[valid_idx] = p_values_corrected_valid
        # identify outliers
        outliers = np.where(p_values_corrected < 0.05)[0]
        outliers_p = p_values_corrected[outliers]
        outliers_dist = mahalanobis_dist[outliers]
        outliers_alleles = np.array(alleles)[outliers]
        return median_alleles, outliers, outliers_p, outliers_dist, outliers_alleles
    except:
        return 'NA', 'NA', 'NA', 'NA', 'NA' 

# Function that return the shorter or longer allele given a form of Allele1|Allele2
def giveAllele(allele, type):
    try:
        allele_list = [int(x) for x in str(allele).split('|')]
        if type == 'short':
            return min(allele_list)
        elif type == 'long':
            return max(allele_list)
        elif type == 'sum':
            return(sum(allele_list))
    except:
        return np.nan

# Function for case-control analysis
def casecontrol_analysis(inp_vcf, region, labels_dic, cpu, out_dir, out_name):
    # initialize output file
    fout = open('%s/%s' %(out_dir, out_name), 'w')
    fout.write('REGION\tBETA_SHORT\tSE_SHORT\tLOW_CI_SHORT\tUP_CI_SHORT\tP_SHORT\tBETA_LONG\tSE_LONG\tLOW_CI_LONG\tUP_CI_LONG\tP_LONG\tBETA_JOIN\tSE_JOIN\tLOW_CI_JOIN\tUP_CI_JOIN\tP_JOIN\n')
    fout.close()
    # list for sample names
    sample_names = []
    # iterate over regions
    with gzip.open(inp_vcf) as finp:
        for line in finp:
            if line.startswith(b'##'):
                pass
            elif line.startswith(b'#CHROM'):
                line = line.rstrip().split()
                sample_names = [x.decode('utf-8') for x in line[9::]]
                # check if we have the phenotype for these
                check_values_in_dict(sample_names, labels_dic)
            else:
                line = line.rstrip().split()
                chrom, pos, id, ref, format = [x.decode('utf-8') for x in [line[0], line[1], line[2], line[3], line[8]]]
                # check input file is from treat
                if format != 'QC;GT;GT_LEN;MOTIF;CN;CN_REF;DP':
                    print('!! VCF was likely not produced with TREAT. Please check.')
                    sys.exit(1)
                else:
                    # check if 1 region need to be processed or all of them
                    if (isinstance(region, list) == True and id in region) or (isinstance(region, list) == False):
                        # extract genotypes per sample
                        allele_sizes = [x.decode('utf-8').split(';')[2] for x in line[9::]] + [len(ref)]
                        # extract list of shorter alleles
                        short_alleles = [giveAllele(x, 'short') for x in allele_sizes]
                        # extract list of longer alleles
                        long_alleles = [giveAllele(x, 'long') for x in allele_sizes]
                        # extract list of sum of both alleles
                        join_alleles = [giveAllele(x, 'sum') for x in allele_sizes]
                        # combine in a dataframe
                        df_alleles = createDF(labels_dic, sample_names, short_alleles, long_alleles, join_alleles)
                        # run association for short
                        res_short = logit(df_alleles, 'Short')
                        # run association for long
                        res_long = logit(df_alleles, 'Long')
                        # run association for join
                        res_join = logit(df_alleles, 'Join')
                        # combine results
                        comb_res = [[id] + res_short + res_long + res_join]
                        comb_df = pd.DataFrame(comb_res, columns = ['REGION', 'BETA_SHORT', 'SE_SHORT', 'LOW_CI_SHORT', 'UP_CI_SHORT', 'P_SHORT', 'BETA_LONG', 'SE_LONG', 'LOW_CI_LONG', 'UP_CI_LONG', 'P_LONG', 'BETA_JOIN', 'SE_JOIN', 'LOW_CI_JOIN', 'UP_CI_JOIN', 'P_JOIN'])
                        # append to file
                        comb_df.to_csv('%s/%s' %(out_dir, out_name), mode = 'a', header = False, index = False, sep = "\t")
                    else:
                        pass
    return'** Analysis completed. Cheers!'

# Function to do logistic regression
def logit(df_alleles, type):
    try:
        # short allele
        model = Logit.from_formula('Label ~ %s' %(type), data = df_alleles)
        result = model.fit(disp=0)
        beta = result.params[1]
        se = result.bse[1]
        lowci = result.conf_int()[0][1]
        upci = result.conf_int()[1][1]
        pval = result.pvalues[1]
        return [beta, se, lowci, upci, pval]
    except:
        return ['NA', 'NA', 'NA', 'NA', 'NA']

# Function to create dataframe with values of interest for the logistic regression
def createDF(labels_dic, sample_names, short_alleles, long_alleles, join_alleles):
    temp = []
    for label, samples in labels_dic.items():
        for sample in samples:
            if sample in sample_names:
                temp_short = short_alleles[sample_names.index(sample)]
                temp_long = long_alleles[sample_names.index(sample)]
                temp_join = join_alleles[sample_names.index(sample)]
                temp.append({'Sample' : sample, 'Label' : int(label), 'Short' : temp_short, 'Long' : temp_long, 'Join' : temp_join})
    df = pd.DataFrame(temp)
    return(df)

# Function to check overlap between a list and a dictionary value
def check_values_in_dict(sample_names, labels_dic):
    # Extract dictionary values
    labels_samples = list(labels_dic.values())[0] + list(labels_dic.values())[1]
    # Check each value in the list against dictionary values
    if all(sample_name in labels_samples for sample_name in sample_names):
        return '** All samples in VCF have a phenotype.'
    else:
        return '** Not all samples in VCF have a phenotype. Those without phenotype will be excluded.'

# Main
###########################################################
# Parse input arguments and set up for running
inp_vcf, analysis, labels, out_dir, out_name, region, madThr, cpu = args.vcf, args.analysis, args.labels, args.outDir, args.outName, args.region, int(args.madThr), args.cpu
if analysis == 'case-control' and labels == 'None':
    print('!! Case-control analysis was chosen, but no case-control labels were given. Exiting.\n')
    sys.exit(1)
elif analysis == 'case-control' and labels != 'None':
    print('** Outlier analysis was chosen')
    # check output directory
    out_dir = checkOutDir(out_dir)
    # check input vcf
    inp_vcf = checkVCF(inp_vcf)
    # check regions
    regions = checkRegion(region)
    # check labels
    labels_dic = checkLabels(labels)
    # check name
    if out_name == 'treat_analysis_output.txt':
        out_name = 'treat_casecontrol_analysis.txt'
    # case-control analysis
    casecontrol_analysis(inp_vcf, regions, labels_dic, cpu, out_dir, out_name)
elif analysis == 'outlier':
    print('** Outlier analysis was chosen')
    # check output directory
    out_dir = checkOutDir(out_dir)
    # check input vcf
    inp_vcf = checkVCF(inp_vcf)
    # check regions
    regions = checkRegion(region)
    # do outlier analysis
    outlier_analysis(inp_vcf, regions, madThr, cpu, out_dir, out_name)
else:
    print('!! Something went wrong with your analysis. Check again your input parameters. Exiting.')
    sys.exit(1)
###########################################################
