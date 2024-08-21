# This script allows to merge VCF files of different individual samples together

# Libraries
print('* Loading libraries')
import pandas as pd
import os
import re
import gzip
import io
import sys
import time
import numpy as np

# Functions
### Functions to check input arguments
# Check directory
def checkOutDir(outDir):
    if outDir == './':
        return("** Output directory valid. Will output files in the current directory")
    else:
        # see whether directory exists
        if os.path.isdir(outDir) == True:
            return("** Output directory valid. Will output files %s" %(outDir))
        else:
            # create directory and put files there
            os.system('mkdir %s' %(outDir))
            return("** Output directory doesn't exist. Will create and will output files %s" %(outDir))

# Check if VCF files exist
def checkVCF(vcf):
    # then check them one by one if they exist
    valid_files = []
    non_valid_files = []
    for f in vcf:
        # check if file exist
        if os.path.exists(f) == True:
            # if file exists, check if it was produced by TREAT
            cmd_cat = 'cat' if '.gz' not in f else 'zcat'
            cmd = '%s %s | head | grep "TREAT"' %(cmd_cat, f)
            res = [x.rstrip() for x in os.popen(cmd)]
            if len(res) >0:
                valid_files.append(f)
            else:
                non_valid_files.append(f)
        else:
            non_valid_files.append(f)
    # summary
    if len(non_valid_files) == 0:
        return('** %s valid VCF provided.' %(len(vcf)))
    else:
        print('** %s non valid VCF provided: %s.\nExecution halted.' %(len(non_valid_files), ' '.join(non_valid_files)))
        sys.exit(1)  # Exit the script with a non-zero status code

### Functions to write output
# function to write outputs
def writeOutputs(output_fname, vcf_data):
    # write header of vcf file
    cmd = 'grep "##" %s > %s' %(vcf[0], output_fname) if '.gz' not in vcf[0] else 'zgrep "##" %s > %s' %(vcf[0], output_fname)
    os.system(cmd)
    # then write the dataframe
    with open(output_fname, mode='a') as file:
        vcf_data.to_csv(file, header=True, index=False, sep='\t')
    # then compress it
    os.system('gzip %s' %(output_fname))
    return '** Output produced --> %s' %(output_fname)

# Define a function to skip lines starting with '##'
def skip_comment_lines(file):
    if '.gz' in file:
        with gzip.open(file, 'rt') as f:
            lines = [line.rstrip().split('\t') for line in f if not line.startswith('##')]
    else:
        with open(file, 'r') as f:
            lines = [line for line in f if not line.startswith('##')]
    return pd.DataFrame(lines[1::], columns = lines[0])

# Combine 2 dataframes resembling the VCF
def combine_opt(df1, df2):
    # first merge df1 (entire) with df2 (sample-specific)
    df2_subset = df2.iloc[:, [2, 4] + list(range(9, df2.shape[1]))].copy()
    df2_subset.rename(columns={df2_subset.columns[1]: 'ALT_2'}, inplace=True)
    merged_df = pd.merge(df1, df2_subset, on='ID', how='outer')
    # now that it's merged, i need to fix the alleles
    for i in range(merged_df.shape[0]):
        percentage = round(i/merged_df.shape[0]*100, 2)
        print(f'{percentage:.2f}%', end='\r')
        # get alleles from df1 and df2
        df1_al = merged_df.iat[i, merged_df.columns.get_loc('ALT')]
        df2_al = merged_df.iat[i, merged_df.columns.get_loc('ALT_2')]
        # check nan
        if not pd.isna(df1_al) and not pd.isna(df2_al):
            df1_al = df1_al.split(',')
            df2_al = df2_al.split(',')
            # define an allele mapping dictionary
            al_map = {0 : 0}
            # iterate through df2 alleles
            for al in df2_al:
                if al in df1_al:
                    # if this is a known allele, get the index
                    al_index = df1_al.index(al)
                    al_map[df2_al.index(al) + 1] = al_index + 1
                else:
                    # otherwise, save the mapping and update the allele list
                    df1_al.append(al)
                    al_index = df1_al.index(al)
                    al_map[df2_al.index(al) + 1] = al_index + 1
            # update alleles
            merged_df.iat[i, 4] = ','.join(df1_al)
            # then need to update the genotypes
            for col in range(df1.shape[1] + 1, merged_df.shape[1]):
                merged_df.iat[i, col] = replace_field(arr = merged_df.iat[i, col], replace_dict = al_map)
        elif pd.isna(df1_al) and not pd.isna(df2_al):
            # here all samples from df1 should be NA;NA;NA;NA;NA;NA;NA
            for col in range(9, df1.shape[1]):
                merged_df.iat[i, col] = 'NA;NA;NA;NA;NA;NA;NA'
        elif not pd.isna(df1_al) and pd.isna(df2_al):
            # here all samples from df2 should be NA;NA;NA;NA;NA;NA;NA
            for col in range(df1.shape[1] + 1, merged_df.shape[1]):
                merged_df.iat[i, col] = 'NA;NA;NA;NA;NA;NA;NA'
    merged_df = merged_df.drop('ALT_2', axis = 1)
    return merged_df

# Function to replace array entry based on a dictionary
def replace_field(arr, replace_dict):
    # Split the string by semicolons
    fields = arr.split(';')
    # Split the second field by pipe and replace using the dictionary
    field_2 = fields[1].split('|')
    field_2 = [str(replace_dict[int(x)]) for x in field_2]
    # Join the modified field back with pipe
    fields[1] = '|'.join(field_2)
    # Join all fields back into a single string
    modified_string = ';'.join(fields)
    return modified_string

# Main
# Read arguments
vcf, outDir, outName = sys.argv[1::]
vcf = vcf.split(',')

# 1. Check arguments: VCFs and output directory
print('* Analysis started')
ts_total = time.time()
# 1.1 Check output directory
print(checkOutDir(outDir))
# 1.2 Check VCF files
print(checkVCF(vcf))

# 2. Read the first VCF
combined_df = skip_comment_lines(vcf[0])
# then iteratively read the others and merge them with the first
for f in vcf[1::]:
    print('** combining VCF %s/%s' %(vcf.index(f) + 1, len(vcf)))
    tmp = skip_comment_lines(f)
    # combine the 2 dfs
    combined_df = combine_opt(df1 = combined_df, df2 = tmp)

# 3. write output file
output_fname = '%s/%s' %(outDir, outName)
output_fname = output_fname.replace('.gz', '')
if '.vcf' not in output_fname:
    output_fname = output_fname + '.vcf'
writeOutputs(output_fname, combined_df)

te_total = time.time()
time_total = te_total - ts_total
print('\n* VCF combined in %s seconds. Ciao!\t\t\t\t\t\t\t\t' %(round(time_total, 0)))

