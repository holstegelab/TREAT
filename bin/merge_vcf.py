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
def combine_2_df(df1, df2):
    common_df = pd.DataFrame()
    # iterates over the rows of dataframe 2 (not the merged)
    for index, row in df2.iterrows():
        # take regions
        region2 = row['ID']
        # create dictionary with alleles mapping
        al_map = {}
        all_alt = []
        # check if region is in merged df1
        if region2 in list(df1['ID']):
            temp_df = df1[df1['ID'] == region2].copy()
            # then extract the alleles from both
            alt2 = row['ALT'].split(',')
            all_alt = [x.split(',') for x in temp_df['ALT']][0]
            # check if alleles are the same, otherwise add them to the merged df
            for al in alt2:
                if al in all_alt:
                    al_index = all_alt.index(al) + 1
                    al_map[alt2.index(al) + 1] = al_index
                else:
                    all_alt.append(al)
                    al_index = len(all_alt)
                    al_map[alt2.index(al) + 1] = al_index
            # after checking the alleles, we need to reformat the samples
            temp_df['ALT'] = ','.join(all_alt)
            for col in df2.iloc[:, 9:]:
                tmp_value = row[col]
                # update the genotype based on the mapping dictionary
                tmp_mod = np.array([replace_field(arr = tmp_value, replace_dict = al_map)])
                temp_df[col] = tmp_mod
            # finally add to the common df
            common_df = pd.concat([common_df, temp_df])
        else:
            temp_df = pd.DataFrame([row.values[0:9]], columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'])
            # add all samples from df1
            for col in df1.iloc[:, 9:]:
                temp_df[col] = 'NA;NA;NA;NA;NA;NA;NA'
            # then add all samples from df2
            for col in range(9, len(row)):
                temp_df[df2.columns[col]] = row[col]
            # finally add to the common df
            common_df = pd.concat([common_df, temp_df])
    # then we need to add the rows of dataframe 1 not in 2
    for index, row in df1.iterrows():
        region1 = row['ID']
        if region1 not in list(df2['ID']):
            temp_df = pd.DataFrame([row.values], columns = df1.columns)
            # then add all samples from df2
            for col in df2.iloc[:, 9:]:
                temp_df[col] = 'NA;NA;NA;NA;NA;NA;NA'
            # finally add to the common df
            common_df = pd.concat([common_df, temp_df])
    return common_df

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
    combined_df = combine_2_df(df1 = combined_df, df2 = tmp)

# 3. write output file
output_fname = '%s/%s' %(outDir, outName)
output_fname = output_fname.replace('.gz', '')
if '.vcf' not in output_fname:
    output_fname = output_fname + '.vcf'
writeOutputs(output_fname, combined_df)

te_total = time.time()
time_total = te_total - ts_total
print('\n* VCF combined in %s seconds. Ciao!\t\t\t\t\t\t\t\t' %(round(time_total, 0)))

