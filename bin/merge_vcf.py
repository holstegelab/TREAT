# This script allows to merge VCF files of different individual samples together

# Libraries
print('* Loading libraries')
import pandas as pd
import os
import re
import sys
import time

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
# function to write header of vcf file
def writeVCFheader(output_fname):
    # open file
    outf = open(output_fname, 'w')
    # write header
    outf.write('##fileformat=VCFv4.2\n##INFO=<ID=REFERENCE_INFO,Number=2,Type=String,Description="Motif observed in the reference genome (GRCh38), and relative number of motif repetitions."\n##FORMAT=<ID=QC,Number=1,Type=String,Description="Quality summary of TREAT genotyping. PASS_BOTH: genotype agreed between reads-spanning and assembly. PASS_RSP: genotype from reads-spanning. PASS_ASM: genotype from assembly."\n##FORMAT=<ID=GT,Number=2,Type=String,Description="Phased size of the tandem repeats. H1_size | H2_size"\n##FORMAT=<ID=MOTIF,Number=2,Type=String,Description="Phased consensus motif found in the sample. H1_motif | H2_motif"\n##FORMAT=<ID=CN,Number=2,Type=String,Description="Phased number of copies of the motif found in the sample. H1_copies | H2_copies"\n##FORMAT=<ID=CN_REF,Number=2,Type=String,Description="Phased estimation of the reference motif as found in the sample. H1_motif_ref | H2_motif_ref"\n##FORMAT=<ID=DP,Number=1,Type=String,Description="Phased depth found in the sample. H1_depth | H2_depth"\n')
    outf.close()
    return

# function to write outputs
def writeOutputs(output_fname, vcf_data):
    # write header of vcf file
    writeVCFheader(output_fname)
    with open(output_fname, mode='a') as file:
        vcf_data.to_csv(file, header=True, index=False, sep='\t')
    # then compress it
    os.system('gzip %s' %(output_fname))
    return

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
first = pd.read_csv(vcf[0], sep="\t", skiprows=8)
# then iteratively read the others and merge them with the first
for f in vcf[1::]:
    tmp = pd.read_csv(f, sep="\t", skiprows=8)
    # merge common ids
    tmp_sb = tmp.iloc[:, [2] + list(range(9, len(tmp.columns)))]
    merged_df = first.merge(tmp_sb, on='ID', how='inner')
    # find common elements in "ID" column
    common_ids = first['ID'].isin(merged_df['ID'])
    common_ids_df2 = tmp['ID'].isin(merged_df['ID'])
    # find elements unique to df1
    unique_to_df1 = first[~common_ids].copy()
    # add NAs
    samples_df2 = list(tmp.columns[9::])
    for newsam in samples_df2:
        unique_to_df1[newsam] = 'NA;NA|NA;NA|NA;NA|NA;NA|NA;NA|NA'
    # the the other way around - find elements unique to df2
    unique_to_df2 = tmp[~common_ids_df2].copy()
    # add NAs
    samples_df1 = list(first.columns[9::])
    for newsam in samples_df1:
        unique_to_df2[newsam] = 'NA;NA|NA;NA|NA;NA|NA;NA|NA;NA|NA'
    # combine all together
    complete_df = pd.concat([merged_df, unique_to_df1, unique_to_df2], ignore_index=True)
    first = complete_df

# 3. write output file
output_fname = '%s/%s' %(outDir, outName)
output_fname = output_fname.replace('.gz', '')
writeOutputs(output_fname, first)

te_total = time.time()
time_total = te_total - ts_total
print('\n* VCF combined in %s seconds. Ciao!\t\t\t\t\t\t\t\t' %(round(time_total, 0)))

