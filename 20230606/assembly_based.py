#!/usr/bin/env python3

# This script manages the read-based analysis

# Libraries
print('**** Loading libraries')
import sys
import os
import pysam
from functools import partial
import multiprocessing
import pandas as pd
import re
import time

# Functions
# Read bed file
def readBed(bed_dir):
    bed = {}
    count_reg = 0
    with open(bed_dir) as finp:
        for line in finp:
            if line.startswith('#'):
                pass
            else:
                line = line.rstrip().split()
                if len(line) >= 3:
                    chrom, start, end = line[0:3]
                    region_id = chrom + ':' + start + '-' + end
                    count_reg += 1
                    if 'chr' not in chrom:
                        chrom = 'chr' + str(chrom)
                    if chrom in bed.keys():
                        bed[chrom].append([start, end, region_id])
                    else:
                        bed[chrom] = [[start, end, region_id]]
    print('**** Found %s regions in %s chromosomes' %(count_reg, len(bed)))
    return bed, count_reg

# Check directory
def checkOutDir(out_dir):
    if out_dir[-1] == '/':
        out_dir = out_dir[:-1]
    if os.path.isdir(out_dir) == False:
        os.system('mkdir %s' %(out_dir))
        return("**** Output directory not found, will create.")
    else:
        return("**** Output directory found, will add outputs there.")

# Check bam file(s)
def checkBAM(bam_dir):
    if bam_dir[-1] == '/':
        bam_dir = bam_dir[:-1]
    if os.path.isdir(bam_dir) == True:              # in case a directory was submitted
        all_bams = [x.rstrip()for x in list(os.popen('ls %s/*bam' %(bam_dir)))]
        print("**** Found directory with %s bam" %(len(all_bams)))
    elif os.path.isfile(bam_dir) == True:           # in case is a single bam file
        print("**** Found single bam")
        all_bams = [bam_dir]
    elif ',' in bam_dir:                            # in case there is a comma-separated list of bams
        all_bams = bam_dir.split(',')
        print("**** Found %s bam" %(len(all_bams)))
    return all_bams

# Function to make assembly with otter and produce fasta files suitable for TRF
def assembly_otter(s, output_directory, ref_fasta, bed_file, number_threads):
    # define output name with the right directory
    outname = s.split('/')[-1].replace('.bam', '.fa')
    # run otter
    cmd = '/project/holstegelab/Software/nicco/bin/otter/build/otter assemble -b %s -r %s %s -t %s > %s/otter_local_asm/%s' %(bed_file, ref_fasta, s, number_threads, output_directory, outname)
    os.system(cmd)
    # collect otter sequences
    f_open = [x.rstrip() for x in open('%s/otter_local_asm/%s' %(output_directory, outname), 'r').readlines()]
    # define padding size and output results list
    padd_size = 50
    # also define the new fasta output containing only the repetitive sequence for TRF
    trf_input = open('%s/otter_local_asm/%s' %(output_directory, outname.replace('.fa', '_trf.fa')), 'w')
    tmp_res = []  
    prev_region = ''
    for x in f_open:
        if x.startswith('>'):
            region = x.split()[0].replace('>', '')
            read_id = x.replace(' ', '_') + '_1' if region != prev_region else x.replace(' ', '_') + '_2'
            prev_region = region
        else:
            seq = x; seq_len_with_paddings = len(seq); seq_no_paddings = x[(padd_size-1):len(x)-(padd_size+1)]; seq_len_no_paddings = len(seq_no_paddings)
            # save hit at this point
            tmp_res.append([outname.replace('.fa', ''), region, read_id.replace('>', ''), seq, seq_len_with_paddings, seq_no_paddings, seq_len_no_paddings])
            # then write the sequence without paddings
            trf_input.write('>%s;%s;%s\n' %(region, outname.replace('.fa', ''), read_id.replace('>', '')))
            trf_input.write('%s\n' %(seq_no_paddings))
    trf_input.close()
    return tmp_res

# Function to write fasta files for TRF
def writeFastaTRF(all_seqs, fasta_name):
    # define container for fasta outputs
    fasta_outputs = []
    # open file and write things
    with open(fasta_name, 'w') as outFile:
        for region in all_seqs:
            outFile.write('>%s;%s;%s\n%s\n' %(region[1], region[0], region[2], region[-4]))
    outFile.close()

# Run TRF given a sequence
def run_trf(index, all_fasta, distances, type):
    # then run tandem repeat finder
    cmd = 'trf4.10.0-rc.2.linux64.exe %s 2 7 7 80 10 50 200 -ngs -h' %(all_fasta[index])
    trf = [x for x in os.popen(cmd).read().split('\n') if x != '']
    # loop on trf results and save them into a list of lists
    x = 0; trf_matches = []
    sample_name = re.sub(r'^a[a-z]\.tmp_', '', os.path.basename(all_fasta[index])).replace('.fa', '')
    if type == 'otter':
        sample_name = sample_name.replace('_trf', '')
        if sample_name.split('_')[0] == 'reference':
            sample_name = 'reference'
    while x < len(trf):
        # check if the line is the header of an entry
        if trf[x].startswith('@'):
            # if so, save the corresponding information
            region, sample, read_id = trf[x].split(';')
            region = region.replace('@', '')
            x += 1
        while x < len(trf) and not trf[x].startswith('@'):
            tmp_trf_match = [read_id + '_' + region, 'NA'] + trf[x].split()
            trf_matches.append(tmp_trf_match)
            x += 1
    # finally create pandas df and assign column names
    if len(trf_matches) == 0:
        trf_matches = [['NA' for i in range(19)]] 
    df = pd.DataFrame(trf_matches)
    df.columns = ['ID', 'EXPECTED_MOTIF', 'START_TRF', 'END_TRF', 'LENGTH_MOTIF_TRF', 'COPIES_TRF', 'TRF_CONSENSUS_SIZE', 'TRF_PERC_MATCH', 'TRF_PERC_INDEL', 'TRF_SCORE', 'TRF_A_PERC', 'TRF_C_PERC', 'TRF_G_PERC', 'TRF_T_PERC', 'TRF_ENTROPY', 'TRF_MOTIF', 'TRF_REPEAT_SEQUENCE', 'TRF_PADDING_BEFORE', 'TRF_PADDING_AFTER']
    # finally, we need to add the reads where trf didn't find any motif
    if type == 'reads':
        df_seqs = pd.DataFrame(distances[index][0])
        df_seqs.columns = ['SAMPLE_NAME', 'REGION', 'READ_NAME', 'PASSES', 'READ_QUALITY', 'MAPPING_CONSENSUS', 'SEQUENCE_FOR_TRF', 'SEQUENCE_WITH_PADDINGS', 'LEN_SEQUENCE_FOR_TRF', 'LEN_SEQUENCE_WITH_PADDINGS']
    else:
        df_seqs = pd.DataFrame(distances[index])
        if sample_name == 'reference':
            df_seqs.columns = ['SAMPLE_NAME', 'REGION', 'READ_NAME', 'PASSES', 'READ_QUALITY', 'MAPPING_CONSENSUS', 'SEQUENCE_FOR_TRF', 'SEQUENCE_WITH_PADDINGS', 'LEN_SEQUENCE_FOR_TRF', 'LEN_SEQUENCE_WITH_PADDINGS']
            df_seqs['HAPLOTAG'] = 1
        else:
            df_seqs.columns = ['SAMPLE_NAME', 'REGION', 'READ_NAME', 'SEQUENCE_WITH_PADDINGS', 'LEN_SEQUENCE_WITH_PADDINGS', 'SEQUENCE_FOR_TRF', 'LEN_SEQUENCE_FOR_TRF']
            # add haplotag
            df_seqs['HAPLOTAG'] = df_seqs['READ_NAME'].str.split('_').str[-1]
            # add other columns and put NA
            df_seqs['PASSES'] = 'NA'; df_seqs['READ_QUALITY'] = 'NA'; df_seqs['MAPPING_CONSENSUS'] = 'NA'
    # add id
    df_seqs['ID'] = df_seqs['READ_NAME'].str.cat(df_seqs['REGION'], sep='_')
    # merge trf dataframe and reads dataframes
    complete_df = pd.merge(df_seqs, df, left_on = 'ID', right_on = 'ID', how = 'outer')
    return complete_df

# Otter pipeline
def otterPipeline(outDir, cpu, ref, bed_dir, inBam, count_reg):
    # create directory for outputs
    os.system('mkdir %s/otter_local_asm' %(outDir))
    # run local assembly in multiprocessing
    pool = multiprocessing.Pool(processes=cpu)
    otter_fun = partial(assembly_otter, output_directory = outDir, ref_fasta = ref, bed_file = bed_dir, number_threads = cpu)
    extract_results = pool.map(otter_fun, inBam)
    # do the same on the reference genome
    temp_beds = splitBed(bed_dir, cpu, outDir, count_reg)
    pool = multiprocessing.Pool(processes=cpu)
    extract_fun = partial(measureDistance_reference, window = window, ref = ref, output_directory = outDir)
    extract_results_ref = pool.map(extract_fun, temp_beds)
    all_fasta_ref = [outer_list[1] for outer_list in extract_results_ref]
    extract_results_ref = [outer_list[0] for outer_list in extract_results_ref]
    # combine reference and local assembly
    extract_results.extend(extract_results_ref)
    # TRF on local assemblies and reference
    trf_out_dir = '%s/trf_reads' %(outDir); os.system('mkdir %s' %(trf_out_dir));
    # list fasta files to do analysis on
    reads_fasta = [x.rstrip() for x in os.popen('ls %s/otter_local_asm/*_trf.fa' %(outDir))]
    # add the reference as well
    reads_fasta.extend(all_fasta_ref)
    # Run TRF in multiprocessing for each sample
    pool = multiprocessing.Pool(processes=cpu)
    trf_fun = partial(run_trf, all_fasta = reads_fasta, distances = extract_results, type = 'otter')
    index_fasta = [x for x in range(len(reads_fasta))]
    trf_results = pool.map(trf_fun, index_fasta)
    # Combine df from different samples together
    df_trf_combined = pd.concat(trf_results)
    return df_trf_combined

# Measure the distance in the reference genome
def measureDistance_reference(bed_file, window, ref, output_directory):    
    # sequence with paddings
    awk_command = """awk '{print $1":"$2-%s"-"$3+%s}' %s > %s_reformatted.txt""" %(window, window, bed_file, bed_file)
    os.system(awk_command)
    sequence_in_reference_with_padding = [x.rstrip() for x in list(os.popen('samtools faidx -r %s_reformatted.txt %s' %(bed_file, ref)))]        # sequence without padding
    # then store these results
    distances = []
    i = 0
    total_sequence = ''
    while i<len(sequence_in_reference_with_padding):
        x = sequence_in_reference_with_padding[i]
        if x.startswith('>'):
            if total_sequence != '':
                sequence_with_paddings, sequence = total_sequence, total_sequence[window:-window]
                total_sequence = ''
                distances.append(['reference', region, 'reference', 'NA', 'NA', 'NA', sequence, sequence_with_paddings, len(sequence), len(sequence_with_paddings)])
            chrom = x.replace('>', '').split(':')[0]
            start = int(x.replace('>', '').split(':')[1].split('-')[0])
            end = int(x.replace('>', '').split(':')[1].split('-')[1])
            region = chrom + ':' + str(start + window) + '-' + str(end - window)
            i += 1
        else:
            total_sequence = total_sequence + x
            i += 1
    # add last element
    sequence_with_paddings, sequence = total_sequence, total_sequence[window:-window]
    distances.append(['reference', region, 'reference', 'NA', 'NA', 'NA', sequence, sequence_with_paddings, len(sequence), len(sequence_with_paddings)])
    # then we write the fasta
    outfasta = '%s/reference_%s.fa' %(output_directory, bed_file.split('.')[-1])
    writeFastaTRF(distances, outfasta)
    return distances, outfasta

# Function to split bed file in n bed files
def splitBed(bed_dir, n, outDir, count_reg):
    number_lines_per_file = int(count_reg/n)
    cmd = 'split -l %s %s %s/tmp_bed.' %(number_lines_per_file, bed_dir, outDir)
    os.system(cmd)
    # then read the obtained temporary beds
    cmd = 'ls %s/tmp_bed.*' %(outDir)
    tmp_beds = [x.rstrip() for x in list(os.popen(cmd))]
    return tmp_beds

# Function to remove temporary files
def removeTemp(outDir):
    # list all files
    all_files = [x.rstrip() for x in list(os.popen('ls %s' %(outDir)))]
    all_files = ['%s/%s' %(outDir, x) for x in all_files if 'gz' not in x]
    all_files = [x for x in all_files if 'otter_local_asm' not in x]
    all_files = [x for x in all_files if 'trf_reads' not in x]
    # and remove them
    for x in all_files:
        os.remove(x)
    # then remove the folders
    os.system('rm -rf %s/otter_local_asm' %(outDir))
    os.system('rm -rf %s/trf_reads' %(outDir))

# Main
# Read arguments and make small changes
inBam_dir, bed_dir, outDir, ref, window, windowAss, cpu, ploidy, software, HaploDev, minimumSupport, minimumCoverage = sys.argv[1::]
window = int(window); cpu = int(cpu); ploidy = int(ploidy); windowAss = int(windowAss); minimumSupport = int(minimumSupport)

# 1. Check arguments: BED, output directory and BAMs
print('** Analysis started')
ts_total = time.time()
# 1.1 Check output directory
print(checkOutDir(outDir))
# 1.2 Read bed file
bed, count_reg = readBed(bed_dir)
# 1.3 Check BAM files
inBam = checkBAM(inBam_dir)

# 2. Check which software was selected and do things accordingly
if software == 'otter':
    # Run local assembly and TRF
    df_trf_phasing_combined = otterPipeline(outDir, cpu, ref, bed_dir, inBam, count_reg)
    # Save output
    outf = '%s/assembly_trf_phasing.txt.gz' %(outDir)
    df_trf_phasing_combined.to_csv(outf, sep = "\t", index=False, na_rep='NA', compression='gzip')
    print('**** Data combined and outputs are ready')
    # Haplotyping
    file_path = os.path.realpath(__file__)
    file_path = '/'.join(file_path.split('/')[:-1])
    os.system("python %s/call_haplotypes.py %s/assembly_trf_phasing.txt.gz %s %s %s %s %s" %(file_path, outDir, outDir, cpu, HaploDev, 'assembly', minimumSupport))
    # Remove temporary files
    removeTemp(outDir)
    te_total = time.time()
    time_total = te_total - ts_total
    print('\n** Analysis completed in %s seconds. Ciao!                   ' %(round(time_total, 0)))
else:
    print('!! Selected hifiasm but not implemented !! Quitting.')
