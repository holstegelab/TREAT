# LIBRARIES
import sys
import os
import random
import pysam
from functools import partial
import multiprocessing
import pandas as pd
import re
import math
import time
from Bio.Seq import reverse_complement
import numpy as np
import itertools
import scipy.stats as stats
import statistics
from sklearn.cluster import KMeans
import shutil
import warnings
import gzip

### FUNCTIONS TO CHECK DIRECTORIES AND FILES
# Read bed file
def readBed(bed_dir, out_dir):
    # in any case, we should re-write the bed file as there are things to check: tab separated, and intervals should be of at least 1 bp
    fout_name = '%s/correct_bed.bed' %(out_dir)
    fout = open(fout_name, 'w')
    # counter for valid and modified intervals
    counter_valid = 0
    counter_invalid = 0
    bed = {}
    count_reg = 0
    try:
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
                        # add chromosome label if not there
                        if 'chr' not in chrom:
                            chrom = 'chr' + str(chrom)
                        # check size of interval, if 0 or negative, put 1
                        if int(end) - int(start) <= 0:
                            counter_invalid += 1
                        else:
                            counter_valid += 1
                            # add region to the new bed file
                            fout.write('%s\t%s\t%s\n' %(chrom, start, end))
                            # add elements to the global dictionary
                            if chrom in bed.keys():
                                bed[chrom].append([start, end, region_id])
                            else:
                                bed[chrom] = [[start, end, region_id]]
        fout.close()
        print('** BED file: found %s regions in %s chromosomes' %(count_reg, len(bed)))
        if counter_invalid >0:
            print('** Of these, %s are valid intervals, and %s are invalid. Invalid intervals (distance between start and end of 0) have been removed.' %(counter_valid, counter_invalid))
        else:
            print('** All intervals are valid.')
        return bed, count_reg, fout_name
    except:
        print("\n!!! Input BED file is missing or wrongly formatted. Make sure to provide a genuine tab-separated BED file without header.\nExecution halted.")
        sys.exit(1)  # Exit the script with a non-zero status code

# Check directory
def checkOutDir(out_dir):
    if out_dir[-1] == '/':
        out_dir = out_dir[:-1]
    if os.path.isdir(out_dir) == False:
        os.system('mkdir %s' %(out_dir))
        return("** Output directory valid.")
    else:
        # check if directory is empty or not
        if any(os.scandir(out_dir)):
            print("\n!!! Output directory you indicated is not empty. This may cause issues and unexpected results. Please indicate a new folder instead.\nExecution halted.")
            sys.exit(1)  # Exit the script with a non-zero status code
        else:
            return("** Output directory exists but empty. Will add outputs there.")

# Check bam file(s)
def checkBAM(bam_dir):
    try:
        if bam_dir[-1] == '/':
            bam_dir = bam_dir[:-1]
        if os.path.isdir(bam_dir) == True:              # in case a directory was submitted
            all_bams = [x.rstrip()for x in list(os.popen('ls %s/*bam' %(bam_dir)))]
            print("** BAM file(s): found directory with %s bam" %(len(all_bams)))
        elif os.path.isfile(bam_dir) == True:           # in case is a single bam file
            print("** BAM file(s): found single bam")
            all_bams = [bam_dir]
        elif ',' in bam_dir:                            # in case there is a comma-separated list of bams
            all_bams = bam_dir.split(',')
            print("** BAM file(s): found %s bam" %(len(all_bams)))
        return all_bams
    except:
        print("\n!!! Input BAM file is missing or wrongly formatted. Make sure to provide a genuine BAM file.\nExecution halted.")
        sys.exit(1)  # Exit the script with a non-zero status code

# Function to create Log file -- Reads analysis
def createLogReads(inBam, bed_dir, outDir, ref, window, cpu, phasingData, mappingSNP, HaploDev, minimumSupport, minimumCoverage):
    foutname = open('%s/treat_run.log' %(outDir), 'w')
    foutname.write('Read-based analysis selected\n')
    foutname.write('** Required argument:\n')
    foutname.write("\tInput BAM file(s): %s\n" %(inBam))
    foutname.write("\tInput BED file: %s\n" %(bed_dir))
    foutname.write("\tOutput directory: %s\n" %(outDir))
    foutname.write("\tReference genome: %s\n" %(ref))
    foutname.write('** Optional arguments:\n')
    foutname.write("\tWindow: %s\n" %(window))
    foutname.write("\tNumber of threads: %s\n" %(cpu))
    foutname.write("\tPhasing data: %s\n" %(phasingData))
    foutname.write("\tPhasing data IDs: %s\n" %(mappingSNP))
    foutname.write("\tHaplotyping deviation: %s\n" %(HaploDev))
    foutname.write("\tMinimum supporting reads: %s\n" %(minimumSupport))
    foutname.write("\tMinimum coverage: %s\n" %(minimumCoverage))
    foutname.write("\n")
    foutname.write('Effective command line:\nTREAT.py reads -i %s -b %s -o %s -r %s -w %s -t %s -p %s -m %s -d %s -minSup %s -minCov %s\n' %(inBam, bed_dir, outDir, ref, window, cpu, phasingData, mappingSNP, HaploDev, minimumSupport, minimumCoverage))
    foutname.close()
    print('** Log file written to %s/treat_run.log' %(outDir))
    return foutname

# Function to create Log file -- Assembly analysis
def createLogAsm(inBam, bed_dir, outDir, ref, window, cpu, windowAss, ploidy, software, HaploDev, minimumSupport, minimumCoverage):
    foutname = open('%s/treat_run.log' %(outDir), 'w')
    foutname.write('Assembly-based analysis selected\n')
    foutname.write('** Required argument:\n')
    foutname.write("\tInput BAM file(s): %s\n" %(inBam))
    foutname.write("\tInput BED file: %s\n" %(bed_dir))
    foutname.write("\tOutput directory: %s\n" %(outDir))
    foutname.write("\tReference genome: %s\n" %(ref))
    foutname.write('** Optional arguments:\n')
    foutname.write("\tWindow: %s\n" %(window))
    foutname.write("\tWindow for assembly: %s\n" %(windowAss))
    foutname.write("\tNumber of threads: %s\n" %(cpu))
    foutname.write("\tPloidy: %s\n" %(ploidy))
    foutname.write("\tAssembler: %s\n" %(software))
    foutname.write("\tHaplotyping deviation: %s\n" %(HaploDev))
    foutname.write("\tMinimum supporting reads: %s\n" %(minimumSupport))
    foutname.write("\tMinimum coverage: %s\n" %(minimumCoverage))
    foutname.write("\n")
    foutname.write('Effective command line:\nTREAT.py assembly -i %s -b %s -o %s -r %s -w %s -wAss %s -t %s -s %s -p %s -d %s -minSup %s -minCov %s\n' %(inBam, bed_dir, outDir, ref, window, windowAss, cpu, software, ploidy, HaploDev, minimumSupport, minimumCoverage))
    foutname.close()
    print('** Log file written to %s/treat_run.log' %(outDir))
    return foutname

### FUNCTIONS FOR OTTER-BASED ASSEMBLY
# Function to make assembly with otter and produce fasta files suitable for TRF
def assembly_otter(s, output_directory, ref_fasta, bed_file, number_threads, windowAss):
    # define output name with the right directory
    outname = s.split('/')[-1].replace('.bam', '.fa')
    # run otter
    cmd = 'otter assemble -b %s -r %s %s -t %s -o %s > %s/otter_local_asm/%s' %(bed_file, ref_fasta, s, number_threads, windowAss, output_directory, outname)
    os.system(cmd)
    # collect otter sequences
    f_open = [x.rstrip() for x in open('%s/otter_local_asm/%s' %(output_directory, outname), 'r').readlines()]
    # define padding size and output results list
    padd_size = windowAss
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
def run_trf_otter(index, all_fasta, distances, type):
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
def otterPipeline(outDir, cpu, ref, bed_dir, inBam, count_reg, windowAss, window):
    print('** Assembler: otter')
    # create directory for outputs
    os.system('mkdir %s/otter_local_asm' %(outDir))
    # run local assembly in multiprocessing
    pool = multiprocessing.Pool(processes=cpu)
    otter_fun = partial(assembly_otter, output_directory = outDir, ref_fasta = ref, bed_file = bed_dir, number_threads = cpu, windowAss = windowAss)
    extract_results = pool.map(otter_fun, inBam)
    pool.close()
    # do the same on the reference genome
    temp_beds = splitBed(bed_dir, cpu, outDir, count_reg)
    pool = multiprocessing.Pool(processes=cpu)
    extract_fun = partial(measureDistance_reference, window = window, ref = ref, output_directory = outDir)
    extract_results_ref = pool.map(extract_fun, temp_beds)
    pool.close()
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
    trf_fun = partial(run_trf_otter, all_fasta = reads_fasta, distances = extract_results, type = 'otter')
    index_fasta = [x for x in range(len(reads_fasta))]
    trf_results = pool.map(trf_fun, index_fasta)
    pool.close()
    #xx = run_trf(index = 0, all_fasta = reads_fasta, distances = extract_results, type = 'otter')
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
    number_lines_per_file = math.ceil(count_reg/n)
    cmd = 'split -l %s %s %s/tmp_bed.' %(number_lines_per_file, bed_dir, outDir)
    os.system(cmd)
    # then read the obtained temporary beds
    cmd = 'ls %s/tmp_bed.*' %(outDir)
    tmp_beds = [x.rstrip() for x in list(os.popen(cmd))]
    return tmp_beds

### FUNCTIONS FOR HIFIASM-BASED ASSEMBLY
def hifiasmPipeline(outDir, cpu, ref, bed_dir, inBam, count_reg, ploidy):
    print('** Assembler: hifiasm')
    # extract reads for assembly
    print("** Extract reads")
    temp_fasta, temp_beds = extractRead(inBam, bed_dir, outDir, cpu, count_reg, 'fasta')
    # do assembly and conversion of gfa to fasta
    print("** Start with assembly: 1 sample at the time")
    print("** Parameters in use for hifiasm: -n 2 -r 3 -N 200")
    assembly_dir = hifiasmAssembly(outDir, cpu, ref, temp_fasta, inBam, ploidy)
    # clean contigs
    assembly_dir_clean = cleanContigs(assembly_dir, cpu)
    # align contigs
    print("** Start with alignment: 1 sample at the time")
    print("** Parameters in use for minimap2: -aYx asm10")
    assembly_aligned = alignContigs(assembly_dir_clean, cpu, ref)
    # then we re-extract sequences of interest like in the read-based method
    print("** Extract sequences of interest")
    temp_bams, temp_beds = extractRead(assembly_aligned, temp_beds, outDir, cpu, count_reg, 'bam')
    # get sequences
    pool = multiprocessing.Pool(processes=cpu)
    extract_fun = partial(distributeExtraction, bed = bed, window = window)
    extract_results = pool.map(extract_fun, temp_bams)
    pool.close()
    print('** Exact SV intervals extracted')
    all_fasta = [outer_list[1] for outer_list in extract_results]
    # then do the same on the reference genome
    pool = multiprocessing.Pool(processes=cpu)
    extract_fun = partial(measureDistance_reference, window = window, ref = ref, output_directory = outDir)
    extract_results_ref = pool.map(extract_fun, temp_beds)
    pool.close()
    all_fasta_ref = [outer_list[1] for outer_list in extract_results_ref]
    print('** Exact SV intervals from reference extracted')
    # combine reference with other samples
    extract_results.extend(extract_results_ref)
    all_fasta.extend(all_fasta_ref)
    # trf
    # run TRF in multiprocessing for each sample
    print("** TRF running")
    pool = multiprocessing.Pool(processes=cpu)
    trf_fun = partial(run_trf, all_fasta = all_fasta, distances = extract_results, type = 'reads')
    index_fasta = [x for x in range(len(all_fasta))]
    trf_results = pool.map(trf_fun, index_fasta)
    pool.close()
    # combine df from different samples together
    df_trf_combined = pd.concat(trf_results)
    print('** TRF done on all reads and samples')
    # assign haplotags
    print("** Assigning haplotag based on assembly")
    df_trf_phasing_combined = assignHaplotags(df_trf_combined)
    return df_trf_phasing_combined

# Function to assign haplotags
def assignHaplotags(df_trf_combined):
    # find samples and regions
    all_samples = list(set(list(df_trf_combined['SAMPLE_NAME'])))
    all_regions = list(set(list(df_trf_combined['REGION'])))
    newdf = pd.DataFrame()
    # iterate over samples and regions
    for s in all_samples:
        for r in all_regions:
            # get data
            tmp = df_trf_combined[(df_trf_combined['SAMPLE_NAME'] == s) & (df_trf_combined['REGION'] == r)].copy()
            # Apply the function to each row to set the HAPLOTAG value
            tmp['HAPLOTAG'] = tmp.apply(set_hap_value, axis=1) if len(tmp) >0 else None
            # then add to new dataframe
            newdf = newdf.append(tmp)
    return newdf

# Set the HAPLOTAG value based on the contig name
def set_hap_value(row):
    if row['SAMPLE_NAME'] == 'reference':
        return 1
    elif 'h1' in row['READ_NAME']:
        return 1
    elif 'h2' in row['READ_NAME']:
        return 2
    else:
        return None

# Function to execute the sequence extraction in multiple processors
def distributeExtraction(x, bed, window):
    # container for results
    tmp_results = []
    # get sample name
    sample_name = re.sub(r'^a[a-z]\.tmp_', '', os.path.basename(x)).replace('.bam', '')
    # loop over reads with pysam
    with pysam.AlignmentFile(x, 'rb', check_sq=False) as bamfile:
        for read in bamfile:
            # extract read information
            ref_chrom, ref_start, ref_end, query_name, query_sequence, cigartuples, tags, is_secondary, is_supplementary, cigarstring = read.reference_name, int(read.reference_start), int(read.reference_end), read.query_name, read.query_sequence, read.cigartuples, read.tags, read.is_supplementary, read.is_secondary, read.cigarstring
            # check how many regions we overlap with this read
            regions_overlapping = checkIntervals(bed, ref_chrom, ref_start, ref_end, window)
            # then get the sequence of the read in the interval
            regions_overlapping_info = getSequenceInterval(regions_overlapping, tags, is_secondary, is_supplementary, query_name, query_sequence, window, ref_start, ref_end, cigartuples, sample_name)       
            # add to results
            for lst in regions_overlapping_info:
                tmp_results.append(lst)
    # name of the fasta output
    fasta_name = x.replace('.bam', '.fa')
    # finally write fasta files
    writeFastaTRF(tmp_results, fasta_name)
    return tmp_results, fasta_name

# Check how many intervals a sequence is included in
def checkIntervals(bed, chrom, start, end, window):
    # define output sublist
    sublist = []
    # iterate over bed file
    for interval in bed[chrom]:
        if (int(interval[0]) - window) >= start and (int(interval[1]) + window) <= end:
            sublist.append(interval[-1])
    return sublist

# Function to parse the CIGAR string
def findPositionOfInterestWhile(cigar, region_start, region_end, ref_start, ref_end, window):
    # define start and ending positions of interest, with and without padding
    positions_of_interest_with_padding = (region_start - window) - ref_start
    positions_of_interest_with_padding_end = (region_end + window) - ref_start
    positions_of_interest = region_start - ref_start
    positions_of_interest_end = region_end - ref_start
    # define counter of the reference and the raw sequences
    counter_ref = 0; counter_raw = 0; counter_raw_padd = 0
    # define positions of interest
    pos_interest = 0; pos_interest_padd = 0; pos_interest_end = 0; pos_interest_padd_end = 0
    # make a list of 1 cigar element per position
    cigar_per_base = [x for cse in cigar for x in [cse[0]] * cse[1]]
    # Then loop on this list
    i = 0; run = True
    while (run == True) and (i < len(cigar_per_base)):
        x = cigar_per_base[i]
        # Parse cigar types: 
        if (x == 7) or (x == 0):   # 7 --> =
            counter_raw += 1; counter_ref += 1; counter_raw_padd += 1
        elif x == 8:   # 8 --> X
            counter_raw += 1; counter_ref += 1; counter_raw_padd += 1
        elif x == 1:   # 1 --> I
            counter_raw +=1; counter_raw_padd += 1
        elif x == 2:   # 2 --> D
            counter_ref += 1
        elif x == 4:  # 4 --> S
            counter_raw += 1; counter_raw_padd += 1
        elif x == 5:  # 5 --> H
            counter_raw += 1; counter_raw_padd += 1
            print("!!! Alignments are hard clipped. Impossible to take actual sequence!")
        else:
            print("!!! Unknown term in cigar string --> %s" % (x))
            break
        # Then check if we reached the start/end position without padding
        if pos_interest == 0 and counter_ref == (positions_of_interest-1):
            pos_interest = counter_raw
        if pos_interest_end == 0 and counter_ref == (positions_of_interest_end-1):
            pos_interest_end = counter_raw
        # Then check if we reached the start/end position with padding
        if pos_interest_padd == 0 and counter_ref == (positions_of_interest_with_padding-1):
            pos_interest_padd = counter_raw_padd
        if pos_interest_padd_end == 0 and counter_ref == (positions_of_interest_with_padding_end-1):
            pos_interest_padd_end = counter_raw_padd
        # Finally check if we need to loop again
        if 0 in [pos_interest, pos_interest_end, pos_interest_padd, pos_interest_padd_end]:
            run = True; i += 1
        else:
            run = False
    return pos_interest, pos_interest_end, pos_interest_padd, pos_interest_padd_end

# Extract sequence given interval and read looking at CIGAR
def getSequenceInterval(regions_overlapping, tags, is_secondary, is_supplementary, query_name, query_sequence, window, ref_start, ref_end, cigartuples, sample_name):
    # define container for the information
    info_reads = []
    # extract tags from read
    info = tags; np, rq, mc = 'NA', 'NA', 'NA'
    for x in info:
        if x[0] == "np":
            np = "NP:%s" %(x[1])
        elif x[0] == "rq":
            rq = "RQ:%s" %(x[1])
        elif x[0] == "mc":
            mc = "MC:%s" %(x[1])
    # iterate over the regions encompassed by the read
    for region in regions_overlapping:
        # extract region stats
        chrom, interval = region.split(':')
        start, end = [int(x) for x in interval.split('-')]
        # exclude secondary and supplementary alignments
        if not is_secondary and not is_supplementary:
            # extract read name
            read_name = query_name
            # look into CIGAR to find positions
            pos_interest, pos_interest_end, pos_interest_padd, pos_interest_padd_end = findPositionOfInterestWhile(cigartuples, start, end, ref_start, ref_end, window)
            # then extract sequence
            sequence_interest = str(query_sequence)[pos_interest : pos_interest_end]
            sequence_interest_len = len(sequence_interest)
            sequence_interest_with_padding = str(query_sequence)[pos_interest_padd : pos_interest_padd_end]
            sequence_interest_with_padding_len = len(sequence_interest_with_padding)
            # save info
            info_reads.append([sample_name, region, query_name, np, rq, mc, sequence_interest, sequence_interest_with_padding, sequence_interest_len, sequence_interest_with_padding_len])
        else:
            info_reads.append([sample_name, region, query_name, np, rq, mc, 'NA', 'NA', 'NA', 'NA'])
    return info_reads

# Alignment of contigs
def alignContigs(assembly_dir_clean, cpu, ref):
    alignments = []
    # iterate over files
    for asm in assembly_dir_clean:
        print('*** Alignment of %s\t\t\t\t\t\t\t\t\t\t' %(os.path.basename(asm)), end = '\r')
        asm_out = asm.replace('.fasta', '.ref.bam')
        # do alignment
        cmd = "minimap2 -aYx asm10 %s -t %s %s | samtools sort - -@ %s -O bam -o %s" %(ref, cpu, asm, cpu, asm_out); os.system(cmd)
        # index
        cmd = 'samtools index %s' %(asm_out); os.system(cmd)
        alignments.append(asm_out)
    return alignments

# Remove duplicated contigs
def removeDuplicatedContigs(f):
    # read fasta file and create output fasta
    fasta_sequences = SeqIO.parse(open(f), 'fasta')
    # write unique sequences to fasta
    fasta_seqs = []
    with open(f.replace('.fasta', '_cleaned.fasta'), 'w') as fout:
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            if sequence not in fasta_seqs:
                fasta_seqs.append(sequence)
                fout.write('>%s\n%s\n' %(name, sequence))
    fout.close()
    # and remove the original fasta
    os.system('rm %s' %(f))
    return f.replace('.fasta', '_cleaned.fasta')

# Clean assembled contigs
def cleanContigs(assembly_dir, cpu):
    # list the haplotype-aware assemblies in the output folder
    flist = [x.rstrip() for x in os.popen('ls %s*fasta' %(assembly_dir))]
    # run in multiprocessing
    pool = multiprocessing.Pool(processes=cpu)
    extract_results = pool.map(removeDuplicatedContigs, flist)
    pool.close()
    return extract_results

# Run hifiasm assembly
def hifiasmAssembly(outDir, cpu, ref, temp_fasta, inBam, ploidy):
    # create directory for storing assemblies
    os.system('mkdir %s/hifiasm' %(outDir))
    # iterate over the samples
    for s in inBam:
        # get sample name
        s_name = os.path.basename(s).replace('.bam', '')
        print('*** Local Assembly of %s\t\t\t\t\t\t\t\t\t\t' %(s_name), end = '\r')
        # identify the relative fasta files
        fasta_interest = [x for x in temp_fasta if s_name in x]
        # run assembly
        cmd_assembly = 'hifiasm -o %s/hifiasm/%s -t %s --n-hap %s -n 2 -r 3 -N 200 %s >/dev/null 2>&1' %(outDir, s_name, cpu, ploidy, ' '.join(fasta_interest))
        try:
            os.system(cmd_assembly)
            # convert hap1, hap2 and primary contigs to fasta
            subprocess.run("gfatools gfa2fa %s/hifiasm/%s.bp.hap1.p_ctg.gfa > %s/hifiasm/%s_haps.fasta" %(outDir, s_name, outDir, s_name), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            subprocess.run("gfatools gfa2fa %s/hifiasm/%s.bp.hap2.p_ctg.gfa >> %s/hifiasm/%s_haps.fasta" %(outDir, s_name, outDir, s_name), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            subprocess.run("gfatools gfa2fa %s/hifiasm/%s.bp.p_ctg.gfa > %s/hifiasm/%s_primary.fasta" %(outDir, s_name, outDir, s_name), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            # clean temporary files
            subprocess.run("rm %s/hifiasm/%s.*" %(outDir, s_name), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            print('*** Assembly of %s done\t\t\t\t\t\t\t\t\t\t' %(s_name), end = '\r')
        except:
            print('*** Assembly of %s failed. Moving to next sample\t\t\t\t\t\t\t\t\t\t' %(s_name))
    return '%s/hifiasm/' %(outDir)

# Function to extract reads using multiple processors and convert to fasta
def samtoolsExtract(x, bam, out_dir, temp_name, type):
    # combine temporary name with the splitted bed name
    bed_ext = x.split('_bed.')[-1]
    temp_name = '%s/%s.%s' %(out_dir, bed_ext, temp_name)
    if type == 'fasta':
        temp_name = temp_name.replace('.bam', '.fa')
    # define command for the extraction
    cmd = 'samtools view -M -b -L %s %s | samtools fasta > %s' %(x, bam, temp_name) if type == 'fasta' else 'samtools view -M -b -L %s %s > %s' %(x, bam, temp_name)
    os.system(cmd)
    # and index if it's a bam file
    if type == 'bam':
        cmd = 'samtools index %s' %(temp_name); os.system(cmd)
    return temp_name

# Extract reads of interest to temporary FASTA files
def extractRead(bam_dir, bed_dir, out_dir, cpu, count_reg, type):
    # first split the bed files in n smaller bed depending on the cpu number
    split_beds = splitBed(bed_dir, cpu, out_dir, count_reg) if type == 'fasta' else bed_dir
    # make list of temporary outputs
    temp_bams = []
    for bam in bam_dir:
        # define temporary name
        temp_name = 'tmp_' + os.path.basename(bam)
        # define command to extract in MP for each of the splitted bed files
        pool = multiprocessing.Pool(processes=cpu)
        extract_fun = partial(samtoolsExtract, bam = bam, out_dir = out_dir, temp_name = temp_name, type = type)
        extract_results = pool.map(extract_fun, split_beds)
        pool.close()
        # append temporary names
        temp_bams.extend(extract_results)
    return temp_bams, split_beds

# Run TRF given a sequence
def run_trf(index, all_fasta, distances, type):
    # then run tandem repeat finder
    cmd = 'trf4.10.0-rc.2.linux64.exe %s 2 7 7 80 10 50 200 -ngs -h' %(all_fasta[index])
    trf = [x for x in os.popen(cmd).read().split('\n') if x != '']
    # loop on trf results and save them into a list of lists
    x = 0; trf_matches = []
    sample_name = re.sub(r'^a[a-z]\.tmp_', '', os.path.basename(all_fasta[index])).replace('.fa', '')
    while x < len(trf):
        # check if the line is the header of an entry
        if trf[x].startswith('@'):
            # if so, save the corresponding information depending on the type of input
            if type != 'otter' or sample_name == 'reference':
                region, sample, read_id = trf[x].split(';')
            else:
                read_id, region, seq_size_with_padding, seq_size = trf[x].split(';')
                read_id = '@>' + read_id
            x += 1
        while x < len(trf) and not trf[x].startswith('@'):
            tmp_trf_match = [read_id + '_' + region.replace('@', ''), 'NA'] + trf[x].split()
            trf_matches.append(tmp_trf_match)
            x += 1
    # finally create pandas df and assign column names
    if len(trf_matches) == 0:
        trf_matches = [['NA' for i in range(19)]] 
    df = pd.DataFrame(trf_matches)
    df.columns = ['ID', 'EXPECTED_MOTIF', 'START_TRF', 'END_TRF', 'LENGTH_MOTIF_TRF', 'COPIES_TRF', 'TRF_CONSENSUS_SIZE', 'TRF_PERC_MATCH', 'TRF_PERC_INDEL', 'TRF_SCORE', 'TRF_A_PERC', 'TRF_C_PERC', 'TRF_G_PERC', 'TRF_T_PERC', 'TRF_ENTROPY', 'TRF_MOTIF', 'TRF_REPEAT_SEQUENCE', 'TRF_PADDING_BEFORE', 'TRF_PADDING_AFTER']
    # finally, we need to add the reads where trf didn't find any motif
    if type != 'otter':
        df_seqs = pd.DataFrame(distances[index][0])
        df_seqs.columns = ['SAMPLE_NAME', 'REGION', 'READ_NAME', 'PASSES', 'READ_QUALITY', 'MAPPING_CONSENSUS', 'SEQUENCE_FOR_TRF', 'SEQUENCE_WITH_PADDINGS', 'LEN_SEQUENCE_FOR_TRF', 'LEN_SEQUENCE_WITH_PADDINGS']
    else:
        if sample_interest == 'reference__rawReads.fasta':
            sample_interest = 'reference'
        distances_sample = distances[sample_interest]
        distances_sample_df = pd.DataFrame(distances_sample)
        distances_sample_df.columns = ['REGION', 'READ_NAME', 'SEQUENCE_WITH_PADDINGS', 'LEN_SEQUENCE_WITH_PADDINGS', 'SEQUENCE_FOR_TRF', 'LEN_SEQUENCE_FOR_TRF']
        # make same identifier
        #distances_sample_df['ID'] = distances_sample_df['REGION'] + '_' + #distances_sample_df['LEN_SEQUENCE_WITH_PADDINGS'].astype(str) + '_' + #distances_sample_df['LEN_SEQUENCE_FOR_TRF'].astype(str)
        # add other columns and put NA
        distances_sample_df['PASSES'] = 'NA'; distances_sample_df['READ_QUALITY'] = 'NA'; distances_sample_df['MAPPING_CONSENSUS'] = 'NA'; distances_sample_df['WINDOW'] = 50;
    # add id
    df_seqs['ID'] = df_seqs['READ_NAME'].str.cat(df_seqs['REGION'], sep='_')
    # merge trf dataframe and reads dataframes
    complete_df = pd.merge(df_seqs, df, left_on = 'ID', right_on = 'ID', how = 'outer')
    return complete_df

### FUNCTIONS FOR CLEANING
# Function to remove temporary files
def removeTemp(outDir):
    # list all files
    all_files = [x.rstrip() for x in list(os.popen('ls %s' %(outDir)))]
    all_files = ['%s/%s' %(outDir, x) for x in all_files if 'gz' not in x]
    all_files = [x for x in all_files if 'otter_local_asm' not in x]
    all_files = [x for x in all_files if 'trf_reads' not in x]
    all_files = [x for x in all_files if 'log' not in x]
    # and remove them
    for x in all_files:
        if os.path.isfile(x):
            os.remove(x)
    # then remove the folders
    #os.system('rm -rf %s/otter_local_asm' %(outDir))
    os.system('rm -rf %s/trf_reads' %(outDir))

### FUNCTIONS FOR HAPLOTYPING
# main function that guides haplotyping
def haplotyping_steps(data, n_cpu, thr_mad, min_support, type, outDir):
    # STEP 1 IS TO ADJUST THE DATA BEFORE WE START
    data['START_TRF'] = pd.to_numeric(data['START_TRF'], errors='coerce')
    data['END_TRF'] = pd.to_numeric(data['END_TRF'], errors='coerce')
    data['LEN_SEQUENCE_FOR_TRF'] = pd.to_numeric(data['LEN_SEQUENCE_FOR_TRF'], errors='coerce')
    data['TRF_SCORE'] = pd.to_numeric(data['TRF_SCORE'], errors='coerce')
    # STEP 2 IS TO ADJUST THE MOTIFS IN THE DATA
    print('** Adjust motifs')
    data['TRF_MOTIF'] = data['TRF_MOTIF'].replace("NA", np.nan)
    all_motifs = data['TRF_MOTIF'].dropna().unique()
    main_motifs = [permutMotif(motif) for motif in all_motifs]
    motifs_df = pd.DataFrame({'motif' : all_motifs, 'UNIFORM_MOTIF' : main_motifs})
    data = pd.merge(data, motifs_df, left_on='TRF_MOTIF', right_on='motif', how='left')
    data['UNIQUE_NAME'] = data.apply(lambda row: str(row['READ_NAME']) + '___' + str(row['SAMPLE_NAME']) + '___' + str(row['REGION']) + '___' + str(row['LEN_SEQUENCE_FOR_TRF']), axis = 1)
    # STEP 3 IS TO ADJUST THE MOTIFS IN THE REFERENCE DATA
    print('** Reference motifs                                     ')
    ref = data[data['SAMPLE_NAME'] == 'reference'].copy()
    all_regions = list(ref['REGION'].dropna().unique())
    intervals = prepareIntervals(all_regions)
    ref['HAPLOTAG'] = 1; ref['POLISHED_HAPLO'] = ref['LEN_SEQUENCE_FOR_TRF']
    all_regions = list(ref['REGION'].dropna().unique())
    pool = multiprocessing.Pool(processes=n_cpu)
    motif_fun = partial(referenceMotifs, ref = ref, intervals = intervals)
    motif_res = pool.map(motif_fun, all_regions)
    pool.close()
    # combine dictionaries
    reference_motif_dic = {k: v for d in motif_res for k, v in d.items()}
    # STEP 4 IS TO ADD A UNIQUE ID AND SPLIT DUPLICATES BEFORE HAPLOTYPING
    data = data[data['SAMPLE_NAME'] != 'reference']
    data_nodup = data.drop_duplicates(subset = 'UNIQUE_NAME')
    dup_df = data[data.duplicated(subset = 'UNIQUE_NAME', keep=False)]
    # STEP 6 IS HAPLOTYPING BASED ON THE SIZES
    print('** Genotyping                                         ')
    all_samples = data_nodup['SAMPLE_NAME'].dropna().unique()
    all_regions = list(data_nodup['REGION'].dropna().unique())
    intervals = prepareIntervals(all_regions)
    sample_res = []
    for s in all_samples:
        print('**** %s                      ' %(s))
        sbs = data_nodup[(data_nodup['SAMPLE_NAME'] == s)]
        sbs_dups = dup_df[(dup_df['SAMPLE_NAME'] == s)]
        # create a dictionary of lists of lists of the rows, grouped by the 'group' column
        grouped_rows = {k: [list(row) for _, row in v.iterrows()] for k, v in sbs.groupby('REGION')}
        grouped_rows_dups = {k: [list(row) for _, row in v.iterrows()] for k, v in sbs_dups.groupby('REGION')}
        # Create a list of lists of lists from the dictionary
        list_of_lists_of_lists = [group_rows for group_rows in grouped_rows.values()]
        # idea is to create 2 lists with the same length and order for the reads and duplicated reads info
        list_of_lists_of_lists_dups = []
        for x in grouped_rows.keys():
            if x in grouped_rows_dups.keys():
                list_of_lists_of_lists_dups.append(grouped_rows_dups[x])
            else:
                list_of_lists_of_lists_dups.append([])
        # pair non-duplicated and duplicated for parallelization, assuming list_of_lists_of_lists and list_of_lists_of_lists_dups have the same length
        list_pairs = zip(list_of_lists_of_lists, list_of_lists_of_lists_dups)
        #for r in list(set(list(sbs['REGION']))):
        #    yy = [] if r not in grouped_rows_dups.keys() else grouped_rows_dups[r]
        #    xx = grouped_rows[r]
        #    haplo_results = haplotyping([xx, yy], s, thr_mad, type, reference_motif_dic, intervals, min_support)
        #    sample_res.append(haplo_results)
        pool = multiprocessing.Pool(processes=n_cpu)
        haplo_fun = partial(haplotyping, s = s, thr_mad = thr_mad, type = type, reference_motif_dic = reference_motif_dic, intervals = intervals, min_support = min_support)
        # use list_of_lists_of_lists below instead of list_pairs to restore
        haplo_results = pool.map(haplo_fun, list_pairs)
        pool.close()
        sample_res.append(haplo_results)
    # STEP 7 IS TO COMPOSE THE OUTPUTS: VCF AND SEQUENCES
    df_vcf = pd.DataFrame([x[0] for x in sample_res[0]], columns=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', all_samples[0]])
    df_seq = pd.DataFrame([x for y in sample_res[0] for x in y[1]], columns=['READ_NAME', 'HAPLOTAG', 'REGION', 'PASSES', 'READ_QUALITY', 'LEN_SEQUENCE_FOR_TRF', 'START_TRF', 'END_TRF', 'type', 'SAMPLE_NAME', 'POLISHED_HAPLO', 'DEPTH', 'CONSENSUS_MOTIF', 'CONSENSUS_MOTIF_COPIES', 'MOTIF_REF', 'REFERENCE_MOTIF_COPIES', 'SEQUENCE_WITH_PADDINGS', 'SEQUENCE_FOR_TRF'])
    # and the raw output
    if type == 'reads':
        raw_seq_list = []
        for sample in sample_res:
            for region in sample:
                if isinstance(region[-1], pd.DataFrame):
                    raw_seq_list.append(region[-1])
        # combine all dataframes and write as output
        raw_seq_df = pd.concat(raw_seq_list, ignore_index=True)
        #raw_seq_df.to_csv('%s/sample.raw.txt.gz' %(outDir), sep="\t", compression = 'gzip', header=True, index=False)
    # if there are more samples, write them as well
    if len(sample_res) >1:
        for i in range(1, len(sample_res)):
            tmp = pd.DataFrame([[x[0][2], x[0][-1]] for x in sample_res[i]], columns=['ID', all_samples[i]])
            # merge with main df
            df_vcf = pd.merge(df_vcf, tmp, on='ID', how='outer')
            # then compose the seq dataframe
            tmp = pd.DataFrame([x for y in sample_res[i] for x in y[1]], columns=['READ_NAME', 'HAPLOTAG', 'REGION', 'PASSES', 'READ_QUALITY', 'LEN_SEQUENCE_FOR_TRF', 'START_TRF', 'END_TRF', 'type', 'SAMPLE_NAME', 'POLISHED_HAPLO', 'DEPTH', 'CONSENSUS_MOTIF', 'CONSENSUS_MOTIF_COPIES', 'MOTIF_REF', 'REFERENCE_MOTIF_COPIES', 'SEQUENCE_WITH_PADDINGS', 'SEQUENCE_FOR_TRF'])
            # add to main df
            df_seq = pd.concat([df_seq, tmp], axis=0)
    # STEP 8 IS TO WRITE THE OUTPUTS
    print('** Producing outputs: VCF file and table with sequences                        ')
    seq_file = '%s/sample.seq.txt.gz' %(outDir)
    vcf_file = '%s/sample.vcf' %(outDir)
    writeOutputs(df_vcf, df_seq, seq_file, vcf_file, all_samples)
    return('Haplotyping analysis done!')

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

# function to prepare intervals
def prepareIntervals(all_regions):
    intervals = []
    for percentile in np.arange(0.05, 1, 0.05).tolist():
        index = int(len(all_regions) * percentile)
        intervals.append(all_regions[index])
    return intervals

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
    # check if the best motif is NA: if so, then try to use any motif that's in there
    if pd.isna(best_motif):
        # extract all motifs non NA motifs
        all_motifs = [value for value in list(set(haplo_data['UNIFORM_MOTIF'])) if not (isinstance(value, float) and math.isnan(value))]
        if len(all_motifs) >0:
            # consider the longest non-NA motif as the best motif
            sorted_list = sorted(all_motifs, key=len, reverse=True)
            best_motif = sorted_list[0]
            # reorder dataframe based on this
            haplo_data = haplo_data.sort_values(by='UNIFORM_MOTIF', ascending=(haplo_data['UNIFORM_MOTIF'] != best_motif).any())
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

# function to guide haplotyping
def haplotyping(pair, s, thr_mad, type, reference_motif_dic, intervals, min_support):
    # recover information for the reads and duplicates -- comment this and add dup_df as argument for the function to restore to previous, also look few lines below
    x, y = pair
    # define columns based on the data type
    columns = ['SAMPLE_NAME', 'REGION', 'READ_NAME', 'PASSES', 'READ_QUALITY', 'MAPPING_CONSENSUS', 'SEQUENCE_FOR_TRF', 'SEQUENCE_WITH_PADDINGS', 'LEN_SEQUENCE_FOR_TRF', 'LEN_SEQUENCE_WITH_PADDINGS', 'ID', 'EXPECTED_MOTIF', 'START_TRF', 'END_TRF', 'LENGTH_MOTIF_TRF', 'COPIES_TRF', 'TRF_CONSENSUS_SIZE', 'TRF_PERC_MATCH', 'TRF_PERC_INDEL', 'TRF_SCORE', 'TRF_A_PERC', 'TRF_C_PERC', 'TRF_G_PERC', 'TRF_T_PERC', 'TRF_ENTROPY', 'TRF_MOTIF', 'TRF_REPEAT_SEQUENCE', 'TRF_PADDING_BEFORE', 'TRF_PADDING_AFTER', 'HAPLOTAG', 'motif', 'UNIFORM_MOTIF', 'UNIQUE_NAME'] if type in ['reads', 'hifiasm'] else ['SAMPLE_NAME', 'REGION', 'READ_NAME', 'SEQUENCE_WITH_PADDINGS', 'LEN_SEQUENCE_WITH_PADDINGS', 'SEQUENCE_FOR_TRF', 'LEN_SEQUENCE_FOR_TRF', 'HAPLOTAG', 'PASSES', 'READ_QUALITY', 'MAPPING_CONSENSUS', 'ID', 'EXPECTED_MOTIF', 'START_TRF', 'END_TRF', 'LENGTH_MOTIF_TRF', 'COPIES_TRF', 'TRF_CONSENSUS_SIZE', 'TRF_PERC_MATCH', 'TRF_PERC_INDEL', 'TRF_SCORE', 'TRF_A_PERC', 'TRF_C_PERC', 'TRF_G_PERC', 'TRF_T_PERC', 'TRF_ENTROPY', 'TRF_MOTIF', 'TRF_REPEAT_SEQUENCE', 'TRF_PADDING_BEFORE', 'TRF_PADDING_AFTER', 'motif', 'UNIFORM_MOTIF', 'UNIQUE_NAME']
    # data of interest to dataframe
    sbs = pd.DataFrame(x, columns=columns)
    #print(list(set(list(sbs['REGION'])))[0])
    # comment out this below to restore
    dup_df = pd.DataFrame(y, columns=columns)
    # extract region
    r = list(sbs['REGION'].unique())[0]
    # exclude nas
    sbs = sbs.dropna(subset=['LEN_SEQUENCE_FOR_TRF'])
    # check if there are rows
    if sbs.shape[0] >0:
        if type in ['otter', 'hifiasm']:
            # identify haplotypes
            phased_sbs = assemblyBased_size(sbs, r)
            # polish haplotypes
            pol_sbs = polishHaplo_asm(phased_sbs, r)
            # add duplicates
            all_sbs = addDups(pol_sbs, dup_df, s, r, type)
            # finally look at the motif
            final_sbs_h1 = sampleMotifs(r, all_sbs, reference_motif_dic, 1, type)
            final_sbs_h2 = sampleMotifs(r, all_sbs, reference_motif_dic, 2, type)
            final_sbs = pd.concat([final_sbs_h1, final_sbs_h2], axis=0)
            # prepare for file writing
            tmp_vcf, tmp_seq = prepareOutputs(final_sbs, reference_motif_dic, r, type, 'None')
            all_sbs = []
            return tmp_vcf, tmp_seq, all_sbs
        elif type == 'reads':
            # check minimum support: minimum support is for alleles --> 2*min_support is the total minimum coverage required for autosomal regions. For sex-regions, we will use min_support directly
            # find chromosome to adapt coverage
            chrom = r.split(':')[0]
            minimum_coverage = min_support if chrom in ['chrY', 'Y'] else min_support*2
            # check if there is minimum support
            if int(sbs.shape[0]) >= int(minimum_coverage):
                # identify haplotypes
                warnings.filterwarnings("ignore", category=RuntimeWarning)
                pol_sbs = readBased_size(sbs, r, chrom, min_support, thr_mad)
                warnings.resetwarnings()
                # add duplicates
                all_sbs = addDups(pol_sbs, dup_df, s, r, type)
                # finally look at the motif
                final_sbs_h1, depth_h1 = sampleMotifs(r, all_sbs, reference_motif_dic, 0, type)
                final_sbs_h2, depth_h2 = sampleMotifs(r, all_sbs, reference_motif_dic, 1, type)
                final_sbs = pd.concat([final_sbs_h1, final_sbs_h2], axis=0)
                # prepare for file writing
                tmp_vcf, tmp_seq = prepareOutputs(final_sbs, reference_motif_dic, r, type, [depth_h1, depth_h2])
                return tmp_vcf, tmp_seq, all_sbs
            else:
                tmp_vcf = [chrom, r.split(':')[1].split('-')[0], r, reference_motif_dic[r][1], '.', '.', 'LOW_COVERAGE', '%s;%s' %(reference_motif_dic[r][0], reference_motif_dic[r][2]), 'QC;GT;MOTIF;CN;CN_REF;DP', 'QC_ISSUE;NA|NA;NA|NA;NA|NA;NA|NA;%s|0' %(sbs.shape[0])]
                tmp_seq = []
                all_sbs = []
                return tmp_vcf, tmp_seq, all_sbs
    else:
        chrom = r.split(':')[0]
        tmp_vcf = [chrom, r.split(':')[1].split('-')[0], r, reference_motif_dic[r][1], '.', '.', 'LOW_COVERAGE', '%s;%s' %(reference_motif_dic[r][0], reference_motif_dic[r][2]), 'QC;GT;MOTIF;CN;CN_REF;DP', 'QC_ISSUE;NA|NA;NA|NA;NA|NA;NA|NA;%s|0' %(sbs.shape[0])]
        tmp_seq = []
        all_sbs = []
        return tmp_vcf, tmp_seq, all_sbs

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
    sbs_dups = sbs_dups.dropna(subset=['LEN_SEQUENCE_FOR_TRF']).copy()
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
                haplo.append(tmp_haplo-1); size.append(tmp_size)
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

# function to make data for vcf writing
def prepareOutputs(final_sbs, reference_motif_dic, r, type, depths):
    # prepare data for VCF
    chrom, start, end = [r.split(':')[0]] + r.split(':')[-1].split('-')
    # check if 'chr' is in both the dictionary and the region of interest
    if 'chr' not in list(reference_motif_dic.keys())[0]:
        reference_motif_dic = {'chr' + key: value for key, value in reference_motif_dic.items()}
    if r in reference_motif_dic.keys():
        ref_motif, ref_len, ref_copies = reference_motif_dic[r]
    else:
        ref_motif, ref_len, ref_copies = 'NA', int(end) - int(start), 'NA'
    info_field = '%s;%s' %(ref_motif, ref_copies)
    format_field = 'QC;GT;MOTIF;CN;CN_REF;DP'
    sam_gt = '|'.join(list(final_sbs['POLISHED_HAPLO'].map(str))) if final_sbs.shape[0] >1 else list(final_sbs["POLISHED_HAPLO"].map(str).apply(lambda x: x + "|" + x))[0]
    sam_mot = '|'.join(list(final_sbs['CONSENSUS_MOTIF'].map(str))) if final_sbs.shape[0] >1 else list(final_sbs["CONSENSUS_MOTIF"].map(str).apply(lambda x: x + "|" + x))[0]
    sam_cop = '|'.join(list(final_sbs['CONSENSUS_MOTIF_COPIES'].map(str))) if final_sbs.shape[0] >1 else list(final_sbs["CONSENSUS_MOTIF_COPIES"].map(str).apply(lambda x: x + "|" + x))[0]
    sam_cop_ref = '|'.join(list(final_sbs['REFERENCE_MOTIF_COPIES'].map(str))) if final_sbs.shape[0] >1 else list(final_sbs["REFERENCE_MOTIF_COPIES"].map(str).apply(lambda x: x + "|" + x))[0]
    if type == 'reads':
        sam_depth = '|'.join([str(x) for x in depths])
    elif type == 'hifiasm':
        sam_depth = 'Na|Na'
    elif type == 'otter':
        sam_depth = '|'.join(list(final_sbs["READ_NAME"].str.split("_").str[1])) if final_sbs.shape[0] >1 else list(final_sbs["READ_NAME"].str.split("_").str[1].apply(lambda x: x + "|" + x))[0]
    sample_field = '%s;%s;%s;%s;%s;%s' %('PASS_ASM', sam_gt, sam_mot, sam_cop, sam_cop_ref, sam_depth)
    res_vcf = [chrom, start, r, ref_len, '.', '.', 'PASS', info_field, format_field, sample_field]
    # then prepare the sequence output
    if type == 'reads':
        final_sbs['DEPTH'] = [x for x in depths if x != 0]
        final_sbs['type'] = type
    elif type == 'hifiasm':
        final_sbs['DEPTH'] = ['NA' for x in range(final_sbs.shape[0])]
    elif type == 'otter':
        final_sbs['DEPTH'] = final_sbs["READ_NAME"].str.split("_").str[1]
    final_sbs['MOTIF_REF'] = ref_motif
    keep_cols = ['READ_NAME', 'HAPLOTAG', 'REGION', 'PASSES', 'READ_QUALITY', 'LEN_SEQUENCE_FOR_TRF', 'START_TRF', 'END_TRF', 'type', 'SAMPLE_NAME', 'POLISHED_HAPLO', 'DEPTH', 'CONSENSUS_MOTIF', 'CONSENSUS_MOTIF_COPIES', 'MOTIF_REF', 'REFERENCE_MOTIF_COPIES', 'SEQUENCE_WITH_PADDINGS', 'SEQUENCE_FOR_TRF']
    # subset dataframe to keep only selected columns
    df_subset = final_sbs.loc[:, keep_cols]
    res_seq = df_subset.values.tolist()    
    return res_vcf, res_seq

# function to write outputs
def writeOutputs(df_vcf, df_seq, seq_file, vcf_file, all_samples):
    # write header of vcf file
    writeVCFheader(vcf_file, all_samples)
    with open(vcf_file, mode='a') as file:
        df_vcf.to_csv(file, header=True, index=False, sep='\t')
    # then compress it
    os.system('gzip %s' %(vcf_file))
    # write sequence file
    #df_seq.to_csv(seq_file, header=True, index=False, sep='\t', compression = 'gzip')
    return

# function to write header of vcf file
def writeVCFheader(vcf_file, samples):
    # open file
    outf = open(vcf_file, 'w')
    # write header
    outf.write('##fileformat=VCFv4.2\n##INFO=<ID=REFERENCE_INFO,Number=2,Type=String,Description="Motif observed in the reference genome (GRCh38), and relative number of motif repetitions."\n##FORMAT=<ID=QC,Number=1,Type=String,Description="Quality summary of TREAT genotyping. PASS_BOTH: genotype agreed between reads-spanning and assembly. PASS_RSP: genotype from reads-spanning. PASS_ASM: genotype from assembly."\n##FORMAT=<ID=GT,Number=2,Type=String,Description="Phased size of the tandem repeats. H1_size | H2_size"\n##FORMAT=<ID=MOTIF,Number=2,Type=String,Description="Phased consensus motif found in the sample. H1_motif | H2_motif"\n##FORMAT=<ID=CN,Number=2,Type=String,Description="Phased number of copies of the motif found in the sample. H1_copies | H2_copies"\n##FORMAT=<ID=CN_REF,Number=2,Type=String,Description="Phased estimation of the reference motif as found in the sample. H1_motif_ref | H2_motif_ref"\n##FORMAT=<ID=DP,Number=1,Type=String,Description="Phased depth found in the sample. H1_depth | H2_depth"\n')
    outf.close()
    return    

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

