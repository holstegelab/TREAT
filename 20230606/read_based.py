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

# Function to extract reads using multiple processors
def samtoolsExtract(x, bam, out_dir, temp_name):
    # combine temporary name with the splitted bed name
    bed_ext = x.split('_bed.')[-1]
    temp_name = '%s/%s.%s' %(out_dir, bed_ext, temp_name)
    # define command for the extraction
    cmd = 'samtools view -M -b -L %s %s > %s' %(x, bam, temp_name)
    os.system(cmd)
    # and index
    cmd = 'samtools index %s' %(temp_name)
    os.system(cmd)
    return temp_name

# Function to split bed file in n bed files
def splitBed(bed_dir, n, outDir, count_reg):
    number_lines_per_file = int(count_reg/n)
    cmd = 'split -l %s %s %s/tmp_bed.' %(number_lines_per_file, bed_dir, outDir)
    os.system(cmd)
    # then read the obtained temporary beds
    cmd = 'ls %s/tmp_bed.*' %(outDir)
    tmp_beds = [x.rstrip() for x in list(os.popen(cmd))]
    return tmp_beds

# Extract reads of interest to temporary BAM files
def extractRead(bam_dir, bed_dir, out_dir, cpu, count_reg):
    # first split the bed files in n smaller bed depending on the cpu number
    split_beds = splitBed(bed_dir, cpu, out_dir, count_reg)
    # make list of temporary outputs
    temp_bams = []
    for bam in bam_dir:
        # define temporary name
        temp_name = 'tmp_' + os.path.basename(bam)
        # define command to extract in MP for each of the splitted bed files
        pool = multiprocessing.Pool(processes=cpu)
        extract_fun = partial(samtoolsExtract, bam = bam, out_dir = out_dir, temp_name = temp_name)
        extract_results = pool.map(extract_fun, split_beds)
        # append temporary names
        temp_bams.extend(extract_results)
        print('**** Reads overlapping BED file extracted')
    return temp_bams, split_beds

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

# Function to execute the sequence extraction in multiple processors
def distributeExtraction(x, bed, window):
    # container for results
    tmp_results = []
    # get sample name
    sample_name = re.sub(r'^a[az]\.tmp_', '', os.path.basename(x)).replace('.bam', '')
    # loop over reads with pysam
    with pysam.AlignmentFile(x, 'rb', check_sq=False) as bamfile:
        for read in bamfile:
            # extract read information
            ref_chrom, ref_start, ref_end, query_name, query_sequence, cigartuples, tags, is_secondary, is_supplementary = read.reference_name, int(read.reference_start), int(read.reference_end), read.query_name, read.query_sequence, read.cigartuples, read.tags, read.is_supplementary, read.is_secondary
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

# Function to write fasta files for TRF
def writeFastaTRF(all_seqs, fasta_name):
    # define container for fasta outputs
    fasta_outputs = []
    # open file and write things
    with open(fasta_name, 'w') as outFile:
        for region in all_seqs:
            outFile.write('>%s;%s;%s\n%s\n' %(region[1], region[0], region[2], region[-4]))
    outFile.close()

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

# Run TRF given a sequence
def run_trf(index, all_fasta, distances, type):
    # then run tandem repeat finder
    cmd = 'trf4.10.0-rc.2.linux64.exe %s 2 7 7 80 10 50 200 -ngs -h' %(all_fasta[index])
    trf = [x for x in os.popen(cmd).read().split('\n') if x != '']
    # loop on trf results and save them into a list of lists
    x = 0; trf_matches = []
    sample_name = re.sub(r'^a[az]\.tmp_', '', os.path.basename(all_fasta[index])).replace('.fa', '')
    while x < len(trf):
        # check if the line is the header of an entry
        if trf[x].startswith('@'):
            # if so, save the corresponding information depending on the type of input
            if type != 'otter' or sample_interest == 'reference':
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

# Main
# Read arguments and make small changes
inBam_dir, bed_dir, outDir, ref, window, cpu, phasingData, mappingSNP, HaploDev, minimumSupport, minimumCoverage = sys.argv[1::]
window = int(window); cpu = int(cpu)

# 1. Check arguments: BED, output directory and BAMs
print('** Analysis started')
ts_total = time.time()
# 1.1 Check output directory
print(checkOutDir(outDir))
# 1.2 Read bed file
bed, count_reg = readBed(bed_dir)
# 1.3 Check BAM files
inBam = checkBAM(inBam_dir)

# 2. Extract sequence of interest
ts = time.time()
# 2.1 Extract reads using samtools
temp_bams, temp_beds = extractRead(inBam, bed_dir, outDir, cpu, count_reg)
# 2.2 Parse output and get sequences
pool = multiprocessing.Pool(processes=cpu)
extract_fun = partial(distributeExtraction, bed = bed, window = window)
extract_results = pool.map(extract_fun, temp_bams)
all_fasta = [outer_list[1] for outer_list in extract_results]
# 2.3 Then do the same on the reference genome
pool = multiprocessing.Pool(processes=cpu)
extract_fun = partial(measureDistance_reference, window = window, ref = ref, output_directory = outDir)
extract_results_ref = pool.map(extract_fun, temp_beds)
all_fasta_ref = [outer_list[1] for outer_list in extract_results]
# 2.5 combine reference with other samples
extract_results.extend(extract_results_ref)
all_fasta.extend(all_fasta_ref)
te = time.time()
time_extraction = te-ts
print('**** Operation took %s seconds                                 ' %(round(time_extraction, 0)))

# 3. TRF
ts = time.time()
# 3.1 Run TRF in multiprocessing for each sample
pool = multiprocessing.Pool(processes=cpu)
trf_fun = partial(run_trf, all_fasta = all_fasta, distances = extract_results, type = 'reads')
index_fasta = [x for x in range(len(all_fasta))]
trf_results = pool.map(trf_fun, index_fasta)
# 3.2 combine df from different samples together
df_trf_combined = pd.concat(trf_results)
print('**** TRF done on all reads and samples')
te = time.time()
time_trf = te-ts
print('**** Operation took %s seconds                                 ' %(round(time_trf, 0)))

# 4. Phasing and haplotagging and combine with sequences
ts = time.time()
# 4.1 Check whether we need to do this
if phasingData == 'None':
    print('**** Phasing NOT selected (not specified any SNP data)')
    combined_haplotags_df = pd.DataFrame(columns=['READ_NAME', 'HAPLOTAG'])
else:
    # check this
    print('**** Phasing and haplotagging')
    os.system('mkdir %s/phasing' %(output_directory))
    print('**** finding SNPs for phasing')
    snps_for_phasing = find_SNPs_Samples_plink(snp_dir, output_directory, snp_data_ids, reads_bam)
    print('**** based on mapping ids provided, phasing will be done on %s/%s samples' %(len(snps_for_phasing), len(reads_bam)))
    print('**** start phasing                                  ', end = '\r')
    pool = multiprocessing.Pool(processes=number_threads)
    phasing_fun = partial(phase_reads_MP, output_directory = output_directory, SNPs_data_directory = snp_dir, ref_path = ref_fasta, bed_file = bed_file, window = window_for_phasing)
    ts = time.time()
    phasing_results = pool.map(phasing_fun, snps_for_phasing)
    te = time.time()
    time_phasing = te-ts
    print('**** phasing done in %s seconds                                       ' %(round(time_phasing, 0)))
    combined_haplotags = sum(phasing_results, [])
    combined_haplotags_df = pd.DataFrame(combined_haplotags)
    combined_haplotags_df.columns = ['READ_NAME', 'HAPLOTAG']
te = time.time()
time_phasing = te-ts
# 4.2 Combine with sequences
df_trf_phasing_combined = pd.merge(df_trf_combined, combined_haplotags_df, left_on = 'READ_NAME', right_on = 'READ_NAME', how = 'outer')
print('**** Operation took %s seconds                                 ' %(round(time_phasing, 0)))

# 5. Output
ts = time.time()
# 5.2 Output file for haplotyping
outf = '%s/spanning_reads_trf_phasing.txt' %(outDir)
df_trf_phasing_combined.to_csv(outf, sep = "\t", index=False, na_rep='NA', compression='gzip')
print('**** Data combined and outputs are ready')
te = time.time()
time_write = te-ts
print('**** Operation took %s seconds                                 ' %(round(time_write, 0)))

# 6. Haplotyping
ts = time.time()
file_path = os.path.realpath(__file__)
file_path = '/'.join(file_path.split('/')[:-1])
main_script = '%s/call_haplotypes.py %s/spanning_reads_trf_phasing.txt %s %s %s %s %s ' %(file_path, outDir, outDir, cpu, HaploDev, 'reads', minimumSupport)
print('\n', main_script, '\n')
os.system('%s/call_haplotypes.py %s/spanning_reads_trf_phasing.txt %s %s %s %s %s ' %(file_path, outDir, outDir, cpu, HaploDev, 'reads', minimumSupport))
te = time.time()
time_write = te-ts
print('**** Operation took %s seconds                                 ' %(round(time_write, 0)))
te_total = time.time()
time_total = te_total - ts_total
print('\n** Analysis completed in %s seconds. Ciao!                   ' %(round(time_total, 0)))