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
import warnings
import gzip
import pyfastx
import pyfaidx
import pytrf

### FUNCTIONS TO CHECK DIRECTORIES AND FILES
# Function to read bed file - OK
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

# Function to check directory - OK
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

# Function to check bam file(s) - OK
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

# Function to create Log file - OK
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
# Function to make assembly with otter and produce fasta files suitable for TRF - OK
def assembly_otter_opt(s, output_directory, ref_fasta, bed_file, number_threads, windowAss):
    # define output name with the right directory
    outname = s.split('/')[-1].replace('.bam', '.fa')
    # run otter -- -l was for spanning only
    cmd = 'otter assemble -c 150 --fasta -b %s -r %s -R %s %s -t %s -o %s > %s/otter_local_asm/%s' %(bed_file, ref_fasta, outname, s, number_threads, windowAss, output_directory, outname)
    os.system(cmd)
    # adjust otter sequences
    return '%s/otter_local_asm/%s' %(output_directory, outname)

# Function to write fasta files for TRF - OK
def writeFastaTRF(all_seqs, fasta_name):
    # define container for fasta outputs
    fasta_outputs = []
    # open file and write things
    with open(fasta_name, 'w') as outFile:
        for region in all_seqs:
            outFile.write('>%s;%s;%s\n%s\n' %(region[1], region[0], region[2], region[-4]))
    outFile.close()

# Function to run pytrf given a sequence on the reference genome - OK
def run_trf_ref_opt(x):
    res = []
    for k in x:
        name, seq = k[1], k[-4]
        # run approximate TRFinder
        temp = [list(i) for i in pytrf.ATRFinder(name, seq, min_motif_size=1, max_motif_size=100).as_list()]
        if len(temp) == 0:
            # if there are no hits, lower parameters and try again
            temp = [list(i) for i in pytrf.ATRFinder(name, seq, min_motif_size=1, max_motif_size=100, min_seed_repeat=2).as_list()]
            # if there are still no results, decrease parameters even lower
            if len(temp) == 0:
                temp = [list(i) for i in pytrf.ATRFinder(name, seq, min_motif_size=1, max_motif_size=100, min_seed_repeat=2, min_seed_length=8).as_list()]
        if len(temp) >0:
            for x in temp:
                x.append(seq)
                x.append(len(seq))
                x.append('REFERENCE')
                res.append(x)
        else:
            temp = [name, 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', seq, len(seq), 'REFERENCE']
            res.append(temp)
    return res

# Function to run pyTRF on otter assemblies - OK
def run_trf_asm_opt(x, w):
    sample_name = os.path.basename(x).replace('.fa', '')
    res = []
    fa = pyfastx.Fastx(x, uppercase=True)
    for name, seq in fa:
        if w == 0:
            tmp = [list(i) for i in pytrf.ATRFinder(name, seq, min_motif_size = 1, max_motif_size=100).as_list()]
            # if there are no hits, lower parameters and try again
            if len(tmp) == 0:
                tmp = [list(i) for i in pytrf.ATRFinder(name, seq, min_motif_size = 1, max_motif_size=100, min_seed_repeat=2).as_list()]
                # if there are still no results, decrease parameters even lower
                if len(tmp) == 0:
                    tmp = [list(i) for i in pytrf.ATRFinder(name, seq, min_motif_size = 1, max_motif_size=100, min_seed_repeat=2, min_seed_length=8).as_list()]
        else:
            tmp = [list(i) for i in pytrf.ATRFinder(name, seq[(w-1):-w], min_motif_size = 1, max_motif_size=100).as_list()]
            # if there are no hits, lower parameters and try again
            if len(tmp) == 0:
                tmp = [list(i) for i in pytrf.ATRFinder(name, seq[(w-1):-w], min_motif_size = 1, max_motif_size=100, min_seed_repeat=2).as_list()]
                # if there are still no results, decrease parameters even lower
                if len(tmp) == 0:
                    tmp = [list(i) for i in pytrf.ATRFinder(name, seq[(w-1):-w], min_motif_size = 1, max_motif_size=100, min_seed_repeat=2, min_seed_length=8).as_list()]
        if len(tmp) >0:
            for k in tmp:
                k.append(seq[(w-1):-w])
                k.append(len(seq[(w-1):-w]))
                k.append(sample_name)
                k.append(name.split('#')[2])
                k.append(name.split('#')[1])
                k.append(name.split('#')[-3].split(':')[-1])
                k.append(name.split('#')[-2].split(':')[-1])
                res.append(k)
        else:
            tmp = [name, 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', seq[(w-1):-w], len(seq[(w-1):-w]), sample_name, name.split('#')[2], name.split('#')[1], 'NA', 'NA']
            res.append(tmp)
    return res

# Function for otter pipeline - OK
def otterPipeline_opt(outDir, cpu, ref, bed_dir, inBam, count_reg, windowAss, window, bed):
    print('** Assembler: otter')
    # create directory for outputs
    os.system('mkdir %s/otter_local_asm' %(outDir))
    # run local assembly in multiprocessing -- optimized
    otter_start_time = time.time()
    pool = multiprocessing.Pool(processes=cpu)
    otter_fun = partial(assembly_otter_opt, output_directory = outDir, ref_fasta = ref, bed_file = bed_dir, number_threads = cpu, windowAss = windowAss)
    extract_results = pool.map(otter_fun, inBam)
    pool.close()
    otter_end_time = time.time()
    time_otter = otter_end_time - otter_start_time
    print('*** Otter took %s seconds\t\t\t\t\t\t\t\t\t\t\t\t' %(round(time_otter, 0)))
    annot_start_time = time.time()
    # do the same on the reference genome -- optimized
    all_regions = [entry[2] for chromosome in bed for entry in bed[chromosome]]
    # divide into n lists based on the number of regions
    regions_list = [all_regions[i * (len(all_regions) // cpu) + min(i, len(all_regions) % cpu):(i + 1) * (len(all_regions) // cpu) + min(i + 1, len(all_regions) % cpu)] for i in range(cpu)]
    pool = multiprocessing.Pool(processes=cpu)
    extract_fun = partial(measureDistance_reference_opt, ref = ref, w = window)
    extract_results_ref = pool.map(extract_fun, regions_list)
    pool.close()
    # run trf on reference
    pool = multiprocessing.Pool(processes=cpu)
    trf_ref = pool.map(run_trf_ref_opt, extract_results_ref)
    pool.close()
    # run trf on the assemblies
    pool = multiprocessing.Pool(processes=cpu)
    trf_asm = partial(run_trf_asm_opt, w = window)
    trf_asm_res = pool.map(trf_asm, extract_results)
    pool.close()
    # Combine df from different samples together
    # flatten the lists first
    flattened_ref = [sublist for sublist_list in trf_ref for sublist in sublist_list]
    flattened_asm = [sublist for sublist_list in trf_asm_res for sublist in sublist_list]
    # make dataframes
    df_asm = pd.DataFrame(flattened_asm, columns = ['ID', 'SEED_START', 'SEED_END', 'MOTIF', 'MOTIF_SIZE', 'SEED_REPEAT', 'ATR_START', 'ATR_END', 'ATR_REPEAT', 'ATR_SIZE', 'MATCHES', 'SUBSTITUTIONS', 'INSERTIONS', 'DELETIONS', 'IDENTITY', 'SEQUENCE', 'SEQUENCE_LEN', 'SAMPLE', 'HAPLOTYPE', 'REGION', 'TOTAL_COVERAGE', 'COVERAGE_HAPLO'])
    df_ref = pd.DataFrame(flattened_ref, columns = ['REGION', 'SEED_START', 'SEED_END', 'MOTIF', 'MOTIF_SIZE', 'SEED_REPEAT', 'ATR_START', 'ATR_END', 'ATR_REPEAT', 'ATR_SIZE', 'MATCHES', 'SUBSTITUTIONS', 'INSERTIONS', 'DELETIONS', 'IDENTITY', 'SEQUENCE', 'SEQUENCE_LEN', 'SAMPLE'])
    # add HAPLOTYPE to reference -- 1
    df_ref['HAPLOTYPE'] = 1
    # finally concatenate with reference info
    df_all = pd.concat([df_asm, df_ref])
    annot_end_time = time.time()
    time_annot = annot_end_time - annot_start_time
    print('*** Annotation took %s seconds\t\t\t\t\t\t\t\t\t\t\t\t' %(round(time_annot, 0)))
    return df_all

# Measure the distance in the reference genome - OK
def measureDistance_reference_opt(x, ref, w):
    distances = []
    # open the fasta
    fa = pyfaidx.Fasta(ref)
    # x is a list of regions -- iterate through these regions
    for i in x:
        chrom, start, end = i.split(':')[0], int(i.split(':')[-1].split('-')[0]), int(i.split(':')[-1].split('-')[1])
        tmp = fa[chrom][(start-1) : end]
        distances.append(['reference', i, 'reference', 'NA', 'NA', 'NA', str(tmp).upper(), str(tmp).upper(), len(tmp), len(tmp)])
    return distances

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
# Function that guides haplotyping - OK
def haplotyping_steps_opt(data, n_cpu, thr_mad, min_support, type, outDir, inBam):
    # Adjust the motifs in the data to find a uniform representation
    #data['UNIQUE_NAME'] = data.apply(lambda row: str(row['SAMPLE']) + '___' + str(row['REGION']) + '___' + str(row['HAPLOTYPE']), axis = 1)
    motif_start_time = time.time()
    data['MOTIF'] = data['MOTIF'].replace("NA", np.nan)
    all_motifs = data['MOTIF'].dropna().unique()
    main_motifs = [permutMotif(motif) for motif in all_motifs]
    motifs_df = pd.DataFrame({'motif' : all_motifs, 'UNIFORM_MOTIF' : main_motifs})
    data = pd.merge(data, motifs_df, left_on='MOTIF', right_on='motif', how='left')
    # Check motifs of the reference
    ref = data[data['SAMPLE'] == 'REFERENCE'].copy()
    ref['POLISHED_HAPLO'] = ref['SEQUENCE_LEN']
    ref_ok = ref.drop_duplicates(subset='REGION', keep=False).copy()
    ref_ok_dic = {row['REGION']: [row['MOTIF'], row['ATR_REPEAT'], row['SEQUENCE_LEN'], row['SEQUENCE']] for _, row in ref_ok.iterrows()}
    ref_tocheck = ref[ref.duplicated(subset='REGION', keep=False)].copy()
    # Only adjust motifs that need to be adjusted
    all_regions = list(ref_tocheck['REGION'].dropna().unique())
    pool = multiprocessing.Pool(processes=n_cpu)
    motif_fun = partial(referenceMotifs_opt, ref = ref_tocheck)
    motif_res = pool.map(motif_fun, all_regions)
    pool.close()
    # combine dictionaries
    reference_motif_dic = {k: v for d in motif_res for k, v in d.items()}
    reference_motif_dic.update(ref_ok_dic)
    # Then move to the sample(s)
    data_sample = data[data['SAMPLE'] != 'REFERENCE'].copy()
    data_sample['POLISHED_HAPLO'] = data_sample['SEQUENCE_LEN']
    # Only adjust motifs that need to be adjusted
    data_sample_ok = data_sample.drop_duplicates(subset='ID', keep=False).copy()
    data_sample_ok['CONSENSUS_MOTIF'] = data_sample_ok['UNIFORM_MOTIF']
    data_sample_ok['CONSENSUS_MOTIF_COPIES'] = data_sample_ok['ATR_REPEAT']
    # fix those that need to be fixed
    data_sample_tocheck = data_sample[data_sample.duplicated(subset='ID', keep=False)].copy()
    ids_to_fix = list(data_sample_tocheck['ID'].dropna().unique())
    pool = multiprocessing.Pool(processes=n_cpu)
    motif_fun = partial(sampleMotifs_opt, df = data_sample_tocheck)
    motif_res = pool.map(motif_fun, ids_to_fix)
    pool.close()
    motif_end_time = time.time()
    time_motif = motif_end_time - motif_start_time
    print('*** Motif merge took %s seconds\t\t\t\t\t\t\t\t\t\t\t\t' %(round(time_motif, 0)))
    prepare_start_time = time.time()
    # generate a final dataframe
    try:
        data_temp = pd.concat(motif_res, ignore_index=True)
        data_final = pd.concat([data_sample_ok, data_temp])
    except:
        data_final = data_sample_ok
    # prepare data for output
    all_samples = list(set(list(data_final['SAMPLE'])))
    all_regions = list(data_final['REGION'].dropna().unique())
    # divide in n chunks where n is th enumber of cores
    prepare_start_time = time.time()
    chunk_size = math.ceil(len(all_regions) / 4)
    chunks = [all_regions[i * chunk_size:(i + 1) * chunk_size] for i in range(n_cpu)]
    pool = multiprocessing.Pool(processes=n_cpu)
    prep_fun = partial(prepareOutputs_opt, final_sbs = data_final, reference_motif_dic = reference_motif_dic, all_samples = all_samples)
    vcf = pool.map(prep_fun, chunks)
    pool.close()
    prepare_end_time = time.time()
    time_prepare = prepare_end_time - prepare_start_time
    print('*** preparation took %s seconds\t\t\t\t\t\t\t\t\t\t\t\t' %(round(time_prepare, 0)))
    # create dataframe
    write_start_time = time.time()
    flattened_vcf = [item for sublist in vcf for item in sublist]
    df_vcf = pd.DataFrame(flattened_vcf, columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + all_samples)
    # write vcf
    vcf_file = '%s/sample.vcf' %(outDir)
    writeOutputs_opt(df_vcf, vcf_file, all_samples, inBam)
    #write_done = writeOutDirect(outDir, data_final, reference_motif_dic, all_samples, all_regions)
    write_end_time = time.time()
    time_write = write_end_time - write_start_time
    print('*** Writing took %s seconds\t\t\t\t\t\t\t\t\t\t\t\t' %(round(time_write, 0)))
    return('Haplotyping analysis done!')

# Function to make permutations - OK
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

# Function to look at reference motifs - OK
def referenceMotifs_opt(r, ref):
    # subset of reference data
    sbs = ref[ref['REGION'] == r].copy()
    # if there's only 1 motif, we are done
    if sbs.shape[0] == 1:
        sbs['CONSENSUS_MOTIF'] = sbs['UNIFORM_MOTIF']
        sbs['CONSENSUS_MOTIF_COPIES'] = sbs['ATR_REPEAT']
    elif sbs.shape[0] >1:
        # first align motifs
        sbs = sbs[sbs['HAPLOTYPE'] == 1].copy()
        sbs = motif_generalization_opt(sbs, r)
    # in the end, take only what we need to bring along
    tmp_dic = {row["REGION"]: [row["CONSENSUS_MOTIF"], row["POLISHED_HAPLO"], row['CONSENSUS_MOTIF_COPIES'], row['SEQUENCE']] for _, row in sbs.iterrows()}
    return tmp_dic

# Function to generate consensus motif using majority rule - OK
def motif_generalization_opt(haplo_data, r):
    # calculate fraction of sequence covered
    try:
        haplo_data['COVERAGE_TR'] = (haplo_data['ATR_END'] - haplo_data['ATR_START'] + 1) / haplo_data['SEQUENCE_LEN']
        haplo_data = haplo_data.sort_values('COVERAGE_TR', ascending=False).copy()
        # if >95% of the sequence is covered, stop
        high_coverage = haplo_data[haplo_data['COVERAGE_TR'] >0.95]
        if high_coverage.shape[0] == 1:
            best_motif = list(high_coverage['UNIFORM_MOTIF'])[0]
            copies = list(high_coverage['ATR_REPEAT'])[0]
            start = list(high_coverage['ATR_START'])[0]
            end = list(high_coverage['ATR_END'])[0]
        elif high_coverage.shape[0] >1:
            high_coverage['MOTIF_SCORE'] = high_coverage['COVERAGE_TR'] * high_coverage['IDENTITY']
            high_coverage = high_coverage.sort_values('MOTIF_SCORE', ascending=False)
            best_motif = high_coverage['UNIFORM_MOTIF'].iloc[0]
            copies = high_coverage['ATR_REPEAT'].iloc[0]
            start = high_coverage['ATR_START'].iloc[0]
            end = high_coverage['ATR_END'].iloc[0]
        elif haplo_data['motif'].isna().all():
            best_motif, copies, indexes, start, end = 'NA', 'NA', 'NA', 'NA', 'NA'
        else:
            best_motif, copies, indexes = combineMotifs_opt(haplo_data, r)
            start, end = indexes[0], indexes[-1]
    except:
        best_motif, copies, indexes, start, end = 'NA', 'NA', 'NA', 'NA', 'NA'
    # then combine with haplotype data
    haplo_data = haplo_data.iloc[[0]]
    haplo_data['CONSENSUS_MOTIF'] = best_motif
    haplo_data['CONSENSUS_MOTIF_COPIES'] = copies
    haplo_data['ATR_START'] = start
    haplo_data['ATR_END'] = end
    return haplo_data

# Function to combine motifs - OK
def combineMotifs_opt(haplo_data, r):
    # take the most representative motif
    haplo_data['MOTIF_SCORE'] = haplo_data['COVERAGE_TR'] * haplo_data['IDENTITY']
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
    best_motif_copies = haplo_data['ATR_REPEAT'].iloc[0]
    best_motif_index = range(int(haplo_data['ATR_START'].iloc[0]), int(haplo_data['ATR_END'].iloc[0]))
    # loop on the other motifs
    for i in range(1, haplo_data.shape[0]):
        try:
            # take corresponding motif and range
            tmp_motif = haplo_data['UNIFORM_MOTIF'].iloc[i]
            tmp_motif_copies = haplo_data['ATR_REPEAT'].iloc[i]
            tmp_motif_index = range(int(haplo_data['ATR_START'].iloc[i]), int(haplo_data['ATR_END'].iloc[i]))
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

# Function to polish haplotypes - OK
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
        # for some reasons, sometimes there are 2 haplotags and 3 contigs, clean it here
        pol_sbs = phased_sbs.drop_duplicates(subset=['HAPLOTAG'])
        pol_sbs = pol_sbs.copy()
        pol_sbs['POLISHED_HAPLO'] = pol_sbs['LEN_SEQUENCE_FOR_TRF']
    return pol_sbs

# Function to look at the motif of the samples - OK
def sampleMotifs_opt(r, df):
    # subset of data
    sb = df[df['ID'] == r].copy()
    sb_merged = motif_generalization_opt(sb, r)
    return sb_merged

# Function to make data for vcf writing - OK
def prepareOutputs_opt(chunk, final_sbs, reference_motif_dic, all_samples):
    # Ensure 'chr' is in the dictionary keys if needed
    if 'chr' not in list(reference_motif_dic.keys())[0]:
        reference_motif_dic = {'chr' + key: value for key, value in reference_motif_dic.items()}
    # Convert 'HAPLOTYPE' to numeric once
    final_sbs['HAPLOTYPE'] = pd.to_numeric(final_sbs['HAPLOTYPE'], errors='coerce')
    res_vcf = []
    columns_of_interest = ['SAMPLE', 'REGION', 'HAPLOTYPE', 'SEQUENCE', 'POLISHED_HAPLO', 'CONSENSUS_MOTIF', 'CONSENSUS_MOTIF_COPIES', 'COVERAGE_HAPLO']
    # Subset based on the chunks
    final_sbs = final_sbs[final_sbs['REGION'].isin(chunk)][columns_of_interest].copy()
    # Precompute static values
    default_sample_field = 'PASS;.|.;.|.;.|.;.|.;.|.;.|.'
    format_field = 'QC;GT;GT_LEN;MOTIF;CN;CN_REF;DP'
    for r in chunk:
        # Prepare data for VCF
        chrom, start, end = [r.split(':')[0]] + r.split(':')[-1].split('-')
        # Restrict to data of interest
        sbs = final_sbs.loc[final_sbs['REGION'] == r].copy()
        # Check if region is in reference_motif_dic
        if r in reference_motif_dic:
            ref_motif, ref_copies, ref_len, ref_seq = reference_motif_dic[r]
            try:
                motif_len = len(ref_motif.replace('+', ''))
                sbs.loc[:, 'REFERENCE_MOTIF_COPIES'] = sbs['POLISHED_HAPLO'] / motif_len
            except:
                sbs.loc[:, 'REFERENCE_MOTIF_COPIES'] = 'NA'
        else:
            ref_motif, ref_len, ref_copies = 'NA', int(end) - int(start), 'NA'
            sbs['REFERENCE_MOTIF_COPIES'] = 'NA'
            ref_seq = 'N' * (int(end) - int(start))  # Assuming ref_seq as 'N' for missing ref_seq
        info_field = f'{ref_motif};{ref_copies};{ref_len}'
        alt_seq = []
        sample_fields = []
        for s in all_samples:
            sbs_s = sbs[sbs['SAMPLE'] == s]
            if sbs_s.empty or sbs_s['HAPLOTYPE'].max() > 2:
                sample_fields.append(default_sample_field)
            else:
                gt, alt_seq = manageSequence(list(sbs_s['SEQUENCE']), alt_seq, ref_seq)
                def prepare_field(col):
                    return '|'.join(map(str, sbs_s[col])) if sbs_s.shape[0] > 1 else f'{sbs_s[col].iloc[0]}|{sbs_s[col].iloc[0]}'
                sam_gt = prepare_field('POLISHED_HAPLO')
                sam_mot = prepare_field('CONSENSUS_MOTIF')
                sam_cop = prepare_field('CONSENSUS_MOTIF_COPIES')
                sam_cop_ref = prepare_field('REFERENCE_MOTIF_COPIES')
                sam_depth = prepare_field('COVERAGE_HAPLO')
                sample_field = f'PASS;{gt};{sam_gt};{sam_mot};{sam_cop};{sam_cop_ref};{sam_depth}'
                sample_fields.append(sample_field)
        alt_seq = '.' if not alt_seq else ','.join(alt_seq)
        tmp_vcf = [chrom, start, r, ref_seq, alt_seq, '.', '.', info_field, format_field] + sample_fields
        res_vcf.append(tmp_vcf)
    return res_vcf

# Function to manage reference and alternative sequences - OK
def manageSequence(sequences, alt, ref_seq):
    gt = []
    for seq in sequences:
        # check if the same as the reference
        if seq == ref_seq:
            gt.append('0')
        else:
            # if not reference, check if already in alt list
            seq = 'N' if seq == '' else seq
            if seq in alt:
                gt.append(str(alt.index(seq)))
            else:
                alt.append(seq)
                gt.append(str(alt.index(seq) + 1))
    # then check for homozygous
    if len(gt) == 1:
        gt = gt + gt
    gt = '|'.join(gt)
    return gt, alt

# Function to write outputs - OK
def writeOutputs_opt(df_vcf, vcf_file, all_samples, inBam):
    # write header of vcf file
    writeVCFheader(vcf_file, all_samples, inBam)
    with open(vcf_file, mode='a') as file:
        df_vcf.to_csv(file, header=True, index=False, sep='\t')
    # then compress it
    os.system('gzip %s' %(vcf_file))
    # write sequence file
    #df_seq.to_csv(seq_file, header=True, index=False, sep='\t', compression = 'gzip')
    return

# Function to write header of vcf file - OK
def writeVCFheader(vcf_file, samples, inBam):
    # open file
    outf = open(vcf_file, 'w')
    # write header
    outf.write('##fileformat=VCFv4.2\n##INFO=<ID=REFERENCE_INFO,Number=2,Type=String,Description="Motif observed in the reference genome (GRCh38), and relative number of motif repetitions."\n##FORMAT=<ID=QC,Number=1,Type=String,Description="Quality summary of TREAT genotyping. PASS: passed quality filter."\n##FORMAT=<ID=GT,Number=2,Type=String,Description="Phased genotype of the tandem repeats. H1_genotype | H2_genotype"\n##FORMAT=<ID=GT_LEN,Number=2,Type=Number,Description="Phased size of the tandem repeat genotypes. H1_size | H2_size"\n##FORMAT=<ID=MOTIF,Number=2,Type=String,Description="Phased consensus motif found in the sample. H1_motif | H2_motif"\n##FORMAT=<ID=CN,Number=2,Type=String,Description="Phased number of copies of the motif found in the sample. H1_copies | H2_copies"\n##FORMAT=<ID=CN_REF,Number=2,Type=String,Description="Phased estimation of the reference motif as found in the sample. H1_motif_ref | H2_motif_ref"\n##FORMAT=<ID=DP,Number=1,Type=String,Description="Phased depth found of the tandem repeat. H1_depth | H2_depth"\n')
    # need to add the contig information
    contig_info = '\n'.join([convert_sq_to_contig(x.rstrip())for x in os.popen('samtools view -H %s' %(inBam[0])) if '@SQ' in x])
    outf.write('%s\n' %(contig_info))
    outf.close()
    return

# Function to convert sequence to contig - OK
def convert_sq_to_contig(sq_string):
    match = re.match(r'@SQ\tSN:(\S+)\tLN:(\d+)', sq_string)
    if not match:
        raise ValueError("The input string does not match the expected format")
    return f'##contig=<ID={match.group(1)},length={match.group(2)}>'

# Function to find sequences of consecutive integers - OK
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

