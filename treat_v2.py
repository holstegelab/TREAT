#!/usr/bin/env python

#############################################################
# SET OF TOOLS THAT CAN BE USED TO MANAGE AND ANALYZE DATA  #
# FROM PACBIO OR OTHER SEQUENCING PLATFORMS.                #
#############################################################

## Libraries
import os; from inspect import getsourcefile; from os.path import abspath; import sys; import pathlib; import time; import itertools
from functions_treat_v2 import *; import argparse; import os.path; import pysam; from Bio.Seq import Seq
from Bio import SeqIO; from random import random; import gzip; import subprocess; import json; import pandas as pd
from datetime import datetime; from functools import partial; from itertools import repeat; import multiprocessing

## Main
## Start time
start_time = time.time()
## Define Arguments
parser = argparse.ArgumentParser(description = 'Find information about a specific region/tandem repeat')

## Main arguments
# analysis type
parser.add_argument('--analysis-type', dest = 'analysis_type', type = str, help = 'Type of analysis to perform [reads_spanning_trf / assembly_trf / haplotyping / extract_snps (beta) / annotate_snps (beta) / complete (beta) / extract_annotate (beta) / genotype_snps_pacbio (beta) / coverage_profile (beta) / extract_raw_reads (beta) / realign (beta)]. See docs for further information.', required = True)
# bed file
parser.add_argument('--bed', dest = 'bed_dir', type = str, help = '.bed file containing the region(s) to look. Header is not required but if present, it must starts with #.', required = False, default = 'None')
# bam file
parser.add_argument('--bam-dir', dest = 'bam_dir', type = str, help = 'Directory of bam file(s). If a directory is provided, will use all .bam in the directory. If a single .bam file is provided, will use that file.', required = False, default = 'None')
# output directory
parser.add_argument('--out-dir', dest = 'out_dir', type = str, help = 'Directory where to place output files. If the directory exists, will place files in, otherwise will create the folder and place results in.', required = False, default = 'None')
# reference genome
parser.add_argument('--ref', dest = 'ref', type = str, help = 'Path to reference genome data.', required = False, default = 'None')
# flanking window to use
parser.add_argument('--window', dest = 'window', type = int, help = 'Integer. Will use this value to take reads surrounding the provided .bed file.', required = False, default = 10)
# whether to store temporary files (reads aligning to bed file)
parser.add_argument('--store-temp', dest = 'store_temp', type = str, help = 'Boolean (True/False). If True, will store all temporary .bam and .fasta files.', required = False, default = 'False')
# number of threads overall
parser.add_argument('--thread', dest = 'thread', type = int, help = 'Number of parallel threads to be used', required = False, default = 1)

## Assembly arguments
# numer of threads for assembly and alignment
parser.add_argument('--thread-asm-aln', dest = 'thread_asm_aln', type = str, help = 'Number of parallel threads to be used during assembly and alignment', required = False, default = "default")
# ploidy
parser.add_argument('--assembly-ploidy', dest = 'ass_ploidy', type = int, help = 'Ploidy to be used for the local assembly procedure. Default value is 2 (for diploid organisms).', required = False, default = 2)
# assembly type
parser.add_argument('--assembly-type', dest = 'ass_type', type = str, help = 'Type of local assembly to perform. By default, each .bam will result in an assembly. If you prefer to use multiple .bam files for an assembly, please submit a file with no header and 2 columns: the first column should report, for each line, a comma-separated list of .bam files to combine in the assembly. The second column, for each line, should report the output prefix of the final assembly for each group.', required = False, default = 'asm_per_bam')
parser.add_argument('--assembly-tool', dest = 'ass_tool', type = str, help = '[hifiasm / otter]: the tool to use for local assembly.', required = False, default = 'hifiasm')
# window for assembly
parser.add_argument('--window-asm', dest = 'window_asm', type = int, help = 'Window to use to recruit reads for assembly.', required = False, default = 5000)

## Phasing arguments
# path to snp data
parser.add_argument('--snp-data', dest = 'snp_dir', type = str, help = 'If phasing is selected, please add here the path to SNP data (in PLINK2 format). This information is necessary for phasing.', required = False, default = 'False')
# path to mapping file between snp data and read
parser.add_argument('--snp-data-ids', dest = 'snp_data_ids', type = str, help = 'Please submit here a 2-column file with GWAS ID and ID in sequencing data. If not provided, will assume the IDs are the same.', required = False, default = 'False')

## Haplotyping analysis
# path to trf file for haplotyping
parser.add_argument('--trf', dest = 'trf_file', type = str, help = 'If analysis type was haplotyping, please input here the path to the TRF output file.', required = False, default = 'None')
# path to phasing file for haplotyping
parser.add_argument('--phase', dest = 'phase_file', type = str, help = 'If analysis type was haplotyping, please input here the path to the PHASING output file.', required = False, default = 'None')
# path to assembly file for haplotyping
parser.add_argument('--asm', dest = 'asm_file', type = str, help = 'If analysis type was haplotyping, please input here the path to the TRF output of assembly.', required = False, default = 'None')
# deviation value for haplotyping to call haplotypes
parser.add_argument('--haplotyping-deviation', dest = 'thr_mad', type = float, help = 'During haplotying analysis, median absolute deviation to assign reads to the same allele.', required = False, default = 0.10)

## Other analyses parameters
# path to file including SNPs with annotate or genotype
parser.add_argument('--variant-file', dest = 'variant_file', type = str, help = 'If the analysis_type is annotate_snps or genotype_snps_pacbio, please provide here the path to the file including the SNPs to annotate/extract.', required = False, default = 'None')
# path to directory containing fasta files to be (re)aligned
parser.add_argument('--fasta-dir', dest = 'fasta_dir', type = str, help = 'Directory of fasta file(s). If a directory is provided, will use all .fasta/fa in the directory. If a single .fasta file is provided, will use that file.', required = False, default = 'None')
# whether to polish reads or not using Alex algorithm
parser.add_argument('--polish', dest = 'polish', type = str, help = 'Boolean (True/False). If True, reads containing the region of interest will also be polished.', required = False, default = 'False')
# coverage step for coverage analysis -- coverage analysis to be changed with the samtools
parser.add_argument('--coverage-step', dest = 'step', type = int, help = 'Number of nt based on which the region(s) of interest will be split to calculate coverage.', required = False, default = 500)
# path to file containing read identifiers to extract
parser.add_argument('--reads-ids', dest = 'target_reads', type = str, help = 'If analysis type was extract_raw_reads, please provide here a tab-separated file with two columns: the first column should contain the .bam file to extract reads from, the second column should contain the reads ID to extract, one read ID per line.', required = False, default = 'False')

## Load arguments
args = parser.parse_args()

## Check arguments and throw errors when some require argument is not provided
# Throw error when out_dir is not specified (apart from when analysis type is realign)
if (args.analysis_type != 'realign') and (args.out_dir == 'None'):
    parser.error('!! You should provide the desired output directory if analysis_type is --> %s' %(args.analysis_type))
# Throw error when bam_dir is not specified for a wide range of analyses (except those specified in the list below)
if (args.analysis_type not in ['haplotyping', 'annotate_snps', 'realign_assembly', 'extract_snps', 'extract_annotate', 'realign']) and (args.bam_dir == 'None'):
    parser.error('!! You should provide at least a .bam file if analysis_type is --> %s' %(args.analysis_type))
# Throw error when bed_dir is not specified for a wide range of analyses (except those specified in the list below)
elif (args.analysis_type not in ['haplotyping', 'genotype_snps_pacbio', 'annotate_snps', 'extract_raw_reads', 'realign']) and (args.bed_dir == 'None'):
    parser.error('!! You should provide at least a .bed file if analysis_type is --> %s' %(args.analysis_type))
# Throw error when fasta_dir is not specified and analysis type is realign
if (args.analysis_type == 'realign') and (args.fasta_dir == 'None'):
    parser.error('!! You should provide the input fasta file(s) directory if analysis_type is --> %s' %(args.analysis_type))
# Throw error when target_reads is not specified and analysis type is extract_raw_reads
if (args.analysis_type == 'extract_raw_reads' and args.target_reads == 'False'):
    parser.error('!! You should provide either a text file containing the IDs of interest or a comma-separated list of IDs')
# Throw error when neither trf_file nor asm_trf are not specified and analysis type is haplotyping
if (args.analysis_type == 'haplotyping' and args.trf_file == 'None' and args.asm_file == 'None'):
    parser.error('!! You should provide at least one TRF-output of single-reads or assembly (or both) when the analysis type is haplotyping.')
# Throw error when measure is specified but no reference data is provided
if (args.analysis_type in ['measure', 'trf', 'assembly', 'realign', 'phase_reads'] and args.ref == 'None'):
    parser.error('!! You should provide the path to the reference genome when the analysis type is measure.')
# Throw error when phasing analysis is specific but no SNP data is provided
if (args.analysis_type in ['phase_reads'] and args.snp_dir == 'False'):
    parser.error('!! You should provide the path to the SNP data when the analysis type is phase_reads.')

## Print arguments
# Introduction message
print("\n**********************************************")
print("** Tandem REpeat Annotation Toolkit (Treat) **")
print('******** put together by Niccolo Tesi ********')
print("**********************************************\n")
print("Your settings in use:")
# Main arguments depending on the analysis type
print("** analysis type --> %s" %(args.analysis_type))
# Haplotyping analysis
if args.analysis_type == 'haplotyping':
    print("** trf file --> %s" %(args.trf_file))
    print("** asm file --> %s" %(args.asm_file))
    print("** phasing file --> %s" %(args.phase_file))
    print("** median absolute deviation --> %s" %(args.thr_mad))
# Realignment analysis
elif args.analysis_type == 'realign':
    print("** fasta file(s) --> %s" %(args.fasta_dir))
    print('** reference genome --> %s' %(args.ref))
    print("** output folder is the same as input fasta file(s)")
# Complete analysis
elif args.analysis_type == 'complete':
    print("** bam file(s) --> %s" %(args.bam_dir))
    print("** output folder --> %s" %(args.out_dir))
    print("** bed file --> %s" %(args.bed_dir))
    print("** intermediate files --> %s" %(args.store_temp))
    print('** reference genome --> %s' %(args.ref))
    print('** number of cpus --> %s' %(args.thread))
    print("** window used --> %s" %(args.window))
    print('** polishing --> %s' %(args.polish))
    print("** assembly type --> %s" %(args.ass_type))
    print("** assembly ploidy --> %s" %(args.ass_ploidy))
    print("** window used for assembly --> %s" %(args.window_asm))
    print("** number of cpus for assembly/alignment --> %s" %(args.thread_asm_aln))
    print("** phasing using SNPs in --> %s" %(args.snp_dir))
    print("** window for recruiting SNPs for phasing --> 2500bp")
    print("** mapping IDs between long-read and SNPs using --> %s" %(args.snp_data_ids))
    print("** step for coverage --> %s" %(args.step))
    print("** median absolute deviation --> %s" %(args.thr_mad))
# Coverage analysis
elif args.analysis_type == 'coverage_profile':
    print("** bam file(s) --> %s" %(args.bam_dir))
    print("** output folder --> %s" %(args.out_dir))
    print("** bed file --> %s" %(args.bed_dir))
    print("** intermediate files --> %s" %(args.store_temp))
    print('** reference genome --> %s' %(args.ref))
    print('** number of cpus --> %s' %(args.thread))
    print("** step for coverage --> %s" %(args.step))
# Reads-spanning-TRF analysis
elif args.analysis_type == 'reads_spanning_trf':
    print("** bam file(s) --> %s" %(args.bam_dir))
    print("** output folder --> %s" %(args.out_dir))
    print("** bed file --> %s" %(args.bed_dir))
    print("** intermediate files --> %s" %(args.store_temp))
    print('** reference genome --> %s' %(args.ref))
    print('** number of cpus --> %s' %(args.thread))
    print("** window used --> %s" %(args.window))
    print('** polishing --> %s' %(args.polish))
    print("** phasing using SNPs in --> %s" %(args.snp_dir))
    print("** mapping IDs between long-read and SNPs using --> %s" %(args.snp_data_ids))
    print("** window for recruiting SNPs for phasing --> 2500bp")
    print("** median absolute deviation --> %s" %(args.thr_mad))
# Assembly-TRF analysis
elif args.analysis_type == 'assembly_trf':
    print("** bam file(s) --> %s" %(args.bam_dir))
    print("** output folder --> %s" %(args.out_dir))
    print("** bed file --> %s" %(args.bed_dir))
    print("** intermediate files --> %s" %(args.store_temp))
    print('** reference genome --> %s' %(args.ref))
    print('** number of cpus --> %s' %(args.thread))
    print("** window used for TRF --> %s" %(args.window))
    print('** polishing --> %s' %(args.polish))
    print("** assembly type --> %s" %(args.ass_type))
    print("** assembly ploidy --> %s" %(args.ass_ploidy))
    print("** window used for assembly --> %s" %(args.window_asm))
    print("** number of cpus for assembly/alignment --> %s" %(args.thread_asm_aln))
    print("** median absolute deviation --> %s" %(args.thr_mad))
# Extract SNPs analysis
elif args.analysis_type in ['extract_snps', 'extract_annotate']:
    print("** bam file(s) --> %s" %(args.bam_dir))
    print("** output folder --> %s" %(args.out_dir))
    print('** variant file --> %s' %(args.variant_file))
# Annotate SNPs analysis
elif args.analysis_type == 'annotate_snps':
    print('** variant file --> %s' %(args.variant_file))
# SNP calling analysis
elif args.analysis_type == 'genotype_snps_pacbio':
    print("** bam file(s) --> %s" %(args.bam_dir))
    print('** variant file --> %s' %(args.variant_file))
# Put some space before actual analysis
print('********************\n')

## Store arguments
bed_file, anal_type, var_file, bam_directory, output_directory, store_temporary = args.bed_dir, args.analysis_type, args.variant_file, args.bam_dir, args.out_dir, args.store_temp
window_size, assembly_type, assembly_ploidy, number_threads, polishing, snp_dir = args.window, args.ass_type, args.ass_ploidy, args.thread, args.polish, args.snp_dir
snp_data_ids, step, target_reads, fasta_dir, trf_file, phase_file = args.snp_data_ids, args.step, args.target_reads, args.fasta_dir, args.trf_file, args.phase_file
asm_file, ref_fasta, thr_mad, number_threads_asm_aln, window_asm, ass_tool = args.asm_file, args.ref, args.thr_mad, args.thread_asm_aln, args.window_asm, args.ass_tool
# window for phasing should be larger
window_for_phasing = 2500

## Store arguments (for debugging only)
#bed_file, anal_type, var_file, bam_directory, output_directory, store_temporary, window_size, assembly_type, assembly_ploidy, number_threads, polishing, snp_dir, snp_data_ids, step, target_reads = '/project/holstegelab/Share/nicco/workspaces/20211013_target_approach/automatic_pipeline/DMPK/dmpk.bed', '', '', '/project/holstegelab/Share/nicco/workspaces/20211013_target_approach/automatic_pipeline/DMPK/extract_reads/case/DNA15-20132_2.haplotagged_step1_hifi.bam', '/project/holstegelab/Share/nicco/workspaces/20211013_target_approach/automatic_pipeline/DMPK/phase_reads/case', 'True', 100000, '', 2, 4, 'True', '/project/holstegelab/Share/pacbio/radbound_rfc1_cases/DNA15-20132-DMPK/Analyzed/GRCh38_20220408/SNVCalling_20220411_deepvariant//DNA15-20132_2.phased.pvar', '/project/holstegelab/Share/nicco/workspaces/20211013_target_approach/automatic_pipeline/DMPK/phase_reads/case/map.txt', 500, '/project/holstegelab/Share/nicco/workspaces/20211013_target_approach/automatic_pipeline/RFC1/subreads/c7_all_reads.txt'

## Check whether output directory exists otherwise create it
print("** checking directories")
print(checkOutDir(output_directory, anal_type))

## Add log file with run information in the output folder
print(addLogRun(bed_file, anal_type, var_file, bam_directory, output_directory, store_temporary, window_size, assembly_type, assembly_ploidy, number_threads, polishing, snp_dir, snp_data_ids, step, target_reads, ref_fasta, fasta_dir, trf_file, phase_file, asm_file))

## Different analyses mode
if anal_type == 'extract_snps':
    print("** reading bed file")
    bed = readBed(args.bed_dir)
    snps_topmed = extractSNPs(bed, args.window)
    writeSNPs(args.out_dir, snps_topmed)
elif anal_type == 'annotate_snps':
    input_variants, variant_info = readSNPs(args.variant_file)
    print('** input variants correctly read')
    snps_topmed_coding = annotateSNPs(input_variants, args.out_dir, args.analysis_type)
    writeSNPsAnnotation(args.out_dir, snps_topmed_coding)
elif anal_type == 'extract_annotate':
    print("** reading bed file")
    bed = readBed(args.bed_dir)
    snps_topmed = extractSNPs(bed, args.window)
    writeSNPs(args.out_dir, snps_topmed)
    snps_topmed_coding = annotateSNPs(snps_topmed, args.out_dir, args.analysis_type)
    writeSNPsAnnotation(args.out_dir, snps_topmed_coding)
elif anal_type == 'genotype_snps_pacbio':
    input_variants, variant_info = readSNPs(args.variant_file)
    print('** input variants correctly read')
    print("** checking bam files")
    all_bams = checkBAM(args.bam_dir)
    pacbio_genotypes = extractSNPsFromBam(input_variants, all_bams, variant_info)
    print("\n** calling genotypes -- using minimun coverage of 5 reads")
    calledGenotypes = callGenotypes(pacbio_genotypes, 6, 0.05)
    print("** creating output")
    writeGenotypes(args.out_dir, calledGenotypes)
elif anal_type == 'coverage_profile':
    print("** reading bed file")
    bed = readBed(bed_file)
    print("** checking bam files")
    all_bams = checkBAM(bam_directory)
    print("** generating coverage profile")
    pool = multiprocessing.Pool(processes=number_threads)
    coverage_fun = partial(generateCoverageProfile_MP, bed = bed, all_bams = all_bams, window_size = window_size, step = step, output_directory = output_directory)
    coverage_results = pool.map(coverage_fun, all_bams)
    print('**** done with coverage profiles                                     ')
    # combine results and output files
    coverage_info = {k:v for element in coverage_results for k,v in element.items()}
    outname = open('%s/coverage_profiles.bed' %(output_directory), 'w')
    outname.write('CHROM\tSTART_POS\tEND_POS\tSAMPLE\tCOVERAGE\tREGION_ID\n')
    for s in coverage_info.keys():
        for x in coverage_info[s]:
            outname.write('%s\t%s\t%s\t%s\t%s\t%s\n' %(x[0], x[1], x[2], x[3], x[4], x[5]))
    outname.close()
elif anal_type == 'complete':
    # 1. read bed regions
    print("**** reading bed file")
    bed_regions = readBed(bed_file)

    # 2. check bam files
    print("** checking bam files")
    all_bam_files = checkBAM(bam_directory)
    
    # 3. read motif
    motif = readMotif(bed_file)
    
    # 4. extract reads mapping to regions of interest in bed -- multiprocessing
    print("\n** 1. extracting reads")
    os.system('mkdir %s/raw_reads' %(output_directory))
    pool = multiprocessing.Pool(processes=number_threads)
    extract_fun = partial(extractReads_MP, bed = bed_regions, out_dir = '%s/raw_reads' %(output_directory), window = window_size, all_bams = all_bam_files)
    extract_results = pool.map(extract_fun, all_bam_files); print('**** read extraction done!                                         ')
    
    # 5. combine results together
    reads_bam = [x[0] for x in extract_results]; reads_fasta = [x[1] for x in extract_results]
    
    # 6. calculate size of the regions of interest -- multiprocessing
    print('** 2. calculate size of the regions of interest')
    pool = multiprocessing.Pool(processes=number_threads)
    measure_fun = partial(measureDistance_MP, bed = bed_regions, window = 1)
    extract_results = pool.map(measure_fun, reads_bam)
    print('**** done measuring reference                                     ', end = '\r')
    dist_reference = measureDistance_reference(bed_regions, 1, ref_fasta); print('**** read measurement done!                                         ')
    extract_results.append(dist_reference)
    # combine results
    distances = {k:v for element in extract_results for k,v in element.items()}
    
    # 7. TRF on single-reads
    print('** 3. tandem repeat finder on the single-reads')
    all_bams = list(distances.keys()); trf_out_dir = '%s/trf_reads' %(output_directory); os.system('mkdir %s' %(trf_out_dir));
    pool = multiprocessing.Pool(processes=number_threads)
    trf_fun = partial(trf_MP, out_dir = trf_out_dir, motif = motif, polished = 'False', distances = distances)
    trf_results = pool.map(trf_fun, all_bams)
    # combine df from different samples together
    df_trf = pd.concat(trf_results)
    print('**** done running TRF                                     ')

    # 8. combine these results and make output
    outf = output_directory + '/trf_reads/measures_spanning_reads_and_trf.txt'
    df_trf.to_csv(outf, sep = "\t", index=False)

    # 9. phasing and haplotagging -- multiprocessing
    if snp_dir == 'False':
        print('** 4. phasing and haplotagging NOT selected (not specified any SNP data)')
    else:
        print('** 4. phasing and haplotagging')
        os.system('mkdir %s/phasing' %(output_directory))
        print('**** finding SNPs for phasing')
        snps_for_phasing, snps_to_keep, SNPs_data_path = findSNPs_gwas(snp_dir, bed_regions, window_for_phasing)
        print('**** start phasing                                  ', end = '\r')
        pool = multiprocessing.Pool(processes=number_threads)
        phasing_fun = partial(phase_reads_MP, reads_bam = reads_bam, snps_to_keep = snps_to_keep, output_directory = '%s/phasing' %(output_directory), SNPs_data_directory = SNPs_data_path, snp_data_ids = snp_data_ids)
        phasing_results = pool.map(phasing_fun, reads_bam)
        print('**** done with phasing                                     ')

    # 10. combine results and output files if phasing was selected
    if snp_dir != 'False':
        phasing_info = {k:v for element in phasing_results for k,v in element.items()}
        fout = open('%s/phasing/haplotags_reads.txt' %(output_directory), 'w')
        header = 'SAMPLE\tREAD_ID\tHAPLOTYPE\n'
        fout.write(header)
        for sample in phasing_info.keys():
            for read in phasing_info[sample]:
                if read == 'NA':
                    to_write = '%s\t%s\t%s\n' %(sample, 'NA', 'NA')
                else:
                    to_write = '%s\t%s\t%s\n' %(sample, read, phasing_info[sample][read])
                fout.write(to_write)
        fout.close()
    
    # 11. assembly
    print('** 5. assembly')
    os.system('mkdir %s/assembly' %(output_directory))
    strategy = AsmStrategy(assembly_type, reads_fasta, '%s/assembly' %(output_directory))
    # decide how many assembly in parallel to run: if there's only 1 genome to do, give all the cpus that were submitted (also true for alignment)
    if len(reads_fasta) == 1:
        threads_per_asm_aln = number_threads; parallel_assemblies = 1; parallel_alignment = 1; all_samples = list(strategy.keys())
    # otherwise (if there are multiple samples to process), if the user requested a specif number of cpu use that, otherwise use 4 cpus for assembly and alignment
    else:
        if number_threads_asm_aln == 'default':
            threads_per_asm_aln = 4; parallel_assemblies = int(number_threads / threads_per_asm_aln); parallel_alignment = parallel_assemblies
            threads_per_asm_aln = 4 if parallel_assemblies >1 else number_threads
            all_samples = list(strategy.keys())
        else:
            threads_per_asm_aln = int(number_threads_asm_aln); parallel_assemblies = int(number_threads / threads_per_asm_aln); all_samples = list(strategy.keys())
        if parallel_assemblies == 0:
            parallel_assemblies = 1; parallel_alignment = 1
    pool = multiprocessing.Pool(processes=parallel_assemblies)
    assembly_fun = partial(localAssembly_MP, strategy = strategy, out_dir = '%s/assembly' %(output_directory), ploidy = assembly_ploidy, thread = threads_per_asm_aln)
    assembly_results = pool.map(assembly_fun, all_samples)
    print('**** done with assembly                                     ')

    # 12. clean contigs
    print('** 6. clean assembled contigs')
    out_dir = cleanContigs('%s/assembly' %(output_directory))

    # 13. alignment
    print('** 7. align contigs')
    pool = multiprocessing.Pool(processes=parallel_alignment)
    align_fun = partial(alignAssembly_MP, outname_list = assembly_results, thread = threads_per_asm_aln, reference = ref_fasta)
    align_results = pool.map(align_fun, assembly_results)
    print('**** done with contig alignment                                     ')

    # 14. measure distance on contigs
    print('** 8. calculate size of the regions of interest in contigs')
    # list all files that should be processed
    haps_to_process = ['%s_haps_hg38.bam' %(x) for x in assembly_results]; prim_to_process = ['%s_p_ctg_hg38.bam' %(x) for x in assembly_results]; files_to_process = haps_to_process + prim_to_process
    pool = multiprocessing.Pool(processes=number_threads)
    measure_fun = partial(measureDistance_MP, bed = bed_regions, window = 1)
    extract_results = pool.map(measure_fun, files_to_process)
    print('**** read measurement done!                                         ')
    
    # 15. combine results together and make output
    distances = {k:v for element in extract_results for k,v in element.items()}; os.system('mkdir %s/trf_assembly' %(output_directory))

    # 16. tandem repeat finder on assembled contigs
    print('** 9. tandem repeat finder on contigs')
    all_bams = list(distances.keys())
    motif = readMotif(bed_file)
    pool = multiprocessing.Pool(processes=number_threads)
    trf_fun = partial(trf_MP, out_dir = '%s/trf_assembly' %(output_directory), motif = motif, polished = 'False', distances = distances)
    trf_results = pool.map(trf_fun, all_bams)
    df_trf = pd.concat(trf_results)
    print('**** done running TRF on assemblies                                     ')

    # 17. make output
    outf = output_directory + '/trf_assembly/measures_spanning_reads_and_trf.txt'
    df_trf.to_csv(outf, sep = "\t", index=False)

    # 18. coverage profile
    print("** 10. generating coverage profile")
    os.system('mkdir %s/coverage' %(output_directory))
    pool = multiprocessing.Pool(processes=number_threads)
    coverage_fun = partial(generateCoverageProfile_MP, bed = bed_regions, all_bams = reads_bam, window_size = window_size, step = step, output_directory = '%s/coverage' %(output_directory))
    coverage_results = pool.map(coverage_fun, reads_bam)
    print('**** done with coverage profiles                                     ')

    # 19. combine results and output files
    coverage_info = {k:v for element in coverage_results for k,v in element.items()}
    outname = open('%s/coverage/coverage_profiles.bed' %(output_directory), 'w')
    outname.write('CHROM\tSTART_POS\tEND_POS\tSAMPLE\tCOVERAGE\tREGION_ID\n')
    for s in coverage_info.keys():
        for x in coverage_info[s]:
            outname.write('%s\t%s\t%s\t%s\t%s\t%s\n' %(x[0], x[1], x[2], x[3], x[4], x[5]))
    outname.close()

    # 20. haplotyping
    print("** 11. haplotype calling and reads-spanning vs. assembly comparison")
    file_path = os.path.realpath(__file__)
    file_path = '/'.join(file_path.split('/')[:-1])
    if snp_dir == 'False':
        os.system("Rscript %s/call_haplotypes.R --reads_spanning %s/trf_reads/measures_spanning_reads_and_trf.txt --asm %s/trf_assembly/measures_spanning_reads_and_trf.txt --out %s/haplotyping --cpu %s" %(file_path, output_directory, output_directory, output_directory, output_directory, number_threads))
    else:
        os.system("Rscript %s/call_haplotypes.R --reads_spanning %s/trf_reads/measures_spanning_reads_and_trf.txt --phase %s/phasing/haplotags_reads.txt --asm %s/trf_assembly/measures_spanning_reads_and_trf.txt --out %s/haplotyping --cpu %s" %(file_path, output_directory, output_directory, output_directory, output_directory, number_threads))
elif anal_type == 'realign':
    filestorealign = getFilesToRealign(fasta_dir)
    print("** aligning fasta files")
    output_aligned = realignFasta(filestorealign, number_threads, ref_fasta)
elif anal_type == 'haplotyping':
    print('** haplotyping')
    file_path = os.path.realpath(__file__)
    file_path = '/'.join(file_path.split('/')[:-1])
    os.system('Rscript %s/call_haplotypes.R --reads_spanning %s --phase %s --asm %s --out %s --deviation %s' %(file_path, trf_file, phase_file, asm_file, output_directory, thr_mad))
    #if snp_dir == 'False':
    #    os.system("Rscript %s/call_haplotypes.R --reads_spanning %s --asm %s --out %s" %(file_path, trf_file, asm_file, output_directory))
    #else:
    #    os.system("Rscript %s/call_haplotypes.R --reads_spanning %s --phase %s --asm %s --out %s" %(file_path, trf_file, phase_file, asm_file, output_directory))
elif anal_type == 'reads_spanning_trf':
    # 1. read bed regions
    print("**** reading bed file")
    bed_regions = readBed(bed_file)

    # 2. check bam files
    print("** checking bam files")
    all_bam_files = checkBAM(bam_directory)
    
    # 3. read motif
    motif = readMotif(bed_file)
    
    # 4. extract reads mapping to regions of interest in bed -- multiprocessing
    print("\n** 1. extracting reads")
    os.system('mkdir %s/raw_reads' %(output_directory))
    pool = multiprocessing.Pool(processes=number_threads)
    extract_fun = partial(extract_sequences_nowrite_nosort, bed = bed_regions, out_dir = '%s/raw_reads' %(output_directory), window = window_size)
    ts = time.time()
    extract_results = pool.map(extract_fun, all_bam_files)
    te = time.time()
    time_extraction = te-ts
    print('**** read extraction done in %s seconds                                 ' %(round(time_extraction, 0)))
    # combine results together: this will give a list for each sample
    reads_info = [x[-2] for x in extract_results]; reads_bam = [x[0] for x in extract_results]; reads_fasta = [x[1] for x in extract_results]; reads_fasta_padding = [x[2] for x in extract_results]; reads_ids = [x[-1] for x in extract_results]
    # separately do the reference genome
    dist_reference, reads_ids_reference, fasta_ref, fasta_withPadd_ref = measureDistance_reference(bed_regions, window_size, ref_fasta, output_directory)
    # combine reads_ids into a single dictionary
    reads_ids.append(reads_ids_reference)
    reads_ids_combined = {}
    for dic in reads_ids:
        for key, value in dic.items():
            reads_ids_combined[key] = value
    # also combine the fasta files and the reads info and the bam
    reads_fasta.append(fasta_ref)
    reads_info.append(dist_reference)
    all_bam_files.append('reference')

    # 5. TRF on single-reads
    print('** 3. tandem repeat finder on the single-reads')
    trf_out_dir = '%s/trf_reads' %(output_directory); os.system('mkdir %s' %(trf_out_dir));
    pool = multiprocessing.Pool(processes=number_threads)
    trf_fun = partial(run_trf, out_dir = trf_out_dir, motif = motif, polished = 'False', distances = reads_info, reads_ids_combined = reads_ids_combined, all_bam_files = all_bam_files)
    ts = time.time()
    trf_results = pool.map(trf_fun, reads_fasta)
    te = time.time()
    time_trf = te-ts
    # combine df from different samples together
    df_trf_combined = pd.concat(trf_results)
    print('**** TRF estimation done in %s seconds                                 ' %(round(time_trf, 0)))
    os.system('rm -rf %s/trf_reads' %(output_directory))

    # 6. phasing and haplotagging -- multiprocessing
    if snp_dir == 'False':
        print('** 4. phasing and haplotagging NOT selected (not specified any SNP data)')
        combined_haplotags_df = pd.DataFrame(columns=['READ_NAME', 'HAPLOTAG'])
    else:
        print('** 4. phasing and haplotagging')
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

    # 7. combine these results and make output
    df_trf_phasing_combined = pd.merge(df_trf_combined, combined_haplotags_df, left_on = 'READ_NAME', right_on = 'READ_NAME', how = 'outer')
    outf = '%s/spanning_reads_trf_phasing.txt' %(output_directory)
    df_trf_phasing_combined.to_csv(outf, sep = "\t", index=False, na_rep='NA')

    # 11. haplotyping
    print("** 11. haplotype calling on reads-spanning and eventually phasing")
    file_path = os.path.realpath(__file__)
    file_path = '/'.join(file_path.split('/')[:-1])
    if snp_dir == 'False':
        os.system("Rscript %s/call_haplotypes_v2.R --reads_spanning %s/spanning_reads_trf_phasing.txt --out %s --cpu %s" %(file_path, output_directory, output_directory, number_threads))
    else:
        os.system("Rscript %s/call_haplotypes_v2.R --reads_spanning %s/spanning_reads_trf_phasing.txt --out %s --cpu %s" %(file_path, output_directory, output_directory, number_threads))
elif anal_type == 'assembly_trf':
    # 1. read bed regions
    print("**** reading bed file")
    bed_regions = readBed(bed_file)

    # 2. check bam files
    print("** checking bam files")
    all_bam_files = checkBAM(bam_directory)
    
    # 3. read motif
    motif = readMotif(bed_file)
    
    # 4. extract reads mapping to regions of interest in bed -- multiprocessing
    if ass_tool == 'hifiasm':
        print("\n** 1. extracting reads")
        os.system('mkdir %s/raw_reads' %(output_directory))
        pool = multiprocessing.Pool(processes=number_threads)
        extract_fun = partial(extractReadsAssembly, bed = bed_regions, out_dir = '%s/raw_reads' %(output_directory), window = window_asm, all_bams = all_bam_files)
        extract_results = pool.map(extract_fun, all_bam_files); print('**** read extraction done!                                         ')
        
        # 5. combine results together
        reads_fasta = [x for x in extract_results]

        # 6. assembly
        print('** 5. assembly')
        os.system('mkdir %s/assembly' %(output_directory))
        strategy = AsmStrategy(assembly_type, reads_fasta, '%s/assembly' %(output_directory))
        # decide how many assembly in parallel to run: if there's only 1 genome to do, give all the cpus that were submitted (also true for alignment)
        if len(reads_fasta) == 1:
            threads_per_asm_aln = number_threads; parallel_assemblies = 1; parallel_alignment = 1; all_samples = list(strategy.keys())
        # otherwise (if there are multiple samples to process), if the user requested a specif number of cpu use that, otherwise use 4 cpus for assembly and alignment
        else:
            if number_threads_asm_aln == 'default':
                threads_per_asm_aln = 4; parallel_assemblies = int(number_threads / threads_per_asm_aln); parallel_alignment = parallel_assemblies
                threads_per_asm_aln = 4 if parallel_assemblies >1 else number_threads
                all_samples = list(strategy.keys())
            else:
                threads_per_asm_aln = int(number_threads_asm_aln); parallel_assemblies = int(number_threads / threads_per_asm_aln); all_samples = list(strategy.keys())
            if parallel_assemblies == 0:
                parallel_assemblies = 1; parallel_alignment = 1 
        pool = multiprocessing.Pool(processes=parallel_assemblies)
        assembly_fun = partial(localAssembly_MP, strategy = strategy, out_dir = '%s/assembly' %(output_directory), ploidy = assembly_ploidy, thread = threads_per_asm_aln)
        assembly_results = pool.map(assembly_fun, all_samples)
        print('**** done with assembly                                     ')

        # 7. clean contigs
        print('** 6. clean assembled contigs')
        out_dir = cleanContigs('%s/assembly' %(output_directory))

        # 8. alignment
        print('** 7. align contigs')
        threads_per_aln = 4; parallel_alignment = int(number_threads / threads_per_aln)
        pool = multiprocessing.Pool(processes=parallel_alignment)
        align_fun = partial(alignAssembly_MP, outname_list = assembly_results, thread = threads_per_asm_aln, reference = ref_fasta)
        align_results = pool.map(align_fun, assembly_results)
        print('**** done with contig alignment                                     ')

        # 4. extract reads mapping to regions of interest in bed -- multiprocessing
        # list all files that should be processed
        haps_to_process = ['%s_haps_aln.bam' %(x) for x in assembly_results]; prim_to_process = ['%s_p_ctg_aln.bam' %(x) for x in assembly_results]; files_to_process = haps_to_process + prim_to_process
        print("\n** 8. extracting sequence of interest from aligned contigs")
        pool = multiprocessing.Pool(processes=number_threads)
        extract_fun = partial(extract_sequences_nowrite_nosort, bed = bed_regions, out_dir = '%s/raw_reads' %(output_directory), window = window_size)
        ts = time.time()
        extract_results = pool.map(extract_fun, files_to_process)
        te = time.time()
        time_extraction = te-ts
        print('**** sequence extraction done in %s seconds                                 ' %(round(time_extraction, 0)))
        # combine results together: this will give a list for each sample
        reads_info = [x[-2] for x in extract_results]; reads_bam = [x[0] for x in extract_results]; reads_fasta = [x[1] for x in extract_results]; reads_fasta_padding = [x[2] for x in extract_results]; reads_ids = [x[-1] for x in extract_results]
        # separately do the reference genome
        dist_reference, reads_ids_reference, fasta_ref, fasta_withPadd_ref = measureDistance_reference(bed_regions, window_size, ref_fasta, output_directory)
        # combine reads_ids into a single dictionary
        reads_ids.append(reads_ids_reference)
        reads_ids_combined = {}
        for dic in reads_ids:
            for key, value in dic.items():
                reads_ids_combined[key] = value
        # also combine the fasta files and the reads info and the bam
        reads_fasta.append(fasta_ref)
        reads_info.append(dist_reference)
        files_to_process.append('reference')

        # 5. TRF on single-reads
        print('** 3. tandem repeat finder on the single-reads')
        trf_out_dir = '%s/trf_reads' %(output_directory); os.system('mkdir %s' %(trf_out_dir));
        pool = multiprocessing.Pool(processes=number_threads)
        trf_fun = partial(run_trf, out_dir = trf_out_dir, motif = motif, polished = 'False', distances = reads_info, reads_ids_combined = reads_ids_combined, all_bam_files = files_to_process)
        ts = time.time()
        trf_results = pool.map(trf_fun, reads_fasta)
        te = time.time()
        time_trf = te-ts
        # combine df from different samples together
        df_trf_combined = pd.concat(trf_results)
        print('**** TRF estimation done in %s seconds                                 ' %(round(time_trf, 0)))
        os.system('rm -rf %s/trf_reads' %(output_directory))
        # add haplotag column
        combined_haplotags_df = pd.DataFrame(columns=['READ_NAME', 'HAPLOTAG'])
        
        # 6. combine these results and make output
        df_trf_phasing_combined = pd.merge(df_trf_combined, combined_haplotags_df, left_on = 'READ_NAME', right_on = 'READ_NAME', how = 'outer')
        outf = '%s/assembly_trf_phasing.txt' %(output_directory)
        df_trf_phasing_combined.to_csv(outf, sep = "\t", index=False, na_rep='NA')

        # 11. haplotyping
        print("** 11. haplotype calling on reads-spanning and eventually phasing")
        file_path = os.path.realpath(__file__)
        file_path = '/'.join(file_path.split('/')[:-1])
        os.system("Rscript %s/call_haplotypes_v2.R --asm %s/assembly_trf_phasing.txt --out %s --cpu %s" %(file_path, output_directory, output_directory, number_threads))
    else:
        print("\n** 1. run otter")
        os.system('mkdir %s/otter_local_asm' %(output_directory))
        # run otter for all samples and regions
        for s in all_bam_files:
            outname = s.split('/')[-1].replace('.bam', '.fa')
            cmd = '/project/holstegelab/Software/nicco/bin/otter/build/otter assemble -b %s -r %s %s > %s/otter_local_asm/%s' %(bed_file, ref_fasta, s, output_directory, outname)
            os.system(cmd)
        # collect otter sequences
        otter_files = [x.rstrip() for x in os.popen('ls %s/otter_local_asm/' %(output_directory))]
        res = {}
        for f in otter_files:
            f_open = [x.rstrip() for x in open('%s/otter_local_asm/%s' %(output_directory, f), 'r').readlines()]
            for x in f_open:
                if x.startswith('>'):
                    region = x.split()[0].replace('>', '')
                else:
                    seq = x; seq_len = len(seq)
                    # save hit at this point
                    if f in res.keys():
                        res[f].append([region, seq, seq_len])
                    else:
                        res[f] = [[region, seq, seq_len]]
        # at the end, save the output
        outf = open('%s/otter_local_asm/otter_sizes.txt' %(output_directory), 'w')
        for s in res.keys():
            for r in res[s]:
                outf.write('%s\t%s\t%s\t%s\n' %(s, r[0], r[1], r[2]))
        outf.close()

## Check whether to keep or delete temporary files
if (store_temporary == 'False' and anal_type not in ['reads_spanning_trf', 'assembly_trf', 'haplotyping', 'extract_reads', 'coverage_profile', 'complete']):
    print(cleanTemp(output_directory, anal_type))

## Final time
end_time = time.time()
total_time = end_time - start_time

## Final message 
print('\n** run complete in %s seconds! all results are correctly stored. \ngoing to sleep now \nciao!' %(round(total_time)))
