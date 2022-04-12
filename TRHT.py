#############################################################
# SET OF TOOLS THAT CAN BE USED TO MANAGE AND ANALYZE DATA  #
# FROM PACBIO OR OTHER SEQUENCING PLATFORMS.                #
#############################################################

## Libraries
import os
import sys
from functions_trht import *
import argparse
import os.path
import pysam
from Bio.Seq import Seq
from Bio import SeqIO
from random import random
import gzip
import subprocess
import json
import pandas as pd
from datetime import datetime

## Main
## Define Arguments
parser = argparse.ArgumentParser(description = 'Find information about a specific region/tandem repeat')
parser.add_argument('--bed', dest = 'bed_dir', type = str, help = '.bed file containing the region(s) to look. Header is not required but if present, it must starts with #.', required = False, default = 'None')
parser.add_argument('--analysis-type', dest = 'analysis_type', type = str, help = 'Type of analysis to perform [extract_reads / measure / trf / assembly / extract_snps / annotate_snps / extract_annotate / genotype_snps_pacbio / phase_reads / coverage_profile / extract_raw_reads]. See docs for further information.', required = True)
parser.add_argument('--variant-file', dest = 'variant_file', type = str, help = 'If the analysis_type is annotate_snps or genotype_snps_pacbio, please provide here the path to the file including the SNPs to annotate/extract.', required = False, default = 'None')
parser.add_argument('--bam-dir', dest = 'bam_dir', type = str, help = 'Directory of bam file(s). If a directory is provided, will use all .bam in the directory. If a single .bam file is provided, will use that file.', required = False, default = 'None')
parser.add_argument('--out-dir', dest = 'out_dir', type = str, help = 'Directory where to place output files. If the directory exists, will place files in, otherwise will create the folder and place results in.', required = True)
parser.add_argument('--store-temp', dest = 'store_temp', type = str, help = 'Boolean (True/False). If True, will store all temporary .bam and .fasta files.', required = False, default = 'False')
parser.add_argument('--window', dest = 'window', type = int, help = 'Integer. Will use this value to take reads surrounding the provided .bed file.', required = False, default = 1000)
parser.add_argument('--assembly-type', dest = 'ass_type', type = str, help = 'Type of local assembly to perform. By default, each .bam will result in an assembly. If you prefer to use multiple .bam files for an assembly, please submit a file with no header and 2 columns: the first column should report, for each line, a comma-separated list of .bam files to combine in the assembly. The second column, for each line, should report the output prefix of the final assembly for each group.', required = False, default = 'asm_per_bam')
parser.add_argument('--assembly-ploidy', dest = 'ass_ploidy', type = int, help = 'Ploidy to be used for the local assembly procedure. Default value is 2 (for diploid organisms).', required = False, default = 2)
parser.add_argument('--thread', dest = 'thread', type = int, help = 'Number of parallel threads to be used during assembly.', required = False, default = 4)
parser.add_argument('--polish', dest = 'polish', type = str, help = 'Boolean (True/False). If True, reads containing the region of interest will also be polished.', required = False, default = 'False')
parser.add_argument('--snp-data', dest = 'snp_dir', type = str, help = 'If assembly is selected, please add here the path to SNP data (in PLINK2 format). This information is necessary for phasing.', required = False, default = 'False')
parser.add_argument('--snp-data-ids', dest = 'snp_data_ids', type = str, help = 'Please submit here a 2-column file with GWAS ID and ID in sequencing data. If not provided, will assume the IDs are the same.', required = False, default = 'False')
parser.add_argument('--coverage-step', dest = 'step', type = int, help = 'Number of nt based on which the region(s) of interest will be split to calculate coverage.', required = False, default = 500)
parser.add_argument('--reads-ids', dest = 'target_reads', type = int, help = 'If analysis type was extract_raw_reads, please provide here a tab-separated file with two columns: the first column should contain the .bam file to extract reads from, the second column should contain the reads ID to extract, one read ID per line.', required = False, default = 'False')

args = parser.parse_args()
# Check arguments
if (args.analysis_type not in ['annotate_snps', 'extract_snps', 'extract_annotate', 'extract_raw_reads']) and (args.bam_dir == 'None'):
    parser.error('!! You should provide at least a .bam file if analysis_type is --> %s' %(args.analysis_type))
elif (args.analysis_type not in ['genotype_snps_pacbio', 'annotate_snps', 'extract_raw_reads']) and (args.bed_dir == 'None'):
    parser.error('!! You should provide at least a .bed file if analysis_type is --> %s' %(args.analysis_type))
if (args.analysis_type == 'extract_raw_reads' and args.target_reads == 'False'):
    parser.error('!! Your should provide either a text file containing the IDs of interest or a comma-separated list of IDs')

# Print arguments
print("\n** Tandem Repeat Haplotyping Toolkit (TRHT) **\n")
print("********************\n** Your settings:")
print("** bed file --> %s" %(args.bed_dir))
print("** bam file(s) --> %s" %(args.bam_dir))
print("** output folder --> %s" %(args.out_dir))
print("** analysis type --> %s" %(args.analysis_type))
if args.analysis_type == 'phase_reads':
    print("** phasing using SNPs in --> %s" %(args.snp_dir))
    print("** mapping IDs between long-read and SNPs using --> %s" %(args.snp_data_ids))
print('** polishing --> %s' %(args.polish))
if args.analysis_type == 'assembly':
    print("** assembly type --> %s" %(args.ass_type))
if args.analysis_type == 'coverage_profile':
    print("** step --> %s" %(args.step))
if args.analysis_type == 'extract_raw_reads':
    print('** read IDs to extract --> %s' %(args.target_reads))
print("** intermediate files --> %s" %(args.store_temp))
print('** variant file --> %s' %(args.variant_file))
print("** window used --> %s\n********************\n" %(args.window))

# Store arguments
bed_file, anal_type, var_file, bam_directory, output_directory, store_temporary, window_size, assembly_type, assembly_ploidy, number_threads, polishing, snp_dir, snp_data_ids, step, target_reads = args.bed_dir, args.analysis_type, args.variant_file, args.bam_dir, args.out_dir, args.store_temp, args.window, args.ass_type, args.ass_ploidy, args.thread, args.polish, args.snp_dir, args.snp_data_ids, args.step, args.target_reads
# Store arguments (for debugging only)
#bed_file, anal_type, var_file, bam_directory, output_directory, store_temporary, window_size, assembly_type, assembly_ploidy, number_threads, polishing, snp_dir, snp_data_ids, step = '../bed_files/RFC1.bed', '', '', '/project/holstegelab/Share/pacbio/blood_brain_child/hifi/c7_child_merged_hifi.bam', 'coverage_profile/blood_brain_child', 'True', 250000, '', 2, 4, 'True', '/project/holstegelab/Share/gwas_array/Imputed_data_100plus/chrAll_imputed_100plus.pvar', '/project/holstegelab/Share/nicco/workspaces/20211013_target_approach/automatic_pipeline/RFC1/phenotypes/map_pabio_gwas_allSamples.txt', 500

## Check whether output directory exists otherwise create it
print("** checking directories")
print(checkOutDir(output_directory))

## Extract reads mapping to the location of interest
if anal_type == 'extract_reads':
    print("** reading bed file")
    bed = readBed(bed_file)
    print("** checking bam files")
    all_bams = checkBAM(bam_directory)
    reads_bam, reads_fasta = extractReads(bed, all_bams, output_directory, window_size)
    store_temporary = 'True'
elif anal_type == 'measure':    
    print("** reading bed file")
    bed = readBed(bed_file)
    print("** checking bam files")
    all_bams = checkBAM(bam_directory)
    reads_bam, reads_fasta = extractReads(bed, all_bams, output_directory, window_size)
    print('\n** calculate size of the regions of interest')
    distances = measureDistance(bed, reads_bam, 25, output_directory)
    if polishing == 'True':
        print('** polishing reads of interest now!')
        distances_polished = polishReads(bed, distances, output_directory)          # check from here!
        os.system('rm ' + output_directory + '/measures_spanning_reads.txt')
elif anal_type == 'trf':
    print("** reading bed file")
    bed = readBed(bed_file)
    print("** checking bam files")
    all_bams = checkBAM(bam_directory)
    motif = readMotif(bed_file)
    reads_bam, reads_fasta = extractReads(bed, all_bams, output_directory, window_size)
    print('\n** calculate size of the regions of interest')
    distances = measureDistance(bed, reads_bam, 25, output_directory)
    if polishing == 'True':
        print('** polishing reads of interest now!')
        distances_polished = polishReads(bed, distances, output_directory)   
    print('** tandem repeat finder')
    trf_info = trf(distances, output_directory, motif, polished = 'False')
    if polishing == 'True':
        print('** tandem repeat finder on corrected reads')
        trf_info_polished = trf(distances_polished, output_directory, motif, polished = 'True')
        # checked until here
elif anal_type == 'assembly':
    print("** reading bed file")
    bed = readBed(bed_file)
    print("** checking bam files")
    all_bams = checkBAM(bam_directory)
    reads_bam, reads_fasta = extractReads(bed, all_bams, output_directory, window_size)
    strategy = AsmStrategy(assembly_type, reads_fasta, output_directory)
    print('** running local assembly')
    outname_list = localAssembly(strategy, output_directory, assembly_ploidy, number_threads)
    print('** cleaning assembled contigs')
    out_dir = cleanContigs(output_directory)
    print('\n** aligning assembled contigs')
    alignAssembly(outname_list, number_threads)
elif anal_type == 'phase_reads':
    print("** reading bed file")
    bed = readBed(bed_file)
    print("** checking bam files")
    all_bams = checkBAM(bam_directory)
    reads_bam, reads_fasta = extractReads(bed, all_bams, output_directory, window_size)
    print('\n** find SNPs for phasing')
    snps_for_phasing, snps_to_keep, SNPs_data_path = findSNPs_gwas(snp_dir, bed, window_size)
    print('** phasing')
    phasing_info = phase_reads(reads_bam, snps_to_keep, output_directory, SNPs_data_path, snp_data_ids)
elif anal_type == 'extract_snps':
    print("** reading bed file")
    bed = readBed(args.bed_dir)
    snps_topmed = extractSNPs(bed, args.window)
    writeSNPs(args.out_dir, snps_topmed)
elif anal_type == 'annotate_snps':
    input_variants = readSNPs(args.variant_file)
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
    coverage_profile = generateCoverageProfile(bed, all_bams, window_size, step, output_directory)   
elif anal_type == 'extract_raw_reads':
    print('** reading target reads and bam file(s)')
    target = readTarget(target_reads)

if (store_temporary == 'False' and analysis_type != 'coverage_profile'):
    print(cleanTemp(args.out_dir, args.analysis_type))
print('\n** run complete! all results are correctly stored. \ngoing to sleep now \nciao!')
