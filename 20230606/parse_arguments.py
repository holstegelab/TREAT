#!/usr/bin/env python

# This script manages the arguments provided along with the script

# Libraries
import argparse

# Define the parser
parser = argparse.ArgumentParser(description='TREAT: Tandem REpeat Haplotyping Toolkit')

# Define the positional argument for the type of run
parser.add_argument('run_type', help='Type of run', choices=['reads', 'assembly'])

# Create subparsers for each run type
subparsers = parse.add_subparsers(title='subcommands_reads')

# Define the subparser for reads analysis
parser_type1 = subparsers.add_parser('reads', help='Reads analysis: help')
# required arguments
# bed file
parser_type1.add_argument('-b', '--bed', required=True, help='BED file with the regions(s) to look at. Header is not required but if present, it must start with #.')
# input bam file(s)
parser_type1.add_argument('-i', '--inBam', required=True, help='BAM file to be used as input. A directory can be provided, in which case all BAM files in the directory will be used.')
# output directory
parser_type1.add_argument('-o', '--outDir', required=True, help='Output directory where to place outputs. If the directory exists, will add files there, otherwise the directory will be created.')
# reference genome
parser_type1.add_argument('-r', '--ref', required=True, help='Path to reference genome data in FASTA format. Reference needs to be indexed.')
# optional arguments
# window around
parser_type1.add_argument('-w', '--window', type = int, help = 'Integer. Will extend the regions defined in the BED file by this value upstream and downstreatm.', required = False, default = 10)
# number of threads
parser_type1.add_argument('-t', dest = '--cpu', type = int, help = 'Number of parallel threads to be used.', required = False, default = 2)
# phasing: path to snp data
parser_type1.add_argument('-p', '--phasingData', type = str, help = 'The path to SNP data (in PLINK2 format).', required = False, default = 'None')
# phasing: path to mapping file between snp data and sequencing data
parser_type1.add_argument('-m', '--mappingSNP', type = str, help = 'Path to a 2-column file with SNP IDs and sequencing IDs. If not provided, will assume the IDs are the same.', required = False, default = 'None')
# haplotyping: deviation
parser_type1.add_argument('-d', '--HaploDev', type = float, help = 'During haplotying analysis, median absolute deviation to assign reads to the same allele.', required = False, default = 0.10)
# haplotyping: minimum supporting read number
parser_type1.add_argument('-minSup', '--minimumSupport', type = int, help = 'During haplotying, minimum number of reads supporting each haplotyping.', required = False, default = 2)
# haplotyping: minimum coverage
parser_type1.add_argument('-minCov', '--minimumCoverage', type = int, help = 'During haplotying, minimum number of total reads necessary for calling.', required = False, default = 5)

# Define the subparser for assembly analysis
parser_type2 = subparsers.add_parser('assembly', help='Assembly analysis: help')
# required arguments
# bed file
parser_type2.add_argument('-b', '--bed', required=True, help='BED file with the regions(s) to look at. Header is not required but if present, it must start with #.')
# input bam file(s)
parser_type2.add_argument('-i', '--inBam', required=True, help='BAM file to be used as input. A directory can be provided, in which case all BAM files in the directory will be used.')
# output directory
parser_type2.add_argument('-o', '--outDir', help='Output directory where to place outputs. If the directory exists, will add files there, otherwise the directory will be created.')
# reference genome
parser_type2.add_argument('-r', '--ref', required=True, help='Path to reference genome data in FASTA format. Reference needs to be indexed.')
# optional arguments
# window around
parser_type2.add_argument('-w', '--window', type = int, help = 'Integer. Will extend the regions defined in the BED file by this value upstream and downstreatm.', required = False, default = 10)
# number of threads
parser_type2.add_argument('-t', dest = '--cpu', type = int, help = 'Number of parallel threads to be used.', required = False, default = 2)
# ploidy
parser_type2.add_argument('-p', dest = '--ploidy', type = int, help = 'Expected ploidy to be used for the local assembly procedure.', required = False, default = 2)
# assembly software
parser_type2.add_argument('-s', '--software', type = str, help = '[hifiasm / otter]: the tool to use for local assembly.', required = False, default = 'hifiasm')
# window for assembly
parser_type2.add_argument('-wAss', '--windowAssembly', type = int, help = 'Window to use to recruit reads for assembly.', required = False, default = 2500)
# haplotyping: deviation
parser_type2.add_argument('-d', '--HaploDev', type = float, help = 'During haplotying analysis, median absolute deviation to assign reads to the same allele.', required = False, default = 0.10)
# haplotyping: minimum supporting read number
parser_type2.add_argument('-minSup', '--minimumSupport', type = int, help = 'During haplotying, minimum number of reads supporting each haplotyping.', required = False, default = 2)
# haplotyping: minimum coverage
parser_type2.add_argument('-minCov', '--minimumCoverage', type = int, help = 'During haplotying, minimum number of total reads necessary for calling.', required = False, default = 5)

# Parse the arguments
args = parser.parse_args()

# Print the values of the arguments
print(f"Run type: {args.run_type}")

if args.run_type == 'reads':
    print('Reads analysis selected:')
    print('Required argument:')
    print(f"\tInput BAM file(s): {args.inBam}")
    print(f"\tInput BED file: {args.bed}")
    print(f"\tOutput directory: {args.outDir}")
    print(f"\tReference genome: {args.ref}")
    print('Optional arguments:')
    print(f"\tWindow: {args.window}")
    print(f"\tNumber of threads: {args.cpu}")
    print(f"\tPhasing data: {args.phasingData}")
    print(f"\tPhasing data IDs: {args.mappingSNP}")
    print(f"\tHaplotyping deviation: {args.HaploDev}")
    print(f"\tMinimum supporting reads: {args.minimumSupport}")
    print(f"\tMinimum coverage: {args.minimumCoverage}")
elif args.run_type == 'assembly':
    print('Assembly analysis selected:')
    print('Required argument:')
    print(f"\tInput BAM file(s): {args.inBam}")
    print(f"\tInput BED file: {args.bed}")
    print(f"\tOutput directory: {args.outDir}")
    print(f"\tReference genome: {args.ref}")
    print('Optional arguments:')
    print(f"\tWindow: {args.window}")
    print(f"\tWindow for assembly: {args.windowAssembly}")
    print(f"\tNumber of threads: {args.cpu}")
    print(f"\tAssembly ploidy: {args.ploidy}")
    print(f"\tAssembly software: {args.software}")
    print(f"\tHaplotyping deviation: {args.HaploDev}")
    print(f"\tMinimum supporting reads: {args.minimumSupport}")
    print(f"\tMinimum coverage: {args.minimumCoverage}")
