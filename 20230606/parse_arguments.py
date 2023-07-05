#!/usr/bin/env python3

# This script manages the arguments provided along with the script

# Libraries
import argparse
import os

# Define the parser
parser = argparse.ArgumentParser(description='TREAT: Tandem REpeat Haplotyping Toolkit')

# Create subparsers
subp = parser.add_subparsers(dest='cmd', description='Analysis types:')
# read-based analysis
readAnal = subp.add_parser('reads', help='Read-based analysis', description='Analysis based on reads spanning the regions of interest.')
# assembly-based analysis
asseAnal = subp.add_parser('assembly', help='Assembly-based analysis', description='Analysis based on local assembly of reads spanning the regions of interest.')

# Define the arguments for reads analysis
# required arguments
# bed file
readAnal.add_argument('-b', '--bed', required=True, help='BED file with the regions(s) to look at. Header is not required but if present, it must start with #.')
# input bam file(s)
readAnal.add_argument('-i', '--inBam', required=True, help='BAM file to be used as input. A directory can be provided, in which case all BAM files in the directory will be used.')
# output directory
readAnal.add_argument('-o', '--outDir', required=True, help='Output directory where to place outputs. If the directory exists, will add files there, otherwise the directory will be created.')
# reference genome
readAnal.add_argument('-r', '--ref', required=True, help='Path to reference genome data in FASTA format. Reference needs to be indexed.')
# optional arguments
# window around
readAnal.add_argument('-w', '--window', type = int, help = 'Integer. Will extend the regions defined in the BED file by this value upstream and downstreatm.', required = False, default = 10)
# number of threads
readAnal.add_argument('-t', '--cpu', type = int, help = 'Number of parallel threads to be used.', required = False, default = 2)
# phasing: path to snp data
readAnal.add_argument('-p', '--phasingData', type = str, help = 'The path to SNP data (in PLINK2 format).', required = False, default = 'None')
# phasing: path to mapping file between snp data and sequencing data
readAnal.add_argument('-m', '--mappingSNP', type = str, help = 'Path to a 2-column file with SNP IDs and sequencing IDs. If not provided, will assume the IDs are the same.', required = False, default = 'None')
# haplotyping: deviation
readAnal.add_argument('-d', '--HaploDev', type = float, help = 'During haplotying analysis, median absolute deviation to assign reads to the same allele.', required = False, default = 0.10)
# haplotyping: minimum supporting read number
readAnal.add_argument('-minSup', '--minimumSupport', type = int, help = 'During haplotying, minimum number of reads supporting each haplotyping.', required = False, default = 2)
# haplotyping: minimum coverage
readAnal.add_argument('-minCov', '--minimumCoverage', type = int, help = 'During haplotying, minimum number of total reads necessary for calling.', required = False, default = 5)

# Define the arguments for assembly analysis
# required arguments
# bed file
asseAnal.add_argument('-b', '--bed', required=True, help='BED file with the regions(s) to look at. Header is not required but if present, it must start with #.')
# input bam file(s)
asseAnal.add_argument('-i', '--inBam', required=True, help='BAM file to be used as input. A directory can be provided, in which case all BAM files in the directory will be used.')
# output directory
asseAnal.add_argument('-o', '--outDir', help='Output directory where to place outputs. If the directory exists, will add files there, otherwise the directory will be created.')
# reference genome
asseAnal.add_argument('-r', '--ref', required=True, help='Path to reference genome data in FASTA format. Reference needs to be indexed.')
# optional arguments
# window around
asseAnal.add_argument('-w', '--window', type = int, help = 'Integer. Will extend the regions defined in the BED file by this value upstream and downstreatm.', required = False, default = 10)
# number of threads
asseAnal.add_argument('-t', '--cpu', type = int, help = 'Number of parallel threads to be used.', required = False, default = 2)
# haplotyping: deviation
asseAnal.add_argument('-d', '--HaploDev', type = float, help = 'During haplotying analysis, median absolute deviation to assign reads to the same allele.', required = False, default = 0.10)
# haplotyping: minimum supporting read number
asseAnal.add_argument('-minSup', '--minimumSupport', type = int, help = 'During haplotying, minimum number of reads supporting each haplotyping.', required = False, default = 2)
# haplotyping: minimum coverage
asseAnal.add_argument('-minCov', '--minimumCoverage', type = int, help = 'During haplotying, minimum number of total reads necessary for calling.', required = False, default = 5)

# Parse the arguments
args = parser.parse_args()

# Flag to run or not the main script
RUN = False

# Print message
if args.cmd == 'reads':
    print('Read-based analysis selected')
    print('** Required argument:')
    print("   Input BAM file(s): ", args.inBam)
    print("   Input BED file: ", args.bed)
    print("   Output directory: ", args.outDir)
    print("   Reference genome: ", args.ref)
    print('** Optional arguments:')
    print("   Window: ", args.window)
    print("   Number of threads: ", args.cpu)
    print("   Phasing data: ", args.phasingData)
    print("   Phasing data IDs: ", args.mappingSNP)
    print("   Haplotyping deviation: ", args.HaploDev)
    print("   Minimum supporting reads: ", args.minimumSupport)
    print("   Minimum coverage: ", args.minimumCoverage)
    print("\n")
    # set flag to true
    RUN = True
    # define script to run and arguments
    script_path = 'read_based.py'
    arguments = [args.inBam, args.bed, args.outDir, args.ref, str(args.window), str(args.cpu), args.phasingData, args.mappingSNP, str(args.HaploDev), str(args.minimumSupport), str(args.minimumCoverage)]
elif args.cmd == 'assembly':
    print('Assembly-based analysis selected')
    print('** Required argument:')
    print("   Input BAM file(s): ", args.inBam)
    print("   Input BED file: ", args.bed)
    print("   Output directory: ", args.outDir)
    print("   Reference genome: ", args.ref)
    print('** Optional arguments:')
    print("   Window: ", args.window)
    print("   Window for assembly: ", args.windowAssembly)
    print("   Number of threads: ", args.cpu)
    print("   Assembly ploidy: ", args.ploidy)
    print("   Assembly software: ", args.software)
    print("   Haplotyping deviation: ", args.HaploDev)
    print("   Minimum supporting reads: ", args.minimumSupport)
    print("   Minimum coverage: ", args.minimumCoverage)
    print("\n")
    # set flag to true
    RUN = True
    # define script to run and arguments
    script_path = 'assembly.py'
    arguments = [args.inBam, args.bed, args.outDir, args.ref, str(args.window), str(args.windowAssembly), str(args.cpu), str(args.ploidy), args.software, str(args.HaploDev), str(args.minimumSupport), str(args.minimumCoverage)]
else:
    print('!! Invalid run_type. Quitting.')
    print("\n")

# If all arguments are good, run the main script
if RUN == True:
    # Run the script
    main_script = './%s %s' %(script_path, ' '.join(arguments))
    print(main_script)