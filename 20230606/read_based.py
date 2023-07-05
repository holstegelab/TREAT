#!/usr/bin/env python3

# This script manages the read-based analysis

# Libraries
import sys
import os

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
    return bed

# Check directory
def checkOutDir(out_dir):
    if out_dir[-1] == '/':
        out_dir = out_dir[:-1]
    if os.path.isdir(out_dir) == False:
        os.system('mkdir %s' %(out_dir))
        return("**** output directory not found, will create.")
    else:
        return("**** output directory found, will add outputs there.")

# Check bam file(s)
def checkBAM(bam_dir):
    if bam_dir[-1] == '/':
        bam_dir = bam_dir[:-1]
    if os.path.isdir(bam_dir) == True:              # in case a directory was submitted
        all_bams = [x.rstrip()for x in list(os.popen('ls %s/*bam' %(bam_dir)))]
        print("**** found directory with %s bam" %(len(all_bams)))
    elif os.path.isfile(bam_dir) == True:           # in case is a single bam file
        print("**** found single bam")
        all_bams = [bam_dir]
    elif ',' in bam_dir:                            # in case there is a comma-separated list of bams
        all_bams = bam_dir.split(',')
        print("**** found %s bam" %(len(all_bams)))
    return all_bams

# Main
# Read arguments
inBam_dir, bed_dir, outDir, ref, window, cpu, phasingData, mappingSNP, HaploDev, minimumSupport, minimumCoverage = sys.argv[1::]

# Start analysis
print('** Analysis started')

# 1. Check output directory
print(checkOutDir(outDir))

# 2. Read bed file
bed = readBed(bed_dir)

# 3. Check BAM files
inBam = checkBAM(inBam_dir)

