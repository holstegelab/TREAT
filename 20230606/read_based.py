#!/usr/bin/env python3

# This script manages the read-based analysis

# Libraries
import sys

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

# Main
# Read arguments
inBam, bed_dir, outDir, ref, window, cpu, phasingData, mappingSNP, HaploDev, minimumSupport, minimumCoverage = sys.argv[1::]

# Start analysis
print('** Analysis started')
# Read bed file
bed = readBed(bed_dir)

