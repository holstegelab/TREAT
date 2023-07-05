#!/usr/bin/env python3

# This script manages the read-based analysis

# Libraries
import sys

# Functions

# Main
# Read arguments
inBam, bed, outDir, ref, window, cpu, phasingData, mappingSNP, HaploDev, minimumSupport, minimumCoverage = sys.argv
print(inBam, bed, outDir)