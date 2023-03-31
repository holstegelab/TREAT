## Packages
import os
import sys
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
from contextlib import contextmanager
from inspect import getsourcefile
from os.path import abspath

## Functions
# Function to suppress standard output
@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:  
            yield
        finally:
            sys.stdout = old_stdout

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
    print('**** found %s regions in %s chromosomes' %(count_reg, len(bed)))
    return bed

# Read bed file with motif
def readMotif(bed_dir):
    motif = {}
    with open(bed_dir) as finp:
        for line in finp:
            if line.startswith('#'):
                pass
            else:
                line = line.rstrip().split()
                if len(line) == 4:
                    chrom, start, end, motif_type = line[0:4]
                else:
                    chrom, start, end = line[0:3]
                    motif_type = "NA"
                region_id = chrom + ':' + start + '-' + end
                motif[region_id] = motif_type
    return motif

# Check directory
def checkOutDir(out_dir, analysis_type):
    if analysis_type == 'realign':
        return("**** output directory will be the same as input")
    else:
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

# Extract sequence of interest precisely through CIGAR string
def findPositionOfInterest(cigar, positions_of_interest, positions_of_interest_end, positions_of_interest_with_padding, positions_of_interest_with_padding_end):
   # define counter of the reference and the raw sequences
   counter_ref = 0
   counter_raw = 0
   counter_raw_padd = 0
   # define positions of interest
   pos_interest = 0
   pos_interest_padd = 0
   pos_interest_end = 0
   pos_interest_padd_end = 0
   # make a list of 1 cigar element per position
   cigar_per_base = [x for cse in cigar for x in [cse[0]] * cse[1]]
   # Then loop on this list
   for x in cigar_per_base:
      # Parse cigar types: 
      if (x == 7) or (x == 0):   # 7 --> =
         counter_raw += 1
         counter_ref += 1
         counter_raw_padd += 1
      elif x == 8:   # 8 --> X
         counter_raw += 1
         counter_ref += 1
         counter_raw_padd += 1
      elif x == 1:   # 1 --> I
         counter_raw +=1
         counter_raw_padd += 1
      elif x == 2:   # 2 --> D
         counter_ref += 1
      elif x == 4:  # 4 --> S
         counter_raw += 1
         counter_raw_padd += 1
      elif x == 5:  # 5 --> H
         counter_raw += 1
         counter_raw_padd += 1
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
   return pos_interest, pos_interest_end, pos_interest_padd, pos_interest_padd_end
# Alternative version of the above to use a while loop: the advantage should be if the repeat is short, it doesn't need to go through all cigar elements
def findPositionOfInterestWhile(cigar, positions_of_interest, positions_of_interest_end, positions_of_interest_with_padding, positions_of_interest_with_padding_end):
   # define counter of the reference and the raw sequences
   counter_ref = 0
   counter_raw = 0
   counter_raw_padd = 0
   # define positions of interest
   pos_interest = 0
   pos_interest_padd = 0
   pos_interest_end = 0
   pos_interest_padd_end = 0
   # make a list of 1 cigar element per position
   cigar_per_base = [x for cse in cigar for x in [cse[0]] * cse[1]]
   # Then loop on this list
   i = 0; run = True
   while (run == True) and (i < len(cigar_per_base)):
      x = cigar_per_base[i]
      # Parse cigar types: 
      if (x == 7) or (x == 0):   # 7 --> =
         counter_raw += 1
         counter_ref += 1
         counter_raw_padd += 1
      elif x == 8:   # 8 --> X
         counter_raw += 1
         counter_ref += 1
         counter_raw_padd += 1
      elif x == 1:   # 1 --> I
         counter_raw +=1
         counter_raw_padd += 1
      elif x == 2:   # 2 --> D
         counter_ref += 1
      elif x == 4:  # 4 --> S
         counter_raw += 1
         counter_raw_padd += 1
      elif x == 5:  # 5 --> H
         counter_raw += 1
         counter_raw_padd += 1
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

# Extract reads mapping the the location of interest
def extractReads(all_bams, bed, out_dir, window):
    all_out_bam = []
    all_out_fasta = []
    for bam in all_bams:
        print("**** %s/%s bam done" %(all_bams.index(bam) + 1, len(all_bams)), end = '\r')
        fasta_seqs = {}
        outname = out_dir + '/' + bam.split('/')[-1][:-4] + '__rawReads.bam'
        outname_fasta = out_dir + '/' + bam.split('/')[-1][:-4] + '__rawReads.fasta'
        inBam = pysam.AlignmentFile(bam, 'rb', check_sq=False)
        outBam = pysam.AlignmentFile(outname, "wb", template=inBam)
        for chrom in bed.keys():
            for region in bed[chrom]:
                start, end = int(region[0]) - window, int(region[1]) + window
                for read in inBam.fetch(chrom, start, end):
                    outBam.write(read)                              # write to new bam file
                    sequence_interest = str(read.query_sequence)    # save fasta sequences
                    read_name = read.query_name
                    fasta_seqs[read_name] = sequence_interest
        outBam.close()
        inBam.close()
        outf = open(outname_fasta, 'w')
        for read in fasta_seqs.keys():
            outf.write('>%s\n%s\n' %(read, fasta_seqs[read]))
        outf.close()
        all_out_bam.append(outname)
        all_out_fasta.append(outname_fasta)
        os.system('samtools sort %s > %s' %(outname, outname + '_sorted'))
        os.system('mv %s %s' %(outname + '_sorted', outname))
        os.system('samtools index %s' %(outname))
    return all_out_bam, all_out_fasta

# Extract reads mapping the the location of interest
def extractReads_MP(bam, bed, out_dir, window, all_bams):
    # dictionary of fasta sequences
    fasta_seqs = {}
    # set object to remember seen reads
    read_ids_list = set()
    # define output names of bam and fasta files
    outname = out_dir + '/' + os.path.basename(bam)[:-4] + '__rawReads.bam'
    outname_fasta = out_dir + '/' + os.path.basename(bam)[:-4] + '__rawReads.fasta'
    # open bam files
    with open(outname_fasta, 'w') as outf:
        with pysam.AlignmentFile(bam, 'rb', check_sq=False) as inBam:
            with pysam.AlignmentFile(outname, "wb", template=inBam) as outBam:
                for chrom in bed.keys():
                    for region in bed[chrom]:
                        start, end = int(region[0]) - window, int(region[1]) + window
                        for read in inBam.fetch(chrom, start, end):
                            read_name = read.query_name
                            if read_name not in read_ids_list:
                                read_ids_list.add(read_name)
                                outf.write('>%s\n%s\n' %(read_name, str(read.query_sequence)))
                                outBam.write(read)                              # write to new bam file
    print("**** done with %s          " %(bam.split('/')[-1]), end = '\r')
    return outname, outname_fasta

# Extract reads mapping to the location of interest for assembly
def extractReadsAssembly(bam, bed, out_dir, window, all_bams):
    # set object to remember seen reads
    read_ids_list = set()
    # define output names of fasta files
    outname_fasta = out_dir + '/' + os.path.basename(bam)[:-4] + '__rawReads.fasta'
    # open bam files
    with open(outname_fasta, 'w') as outf:
        with pysam.AlignmentFile(bam, 'rb', check_sq=False) as inBam:
            for chrom in bed.keys():
                for region in bed[chrom]:
                    start, end = int(region[0]) - window, int(region[1]) + window
                    for read in inBam.fetch(chrom, start, end):
                        read_name = read.query_name
                        if read_name not in read_ids_list:
                            read_ids_list.add(read_name)
                            outf.write('>%s\n%s\n' %(read_name, str(read.query_sequence)))
    print("**** done with %s          " %(bam.split('/')[-1]), end = '\r')
    return outname_fasta

# Polish reads using Alex's script
def polishReads(bed, distances, out_dir):
    # first find the unique region ids
    unique_regions = []
    for chrom in bed.keys():
        for region in bed[chrom]:
            if region[-1] not in unique_regions:
                unique_regions.append(region[-1])
    distances_corrected = {}
    # main loop across samples
    for sample in distances.keys():
        if sample != 'reference':
            # create input fasta for polishing script
            tmp_fname = out_dir + '/' + sample[:-4] + '_inputPolish.fasta'
            tmp_f = open(tmp_fname, 'w')
            sample_reads = distances[sample]
            # order the reads based on region id and sequence length (longer sequences first)
            for region in unique_regions:
                ordered_reads = {}
                for read in sample_reads:
                    if read[0] == region:
                        key = region + ' ' + read[1]
                        ordered_reads[key] = read[-2]
                order_of_reads = sorted(ordered_reads, key=ordered_reads.get, reverse=True)
                for read_in_order in order_of_reads:
                    read_name = read_in_order.split(' ')[1]
                    for read in sample_reads:
                        if read[1] == read_name:
                            to_write = '>' + read_in_order + '\n' + read[-3] + '\n'
                            tmp_f.write(to_write)
            tmp_f.close()
            corrected_reads = os.popen('urchin correct -m 3 ' + tmp_fname).read().split('\n')
            corrected_reads = corrected_reads[:-1]
            corrected_reads_id = [x for x in corrected_reads if x.startswith('>')]
            corrected_reads_seq = [x for x in corrected_reads if not x.startswith('>')]
            for i in range(len(corrected_reads_id)):
                tmp = corrected_reads_id[i].replace('>', '')
                region, read_name = tmp.split(' ')
                sequence = corrected_reads_seq[i]
                for read in sample_reads:
                    if read[1] == read_name:
                        info = [region, read_name, read[2], read[3], read[4], read[5], read[6], read[7], sequence, len(read[5]) - len(sequence)]
                        if sample in distances_corrected.keys():
                            distances_corrected[sample].append(info)
                        else:
                            distances_corrected[sample] = [info]        
        else:
            sample_reads = distances[sample]
            for x in sample_reads:
                tmp = x[:]
                tmp.append('NA')
                tmp.append('NA')
                if sample in distances_corrected.keys():
                    distances_corrected[sample].append(tmp)
                else:
                    distances_corrected[sample] = [tmp]
    # finally make 1 output table
    outf = open(out_dir + '/measures_spanning_reads_after_correction.txt', 'w')
    outf.write('REGION\tSAMPLE_NAME\tREAD_NAME\tPASSES\tREAD_QUALITY\tMAPPING_CONSENSUS\tSEQUENCE_WITH_WINDOW\tLENGTH_SEQUENCE\tPADDING_SIZE\tPOLISHED_SEQUENCE\tLENGTH_POLISHED_SEQUENCE\tSIZE_DIFF\n')
    for x in distances_corrected.keys():
        for read in distances_corrected[x]:
            outf.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(read[0], x, read[1], read[2], read[3], read[4], read[5], read[6], read[7], read[8], len(read[8]), read[9]))
    outf.close()
    return(distances_corrected)

# Given two position, measure the size in the reads
def measureDistance(bed, reads_bam, window, out_dir):
    distances = {}
    for bam in reads_bam:
        bam_name = bam.split('/')[-1]
        inBam = pysam.AlignmentFile(bam, 'rb', check_sq=False)
        for chrom in bed.keys():
            for region in bed[chrom]:
                region_id = chrom + ':' + region[0] + '-' + region[1]
                start, end = int(region[0]) - window, int(region[1]) + window
                for read in inBam.fetch(chrom, start, end):
                    ref_start = int(read.reference_start)           # take the reference start position
                    ref_end = int(read.reference_end)               # take the reference end position
                    if (ref_end != "NA") and (int(ref_start) <= start) and (int(ref_end) >= end):
                        start_poi_dist = start - int(ref_start)         # calculate distance in the aligned sequence between start and positions of interest
                        end_poi_dist = end - int(ref_start)
                        cigar = read.cigartuples                        # take cigar string
                        poi_start = findPositionOfInterest(cigar, start_poi_dist)                             # find start position of interest
                        poi_end = findPositionOfInterest(cigar, end_poi_dist)                                 # find end position of interest
                        padding_before = str(read.query_sequence)[0 : poi_start]                                # find padding sequence before sequence of interest
                        padding_after = str(read.query_sequence)[poi_end :]                                     # find padding sequence after sequence of interest   
                        sequence_interest = str(read.query_sequence)[poi_start : poi_end]                       # find sequence of interest
                        read_name = str(read.query_name)
                        info = read.tags                                                                        # read tags
                        np, rq, mc = 'NA', 'NA', 'NA'
                        for x in info:
                            if x[0] == "np":
                                np = "NP:%s" %(x[1])
                            elif x[0] == "rq":
                                rq = "RQ:%s" %(x[1])
                            elif x[0] == "mc":
                                mc = "MC:%s" %(x[1])
                        sec_aln = str(read.is_secondary)
                        sup_aln = str(read.is_supplementary)
                        if bam_name in distances.keys():
                            distances[bam_name].append([region_id, read_name, np, rq, mc, sequence_interest, len(sequence_interest), window])
                        else:
                            distances[bam_name] = [[region_id, read_name, np, rq, mc, sequence_interest, len(sequence_interest), window]]
        inBam.close()
    for chrom in bed.keys():
        for region in bed[chrom]:
            region_id = chrom + ':' + region[0] + '-' + region[1]
            start, end = int(region[0]) - window, int(region[1]) + window
            sequence_in_reference = [x.rstrip() for x in list(os.popen('samtools faidx /project/holstegelab/Software/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa %s:%s-%s' %(chrom, start, end)))]
            seq_merged = ''.join(sequence_in_reference[1:])
            if 'reference' in distances.keys():
                distances['reference'].append([region_id, 'NA', 'NA', 'NA', 'NA', seq_merged, len(seq_merged), window])
            else:
                distances['reference'] = [[region_id, 'NA', 'NA', 'NA', 'NA', seq_merged, len(seq_merged), window]]
    # finally make 1 output table
    outf = open(out_dir + '/measures_spanning_reads.txt', 'w')
    outf.write('REGION\tSAMPLE_NAME\tREAD_NAME\tPASSES\tREAD_QUALITY\tMAPPING_CONSENSUS\tSEQUENCE_WITH_WINDOW\tLENGTH_SEQUENCE\tPADDING_SIZE\n')
    for x in distances.keys():
        for read in distances[x]:
            outf.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(read[0], x, read[1], read[2], read[3], read[4], read[5], read[6], read[7]))
    outf.close()
    return distances

# Given two position, measure the size in the reads
def measureDistance_MP(reads_bam, bed, window):
    distances = {}
    bam_name = reads_bam.split('/')[-1]
    inBam = pysam.AlignmentFile(reads_bam, 'rb', check_sq=False)
    for chrom in bed.keys():
        for region in bed[chrom]:
            region_id = chrom + ':' + region[0] + '-' + region[1]
            start, end = int(region[0]) - window, int(region[1]) + window
            for read in inBam.fetch(chrom, start, end):
                ref_start = int(read.reference_start)           # take the reference start position
                ref_end = int(read.reference_end)               # take the reference end position
                if (ref_end != "NA") and (int(ref_start) <= start) and (int(ref_end) >= end):
                    start_poi_dist_with_padding = start - int(ref_start)         # calculate distance in the aligned sequence between start and positions of interest: WITH PADDING
                    start_poi_dist = (start + window) - int(ref_start)         # calculate distance in the aligned sequence between start and positions of interest: WITHOUT PADDING
                    end_poi_dist_with_padding = end - int(ref_start)
                    end_poi_dist = (end - window) - int(ref_start)
                    cigar = read.cigartuples                        # take cigar string
                    poi_start_with_padding = findPositionOfInterest(cigar, start_poi_dist_with_padding)                             # find start position of interest WITH PADDING
                    poi_start = findPositionOfInterest(cigar, start_poi_dist)                             # find start position of interest WITHOUT PADDING
                    poi_end_with_padding = findPositionOfInterest(cigar, end_poi_dist_with_padding)                                 # find end position of interest
                    poi_end = findPositionOfInterest(cigar, end_poi_dist)                                 # find end position of interest
                    padding_before = str(read.query_sequence)[0 : poi_start]                                # find padding sequence before sequence of interest
                    padding_after = str(read.query_sequence)[poi_end :]                                     # find padding sequence after sequence of interest   
                    sequence_interest = str(read.query_sequence)[poi_start : poi_end]                       # find sequence of interest WITHOUT PADDING
                    sequence_interest_with_padding = str(read.query_sequence)[poi_start_with_padding : poi_end_with_padding]                       # find sequence of interest WITHOUT PADDING
                    read_name = str(read.query_name)
                    info = read.tags                                                                        # read tags
                    np, rq, mc = 'NA', 'NA', 'NA'
                    for x in info:
                        if x[0] == "np":
                            np = "NP:%s" %(x[1])
                        elif x[0] == "rq":
                            rq = "RQ:%s" %(x[1])
                        elif x[0] == "mc":
                            mc = "MC:%s" %(x[1])
                    sec_aln = str(read.is_secondary)
                    sup_aln = str(read.is_supplementary)
                    if bam_name in distances.keys():
                        distances[bam_name].append([region_id, read_name, np, rq, mc, sequence_interest, sequence_interest_with_padding, len(sequence_interest), len(sequence_interest_with_padding), window])
                    else:
                        distances[bam_name] = [[region_id, read_name, np, rq, mc, sequence_interest, sequence_interest_with_padding, len(sequence_interest), len(sequence_interest_with_padding), window]]
    inBam.close()
    print('**** done measuring %s                   ' %(bam_name), end = '\r')
    return distances

# Function to measure distance
def measureDist(read, chrom, start, end, window):
   # region identifier
   region_id = '%s:%s-%s' %(chrom, start + window, end - window)
   # read identifier
   read_id = str(read.query_name)
   # take reference start position and end position
   ref_start = int(read.reference_start)
   ref_end = int(read.reference_end)
   # check if read is spanning the region of interest 
   if (ref_end != "NA") and (int(ref_start) <= start) and (int(ref_end) >= end):
      # calculate distance between start/end position and position of interest, with and without the paddings on the sides
      start_poi_dist_with_padding = start - int(ref_start)
      end_poi_dist_with_padding = end - int(ref_start)
      start_poi_dist = (start + window) - int(ref_start)
      end_poi_dist = (end - window) - int(ref_start)
      # take cigar string
      cigar = read.cigartuples
      # find start/end positions of interest for the sequence with and witout paddings
      pos_interest, pos_interest_end, pos_interest_padd, pos_interest_padd_end = findPositionOfInterestWhile(cigar, start_poi_dist, end_poi_dist, start_poi_dist_with_padding, end_poi_dist_with_padding)
      # then extract the sequence in the region of interest, with and without padding
      sequence_interest = str(read.query_sequence)[pos_interest : pos_interest_end]
      sequence_interest_len = len(sequence_interest)
      sequence_interest_with_padding = str(read.query_sequence)[pos_interest_padd : pos_interest_padd_end]
      sequence_interest_with_padding_len = len(sequence_interest_with_padding)
      # also look into tags: number of passes, read quality, mapping consensus
      info = read.tags
      np, rq, mc = 'NA', 'NA', 'NA'
      for x in info:
         if x[0] == "np":
            np = "NP:%s" %(x[1])
         elif x[0] == "rq":
            rq = "RQ:%s" %(x[1])
         elif x[0] == "mc":
            mc = "MC:%s" %(x[1])
   else:
      np, rq, mc, sequence_interest, sequence_interest_with_padding, sequence_interest_len, sequence_interest_with_padding_len = "NA", "NA", "NA", "NA", "NA", "NA", "NA"
   return [read_id, region_id, np, rq, mc, sequence_interest, sequence_interest_with_padding, sequence_interest_len, sequence_interest_with_padding_len, window]

# Measure the distance in the reference genome
def measureDistance_reference(bed, window, ref, output_directory):
    distances = []
    reads_ids = {'reference' : []}
    for chrom in bed.keys():
        for region in bed[chrom]:
            region_id = chrom + ':' + region[0] + '-' + region[1]
            start_with_padding, end_with_padding = int(region[0]) - window, int(region[1]) + window           # coordinates with padding
            start, end = int(region[0]), int(region[1])           # coordinates without padding
            sequence_in_reference_with_padding = [x.rstrip() for x in list(os.popen('samtools faidx %s %s:%s-%s' %(ref, chrom, start_with_padding, end_with_padding)))]        # sequence with padding
            sequence_in_reference = [x.rstrip() for x in list(os.popen('samtools faidx %s %s:%s-%s' %(ref, chrom, start, end)))]        # sequence without padding
            seq_merged = ''.join(sequence_in_reference[1:])
            seq_merged_with_padding = ''.join(sequence_in_reference_with_padding[1:])
            distances.append(['reference', region_id, 'NA', 'NA', 'NA', seq_merged, seq_merged_with_padding, len(seq_merged), len(seq_merged_with_padding), window])
            reads_ids['reference'].append(region_id)
    # then we write the fasta
    outfasta = '%s/raw_reads/reference__rawReads.fasta' %(output_directory)
    outf = open(outfasta, 'w')
    outfasta_withPad = '%s/raw_reads/reference__rawReads_withPaddings.fasta' %(output_directory)
    outf_withPad = open(outfasta_withPad, 'w')
    for x in distances:
        read_id, region_id, np, rq, mc, sequence_interest_len, sequence_interest, sequence_interest_with_padding_len, sequence_interest_with_padding = x[0], x[1], x[2], x[3], x[4], x[-3], x[5], x[-2], x[-4]
        outf.write('>%s;%s;%s;%s;%s;%s\n%s\n' %(read_id, region_id, np, rq, mc, sequence_interest_len, sequence_interest))
        outf_withPad.write('>%s;%s;%s;%s;%s;%s;%s\n%s\n' %(read_id, region_id, np, rq, mc, window, sequence_interest_with_padding_len, sequence_interest_with_padding))
    return distances, reads_ids, outfasta, outfasta_withPad

# Run TRF given a sequence
def trf(distances, out_dir, motif, polished):
    trf_info = {}
    all_bams = list(distances.keys())
    for bam in distances.keys():                    # loop on samples, that is bam files
        print("**** %s/%s bam done" %(all_bams.index(bam) + 1, len(all_bams)), end = '\r')
        for read in distances[bam]:                 # then loop on each read, that is, every read on every region
            seq = read[-2] if polished == 'True' else read[-3]
            if bam == 'reference' and polished == 'True':
                seq = read[-5]
            elif bam == 'reference' and polished != 'True':
                seq = read[-3]
            region = read[0]         # extract sequence and region name
            motif_type = motif[region]              # get the correspondin motif
            exact_motifs = motif_type.split(',') if ',' in motif_type else [motif_type]         # find motif and all permutations of that
            perm_motifs, rev_motifs, rev_perm_motifs = [], [], []       
            for m in exact_motifs:
                tmp_perm = [m[x:] + m[:x] for x in range(len(m))]
                tmp_perm.remove(m)
                reverse = str(Seq(m).reverse_complement())
                perm_reverse = [reverse[x:] + reverse[:x] for x in range(len(reverse))]
                perm_reverse.remove(reverse)
                perm_motifs, rev_motifs, rev_perm_motifs = perm_motifs + tmp_perm, rev_motifs + [reverse], rev_perm_motifs + perm_reverse
            tmp_name = out_dir + '/tmp_trf_' + str(random()).replace('.', '') + '.fasta'            # then create fasta file for running trf
            tmp_out = open(tmp_name, 'w')
            tmp_out.write('>%s\n%s\n' %(bam, seq))
            tmp_out.close()
            # run tandem repeat finder specifying as maximum motif length the known one
            cmd = 'trf4.10.0-rc.2.linux64.exe ' + tmp_name + ' 2 7 7 80 10 50 %s -ngs -h' %(len(exact_motifs[0]))
            trf = [x for x in os.popen(cmd).read().split('\n') if x != '']
            motif_len = 'no_motif'
            if len(trf) == 0:
                add_trial = 0
                while len(trf) == 0 and add_trial <3:
                    # if no trf results were found, first try to lower the minimum score to 40 for TRF match
                    if add_trial == 0:
                        cmd = 'trf4.10.0-rc.2.linux64.exe ' + tmp_name + ' 2 7 7 80 10 40 %s -ngs -h' %(len(exact_motifs[0]))
                        trf = [x for x in os.popen(cmd).read().split('\n') if x != '']
                        add_trial += 1
                        motif_len = 'same_length'
                    # if no trf results were found again, first try to lower the minimum score to 30 for TRF match
                    elif add_trial == 1:
                        cmd = 'trf4.10.0-rc.2.linux64.exe ' + tmp_name + ' 2 7 7 80 10 30 %s -ngs -h' %(len(exact_motifs[0]))
                        trf = [x for x in os.popen(cmd).read().split('\n') if x != '']
                        add_trial += 1
                        motif_len = 'same_length'
                    # if again no matches, try to match any motif with higher score
                    elif add_trial == 2:
                        cmd = 'trf4.10.0-rc.2.linux64.exe ' + tmp_name + ' 2 7 7 80 10 50 200 -ngs -h'
                        trf = [x for x in os.popen(cmd).read().split('\n') if x != '']
                        add_trial += 1
                        motif_len = 'different_length'
                    # if again no matches, try to match any motif with looser score
                    else:
                        cmd = 'trf4.10.0-rc.2.linux64.exe ' + tmp_name + ' 2 7 7 80 10 30 200 -ngs -h'
                        trf = [x for x in os.popen(cmd).read().split('\n') if x != '']
                        add_trial += 1
                        motif_len = 'different_length'
            else:
                motif_len = 'same_length'
            # look whether there are hits for the corresponding motif -- if so, save results into dictionary
            if len(trf) > 1:
                for match in trf:
                    if '@' in match:
                        pass
                    elif motif_len == 'same_length':
                        trf_motif = match.split()[13]
                        if trf_motif in exact_motifs:
                            motif_match = "exact_motif"
                        elif trf_motif in perm_motifs:
                            motif_match = "perm_motif"
                        elif trf_motif in rev_motifs:
                            motif_match = "rev_motif"
                        elif trf_motif in rev_perm_motifs:
                            motif_match = "perm_rev_motif"
                        else:
                            motif_match = "different_motif"
                        read.append([match, motif_match])
                    elif motif_len == 'different_length':
                        read.append([match, 'different_length'])
            else:
                read.append(['No matches'])
            if bam in trf_info.keys():
                trf_info[bam].append(read)
            else:
                trf_info[bam] = [read]
            # finally remove temporary file
            os.system('rm ' + tmp_name)          
    # finally make 1 output table
    if polished == 'True':
        outf = open(out_dir + '/measures_spanning_reads_and_trf_polished.txt', 'w')
        outf.write('REGION\tSAMPLE_NAME\tREAD_NAME\tPASSES\tREAD_QUALITY\tMAPPING_CONSENSUS\tSEQUENCE_WITH_WINDOW\tLENGTH_SEQUENCE\tPADDING_SIZE\tPOLISHED_SEQUENCE\tLENGTH_POLISHED\tDIFF_WITH_ORIGINAL\tSTART_TRF\tEND_TRF\tLENGTH_MOTIF_TRF\tCOPIES_TRF\tPC_MATCH_TRF\tPC_INDEL_TRF\tMOTIF_TRF\tPADDING_BEFORE\tSEQUENCE_TRF\tPADDING_AFTER\tMATCH_TYPE\n')
    else:
        outf = open(out_dir + '/measures_spanning_reads_and_trf.txt', 'w')
        outf.write('REGION\tSAMPLE_NAME\tREAD_NAME\tPASSES\tREAD_QUALITY\tMAPPING_CONSENSUS\tSEQUENCE_WITH_WINDOW\tLENGTH_SEQUENCE\tPADDING_SIZE\tSTART_TRF\tEND_TRF\tLENGTH_MOTIF_TRF\tCOPIES_TRF\tPC_MATCH_TRF\tPC_INDEL_TRF\tMOTIF_TRF\tPADDING_BEFORE\tSEQUENCE_TRF\tPADDING_AFTER\tMATCH_TYPE\n')
    for x in trf_info.keys():
        for read in trf_info[x]:
            trf_res = read[8::] if polished != 'True' else read[10::]
            for trf_match in trf_res:
                if isinstance(trf_match, list):
                    if trf_match[0] == 'No matches':
                        start_rep, end_rep, mot_len, cp, cons_size, pc_match, pc_indel, aln_score, pc_a, pc_c, pc_g, pc_t, entropy, motif, sequence_trf, padding_before, padding_after, match_type = ['NA']*17 + ['No matches']
                    else:
                        start_rep, end_rep, mot_len, cp, cons_size, pc_match, pc_indel, aln_score, pc_a, pc_c, pc_g, pc_t, entropy, motif, sequence_trf, padding_before, padding_after, match_type = trf_match[0].split(' ') + [trf_match[1]]
                    if polished == 'True':
                        outf.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(read[0], x, read[1], read[2], read[3], read[4], read[5], read[6], read[7], read[8], len(read[8]), read[9], start_rep, end_rep, mot_len, cp, pc_match, pc_indel, motif, padding_before, sequence_trf, padding_after, match_type))
                    else:
                        outf.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(read[0], x, read[1], read[2], read[3], read[4], read[5], read[6], read[7], start_rep, end_rep, mot_len, cp, pc_match, pc_indel, motif, padding_before, sequence_trf, padding_after, match_type))
                else:
                    if polished == 'True':
                        outf.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(read[0], x, read[1], read[2], read[3], read[4], read[5], read[6], read[7], read[8], len(read[8]), read[9], 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'))
                    else:
                        outf.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(read[0], x, read[1], read[2], read[3], read[4], read[5], read[6], read[7], 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'))
    outf.close()
    #os.system('rm %s/measures_spanning_reads.txt' %(out_dir))
    return trf_info

# Run TRF given a sequence -- MP without adjustments
def trf_MP_old(bam, out_dir, motif, polished, distances):
    trf_info = {}
    for read in distances[bam]:                 # then loop on each read, that is, every read on every region
        seq = read[-2] if polished == 'True' else read[-3]
        if bam == 'reference' and polished == 'True':
            seq = read[-5]
        elif bam == 'reference' and polished != 'True':
            seq = read[-3]
        region = read[0]         # extract sequence and region name
        motif_type = motif[region]              # get the correspondin motif
        exact_motifs = motif_type.split(',') if ',' in motif_type else [motif_type]         # find motif and all permutations of that
        perm_motifs, rev_motifs, rev_perm_motifs = [], [], []       
        for m in exact_motifs:
            tmp_perm = [m[x:] + m[:x] for x in range(len(m))]
            tmp_perm.remove(m)
            reverse = str(Seq(m).reverse_complement())
            perm_reverse = [reverse[x:] + reverse[:x] for x in range(len(reverse))]
            perm_reverse.remove(reverse)
            perm_motifs, rev_motifs, rev_perm_motifs = perm_motifs + tmp_perm, rev_motifs + [reverse], rev_perm_motifs + perm_reverse
        tmp_name = out_dir + '/tmp_trf_' + str(random()).replace('.', '') + '.fasta'            # then create fasta file for running trf
        tmp_out = open(tmp_name, 'w')
        tmp_out.write('>%s\n%s\n' %(bam, seq))
        tmp_out.close()
        # run tandem repeat finder specifying as maximum motif length the known one
        cmd = 'trf4.10.0-rc.2.linux64.exe ' + tmp_name + ' 2 7 7 80 10 50 %s -ngs -h' %(len(exact_motifs[0]))
        trf = [x for x in os.popen(cmd).read().split('\n') if x != '']
        motif_len = 'no_motif'
        if len(trf) == 0:
            add_trial = 0
            while len(trf) == 0 and add_trial <3:
                # if no trf results were found, first try to lower the minimum score to 40 for TRF match
                if add_trial == 0:
                    cmd = 'trf4.10.0-rc.2.linux64.exe ' + tmp_name + ' 2 7 7 80 10 40 %s -ngs -h' %(len(exact_motifs[0]))
                    trf = [x for x in os.popen(cmd).read().split('\n') if x != '']
                    add_trial += 1
                    motif_len = 'same_length'
                # if no trf results were found again, first try to lower the minimum score to 30 for TRF match
                elif add_trial == 1:
                    cmd = 'trf4.10.0-rc.2.linux64.exe ' + tmp_name + ' 2 7 7 80 10 30 %s -ngs -h' %(len(exact_motifs[0]))
                    trf = [x for x in os.popen(cmd).read().split('\n') if x != '']
                    add_trial += 1
                    motif_len = 'same_length'
                # if again no matches, try to match any motif with higher score
                elif add_trial == 2:
                    cmd = 'trf4.10.0-rc.2.linux64.exe ' + tmp_name + ' 2 7 7 80 10 50 200 -ngs -h'
                    trf = [x for x in os.popen(cmd).read().split('\n') if x != '']
                    add_trial += 1
                    motif_len = 'different_length'
                # if again no matches, try to match any motif with looser score
                else:
                    cmd = 'trf4.10.0-rc.2.linux64.exe ' + tmp_name + ' 2 7 7 80 10 30 200 -ngs -h'
                    trf = [x for x in os.popen(cmd).read().split('\n') if x != '']
                    add_trial += 1
                    motif_len = 'different_length'
        else:
            motif_len = 'same_length'
        # look whether there are hits for the corresponding motif -- if so, save results into dictionary
        if len(trf) > 1:
            for match in trf:
                if '@' in match:
                    pass
                elif motif_len == 'same_length':
                    trf_motif = match.split()[13]
                    if trf_motif in exact_motifs:
                        motif_match = "exact_motif"
                    elif trf_motif in perm_motifs:
                        motif_match = "perm_motif"
                    elif trf_motif in rev_motifs:
                        motif_match = "rev_motif"
                    elif trf_motif in rev_perm_motifs:
                        motif_match = "perm_rev_motif"
                    else:
                        motif_match = "different_motif"
                    read.append([match, motif_match])
                elif motif_len == 'different_length':
                    read.append([match, 'different_length'])
        else:
            read.append(['No matches'])
        if bam in trf_info.keys():
            trf_info[bam].append(read)
        else:
            trf_info[bam] = [read]
        # finally remove temporary file
        os.system('rm ' + tmp_name)
    print('**** done TRF on %s                                       ' %(bam), end = '\r')        
    return trf_info

# Run TRF given a sequence
def run_trf(fasta, out_dir, motif, polished, distances, reads_ids_combined, all_bam_files):
   # then run tandem repeat finder
   cmd = 'trf4.10.0-rc.2.linux64.exe %s 2 7 7 80 10 50 200 -ngs -h' %(fasta)
   trf = [x for x in os.popen(cmd).read().split('\n') if x != '']
   # loop on trf results and save them into a list of lists
   x = 0; trf_matches = []
   read_found = []
   while x < len(trf):
      # check if the line is the header of an entry
      if trf[x].startswith('@'):
         # if so, save the corresponding information
         read_id, region, passes, qual, cons, seq_size = trf[x].split(';')
         read_found.append(read_id.replace('@', ''))
         x += 1
         while x < len(trf) and not trf[x].startswith('@'):
            motif_type = motif[region]              # get the correspondin motif
            tmp_trf_match = [read_id.replace('@', '') + '_' + region, motif_type] + trf[x].split()
            trf_matches.append(tmp_trf_match)
            x += 1
   # finally create pandas df and assign column names
   if len(trf_matches) == 0:
      trf_matches = [['NA' for i in range(19)]] 
   df = pd.DataFrame(trf_matches)
   df.columns = ['ID', 'EXPECTED_MOTIF', 'START_TRF', 'END_TRF', 'LENGTH_MOTIF_TRF', 'COPIES_TRF', 'TRF_CONSENSUS_SIZE', 'TRF_PERC_MATCH', 'TRF_PERC_INDEL', 'TRF_SCORE', 'TRF_A_PERC', 'TRF_C_PERC', 'TRF_G_PERC', 'TRF_T_PERC', 'TRF_ENTROPY', 'TRF_MOTIF', 'TRF_REPEAT_SEQUENCE', 'TRF_PADDING_BEFORE', 'TRF_PADDING_AFTER']
   # finally, we need to add the reads where trf didn't find any motif
   sample_interest = os.path.basename(fasta).replace('__rawReads.fasta', '')
   all_samples = [os.path.basename(x).replace('.bam', '') for x in all_bam_files]
   sample_interest_index = all_samples.index(sample_interest)
   distances_sample = distances[sample_interest_index]
   # convert distances to dataframe
   distances_sample_df = pd.DataFrame(distances_sample)
   distances_sample_df.columns = ['READ_NAME', 'REGION', 'PASSES', 'READ_QUALITY', 'MAPPING_CONSENSUS', 'SEQUENCE_FOR_TRF', 'SEQUENCE_WITH_PADDINGS', 'LEN_SEQUENCE_FOR_TRF', 'LEN_SEQUENCE_WITH_PADDINGS', 'WINDOW']
   distances_sample_df['ID'] = distances_sample_df['READ_NAME'].str.cat(distances_sample_df['REGION'], sep='_')
   # add sample name in a column
   distances_sample_df['SAMPLE_NAME'] = sample_interest
   # merge trf dataframe and reads dataframes
   complete_df = pd.merge(distances_sample_df, df, left_on = 'ID', right_on = 'ID', how = 'outer')
   print("**** done with %s                               " %(sample_interest), end = '\r')
   return complete_df

# Read strategy file for assembly
def AsmStrategy(ass_type, reads_fasta, out_dir):
    strategy = {}
    if ass_type != 'asm_per_bam':
        finp = open(ass_type).readlines()
        for line in finp:
            samples_to_combine, asm_name = line.rstrip().split('\t')
            samples_to_combine = samples_to_combine.replace(' ', '').replace('.bam', '').replace('.fasta', '').replace('.fa', '').split(',')
            samples_to_combine = [out_dir + '/' + x + '__rawReads.fasta' for x in samples_to_combine]
            strategy[asm_name] = samples_to_combine
        print('**** file describing assembly with %s groups found' %(len(strategy)))
    else:
        print('**** using default settings for assembly, that is, one per .bam')
        for x in reads_fasta:
            x_name = os.path.basename(x).replace('__rawReads.fasta', '')
            strategy[x_name] = [x]
    return strategy

# Local assembly given the strategy and output directory for the fasta files
def localAssembly(strategy, out_dir, ploidy, thread):
    outnames_list = []
    for group_name in strategy.keys():
        inp_fasta = ' '.join(strategy[group_name])
        outname = '%s/%s_asm' %(out_dir, group_name)
        cmd_assembly = 'hifiasm -o %s -t %s --n-hap %s -n 2 -r 3 -N 200 %s >/dev/null 2>&1' %(outname, thread, ploidy, inp_fasta)
        try:
            os.system(cmd_assembly)
            subprocess.run("gfatools gfa2fa %s.bp.hap1.p_ctg.gfa > %s_haps.fasta" %(outname, outname), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)          # convert gfa of hap 1 to fasta
            subprocess.run("gfatools gfa2fa %s.bp.hap2.p_ctg.gfa >> %s_haps.fasta" %(outname, outname), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)          # convert gfa of hap 2 to fasta
            subprocess.run("mv %s.bp.p_ctg.gfa %s_p_ctg.gfa" %(outname, outname), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)    # rename primary contigs
            subprocess.run("gfatools gfa2fa %s_p_ctg.gfa >> %s_p_ctg.fasta" %(outname, outname), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)          # convert gfa of hap 2 to fasta
            subprocess.run("rm %s.*" %(outname), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)                 # clean non-necessary files
            print('**** %s/%s assemblies done' %(list(strategy.keys()).index(group_name) + 1, len(list(strategy.keys()))), end = '\r')
            outnames_list.append(outname)
        except:
            pass
    return outnames_list

# Local assembly given the strategy and output directory for the fasta files
def localAssembly_MP(group_name, strategy, out_dir, ploidy, thread):
    inp_fasta = ' '.join(strategy[group_name])
    outname = '%s/%s__Assembly' %(out_dir, group_name)
    cmd_assembly = 'hifiasm -o %s -t %s --n-hap %s -n 2 -r 3 -N 200 %s >/dev/null 2>&1' %(outname, thread, ploidy, inp_fasta)
    print('**** assembly done with hifiasm using %s threads. Command line is: %s\n' %(thread, cmd_assembly))    
    #cmd_assembly = 'hifiasm -o %s -t %s --n-hap %s -n 2 -r 3 -N 200 %s' %(outname, thread, ploidy, inp_fasta)
    try:
        os.system(cmd_assembly)
        subprocess.run("gfatools gfa2fa %s.bp.hap1.p_ctg.gfa > %s_haps.fasta" %(outname, outname), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)          # convert gfa of hap 1 to fasta
        subprocess.run("gfatools gfa2fa %s.bp.hap2.p_ctg.gfa >> %s_haps.fasta" %(outname, outname), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)          # convert gfa of hap 2 to fasta
        subprocess.run("mv %s.bp.p_ctg.gfa %s_p_ctg.gfa" %(outname, outname), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)    # rename primary contigs
        subprocess.run("gfatools gfa2fa %s_p_ctg.gfa >> %s_p_ctg.fasta" %(outname, outname), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)          # convert gfa of hap 2 to fasta
        subprocess.run("rm %s.*" %(outname), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)                 # clean non-necessary files
    except:
        pass
    print('**** done assembly for %s                                   ' %(group_name), end = '\r')
    return outname

# Align assembly
def alignAssembly(outname_list, thread, reference):
    if reference == 'chm13':
        ref_hifi = '/project/holstegelab/Share/asalazar/data/chm13/assembly/v2_0/chm13v2.0_hifi.mmi'
        outname_prefix = '_haps_chm13.bam'
        outname_primary_prefix = '_p_ctg_chm13.bam'
    else:
        ref_hifi = '/project/holstegelab/Share/pacbio/resources/h38_ccs.mmi'
        outname_prefix = '_haps_hg38.bam'
        outname_primary_prefix = '_p_ctg_hg38.bam'
    for asm in outname_list:
        # haplotype-aware fastas
        asm_name = asm + '_haps.fasta'
        outname = asm + outname_prefix
        cmd = "pbmm2 align --preset CCS %s -j %s --log-level FATAL --sort %s %s" %(ref_hifi, thread, asm_name, outname)       # command for alignment
        os.system(cmd)
        # primary contigs as well
        asm_name = asm + '_p_ctg.fasta'
        outname_primary = asm + outname_primary_prefix
        cmd = "pbmm2 align --preset CCS %s -j %s --log-level FATAL --sort %s %s" %(ref_hifi, thread, asm_name, outname_primary)       # command for alignment
        os.system(cmd)
        print('**** %s/%s alignment done' %(outname_list.index(asm) + 1, len(outname_list)), end = '\r')
    return(print('\n** alignment done'))

# Align assembly
def alignAssembly_MP(asm, outname_list, thread, reference):
    ref_hifi = reference
    outname_prefix = '_haps_aln.bam'
    outname_primary_prefix = '_p_ctg_aln.bam'
    # haplotype-aware fastas
    asm_name = asm + '_haps.fasta'
    outname = asm + outname_prefix
    #cmd = "pbmm2 align --preset CCS %s -j %s --log-level FATAL --sort %s %s" %(ref_hifi, thread, asm_name, outname)       # command for alignment
    cmd = "minimap2 -aYx asm10 %s -t %s %s | samtools sort - -@ %s -O bam -o %s" %(ref_hifi, thread, asm_name, thread, outname)       # command for alignment
    os.system(cmd)
    # then index
    cmd = 'samtools index %s' %(outname)
    os.system(cmd)
    # primary contigs as well
    asm_name = asm + '_p_ctg.fasta'
    outname_primary = asm + outname_primary_prefix
    #cmd = "pbmm2 align --preset CCS %s -j %s --log-level FATAL --sort %s %s" %(ref_hifi, thread, asm_name, outname_primary)       # command for alignment
    cmd = "minimap2 -ax asm10 %s -t %s %s | samtools sort - -@ %s -O bam -o %s" %(ref_hifi, thread, asm_name, thread, outname_primary)       # command for alignment
    os.system(cmd)
    cmd = 'samtools index %s' %(outname_primary)
    os.system(cmd)
    print('**** done alignment for %s                                            ' %(asm.split('/')[-1]), end = '\r')
    return(outname)

# Clean temporary files
def cleanTemp(out_dir, analysis_type):
    if analysis_type == 'assembly':
        os.system('rm %s/*__rawReads.bam*' %(out_dir))
    elif analysis_type in ['extract_snps', 'annotate_snps', 'extract_annotate', 'realign']:
        pass
    else:
        os.system('rm %s/*__rawReads*bam' %(out_dir))
        os.system('rm %s/*__rawReads*bai' %(out_dir))
        os.system('rm %s/*__rawReads*fasta' %(out_dir))
    return('** intermediate files are removed')

# Extract SNPs from TOPMED of the region of interest
def extractSNPs(bed, window):
    snps_topmed = {}
    for chrom in bed.keys():
        print('**** finding snps in chromosome --> %s' %(chrom))
        for region in bed[chrom]:
            start_pos, end_pos, region_id = int(region[0]), int(region[1]), region[-1]
            plink_file = '/project/holstegelab/Share/gwas_array/TOPMED_BRAVO/PASS.Variants%s.BRAVO_TOPMed_Freeze_8.tab.gz' %(chrom)
            with gzip.open(plink_file) as finp:
                for line in finp:
                    if line.startswith(b'#CHROM'):
                        pass
                    else:
                        line = line.rstrip().split()
                        snp_pos, snp_id, snp_ref, snp_alt, snp_af = int(line[1]), line[2], line[3], line[4], float(line[-1])
                        if snp_pos >= (start_pos - window) and snp_pos <= (end_pos + window):
                            if chrom in snps_topmed.keys():
                                snps_topmed[chrom].append([snp_pos, snp_id, snp_ref, snp_alt, snp_af])
                            else:
                                snps_topmed[chrom] = [[snp_pos, snp_id, snp_ref, snp_alt, snp_af]]
                        elif snp_pos > (end_pos + window):
                            break  
    return snps_topmed 

# Save SNP dataset
def writeSNPs(out_dir, snps_topmed):
    fout = open('%s/SNPs_in_region.tab' %(out_dir), 'w')
    fout.write('#CHROM\tPOS\tID\tREF\tALT\tAF\n')
    for chrom in snps_topmed.keys():
        for snp in snps_topmed[chrom]:
            pos, snp_id, ref, alt, af = snp[0], '%s:%s' %(chrom, str(snp[0])), snp[2].decode('utf-8'), snp[3].decode('utf-8'), snp[-1]
            fout.write('%s\t%s\t%s\t%s\t%s\t%s\n' %(chrom, pos, snp_id, ref, alt, af))
    fout.close()
    return('** variant set correctly saved')

# Save SNP-annotation dataset
def writeSNPsAnnotation(out_dir, snps_topmed_coding):
    fout = open('%s/SNPs_in_region_annotated.tab' %(out_dir), 'w')
    fout.write('#CHROM\tPOS\tSNP_ID\tREF\tALT\tANNOTATION_TYPE\tCONSEQUENCE_DETAIL\tCONSEQUENCE\tGENE_NAME\tPHRED\n')
    for chrom in snps_topmed_coding.keys():
        for snp in snps_topmed_coding[chrom]:
            pos, snp_id, ref, alt, annotype, consdetail, conseq, gene, phred = snp[1], '%s:%s' %(chrom, snp[1]), snp[2], snp[3], snp[4], snp[5], snp[6], snp[7], snp[8]
            fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(chrom, pos, snp_id, ref, alt, annotype, consdetail, conseq, gene, phred))
    fout.close()
    return('** variant set with annotation correctly saved')

# Annotate SNPs with CADD
def annotateSNPs(snps_topmed, out_dir, input_type):
    annotatedSNPs = {}
    for chrom in snps_topmed.keys():
        for snp in snps_topmed[chrom]:
            print("**** %s/%s snps annotated" %(snps_topmed[chrom].index(snp) + 1, len(snps_topmed[chrom])), end = '\r')
            snp_id = '%s:%s' %(chrom, str(snp[0])) if input_type != 'annotate_snps' else '%s:%s' %(chrom, str(snp))
            res_tmp = subprocess.call('curl -i https://cadd.gs.washington.edu/api/v1.0/%s/%s | tail -n +8 > %s/tmp_cadd.json' %('GRCh38-v1.6_inclAnno', snp_id, out_dir), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            res = open('%s/tmp_cadd.json' %(out_dir), 'r')
            os.system('rm %s/tmp_cadd.json' %(out_dir))
            res_df = pd.DataFrame.from_dict(json.load(res))
            if res_df.shape[0] > 0:
                res_sb = res_df[["Chrom", "Pos", "Ref", "Alt", "AnnoType", "ConsDetail", "Consequence", "GeneName", "PHRED"]]
                res_sb = res_sb[res_sb['Consequence'].isin(['NON_SYNONYMOUS', "SYNONYMOUS", "STOP_GAINED", "MISSENSE"])]
                res_sb = res_sb.sort_values(by=['PHRED'], ascending=False)
                if res_sb.shape[0] > 0:
                    record = list(res_sb.iloc[0])
                    if chrom in annotatedSNPs.keys():
                        annotatedSNPs[chrom].append(record)
                    else:
                        annotatedSNPs[chrom] = [record]
    return annotatedSNPs

# Read SNP set
def readSNPs(variant_file):
    if os.path.isfile(variant_file) == True:
        input_variants = {}
        variant_info = {}
        with open(variant_file) as finp:
            for line in finp:
                if line.startswith('#'):
                    pass
                else:
                    line = line.rstrip().split()
                    chrom, pos = line[0], line[1]
                    snp_info = ';'.join(line[2:]) if len(line) > 2 else 'NA'
                    if chrom in input_variants.keys():
                        input_variants[chrom].append(pos)
                        variant_info[chrom].append(snp_info)
                    else:
                        input_variants[chrom] = [pos]
                        variant_info[chrom] = [snp_info]
    else:
        print('!! please provide a valid file including SNPs\n\n')
    return input_variants, variant_info

# Extract SNPs from .bam file
def extractSNPsFromBam(input_variants, bam_files, variant_info):
    pacbio_genotypes = {}
    for bam in bam_files:
        print("**** processing snp of bam %s/%s" %(bam_files.index(bam) + 1, len(bam_files)), end = '\r')
        bam_name = bam.split('/')[-1]
        for chrom in input_variants.keys():
            all_snps_chrom = [int(x) for x in input_variants[chrom]]
            all_snps_info = list(variant_info[chrom])
            n = 300000
            indices = [i + 1 for (x, y, i) in zip(all_snps_chrom, all_snps_chrom[1:], range(len(all_snps_chrom))) if n < abs(x - y)]
            chunks = [all_snps_chrom[start:end] for start, end in zip([0] + indices, indices + [len(all_snps_chrom)])]
            chunks_info = [all_snps_info[start:end] for start, end in zip([0] + indices, indices + [len(all_snps_info)])]
            for i in range(len(chunks)):
                start_pos, end_pos = chunks[i][0], chunks[i][-1]
                tmp = os.popen("pysamstats -c %s -s %s -e %s -u -t variation --fasta /project/holstegelab/Software/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa %s" %(chrom, start_pos, end_pos, bam)).read().split("\n")
                info = [x.split('\t') for x in tmp[1:]]
                info = [x for x in info if len(x) > 1]
                counter_info = 0
                for j in range(len(info)):
                    pos, ref, coverage, a_count, c_count, t_count, g_count  = info[j][1], info[j][2], info[j][3], info[j][13], info[j][15], info[j][17], info[j][19]
                    if int(pos) in all_snps_chrom:
                        snp_type = list(filter(lambda x:pos in x, all_snps_info))[0]
                        snp_id = '%s:%s' %(chrom, pos)
                        snp_info = [snp_id, ref, coverage, 'A:%s' %(a_count), 'C:%s' %(c_count), 'G:%s' %(g_count), 'T:%s' %(t_count), snp_type]
                        counter_info += 1
                        if bam_name in pacbio_genotypes.keys():
                            pacbio_genotypes[bam_name].append(snp_info)
                        else:
                            pacbio_genotypes[bam_name] = [snp_info]
    return pacbio_genotypes

# Call genotypes of pacbio SNPs
def callGenotypes(pacbio_genotypes, min_coverage, maf_alleles):
    called_genotypes = {}
    for sample in pacbio_genotypes.keys():
        for snp in pacbio_genotypes[sample]:
            snp_id, ref, coverage, a_count, c_count, g_count, t_count, snp_type = snp
            if int(coverage) >= min_coverage:
                all_alleles = []
                all_freqs = []
                for allele in [a_count, c_count, g_count, t_count]:
                    if int(allele.split(':')[-1]) > 0:
                        frq = int(allele.split(':')[-1]) / int(coverage)
                        if frq >= maf_alleles:
                            all_alleles.append(allele.split(':')[0])
                            all_freqs.append(str(int(allele.split(':')[-1]) / int(coverage)))
                genotype, filter, freqs, genotype_type = '%s/%s' %(all_alleles[0], all_alleles[0]) if len(all_alleles) == 1 else '/'.join(all_alleles), 'PASS', ';'.join(all_freqs), 'homozygous' if len(all_alleles) == 1 else 'heterozygous'
            else:
                genotype, filter, freqs, genotype_type = 'NA', 'Not_enough_coverage', 'NA', 'NA'
            snp_info = [snp_id, ref, coverage, a_count, c_count, g_count, t_count, genotype, filter, freqs, genotype_type, snp_type]
            if sample in called_genotypes.keys():
                called_genotypes[sample].append(snp_info)
            else:
                called_genotypes[sample] = [snp_info]
    return called_genotypes

# Create output file containing pacbio SNPs
def writeGenotypes(out_dir, calledGenotypes):
    fout = open('%s/SNPs_genotypes.tab' %(out_dir), 'w')
    fout.write('SAMPLE\tSNP_ID\tREF\tCOVERAGE\tINFO_COUNTS\tGENOTYPE\tFILTER\tFREQUENCIES\tGENOTYPE_TYPE\tSNP_INFORMATION\n')
    for sample in calledGenotypes.keys():
        for snp in calledGenotypes[sample]:
            fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(sample, snp[0], snp[1], snp[2], ';'.join(snp[3:7]), snp[7], snp[8], snp[9], snp[10], snp[11]))
    fout.close()
    return('** output created')

# Clean assembly -- remove the same contigs
def cleanContigs(out_dir):
    # list the haplotype-aware assemblies in the output folder
    flist = os.popen('ls ' + out_dir + '/*__Assembly*fasta').read().split('\n'); flist = flist[:-1]
    # loop on the files
    for f in flist:
        # read fasta file
        fasta_sequences = SeqIO.parse(open(f),'fasta')
        # store unique fasta sequences into a list
        unique_sequences, unique_ids = [], []
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            if sequence not in unique_sequences:
                unique_sequences.append(sequence)
                unique_ids.append(name)
        # overwrite the old fasta with the new sequences
        with open(f, 'w') as out_file:
            for i in range(len(unique_ids)):
                to_write = '>' + unique_ids[i] + '\n' + unique_sequences[i] + '\n'
                out_file.write(to_write)
        out_file.close()
    return(out_dir)

# Function to find SNPs from GWAS data in the regions of interest given data of interest
def findSNPs_gwas(SNPs_data_directory, bed, window):
    # read GWAS data depending on data_type
    snps_interest = {}
    all_regions = []
    # if a folder was submitted, then look for the chrAll* file
    if os.path.isdir(SNPs_data_directory):
        inpf_path = os.popen('ls ' + SNPs_data_directory + '/chrAll*pvar').readlines()[0].replace('\n', '')
        inpf = open(inpf_path).readlines()
    else:
        inpf_path = SNPs_data_directory
        inpf = open(SNPs_data_directory).readlines()
    snps_to_keep = []
    chroms_interest = list(bed.keys())
    i = 1
    while i<len(inpf):
        if inpf[i].startswith('#'):
            pass
        else:
            line = inpf[i].rstrip().split()
            chrom_snp, pos, snpid, ref, alt = line[0:5]
            if 'chr' + chrom_snp in chroms_interest:
                for interval in bed['chr' + chrom_snp]:
                    start, end, region_id = int(interval[0]), int(interval[1]), interval[-1]
                    if int(pos) >= (start - window) and int(pos) <= (end + window):
                        snps_to_keep.append(snpid)
                        if region_id in snps_interest.keys():
                            snps_interest[region_id].append([chrom_snp, pos, snpid, ref, alt])
                        else:
                            snps_interest[region_id] = [[chrom_snp, pos, snpid, ref, alt]]
        i += 1
    # finally check if we have results for all genes
    regions_found = list(snps_interest.keys())
    missing_regions = list(set(all_regions) - set(regions_found))
    for x in missing_regions: print('!! No SNPs found in --> %s' %(x))
    return(snps_interest, snps_to_keep, inpf_path)

# Function that uses whatshap to assign a haplotype to reads given a vcf file to take snps -- this uses whatshap haplotag
def phase_reads(reads_bam, snps_to_keep, output_directory, SNPs_data_directory, snp_data_ids):
    map_ids = pd.read_csv(snp_data_ids, sep = " |\t", engine = 'python') if snp_data_ids != 'False' else False
    ref_path = '/project/holstegelab/Software/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa'
    phasing_info = {}
    for f in reads_bam:
        fname = f.split('/')[-1].replace('_step1.bam', '.bam')
        id_gwas = list(map_ids['ID_GWAS'][map_ids['ID_PACBIO'] == fname]) if isinstance(map_ids, pd.DataFrame) == True else [f.split('/')[-1].split('_')[0]]          # find gwas id for the corresponding sample
        # throw error in case no match is found here
        if id_gwas == []:
            parser.error('!! Error while mapping GWAS and Long-read IDs. Please ensure the sample ID are the same or provide a file as described in docs!')
            phasing_info = []
        else:
            # prepare files for phasing
            tmp_snps = open('%s/tmp_snps_interest.txt' %(output_directory), 'w')
            for x in snps_to_keep: to_write = '%s\n' %(x); tmp_snps.write(to_write)
            tmp_snps.close()
            tmp_samples = open('%s/tmp_samples_interest.txt' %(output_directory), 'w')
            header = 'IID\n'; tmp_samples.write(header)
            for x in id_gwas: to_write = '%s\n' %(x); tmp_samples.write(to_write)
            tmp_samples.close()
            # make subset of snps and sample of interest -- one sample at the time
            print('**** %s: extract snps...' %(fname), end = '\r')
            os.system('plink2 --pfile %s --extract %s/tmp_snps_interest.txt --keep %s/tmp_samples_interest.txt --recode vcf --out %s/tmp_genotypes >/dev/null 2>&1' %(SNPs_data_directory[:-5], output_directory, output_directory, output_directory))
            # check if file was created, otherwise skip
            if os.path.isfile('%s/tmp_genotypes.vcf' %(output_directory)):
                # add chr notation for chromosome to vcf
                os.system('bcftools annotate --rename-chrs /project/holstegelab/Software/nicco/bin/TRHT/test_data/chr_name_conv.txt %s/tmp_genotypes.vcf | bgzip > %s/tmp_genotypes.vcf.gz' %(output_directory, output_directory))
                os.system('tabix %s/tmp_genotypes.vcf.gz' %(output_directory))
                # create a phased vcf with whatshap
                print('**** %s: phasing snps...   ' %(fname), end = '\r')
                os.system('whatshap phase -o %s --reference=%s %s/tmp_genotypes.vcf.gz %s --ignore-read-groups --internal-downsampling 5 >/dev/null 2>&1' %(f.replace('.bam', '_phased.vcf.gz'), ref_path, output_directory, f))
                # for haplotag command, need a phased vcf
                os.system('tabix %s' %(f.replace('.bam', '_phased.vcf.gz')))
                print('**** %s: tagging reads...   ' %(fname), end = '\r')
                os.system('whatshap haplotag -o %s/tmp_haplotag.bam --reference=%s %s %s --ignore-read-groups >/dev/null 2>&1' %(output_directory, ref_path, f.replace('.bam', '_phased.vcf.gz'), f))
                # read haplotags
                haplotags = {}
                inBam = pysam.AlignmentFile('%s/tmp_haplotag.bam' %(output_directory), 'rb', check_sq=False)
                for read in inBam:
                    info = read.tags
                    haplo = 'NA'
                    for x in info:
                        if x[0] == 'HP': haplo = x[1]
                    haplotags[read.query_name] = haplo
                # clean environment
                os.system('rm %s/tmp_*' %(output_directory))
                print('**** %s: phasing and haplotagging done!' %(fname))
                phasing_info[f.split('/')[-1]] = haplotags
            else:
                phasing_info[f.split('/')[-1]] = ['NA']
        fout = open('%s/haplotags_reads.txt' %(output_directory), 'w')
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
    return(phasing_info)

# Function that uses whatshap to assign a haplotype to reads given a vcf file to take snps -- this uses whatshap haplotag
def phase_reads_MP(f, output_directory, SNPs_data_directory, ref_path, bed_file, window):
   # extract info
   tmp_rand, tmp_bam, tmp_sample, sample_name = f[0], f[3], f[4], f[2]
   # write vcf for each sample keeping the snps of interest
   vcf_out = '%s/phasing/sample_snps_%s' %(output_directory, tmp_rand)
   cmd = 'plink2 --pfile %s --extract bed1 %s --bed-border-bp %s --keep %s --recode vcf --out %s >/dev/null 2>&1' %(SNPs_data_directory.replace('.pvar', ''), bed_file, window, tmp_sample, vcf_out)
   os.system(cmd)
   # check if file was created, otherwise skip
   if os.path.isfile('%s.vcf' %(vcf_out)):
      # sort bam file
      os.system('samtools sort %s > %s' %(tmp_bam, tmp_bam + '_sorted'))
      os.system('mv %s %s' %(tmp_bam + '_sorted', tmp_bam))
      os.system('samtools index %s' %(tmp_bam))
      # add chr notation for chromosome to vcf
      os.system('bcftools annotate --rename-chrs %s %s.vcf | bgzip > %s.vcf.gz' %('/'.join(abspath(getsourcefile(lambda:0)).split('/')[:-1]) + '/test_data/chr_name_conv.txt', vcf_out, vcf_out))
      # index vcf
      os.system('tabix %s.vcf.gz' %(vcf_out))
      # create a phased vcf with whatshap
      whathap_out = os.path.basename(tmp_bam).replace('__rawReads.bam', '__phased.vcf.gz')
      haplotag_out = os.path.basename(tmp_bam).replace('__rawReads.bam', '__haplotag.bam')
      os.system('whatshap phase -o %s/phasing/%s --reference=%s %s.vcf.gz %s --ignore-read-groups --internal-downsampling 5 >/dev/null 2>&1' %(output_directory, whathap_out, ref_path, vcf_out, tmp_bam))
      # index the phased vcf
      os.system('tabix %s/phasing/%s' %(output_directory, whathap_out))
      # then tag the haplotypes in the bam file
      os.system('whatshap haplotag -o %s/phasing/%s --reference=%s %s/phasing/%s %s --ignore-read-groups --skip-missing-contigs >/dev/null 2>&1' %(output_directory, haplotag_out, ref_path, output_directory, whathap_out, tmp_bam))
      # also index so that everything is ok
      os.system('samtools index %s/phasing/%s >/dev/null 2>&1' %(output_directory, haplotag_out))
      # read haplotags
      haplotags = []
      inBam = pysam.AlignmentFile('%s/phasing/%s' %(output_directory, haplotag_out), 'rb', check_sq=False)
      for read in inBam:
         info = read.tags
         haplo = 'NA'
         for x in info:
            if x[0] == 'HP':
               haplotags.append([read.query_name, x[1]])
      # clean environment
      os.system('rm %s/phasing/*_%s.*' %(output_directory, tmp_rand))
      print('**** done phasing for %s                                     ' %(sample_name), end = '\r')
   else:
      haplotags = [["NA", "NA"]]
   return haplotags

# Function to find SNPs from GWAS data in the regions of interest given data of interest
def find_SNPs_Samples_plink(SNPs_data_directory, output_directory, snp_data_ids, reads_bam):
   # read GWAS data depending on data_type
   snps_interest = {}
   all_regions = []
   # if a folder was submitted, then look for the chrAll* file
   if os.path.isdir(SNPs_data_directory):
      inpf_path = os.popen('ls ' + SNPs_data_directory + '/chrAll*pvar').readlines()[0].replace('\n', '').replace('.pvar', '')
   else:
      inpf_path = SNPs_data_directory.replace('.pvar', '')
   # define output file
   outf = '%s/phasing/snps_samples_selected' %(output_directory)
   # now check samples
   map_ids = pd.read_csv(snp_data_ids, sep = " |\t", engine = 'python') if snp_data_ids != 'False' else False
   # find gwas id for the corresponding sample
   id_gwas = []
   for x in reads_bam:
      fname = os.path.basename(x).replace('__rawReads.bam', '.bam')
      tmp_match = list(map_ids['ID_GWAS'][map_ids['ID_PACBIO'] == fname]) if isinstance(map_ids, pd.DataFrame) == True else [fname.replace('.bam', '')]
      if tmp_match != []:
          id_gwas.append([tmp_match[0], fname, x])
   # generate as many random numbers as the number of samples to phase
   rand_num = [str(random()).replace('.', '') for x in range(len(id_gwas))]
   # loop to write VCF files for each sample keeping SNPs in the interval of the bed file
   samples_phasing_info = []
   for i in range(len(rand_num)):
      # write sample
      samples_path = '%s/phasing/tmp_sample_%s.txt' %(output_directory, rand_num[i])
      tmp_samples = open(samples_path, 'w')
      header = 'IID\n'; tmp_samples.write(header)
      tmp_samples.write('%s\n' %(id_gwas[i][0]))
      tmp_samples.close()
      samples_phasing_info.append([rand_num[i], id_gwas[i][0], id_gwas[i][1], id_gwas[i][2], samples_path])
   return samples_phasing_info

# Extract reads mapping the the location of interest
def extract_sequences_nowrite_nosort(bam, bed, out_dir, window):
   # dictionary of fasta sequences
   reads_info = []
   reads_ids = {}
   reads_ids_list = []
   # define output names of fasta files with and without padding
   outname_fasta = out_dir + '/' + os.path.basename(bam)[:-4] + '__rawReads.fasta'
   outname_fasta_with_padding = out_dir + '/' + os.path.basename(bam)[:-4] + '__rawReads_withPadding.fasta'
   outname_bam = out_dir + '/' + os.path.basename(bam)[:-4] + '__rawReads.bam'
   # open output fasta with padding
   with open(outname_fasta_with_padding, 'w') as outf_with_padding:
      # open output fasta without padding for TRF
      with open(outname_fasta, 'w') as outf:
         # open bam file
         with pysam.AlignmentFile(bam, 'rb', check_sq=False) as inBam:
            # open output bam file
            with pysam.AlignmentFile(outname_bam, "wb", template=inBam) as outBam:
               # loop over the sorted regions
               for chrom in bed.keys():
                  for region in bed[chrom]:
                     # set object to remember seen reads
                     read_ids_list = set()
                     # reads will be extracted if they encompass the repeat + padding
                     start, end = int(region[0]) - window, int(region[1]) + window
                     # loop over the reads that map there
                     for read in inBam.fetch(chrom, start, end):
                        # extract read name for checking duplicated reads
                        read_name = read.query_name
                        # exclude secondary and supplementary alignments and if read was already processed for this region
                        if (not read.is_secondary) and (not read.is_supplementary) and (read_name not in read_ids_list):
                           # add the read name to the processed ones
                           read_ids_list.add(read_name)
                           # process read
                           res = measureDist(read, chrom, start, end, window)
                           read_id, region_id, np, rq, mc, sequence_interest, sequence_interest_with_padding, sequence_interest_len, sequence_interest_with_padding_len, window = res
                           # write raw fasta read if there were valid results
                           if len(res) >0:
                              # save results
                              reads_ids_list.append('%s;%s' %(read_id, region_id))
                              reads_info.append(res)
                              outf.write('>%s;%s;%s;%s;%s;%s\n%s\n' %(read_id, region_id, np, rq, mc, sequence_interest_len, sequence_interest))
                              outf_with_padding.write('>%s;%s;%s;%s;%s;%s;%s\n%s\n' %(read_id, region_id, np, rq, mc, window, sequence_interest_with_padding_len, sequence_interest_with_padding))
                              outBam.write(read)
   # close open files
   outf_with_padding.close()
   outf.close()
   inBam.close()
   outBam.close()
   reads_ids[os.path.basename(bam)[:-4]] = reads_ids_list
   print("**** done with %s                               " %(os.path.basename(bam)), end = '\r')
   return outname_bam, outname_fasta, outname_fasta_with_padding, reads_info, reads_ids

# Function to generate a coverage profile given a bam file and a bed file -- just counting the number of reads
def generateCoverageProfile(bed, all_bams, window_size, step, output_directory):
    outname = open('%s/coverage_profile.bed' %(output_directory), 'w')
    for bam in all_bams:    # main loop across samples
        inBam = pysam.AlignmentFile(bam, 'rb', check_sq=False)
        bam_name = bam.split('/')[-1]
        for chrom in bed:   # loop on chromosomes
            for region in bed[chrom]:   # loop on single regions
                region_profile = []
                start, end, locus = region
                start, end = int(start) - int(window_size/2), int(end) + int(window_size/2)
                intervals = [x for x in range(start, end, step)]
                for x in intervals:         # loop on intervals
                    counter = 0
                    for read in inBam.fetch(chrom, x, x+step):
                        counter += 1
                    outname.write('%s\t%s\t%s\t%s\t%s\n' %(chrom, x, x+step, bam_name, counter))
        print('**** %s/%s profiles done' %(all_bams.index(bam) + 1, len(all_bams)), end = '\r')
        inBam.close()
    outname.close()
    return('Done')

# Function to generate a coverage profile given a bam file and a bed file -- just counting the number of reads
def generateCoverageProfile_MP(bam, bed, all_bams, window_size, step, output_directory):
    cov_info = {}
    inBam = pysam.AlignmentFile(bam, 'rb', check_sq=False)
    bam_name = bam.split('/')[-1]
    for chrom in bed:   # loop on chromosomes
        for region in bed[chrom]:   # loop on single regions
            start, end, locus = region
            region_id = region[-1]
            start, end = int(start) - int(window_size/2), int(end) + int(window_size/2)
            intervals = [x for x in range(start, end, step)]
            for x in intervals:         # loop on intervals
                counter = 0
                for read in inBam.fetch(chrom, x, x+step):
                    counter += 1
                if bam_name in cov_info.keys():
                    cov_info[bam_name].append([chrom, x, x+step, bam_name, counter, region_id])
                else:
                    cov_info[bam_name] = [[chrom, x, x+step, bam_name, counter, region_id]]
    print('**** profiles done for %s                           ' %(bam_name), end = '\r')
    inBam.close()
    return cov_info

# Function to read target file containing the reads of interest and the target .bam file(s)
def findTargetReads(target, all_bams, output_directory):
    # clean targets
    target = [x.rstrip().split('/ccs')[0] for x in target]
    # list of temporary outputs
    tmp_list = []
    # define output bam
    for bam in all_bams:
        counter = 0
        inBam = pysam.AlignmentFile(bam, 'rb', check_sq=False)
        outname = '%s/tmp_subreads_information_%s.bam' %(output_directory, all_bams.index(bam))
        tmp_list.append(outname)
        outBam = pysam.AlignmentFile(outname, "wb", template=inBam)
        for read in inBam:
            counter += 1
            print('**** %s reads processed' %(counter), end = '\r')
            read_name = '/'.join(read.query_name.split('/')[:2])
            if read_name in target:
                outBam.write(read)
        print('\n')
        inBam.close()
        outBam.close()
    # then merge the bam files in tmp_list and remove the temporary
    cmd = 'samtools merge %s/subreads_information.bam %s' %(output_directory, ' '.join(tmp_list))
    os.system(cmd)
    for bam in tmp_list:
        cmd = 'rm %s' %(bam)
        os.system(cmd)
    return('%s/subreads_information.bam' %(output_directory))

# Function to align raw subreads (non-hifi) reads to the reference genome and chm13
def alignRawReads(target_bam, output_directory):
    # set output names
    outname_hg38 = '%s/sub_information_hg38.bam' %(output_directory)
    outname_chm13 = '%s/sub_information_chm13.bam' %(output_directory)
    # set reference genomes
    hg38 = '/project/holstegelab/Share/pacbio/resources/h38_subread.mmi'
    chm13 = '/project/holstegelab/Share/asalazar/data/chm13/assembly/v2_0/chm13v2.0_subreads.mmi'
    # commands for alignment to hg38
    print('**** aligning to GRCh38')
    cmd = "pbmm2 align --preset SUBREAD %s -j %s --log-level FATAL --sort %s %s" %(hg38, thread, target_bam, outname_hg38)
    os.system(cmd)
    # commands for alignment to chm13
    print('**** aligning to chm13')
    cmd = "pbmm2 align --preset SUBREAD %s -j %s --log-level FATAL --sort %s %s" %(hg38, thread, target_bam, outname_chm13)
    os.system(cmd)
    return('Files are aligned')

# Function to create and save log file with run information
def addLogRun(bed_file, anal_type, var_file, bam_directory, output_directory, store_temporary, window_size, assembly_type, assembly_ploidy, number_threads, polishing, snp_dir, snp_data_ids, step, target_reads, reference, fasta_dir, trf_file, phase_file, asm_file):
    if anal_type == 'realign':
        if os.path.isdir(fasta_dir) == True:
            output_directory = fasta_dir
            if fasta_dir[-1] == '/':
                fasta_dir = fasta_dir[:-1]
        else:
            output_directory = '/'.join(fasta_dir.split('/')[:-1])
    logfile = open('%s/treat_run_info.log' %(output_directory), 'w')
    logfile.write('**  Tandem Repeat Haplotyping Toolkit  **\n\n')
    logfile.write('Please find here belo your run parameters:\n')
    logfile.write('**** GENERAL OPTIONS\n')
    logfile.write('** analysis type --> %s\n' %(anal_type))
    logfile.write('** bed file --> %s\n' %(bed_file))
    logfile.write('** bam file(s) --> %s\n' %(bam_directory))
    logfile.write('** output directory --> %s\n' %(output_directory))
    logfile.write('** store temporary files --> %s\n' %(store_temporary))
    logfile.write('** window size --> %s\n' %(window_size))
    logfile.write('** polishing --> %s\n\n' %(polishing))
    logfile.write('**** ASSEMBLY OPTIONS\n')
    logfile.write('** assembly type --> %s\n' %(assembly_type))
    logfile.write('** ploidy --> %s\n' %(assembly_ploidy))
    logfile.write('** n. cpu --> %s\n' %(number_threads))
    logfile.write('** assembly type --> %s\n' %(assembly_type))
    logfile.write('** reference genome for alignment --> %s\n\n' %(reference))
    logfile.write('**** PHASING OPTIONS\n')
    logfile.write('** snp data directory --> %s\n' %(snp_dir))
    logfile.write('** snp data id mapping --> %s\n\n' %(snp_data_ids))
    logfile.write('**** COVERAGE ANALYSIS OPTIONS\n')
    logfile.write('** step --> %s\n\n' %(step))
    logfile.write('**** SNP ANALYSIS OPTIONS\n')
    logfile.write('** SNP data --> %s\n\n' %(var_file))
    logfile.write('**** SUBREADS ANALYSIS OPTIONS\n')
    logfile.write('** Subreads data --> %s\n\n' %(target_reads))
    logfile.write('**** REALIGNMENT\n')
    logfile.write('** fasta file(s) input --> %s\n' %(fasta_dir))
    logfile.write('** reference genome for alignment --> %s\n\n' %(reference))
    logfile.write('**** HAPLOTYPING\n')
    logfile.write('** TRF file(s) on reads-spanning --> %s\n' %(trf_file))
    logfile.write('** TRF file(s) on assembly --> %s\n' %(asm_file))
    logfile.write('** Phasing file(s) input --> %s\n\n' %(phase_file))
    logfile.close()
    return('** log file created')

# Function to find files to realign. Bam files are the input. Check if there are fasta files with the same name and use those
def getFilesToRealign(fasta_dir):
    if fasta_dir[-1] == '/':
        fasta_dir = fasta_dir[:-1]
    if os.path.isdir(fasta_dir) == True:              # in case a directory was submitted
        all_fasta = [x.rstrip() for x in list(os.popen('ls %s/*.fasta' %(fasta_dir)))]
        print("** found directory with %s fasta" %(len(all_fasta)))
    elif os.path.isfile(fasta_dir) == True:           # in case is a single bam file
        print("** found single fasta")
        all_fasta = [all_fasta]
    elif ',' in fasta_dir:                            # in case there is a comma-separated list of bams
        all_fasta = fasta_dir.split(',')
        print("** found %s fasta" %(len(all_fasta)))
    return all_fasta

# Function to realign fasta files given fasta, threads and reference genome, output directory is the same as the input
def realignFasta(fasta_files, threads, reference):
    # define output files
    output_aligned = [x.replace('.fasta', '_aligned_%s.bam' %(reference)) for x in fasta_files]
    # prepare reference fasta depending on selected reference from user
    reference_mmi = '/project/holstegelab/Share/asalazar/data/chm13/assembly/v2_0/chm13v2.0_hifi.mmi' if reference == 'chm13' else '/project/holstegelab/Share/pacbio/resources/h38_ccs.mmi'
    # loop to align fasta
    for i in range(len(fasta_files)):
        os.system("pbmm2 align --preset CCS %s -j %s --log-level FATAL --sort %s %s" %(reference_mmi, threads, fasta_files[i], output_aligned[i]))       # command for alignment
        print('**** %s/%s alignment done' %(fasta_files.index(fasta_files[i]) + 1, len(fasta_files)), end = '\r')
    return output_aligned

# --------------------------------------------------------------------- #
# functions deleted or modified -- put here for reference
def extractReads_MP(bam, bed, out_dir, window, all_bams):
    fasta_seqs = {}
    read_ids_list = []
    outname = out_dir + '/' + bam.split('/')[-1][:-4] + '__rawReads.bam'
    outname_fasta = out_dir + '/' + bam.split('/')[-1][:-4] + '__rawReads.fasta'
    inBam = pysam.AlignmentFile(bam, 'rb', check_sq=False)
    outBam = pysam.AlignmentFile(outname, "wb", template=inBam)
    for chrom in bed.keys():
        for region in bed[chrom]:
            start, end = int(region[0]) - window, int(region[1]) + window
            for read in inBam.fetch(chrom, start, end):
                read_name = read.query_name
                if read_name not in read_ids_list:
                    read_ids_list.append(read_name)
                    sequence_interest = str(read.query_sequence)    # save fasta sequences
                    fasta_seqs[read_name] = sequence_interest
                    outBam.write(read)                              # write to new bam file
    outBam.close()
    inBam.close()
    outf = open(outname_fasta, 'w')
    for read in fasta_seqs.keys():
        outf.write('>%s\n%s\n' %(read, fasta_seqs[read]))
    outf.close()
    os.system('samtools sort %s > %s' %(outname, outname + '_sorted'))
    os.system('mv %s %s' %(outname + '_sorted', outname))
    os.system('samtools index %s' %(outname))
    print("**** done with %s          " %(bam.split('/')[-1]), end = '\r')
    return outname, outname_fasta
