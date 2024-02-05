### FUNCTIONS TO CHECK DIRECTORIES AND FILES
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
    print('** BED file: found %s regions in %s chromosomes' %(count_reg, len(bed)))
    return bed, count_reg

# Check directory
def checkOutDir(out_dir):
    if out_dir[-1] == '/':
        out_dir = out_dir[:-1]
    if os.path.isdir(out_dir) == False:
        os.system('mkdir %s' %(out_dir))
        return("** Output directory not found, will create.")
    else:
        return("** Output directory found, will add outputs there.")

# Check bam file(s)
def checkBAM(bam_dir):
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

# Function to create Log file
def createLog(inBam, bed_dir, outDir, ref, window, cpu, windowAss, ploidy, software, HaploDev, minimumSupport, minimumCoverage):
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
# Function to make assembly with otter and produce fasta files suitable for TRF
def assembly_otter(s, output_directory, ref_fasta, bed_file, number_threads, windowAss):
    # define output name with the right directory
    outname = s.split('/')[-1].replace('.bam', '.fa')
    # run otter
    cmd = 'otter assemble -b %s -r %s %s -t %s -o %s > %s/otter_local_asm/%s' %(bed_file, ref_fasta, s, number_threads, windowAss, output_directory, outname)
    os.system(cmd)
    # collect otter sequences
    f_open = [x.rstrip() for x in open('%s/otter_local_asm/%s' %(output_directory, outname), 'r').readlines()]
    # define padding size and output results list
    padd_size = windowAss
    # also define the new fasta output containing only the repetitive sequence for TRF
    trf_input = open('%s/otter_local_asm/%s' %(output_directory, outname.replace('.fa', '_trf.fa')), 'w')
    tmp_res = []  
    prev_region = ''
    for x in f_open:
        if x.startswith('>'):
            region = x.split()[0].replace('>', '')
            read_id = x.replace(' ', '_') + '_1' if region != prev_region else x.replace(' ', '_') + '_2'
            prev_region = region
        else:
            seq = x; seq_len_with_paddings = len(seq); seq_no_paddings = x[(padd_size-1):len(x)-(padd_size+1)]; seq_len_no_paddings = len(seq_no_paddings)
            # save hit at this point
            tmp_res.append([outname.replace('.fa', ''), region, read_id.replace('>', ''), seq, seq_len_with_paddings, seq_no_paddings, seq_len_no_paddings])
            # then write the sequence without paddings
            trf_input.write('>%s;%s;%s\n' %(region, outname.replace('.fa', ''), read_id.replace('>', '')))
            trf_input.write('%s\n' %(seq_no_paddings))
    trf_input.close()
    return tmp_res

# Function to write fasta files for TRF
def writeFastaTRF(all_seqs, fasta_name):
    # define container for fasta outputs
    fasta_outputs = []
    # open file and write things
    with open(fasta_name, 'w') as outFile:
        for region in all_seqs:
            outFile.write('>%s;%s;%s\n%s\n' %(region[1], region[0], region[2], region[-4]))
    outFile.close()

# Run TRF given a sequence
def run_trf_otter(index, all_fasta, distances, type):
    # then run tandem repeat finder
    cmd = 'trf4.10.0-rc.2.linux64.exe %s 2 7 7 80 10 50 200 -ngs -h' %(all_fasta[index])
    trf = [x for x in os.popen(cmd).read().split('\n') if x != '']
    # loop on trf results and save them into a list of lists
    x = 0; trf_matches = []
    sample_name = re.sub(r'^a[a-z]\.tmp_', '', os.path.basename(all_fasta[index])).replace('.fa', '')
    if type == 'otter':
        sample_name = sample_name.replace('_trf', '')
        if sample_name.split('_')[0] == 'reference':
            sample_name = 'reference'
    while x < len(trf):
        # check if the line is the header of an entry
        if trf[x].startswith('@'):
            # if so, save the corresponding information
            region, sample, read_id = trf[x].split(';')
            region = region.replace('@', '')
            x += 1
        while x < len(trf) and not trf[x].startswith('@'):
            tmp_trf_match = [read_id + '_' + region, 'NA'] + trf[x].split()
            trf_matches.append(tmp_trf_match)
            x += 1
    # finally create pandas df and assign column names
    if len(trf_matches) == 0:
        trf_matches = [['NA' for i in range(19)]] 
    df = pd.DataFrame(trf_matches)
    df.columns = ['ID', 'EXPECTED_MOTIF', 'START_TRF', 'END_TRF', 'LENGTH_MOTIF_TRF', 'COPIES_TRF', 'TRF_CONSENSUS_SIZE', 'TRF_PERC_MATCH', 'TRF_PERC_INDEL', 'TRF_SCORE', 'TRF_A_PERC', 'TRF_C_PERC', 'TRF_G_PERC', 'TRF_T_PERC', 'TRF_ENTROPY', 'TRF_MOTIF', 'TRF_REPEAT_SEQUENCE', 'TRF_PADDING_BEFORE', 'TRF_PADDING_AFTER']
    # finally, we need to add the reads where trf didn't find any motif
    if type == 'reads':
        df_seqs = pd.DataFrame(distances[index][0])
        df_seqs.columns = ['SAMPLE_NAME', 'REGION', 'READ_NAME', 'PASSES', 'READ_QUALITY', 'MAPPING_CONSENSUS', 'SEQUENCE_FOR_TRF', 'SEQUENCE_WITH_PADDINGS', 'LEN_SEQUENCE_FOR_TRF', 'LEN_SEQUENCE_WITH_PADDINGS']
    else:
        df_seqs = pd.DataFrame(distances[index])
        if sample_name == 'reference':
            df_seqs.columns = ['SAMPLE_NAME', 'REGION', 'READ_NAME', 'PASSES', 'READ_QUALITY', 'MAPPING_CONSENSUS', 'SEQUENCE_FOR_TRF', 'SEQUENCE_WITH_PADDINGS', 'LEN_SEQUENCE_FOR_TRF', 'LEN_SEQUENCE_WITH_PADDINGS']
            df_seqs['HAPLOTAG'] = 1
        else:
            df_seqs.columns = ['SAMPLE_NAME', 'REGION', 'READ_NAME', 'SEQUENCE_WITH_PADDINGS', 'LEN_SEQUENCE_WITH_PADDINGS', 'SEQUENCE_FOR_TRF', 'LEN_SEQUENCE_FOR_TRF']
            # add haplotag
            df_seqs['HAPLOTAG'] = df_seqs['READ_NAME'].str.split('_').str[-1]
            # add other columns and put NA
            df_seqs['PASSES'] = 'NA'; df_seqs['READ_QUALITY'] = 'NA'; df_seqs['MAPPING_CONSENSUS'] = 'NA'
    # add id
    df_seqs['ID'] = df_seqs['READ_NAME'].str.cat(df_seqs['REGION'], sep='_')
    # merge trf dataframe and reads dataframes
    complete_df = pd.merge(df_seqs, df, left_on = 'ID', right_on = 'ID', how = 'outer')
    return complete_df

# Otter pipeline
def otterPipeline(outDir, cpu, ref, bed_dir, inBam, count_reg):
    print('** Assembler: otter')
    # create directory for outputs
    os.system('mkdir %s/otter_local_asm' %(outDir))
    # run local assembly in multiprocessing
    pool = multiprocessing.Pool(processes=cpu)
    otter_fun = partial(assembly_otter, output_directory = outDir, ref_fasta = ref, bed_file = bed_dir, number_threads = cpu, windowAss = windowAss)
    extract_results = pool.map(otter_fun, inBam)
    pool.close()
    # do the same on the reference genome
    temp_beds = splitBed(bed_dir, cpu, outDir, count_reg)
    pool = multiprocessing.Pool(processes=cpu)
    extract_fun = partial(measureDistance_reference, window = window, ref = ref, output_directory = outDir)
    extract_results_ref = pool.map(extract_fun, temp_beds)
    pool.close()
    all_fasta_ref = [outer_list[1] for outer_list in extract_results_ref]
    extract_results_ref = [outer_list[0] for outer_list in extract_results_ref]
    # combine reference and local assembly
    extract_results.extend(extract_results_ref)
    # TRF on local assemblies and reference
    trf_out_dir = '%s/trf_reads' %(outDir); os.system('mkdir %s' %(trf_out_dir));
    # list fasta files to do analysis on
    reads_fasta = [x.rstrip() for x in os.popen('ls %s/otter_local_asm/*_trf.fa' %(outDir))]
    # add the reference as well
    reads_fasta.extend(all_fasta_ref)
    # Run TRF in multiprocessing for each sample
    pool = multiprocessing.Pool(processes=cpu)
    trf_fun = partial(run_trf_otter, all_fasta = reads_fasta, distances = extract_results, type = 'otter')
    index_fasta = [x for x in range(len(reads_fasta))]
    trf_results = pool.map(trf_fun, index_fasta)
    pool.close()
    #xx = run_trf(index = 0, all_fasta = reads_fasta, distances = extract_results, type = 'otter')
    # Combine df from different samples together
    df_trf_combined = pd.concat(trf_results)
    return df_trf_combined

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

# Function to split bed file in n bed files
def splitBed(bed_dir, n, outDir, count_reg):
    number_lines_per_file = math.ceil(count_reg/n)
    cmd = 'split -l %s %s %s/tmp_bed.' %(number_lines_per_file, bed_dir, outDir)
    os.system(cmd)
    # then read the obtained temporary beds
    cmd = 'ls %s/tmp_bed.*' %(outDir)
    tmp_beds = [x.rstrip() for x in list(os.popen(cmd))]
    return tmp_beds

### FUNCTIONS FOR HIFIASM-BASED ASSEMBLY
def hifiasmPipeline(outDir, cpu, ref, bed_dir, inBam, count_reg, ploidy):
    print('** Assembler: hifiasm')
    # extract reads for assembly
    print("** Extract reads")
    temp_fasta, temp_beds = extractRead(inBam, bed_dir, outDir, cpu, count_reg, 'fasta')
    # do assembly and conversion of gfa to fasta
    print("** Start with assembly: 1 sample at the time")
    print("** Parameters in use for hifiasm: -n 2 -r 3 -N 200")
    assembly_dir = hifiasmAssembly(outDir, cpu, ref, temp_fasta, inBam, ploidy)
    # clean contigs
    assembly_dir_clean = cleanContigs(assembly_dir, cpu)
    # align contigs
    print("** Start with alignment: 1 sample at the time")
    print("** Parameters in use for minimap2: -aYx asm10")
    assembly_aligned = alignContigs(assembly_dir_clean, cpu, ref)
    # then we re-extract sequences of interest like in the read-based method
    print("** Extract sequences of interest")
    temp_bams, temp_beds = extractRead(assembly_aligned, temp_beds, outDir, cpu, count_reg, 'bam')
    # get sequences
    pool = multiprocessing.Pool(processes=cpu)
    extract_fun = partial(distributeExtraction, bed = bed, window = window)
    extract_results = pool.map(extract_fun, temp_bams)
    pool.close()
    print('** Exact SV intervals extracted')
    all_fasta = [outer_list[1] for outer_list in extract_results]
    # then do the same on the reference genome
    pool = multiprocessing.Pool(processes=cpu)
    extract_fun = partial(measureDistance_reference, window = window, ref = ref, output_directory = outDir)
    extract_results_ref = pool.map(extract_fun, temp_beds)
    pool.close()
    all_fasta_ref = [outer_list[1] for outer_list in extract_results_ref]
    print('** Exact SV intervals from reference extracted')
    # combine reference with other samples
    extract_results.extend(extract_results_ref)
    all_fasta.extend(all_fasta_ref)
    # trf
    # run TRF in multiprocessing for each sample
    print("** TRF running")
    pool = multiprocessing.Pool(processes=cpu)
    trf_fun = partial(run_trf, all_fasta = all_fasta, distances = extract_results, type = 'reads')
    index_fasta = [x for x in range(len(all_fasta))]
    trf_results = pool.map(trf_fun, index_fasta)
    pool.close()
    # combine df from different samples together
    df_trf_combined = pd.concat(trf_results)
    print('** TRF done on all reads and samples')
    # assign haplotags
    print("** Assigning haplotag based on assembly")
    df_trf_phasing_combined = assignHaplotags(df_trf_combined)
    return df_trf_phasing_combined

# Function to assign haplotags
def assignHaplotags(df_trf_combined):
    # find samples and regions
    all_samples = list(set(list(df_trf_combined['SAMPLE_NAME'])))
    all_regions = list(set(list(df_trf_combined['REGION'])))
    newdf = pd.DataFrame()
    # iterate over samples and regions
    for s in all_samples:
        for r in all_regions:
            # get data
            tmp = df_trf_combined[(df_trf_combined['SAMPLE_NAME'] == s) & (df_trf_combined['REGION'] == r)].copy()
            # Apply the function to each row to set the HAPLOTAG value
            tmp['HAPLOTAG'] = tmp.apply(set_hap_value, axis=1) if len(tmp) >0 else None
            # then add to new dataframe
            newdf = newdf.append(tmp)
    return newdf

# Set the HAPLOTAG value based on the contig name
def set_hap_value(row):
    if row['SAMPLE_NAME'] == 'reference':
        return 1
    elif 'h1' in row['READ_NAME']:
        return 1
    elif 'h2' in row['READ_NAME']:
        return 2
    else:
        return None

# Function to execute the sequence extraction in multiple processors
def distributeExtraction(x, bed, window):
    # container for results
    tmp_results = []
    # get sample name
    sample_name = re.sub(r'^a[a-z]\.tmp_', '', os.path.basename(x)).replace('.bam', '')
    # loop over reads with pysam
    with pysam.AlignmentFile(x, 'rb', check_sq=False) as bamfile:
        for read in bamfile:
            # extract read information
            ref_chrom, ref_start, ref_end, query_name, query_sequence, cigartuples, tags, is_secondary, is_supplementary, cigarstring = read.reference_name, int(read.reference_start), int(read.reference_end), read.query_name, read.query_sequence, read.cigartuples, read.tags, read.is_supplementary, read.is_secondary, read.cigarstring
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

# Alignment of contigs
def alignContigs(assembly_dir_clean, cpu, ref):
    alignments = []
    # iterate over files
    for asm in assembly_dir_clean:
        print('*** Alignment of %s\t\t\t\t\t\t\t\t\t\t' %(os.path.basename(asm)), end = '\r')
        asm_out = asm.replace('.fasta', '.ref.bam')
        # do alignment
        cmd = "minimap2 -aYx asm10 %s -t %s %s | samtools sort - -@ %s -O bam -o %s" %(ref, cpu, asm, cpu, asm_out); os.system(cmd)
        # index
        cmd = 'samtools index %s' %(asm_out); os.system(cmd)
        alignments.append(asm_out)
    return alignments

# Remove duplicated contigs
def removeDuplicatedContigs(f):
    # read fasta file and create output fasta
    fasta_sequences = SeqIO.parse(open(f), 'fasta')
    # write unique sequences to fasta
    fasta_seqs = []
    with open(f.replace('.fasta', '_cleaned.fasta'), 'w') as fout:
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            if sequence not in fasta_seqs:
                fasta_seqs.append(sequence)
                fout.write('>%s\n%s\n' %(name, sequence))
    fout.close()
    # and remove the original fasta
    os.system('rm %s' %(f))
    return f.replace('.fasta', '_cleaned.fasta')

# Clean assembled contigs
def cleanContigs(assembly_dir, cpu):
    # list the haplotype-aware assemblies in the output folder
    flist = [x.rstrip() for x in os.popen('ls %s*fasta' %(assembly_dir))]
    # run in multiprocessing
    pool = multiprocessing.Pool(processes=cpu)
    extract_results = pool.map(removeDuplicatedContigs, flist)
    pool.close()
    return extract_results

# Run hifiasm assembly
def hifiasmAssembly(outDir, cpu, ref, temp_fasta, inBam, ploidy):
    # create directory for storing assemblies
    os.system('mkdir %s/hifiasm' %(outDir))
    # iterate over the samples
    for s in inBam:
        # get sample name
        s_name = os.path.basename(s).replace('.bam', '')
        print('*** Local Assembly of %s\t\t\t\t\t\t\t\t\t\t' %(s_name), end = '\r')
        # identify the relative fasta files
        fasta_interest = [x for x in temp_fasta if s_name in x]
        # run assembly
        cmd_assembly = 'hifiasm -o %s/hifiasm/%s -t %s --n-hap %s -n 2 -r 3 -N 200 %s >/dev/null 2>&1' %(outDir, s_name, cpu, ploidy, ' '.join(fasta_interest))
        try:
            os.system(cmd_assembly)
            # convert hap1, hap2 and primary contigs to fasta
            subprocess.run("gfatools gfa2fa %s/hifiasm/%s.bp.hap1.p_ctg.gfa > %s/hifiasm/%s_haps.fasta" %(outDir, s_name, outDir, s_name), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            subprocess.run("gfatools gfa2fa %s/hifiasm/%s.bp.hap2.p_ctg.gfa >> %s/hifiasm/%s_haps.fasta" %(outDir, s_name, outDir, s_name), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            subprocess.run("gfatools gfa2fa %s/hifiasm/%s.bp.p_ctg.gfa > %s/hifiasm/%s_primary.fasta" %(outDir, s_name, outDir, s_name), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            # clean temporary files
            subprocess.run("rm %s/hifiasm/%s.*" %(outDir, s_name), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            print('*** Assembly of %s done\t\t\t\t\t\t\t\t\t\t' %(s_name), end = '\r')
        except:
            print('*** Assembly of %s failed. Moving to next sample\t\t\t\t\t\t\t\t\t\t' %(s_name))
    return '%s/hifiasm/' %(outDir)

# Function to extract reads using multiple processors and convert to fasta
def samtoolsExtract(x, bam, out_dir, temp_name, type):
    # combine temporary name with the splitted bed name
    bed_ext = x.split('_bed.')[-1]
    temp_name = '%s/%s.%s' %(out_dir, bed_ext, temp_name)
    if type == 'fasta':
        temp_name = temp_name.replace('.bam', '.fa')
    # define command for the extraction
    cmd = 'samtools view -M -b -L %s %s | samtools fasta > %s' %(x, bam, temp_name) if type == 'fasta' else 'samtools view -M -b -L %s %s > %s' %(x, bam, temp_name)
    os.system(cmd)
    # and index if it's a bam file
    if type == 'bam':
        cmd = 'samtools index %s' %(temp_name); os.system(cmd)
    return temp_name

# Extract reads of interest to temporary FASTA files
def extractRead(bam_dir, bed_dir, out_dir, cpu, count_reg, type):
    # first split the bed files in n smaller bed depending on the cpu number
    split_beds = splitBed(bed_dir, cpu, out_dir, count_reg) if type == 'fasta' else bed_dir
    # make list of temporary outputs
    temp_bams = []
    for bam in bam_dir:
        # define temporary name
        temp_name = 'tmp_' + os.path.basename(bam)
        # define command to extract in MP for each of the splitted bed files
        pool = multiprocessing.Pool(processes=cpu)
        extract_fun = partial(samtoolsExtract, bam = bam, out_dir = out_dir, temp_name = temp_name, type = type)
        extract_results = pool.map(extract_fun, split_beds)
        pool.close()
        # append temporary names
        temp_bams.extend(extract_results)
    return temp_bams, split_beds

# Run TRF given a sequence
def run_trf(index, all_fasta, distances, type):
    # then run tandem repeat finder
    cmd = 'trf4.10.0-rc.2.linux64.exe %s 2 7 7 80 10 50 200 -ngs -h' %(all_fasta[index])
    trf = [x for x in os.popen(cmd).read().split('\n') if x != '']
    # loop on trf results and save them into a list of lists
    x = 0; trf_matches = []
    sample_name = re.sub(r'^a[a-z]\.tmp_', '', os.path.basename(all_fasta[index])).replace('.fa', '')
    while x < len(trf):
        # check if the line is the header of an entry
        if trf[x].startswith('@'):
            # if so, save the corresponding information depending on the type of input
            if type != 'otter' or sample_name == 'reference':
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
