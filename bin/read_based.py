# This script manages the read-based analysis

# Libraries
print('* Loading libraries')
from functions_read_based import *

# Main
# Read arguments and make small changes
inBam_dir, bed_dir, outDir, ref, window, cpu, phasingData, mappingSNP, HaploDev, minimumSupport, minimumCoverage, rawSequences = sys.argv[1::]
window = int(window); cpu = int(cpu); minimumSupport = int(minimumSupport)
if HaploDev == 'None':
    HaploDev = 0.10
else:
    HaploDev = float(HaploDev)

# 1. Check arguments: BED, output directory and BAMs
print('* Analysis started')
ts_total = time.time()
# 1.1 Check output directory
print(checkOutDir(outDir))
# 1.2 Create Log file
logfile = createLogReads(inBam_dir, bed_dir, outDir, ref, window, cpu, phasingData, mappingSNP, HaploDev, minimumSupport, minimumCoverage)
# 1.3 Read bed file
bed, count_reg, bed_dir = readBed(bed_dir, outDir)
# 1.3 Check BAM files
inBam = checkBAM(inBam_dir)

# 2. Extract sequence of interest
ts = time.time()
# 2.1 Extract reads using samtools
temp_bams, temp_beds = extractRead(inBam, bed_dir, outDir, cpu, count_reg)
# 2.2 Parse output and get sequences
pool = multiprocessing.Pool(processes=cpu)
extract_fun = partial(distributeExtraction, bed = bed, window = window)
extract_results = pool.map(extract_fun, temp_bams)
pool.close()
print('** Exact SV intervals extracted')
all_fasta = [outer_list[1] for outer_list in extract_results]
all_clipping = [outer_list[2] for outer_list in extract_results]
all_clipping_flatten = [item for sublist in all_clipping for item in sublist]
all_clipping_df = pd.DataFrame(all_clipping_flatten, columns=['REGION', 'SAMPLE', 'READ_NAME'])
# 2.3 Then do the same on the reference genome
pool = multiprocessing.Pool(processes=cpu)
extract_fun = partial(measureDistance_reference, window = window, ref = ref, output_directory = outDir)
extract_results_ref = pool.map(extract_fun, temp_beds)
pool.close()
all_fasta_ref = [outer_list[1] for outer_list in extract_results_ref]
print('** Exact SV intervals from reference extracted')
# 2.5 combine reference with other samples
extract_results.extend(extract_results_ref)
all_fasta.extend(all_fasta_ref)
te = time.time()
time_extraction = te-ts
print('** Read extraction took %s seconds\t\t\t\t\t\t\t\t\t\t' %(round(time_extraction, 0)))

# 3. TRF
ts = time.time()
# 3.1 Run TRF in multiprocessing for each sample
pool = multiprocessing.Pool(processes=cpu)
trf_fun = partial(run_trf, all_fasta = all_fasta, type = 'reads')
index_fasta = [x for x in range(len(all_fasta))]
trf_results = pool.map(trf_fun, index_fasta)
pool.close()
# 3.2 combine df from different chunks together
df_trf_combined = combineTRF_res(trf_results, extract_results, all_fasta)
print('** TRF done on all reads and samples')
te = time.time()
time_trf = te-ts
print('*** TRF took %s seconds\t\t\t\t\t\t\t\t\t\t' %(round(time_trf, 0)))

# 4. Phasing and haplotagging and combine with sequences
ts = time.time()
# 4.1 Check whether we need to do this
if phasingData == 'None':
    print('** Phasing NOT selected (not specified any SNP data)')
    combined_haplotags_df = pd.DataFrame(columns=['READ_NAME', 'HAPLOTAG'])
else:
    print('** Phasing and haplotagging with whatshap')
    # create directory for phasing
    ts = time.time()
    os.system('mkdir %s/phasing' %(outDir))
    print('** Phasing started\t\t\t\t\t\t\t\t\t\t\t')
    pool = multiprocessing.Pool(processes=cpu)
    phasing_fun = partial(phase_reads, temp_bams = temp_bams, temp_beds = temp_beds, phasingData = phasingData, mappingSNP = mappingSNP, outDir = outDir, snpWindow = 10000)
    #tmp = phase_reads(5, temp_bams = temp_bams, temp_beds = temp_beds, phasingData = phasingData, mappingSNP = mappingSNP, outDir = outDir, snpWindow = 10000)
    phasing_res = pool.map(phasing_fun, [i for i in range(len(temp_bams))])
    pool.close()
    te = time.time()
    time_phasing = te-ts
    print('** Phasing done in %s seconds\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t' %(round(time_phasing, 0)))
    combined_haplotags = sum(phasing_res, [])
    if combined_haplotags == []:
        combined_haplotags_df = pd.DataFrame(columns=['READ_NAME', 'HAPLOTAG'])
    else:
        combined_haplotags_df = pd.DataFrame(combined_haplotags, columns = ['READ_NAME', 'HAPLOTAG'])
    # combine phased VCF and haplotagged bam files
    combined_data = combine_data_afterPhasing(outDir)
    print('*** Phasing took %s seconds\t\t\t\t\t\t\t\t\t\t' %(round(time_phasing, 0)))
# 4.2 Combine with sequences
df_trf_phasing_combined = pd.merge(df_trf_combined, combined_haplotags_df, left_on = 'READ_NAME', right_on = 'READ_NAME', how = 'outer')
# check with some actual data where phasing is expected

# 5. Do directly the haplotyping so that we save on IO usage
ts = time.time()
df_seq = haplotyping_steps(data = df_trf_phasing_combined, n_cpu = cpu, thr_mad = HaploDev, min_support = minimumSupport, type = 'reads', outDir = outDir, all_clipping_df = all_clipping_df)
te = time.time()
time_write = te-ts
print('*** Operation took %s seconds\t\t\t\t\t\t\t\t\t\t\t\t' %(round(time_write, 0)))
te_total = time.time()
time_total = te_total - ts_total

# 6. Output also the raw data
ts = time.time()
# 6.1 Output file for haplotyping if requested
if rawSequences == 'True':
    outf = '%s/spanning_reads_trf_phasing.txt.gz' %(outDir)
    print('** Writing raw data sequences')
    df_seq.to_csv(outf, sep = " ", index=False, na_rep='NA', compression='gzip')
# 6.2 Removing temporary files
print('** Cleaning')
tmp = removeTemp(outDir)
print('\n* Analysis completed in %s seconds. Ciao!\t\t\t\t\t\t\t\t' %(round(time_total, 0)))
