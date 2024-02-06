# This script manages the assembly-based analysis

# Libraries
print('* Loading libraries')
from functions_assembly_based import *

# Main
# Read arguments and make small changes
inBam_dir, bed_dir, outDir, ref, window, windowAss, cpu, ploidy, software, HaploDev, minimumSupport, minimumCoverage = sys.argv[1::]
window = int(window); cpu = int(cpu); ploidy = int(ploidy); windowAss = int(windowAss); minimumSupport = int(minimumSupport)

# 1. Check arguments: BED, output directory and BAMs
print('** Analysis started')
ts_total = time.time()
# 1.1 Check output directory
print(checkOutDir(outDir))
# 1.2 Create Log file
logfile = createLogAsm(inBam_dir, bed_dir, outDir, ref, window, cpu, windowAss, ploidy, software, HaploDev, minimumSupport, minimumCoverage)
# 1.3 Read bed file
bed, count_reg, bed_dir = readBed(bed_dir, outDir)
# 1.4 Check BAM files
inBam = checkBAM(inBam_dir)

# 2. Check which software was selected and do things accordingly
if software == 'otter':
    # Run local assembly and TRF
    df_trf_phasing_combined = otterPipeline(outDir, cpu, ref, bed_dir, inBam, count_reg, windowAss)
    # Do directly the haplotyping so that we save on IO usage
    ts = time.time()
    print(haplotyping_steps(data = df_trf_phasing_combined, n_cpu = cpu, thr_mad = HaploDev, min_support = minimumSupport, type = 'otter', outDir = outDir))
    te = time.time()
    time_write = te-ts
    print('*** Operation took %s seconds\t\t\t\t\t\t\t\t\t\t\t\t' %(round(time_write, 0)))
    te_total = time.time()
    time_total = te_total - ts_total
    # Remove temporary files
    removeTemp(outDir)
    te_total = time.time()
    time_total = te_total - ts_total
    print('\n** Analysis completed in %s seconds. Ciao!                   ' %(round(time_total, 0)))
elif software == 'hifiasm':
    # Run local assembly with hifiasm
    df_trf_phasing_combined_all = hifiasmPipeline(outDir, cpu, ref, bed_dir, inBam, count_reg, ploidy)
    # Probably good to exclude primary assemblies -- but if not, comment the next line out
    mask = ~df_trf_phasing_combined_all['SAMPLE_NAME'].str.contains('_primary_cleaned')
    df_trf_phasing_combined = df_trf_phasing_combined_all[mask].copy()
    # Save output
    outf = '%s/assembly_trf_phasing.txt.gz' %(outDir)
    df_trf_phasing_combined.to_csv(outf, sep = "\t", index=False, na_rep='NA', compression='gzip')
    print('** Data combined and outputs are ready')
    # Haplotyping
    file_path = os.path.realpath(__file__)
    file_path = '/'.join(file_path.split('/')[:-1])
    os.system("~/.conda/envs/treat/bin/python %s/call_haplotypes.py %s/assembly_trf_phasing.txt.gz %s %s %s %s %s" %(file_path, outDir, outDir, cpu, HaploDev, 'hifiasm', minimumSupport))
    # Remove temporary files
    removeTemp(outDir)
    te_total = time.time()
    time_total = te_total - ts_total
    print('\n** Analysis completed in %s seconds. Ciao!                   ' %(round(time_total, 0)))