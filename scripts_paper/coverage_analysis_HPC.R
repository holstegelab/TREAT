library(stringr)
setwd('/project/holstegelab/Share/nicco/paper_treat/20230802_final_runs/coverages')
bam_path = '/project/holstegelab/Share/pacbio/public_data/hpc/hifi/alignments/mm2/grch38/'
bam_files = system(paste0('ls ', bam_path, '*bam'), intern = T)
bam_files
all_bam_info = data.frame()
    for (bam in bam_files){
        print(bam)
        tmp_res = data.frame(str_split_fixed(system(paste0('asbt cov -g 3088000000 ', bam), intern=T), '\t', 8))
        colnames(tmp_res) = c('Global_median_read_length', 'Global_coverage', 'Mapped_coverage', 'Mapped_reads', 'Alt_coverage', 'Alt_reads', 'Unmapped_coverage', 'Unmapped_reads')
        all_bam_info = rbind(all_bam_info, tmp_res)
    }
ls()
head(all_bam_info)
dim(all_bam_info)
all_bam_info
library(stringr)
library(stringr)
ls()
save.image('coverage_HPC_hifi_hg38.RData')

# 