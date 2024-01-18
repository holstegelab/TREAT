# Libraries and working directories
    library(ggpubr)
    library(ggplot2)
    library(stringr)
    library(data.table)
    setwd('/project/holstegelab/Share/nicco/paper_treat/20230802_final_runs/coverages')

# File paths
    bam_path = '/project/holstegelab/Share/pacbio/public_data/hpc/hifi/alignments/mm2/grch38/'
    bam_files = system(paste0('ls ', bam_path, '*bam'), intern = T)
    
# Iterate over files, calculate coverage and save information
    all_bam_info = data.frame()
        for (bam in bam_files){
            print(bam)
            tmp_res = data.frame(str_split_fixed(system(paste0('asbt cov -g 3088000000 ', bam), intern=T), '\t', 8))
            colnames(tmp_res) = c('Global_median_read_length', 'Global_coverage', 'Mapped_coverage', 'Mapped_reads', 'Alt_coverage', 'Alt_reads', 'Unmapped_coverage', 'Unmapped_reads')
            all_bam_info = rbind(all_bam_info, tmp_res)
        }

# Save environment
    save.image('coverage_HPC_hifi_hg38.RData')

# Re-load environment
    load('coverage_HPC_hifi_hg38.RData')

# Extract sample names
    sample_names = c()
    for (x in bam_files){ tmp = unlist(strsplit(x, '/'))[[length(unlist(strsplit(x, '/')))]]; tmp = str_replace_all(tmp, '.bam', ''); sample_names = c(sample_names, tmp) }

# Get sample information
    sample_info = fread('/project/holstegelab/Share/gwas_array/1000Genome_common_snps/samples_description.txt', h=T, sep=" ")
    table(sample_names %in% sample_info$Sample)
    # 4 samples are not in the list
    sample_names[which(!(sample_names %in% sample_info$Sample))]
    # missing HG002 - HG005 - HG02109 - NA21309
    missing_samples = data.frame(Sample = c('HG02109', 'HG002', 'HG005', 'NA21309'), 'Biosample' = c('ACB', NA, NA, NA), 'Population' = c('AFR', 'ASH', 'Chinese', 'AFR'))
    sample_info = rbind(sample_info, missing_samples)
    table(sample_names %in% sample_info$Sample)
    # find subset of samples of interest
    sample_hpc = sample_info[which(sample_info$Sample %in% sample_names),]
    table(sample_hpc$Population)
    # make main populations
    sample_hpc$SuperPop = NA
    sample_hpc$SuperPop[which(sample_hpc$Population %in% c('AFR', 'Mandinka'))] = 'African'
    sample_hpc$SuperPop[which(sample_hpc$Population %in% c('ASH', 'AMR', 'Rican'))] = 'American'
    sample_hpc$SuperPop[which(sample_hpc$Population %in% c('Chinese', 'Han', 'SAS', 'SW', 'Vietnamese'))] = 'Asian'
    table(sample_hpc$SuperPop, exclude=F)

# Add sample name to bam stats -- assuming the same order which should be the case
    all_bam_info$Sample = sample_names
    all_bam_info_pheno = merge(all_bam_info, sample_hpc, by = 'Sample')

# Adjust variables for the plot
    all_bam_info_pheno$Global_coverage = as.numeric(all_bam_info_pheno$Global_coverage)
    all_bam_info_pheno$Global_median_read_length = as.numeric(all_bam_info_pheno$Global_median_read_length)
    tb = data.frame(table(all_bam_info_pheno$SuperPop)/nrow(all_bam_info))
    tb = tb[order(tb$Freq),]
    tb$cum = cumsum(tb$Freq)
    tb$pos = NA
    for (i in 1:nrow(tb)){ if (i == 1){ tb$pos[i] = tb$Freq[i]/2 } else { tb$pos[i] = tb$cum[i-1] + tb$Freq[i]/2 } }
    tb$Var1 = as.character(tb$Var1)
    tb$Var1[which(tb$Var1 == 'African')] = 'African\n(n=23)'
    tb$Var1[which(tb$Var1 == 'American')] = 'American\n(n=17)'
    tb$Var1[which(tb$Var1 == 'Asian')] = 'Asian\n(n=7)'
    #tb$Var1 = factor(tb$Var1, levels = c('African', 'American', 'Asian'))
    
# Plots -- 3 panels: mean coverage of all samples, mean read-length of all samples, ethnicity
    # coverage
    plot1 = ggplot(all_bam_info_pheno, aes(x = SuperPop, y = Global_coverage, fill = SuperPop)) + geom_violin() + geom_boxplot(width = 0.10) + geom_jitter(stat = 'identity', width = 0.05, alpha = 0.7, size = 2) + xlab('Super Populations') + ylab('Global Coverage') + theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18), legend.position = 'None')
    # read-length
    plot2 = ggplot(all_bam_info_pheno, aes(x = SuperPop, y = Global_median_read_length, fill = SuperPop)) + geom_violin() + geom_boxplot(width = 0.10) + geom_jitter(stat = 'identity', width = 0.05, alpha = 0.7, size = 2) + xlab('Super Populations') + ylab('Median Read Length') + theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18), legend.position = 'None')
    # sample distribution
    plot3 = ggplot() + geom_col(aes(x = 1, y = Freq, fill = Var1), data = tb, color = "black") + geom_text(aes(label = Var1, x = 1, y = pos), data = tb, size = 5) + xlim(0, 2) + labs(x = NULL, y = NULL) + theme(axis.ticks=element_blank(), axis.text=element_blank(), axis.title=element_blank())+ coord_polar(theta = "y") + theme_classic() + theme(axis.text = element_blank(), axis.title = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), legend.position = 'None')
    # combine plots
    plt_combined = ggarrange(plot3, plot1, plot2, nrow=1, ncol=3, labels=c('A', 'B', 'C'), font.label = list(size = 19))
    pdf('figure1_summary_data.pdf', height = 7, width = 18)
    plt_combined
    dev.off()



