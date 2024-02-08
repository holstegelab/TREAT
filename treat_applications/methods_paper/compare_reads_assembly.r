# Script to compare reads-based and assembly-based approached

# Script to compare reads-based and assembly-based methods
# This script accompanies Treat/Otter manuscript

# Libraries
    library(data.table)
    library(stringr)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(ggpubr)

# Functions
    # Function to read VCF, and restrict to region of interest
    readVCF <- function(rs_vcf, region){
        combined_data = data.frame()
        for (vcf in rs_vcf){
            print(paste0('** reading ', vcf))
            # read vcf
            d = fread(vcf, h=T)
            # restrict to region of interest
            if (length(region) == 1 && (region == 'all')){ sub = data.frame(d) } else { sub = data.frame(d[which(d$ID %in% region),]) }
            # check if region existed
            if (nrow(combined_data) == 0){
                combined_data = sub
            } else {
                sub = sub[, c(3, ncol(sub))]
                combined_data = merge(combined_data, sub, by = 'ID')
            }
        }
        return(combined_data)
    }

    # Function to extract sizes
    extractHaploSize = function(vcf){
        # restrict to samples and regions
        tmp = select(vcf, 1, 10:ncol(vcf))
        # split and keep genotypes
        tmp_split <- tmp %>% mutate(across(2:last_col(), ~str_split(., ";"))) %>% mutate(across(2:last_col(), ~sapply(., function(x) x[2])))
        # Reshape the dataframe from wide to long format, and add a column with the sample name
        tmp_split_reshape <- tmp_split %>% gather(key = "Sample", value = "Value", -ID)
        # add short, long and sum of alleles
        tmp_split_reshape <- tmp_split_reshape %>% mutate(Short = sapply(strsplit(Value, "\\|"), function(x) min(as.numeric(x))))
        tmp_split_reshape <- tmp_split_reshape %>% mutate(Long = sapply(strsplit(Value, "\\|"), function(x) max(as.numeric(x))))
        tmp_split_reshape <- tmp_split_reshape %>% mutate(Sum = sapply(strsplit(Value, "\\|"), function(x) sum(as.numeric(x))))
        return(tmp_split_reshape)
    }

# Main
    # 1. define paths
        data_path = '/project/holstegelab/Share/nicco/paper_treat/20240206_final'

    # 2. create directory for comparative analysis
        system(paste0('mkdir ', data_path, '/comparison_reads_asm'))

    # 3. identify the samples's data (individual data)
        reads_vcf = system(paste0("find ", data_path, "/ -name '*vcf*' | grep reads"), intern=T)
        asm_vcf = system(paste0("find ", data_path, "/ -name '*vcf*' | grep asm"), intern=T)

    # 4. iterate through files, read them and combine
        reads_data = readVCF(reads_vcf, 'all')
        asm_data = readVCF(asm_vcf, 'all')

    # 5. reformat: extract genotypes
        reads_gt = extractHaploSize(reads_data)
        asm_gt = extractHaploSize(asm_data)

    # 6. put data together
        reads_gt$Type = 'Read-based'
        asm_gt$Type = 'Assembly-based'

    # 7. add unique ID with the sample and the region
        reads_gt$ID = paste(reads_gt$ID, reads_gt$Sample, sep='_')
        asm_gt$ID = paste(asm_gt$ID, asm_gt$Sample, sep='_')

    # 8. rename columns and combine reads and assembly
        colnames(reads_gt)[2:length(colnames(reads_gt))] = paste0(colnames(reads_gt)[2:length(colnames(reads_gt))], '_Read')
        colnames(asm_gt)[2:length(colnames(asm_gt))] = paste0(colnames(asm_gt)[2:length(colnames(asm_gt))], '_Assembly')
        combined = merge(reads_gt, asm_gt, by = 'ID')

    # 9. correlation
        cor_all_samples_short = cor.test(combined$Short_Read, combined$Short_Assembly, use = 'complete', method = 'spearman')
        cor_all_samples_long = cor(combined$Long_Read, combined$Long_Assembly, use = 'complete', method = 'spearman')
        cor_all_samples_sum = cor(combined$Sum_Read, combined$Sum_Assembly, use = 'complete', method = 'spearman')

    # 10. save workspace
        save.image('20240207_workspace.RData')

    # 11. plot
        # randomly subsample for the plot
        sampled_data <- combined %>% sample_n(100000)
        # then plot
        ggplot(data = sampled_data, aes(x = Short_Read, y = Short_Assembly)) + stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + theme(legend.position='none')


    ggplot(data = combined, aes(x = Short_Read, y = Short_Assembly)) + geom_point(stat = 'identity', alpha = 0.5) + geom_smooth(method = 'lm') + geom_abline(slope = 1, intercept = 0, linetype = 'dashed') + xlab('Read-based') + ylab('Assembly-based') + ggtitle('Short allele')
    ggplot(data = combined, aes(x = Long_Read, y = Long_Assembly)) + geom_point(stat = 'identity', alpha = 0.5) + geom_smooth(method = 'lm') + geom_abline(slope = 1, intercept = 0, linetype = 'dashed') + xlab('Read-based') + ylab('Assembly-based') + ggtitle('Long allele')
    ggplot(data = combined, aes(x = Sum_Read, y = Sum_Assembly)) + geom_point(stat = 'identity', alpha = 0.5) + geom_smooth(method = 'lm') + geom_abline(slope = 1, intercept = 0, linetype = 'dashed') + xlab('Read-based') + ylab('Assembly-based') + ggtitle('Sum of alleles')
    ggarrange(short_plot, long_plot, sum_plot, nrow = 1, ncol = 3)

    # some deviation in the long ones
    combined[which(combined$Long_Assembly >300 & combined$Sample_Assembly == 'HG005'),]
    # extract reads of these and manually look at them
    regions_divergent = str_split_fixed(combined$ID[which(combined$Long_Assembly >300 & combined$Sample_Assembly == 'HG005')], '_', 2)[, 1]
    regions_divergent_df = data.frame(str_split_fixed(regions_divergent, ':', 2)[, 1], str_split_fixed(str_split_fixed(regions_divergent, ':', 2)[, 2], '-', 2)[, 1], str_split_fixed(str_split_fixed(regions_divergent, ':', 2)[, 2], '-', 2)[, 2])
    # extract reads from bam file
    bam = '/project/holstegelab/Share/pacbio/public_data/hpc/hifi/alignments/mm2/grch38/HG005.bam'
    write.table(regions_divergent_df, 'divergent_regions.bed', quote=F, row.names=F, col.names=F)
    cmd = paste0('samtools view -M -b -L divergent_regions.bed ', bam, ' > divergent_regions.bam')
    system(cmd)
    system('samtools index divergent_regions.bam')
    # check assembly and reads
    read[which(read$ID == 'chr1:54510616-54510666'),]
    assembly[which(assembly$ID == 'chr1:54510616-54510666'),]
