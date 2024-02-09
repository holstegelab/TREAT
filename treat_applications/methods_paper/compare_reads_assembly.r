# Script to compare reads-based and assembly-based approached

# Script to compare reads-based and assembly-based methods
# This script accompanies Treat/Otter manuscript

# Libraries
    library(data.table)
    library(stringr)
    library(corrplot)
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

    # Function to do QC on data
    QC_beforePCA = function(combined, type){
        # extract features of interest
        if (type == 'asm'){
            combined_sb = combined[, c('ID', 'Sum_Assembly')]
        } else {
            combined_sb = combined[, c('ID', 'Sum_Read')]
        }
        combined_sb$Sample = stringr::str_split_fixed(combined_sb$ID, '_', 2)[, 2]
        combined_sb$Region = stringr::str_split_fixed(combined_sb$ID, '_', 2)[, 1]
        combined_sb$ID = NULL
        # reshape
        if (type == 'asm'){
            wide_df <- data.frame(pivot_wider(data = combined_sb, names_from = Region, values_from = Sum_Assembly))
        } else {
            wide_df <- data.frame(pivot_wider(data = combined_sb, names_from = Region, values_from = Sum_Read))
        }
        # calculate mean, median and standard deviation
        column_medians <- apply(wide_df, 2, median, na.rm = TRUE)
        column_sd <- apply(wide_df, 2, sd, na.rm = TRUE)
        column_sd_vector = as.vector(column_sd)
        # select regions variables - standard deviation >50%
        threshold_quantile = as.vector(quantile(column_sd_vector, probs=0.50, na.rm=T))
        # select these columns
        combined_sb_variable = wide_df[, c(1, which(column_sd_vector > threshold_quantile))]
        rownames(combined_sb_variable) = combined_sb_variable$Sample
        combined_sb_variable$Sample = NULL
        rownames(wide_df) = wide_df$Sample
        wide_df$Sample = NULL
        # combine full results and the subset of variables
        res = list(combined_sb_variable, wide_df)
        return(res)       
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
        cor_all_samples_short_pea = cor.test(combined$Short_Read, combined$Short_Assembly, use = 'complete', method = 'pearson')
        cor_all_samples_long_pea = cor.test(combined$Long_Read, combined$Long_Assembly, use = 'complete', method = 'pearson')
        cor_all_samples_sum_pea = cor.test(combined$Sum_Read, combined$Sum_Assembly, use = 'complete', method = 'pearson')

    # 10. save workspace
        save.image('20240207_workspace.RData')

    # 11. plots
        # scatterplots
            ggplot(data = combined, aes(x = Short_Read, y = Short_Assembly)) + geom_point(stat = 'identity', alpha = 0.5) + geom_smooth(method = 'lm') + geom_abline(slope = 1, intercept = 0, linetype = 'dashed') + xlab('Read-based') + ylab('Assembly-based') + ggtitle('Short allele')
            ggplot(data = combined, aes(x = Long_Read, y = Long_Assembly)) + geom_point(stat = 'identity', alpha = 0.5) + geom_smooth(method = 'lm') + geom_abline(slope = 1, intercept = 0, linetype = 'dashed') + xlab('Read-based') + ylab('Assembly-based') + ggtitle('Long allele')
            ggplot(data = combined, aes(x = Sum_Read, y = Sum_Assembly)) + geom_point(stat = 'identity', alpha = 0.5) + geom_smooth(method = 'lm') + geom_abline(slope = 1, intercept = 0, linetype = 'dashed') + xlab('Read-based') + ylab('Assembly-based') + ggtitle('Sum of alleles')
            ggarrange(short_plot, long_plot, sum_plot, nrow = 1, ncol = 3)
        # alternatively, a correlation plot could also be ok
            corr_mt_spe = cor(combined[, c('Short_Read', 'Long_Read', 'Sum_Read', 'Short_Assembly', 'Long_Assembly', 'Sum_Assembly')], method = 'spearman', use = 'complete.obs')
            corr_mt_pea = cor(combined[, c('Short_Read', 'Long_Read', 'Sum_Read', 'Short_Assembly', 'Long_Assembly', 'Sum_Assembly')], method = 'pearson', use = 'complete.obs')
            corrplot(corr = corr_mt_pea)
        # maybe, just not show a plot

    # 12. in a subset of cases there is deviation between reads and assembly - check a few examples
        # some deviation in the long ones
            combined$Difference_Reads_Asm = abs(combined$Sum_Assembly - combined$Sum_Read)
            combined = combined[order(-combined$Difference_Reads_Asm),]
            head(combined)
            summary(combined$Difference_Reads_Asm)
            # multiple divergence at chr2:41686495-41686521: the reads analysis seems to be missing an allele -- extract the bam files and check
            system('mkdir /project/holstegelab/Share/nicco/paper_treat/20240206_final/comparison_reads_asm/examples_divergence')
            system('samtools view -b /project/holstegelab/Share/pacbio/public_data/hpc/hifi/alignments/mm2/grch38/HG03098.bam chr2:41686495-41686521 > /project/holstegelab/Share/nicco/paper_treat/20240206_final/comparison_reads_asm/examples_divergence/HG03098_divergence_example.bam')
            system('samtools index /project/holstegelab/Share/nicco/paper_treat/20240206_final/comparison_reads_asm/examples_divergence/HG03098_divergence_example.bam')
            # in conclusion, it looks heterozygous with an insertion, but there's only 1 read spanning the large insertion, the others are left clipped. Maybe in similar cases, detect that a problem may arise and report this in the VCF file.
    # 13. pca with population information
        # extract samples information
            samples_info = fread('/project/holstegelab/Share/nicco/paper_treat/20240206_final/population_analysis/samples_info_hprc.txt', h=T, stringsAsFactors=F)
        # we should do some qc on the data - use function for it
            combined_data_forPCA_res = QC_beforePCA(combined, 'reads')
            subset_variable = combined_data_forPCA_res[[1]]
            all_variables = combined_data_forPCA_res[[2]]
        # do pca
            pca = prcomp(subset_variable))
            pca_all = prcomp(all_variables)
        # add samples information
            pca_components = pca$x
            pca_components = merge(pca_components, samples_info, by.x = 'row.names', by.y = 'Sample')
            pca_components_all = pca_all$x
            pca_components_all = merge(pca_components_all, samples_info, by.x = 'row.names', by.y = 'Sample')
        # plot
            ggplot(data = pca_components, aes(x = PC1, y = PC2, color = Superpopulation)) + geom_point(stat = 'identity', size = 5, alpha = 0.7) + scale_color_manual(values = RColorBrewer::brewer.pal(n = 10, name = "Dark"))
            ggplot(data = pca_components_all, aes(x = PC1, y = PC2, color = Superpopulation)) + geom_point(stat = 'identity', size = 5, alpha = 0.7) + scale_color_manual(values = RColorBrewer::brewer.pal(n = 10, name = "Paired"))
