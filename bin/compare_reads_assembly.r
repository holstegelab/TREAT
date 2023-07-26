# Script to compare reads-based and assembly-based approached

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
        # read vcf
        d = fread(rs_vcf, h=T)
        # restrict to region of interest
        if (length(region) == 1 && (region == 'all')){ sub = d } else { sub = d[which(d$ID %in% region),] }
        # check if region existed
        if (nrow(sub) >0){
            return(sub)
        } else {
            stop('The region you provided does not exists! Stopping.')
        }
    }

    # Function to extract sizes
    extractHaploSize = function(vcf){
        # restrict to samples and regions
        tmp = select(vcf, 3, 10:ncol(vcf))
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
    # 1. read read-based and assembly-based VCF files
    read = readVCF('reads/sample.vcf.gz', 'all')
    assembly = readVCF('assembly/sample.vcf.gz', 'all')

    # 2. reformat: extract genotypes
    read_gt = extractHaploSize(read)
    assembly_gt = extractHaploSize(assembly)

    # 3. put data together
    read_gt$Type = 'Read-based'
    assembly_gt$Type = 'Assembly-based'
    # add unique ID with the sample
    read_gt$ID = paste(read_gt$ID, read_gt$Sample, sep='_')
    assembly_gt$ID = paste(assembly_gt$ID, assembly_gt$Sample, sep='_')
    # rename columns
    colnames(read_gt)[2:length(colnames(read_gt))] = paste0(colnames(read_gt)[2:length(colnames(read_gt))], '_Read')
    colnames(assembly_gt)[2:length(colnames(assembly_gt))] = paste0(colnames(assembly_gt)[2:length(colnames(assembly_gt))], '_Assembly')
    combined = merge(read_gt, assembly_gt, by = 'ID')

    # 4. correlation
    cor(combined$Short_Read[which(combined$Sample_Read == 'HG005')], combined$Short_Assembly[which(combined$Sample_Assembly == 'HG005')], use = 'complete', method = 'spearman')
    cor(combined$Long_Read[which(combined$Sample_Read == 'HG005')], combined$Long_Assembly[which(combined$Sample_Assembly == 'HG005')], use = 'complete', method = 'spearman')
    cor(combined$Sum_Read, combined$Sum_Assembly, use = 'complete', method = 'spearman')

    # 5. plot
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


while (TRUE) {
  print("ok")
  Sys.sleep(20)
}