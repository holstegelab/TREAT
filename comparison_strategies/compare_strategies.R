# Script to compare the phasing and non-phasing and the assembly based method

# Libraries
    library(data.table)
    library(stringr)
    library(ggplot2)
    library(ggpubr)

# Functions
    # Function to extract sizes
    extractHaploSize = function(vcf){
        # restrict to sample's columns
        sub = vcf[, 10:ncol(vcf)]
        # take the transpose of this
        sub_trans = t(sub)
        # split columns
        sub_trans_df = data.frame(str_split_fixed(sub_trans[, 1], ';', 5))
        sub_trans_df$sample = rownames(sub_trans)
        # rename columns
        colnames(sub_trans_df) = c('qc', 'haplo', 'motifs', 'copies', 'copies_reference', 'sample')
        # add phased values for short and long haplo
        sub_trans_df$haplo_h1 = as.numeric(str_split_fixed(sub_trans_df$haplo, '\\|', 2)[, 1])
        sub_trans_df$haplo_h2 = as.numeric(str_split_fixed(sub_trans_df$haplo, '\\|', 2)[, 2])
        # define short and long alleles and link to the motifs and number of copies
        sub_trans_df$short_allele = NA; sub_trans_df$long_allele = NA; sub_trans_df$short_allele_motif = NA; sub_trans_df$long_allele_motif = NA; sub_trans_df$short_allele_cn = NA; sub_trans_df$long_allele_cn = NA; sub_trans_df$short_allele_cn_ref = NA; sub_trans_df$long_allele_cn_ref = NA
        for (x in 1:nrow(sub_trans_df)){
            if (!is.na(sub_trans_df$haplo_h1[x])){
                if (sub_trans_df$haplo_h1[x] > sub_trans_df$haplo_h2[x]){
                    # haplotype size
                    sub_trans_df$short_allele[x] = sub_trans_df$haplo_h2[x]; sub_trans_df$long_allele[x] = sub_trans_df$haplo_h1[x]
                    # motifs
                    sub_trans_df$short_allele_motif[x] = str_split_fixed(sub_trans_df$motifs[x], '\\|', 2)[, 2]; sub_trans_df$long_allele_motif[x] = str_split_fixed(sub_trans_df$motifs[x], '\\|', 2)[, 1]
                    # copy number
                    sub_trans_df$short_allele_cn[x] = as.numeric(str_split_fixed(sub_trans_df$copies[x], '\\|', 2)[, 2]); sub_trans_df$long_allele_cn[x] = str_split_fixed(sub_trans_df$copies[x], '\\|', 2)[, 1]
                    # copy number of reference
                    sub_trans_df$short_allele_cn_ref[x] = as.numeric(str_split_fixed(sub_trans_df$copies_reference[x], '\\|', 2)[, 2]); sub_trans_df$long_allele_cn_ref[x] = str_split_fixed(sub_trans_df$copies_reference[x], '\\|', 2)[, 1]
                } else {
                    # haplotype size
                    sub_trans_df$short_allele[x] = sub_trans_df$haplo_h1[x]; sub_trans_df$long_allele[x] = sub_trans_df$haplo_h2[x]
                    # motifs
                    sub_trans_df$short_allele_motif[x] = str_split_fixed(sub_trans_df$motifs[x], '\\|', 2)[, 1]; sub_trans_df$long_allele_motif[x] = str_split_fixed(sub_trans_df$motifs[x], '\\|', 2)[, 2]
                    # copy number
                    sub_trans_df$short_allele_cn[x] = as.numeric(str_split_fixed(sub_trans_df$copies[x], '\\|', 2)[, 1]); sub_trans_df$long_allele_cn[x] = str_split_fixed(sub_trans_df$copies[x], '\\|', 2)[, 2]
                    # copy number of reference
                    sub_trans_df$short_allele_cn_ref[x] = as.numeric(str_split_fixed(sub_trans_df$copies_reference[x], '\\|', 2)[, 1]); sub_trans_df$long_allele_cn_ref[x] = str_split_fixed(sub_trans_df$copies_reference[x], '\\|', 2)[, 2]
                }
            } else {
                # haplotype size
                sub_trans_df$short_allele[x] = NA; sub_trans_df$long_allele[x] = NA
                # motifs
                sub_trans_df$short_allele_motif[x] = NA; sub_trans_df$long_allele_motif[x] = NA
                # copy number
                sub_trans_df$short_allele_cn[x] = NA; sub_trans_df$long_allele_cn[x] = NA
                # copy number of reference
                sub_trans_df$short_allele_cn_ref[x] = NA; sub_trans_df$long_allele_cn_ref[x] = NA
            }
        }
        # subset of the fields of interest
        sub_trans_df = sub_trans_df[, c('sample', 'short_allele', 'long_allele', 'short_allele_motif', 'long_allele_motif', 'short_allele_cn', 'long_allele_cn', 'short_allele_cn_ref', 'long_allele_cn_ref')]
        return(sub_trans_df)
        }
# Main
    # Define input directory
    inp_dir = '/project/holstegelab/Share/nicco/paper_treat/20230404_hpc_treat'

    # Read treat output with phasing
    treat_phasing = fread(file.path(inp_dir, 'treat_phasing/samples_genotypes.vcf'), h=T, stringsAsFactors=F)

    # Read treat output without phasing
    treat = fread(file.path(inp_dir, 'treat_no_phasing/samples_genotypes.vcf'), h=T, stringsAsFactors=F)

    # Read treat output with asm


    # Extract short allele and long alleles
    all_regions = unique(c(treat_phasing$ID, treat$ID)) 
    # Loop across regions and extract information
    info_phasing = data.frame()
    for (r in all_regions){
        # extract info
        tmp_phasing = extractHaploSize(treat_phasing[which(treat_phasing$ID == r),])
        tmp = extractHaploSize(treat[which(treat$ID == r),])
        # change colnames
        colnames(tmp_phasing)[2:ncol(tmp_phasing)] = paste(colnames(tmp_phasing)[2:ncol(tmp_phasing)], 'phasing', sep='_')
        # combine and add to dataframe
        tmp_combined = merge(tmp, tmp_phasing, by = 'sample', all=T)
        # add region
        tmp_combined$region = r
        # then add to dataframe
        info_phasing = rbind(info_phasing, tmp_combined)
    }

    # Basica correlation
    # Short allele
    cor(info_phasing$short_allele, info_phasing$short_allele_phasing, use='complete.obs')
    # Long allele
    cor(info_phasing$long_allele, info_phasing$long_allele_phasing, use='complete.obs')
    # plot
    short_plot = ggplot(info_phasing, aes(x = short_allele, y = short_allele_phasing)) + geom_point(stat = 'identity') + geom_smooth() + xlab('Short allele (No phasing)') + ylab('Short allele (Phasing)') + ggtitle('Short Allele')
    long_plot = ggplot(info_phasing, aes(x = long_allele, y = long_allele_phasing)) + geom_point(stat = 'identity') + geom_smooth() + xlab('Long allele (No phasing)') + ylab('Long allele (Phasing)') + ggtitle('Long Allele')
    combined_plot = ggarrange(plotlist=list(short_plot, long_plot), nrow=1)