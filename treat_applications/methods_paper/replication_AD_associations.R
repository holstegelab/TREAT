
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

# Also try to replicate the STR from the ADSP paper
ad_vcf = data.table::fread('/project/holstegelab/Share/nicco/paper_treat/20240206_final/replication_adsp/ad_asm/sample.vcf.gz', h=T, stringsAsFactors=F, sep="\t")
chc_vcf = data.table::fread('/project/holstegelab/Share/nicco/paper_treat/20240206_final/replication_adsp/chc_asm/sample.vcf.gz', h=T, stringsAsFactors=F, sep="\t")

res_list = list()
all_regions = unique(ad_vcf$ID)
for (r in all_regions){
  # Extract sizes
  tmp_ad_vcf = extractHaploSize(ad_vcf[which(ad_vcf$ID == r),])
  tmp_chc_vcf = extractHaploSize(chc_vcf[which(chc_vcf$ID == r),])
  # combine sizes
  tmp_ad_vcf$combined = tmp_ad_vcf$short_allele + tmp_ad_vcf$long_allele
  tmp_chc_vcf$combined = tmp_chc_vcf$short_allele + tmp_chc_vcf$long_allele
  # combine
  all = rbind(tmp_ad_vcf[, c('sample', 'combined', 'short_allele', 'long_allele')], tmp_chc_vcf[, c('sample', 'combined', 'short_allele', 'long_allele')])
  all$pheno = c(rep(1, nrow(tmp_ad_vcf)), rep(0, nrow(tmp_chc_vcf)))
  # model
  model_asm_combined = glm(pheno ~ combined, data = all, family = 'binomial')
  model_asm_short = glm(pheno ~ short_allele, data = all, family = 'binomial')
  model_asm_long = glm(pheno ~ long_allele, data = all, family = 'binomial')
  cat('** REGION: ', r, '\n')
  cat('** MODEL Combined asm beta: ', summary(model_asm_combined)$coefficient[2, 1], '\n')
  cat('** MODEL asm p: ', summary(model_asm_combined)$coefficient[2, 4], '\n')
  cat('** MODEL Short asm beta: ', summary(model_asm_short)$coefficient[2, 1], '\n')
  cat('** MODEL asm p: ', summary(model_asm_short)$coefficient[2, 4], '\n')
  cat('** MODEL Long asm beta: ', summary(model_asm_long)$coefficient[2, 1], '\n')
  cat('** MODEL asm p: ', summary(model_asm_long)$coefficient[2, 4], '\n\n\n\n')
  tmp = list(all, model_asm_combined, model_asm_long, model_asm_short)
  res_list[[length(res_list) + 1]] = tmp
}




