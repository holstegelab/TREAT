# open the vcf file and extract allele sizes

# libraries
import os
import pandas as pd
#from collections import Counter
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
#from scipy.stats import binom_test
#import matplotlib.pyplot as plt

# main
# find files on which to work
file_dir = '/project/holstegelab/Share/nicco/paper_treat/20240206_final'
filelist_asm = [x.rstrip() for x in list(os.popen('find %s -name "*vcf*" | grep "_asm"' %(file_dir)))]

# iterate across files
all_df = pd.DataFrame()
for sample in filelist_asm:
    print(sample)
    # open file
    df_sample = pd.read_csv(sample, comment='#', sep='\t', header=None)
    # grep header from the file
    header = list(os.popen('zgrep -w "%s" %s' %('#CHROM', sample)))[0].rstrip().split('\t')
    # extract sample name
    sample_name = header[-1]
    # assign header
    df_sample.columns = header
    # extract allele sizes
    df_sample['Alleles'] = df_sample[sample_name].str.split(';').str[1]
    df_sample['Depth'] = df_sample[sample_name].str.split(';').str[-1]
    df_sample = df_sample[df_sample['Alleles'] != 'NA|NA'].copy()
    df_sample['Alleles_num'] = df_sample['Alleles'].str.split('|').apply(lambda x: [float(i) for i in x])
    df_sample[sample_name + '_join'] = df_sample['Alleles_num'].apply(sum)
    df_sample_sb = df_sample[['ID', sample_name + '_join']]
    if all_df.shape[0] == 0:
        all_df = df_sample[['#CHROM', 'POS', 'ID', 'REF', sample_name + 'join']]
    else:
        all_df = pd.merge(all_df, df_sample_sb, left_on = 'ID', right_on = 'ID')

# with a unique dataset, check the size
all_df.shape
# then add coefficient of variation to each variable
all_df['CV'] = (all_df.std(axis=1) / all_df.mean(axis=1)) * 100
# then take top x% with highest coefficient of variation
sorted_cv = sorted(list(all_df['CV']), reverse=True)
threshold = sorted_cv[int(len(sorted_cv) * 0.20)]
# subset the original df
subset_df = all_df[all_df['CV'] >= threshold]
subset_df = subset_df.drop('CV', axis=1)
# transpose
sample_names = list(subset_df.columns)[1::]
subset_df_t = subset_df.T
subset_df_t.columns = subset_df_t.iloc[0]
subset_df_t = subset_df_t.iloc[1:].reset_index(drop=True)
subset_df_t.index = sample_names
subset_df_t.head()
# pca
scaler = StandardScaler()
# Normalize the dataframe using the StandardScaler object
subset_df_t_norm = scaler.fit_transform(subset_df_t)
# Create a PCA object with 2 components
pca = PCA(n_components=2)
# Fit and transform the dataframe using PCA
pca_df = pca.fit_transform(subset_df_t_norm)
# Create a new dataframe with the PCA results
pca_df = pd.DataFrame(data=pca_df, columns=['PC1', 'PC2'])
pca_df.index = sample_names
pca_df.to_csv('20240217_pca_results_QC_normalized.csv')
explained_variance_pc1 = pca.explained_variance_ratio_[0]
explained_variance_pc2 = pca.explained_variance_ratio_[1]

### plot in R
d = data.table::fread('20240217_pca_results_QC_normalized.csv', h=T, sep=",")
d$Sample = stringr::str_split_fixed(d$V1, '_', 2)[, 1]
info = data.table::fread('samples_info_hprc.txt', h=T, fill=T)
d_info = merge(d, info, by = 'Sample')
dim(d_info)
library(ggplot2)
pdf('20240217_population_plot.pdf', height=7, width=7)
ggplot(data = d_info, aes(x = PC1, y = PC2, color = Superpopulation)) + geom_point(stat = 'identity', size = 5, alpha = 0.7) + scale_color_manual(values = RColorBrewer::brewer.pal(n = 10, name = "Paired"))
dev.off()

### in R, also check the snps
kg = read.table('/project/holstegelab/Share/gwas_array/1000Genome_common_snps/additional_samples/chrAll_1KG.psam', h=F, stringsAsFactors=F)
hpc = read.table('samples_info_hprc.txt', h=T)
common_samples = kg[which(kg$V1 %in% hpc$Sample),]
# 40 are present, 7 are missing
write.table(common_samples$V1, 'common_samples_KG_hpc.txt', quote=F, row.names=F)
system('plink2 --pfile /project/holstegelab/Share/gwas_array/1000Genome_common_snps/additional_samples/chrAll_1KG --keep common_samples_KG_hpc.txt --maf 0.10 --thin-count 150000 --export A --out KG_SNPs_30k_random')

# in Python, do PCA in the same way as with the TR
df = pd.read_csv('KG_SNPs_30k_random.raw', sep='\t')
df = df.drop(columns=['FID'])
df = df.drop(columns=['PAT'])
df = df.drop(columns=['MAT'])
df = df.drop(columns=['PHENOTYPE'])
df = df.drop(columns=['SEX'])
df.set_index('IID', inplace=True)
# pca
scaler = StandardScaler()
# Normalize the dataframe using the StandardScaler object
subset_df_t_norm = scaler.fit_transform(df)
# Create a PCA object with 2 components
pca = PCA(n_components=2)
# Fit and transform the dataframe using PCA
pca_df = pca.fit_transform(subset_df_t_norm)
# Create a new dataframe with the PCA results
pca_df = pd.DataFrame(data=pca_df, columns=['PC1', 'PC2'])
pca_df.index = df.index
pca_df.to_csv('20240313_pca_results_QC_normalized_SNPs.csv')
explained_variance_pc1 = pca.explained_variance_ratio_[0]
explained_variance_pc2 = pca.explained_variance_ratio_[1]

### plot in R
d = data.table::fread('20240313_pca_results_QC_normalized_SNPs.csv', h=T, sep=",")
d$Sample = stringr::str_split_fixed(d$V1, '_', 2)[, 1]
info = data.table::fread('samples_info_hprc.txt', h=T, fill=T)
d_info = merge(d, info, by.x = 'IID', by.y = 'Sample')
dim(d_info)
library(ggplot2)
pdf('20240313_population_plot_SNPs.pdf', height=7, width=7)
ggplot(data = d_info, aes(x = PC1, y = PC2, color = Superpopulation)) + geom_point(stat = 'identity', size = 5, alpha = 0.7) + scale_color_manual(values = RColorBrewer::brewer.pal(n = 10, name = "Paired")) + xlab('PC1 (14%)') + ylab('PC2 (4%)') + ggtitle('PCA of 30544 random common SNPs (40 overlapping individuals)')
dev.off()