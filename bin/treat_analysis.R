# Libraries: check if the required packages are installed, and if not, install them
  packages <- c("data.table", "stringr", "argparse")
  check_install_packages <- function(packages) {
      installed_packages <- installed.packages()
      for (package in packages){
        if (package %in% installed_packages){
          suppressPackageStartupMessages(library(package, character.only = TRUE))
        } else {
          install.packages(package)
          suppressPackageStartupMessages(library(package, character.only = TRUE))
        }
      }
    }
  check_install_packages(packages)


# Functions
    # Pipeline to analyze repeats
    analysisPipeline <- function(rs_vcf, out_dir, out_name, region, mad_thr){
      cat('\n')
      # Disable warnings
      defaultW <- getOption("warn"); options(warn = -1)
      # Read VCF
      vcf = data.frame()
      for (v in rs_vcf){
        tmp = readVCF(v, region)
        if (nrow(vcf) == 0){
          vcf = tmp
        } else {
          columns_to_select = c('ID', colnames(tmp)[10:ncol(tmp)])
          tmp = tmp[, ..columns_to_select]
          vcf = merge(vcf, tmp, by = 'ID', all=T)
        }
      }
      # Identify the regions to be plotted
      all_regions = unique(vcf$ID)
      # Create an empty dataframe for the results
      outlier_results = data.frame()
      # The iterate over all regions
      for (r in all_regions){
        # Extract sizes
        vcf_info = extractHaploSize(vcf[which(vcf$ID == r),])
        # Extract Reference and add it to the data
        vcf_info_withRef = extractReference(vcf, r, vcf_info)
        # Outlier analysis
        vcf_info_withRef = outlierAnal(vcf_info_withRef, r, mad_thr = mad_thr)
        # Add to dataframe
        outlier_results = rbind(outlier_results, vcf_info_withRef)
      }
      # Define the path to the output file
      output_file = file.path(out_dir, out_name)
      write.table(outlier_results, output_file, sep="\t", quote=F, row.names=F)
      cat('\n')
    }

    # Single functions
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

    # Function to print summary of the run
    summaryRun <- function(rs_vcf, out_dir, out_name, region, mad_thr){
        cat('\n********************')
        cat('\n** TREAT analysis **')
        cat('\n*** Arguments:')
        cat(paste0('\n*** Input VCF: ', paste(rs_vcf, collapse = ',')))
        cat(paste0('\n*** Region(s): ', paste(region, collapse = ', ')))
        cat(paste0('\n*** MAD threshold: ', mad_thr))
        cat(paste0('\n*** Output directory: ', out_dir))
        cat(paste0('\n*** Output name: ', out_name))
        return('\n*** Analysis started!')
    }

    # Function to extract sizes
    extractHaploSize = function(vcf){
      # restrict to sample's columns
      sub = vcf[, 10:ncol(vcf)]
      # take the transpose of this
      sub_trans = t(sub)
      # split columns
      sub_trans_df = data.frame(str_split_fixed(sub_trans[, 1], ';', 6))
      sub_trans_df$sample = rownames(sub_trans)
      # rename columns
      colnames(sub_trans_df) = c('qc', 'haplo', 'motifs', 'copies', 'copies_reference', 'depth', 'sample')
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
            sub_trans_df$short_allele_cn[x] = as.numeric(str_split_fixed(sub_trans_df$copies[x], '\\|', 2)[, 2]); sub_trans_df$long_allele_cn[x] = as.numeric(str_split_fixed(sub_trans_df$copies[x], '\\|', 2)[, 1])
            # copy number of reference
            sub_trans_df$short_allele_cn_ref[x] = as.numeric(str_split_fixed(sub_trans_df$copies_reference[x], '\\|', 2)[, 2]); sub_trans_df$long_allele_cn_ref[x] = as.numeric(str_split_fixed(sub_trans_df$copies_reference[x], '\\|', 2)[, 1])
          } else {
            # haplotype size
            sub_trans_df$short_allele[x] = sub_trans_df$haplo_h1[x]; sub_trans_df$long_allele[x] = sub_trans_df$haplo_h2[x]
            # motifs
            sub_trans_df$short_allele_motif[x] = str_split_fixed(sub_trans_df$motifs[x], '\\|', 2)[, 1]; sub_trans_df$long_allele_motif[x] = str_split_fixed(sub_trans_df$motifs[x], '\\|', 2)[, 2]
            # copy number
            sub_trans_df$short_allele_cn[x] = as.numeric(str_split_fixed(sub_trans_df$copies[x], '\\|', 2)[, 1]); sub_trans_df$long_allele_cn[x] = as.numeric(str_split_fixed(sub_trans_df$copies[x], '\\|', 2)[, 2])
            # copy number of reference
            sub_trans_df$short_allele_cn_ref[x] = as.numeric(str_split_fixed(sub_trans_df$copies_reference[x], '\\|', 2)[, 1]); sub_trans_df$long_allele_cn_ref[x] = as.numeric(str_split_fixed(sub_trans_df$copies_reference[x], '\\|', 2)[, 2])
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

    # Function to extract reference sizes
    extractReference = function(vcf, r, vcf_info){
      # get reference size
      ref_size = as.numeric(vcf$REF[which(vcf$ID == r)])
      # then get reference info
      ref_motif = str_split_fixed(vcf$INFO[which(vcf$ID == r)], ';', 2)[, 1]
      ref_copies = str_split_fixed(vcf$INFO[which(vcf$ID == r)], ';', 2)[, 2]
      # create dataframe
      ref_df = data.frame(sample = 'GRCh38', short_allele = ref_size, long_allele = ref_size, short_allele_motif = ref_motif, long_allele_motif = ref_motif, short_allele_cn = ref_copies, long_allele_cn = ref_copies, short_allele_cn_ref = ref_copies, long_allele_cn_ref = ref_copies)
      # join with sample's data
      vcf_info = rbind(vcf_info, ref_df)
      return(vcf_info)
    }

    # Function to find outliers based on distance
    findOutliers_distance = function(size, mad_thr){
      mean_size <- mean(size, na.rm=T)
      covariance_matrix <- var(size, na.rm=T)
      # Compute Mahalanobis distance
      mahalanobis_dist <- abs(size - mean_size) / sqrt(covariance_matrix)
      # Identify outliers
      outliers <- which(mahalanobis_dist > mad_thr)
      # Compute p-values using chi-squared distribution
      p_values <- 1 - pchisq(mahalanobis_dist^2, df = 1)
      outliers <- which(p_values < 0.05)
      return(list(outliers, p_values))
    }

    # Function to find thresholds for outlier analysis
    findOutliers = function(size, mad_thr){
        # Compute the modified Z-scores
        med <- median(size, na.rm=TRUE)
        mad <- mad(size, constant = 1.4826, na.rm=TRUE)
        modified_z <- abs((size - med) / mad)
        # Identify the outliers
        outliers = which(modified_z > mad_thr)
        # Calculate the p-value using modified z-score method
        pvalues <- sapply(size[outliers], function(x) {
            mod_z <- (x - med) / (1.4826 * mad)
            p <- 2 * pnorm(-abs(mod_z))
            return(p)
        })
        return(list(outliers, pvalues))
    }

    # Function to run the outlier analysis
    outlierAnal = function(vcf_info_withRef, region, mad_thr){
        # score outliers based on the short allele
        #short_allele_res = findOutliers(vcf_info_withRef$short_allele, mad_thr)
        short_allele_res = findOutliers_distance(vcf_info_withRef$short_allele, mad_thr)
        # score outliers based on the long allele
        #long_allele_res = findOutliers(vcf_info_withRef$long_allele, mad_thr)
        long_allele_res = findOutliers_distance(vcf_info_withRef$long_allele, mad_thr)
        # score outliers based on the sum of the alleles
        #sum_allele_res = findOutliers(vcf_info_withRef$short_allele + vcf_info_withRef$long_allele, mad_thr)
        sum_allele_res = findOutliers_distance(vcf_info_withRef$short_allele + vcf_info_withRef$long_allele, mad_thr)
        # put results in the dataframe for short allele
        vcf_info_withRef$outliers_short_allele = ifelse(rownames(vcf_info_withRef) %in% short_allele_res[[1]], 'yes', NA)
        vcf_info_withRef$outliers_short_allele_pvalue = ifelse(rownames(vcf_info_withRef) %in% short_allele_res[[1]], short_allele_res[[2]][match(rownames(vcf_info_withRef), as.character(short_allele_res[[1]]))], NA)
        # put results in the dataframe for long allele
        vcf_info_withRef$outliers_long_allele = ifelse(rownames(vcf_info_withRef) %in% long_allele_res[[1]], 'yes', NA)
        vcf_info_withRef$outliers_long_allele_pvalue = ifelse(rownames(vcf_info_withRef) %in% long_allele_res[[1]], long_allele_res[[2]][match(rownames(vcf_info_withRef), as.character(long_allele_res[[1]]))], NA)
        # put results in the dataframe for the sum of the alleles
        vcf_info_withRef$outliers_sum_allele = ifelse(rownames(vcf_info_withRef) %in% sum_allele_res[[1]], 'yes', NA)
        vcf_info_withRef$outliers_sum_allele_pvalue = ifelse(rownames(vcf_info_withRef) %in% sum_allele_res[[1]], sum_allele_res[[2]][match(rownames(vcf_info_withRef), as.character(sum_allele_res[[1]]))], NA)
        # Remove NAs and add region
        vcf_info_withRef <- vcf_info_withRef[!(is.na(vcf_info_withRef$outliers_short_allele) & is.na(vcf_info_withRef$outliers_long_allele) & is.na(vcf_info_withRef$outliers_sum_allele)), ]
        if (nrow(vcf_info_withRef) >0){ vcf_info_withRef$region = region; vcf_info_withRef$mad_threshold = mad_thr }
        return(vcf_info_withRef)
    }

# Arguments definition
    parser <- ArgumentParser()
    # add arguments: --reads_spannning is the VCF file of the output of read_spanning_analysis
    parser$add_argument("--vcf", default = 'None', help = "VCF file output of TREAT. Multiple files should be comma-separated.")
    # add arguments: --out is the output directory
    parser$add_argument("--out", default = './', help = "Output directory where output will be placed. Default is the current directory.")
    # add arguments: --outname is the name of the output file
    parser$add_argument("--outname", default = 'treat_analysis_output.txt', help = "Name of the plot. Default name is treat_analysis_output.txt")
    # add arguments: --region is the name of the region to plot
    parser$add_argument("--region", default = 'all', help = "Name of the region to plot. By default, all regions will be analyzed.")
    # add arguments: --madThr is the value to call outliers
    parser$add_argument("--madThr", default = 3, help = "Median Absolute deviation value used to call outliers.")

# Read arguments
    args <- parser$parse_args()
    # vcf of read-spanning analysis
    rs_vcf = args$vcf; rs_vcf = unlist(strsplit(rs_vcf, ','))
    # rs_vcf = 'samples_genotypes.vcf'
    # output directory
    out_dir = args$out
    # out_dir = './'
    # output name
    out_name = args$outname
    # out_name = 'treat_analysis_output.txt'
    # region
    region = args$region
    # region = 'all'
    # madThr
    mad_thr = as.numeric(args$madThr)

# Check arguments and stop if no input data is provided
    run = 'false'
    if ((length(rs_vcf) == 1) && (rs_vcf == 'None')){ stop("Input error: Missing input file(s)!!") } else if (!is.numeric(mad_thr)){ "The --madThr should be numeric" } else { run = 'true'}
# Main
if (run == 'true'){
    # Split regions if it is not all
    if (region != 'all'){ region = unlist(strsplit(region, ',')) }
    # Create directory if not present
    if (!dir.exists(out_dir)){ system(paste0('mkdir ', out_dir)) }
    # Print summary of the run
    x = summaryRun(rs_vcf, out_dir, out_name, region, mad_thr)

    # Pipeline to plot repeats
    analysisPipeline(rs_vcf, out_dir, out_name, region, mad_thr)
  }
