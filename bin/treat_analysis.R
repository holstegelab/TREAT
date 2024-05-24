# Libraries: check if the required packages are installed, and if not, install them
  packages <- c("data.table", "stringr", "argparse", "parallel")
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
      # Identify the regions to be analysed
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
      cat(paste0('\n** Analysis finished: output table is ', output_file))
      cat('\n\n')
    }

    # Function for the case-control analysis
    caseControlPipeline <- function(rs_vcf, out_dir, out_name, region, cc_labels, cpu){
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
      # Compare labels with sample's names
      vcf_sample_labels = compareLabelsVCF(cc_labels, vcf)
      # Identify the regions to be analysed
      all_regions = unique(vcf$ID)
      # Do association
      assoc_results = rbindlist(mclapply(all_regions, casecontrolAssoc, vcf = vcf, vcf_sample_labels = vcf_sample_labels, mc.cores = cpu))
      # Define the path to the output file
      output_file = file.path(out_dir, out_name)
      write.table(assoc_results, output_file, sep="\t", quote=F, row.names=F)
      cat(paste0('\n** Analysis finished: output table is ', output_file))
      cat('\n\n')
    }

    # Single functions
    # Function to do case-control association 
    casecontrolAssoc <- function(r, vcf, vcf_sample_labels){
      cat(paste0('**** Analysis of ', r, '\n'))
      tryCatch({
        # Extract sizes
        vcf_info = extractHaploSize(vcf[which(vcf$ID == r),])
        # Add the sum of alleles
        vcf_info$join_allele = vcf_info$short_allele + vcf_info$long_allele
        # Add phenotype
        vcf_info_ph = merge(vcf_info, vcf_sample_labels, by.x = 'sample', by.y = 'V1')
        colnames(vcf_info_ph)[which(colnames(vcf_info_ph) == 'V2')] = 'pheno'
        # Model: short, long and joint alleles
        m_short = glm(pheno ~ short_allele, data = vcf_info_ph, family = 'binomial')
        m_long = glm(pheno ~ long_allele, data = vcf_info_ph, family = 'binomial')
        m_joint = glm(pheno ~ join_allele, data = vcf_info_ph, family = 'binomial')
        # Combine results in a table
        tmp = data.frame(region = rep(r, 3), n_cases = rep(nrow(vcf_info_ph[which(vcf_info_ph$pheno == 1),]), 3), n_controls = rep(nrow(vcf_info_ph[which(vcf_info_ph$pheno == 0),]), 3), model = c('short', 'long', 'joint'), 
          beta = c(summary(m_short)$coefficient[2, 'Estimate'], summary(m_long)$coefficient[2, 'Estimate'], summary(m_joint)$coefficient[2, 'Estimate']), 
          se = c(summary(m_short)$coefficient[2, 'Std. Error'], summary(m_long)$coefficient[2, 'Std. Error'], summary(m_joint)$coefficient[2, 'Std. Error']), 
          pvalue = c(summary(m_short)$coefficient[2, 'Pr(>|z|)'], summary(m_long)$coefficient[2, 'Pr(>|z|)'], summary(m_joint)$coefficient[2, 'Pr(>|z|)']))
      }, error = function(e) {
        tmp = data.frame(region = rep(r, 3), n_cases = rep(NA, 3), n_controls = rep(NA, 3), model = c('short', 'long', 'joint'), beta = rep(NA, 3), se = rep(NA, 3), pvalue = rep(NA, 3))
      })
      return(tmp)
    }

    # Function to check labels for the case-control analysis
    checkLabels <- function(labels){
      cat('\n\n')
      # Disable warnings
      defaultW <- getOption("warn"); options(warn = -1)
      tryCatch({
        # Attempt to read the file
        cc_labels <- read.table(labels, h=F)
        # check case-control label
        if (is.numeric(cc_labels$V2)){
          cc_labels$V2[which(cc_labels$V2 == min(cc_labels$V2))] = 0
          cc_labels$V2[which(cc_labels$V2 == max(cc_labels$V2))] = 1
          if (length(unique(cc_labels$V2)) == 2){
            cat('** Case-control label is numeric and binary\n')
            cat(paste0('** There are ', nrow(cc_labels[which(cc_labels$V2 == 1),]), ' cases, and ', nrow(cc_labels[which(cc_labels$V2 == 0),]), ' controls\n'))
            return(cc_labels)
          } else if (length(unique(cc_labels$V2)) == 1){
            cat('** Case-control label is NOT binary: only 1 label found\n')
            cat('** Stopping analysis! Please check the case-control labels, and make sure they are binary\n')
            return(NULL)
          } else if (length(unique(cc_labels$V2)) >2){
            cat('** Case-control label is NOT binary: >2 labels found\n')
            cat('** Stopping analysis! Please check the case-control labels, and make sure they are binary\n')
            return(NULL)
          }
        } else {
          if (length(unique(cc_labels$V2)) == 2){
            lbs = unique(cc_labels$V2)
            cc_labels$V2[which(cc_labels$V2 == lbs[1])] = 0
            cc_labels$V2[which(cc_labels$V2 == lbs[2])] = 1
            cat('** Case-control label is not numeric but it is binary\n')
            cat('** Converting ', lbs[1], ' to controls, and ', lbs[2], ' to cases\n')
            cat(paste0('** In total, there are ', nrow(cc_labels[which(cc_labels$V2 == 1),]), ' cases, and ', nrow(labels_vcf[which(labels_vcf$V2 == 0),]), ' controls\n'))
            return(cc_labels)
          } else if (length(unique(cc_labels$V2)) == 1){
            cat('** Case-control label is NOT binary: only 1 label found\n')
            cat('** Stopping analysis! Please check the case-control labels, and make sure they are binary\n')
            return(NULL)
          } else if (length(unique(cc_labels$V2)) >2){
            cat('** Case-control label is NOT binary: >2 labels found\n')
            cat('** Stopping analysis! Please check the case-control labels, and make sure they are binary\n')
            return(NULL)
          }
        }
      }, error = function(e) {
        # Handle the error if the file doesn't exist
        return(NULL)
      })
    }

    # Function to compare labels to VCF sample names
    compareLabelsVCF <- function(labels, vcf){
      samples_vcf = colnames(vcf)[10:ncol(vcf)]
      # check overlap
      labels_vcf = cc_labels[which(cc_labels$V1 %in% samples_vcf),]
      cat(paste0('** The VCF file includes ', length(samples_vcf), ' samples. Of these, ', nrow(labels_vcf), '/', length(samples_vcf), ' have matching case-control label\n'))
      cat('** Data looks good! starting with association using logistic regression model.\n')
      return(labels_vcf)      
    }

    # Function to read VCF, and restrict to region of interest
    readVCF <- function(rs_vcf, region){
      # read vcf
      d = fread(rs_vcf, h=T, skip = '##', sep="\t")
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
    summaryRun <- function(rs_vcf, out_dir, out_name, region, mad_thr, analysis, labels){
        cat('\n********************')
        cat('\n** TREAT analysis **')
        cat('\n*** Arguments:')
        cat(paste0('\n*** Analysis: ', paste(analysis, collapse = ',')))
        cat(paste0('\n*** Input VCF: ', paste(rs_vcf, collapse = ',')))
        cat(paste0('\n*** Region(s): ', paste(region, collapse = ', ')))
        cat(paste0('\n*** Output directory: ', out_dir))
        cat(paste0('\n*** Output name: ', out_name))
        if (analysis == 'outlier'){
          cat(paste0('\n*** MAD threshold: ', mad_thr))
        } else if (analysis == 'case-control'){
          cat(paste0('\n*** Case-control labels: ', labels))
        }
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
    # add arguments: --analysis is type of analysis to perform: [outlier / case-control]
    parser$add_argument("--analysis", default = 'None', help = "Type of analysis to perform: [outlier / case-control]")
    # add arguments: --labels is the file containing the labels of the samples for comparison
    parser$add_argument("--labels", default = 'None', help = "File including the case-control labels for each sample. This should be a tab-separated file with 2 columns and no header, reporting sample name and a binary label.")
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
    # add arguments: --cpu is the number of cpus to use
    parser$add_argument("--cpu", default = 1, help = "Number of CPUs to use (Default: 1).")

# Read arguments
    args <- parser$parse_args()
    # analysis type
    analysis = args$analysis
    # eventually the case-control labels
    labels = args$labels
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
    # cpu
    cpu = as.numeric(args$cpu)
  
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
    x = summaryRun(rs_vcf, out_dir, out_name, region, mad_thr, analysis, labels)

    # Pipeline to do outlier analysis
    if (analysis == 'outlier'){
      # run outlier analysis pipeline
      analysisPipeline(rs_vcf, out_dir, out_name, region, mad_thr)
    } else if (analysis == 'case-control'){
      # check the labels file
      cc_labels = checkLabels(labels)
      if (is.null(cc_labels)){
        stop("Error: The labels file does not exist or is unreadable. This should be a tab-separated file with 2 columns and no header, reporting sample name and a binary label.")
      } else {
        # then run the case-control analysis if label file is ok
        caseControlPipeline(rs_vcf, out_dir, out_name, region, cc_labels, cpu)
      }
    }
  }
