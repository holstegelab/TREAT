# Libraries: check if the required packages are installed, and if not, install them
  packages <- c("data.table", "stringr", "argparse", "ggplot2", "berryFunctions", "dendextend", "dplyr")
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
  # Pipelines: set of functions
    # Pipeline to plot repeats
    plotRepeatsComplete <- function(rs_vcf, out_dir, out_name, plotFormat, custom_colors, region, path){
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
      # The iterate over all regions
      for (r in all_regions){
        # Extract sizes
        vcf_info = extractHaploSize(vcf[which(vcf$ID == r),])
        # Extract Reference and add it to the data
        vcf_info_withRef = extractReference(vcf, r, vcf_info)
        # Clustering
        clustering_info = cluster_TR(vcf_info_withRef)
        # Plot name
        plt_name = paste0(out_name, '_', r, '.', plotFormat)
        plt_name_af = paste0(out_name, '_', r, '_Freq.', plotFormat)
        plotname = file.path(out_dir, plt_name)
        plotname_af = file.path(out_dir, plt_name_af)
        # Plot tandem repeat
        pdf(plotname, height = 10, width = 12)
        plotComplete(vcf_info_withRef, clustering_info, custom_colors = custom_colors, region = r, path = path)
        dev.off()
        # Plot allele frequency
        pdf(plotname_af, height = 7, width = 12)
        print(plotAlleleFrequency(vcf_info_withRef, custom_colors = custom_colors, region = r))
        dev.off()
        cat('\nPlots done --> ', plotname, ' + ', plotname_af)
      }
      cat('\n')
    }

  # Single functions used in the pipelines
    # Function to read VCF, and restrict to region of interest
    readVCF <- function(rs_vcf, region){
      # read vcf
      d = fread(rs_vcf, h=T, sep="\t")
      # restrict to region of interest
      if (length(region) == 1 && (region == 'all')){ sub = d } else { sub = d[which(d$ID %in% region),] }
      # check if region existed
      if (nrow(sub) >0){
        return(sub)
      } else {
        stop('The region you provided does not exists! Stopping.')
      }
    }

    # Function to plot allele frequency
    plotAlleleFrequency <- function(vcf_info_withRef, custom_colors, region){
        # get all alleles
        all_alleles = na.omit(c(vcf_info_withRef$short_allele, vcf_info_withRef$long_allele))
        # bin size based on fixed value of 20bp
        bin_size_fixed <- 10
        # create bins based on the fixed values
        intervals = seq(from = min(all_alleles), to = max(all_alleles), by = 2*bin_size_fixed)
        intervals = c(intervals, intervals[length(intervals)] + 20)
        # cut in bins
        bins <- data.frame(table(cut(all_alleles, breaks = intervals, include.lowest = TRUE)))
        bins$Size = intervals[1:length(intervals)-1] + 10
        # clean the empty bins
        bins = bins[which(bins$Freq >0),]
        # calculate frequency
        bins$fraction = bins$Freq / sum(bins$Freq)
        # adjust label for the bin
        bins$Var1 = str_replace_all(bins$Var1, c(',' = '-', '\\(' = '', '\\]' = '', '\\[' = ''))
        # sort by frequency
        bins = bins[order(bins$Size),]
        bins$Var1 = factor(bins$Var1, levels = bins$Var1)
        bins$x_pos = seq(1, nrow(bins))
        # add labels
        bins$label = paste0(bins$Size-10, '-', bins$Size+10)
        # plot
        plt = ggplot(bins, aes(x = x_pos, y = fraction, fill=Size)) + geom_bar(stat='identity') + scale_x_continuous(breaks = bins$x_pos, labels = bins$label) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), axis.text.y = element_text(size=14), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey90"), plot.title = element_text(size = 18), legend.text = element_text(size = 14), legend.margin = margin(t = 0, r = 5, b = 0, l = 5), legend.key.size = unit(1, "cm"), legend.key.width = unit(0.5, "cm"), legend.key.height = unit(1, "cm"), legend.title = element_text(size = 14)) + xlab('Bin size in base pairs') + ylab('Allele size frequency') + ggtitle(paste('Allele frequency of', region, 'in', nrow(vcf_info_withRef), 'samples')) + guides() + scale_fill_gradient(low = 'deepskyblue2', high = 'red', limits = c(min(all_alleles), max(all_alleles)))
      return(plt)
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

  # Function to do hierarchical clustering on the TR sizes and eventually motifs
    cluster_TR <- function(vcf_info_withRef){
      # create data for clustering
      sb_data = vcf_info_withRef[, c('short_allele', 'long_allele')]
      labs = vcf_info_withRef$sample
      # convert NAs
      sb_data$short_allele[is.na(sb_data$short_allele)] = -10000
      sb_data$long_allele[is.na(sb_data$long_allele)] = -10000
      # also add a combined size that maybe helps for the clustering
      sb_data$combined = sb_data$short_allele + sb_data$long
      # then do hierarchical clustering -- canberra, binary, minkowski
      hr <- as.dendrogram(hclust(dist(sb_data, method = 'canberra'), method="complete", members = NULL))
      hc = hclust(dist(sb_data, method = 'canberra'), method="complete", members = NULL)
      ordered_labs = labs[hc$order]
      # return dendrogram + ordered labels
      res = list(hr, ordered_labs)
      return(res)
    }

    # Function to plot TR of the samples -- plotted motif is the consensus motif
    plotComplete <- function(vcf_info_withRef, clustering_info, custom_colors = 'None', region = r, path){
      # Find info about position of interest
      chrom = str_split_fixed(region, ':', 2)[, 1]; start = str_split_fixed(str_split_fixed(region, ':', 2)[, 2], '-', 2)[, 1]; stop = str_split_fixed(str_split_fixed(region, ':', 2)[, 2], '-', 2)[, 2]
      
      # Define maximum and minimum allele sizes
      min_len <- floor(min(na.omit(vcf_info_withRef$short_allele)) - min(na.omit(vcf_info_withRef$short_allele))*0.15); if (min_len == 0){ min_len = -1}
      max_len <- ceiling(max(na.omit(vcf_info_withRef$long_allele)) + max(na.omit(vcf_info_withRef$long_allele))*0.15)
        
      # Parse the chromosomal band
      chr_bands <- fread(paste0(path, "/hg38_cytogenetic_bands.txt"), h=T)
      tmp_band <- chr_bands[which(chr_bands$`#chrom` == chrom),]
      # assign colors depending on transcription levels
      tmp_band$col <- "white"; tmp_band$col[which(tmp_band$gieStain == "gpos25")] <- "grey90";tmp_band$col[which(tmp_band$gieStain == "gpos50")] <- "grey60"; tmp_band$col[which(tmp_band$gieStain == "gpos75")] <- "grey40"; tmp_band$col[which(tmp_band$gieStain == "gpos100")] <- "black"; tmp_band$col[which(tmp_band$gieStain == "acen")] <- "deepskyblue3"
          
      # Layout of the plot: chromosome band and location at the top -- main plot with samples' size of TR -- on left side, hirarchical clustering of the samples -- on the right side, motif information
      layout(matrix(c(1,2,2,2,2,2,2,2, 1,3,3,3,3,3,3,3, 1,3,3,3,3,3,3,3, 1,4,4,4,4,4,4,4, 1,4,4,4,4,4,4,4), nrow = 8, ncol = 5, byrow = F))
        
      # First plot: chromosomal band
        # define width of chromosome and graphical parameters
        w = 0.15; par(mar=c(0, 6, 2, 0))
        # background plot
        plot(0,0, pch=16, col="white", xlim=c(0, max(tmp_band$chromEnd) + max(tmp_band$chromEnd)*0.20), ylim=c(0, 1), bty="n", xaxt="none", yaxt="none", ylab="", xlab="", xaxs = "i", yaxs = "i")
        # draw chromosomal bands -- first rounded rectangle until centromere -- need to find position of the centromere
        index_centromere = which(tmp_band$gieStain == "acen")
        if (length(index_centromere) >0){
          # draw chromosomal bands before centromere and rectangle around it
          for (i in 1:(index_centromere[1] - 1)){ rect(xleft = tmp_band$chromStart[i], ybottom = 0.5 - w, xright = tmp_band$chromEnd[i], ytop = 0.5 + w, col = tmp_band$col[i], border = NA) }
          roundedRect(xleft = 0, ybottom = 0.5-w, xright = max(tmp_band$chromEnd[(index_centromere[1])]), ytop = 0.5 + w, lwd=3, rounding = 0.20, xpd=T)
          # draw chromosomal bands after centromere and rectangle around it
          for (i in (index_centromere[2] + 1):nrow(tmp_band)){ rect(xleft = tmp_band$chromStart[i], ybottom = 0.5 - w, xright = tmp_band$chromEnd[i], ytop = 0.5 + w, col = tmp_band$col[i], border = NA) }
          roundedRect(xleft = max(tmp_band$chromEnd[(index_centromere[1])]), ybottom = 0.5-w, xright = max(tmp_band$chromEnd), ytop = 0.5 + w, lwd=3, rounding = 0.20, xpd=T)
          # draw centromere
          points(x = tmp_band$chromStart[index_centromere[2]], y = 0.5, type = "p", bg=tmp_band$col[which(tmp_band$gieStain == "acen")], pch=21, cex=5, lwd=3)
          # highlight region of interest with some padding
          pd = 10000
          rect(xleft = as.numeric(start) - pd, ybottom = 0.5 - w*3/2, xright = as.numeric(stop) + pd, ytop = 0.5 + w*3/2, border = "red", lwd=4)
          # finally the title
          text(x = max(tmp_band$chromEnd)/2, y = 1, labels = paste(chrom, start, stop, sep = " ~ "), cex=3, xpd=T, adj= 0.5, font=4)
        } else {
          text(x = 0, y = 0, labels = "Not mapping a any autosomal chromosome", pos = 4, xpd=T)    
        } 
        
      # Second plot is the clustering of all samples
        par(mar=c(4, 0, 0, 6))
        # run function to get dendrogram
        dendro = clustering_info[[1]]; ordered_labs = clustering_info[[2]]
        # plot
        dendro %>% set("labels_cex", 1.50) %>% set("labels", ordered_labs) %>% plot(yaxt = 'none', horiz = T, center = T, cex.lab = 1.20, ylim = c(0, nrow(vcf_info_withRef) +2))
        
      # Third plot is the TR sizes
        par(mar=c(4, 0, 0, 0))
        # background plot
        plot(x = 0, 0, pch=16, col="white", xlim=c(min_len, max_len), cex.axis = 1.80, cex.lab = 2, ylim = c(0, nrow(vcf_info_withRef) + 2), xlab = "Size of Tandem Repeat", ylab = "", yaxt="none")
        # get nejm colors
        color = c('navy', 'orange', 'red')
        # grid
        for (x in 1:nrow(vcf_info_withRef)){ abline(h = x, lwd=0.4, col='grey80') }
        # add dashed line for the reference
        segments(x0 = vcf_info_withRef$short_allele[which(vcf_info_withRef$sample == 'GRCh38')], y0 = 1, x1 = vcf_info_withRef$short_allele[which(vcf_info_withRef$sample == 'GRCh38')], y1 = nrow(vcf_info_withRef), col = 'black', lty = 2)
        # main loop across samples
        y_value = 1
        gt_info = data.frame()
        for (s in ordered_labs){
          # get samples' data
          sb = vcf_info_withRef[which(vcf_info_withRef$sample == s),]
          # use 2 colors for heterozygous, 1 for homozygous
          gt_type = ifelse(sb$short_allele == sb$long_allele, 'homo', 'hetero')
          # gt_info = rbind(gt_info, data.frame(sample = s, genotype = gt_type))
          # add points for the haplotype sizes
          if (is.na(gt_type)){ text(x = min_len, y = y_value, labels = 'NA', font = 2) } else if (gt_type == 'homo'){ points(x = sb$short_allele, y = y_value, col = color[1], pch = 16, cex = 1.80) } else { points(x = sb$short_allele, y = y_value, col = color[2], pch = 16, cex = 1.80); points(x = sb$long_allele, y = y_value, col = color[3], pch = 16, cex = 1.80) }
          y_value = y_value + 1
        }
        # Legend 
        legend('top', legend = c('Homozygous', 'Allele 1', 'Allele 2', 'Reference'), x.intersp = 0.2, pch = c(16, 16, 16, NA), col = c(color, 'black'), lty = c(NA, NA, NA, 2), bty = 'n', ncol=2, cex = 1.50)
        
      # Fourth plot is about the motifs -- reference motif with attached the estimated number of copies
        par(mar=c(4, 0, 0, 0))
        # to set the plot limits, need to check the motif
        ref_motif = unique(vcf_info_withRef$short_allele_motif[which(vcf_info_withRef$sample == 'GRCh38')])
        xlim_max = 50
        # background plot
        plot(x = 0, 0, pch=16, col="white", xlim=c(0, xlim_max), ylim = c(0, nrow(vcf_info_withRef) + 2), xlab = "", ylab = "", xaxt='none', yaxt="none", bty = 'n')
        # grid
        colors_rects = rep(c('grey80', 'white'), 200)
        for (x in 1:nrow(vcf_info_withRef)){ rect(xleft = 0, ybottom = x-0.5, xright = xlim_max, ytop = x+0.5, col = colors_rects[x], border = NA) }
        # middle line to divide reference-based and sample-based
        segments(x0 = xlim_max/2, y0 = 0.5, x1 = xlim_max/2, y1 = nrow(vcf_info_withRef)+0.5, lwd = 2, col = 'black')
        # text at the top
        text(x = xlim_max/4, y = nrow(vcf_info_withRef) + 1 + nrow(vcf_info_withRef)*0.02, labels = 'Haplotype 1', font = 2, cex = 1.5, xpd=T)
        text(x = xlim_max/2 + xlim_max/4, y = nrow(vcf_info_withRef) + 1 + nrow(vcf_info_withRef)*0.02, labels = 'Haplotype 2', font = 2, cex = 1.5, xpd=T)
        # the size of the motif rectangle depends on the size of the motif
        max_size_motif = max(nchar(vcf_info_withRef$short_allele_motif), nchar(vcf_info_withRef$long_allele_motif), na.rm=T)
        w = ifelse(max_size_motif <10, 0.8, 0.5); sz_bp = ifelse(max_size_motif <10, 0.8, 0.6); wy = 0.5
        # main loop across samples
        y_value = 1
        for (s in ordered_labs){
          # get samples' data
          sb = vcf_info_withRef[which(vcf_info_withRef$sample == s),]
          if (is.na(sb$short_allele_motif)){
            # just report NA
            text(x = 1, y = y_value, labels = 'NA', adj = c(0, 0.5), font = 2); text(x = xlim_max/2+1, y = y_value, labels = 'NA', adj = c(0, 0.5), font = 2)
          } else {
            # define the colors for the bp
            colors_nt = c("A" = 'blue', "C" = 'red', "G" = 'green', "T" = 'yellow')
            # haplotype 1
              if (nchar(sb$short_allele_motif) >15){
                # if the motif is too long (>15bp) report that motif is too large
                text(x = 1, y = y_value, labels = 'Motif too large (>15bp)', adj = c(0, 0.5)); text(x = xlim_max/2-1, y = y_value, labels = paste0('(', round(as.numeric(sb$short_allele_cn), 1), ')'), adj = c(1, 0.5))
              } else {
                # split the letters
                df_bp = data.frame(letter = strsplit(sb$short_allele_motif, "")[[1]], pos = seq(1, length(strsplit(sb$short_allele_motif, "")[[1]])*(w*2), w*2))
                # assign colors of bp
                for (i in 1:nrow(df_bp)){ df_bp$colors[i] = colors_nt[df_bp$letter[i]] }; df_bp$colors[is.na(df_bp$colors)] = 'black'
                # draw rectangles and add text
                for (i in 1:nrow(df_bp)){ rect(xleft = df_bp$pos[i]-w, xright = df_bp$pos[i]+w, ybottom = y_value-wy, ytop = y_value+wy, col = alpha(df_bp$colors[i], 0.6)); text(x = df_bp$pos[i], y = y_value, labels = df_bp$letter[i], cex = sz_bp, font=2) }
                # finally the copy number
                text(x = xlim_max/2-1, y = y_value, labels = paste0('(', round(as.numeric(sb$short_allele_cn), 1), ')'), font = 2, adj = c(1, 0.5))
              }
            # haplotype 2
              if (nchar(sb$long_allele_motif) >15){
                # if the motif is too long (>15bp) report that motif is too large
                text(x = xlim_max/2+1, y = y_value, labels = 'Motif too large (>15bp)', adj = c(0, 0.5)); text(x = xlim_max-1, y = y_value, labels = paste0('(', round(as.numeric(sb$long_allele_cn), 1), ')'), adj = c(1, 0.5))
              } else {
                # split the letters
                df_bp = data.frame(letter = strsplit(sb$long_allele_motif, "")[[1]], pos = seq(1, length(strsplit(sb$long_allele_motif, "")[[1]])*(w*2), w*2))
                # assign colors of bp
                for (i in 1:nrow(df_bp)){ df_bp$colors[i] = colors_nt[df_bp$letter[i]] }; df_bp$colors[is.na(df_bp$colors)] = 'black'
                # draw rectangles and add text
                for (i in 1:nrow(df_bp)){ rect(xleft = xlim_max/2+df_bp$pos[i]-w, xright = xlim_max/2+df_bp$pos[i]+w, ybottom = y_value-wy, ytop = y_value+wy, col = alpha(df_bp$colors[i], 0.6)); text(x = xlim_max/2+df_bp$pos[i], y = y_value, labels = df_bp$letter[i], font = 2, cex = sz_bp) }
                # finally the copy number
                text(x = xlim_max-1, y = y_value, labels = paste0('(', round(as.numeric(sb$long_allele_cn), 1), ')'), font = 2, adj = c(1, 0.5))
              }
          } 
          y_value = y_value + 1
        }
    }

    # Function to print summary of the run
    summaryRun <- function(rs_vcf, out_dir, out_name, plotFormat, custom_colors, region){
      cat('\n********************')
      cat('\n*** TREAT plot *****')
      cat('\n*** Arguments:')
      cat(paste0('\n*** Input VCF: ', paste(rs_vcf, collapse = ',')))
      cat(paste0('\n*** Region(s): ', paste(region, collapse = ', ')))
      cat(paste0('\n*** Output directory: ', out_dir))
      cat(paste0('\n*** Output name: ', out_name))
      cat(paste0('\n*** Plot format: ', plotFormat))
      cat(paste0('\n*** Custom colors: ', custom_colors))
      return('\n*** Analysis started!')
    }

# Arguments definition
  parser <- ArgumentParser()
  # add arguments: --reads_spannning is the VCF file of the output of read_spanning_analysis
  parser$add_argument("--vcf", default = 'None', help = "VCF file output of TREAT. Multiple files should be comma-separated.")
  # add arguments: --out is the output directory
  parser$add_argument("--out", default = './', help = "Output directory where output will be placed. Default is the current directory.")
  # add arguments: --outname is the name of the output file
  parser$add_argument("--outname", default = 'None', help = "Name of the plot. Default values are msa (for multiple sequence alignment) or repeats (for repeat analysis)")
  # add arguments: --region is the name of the region to plot
  parser$add_argument("--region", default = 'None', help = "Name of the region to plot. Can be one region, multiple (comma separated) regions, or all.")
  # add arguments: --plotformat is whether to give png or pdf as output plot. Default is pdf
  parser$add_argument("--plotformat", default = 'pdf', help = "File format of output plot. Choices are png or pdf. Default value is pdf.")
  # add arguments: --customColors accepts a file with 2 columns: sample name (same as in the data) and an additional column. Samples will be colored.
  parser$add_argument("--customColors", default = 'None', help = "Custom file for coloring options. Accepts a file with 2 columns, i.e sample name and a grouping variable.")
  # add arguments: --path gives the path to the rscript
  parser$add_argument("--path", default = 'None', help = "Path to the Rscript")

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
  # out_name = 'None'
  # plot format
  plotFormat = args$plotformat
  # plotFormat = 'pdf'
  # custom colors
  custom_colors = args$customColors
  # custom_colors = 'None'
  # region
  region = args$region
  # region = 'chr4:39348425-39348483'
  path = args$path

# Check arguments
  # stop if no input data is provided
  run = 'false'
  if ((length(rs_vcf) == 1) && (rs_vcf == 'None')){ stop("Input error: Missing input file(s)!!") } else if (region == 'None'){ stop("No region to plot selected!") } else { run = 'true'}

# Main
  if (run == 'true'){
    # Define the automatic variables: output name
    if (out_name == 'None'){ out_name = 'repeats' }
    # Create directory if not present
    if (!dir.exists(out_dir)){ system(paste0('mkdir ', out_dir)) }
    # Split regions
    region = unlist(strsplit(region, ','))
    # Print summary of the run
    x = summaryRun(rs_vcf, out_dir, out_name, plotFormat, custom_colors, region)

    # Pipeline to plot repeats
    plotRepeatsComplete(rs_vcf, out_dir, out_name, plotFormat, custom_colors, region, path)
  }
