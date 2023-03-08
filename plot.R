library(data.table)
library(stringr)
library(ggplot2)
library(openxlsx)
library(berryFunctions)
library(dendextend)
library(viridis)
library(ggsci)
library(ggpubr)
library(parallel)

# Function to do hierarchical clustering on the TR sizes and eventually motifs
  cluster_TR <- function(sb_rs){
    # create data for clustering
    sb_data = sb_rs[, c('H1_HAPLO_SIZE', 'H2_HAPLO_SIZE')]
    labs = str_split_fixed(sb_rs$SAMPLE, '__', 2)[, 1]
    # fix NAs (homozygous calls)
    sb_data$H2_HAPLO_SIZE <- ifelse(is.na(sb_data$H2_HAPLO_SIZE), sb_data$H1_HAPLO_SIZE, sb_data$H2_HAPLO_SIZE)
    # then do hierarchical clustering
    hr <- as.dendrogram(hclust(dist(sb_data, method = 'euclidean'), method="ward.D2", members = NULL))
    hc = hclust(dist(sb_data, method = 'euclidean'), method="ward.D2", members = NULL)
    ordered_labs = labs[hc$order]
    # return dendrogram + ordered labels
    res = list(hr, ordered_labs)
    return(res)
  }

  # Function to plot TR of the samples -- plotted motif is the consensus motif
  TRplot_consensus_v2 <- function(all_rs_info, all_asm_info, region){
    # 1. Find info about position of interest
    chrom = str_split_fixed(region, ':', 2)[, 1]
    start = str_split_fixed(str_split_fixed(region, ':', 2)[, 2], '-', 2)[, 1]
    stop = str_split_fixed(str_split_fixed(region, ':', 2)[, 2], '-', 2)[, 2]
    gene = unique(all_rs_info$Gene[which(all_rs_info$REGION == region)])
    disease = unique(all_rs_info$Disease[which(all_rs_info$REGION == region)])
    cp_normal = unique(all_rs_info$`Repeat.number.(normal)`[which(all_rs_info$REGION == region)])
    cp_risk = unique(all_rs_info$`Repeat.number.(pre-mutation)`[which(all_rs_info$REGION == region)])
    cp_disease = unique(all_rs_info$`Repeat.number.(disease)`[which(all_rs_info$REGION == region)])
    
    # 2. add information about the HPC or the families and restrict to region and samples of interest
    sb_rs = all_rs_info[which(all_rs_info$REGION == region),]
    sb_asm = all_asm_info[which(all_asm_info$REGION == region),]
    # check if there are results otherwise stop
    if (nrow(sb_rs) == 0 & nrow(sb_asm) == 0){
      print("!!! There are no res -- impossible to draw plot for this SV !!!")
    } else {
      # 3. define maximum and minimum
      all_sizes = c()
      # reads-spanning
      for (i in 1:nrow(sb_rs)){
        tmp_reads_h1 = as.numeric(sb_rs$H1_HAPLO_SIZE[i])
        tmp_reads_h2 = as.numeric(sb_rs$H2_HAPLO_SIZE[i])
        all_sizes = c(all_sizes, tmp_reads_h1, tmp_reads_h2)
      }
      # assembly
      for (i in 1:nrow(sb_asm)){
        tmp_asm_h1 = as.numeric(sb_asm$H1_HAPLO_SIZE[i])
        tmp_asm_h2 = as.numeric(sb_asm$H2_HAPLO_SIZE[i])
        all_sizes = c(all_sizes, tmp_asm_h1, tmp_asm_h2)
      }
      
      # 4. Identify minimum and maximum for the plot window
      min_len <- floor(min(na.omit(all_sizes)) - min(na.omit(all_sizes))*0.15)
      if (min_len == 0){ min_len = -1}
      max_len <- ceiling(max(na.omit(all_sizes)) + max(na.omit(all_sizes))*0.15)
      
      # 5. Parse the chromosomal band
      chr_bands <- fread("/Users/nicco/Library/Mobile Documents/com~apple~CloudDocs/Downloads_shared/hg38_cytogenetic_bands.txt", h=T)
      tmp_band <- chr_bands[which(chr_bands$`#chrom` == chrom),]
      # assign colors depending on transcription levels
      tmp_band$col <- "white"; tmp_band$col[which(tmp_band$gieStain == "gpos25")] <- "grey90";tmp_band$col[which(tmp_band$gieStain == "gpos50")] <- "grey60"; tmp_band$col[which(tmp_band$gieStain == "gpos75")] <- "grey40"; tmp_band$col[which(tmp_band$gieStain == "gpos100")] <- "black"; tmp_band$col[which(tmp_band$gieStain == "acen")] <- "deepskyblue3"
        
      # 6. Layout of the plot
      # chromosome band and location at the top
      # main plot with samples' size of TR
      # on left side, hirarchical clustering of the samples
      # on the right side, motif information
      # bottom is for the repeat and disease information
      layout(matrix(c(1,2,2,2,2,2,2,2,5, 1,3,3,3,3,3,3,3,5, 1,3,3,3,3,3,3,3,5, 1,4,4,4,4,4,4,4,5, 1,4,4,4,4,4,4,4,5), nrow = 9, ncol = 5, byrow = F))
      
      # 7. First plot: chromosomal band
      # define width of chromosome
      w = 0.15
      # graphical parameters
      par(mar=c(0, 6, 2, 0))
      #par(mar=c(1, 7, 3.5, 0))
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
        text(x = max(tmp_band$chromEnd)/2, y = 1, labels = paste(chrom, start, stop, gene, sep = " ~ "), cex=3, xpd=T, adj= 0.5, font=4)
      } else {
        text(x = 0, y = 0, labels = "Not mapping a any autosomal chromosome", pos = 4, xpd=T)    
      } 
      
      # 8. Second plot is the clustering of all samples
      par(mar=c(4, 0, 0, 6))
      # run function to get dendrogram
      sb_reference = sb_rs[which(sb_rs$SAMPLE == 'reference'),]
      sb_rs = sb_rs[which(sb_rs$SAMPLE != 'reference'),]
      res_hclust = cluster_TR(sb_rs)
      dendro = res_hclust[[1]]; ordered_labs = res_hclust[[2]]
      # plot
      dendro %>% set("labels_cex", 1.50) %>% set("labels", ordered_labs) %>% plot(yaxt = 'none', horiz = T, center = T, cex.lab = 1.20)
      
      # 9. Third plot is the TR sizes
      par(mar=c(4, 0, 0, 0))
      # background plot
      plot(x = 0, 0, pch=16, col="white", xlim=c(min_len, max_len), cex.axis = 1.80, cex.lab = 2, ylim = c(0, nrow(sb_rs) + 1), xlab = "Size of Tandem Repeat", ylab = "", yaxt="none")
      # get nejm colors
      color = get_palette(palette = 'nejm', 3)
      # grid
      for (x in 1:nrow(sb_rs)){ abline(h = x, lwd=0.4, col='grey80') }
      # add dashed line for the reference
      segments(x0 = sb_reference$H1_HAPLO_SIZE, y0 = 1, x1 = sb_reference$H1_HAPLO_SIZE, y1 = nrow(sb_rs), col = 'black', lty = 2)
      # add dashed line for the risk-copy number
      risk_cn = as.numeric(str_replace_all(str_replace_all(str_split_fixed(sb_reference$`Repeat.number.(disease)`, '-', 2)[, 1], '>', ''), '<', ''))*nchar(sb_reference$Motif)
      segments(x0 = as.numeric(risk_cn), y0 = 1, x1 = as.numeric(risk_cn), y1 = nrow(sb_rs), col = 'red', lty = 2)
      # main loop across samples -- keep the genotype information
      y_value = 1
      gt_info = data.frame()
      for (s in ordered_labs){
        # get samples' data
        sb = sb_rs[grep(s, sb_rs$SAMPLE),]
        # use 2 colors for heterozygous, 1 for homozygous
        min_hp = min(na.omit(as.numeric(sb$H1_HAPLO_SIZE)), na.omit(as.numeric(sb$H2_HAPLO_SIZE)))
        max_hp = max(na.omit(as.numeric(sb$H1_HAPLO_SIZE)), na.omit(as.numeric(sb$H2_HAPLO_SIZE)))
        gt_type = ifelse(min_hp == max_hp, 'homo', 'hetero')
        gt_info = rbind(gt_info, data.frame(sample = s, genotype = gt_type))
        if (gt_type == 'homo'){
          points(x = min_hp, y = y_value, col = color[1], pch = 16, cex = 1.80)
        } else {
          points(x = min_hp, y = y_value, col = color[2], pch = 16, cex = 1.80); points(x = max_hp, y = y_value, col = color[3], pch = 16, cex = 1.80)
        }
        y_value = y_value + 1
      }
      # Legend 
      legend('top', legend = c('Homozygous', 'Allele 1', 'Allele 2', 'Reference', 'Risk'), x.intersp = 0.2, pch = c(16, 16, 16, NA, NA), col = c(color, 'black', 'red'), lty = c(NA, NA, NA, 2, 2), bty = 'n', ncol = 5)
      
      # 10. Fourth plot is about the motifs -- reference motif with attached the estimated number of copies
      par(mar=c(4, 0, 0, 0))
      # to set the plot limits, need to check the motif
      ref_motif = unique(sb_rs$REFERENCE_MOTIF)
      xlim_max = 100
      # background plot
      plot(x = 0, 0, pch=16, col="white", xlim=c(0, xlim_max), ylim = c(0, nrow(sb_rs) + 1), xlab = "", ylab = "", xaxt='none', yaxt="none", bty = 'n')
      # grid
      colors_rects = rep(c('grey80', 'white'), 200)
      for (x in 1:nrow(sb_rs)){ rect(xleft = 1, ybottom = x-0.5, xright = 100, ytop = x+0.5, col = colors_rects[x], border = NA) }
      # middle line to divide reference-based and sample-based
      segments(x0 = 50, y0 = 0.5, x1 = 50, y1 = nrow(sb_rs)+0.5, lwd = 2, col = 'black')
      # text at the top
      text(x = 25, y = nrow(sb_rs) + 1 + nrow(sb_rs)*0.015, labels = 'Reference-based', font = 2, cex = 1.5, xpd=T)
      text(x = 75, y = nrow(sb_rs) + 1 + nrow(sb_rs)*0.015, labels = 'Sample-based', font = 2, cex = 1.5, xpd=T)
      # main loop across samples
      y_value = 1
      for (s in ordered_labs){
        # get samples' data
        sb = sb_rs[grep(s, sb_rs$SAMPLE),]
        #colors_nt = c("A" = 'blue', "C" = 'red', "G" = 'green', "T" = 'yellow')
        # reference-based copy numbers based on genotype
        h1_cn_ref = round(sb$H1_REFERENCE_MOTIF_CN, 2); h2_cn_ref = round(sb$H2_REFERENCE_MOTIF_CN, 2)
        h1_cn_con = round(sb$H1_CONSENSUS, 2); h2_cn_con = round(sb$H2_CONSENSUS, 2)
        if (is.na(h2_cn_ref)){ 
          #df_ref = data.frame(letter = c('(', strsplit(ref_motif, "")[[1]], ')', h1_cn_ref))
          #df_ref$pos = seq(2, nrow(df_ref)*2, 2)
          #for (i in 1:nrow(df_ref)){ df_ref$colors[i] = colors_nt[df_ref$letter[i]] }; df_ref$colors[is.na(df_ref$colors)] = 'black'
          #for (i in 1:nrow(df_ref)){ text(x = df_ref$pos[i], y = y_value, labels = df_ref$letter[i], col = df_ref$colors[i], font = 2, adj = 0) }
          labe_h1 = paste0('(', ref_motif, ')', h1_cn_ref, collapse = '')
          text(x = 25, y = y_value, labels = labe_h1, adj=0.5, col = color[1])
          labe_h1 = paste0('(', sb$H1_CONSENSUS_MOTIF, ')', h1_cn_con, collapse = '')
          text(x = 75, y = y_value, labels = labe_h1, adj=0.5, col = color[1])
        } else {
          #df_ref = data.frame(letter = c('(', strsplit(ref_motif, "")[[1]], ')', min(h1_cn_ref, h2_cn_ref)))
          #df_ref$letter2 = c('(', strsplit(ref_motif, "")[[1]], ')', max(h1_cn_ref, h2_cn_ref))
          #df_ref$pos = seq(2, nrow(df_ref)*2, 2)
          #df_ref$pos2 = seq(26, 24+nrow(df_ref)*2, 2)
          # add colors
          #for (i in 1:nrow(df_ref)){ df_ref$colors[i] = colors_nt[df_ref$letter[i]] }; df_ref$colors[is.na(df_ref$colors)] = 'black'
          #for (i in 1:nrow(df_ref)){ text(x = df_ref$pos[i], y = y_value, labels = df_ref$letter[i], col = df_ref$colors[i], font = 2, adj = 0) }
          #for (i in 1:nrow(df_ref)){ text(x = df_ref$pos2[i], y = y_value, labels = df_ref$letter2[i], col = df_ref$colors[i], font = 2, adj = 0) }
          labe_h1 = paste0('(', ref_motif, ')', min(h1_cn_ref, h2_cn_ref), collapse = ''); labe_h2 = paste('(', ref_motif, ')', max(h1_cn_ref, h2_cn_ref), collapse = '')
          text(x = 2, y = y_value, labels = labe_h1, adj=0, col = color[2])
          text(x = 48, y = y_value, labels = labe_h2, adj=1, col = color[3])
          # check which motif was larger
          if (h1_cn_con > h2_cn_con){
            labe_h1 = paste0('(', sb$H2_CONSENSUS_MOTIF, ')', h2_cn_con, collapse = ''); labe_h2 = paste0('(', sb$H1_CONSENSUS_MOTIF, ')', h1_cn_con, collapse = '')
          } else {
            labe_h1 = paste0('(', sb$H1_CONSENSUS_MOTIF, ')', h1_cn_con, collapse = ''); labe_h2 = paste0('(', sb$H2_CONSENSUS_MOTIF, ')', h2_cn_con, collapse = '')
          }
          text(x = 52, y = y_value, labels = labe_h1, adj=0, col = color[2])
          text(x = 98, y = y_value, labels = labe_h2, adj=1, col = color[3])
        }
        y_value = y_value + 1
      }
      # explain the *
      text(x = 1, y = 0-nrow(sb_rs)*0.01, xpd = T, labels = '*: different motif across alleles', pos = 4, cex = 1)
      
      # 10. Fifth plot is about information on the repeat and disease
      par(mar=c(2, 12, 2, 18))
      plot(0, 0, pch=16, col="white", xlim=c(0,1), ylim=c(0, 1), bty="n", xaxt="none", yaxt="none", xlab="", ylab="")
      # additional info
      text(x = 0, y = 0.9, labels = paste0("Gene = ", gene), pos = 4, offset = 1, xpd = T, cex = 1.25, font = 2)
      text(x = 0.75, y = 0.9, labels = paste0("Motif = ", sb_reference$Motif), pos = 4, offset = 1, xpd = T, cex = 1.25, font = 2)
      text(x = 0.75, y = 0.6, labels = paste0("Disease = ", disease), pos = 4, offset = 1, xpd = T, cex = 1.25, font = 2)
      text(x = 0, y = 0.6, labels = paste0("# Copies Normal = ", cp_normal), pos = 4, offset = 1, xpd = T, cex = 1.25, font = 2)
      text(x = 0.75, y = 0.3, labels = paste0("# Copies Risk = ", cp_risk), pos = 4, offset = 1, xpd = T, cex = 1.25, font = 2)
      text(x = 0, y = 0.3, labels = paste0("# Copies Disease = ", cp_disease), pos = 4, offset = 1, xpd = T, cex = 1.25, font = 2)
    }
  }

all = rbind(ad, chc)
TRplot(all_rs_info = all, all_asm_info = data.frame(), region = "chr19:524699-525186")
