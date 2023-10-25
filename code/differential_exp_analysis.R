make_deseq_obj <- function(count_df, metadata, variable = "Condition", levels = c("control", "mutant")){
  library(tidyverse)
  library(DESeq2)
  metadata$variable <- factor(metadata$variable, levels=levels)
  dds <- DESeqDataSetFromMatrix(countData = count_df,
                                colData = metadata,
                                design = ~ variable)
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  return(dds)
}

est_size_factors <- function(dds_df){
  dds <- estimateSizeFactors(dds_df)
  counts_norm <- counts(dds, normalized=TRUE)
  return(dds)
}


gene_check_plot <- function(dds_df, variable = "Condition", BF9343 = "gene-BF9343_0242", gene_symbol = "gyrB"){
  library(DESeq2)
  library(ggplot2)
  cnts <- counts(df, normalized=TRUE)
  colnames(cnts) <- colData(df)$variable
  df1 <- data.frame(log2(cnts[BF9343,]), colData(df)$variable)
  colnames(df1) <- c(paste0("log2_gene"), variable)
  p1<- ggplot(df1, aes(variable, log2_gene)) + 
    geom_jitter(aes(color = variable)) + 
    ggtitle(paste0(gene_symbol), " - Log2 Normalized counts")
  return(p1)
}

pca_plot <- function(dds_df, metadata, feature_num = 500, variable = "Condition"){
  library(DESeq2)
  library(tidyverse)
  rld <- rlog(dds_df, blind = FALSE)
  var <- rev(rowVars(assay(rld))[order(rowVars(assay(rld)))])

  rv <- rowVars(assay(rld))
  select <- order(rv, decreasing = TRUE)[seq_len(min(feature_num, length(rv)))]
  pca <- prcomp(t(assay(rld)[select, ]))
  
  pca_df <- as.data.frame(pca$x)
  pca_df$variable <- dds@metadata$variable
  pca_df$sample_ids <- colnames(dds)
  
  pca_df$col <- NA
  for(i in 1:length(levels(pca_df$variable))){
    ind1 <- which(pca_df$variable == levels(pca_df$variable)[i])
    pca_df$col[ind1] <- i
  }
  
  PCA_plot <- plot(pca_df[, 1], pca_df[, 2], 
       xlab = paste0("PC1 (", (round(percentVar[1], digits=3)*100), "% variance)"), 
       ylab = paste0("PC2 (", (round(percentVar[2], digits=3)*100), "% variance)"),
       main=paste0("PC1 vs PC2 for ", var_feature_n, " most variable genes"),
       pch=16, cex=1.35, cex.lab=1.3, cex.axis = 1.15, las=1, 
       panel.first = grid(),
       col=pca_df$col)
  text((pca_df[, 2])~(pca_df[, 1]), labels = pca_df$sample_ids, cex=0.5, font=2, pos=4)
  return(PCA_plot)
}

run_diff_expression <- function(dds_df, variable = "Condition", control = "control", treatment = "mutant"){
  library(tidyverse)
  library(DESeq2)
  dds <- DESeq(dds_df)
  res <- results(dds, alpha = 0.05, 
                 contrast = c(variable, control, treatment), 
                 lfcThreshold = 0)
  res_ord <- res[order(res$padj),]
  res_ord <- res_ord[!is.na(res_ord$padj),]
  return(res_ord)
}

make_volcano <- function(res_df, alpha = 0.05, fc_cutoff = 1){
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  res_tmp <- as.data.frame(res_df)
  res_tmp$cols <- c()
  res_tmp$cols <- NA
  for(i in 1:nrow(res_tmp)){
    if(is.na(res_tmp$pvalue[i])){
      res_tmp$cols[i] <- NA
    }
    else if(res_tmp$pvalue[i]<=alpha & res_tmp$log2FoldChange[i] > fc_cutoff){
      res_tmp$cols[i] <- "indianred"
    } 
    else if(res_tmp$pvalue[i]<=alpha & res_tmp$log2FoldChange[i] < -fc_cutoff){
      res_tmp$cols[i] <- "indianred"
    } 
    else if(res_tmp$pvalue[i]<=alpha & res_tmp$log2FoldChange[i]>-fc_cutoff & res_tmp$log2FoldChange[i]<fc_cutoff){
      res_tmp$cols[i] <- "cornflowerblue"
    } 
    else if(res_tmp$pvalue[i]>alpha & res_tmp$log2FoldChange[i] > fc_cutoff){
      res_tmp$cols[i] <- "gray47" 
    }
    else if(res_tmp$pvalue[i]>alpha & res_tmp$log2FoldChange[i] < -fc_cutoff){
      res_tmp$cols[i] <- "gray47" 
    }
    else if(res_tmp$pvalue[i]>alpha & res_tmp$log2FoldChange[i] < fc_cutoff){
      res_tmp$cols[i] <- "gray10" 
    }
  }
  
  
  res_tmp$gene <- rownames(res_tmp)

  p = ggplot(res_tmp, aes(log2FoldChange, -log10(pvalue))) + 
    geom_point(aes(col=col), alpha = 0.5, size =2, colour = res_tmp$cols, fill = res_tmp$cols)  + 
    xlab("Log2 fold change") + ylab("-log10 P-value") +
    ylim(0, 12) + 
    xlim(-5, 5) +
    geom_hline(yintercept = -log10(alpha), color = "black", linetype = "dashed", size = 0.4) + 
    theme(legend.key = element_blank()) + 
    theme_classic()
  
  p2 <- p + 
    geom_label_repel(data = subset(res_tmp, log2FoldChange > 1 & pvalue < alpha), aes(label = gene), 
                     box.padding   = 0.35,
                     nudge_x = 0.1,
                     nudge_y = 0.1,
                     label.size = 0.1,
                     segment.size = 0.3,
                     segment.color = 'grey50', 
                     size = 3) +
    geom_label_repel(data = subset(res_tmp, log2FoldChange < -1 & pvalue < alpha), aes(label = gene), 
                     box.padding   = 0.35,
                     nudge_x = -0.1,
                     nudge_y = 0.1,
                     label.size = 0.1,
                     segment.size = 0.3,
                     segment.color = 'grey50', 
                     size = 3) +
    geom_vline(xintercept = fc_cutoff, colour = "black", linetype="dotted") + 
    geom_vline(xintercept = -fc_cutoff, colour = "black", linetype="dotted")
  return(p2)  
}

make_clustering_heatmap <- function(dds_df, metadata_df, res_ord_df, variable = "Condition", control = "control", treatment = "mutant"){
  library(tidyverse)
  library(DESeq2)
  library(RColorBrewer)
  library(ComplexHeatmap)
  
  FDR_05 <- res_ord_df[res_ord_df$padj<0.05,]
  
  rld <- rlog(dds_df, blind = FALSE)
  ind_to_keep <- c(which(metadata_df(rld)$variable==control), which(metadata_df(rld)$variable==treatment))
  
  mat1 <- assay(rld)[rownames(FDR_05), ind_to_keep]
  mat_scaled = t(apply(mat1, 1, scale))
  
  col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
  cols1 <- brewer.pal(11, "Paired")
  cols2 <- brewer.pal(9, "Greens")
  
  meta_sub <- metadata_df(dds_df)[ind_to_keep, ]
  
  ha1 = HeatmapAnnotation(Group = meta_sub$variable, 
                          col = list(Group = c(control = cols1[1],
                                               treatment = cols1[2])), 
                          show_legend = TRUE)
  
  ha = columnAnnotation(x = anno_text(meta_sub$variable, 
                                      which="column", rot = 45, 
                                      gp = gpar(fontsize = 10)))
  
  ht1 = Heatmap(mat_scaled, name = "Expression", col = col, 
                top_annotation = c(ha1), 
                bottom_annotation = c(ha),
                show_row_names = FALSE)
  
  return(draw(ht1, row_title = "Genes", column_title = "Hierachical clustering of DEGs (padj<0.05)"))
}



