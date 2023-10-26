make_deseq_obj <- function(count_df, metadata, variable = "Condition", levels = c("control", "mutant")){
  library(tidyverse)
  library(DESeq2)
  metadata[[variable]] <- factor(metadata[[variable]], levels=levels)
  design_formula <- as.formula(paste0("~", variable))
  dds <- DESeqDataSetFromMatrix(countData = count_df,
                                colData = metadata,
                                design = design_formula)
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
  library(ggplot2)
  rld <- rlog(dds_df, blind = FALSE)
  rv <- rowVars(assay(rld))
  select <- order(rv, decreasing = TRUE)[seq_len(min(feature_num, length(rv)))]
  pca <- prcomp(t(assay(rld)[select, ]))

  var <- (pca$sdev^2) / sum(pca$sdev^2)
  
  pca_df <- as.data.frame(pca$x)
  pca_df[variable] <- metadata[variable]
  pca_df$sample_ids <- colnames(dds_df)
  
  unique_levels <- unique(pca_df[, variable])
  color_palette <- rainbow(length(unique_levels))
  level_to_color <- setNames(color_palette, unique_levels)
  pca_df$col <- level_to_color[pca_df[, variable]]
  # pca_df$col <- NA
  # for(i in 1:length(levels(pca_df[,variable]))){
  #   ind1 <- which(pca_df[,variable] == levels(pca_df[,variable])[i])
  #   pca_df$col[ind1] <- i
  # }
  col_var <- variable
  
  # Create a new environment to store col_var
  env <- new.env()
  assign("col_var", col_var, envir = env)
  
  # Create a PCA plot using base R plotting functions
  PCA_plot <- plot(pca_df$PC1, pca_df$PC2, 
       col = pca_df$col,
       xlab = paste0("PC1 (", (round(var[1], digits = 3) * 100), "% variance)"), 
       ylab = paste0("PC2 (", (round(var[2], digits = 3) * 100), "% variance)"),
       main = paste0("PC1 vs PC2 for ", feature_num, " most variable genes"))
  
  # Add labels for sample_ids
  text(pca_df$PC1, pca_df$PC2, labels = pca_df$sample_ids, cex = 0.5, font = 2, pos = 4)
  
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

make_volcano <- function(res_df, alpha = 0.05, fc_cutoff = 1, xlim = 5, ylim = 10){
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
    ylim(0, ylim) + 
    xlim(-xlim, xlim) +
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
  library(dplyr)
  library(tidyverse)
  library(DESeq2)
  library(RColorBrewer)
  library(circlize)
  library(ComplexHeatmap)
  
  FDR_05 <- res_ord_df[res_ord_df$padj<=0.05,]

  rld <- rlog(dds_df, blind = FALSE)
  ind_to_keep <- c(which(metadata_df[[variable]] == control), 
                   which(metadata_df[[variable]] == treatment))

  mat1 <- assay(rld)[rownames(FDR_05), ind_to_keep]

  mat_scaled = t(apply(mat1, 1, scale))

  col = colorRamp2(c(-3, 0, 3), c("#1F68B8", "white", "#FF7F20"))

  meta_sub <- metadata_df[ind_to_keep, ]

  ha1 = HeatmapAnnotation(Group = meta_sub[[variable]],
                          col = list(Group = setNames(c("purple2", "goldenrod2"), c(control, treatment))),
                          show_legend = TRUE)

  ha = columnAnnotation(x = anno_text(meta_sub$Sample_Name,
                                      which="column", rot = 45,
                                      gp = gpar(fontsize = 10)))

  ht1 = Heatmap(mat_scaled, name = "Expression", col = col,
                top_annotation = c(ha1),
                bottom_annotation = c(ha),
                show_row_names = FALSE)

  return(draw(ht1, row_title = "Genes", column_title = "Hierachical clustering of DEGs (padj<=0.05)"))
}


