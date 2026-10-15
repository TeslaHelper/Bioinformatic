# ==============================================================================
# Independent R/DESeq2 verification of GSE245658 hub gene analysis
# Dataset: Wang et al. 2024 Nat Commun — meningeal LEC RNA-seq
# Samples: WT (n=4), 5xFAD (n=3), 5xFAD + 808nm light (n=3)
# Contrasts: (1) light vs no-light within 5xFAD; (2) 5xFAD vs WT (disease context)
# ==============================================================================

suppressMessages({
  library(DESeq2)
  library(pheatmap)
  library(RColorBrewer)
  library(ggplot2)
  library(dplyr)
})
set.seed(42)
outdir <- "/home/claude/paper/revision3/deseq2_out"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Load raw expression file
# ------------------------------------------------------------------------------
cat("=== Loading data ===\n")
raw <- read.table("/tmp/GSE245658_all_genes_expression_annot_WT_txt",
                  header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                  quote = "", fill = TRUE, comment.char = "")

count_cols <- c("WT1_LEC_count","WT2_LEC_count","WT3_LEC_count","WT4_LEC_count",
                "5xFAD1_LEC_count","5xFAD2_LEC_count","5xFAD4_LEC_count",
                "X5xFAD.L1_LEC_count","X5xFAD.L2_LEC_count","X5xFAD.L3_LEC_count")
# Handle column name sanitisation (R replaces special chars with dots)
actual_cols <- colnames(raw)[grepl("count", colnames(raw))]
cat("Count columns found:", paste(actual_cols, collapse=", "), "\n")
count_cols <- actual_cols

# Build count matrix (round fractional RSEM expected counts to integers)
counts_raw <- raw[, count_cols]
counts_mat  <- round(as.matrix(counts_raw))
rownames(counts_mat) <- raw$id
cat("Raw matrix dimensions:", nrow(counts_mat), "genes x", ncol(counts_mat), "samples\n")

# Rename columns for clarity
n_cols <- ncol(counts_mat)
colnames(counts_mat) <- c("WT_1","WT_2","WT_3","WT_4",
                           "FAD_1","FAD_2","FAD_3",
                           "FAD_L_1","FAD_L_2","FAD_L_3")[1:n_cols]

# Add gene symbol for lookup
sym_map <- setNames(raw$Symbol, raw$id)

# ------------------------------------------------------------------------------
# 2. Sample metadata
# ------------------------------------------------------------------------------
coldata <- data.frame(
  sample    = colnames(counts_mat),
  condition = factor(c("WT","WT","WT","WT","FAD","FAD","FAD","FAD_L","FAD_L","FAD_L")[1:n_cols],
                     levels = c("WT","FAD","FAD_L")),
  row.names = colnames(counts_mat)
)
cat("Sample design:\n"); print(coldata)

# ------------------------------------------------------------------------------
# 3. Low-count filter (keep genes with ≥10 summed counts across all samples)
# ------------------------------------------------------------------------------
keep <- rowSums(counts_mat) >= 10
counts_filtered <- counts_mat[keep, ]
cat("Genes after low-count filter (sum>=10):", nrow(counts_filtered), "\n\n")

# ------------------------------------------------------------------------------
# 4. DESeq2 — single model for all three groups, extract two contrasts
# ------------------------------------------------------------------------------
cat("=== Running DESeq2 ===\n")
dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
                               colData   = coldata,
                               design    = ~ condition)
dds <- DESeq(dds, quiet = TRUE)
cat("Model converged.\n")
cat("Size factors:", round(sizeFactors(dds), 3), "\n\n")

# Contrast 1: FAD+Light vs FAD (primary — PBM light effect)
res_light <- results(dds, contrast = c("condition","FAD_L","FAD"), alpha = 0.05)
# Normal shrinkage: need coef name for apeglm; use normal type for contrast-based
res_light_shrunk <- lfcShrink(dds, contrast = c("condition","FAD_L","FAD"),
                               type = "normal", quiet = TRUE)

# Contrast 2: FAD vs WT (disease context)
res_disease <- results(dds, contrast = c("condition","FAD","WT"), alpha = 0.05)
res_disease_shrunk <- lfcShrink(dds, contrast = c("condition","FAD","WT"),
                                 type = "normal", quiet = TRUE)

# How many genome-wide significant?
n_light_sig  <- sum(res_light$padj  < 0.05, na.rm = TRUE)
n_disease_sig <- sum(res_disease$padj < 0.05, na.rm = TRUE)
cat(sprintf("Genome-wide DEGs (padj<0.05): light_vs_FAD=%d  FAD_vs_WT=%d\n\n",
            n_light_sig, n_disease_sig))

# Save full results
full_res <- as.data.frame(res_light_shrunk)
full_res$Symbol <- sym_map[rownames(full_res)]
full_res$light_baseMean <- res_light$baseMean
full_res_light <- full_res[order(full_res$pvalue, na.last = TRUE), ]
write.csv(full_res_light, file.path(outdir, "DESeq2_full_light_vs_FAD.csv"))

# ------------------------------------------------------------------------------
# 5. Hub gene extraction
# ------------------------------------------------------------------------------
hub_mouse <- c("Actb","Ctnnb1","Il6","Tgfb1","Gapdh","Tnf","Trp53",
               "Hif1a","Kdr","Pten","Hsp90aa1","Stat3","Ins1","Ins2","Src","Alb")
hub_human  <- c("ACTB","CTNNB1","IL6","TGFB1","GAPDH","TNF","TP53",
                "HIF1A","KDR","PTEN","HSP90AA1","STAT3","INS","INS","SRC","ALB")
hub_gs     <- c(94.8,89.3,82.6,81.7,77.4,77.1,72.5,70.7,63.7,60.1,57.7,56.2,52.7,52.7,51.4,44.3)

fpkm_cols <- colnames(raw)[grepl("fpkm", colnames(raw))]  # 10 cols
fpkm_wt_cols  <- fpkm_cols[1:4]
fpkm_fad_cols <- fpkm_cols[5:7]
fpkm_l_cols   <- fpkm_cols[8:10]

hub_results <- lapply(seq_along(hub_mouse), function(i) {
  gene <- hub_mouse[i]
  eid  <- names(sym_map)[sym_map == gene]
  eid  <- eid[eid %in% rownames(res_light_shrunk)]

  base <- data.frame(human=hub_human[i], mouse=gene, GoldenScore=hub_gs[i],
                      detected=FALSE, log2FC_light=NA, lfcSE_light=NA, padj_light=NA,
                      log2FC_disease=NA, padj_disease=NA,
                      meanFPKM_WT=NA, meanFPKM_FAD=NA, meanFPKM_FAD_L=NA,
                      stringsAsFactors=FALSE)
  if (length(eid) == 0) return(base)

  eid <- eid[1]
  fpkm_row <- raw[raw$id == eid, fpkm_cols, drop=FALSE]

  base$detected        <- TRUE
  base$log2FC_light    <- round(as.numeric(res_light_shrunk[eid,"log2FoldChange"]),3)
  base$lfcSE_light     <- round(as.numeric(res_light_shrunk[eid,"lfcSE"]),3)
  base$padj_light      <- round(as.numeric(res_light[eid,"padj"]),4)
  base$log2FC_disease  <- round(as.numeric(res_disease_shrunk[eid,"log2FoldChange"]),3)
  base$padj_disease    <- round(as.numeric(res_disease[eid,"padj"]),4)
  base$meanFPKM_WT     <- round(mean(as.numeric(fpkm_row[1,fpkm_wt_cols])),2)
  base$meanFPKM_FAD    <- round(mean(as.numeric(fpkm_row[1,fpkm_fad_cols])),2)
  base$meanFPKM_FAD_L  <- round(mean(as.numeric(fpkm_row[1,fpkm_l_cols])),2)
  base
})

hub_df <- do.call(rbind, hub_results)
write.csv(hub_df, file.path(outdir, "hub_gene_deseq2_verified.csv"), row.names = FALSE)

cat("=== Hub Gene Results (DESeq2 R verification) ===\n")
print(hub_df[, c("human","GoldenScore","detected","log2FC_light","padj_light",
                  "log2FC_disease","meanFPKM_WT","meanFPKM_FAD","meanFPKM_FAD_L")],
      row.names = FALSE)

# IL6 raw counts (key finding)
il6_eid <- names(sym_map)[sym_map == "Il6"]
il6_eid <- il6_eid[il6_eid %in% rownames(counts_filtered)]
if (length(il6_eid)) {
  cat("\n=== IL6 raw counts across samples ===\n")
  print(counts_filtered[il6_eid[1], ])
}

# ------------------------------------------------------------------------------
# 6. PCA Plot — VST transformed
# ------------------------------------------------------------------------------
cat("\n=== Generating PCA plot ===\n")
vsd <- vst(dds, blind = TRUE)

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

condition_colors <- c("WT" = "#3C5488", "FAD" = "#E64B35", "FAD_L" = "#00A087")
condition_labels <- c("WT" = "Wild-type (n=4)",
                       "FAD" = "5xFAD, untreated (n=3)",
                       "FAD_L" = "5xFAD + NIR light (n=3)")

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 5.5, alpha = 0.9) +
  geom_text(nudge_y = 1.2, size = 3.5, fontface = "bold") +
  stat_ellipse(aes(group = condition), level = 0.8, linetype = "dashed", linewidth = 0.7) +
  scale_color_manual(values = condition_colors, labels = condition_labels,
                     name = "Condition") +
  labs(title = "PCA of mLEC transcriptomes (VST-normalised)",
       subtitle = sprintf("GSE245658 (Wang et al., 2024) — %d genes after low-count filter",
                          nrow(counts_filtered)),
       x = sprintf("PC1 (%d%% variance)", pct_var[1]),
       y = sprintf("PC2 (%d%% variance)", pct_var[2])) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom",
        plot.title    = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, color = "grey40"))

ggsave(file.path(outdir, "Figure_PCA_mLEC.png"), p_pca,
       width = 8, height = 6.5, dpi = 300, bg = "white")
ggsave(file.path(outdir, "Figure_PCA_mLEC.pdf"), p_pca,
       width = 8, height = 6.5, bg = "white")
cat("PCA saved.\n")

# ------------------------------------------------------------------------------
# 7. Volcano Plot — light vs FAD, hub genes labelled
# ------------------------------------------------------------------------------
cat("=== Generating Volcano Plot ===\n")
volc_df <- as.data.frame(res_light)
volc_df$Symbol <- sym_map[rownames(volc_df)]
volc_df <- volc_df[!is.na(volc_df$padj) & !is.na(volc_df$log2FoldChange), ]
volc_df$neglog10padj <- -log10(volc_df$padj)
volc_df$sig <- ifelse(volc_df$padj < 0.05 & abs(volc_df$log2FoldChange) > 1,
                       "Significant (padj<0.05, |log2FC|>1)", "Not significant")

# Hub gene overlay (shrunk estimates + raw p-values)
hub_volc <- hub_df[hub_df$detected == TRUE, ]
hub_volc$neglog10padj_raw <- -log10(
  as.numeric(res_light[names(sym_map)[sym_map %in% hub_volc$mouse], "padj"])
)
hub_eid <- sapply(hub_volc$mouse, function(g) {
  e <- names(sym_map)[sym_map == g]; e[e %in% rownames(res_light)][1]
})
hub_volc$log2FC_raw  <- as.numeric(res_light[hub_eid, "log2FoldChange"])
hub_volc$padj_raw    <- as.numeric(res_light[hub_eid, "padj"])
hub_volc$nlp_raw     <- pmin(-log10(hub_volc$padj_raw), 6)

p_volc <- ggplot(volc_df, aes(x = log2FoldChange, y = neglog10padj)) +
  geom_point(aes(color = sig), alpha = 0.3, size = 0.9) +
  scale_color_manual(values = c("Significant (padj<0.05, |log2FC|>1)" = "#E64B35",
                                  "Not significant" = "#AAAAAA"),
                      name = NULL) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#666666", linewidth = 0.6) +
  geom_vline(xintercept = c(-1, 1),   linetype = "dashed", color = "#666666", linewidth = 0.6) +
  # Hub gene points
  geom_point(data = hub_volc, aes(x = log2FC_raw, y = nlp_raw),
             color = "#3C5488", size = 4.5, shape = 18, alpha = 0.9) +
  ggrepel::geom_label_repel(data = hub_volc,
    aes(x = log2FC_raw, y = nlp_raw, label = human),
    size = 3.2, fontface = "bold", color = "#3C5488",
    fill = alpha("white", 0.85), box.padding = 0.5,
    max.overlaps = 20, segment.color = "#3C5488", segment.size = 0.4) +
  labs(title = "Volcano plot: light vs. untreated (5xFAD mLECs)",
       subtitle = "Blue diamonds = PBM–MLV network hub genes\nDashed lines: padj = 0.05 and |log2FC| = 1",
       x = "log2 fold change (5xFAD+NIR light vs. 5xFAD)",
       y = expression(-log[10](p.adjusted))) +
  xlim(-6, 6) + ylim(0, 6.5) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom",
        plot.title    = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10.5, color = "grey40"))

# Need ggrepel
tryCatch({
  library(ggrepel)
  ggsave(file.path(outdir, "Figure_Volcano_mLEC.png"), p_volc,
         width = 9, height = 7, dpi = 300, bg = "white")
  ggsave(file.path(outdir, "Figure_Volcano_mLEC.pdf"), p_volc,
         width = 9, height = 7, bg = "white")
  cat("Volcano plot saved.\n")
}, error = function(e) {
  cat("ggrepel not available, saving without labels:", conditionMessage(e), "\n")
  p_volc2 <- p_volc + geom_text(data = hub_volc,
    aes(x = log2FC_raw, y = nlp_raw, label = human),
    size = 3.2, fontface = "bold", color = "#3C5488", vjust = -1)
  ggsave(file.path(outdir, "Figure_Volcano_mLEC.png"), p_volc2,
         width = 9, height = 7, dpi = 300, bg = "white")
  ggsave(file.path(outdir, "Figure_Volcano_mLEC.pdf"), p_volc2,
         width = 9, height = 7, bg = "white")
  cat("Volcano saved (no repel).\n")
})

# ------------------------------------------------------------------------------
# 8. Hub Gene Heatmap — VST counts, three conditions
# ------------------------------------------------------------------------------
cat("=== Generating Hub Gene Heatmap ===\n")
hub_eids_detected <- sapply(hub_mouse, function(g) {
  e <- names(sym_map)[sym_map == g]; e <- e[e %in% rownames(vsd)]; if(length(e)) e[1] else NA
})
hub_eids_detected <- hub_eids_detected[!is.na(hub_eids_detected)]

vst_hub <- assay(vsd)[hub_eids_detected, ]
rownames(vst_hub) <- sym_map[rownames(vst_hub)]

# Scale across rows (z-score per gene) for interpretable heatmap
vst_hub_z <- t(scale(t(vst_hub)))

anno_col <- data.frame(
  Condition = coldata$condition,
  row.names = coldata$sample
)
anno_colors <- list(
  Condition = c("WT" = "#3C5488", "FAD" = "#E64B35", "FAD_L" = "#00A087")
)

png(file.path(outdir, "Figure_Heatmap_HubGenes.png"), width = 2400, height = 2000, res = 250)
pheatmap(vst_hub_z,
         annotation_col   = anno_col,
         annotation_colors = anno_colors,
         cluster_cols     = FALSE,
         cluster_rows     = TRUE,
         show_colnames    = TRUE,
         color            = colorRampPalette(c("#3C5488","white","#E64B35"))(100),
         breaks           = seq(-2.5, 2.5, length.out = 101),
         main             = "Hub gene expression across mLEC conditions\n(VST-normalized, z-scored per gene; GSE245658)",
         fontsize         = 11,
         border_color     = "grey90")
dev.off()
cat("Heatmap saved.\n")

cat("\n=== ALL DESEQ2 ANALYSES COMPLETE ===\n")
cat("Output directory:", outdir, "\n")
