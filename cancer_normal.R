# Differential Gene Expression Analysis Pipeline
# Author: JahanZaib
# Dataset: TCGA-BRCA (Breast Cancer vs Normal)

# Load required libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(TCGAbiolinks)
library(SummarizedExperiment)

# =============================================================================
# 1. DATA ACQUISITION
# =============================================================================

# Load from downloaded GDC files directory
tcga_data <- GDCprepare(
  query = NULL,  # No query needed if loading from local files
  save = FALSE,
  directory = "path/to/your/GDCdata/",  # Your local GDC download directory
  summarizedExperiment = TRUE
)

# Extract count matrix and sample information
count_matrix <- assay(tcga_data)
sample_info <- colData(tcga_data)

# Clean sample information
sample_info$condition <- ifelse(
  sample_info$sample_type == "Primary Tumor", 
  "Tumor", 
  "Normal"
)

# Filter for protein-coding genes and remove low counts
keep_genes <- rowSums(count_matrix >= 10) >= ncol(count_matrix) * 0.1
count_matrix_filtered <- count_matrix[keep_genes, ]

print(paste("Genes after filtering:", nrow(count_matrix_filtered)))
print(paste("Samples:", ncol(count_matrix_filtered)))

# =============================================================================
# 2. DIFFERENTIAL EXPRESSION ANALYSIS WITH DESeq2
# =============================================================================

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_filtered,
  colData = sample_info,
  design = ~ condition
)

# Set reference level
dds$condition <- relevel(dds$condition, ref = "Normal")

# Run DESeq2 analysis
dds <- DESeq(dds)

# Get results
res <- results(dds, contrast = c("condition", "Tumor", "Normal"))

# Summary of results
summary(res, alpha = 0.05)

# Convert to data frame and add gene symbols
res_df <- as.data.frame(res)
res_df$gene_symbol <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]

# Filter significant genes
sig_genes <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)
print(paste("Significantly differentially expressed genes:", nrow(sig_genes)))

# =============================================================================
# 3. DATA VISUALIZATION
# =============================================================================

# 3.1 Volcano Plot
volcano_plot <- EnhancedVolcano(
  res_df,
  lab = res_df$gene_symbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'Volcano Plot: Tumor vs Normal',
  subtitle = 'Breast Cancer TCGA Data',
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2.0,
  labSize = 4.0,
  maxoverlapsConnectors = 10
)

print(volcano_plot)
ggsave("volcano_plot.png", plot = volcano_plot, width = 12, height = 8, dpi = 300)

# 3.2 MA Plot
plotMA(res, alpha = 0.05, main = "MA Plot: Tumor vs Normal")

# 3.3 Heatmap of top 50 genes
# Get variance stabilized data
vsd <- vst(dds, blind = FALSE)

# Select top 50 most variable genes
top_genes <- order(rowVars(assay(vsd)), decreasing = TRUE)[1:50]
heatmap_data <- assay(vsd)[top_genes, ]

# Create annotation for samples
annotation_col <- data.frame(
  Condition = sample_info$condition
)
rownames(annotation_col) <- colnames(heatmap_data)

# Generate heatmap
pheatmap(
  heatmap_data,
  annotation_col = annotation_col,
  scale = "row",
  show_rownames = TRUE,
  show_colnames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Top 50 Most Variable Genes",
  filename = "heatmap_top50.png",
  width = 12,
  height = 10
)

# 3.4 PCA Plot
plotPCA(vsd, intgroup = "condition") +
  ggtitle("PCA Plot: Tumor vs Normal") +
  theme_minimal()
ggsave("pca_plot.png", width = 10, height = 6, dpi = 300)

# =============================================================================
# 4. FUNCTIONAL ENRICHMENT ANALYSIS
# =============================================================================

# Get significant upregulated and downregulated genes
up_genes <- subset(sig_genes, log2FoldChange > 1)$gene_symbol
down_genes <- subset(sig_genes, log2FoldChange < -1)$gene_symbol

# Convert gene symbols to Entrez IDs
up_genes_clean <- gsub("\\..*", "", up_genes)
head(up_genes_clean, 10)
up_genes_entrez <- clusterProfiler::bitr(up_genes_clean, 
                                         fromType = "ENSEMBL", 
                                         toType = "ENTREZID", 
                                         OrgDb = org.Hs.eg.db)
down_genes_clean <- gsub("\\..*", "", down_genes)
head(down_genes_clean, 10)
down_genes_entrez <- clusterProfiler::bitr(down_genes_clean, 
                                         fromType = "ENSEMBL", 
                                         toType = "ENTREZID", 
                                         OrgDb = org.Hs.eg.db)

# GO enrichment analysis for upregulated genes
go_up <- enrichGO(
  gene = up_genes_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

# GO enrichment analysis for downregulated genes
go_down <- enrichGO(
  gene = down_genes_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

# Visualize GO results
if(nrow(go_up@result) > 0) {
  dotplot(go_up, showCategory = 15) + ggtitle("GO Enrichment: Upregulated Genes")
  ggsave("go_upregulated.png", width = 12, height = 8, dpi = 300)
}

if(nrow(go_down@result) > 0) {
  dotplot(go_down, showCategory = 15) + ggtitle("GO Enrichment: Downregulated Genes")
  ggsave("go_downregulated.png", width = 12, height = 8, dpi = 300)
}

# KEGG pathway analysis
kegg_up <- enrichKEGG(
  gene = up_genes_entrez$ENTREZID,
  organism = 'hsa',
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

kegg_down <- enrichKEGG(
  gene = down_genes_entrez$ENTREZID,
  organism = 'hsa',
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

png("kegg_upregulated_dotplot.png", width = 1200, height = 800, res = 150)
dotplot(kegg_up, showCategory = 15, title = "KEGG Upregulated Pathways")
dev.off()

png("kegg_downregulated_dotplot.png", width = 1200, height = 800, res = 150)
dotplot(kegg_down, showCategory = 15, title = "KEGG Downregulated Pathways")
dev.off()

# Filter for significant results (p.adjust < 0.05)
kegg_up_sig <- kegg_up@result[kegg_up@result$p.adjust < 0.05, ]
kegg_down_sig <- kegg_down@result[kegg_down@result$p.adjust < 0.05, ]

# Save only significant results
write.csv(kegg_up_sig, "kegg_upregulated_significant.csv", row.names = FALSE)
write.csv(kegg_down_sig, "kegg_downregulated_significant.csv", row.names = FALSE)

# Create a comprehensive summary
summary_data <- data.frame(
  Direction = c(rep("Up", nrow(kegg_up@result)), rep("Down", nrow(kegg_down@result))),
  rbind(kegg_up@result, kegg_down@result)
)

write.csv(summary_data, "KEGG_enrichment_summary.csv", row.names = FALSE)

# =============================================================================
# 5. SAVE RESULTS
# =============================================================================

# Save significant genes
write.csv(sig_genes, "significant_genes.csv", row.names = FALSE)

# Save all results
write.csv(res_df, "all_differential_expression_results.csv", row.names = FALSE)

# Save top 100 genes for further analysis
top_100 <- head(sig_genes, 100)
write.csv(top_100, "top_100_significant_genes.csv", row.names = FALSE)

# Save GO enrichment results
if(nrow(go_up@result) > 0) {
  write.csv(go_up@result, "GO_enrichment_upregulated.csv", row.names = FALSE)
}
if(nrow(go_down@result) > 0) {
  write.csv(go_down@result, "GO_enrichment_downregulated.csv", row.names = FALSE)
}

# =============================================================================
# 6. SUMMARY REPORT
# =============================================================================

# Create summary statistics
cat("=== DIFFERENTIAL GENE EXPRESSION ANALYSIS SUMMARY ===\n")
cat("Dataset: TCGA Breast Cancer (BRCA)\n")
cat("Total samples analyzed:", ncol(count_matrix_filtered), "\n")
cat("Tumor samples:", sum(sample_info$condition == "Tumor"), "\n")
cat("Normal samples:", sum(sample_info$condition == "Normal"), "\n")
cat("Genes analyzed:", nrow(count_matrix_filtered), "\n")
cat("Significantly differentially expressed genes (p.adj < 0.05, |log2FC| > 1):", nrow(sig_genes), "\n")
cat("Upregulated genes:", sum(sig_genes$log2FoldChange > 1), "\n")
cat("Downregulated genes:", sum(sig_genes$log2FoldChange < -1), "\n")

# Print top 10 upregulated and downregulated genes
cat("\n=== TOP 10 UPREGULATED GENES ===\n")
top_up <- head(subset(sig_genes, log2FoldChange > 1)[order(-subset(sig_genes, log2FoldChange > 1)$log2FoldChange), ], 10)
print(top_up[, c("gene_symbol", "log2FoldChange", "padj")])

cat("\n=== TOP 10 DOWNREGULATED GENES ===\n")
top_down <- head(subset(sig_genes, log2FoldChange < -1)[order(subset(sig_genes, log2FoldChange < -1)$log2FoldChange), ], 10)
print(top_down[, c("gene_symbol", "log2FoldChange", "padj")])

# ===============Analysis completed===============================
