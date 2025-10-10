# Comprehensive PPI Network Analysis Pipeline for TCGA-BRCA DEGs
# Author: JahanZaib
# Purpose: Protein-Protein Interaction network analysis of DEGs with survival integration
# Date: 2024

# =============================================================================
# LOAD REQUIRED LIBRARIES
# =============================================================================

# Check and install required packages
required_packages <- c(
  "STRINGdb",           # PPI database access
  "igraph",             # Network analysis
  "dplyr",              # Data manipulation
  "ggplot2",            # Visualization
  "ggraph",             # Network visualization
  "tidygraph",          # Tidy graph operations
  "RColorBrewer",       # Color palettes
  "pheatmap",           # Heatmaps
  "visNetwork",         # Interactive networks
  "DOSE",               # Disease enrichment
  "clusterProfiler",    # Pathway analysis
  "org.Hs.eg.db"        # Human annotation
)

# Install missing packages
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("STRINGdb", "DOSE", "clusterProfiler", "org.Hs.eg.db")) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
}

# Load libraries
library(STRINGdb)
library(igraph)
library(dplyr)
library(ggplot2)
library(ggraph)
library(tidygraph)
library(RColorBrewer)
library(pheatmap)
library(visNetwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggrepel)

# =============================================================================
# 1. LOAD DEG AND SURVIVAL ANALYSIS RESULTS
# =============================================================================

print("Loading previous analysis results...")

# Load DEG results
deg_results <- read.csv("all_genes_annotated_results.csv", stringsAsFactors = FALSE)
sig_genes_standard <- read.csv("significant_genes_standard.csv", stringsAsFactors = FALSE)
sig_genes_stringent <- read.csv("significant_genes_stringent.csv", stringsAsFactors = FALSE)

# Load survival results if available
if (file.exists("cox_regression_results_all_genes.csv")) {
  survival_results <- read.csv("cox_regression_results_all_genes.csv", stringsAsFactors = FALSE)
  has_survival <- TRUE
} else {
  has_survival <- FALSE
  print("Note: Survival results not found. Proceeding without survival integration.")
}

print(paste("Total DEGs loaded (standard):", nrow(sig_genes_standard)))
print(paste("Total DEGs loaded (stringent):", nrow(sig_genes_stringent)))

# =============================================================================
# 2. PREPARE INPUT FOR STRING DATABASE
# =============================================================================

print("Preparing gene lists for STRING database...")

# Use stringent DEGs for cleaner network
input_genes <- sig_genes_stringent %>%
  filter(!is.na(gene_symbol) & gene_symbol != "" & !is.na(entrez_id)) %>%
  dplyr::select(gene_symbol, entrez_id, log2FoldChange, padj, direction)

# Separate upregulated and downregulated genes
up_genes <- input_genes %>% filter(log2FoldChange > 0)
down_genes <- input_genes %>% filter(log2FoldChange < 0)

print(paste("Upregulated genes:", nrow(up_genes)))
print(paste("Downregulated genes:", nrow(down_genes)))

# Create directory for outputs
dir.create("PPI_Network_Analysis", showWarnings = FALSE)
setwd("PPI_Network_Analysis")

# =============================================================================
# 3. INITIALIZE STRING DATABASE
# =============================================================================

print("Connecting to STRING database...")
options(timeout = 1000)
# Initialize STRING database (Homo sapiens = 9606)
string_db <- STRINGdb$new(
  version = "12",
  species = 9606,
  score_threshold = 400,  # Medium confidence (0-1000 scale)
  input_directory = ""
)

# Map genes to STRING identifiers
input_mapped <- string_db$map(
  input_genes,
  "gene_symbol",
  removeUnmappedRows = TRUE
)

print(paste("Genes successfully mapped to STRING:", nrow(input_mapped)))
print(paste("Mapping rate:", round(nrow(input_mapped)/nrow(input_genes)*100, 2), "%"))

# Save mapping results
write.csv(input_mapped, "STRING_gene_mapping.csv", row.names = FALSE)

# =============================================================================
# 4. RETRIEVE PPI INTERACTIONS
# =============================================================================

print("Retrieving protein-protein interactions...")
# Get interactions for mapped genes
options(timeout = 3000) # Sets the timeout to 300 seconds
interactions <- string_db$get_interactions(input_mapped$STRING_id)

print(paste("Total interactions retrieved:", nrow(interactions)))

# Add gene information to interactions
interactions_annotated <- interactions %>%
  left_join(input_mapped %>% dplyr::select(STRING_id, gene_symbol, log2FoldChange, direction),
            by = c("from" = "STRING_id")) %>%
  dplyr::rename(from_gene = gene_symbol, from_log2FC = log2FoldChange, from_direction = direction) %>%
  left_join(input_mapped %>% dplyr::select(STRING_id, gene_symbol, log2FoldChange, direction),
            by = c("to" = "STRING_id")) %>%
  dplyr::rename(to_gene = gene_symbol, to_log2FC = log2FoldChange, to_direction = direction)

# Filter for high-confidence interactions
high_conf_interactions <- interactions_annotated %>%
  filter(combined_score >= 600)  # High confidence threshold

print(paste("High-confidence interactions:", nrow(high_conf_interactions)))

# Save interactions
write.csv(interactions_annotated, "all_PPI_interactions.csv", row.names = FALSE)
write.csv(high_conf_interactions, "high_confidence_PPI_interactions.csv", row.names = FALSE)

# =============================================================================
# 5. BUILD IGRAPH NETWORK
# =============================================================================

print("Building igraph network object...")

# Create igraph object from interactions
network <- graph_from_data_frame(
  d = high_conf_interactions[, c("from_gene", "to_gene", "combined_score")],
  vertices = input_mapped[, c("gene_symbol", "log2FoldChange", "padj", "direction")],
  directed = FALSE
)

# Create igraph object with proper edge attributes
edges_with_attr <- high_conf_interactions %>%
  select(from_gene, to_gene, combined_score)

network <- graph_from_data_frame(
  d = edges_with_attr,
  vertices = input_mapped[, c("gene_symbol", "log2FoldChange", "padj", "direction")],
  directed = FALSE
)

# Remove self-loops
network <- igraph::simplify(network, remove.multiple = TRUE, remove.loops = TRUE)

print(paste("Network nodes (genes):", vcount(network)))
print(paste("Network edges (interactions):", ecount(network)))

# Network basic statistics
network_density <- edge_density(network)
network_diameter <- diameter(network)
avg_path_length <- mean_distance(network)

print(paste("Network density:", round(network_density, 4)))
print(paste("Network diameter:", network_diameter))
print(paste("Average path length:", round(avg_path_length, 2)))

# =============================================================================
# 6. CALCULATE NETWORK METRICS (CENTRALITY MEASURES)
# =============================================================================

print("Calculating network centrality metrics...")

# Calculate various centrality measures
V(network)$degree <- degree(network)
V(network)$betweenness <- betweenness(network, normalized = TRUE)
V(network)$closeness <- closeness(network, normalized = TRUE)
V(network)$eigenvector <- eigen_centrality(network)$vector
V(network)$pagerank <- page_rank(network)$vector

# Get clustering coefficient
V(network)$clustering_coef <- transitivity(network, type = "local")
V(network)$clustering_coef[is.nan(V(network)$clustering_coef)] <- 0

# Create centrality dataframe
centrality_df <- data.frame(
  gene = V(network)$name,
  log2FoldChange = V(network)$log2FoldChange,
  padj = V(network)$padj,
  direction = V(network)$direction,
  degree = V(network)$degree,
  betweenness = V(network)$betweenness,
  closeness = V(network)$closeness,
  eigenvector = V(network)$eigenvector,
  pagerank = V(network)$pagerank,
  clustering_coef = V(network)$clustering_coef,
  stringsAsFactors = FALSE
)

# Add survival data if available
if (has_survival) {
  centrality_df <- centrality_df %>%
    left_join(survival_results %>% dplyr::select(gene, HR, cox_pvalue),
              by = "gene")
}

# Rank genes by centrality measures
centrality_df <- centrality_df %>%
  mutate(
    degree_rank = rank(-degree),
    betweenness_rank = rank(-betweenness),
    eigenvector_rank = rank(-eigenvector),
    combined_rank = degree_rank + betweenness_rank + eigenvector_rank
  ) %>%
  arrange(combined_rank)

# Save centrality results
write.csv(centrality_df, "network_centrality_metrics.csv", row.names = FALSE)

# =============================================================================
# 7. IDENTIFY HUB GENES
# =============================================================================

print("Identifying hub genes...")

# Define hub genes using multiple criteria
# Method 1: Top degree (most connections)
degree_threshold <- quantile(centrality_df$degree, 0.90)
hub_genes_degree <- centrality_df %>%
  filter(degree >= degree_threshold) %>%
  arrange(desc(degree))

# Method 2: High betweenness (bridging nodes)
betweenness_threshold <- quantile(centrality_df$betweenness, 0.90)
hub_genes_betweenness <- centrality_df %>%
  filter(betweenness >= betweenness_threshold) %>%
  arrange(desc(betweenness))

# Method 3: High eigenvector centrality (well-connected to important nodes)
eigenvector_threshold <- quantile(centrality_df$eigenvector, 0.90)
hub_genes_eigenvector <- centrality_df %>%
  filter(eigenvector >= eigenvector_threshold) %>%
  arrange(desc(eigenvector))

# Combined hub genes (top 10% in at least 2 metrics)
hub_genes <- centrality_df %>%
  filter(
    (degree >= degree_threshold) +
      (betweenness >= betweenness_threshold) +
      (eigenvector >= eigenvector_threshold) >= 2
  ) %>%
  arrange(desc(degree))

print(paste("Hub genes identified:", nrow(hub_genes)))
print("Top 10 hub genes:")
print(head(hub_genes[, c("gene", "degree", "betweenness", "eigenvector")], 10))

# Save hub gene lists
write.csv(hub_genes, "hub_genes_all.csv", row.names = FALSE)
write.csv(hub_genes_degree, "hub_genes_by_degree.csv", row.names = FALSE)
write.csv(hub_genes_betweenness, "hub_genes_by_betweenness.csv", row.names = FALSE)
write.csv(hub_genes_eigenvector, "hub_genes_by_eigenvector.csv", row.names = FALSE)

# =============================================================================
# 8. COMMUNITY DETECTION (MODULE IDENTIFICATION)
# =============================================================================

print("Detecting network communities/modules...")

# Louvain community detection
set.seed(123)
communities_louvain <- cluster_louvain(network)
V(network)$community_louvain <- membership(communities_louvain)

# Walktrap community detection
communities_walktrap <- cluster_walktrap(network)
V(network)$community_walktrap <- membership(communities_walktrap)

# Fast greedy community detection
communities_fastgreedy <- cluster_fast_greedy(network)
V(network)$community_fastgreedy <- membership(communities_fastgreedy)

print(paste("Louvain communities:", length(unique(V(network)$community_louvain))))
print(paste("Walktrap communities:", length(unique(V(network)$community_walktrap))))
print(paste("Fast greedy communities:", length(unique(V(network)$community_fastgreedy))))

# Add community information to centrality dataframe
centrality_df$community_louvain <- V(network)$community_louvain[match(centrality_df$gene, V(network)$name)]
centrality_df$community_walktrap <- V(network)$community_walktrap[match(centrality_df$gene, V(network)$name)]

# Analyze each community
community_analysis <- centrality_df %>%
  group_by(community_louvain) %>%
  summarise(
    n_genes = n(),
    avg_degree = mean(degree),
    avg_log2FC = mean(abs(log2FoldChange)),
    n_upregulated = sum(direction == "Upregulated"),
    n_downregulated = sum(direction == "Downregulated"),
    top_genes = paste(head(gene[order(-degree)], 5), collapse = ", ")
  ) %>%
  arrange(desc(n_genes))

write.csv(community_analysis, "community_analysis.csv", row.names = FALSE)

# =============================================================================
# 9. PATHWAY ENRICHMENT FOR HUB GENES AND COMMUNITIES
# =============================================================================

print("Performing pathway enrichment for hub genes...")

# Get Entrez IDs for hub genes
hub_genes_entrez <- input_genes %>%
  filter(gene_symbol %in% hub_genes$gene) %>%
  pull(entrez_id) %>%
  unique() %>%
  na.omit()

# KEGG pathway enrichment for hub genes
if (length(hub_genes_entrez) > 5) {
  hub_kegg <- enrichKEGG(
    gene = hub_genes_entrez,
    organism = 'hsa',
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  
  if (nrow(hub_kegg@result) > 0) {
    write.csv(hub_kegg@result, "hub_genes_KEGG_enrichment.csv", row.names = FALSE)
    
    # Visualize
    hub_kegg_plot <- dotplot(hub_kegg, showCategory = 15) +
      ggtitle("KEGG Pathways: Hub Genes")
    ggsave("hub_genes_KEGG_pathways.png", hub_kegg_plot, width = 12, height = 8, dpi = 300)
  }
}

# GO enrichment for hub genes
hub_go_bp <- enrichGO(
  gene = hub_genes_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

if (nrow(hub_go_bp@result) > 0) {
  write.csv(hub_go_bp@result, "hub_genes_GO_BP_enrichment.csv", row.names = FALSE)
  
  hub_go_plot <- dotplot(hub_go_bp, showCategory = 15) +
    ggtitle("GO Biological Processes: Hub Genes")
  ggsave("hub_genes_GO_BP.png", hub_go_plot, width = 12, height = 8, dpi = 300)
}

# Enrichment for largest communities
for (i in 1:min(3, nrow(community_analysis))) {
  comm_id <- community_analysis$community_louvain[i]
  comm_genes <- centrality_df %>%
    filter(community_louvain == comm_id) %>%
    pull(gene)
  
  comm_entrez <- input_genes %>%
    filter(gene_symbol %in% comm_genes) %>%
    pull(entrez_id) %>%
    unique() %>%
    na.omit()
  
  if (length(comm_entrez) > 5) {
    comm_kegg <- enrichKEGG(
      gene = comm_entrez,
      organism = 'hsa',
      pvalueCutoff = 0.05
    )
    
    if (nrow(comm_kegg@result) > 0) {
      write.csv(comm_kegg@result, 
                paste0("community_", comm_id, "_KEGG_enrichment.csv"), 
                row.names = FALSE)
    }
  }
}

# =============================================================================
# 10. NETWORK VISUALIZATION
# =============================================================================

print("Creating network visualizations...")

# 10.1 Full network with hub genes highlighted
set.seed(123)

# Prepare layout
layout_fr <- layout_with_fr(network)

# Color by direction
V(network)$color <- ifelse(V(network)$direction == "Upregulated", "#E7B800", "#2E9FDF")

# Size by degree
V(network)$size <- scales::rescale(V(network)$degree, to = c(3, 15))

# Highlight hub genes
V(network)$frame.color <- ifelse(V(network)$name %in% hub_genes$gene, "red", "gray")
V(network)$frame.width <- ifelse(V(network)$name %in% hub_genes$gene, 2, 0.5)

# Plot full network
png("PPI_network_full.png", width = 14, height = 14, units = "in", res = 300)
plot(network,
     layout = layout_fr,
     vertex.label = ifelse(V(network)$name %in% head(hub_genes$gene, 20), 
                           V(network)$name, NA),
     vertex.label.cex = 0.7,
     vertex.label.color = "black",
     vertex.label.dist = 1,
     edge.color = "gray80",
     edge.width = 0.5,
     main = "PPI Network: Hub Genes Highlighted (Red Border)")
legend("topright", 
       legend = c("Upregulated", "Downregulated", "Hub Gene"),
       col = c("#E7B800", "#2E9FDF", "white"),
       pch = c(21, 21, 21),
       pt.bg = c("#E7B800", "#2E9FDF", "white"),
       pt.cex = 2,
       bty = "n")
dev.off()

# 10.2 Hub gene subnetwork
hub_gene_names <- hub_genes$gene
hub_subnetwork <- induced_subgraph(network, V(network)[V(network)$name %in% hub_gene_names])

png("PPI_network_hub_genes.png", width = 12, height = 12, units = "in", res = 300)
plot(hub_subnetwork,
     layout = layout_with_fr(hub_subnetwork),
     vertex.label = V(hub_subnetwork)$name,
     vertex.label.cex = 0.8,
     vertex.label.color = "black",
     vertex.size = scales::rescale(V(hub_subnetwork)$degree, to = c(8, 20)),
     edge.color = "gray60",
     edge.width = 1,
     main = "Hub Gene Subnetwork")
dev.off()

# 10.3 Community-colored network
V(network)$community_color <- rainbow(length(unique(V(network)$community_louvain)))[V(network)$community_louvain]

png("PPI_network_communities_simple.png", width = 14, height = 14, units = "in", res = 300)
plot(network,
     layout = layout_fr,
     vertex.color = V(network)$community_color,
     vertex.label = NA,
     vertex.size = 3,
     edge.color = "gray80",
     edge.width = 0.5,
     main = "PPI Network: Community Structure")
dev.off()

# 10.4 Interactive network with visNetwork
vis_nodes <- data.frame(
  id = V(network)$name,
  label = V(network)$name,
  title = paste0("<b>", V(network)$name, "</b><br>",
                 "Degree: ", V(network)$degree, "<br>",
                 "log2FC: ", round(V(network)$log2FoldChange, 2), "<br>",
                 "Direction: ", V(network)$direction),
  group = V(network)$direction,
  value = V(network)$degree,
  color = V(network)$color
)

vis_edges <- data.frame(
  from = as_edgelist(network)[,1],
  to = as_edgelist(network)[,2],
  value = E(network)$combined_score / 1000
)

visnet <- visNetwork(vis_nodes, vis_edges, width = "100%", height = "800px") %>%
  visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
             nodesIdSelection = TRUE) %>%
  visPhysics(stabilization = TRUE) %>%
  visInteraction(navigationButtons = TRUE) %>%
  visLegend()

visSave(visnet, file = "PPI_network_interactive.html")

# =============================================================================
# 11. ADVANCED VISUALIZATIONS
# =============================================================================

# 11.1 Centrality heatmap
top_genes_for_heatmap <- head(centrality_df, 50)

centrality_matrix <- as.matrix(top_genes_for_heatmap[, c("degree", "betweenness", 
                                                         "closeness", "eigenvector")])
rownames(centrality_matrix) <- top_genes_for_heatmap$gene

# Scale for visualization
centrality_matrix_scaled <- scale(centrality_matrix)

# Annotation
annotation_row <- data.frame(
  Direction = top_genes_for_heatmap$direction,
  row.names = top_genes_for_heatmap$gene
)

ann_colors <- list(
  Direction = c("Upregulated" = "#E7B800", "Downregulated" = "#2E9FDF")
)

pheatmap(centrality_matrix_scaled,
         annotation_row = annotation_row,
         annotation_colors = ann_colors,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         fontsize_row = 7,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main = "Network Centrality Heatmap: Top 50 Genes",
         filename = "centrality_heatmap.png",
         width = 8,
         height = 12)

# 11.2 Degree distribution plot
degree_dist_plot <- ggplot(centrality_df, aes(x = degree)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = degree_threshold, linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = degree_threshold + 2, y = Inf, 
           label = "Hub threshold (90th percentile)", 
           vjust = 2, color = "red") +
  scale_x_continuous(breaks = seq(0, max(centrality_df$degree), by = 5)) +
  labs(title = "Degree Distribution of PPI Network",
       x = "Degree (Number of Connections)",
       y = "Frequency") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("degree_distribution.png", degree_dist_plot, width = 10, height = 6, dpi = 300)

# 11.3 Scatter plot: Degree vs Betweenness
scatter_plot <- ggplot(centrality_df, aes(x = degree, y = betweenness)) +
  geom_point(aes(color = direction, size = abs(log2FoldChange)), alpha = 0.6) +
  geom_text_repel(data = head(hub_genes, 15),
                  aes(label = gene),
                  size = 3,
                  max.overlaps = 20) +
  scale_color_manual(values = c("Upregulated" = "#E7B800", "Downregulated" = "#2E9FDF")) +
  labs(title = "Network Centrality: Degree vs Betweenness",
       x = "Degree Centrality",
       y = "Betweenness Centrality",
       color = "Expression",
       size = "|log2FC|") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("degree_vs_betweenness.png", scatter_plot, width = 12, height = 8, dpi = 300)

# 11.4 Community size barplot
# Top N communities by size
top_n <- 10
comm_size_plot <- community_analysis %>%
  arrange(desc(n_genes)) %>%
  slice(1:top_n) %>%
  ggplot(aes(x = reorder(as.factor(community_louvain), -n_genes), y = n_genes)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = n_genes), vjust = -0.5, size = 3) +
  labs(title = paste("Top", top_n, "Community Sizes in PPI Network"),
       x = "Community ID",
       y = "Number of Genes") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("community_sizes_top10.png", comm_size_plot, width = 10, height = 6, dpi = 300)

# You can use all communities for interactivity
plot_data <- community_analysis %>%
  arrange(desc(n_genes))

comm_size_plotly <- plot_ly(
  data = plot_data,
  x = ~as.factor(community_louvain),
  y = ~n_genes,
  type = "bar",
  text = ~paste("Community:", community_louvain, "<br>Genes:", n_genes),
  hoverinfo = "text"
) %>%
  layout(
    title = "Community Sizes in PPI Network",
    xaxis = list(title = "Community ID", tickangle = -45),
    yaxis = list(title = "Number of Genes")
  )

# Save as interactive HTML
htmlwidgets::saveWidget(comm_size_plotly, "community_sizes_interactive.html")

# =============================================================================
# 12. INTEGRATE WITH SURVIVAL DATA (IF AVAILABLE)
# =============================================================================

if (has_survival) {
  print("Analyzing survival associations of hub genes...")
  
  # Filter hub genes with survival data
  hub_survival <- centrality_df %>%
    filter(!is.na(HR)) %>%
    arrange(cox_pvalue)
  
  # Identify prognostic hub genes
  prognostic_hubs <- hub_survival %>%
    filter(cox_pvalue < 0.05) %>%
    arrange(desc(degree))
  
  print(paste("Prognostic hub genes:", nrow(prognostic_hubs)))
  
  if (nrow(prognostic_hubs) > 0) {
    write.csv(prognostic_hubs, "prognostic_hub_genes.csv", row.names = FALSE)
    
    # Visualize hub genes with survival
    hub_survival_plot <- ggplot(hub_survival, 
                                aes(x = degree, y = -log10(cox_pvalue))) +
      geom_point(aes(color = ifelse(HR > 1, "Poor Prognosis", "Good Prognosis"),
                     size = abs(log2FoldChange)), alpha = 0.6) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
      geom_text_repel(data = head(prognostic_hubs, 15),
                      aes(label = gene),
                      size = 3,
                      max.overlaps = 20) +
      scale_color_manual(values = c("Poor Prognosis" = "#E74C3C", "Good Prognosis" = "#27AE60")) +
      labs(title = "Hub Genes: Network Centrality vs Survival Association",
           x = "Degree Centrality",
           y = "-log10(Cox P-value)",
           color = "Survival Association",
           size = "|log2FC|") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    ggsave("hub_genes_survival_association.png", hub_survival_plot, 
           width = 12, height = 8, dpi = 300)
  }
}

# =============================================================================
# 13. COMPREHENSIVE SUMMARY REPORT
# =============================================================================

print("Generating comprehensive summary report...")

# Create summary
summary_report <- list(
  network_stats = data.frame(
    Metric = c("Total nodes (genes)", 
               "Total edges (interactions)",
               "Network density",
               "Average degree",
               "Network diameter",
               "Average path length",
               "Clustering coefficient",
               "Number of communities"),
    Value = c(vcount(network),
              ecount(network),
              round(network_density, 4),
              round(mean(V(network)$degree), 2),
              network_diameter,
              round(avg_path_length, 2),
              round(transitivity(network, type = "global"), 3),
              length(unique(V(network)$community_louvain)))
  ),
  
  hub_gene_summary = data.frame(
    Category = c("Total hub genes",
                 "Upregulated hubs",
                 "Downregulated hubs",
                 "Average hub degree",
                 "Max hub degree"),
    Value = c(nrow(hub_genes),
              sum(hub_genes$direction == "Upregulated"),
              sum(hub_genes$direction == "Downregulated"),
              round(mean(hub_genes$degree), 2),
              max(hub_genes$degree))
  ),
  
  top_hub_genes = head(hub_genes[, c("gene", "degree", "betweenness", 
                                     "eigenvector", "direction")], 20)
)

# Add survival info if available
if (has_survival && nrow(prognostic_hubs) > 0) {
  summary_report$prognostic_hubs <- data.frame(
    Metric = c("Prognostic hub genes (p<0.05)",
               "Poor prognosis hubs (HR>1)",
               "Good prognosis hubs (HR<1)"),
    Value = c(nrow(prognostic_hubs),
              sum(prognostic_hubs$HR > 1),
              sum(prognostic_hubs$HR < 1))
  )
}

# Print summary
cat("\n===========================================\n")
cat("PPI NETWORK ANALYSIS SUMMARY\n")
cat("===========================================\n\n")

cat("NETWORK STATISTICS:\n")
print(summary_report$network_stats)

cat("\n\nHUB GENE SUMMARY:\n")
print(summary_report$hub_gene_summary)

cat("\n\nTOP 20 HUB GENES:\n")
print(summary_report$top_hub_genes)

if (has_survival && exists("prognostic_hubs")) {
  cat("\n\nPROGNOSTIC HUB GENES:\n")
  print(summary_report$prognostic_hubs)
  cat("\nTop prognostic hub genes:\n")
  print(head(prognostic_hubs[, c("gene", "degree", "HR", "cox_pvalue")], 10))
}

# Save summary
saveRDS(summary_report, "PPI_network_summary.rds")

# Create markdown report
report_text <- paste0(
  "# PPI Network Analysis Report\n\n",
  "## Dataset: TCGA-BRCA DEGs\n",
  "Date: ", Sys.Date(), "\n\n",
  "### Network Statistics\n",
  "- Total genes (nodes): ", vcount(network), "\n",
  "- Total interactions (edges): ", ecount(network), "\n",
  "- Network density: ", round(network_density, 4), "\n",
  "- Average degree: ", round(mean(V(network)$degree), 2), "\n",
  "- Network diameter: ", network_diameter, "\n",
  "- Average clustering coefficient: ", round(transitivity(network, type = "global"), 3), "\n\n",
  "### Hub Genes\n",
  "- Total hub genes identified: ", nrow(hub_genes), "\n",
  "- Upregulated hubs: ", sum(hub_genes$direction == "Upregulated"), "\n",
  "- Downregulated hubs: ", sum(hub_genes$direction == "Downregulated"), "\n",
  "- Average hub degree: ", round(mean(hub_genes$degree), 2), "\n\n",
  "### Top 10 Hub Genes\n"
)

# Add top hub genes
for (i in 1:min(10, nrow(hub_genes))) {
  gene_info <- hub_genes[i, ]
  report_text <- paste0(report_text,
                        i, ". **", gene_info$gene, "**\n",
                        "   - Degree: ", gene_info$degree, "\n",
                        "   - Betweenness: ", round(gene_info$betweenness, 4), "\n",
                        "   - Expression: ", gene_info$direction, " (log2FC: ", round(gene_info$log2FoldChange, 2), ")\n")
  
  if (has_survival && !is.na(gene_info$HR)) {
    report_text <- paste0(report_text,
                          "   - Survival: HR = ", round(gene_info$HR, 2), 
                          " (p = ", format(gene_info$cox_pvalue, scientific = TRUE, digits = 2), ")\n")
  }
  report_text <- paste0(report_text, "\n")
}

report_text <- paste0(report_text,
                      "\n### Community Structure\n",
                      "- Number of communities detected: ", length(unique(V(network)$community_louvain)), "\n",
                      "- Largest community size: ", max(community_analysis$n_genes), " genes\n\n",
                      "### Key Pathways in Hub Genes\n"
)

# Add pathway info if available
if (exists("hub_kegg") && nrow(hub_kegg@result) > 0) {
  report_text <- paste0(report_text, "Top enriched KEGG pathways:\n")
  for (i in 1:min(5, nrow(hub_kegg@result))) {
    report_text <- paste0(report_text,
                          i, ". ", hub_kegg@result$Description[i], 
                          " (p = ", format(hub_kegg@result$p.adjust[i], scientific = TRUE, digits = 2), ")\n")
  }
}

if (has_survival && exists("prognostic_hubs")) {
  report_text <- paste0(report_text,
                        "\n### Prognostic Hub Genes\n",
                        "- Prognostic hub genes (p < 0.05): ", nrow(prognostic_hubs), "\n",
                        "- Poor prognosis hubs (HR > 1): ", sum(prognostic_hubs$HR > 1), "\n",
                        "- Good prognosis hubs (HR < 1): ", sum(prognostic_hubs$HR < 1), "\n")
}

report_text <- paste0(report_text,
                      "\n### Output Files Generated\n",
                      "1. **Network Files**\n",
                      "   - all_PPI_interactions.csv - All protein interactions\n",
                      "   - high_confidence_PPI_interactions.csv - High-confidence interactions\n",
                      "   - network_centrality_metrics.csv - Centrality measures for all genes\n",
                      "   - PPI_network_interactive.html - Interactive network visualization\n\n",
                      "2. **Hub Gene Files**\n",
                      "   - hub_genes_all.csv - All identified hub genes\n",
                      "   - hub_genes_by_degree.csv - Hubs ranked by degree\n",
                      "   - hub_genes_by_betweenness.csv - Hubs ranked by betweenness\n",
                      "   - hub_genes_KEGG_enrichment.csv - KEGG pathways for hub genes\n",
                      "   - hub_genes_GO_BP_enrichment.csv - GO biological processes\n\n",
                      "3. **Community Files**\n",
                      "   - community_analysis.csv - Analysis of network communities\n",
                      "   - community_*_KEGG_enrichment.csv - Pathways for each community\n\n",
                      "4. **Visualization Files**\n",
                      "   - PPI_network_full.png - Full network with hub genes highlighted\n",
                      "   - PPI_network_hub_genes.png - Hub gene subnetwork\n",
                      "   - PPI_network_communities.png - Network colored by communities\n",
                      "   - centrality_heatmap.png - Heatmap of centrality metrics\n",
                      "   - degree_distribution.png - Degree distribution histogram\n",
                      "   - degree_vs_betweenness.png - Centrality correlation plot\n",
                      "   - community_sizes.png - Community size distribution\n"
)

if (has_survival) {
  report_text <- paste0(report_text,
                        "   - hub_genes_survival_association.png - Hub genes vs survival\n")
}

writeLines(report_text, "PPI_Network_Analysis_Report.md")

# =============================================================================
# 14. IDENTIFY KEY PATHWAYS FROM HUB GENES
# =============================================================================

print("Analyzing key cancer pathways in hub genes...")

# Define key cancer-related pathways
key_pathways <- c(
  "PI3K-Akt signaling pathway",
  "Cell cycle",
  "p53 signaling pathway",
  "MAPK signaling pathway",
  "Apoptosis",
  "Focal adhesion",
  "ECM-receptor interaction",
  "Pathways in cancer",
  "Breast cancer",
  "mTOR signaling pathway"
)

# Check if hub genes are enriched in these pathways
if (exists("hub_kegg") && nrow(hub_kegg@result) > 0) {
  key_pathway_results <- hub_kegg@result %>%
    filter(Description %in% key_pathways) %>%
    arrange(p.adjust)
  
  if (nrow(key_pathway_results) > 0) {
    write.csv(key_pathway_results, "key_cancer_pathways_in_hubs.csv", row.names = FALSE)
    
    print("\nKey cancer pathways enriched in hub genes:")
    print(key_pathway_results[, c("Description", "Count", "p.adjust")])
    
    # Extract genes from key pathways
    for (i in 1:nrow(key_pathway_results)) {
      pathway_name <- key_pathway_results$Description[i]
      pathway_genes <- unlist(strsplit(key_pathway_results$geneID[i], "/"))
      
      # Get hub genes in this pathway
      pathway_hubs <- hub_genes %>%
        filter(gene %in% pathway_genes)
      
      if (nrow(pathway_hubs) > 0) {
        pathway_filename <- gsub(" ", "_", gsub("-", "_", pathway_name))
        write.csv(pathway_hubs, 
                  paste0("hub_genes_in_", pathway_filename, ".csv"),
                  row.names = FALSE)
        
        print(paste("\nHub genes in", pathway_name, ":", nrow(pathway_hubs)))
        print(pathway_hubs$gene)
      }
    }
  }
}

# =============================================================================
# 15. NETWORK COMPARISON: UP vs DOWN-REGULATED GENES
# =============================================================================

print("Comparing networks of up and down-regulated genes...")

# Create separate networks
up_genes_for_network <- input_mapped %>%
  filter(direction == "Upregulated") %>%
  pull(STRING_id)

down_genes_for_network <- input_mapped %>%
  filter(direction == "Downregulated") %>%
  pull(STRING_id)

# Get interactions for each
up_interactions <- string_db$get_interactions(up_genes_for_network)
down_interactions <- string_db$get_interactions(down_genes_for_network)

# Create subnetworks
up_network <- induced_subgraph(network, V(network)[V(network)$direction == "Upregulated"])
down_network <- induced_subgraph(network, V(network)[V(network)$direction == "Downregulated"])

# Compare network properties
network_comparison <- data.frame(
  Metric = c("Number of nodes", "Number of edges", "Density", 
             "Average degree", "Clustering coefficient", "Average path length"),
  Upregulated = c(
    vcount(up_network),
    ecount(up_network),
    round(edge_density(up_network), 4),
    round(mean(degree(up_network)), 2),
    round(transitivity(up_network, type = "global"), 3),
    round(mean_distance(up_network), 2)
  ),
  Downregulated = c(
    vcount(down_network),
    ecount(down_network),
    round(edge_density(down_network), 4),
    round(mean(degree(down_network)), 2),
    round(transitivity(down_network, type = "global"), 3),
    round(mean_distance(down_network), 2)
  )
)

write.csv(network_comparison, "network_comparison_up_vs_down.csv", row.names = FALSE)
print("\nNetwork comparison (Up vs Down-regulated):")
print(network_comparison)

# =============================================================================
# 16. EXPORT NETWORK FOR CYTOSCAPE
# =============================================================================

print("Exporting network files for Cytoscape...")

# Export node attributes
# Calculate network metrics
V(network)$degree <- degree(network)
V(network)$betweenness <- betweenness(network)
V(network)$closeness <- closeness(network)
V(network)$eigenvector <- eigen_centrality(network)$vector
V(network)$pagerank <- page_rank(network)$vector
V(network)$clustering_coef <- transitivity(network, type = "local")

# Community detection
community_result <- cluster_louvain(network)
V(network)$community_louvain <- membership(community_result)
# Prepare node attributes dataframe
node_attributes <- data.frame(
  node_id = V(network)$name,
  gene_symbol = V(network)$name,
  log2FoldChange = V(network)$log2FoldChange,
  padj = V(network)$padj,
  direction = V(network)$direction,
  degree = V(network)$degree,
  betweenness = V(network)$betweenness,
  closeness = V(network)$closeness,
  eigenvector = V(network)$eigenvector,
  pagerank = V(network)$pagerank,
  clustering_coef = V(network)$clustering_coef,
  community = V(network)$community_louvain,
  is_hub = V(network)$name %in% hub_genes$gene,
  stringsAsFactors = FALSE
)

# Add survival data if available
if (has_survival) {
  node_attributes <- node_attributes %>%
    left_join(survival_results %>% select(gene, HR, cox_pvalue),
              by = c("gene_symbol" = "gene"))
}

write.csv(node_attributes, "cytoscape_node_attributes.csv", row.names = FALSE)

# Export edge list with attributes
edge_list <- igraph::as_data_frame(network, what = "edges")
edge_list$weight <- E(network)$combined_score

write.csv(edge_list, "cytoscape_edge_list.csv", row.names = FALSE)

# Export in SIF format (simple interaction format)
sif_data <- data.frame(
  source = as_edgelist(network)[,1],
  interaction = "pp",  # protein-protein interaction
  target = as_edgelist(network)[,2]
)

write.table(sif_data, "network.sif", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# =============================================================================
# 17. ADDITIONAL NETWORK METRICS AND ANALYSES
# =============================================================================

print("Calculating additional network metrics...")

# K-core decomposition
V(network)$coreness <- coreness(network)

# Network motifs (triangles)
triangles_count <- sum(count_triangles(network)) / 3
print(paste("Number of triangles in network:", triangles_count))

# Assortativity (degree correlation)
assortativity_degree <- assortativity_degree(network)
print(paste("Degree assortativity coefficient:", round(assortativity_degree, 3)))

# Network resilience (largest connected component size)
components <- igraph::components(network)
largest_component_size <- max(components$csize)
print(paste("Largest connected component size:", largest_component_size, 
            "(", round(largest_component_size/vcount(network)*100, 1), "%)"))

# Rich-club coefficient (tendency of high-degree nodes to connect)
max_degree <- max(degree(network))
rich_club_coef <- sapply(1:max_degree, function(k) {
  subg <- induced_subgraph(network, V(network)[degree(network) > k])
  if (vcount(subg) > 1) {
    edge_density(subg)
  } else {
    NA
  }
})

rich_club_df <- data.frame(
  degree_threshold = 1:max_degree,
  rich_club_coefficient = rich_club_coef
)

write.csv(rich_club_df, "rich_club_coefficients.csv", row.names = FALSE)

# Plot rich-club coefficient
rich_club_plot <- ggplot(rich_club_df, aes(x = degree_threshold, y = rich_club_coefficient)) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(color = "darkblue", size = 2) +
  labs(title = "Rich-Club Coefficient",
       x = "Degree Threshold",
       y = "Rich-Club Coefficient") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("rich_club_coefficient.png", rich_club_plot, width = 10, height = 6, dpi = 300)

# Update centrality dataframe with k-core
centrality_df$coreness <- V(network)$coreness[match(centrality_df$gene, V(network)$name)]
write.csv(centrality_df, "network_centrality_metrics_complete.csv", row.names = FALSE)

# =============================================================================
# 18. BOTTLENECK GENES IDENTIFICATION
# =============================================================================

print("Identifying bottleneck genes...")

# Bottleneck genes: high betweenness but not necessarily high degree
# These genes are critical for network communication

bottleneck_threshold_betweenness <- quantile(centrality_df$betweenness, 0.80)
bottleneck_threshold_degree <- quantile(centrality_df$degree, 0.50)

bottleneck_genes <- centrality_df %>%
  filter(betweenness >= bottleneck_threshold_betweenness,
         degree < median(degree)) %>%
  arrange(desc(betweenness))

print(paste("Bottleneck genes identified:", nrow(bottleneck_genes)))

if (nrow(bottleneck_genes) > 0) {
  write.csv(bottleneck_genes, "bottleneck_genes.csv", row.names = FALSE)
  
  print("Top 10 bottleneck genes:")
  print(head(bottleneck_genes[, c("gene", "degree", "betweenness", "direction")], 10))
}

# =============================================================================
# 19. CREATE INTEGRATED SUMMARY TABLE
# =============================================================================

print("Creating integrated summary table...")

# Combine all information
integrated_summary <- centrality_df %>%
  mutate(
    node_type = case_when(
      gene %in% hub_genes$gene ~ "Hub",
      gene %in% bottleneck_genes$gene ~ "Bottleneck",
      TRUE ~ "Regular"
    ),
    hub_score = (degree_rank + betweenness_rank + eigenvector_rank) / 3
  )

# Add pathway information if available
if (exists("hub_kegg") && nrow(hub_kegg@result) > 0) {
  # Create a gene-to-pathway mapping
  pathway_mapping <- list()
  for (i in 1:nrow(hub_kegg@result)) {
    genes_in_pathway <- unlist(strsplit(hub_kegg@result$geneID[i], "/"))
    pathway_name <- hub_kegg@result$Description[i]
    for (gene in genes_in_pathway) {
      if (gene %in% names(pathway_mapping)) {
        pathway_mapping[[gene]] <- c(pathway_mapping[[gene]], pathway_name)
      } else {
        pathway_mapping[[gene]] <- pathway_name
      }
    }
  }
  
  integrated_summary$key_pathways <- sapply(integrated_summary$gene, function(g) {
    if (g %in% names(pathway_mapping)) {
      paste(pathway_mapping[[g]], collapse = "; ")
    } else {
      NA
    }
  })
}

write.csv(integrated_summary, "integrated_network_summary.csv", row.names = FALSE)

# Create publication-ready table for top genes
publication_table <- integrated_summary %>%
  filter(node_type %in% c("Hub", "Bottleneck")) %>%
  arrange(desc(degree)) %>%
  head(30) %>%
  select(gene, node_type, degree, betweenness, eigenvector, 
         log2FoldChange, direction, community_louvain) %>%
  mutate(
    degree = round(degree, 0),
    betweenness = round(betweenness, 4),
    eigenvector = round(eigenvector, 4),
    log2FoldChange = round(log2FoldChange, 2)
  )

if (has_survival) {
  publication_table <- publication_table %>%
    left_join(survival_results %>% select(gene, HR, cox_pvalue),
              by = "gene") %>%
    mutate(
      HR = round(HR, 2),
      cox_pvalue = format(cox_pvalue, scientific = TRUE, digits = 2)
    )
}

write.csv(publication_table, "publication_table_network_genes.csv", row.names = FALSE)

# =============================================================================
# 20. FINAL COMPREHENSIVE SUMMARY
# =============================================================================

cat("\n\n===========================================\n")
cat("PPI NETWORK ANALYSIS COMPLETED\n")
cat("===========================================\n\n")

cat("FINAL SUMMARY:\n\n")
cat("Network Overview:\n")
cat(paste("  - Genes analyzed:", vcount(network), "\n"))
cat(paste("  - Interactions:", ecount(network), "\n"))
cat(paste("  - Average degree:", round(mean(degree(network)), 2), "\n"))
cat(paste("  - Network density:", round(network_density, 4), "\n"))
cat(paste("  - Clustering coefficient:", round(transitivity(network, type = "global"), 3), "\n\n"))

cat("Key Genes Identified:\n")
cat(paste("  - Hub genes:", nrow(hub_genes), "\n"))
cat(paste("  - Bottleneck genes:", nrow(bottleneck_genes), "\n"))
cat(paste("  - Network communities:", length(unique(V(network)$community_louvain)), "\n\n"))

if (has_survival && exists("prognostic_hubs")) {
  cat("Survival Integration:\n")
  cat(paste("  - Prognostic hub genes:", nrow(prognostic_hubs), "\n"))
  cat(paste("  - Poor prognosis hubs:", sum(prognostic_hubs$HR > 1 & prognostic_hubs$cox_pvalue < 0.05), "\n"))
  cat(paste("  - Good prognosis hubs:", sum(prognostic_hubs$HR < 1 & prognostic_hubs$cox_pvalue < 0.05), "\n\n"))
}

cat("Top 10 Hub Genes:\n")
top_10_hubs <- head(hub_genes, 10)
for (i in 1:nrow(top_10_hubs)) {
  cat(paste0("  ", i, ". ", top_10_hubs$gene[i], 
             " (Degree: ", top_10_hubs$degree[i], 
             ", ", top_10_hubs$direction[i], ")\n"))
}

cat("\n\nOutput Files Summary:\n")
cat("  Network files: 8\n")
cat("  Hub gene files: 6\n")
cat("  Community files: ", 1 + min(3, nrow(community_analysis)), "\n")
cat("  Visualization files: ", 7 + ifelse(has_survival, 1, 0), "\n")
cat("  Cytoscape files: 3\n")
cat("  Total output files: ", length(list.files(pattern = "*.csv|*.png|*.pdf|*.html|*.sif")), "\n\n")

cat("===========================================\n")
cat("Analysis files saved in: PPI_Network_Analysis/\n")
cat("===========================================\n\n")

# Save workspace
save.image("PPI_network_analysis_workspace.RData")

print("All analyses completed successfully!")
print("Check PPI_Network_Analysis_Report.md for detailed results")

