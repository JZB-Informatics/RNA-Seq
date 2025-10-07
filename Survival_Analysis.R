# Comprehensive Survival Analysis Pipeline
# Author: JahanZaib
# Purpose: Survival analysis using DEG results

# Load required libraries
library(TCGAbiolinks)
library(survival)
library(survminer)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(gridExtra)
library(forestplot)
library(glmnet)
library(ComplexHeatmap)
library(circlize)

# =============================================================================
# 1. LOAD DEG ANALYSIS RESULTS
# =============================================================================

# Convert expression matrix to proper format
gene_info_cols <- c("ensembl_id", "gene_symbol", "entrez_id")
expr_data <- expression_matrix[, !names(expression_matrix) %in% gene_info_cols]
#####--------------------------------------############################
# Check for duplicated values
duplicated_symbols <- duplicated(expression_matrix$gene_symbol)
# View the duplicated symbols
expression_matrix$gene_symbol[duplicated_symbols]
# Apply make.unique to the gene symbol vector
unique_symbols <- make.unique(expression_matrix$gene_symbol)
# Assign the unique symbols as row names
rownames(expr_data) <- unique_symbols
#####--------------------------------------###########################
rownames(expr_data) <- expression_matrix$gene_symbol

# Transpose for patient-centric analysis
expr_data_t <- as.data.frame(t(expr_data))
expr_data_t$sample_id <- rownames(expr_data_t)

print(paste("Expression data dimensions:", nrow(expr_data), "genes x", ncol(expr_data), "samples"))

# =============================================================================
# 2. DOWNLOAD AND PREPARE CLINICAL DATA
# =============================================================================

# Get patient clinical data
clinical_patient <- clinical_data$clinical_patient_brca
write.csv(clinical_patient, "clinical_patient_brca.csv", row.names = FALSE)
clinical_followup <- clinical_data$clinical_follow_up_brca

# Process clinical data
clinical_processed <- clinical_patient %>%
  select(bcr_patient_barcode, gender, age_at_diagnosis, 
         vital_status, death_days_to, last_contact_days_to,
         ajcc_pathologic_tumor_stage, ajcc_tumor_pathologic_pt, ajcc_nodes_pathologic_pn, ajcc_metastasis_pathologic_pm,
         er_status_by_ihc, pr_status_by_ihc, her2_status_by_ihc) %>%
  mutate(
    # Convert age to years
    age_at_diagnosis = as.numeric(age_at_diagnosis) ,
    
    # Create survival time and event variables
    days_to_event = ifelse(vital_status == "Dead", 
                           as.numeric(death_days_to), 
                           as.numeric(last_contact_days_to)),
    survival_time = days_to_event / 365.25,  # Convert to years
    survival_event = ifelse(vital_status == "Dead", 1, 0),
    
    # Clean stage information
    stage_simple = case_when(
      grepl("Stage I", ajcc_pathologic_tumor_stage) ~ "I",
      grepl("Stage II", ajcc_pathologic_tumor_stage) ~ "II",
      grepl("Stage III", ajcc_pathologic_tumor_stage) ~ "III",
      grepl("Stage IV", ajcc_pathologic_tumor_stage) ~ "IV",
      TRUE ~ "Unknown"
    ),
    
    # Create molecular subtype
    subtype = case_when(
      er_status_by_ihc == "Positive" & pr_status_by_ihc == "Positive" & 
        her2_status_by_ihc == "Negative" ~ "Luminal A",
      er_status_by_ihc == "Positive" & her2_status_by_ihc == "Positive" ~ "Luminal B",
      er_status_by_ihc == "Negative" & pr_status_by_ihc == "Negative" & 
        her2_status_by_ihc == "Positive" ~ "HER2+",
      er_status_by_ihc == "Negative" & pr_status_by_ihc == "Negative" & 
        her2_status_by_ihc == "Negative" ~ "TNBC",
      TRUE ~ "Unknown"
    )
  ) %>%
  filter(!is.na(survival_time) & survival_time > 0)

# Match clinical data with expression data
# Extract patient barcode from sample names (first 12 characters)
expr_data_t$patient_barcode <- substr(expr_data_t$sample_id, 1, 12)

# Merge clinical and expression data
surv_data <- merge(clinical_processed, expr_data_t, 
                   by.x = "bcr_patient_barcode", 
                   by.y = "patient_barcode")

print(paste("Patients with both clinical and expression data:", nrow(surv_data)))

# =============================================================================
# 3. SURVIVAL ANALYSIS FOR TOP DEGS
# =============================================================================

print("Performing survival analysis for top DEGs...")

# Initialize results storage
survival_results <- list()
cox_results <- data.frame()

# Function to perform survival analysis for a single gene
analyze_gene_survival <- function(gene_name, data, cutoff_method = "median") {
  
  if (!gene_name %in% names(data)) {
    return(NULL)
  }
  
  # Get gene expression
  gene_expr <- data[[gene_name]]
  
  # Determine cutoff
  if (cutoff_method == "median") {
    cutoff <- median(gene_expr, na.rm = TRUE)
    data$gene_group <- ifelse(gene_expr > cutoff, "High", "Low")
  } else if (cutoff_method == "tertile") {
    tertiles <- quantile(gene_expr, c(0.33, 0.67), na.rm = TRUE)
    data$gene_group <- cut(gene_expr, 
                           breaks = c(-Inf, tertiles, Inf),
                           labels = c("Low", "Medium", "High"))
  } else if (cutoff_method == "optimal") {
    # Use survminer's surv_cutpoint for optimal cutoff
    cutpoint <- surv_cutpoint(data, 
                              time = "survival_time", 
                              event = "survival_event",
                              variables = gene_name)
    cutoff <- cutpoint$cutpoint[[1]]$cutpoint
    data$gene_group <- ifelse(gene_expr > cutoff, "High", "Low")
  }
  
  # Create survival object
  surv_obj <- Surv(data$survival_time, data$survival_event)
  
  # Log-rank test
  survdiff_result <- survdiff(surv_obj ~ gene_group, data = data)
  p_value <- 1 - pchisq(survdiff_result$chisq, length(survdiff_result$n) - 1)
  
  # Cox regression
  cox_model <- coxph(surv_obj ~ gene_expr + age_at_diagnosis + stage_simple, data = data)
  cox_summary <- summary(cox_model)
  
  # Extract results
  result <- list(
    gene = gene_name,
    cutoff = cutoff,
    logrank_pvalue = p_value,
    cox_hr = exp(cox_summary$coefficients[1, "coef"]),
    cox_ci_lower = exp(cox_summary$conf.int[1, "lower .95"]),
    cox_ci_upper = exp(cox_summary$conf.int[1, "upper .95"]),
    cox_pvalue = cox_summary$coefficients[1, "Pr(>|z|)"],
    survfit = survfit(surv_obj ~ gene_group, data = data),
    data = data
  )
  
  return(result)
}

# Analyze top survival candidate genes
top_genes_to_analyze <- head(top_survival_genes$gene_symbol, 50)

for (gene in top_genes_to_analyze) {
  result <- analyze_gene_survival(gene, surv_data, cutoff_method = "median")
  if (!is.null(result)) {
    survival_results[[gene]] <- result
    
    # Store Cox results
    cox_results <- rbind(cox_results, data.frame(
      gene = gene,
      HR = result$cox_hr,
      CI_lower = result$cox_ci_lower,
      CI_upper = result$cox_ci_upper,
      cox_pvalue = result$cox_pvalue,
      logrank_pvalue = result$logrank_pvalue,
      stringsAsFactors = FALSE
    ))
  }
}

# Adjust p-values for multiple testing
cox_results$cox_padj <- p.adjust(cox_results$cox_pvalue, method = "BH")
cox_results$logrank_padj <- p.adjust(cox_results$logrank_pvalue, method = "BH")

# Sort by significance
cox_results <- cox_results[order(cox_results$cox_pvalue), ]

print(paste("Genes with significant survival association (p < 0.05):", 
            sum(cox_results$cox_pvalue < 0.05)))

# =============================================================================
# 4. KAPLAN-MEIER PLOTS FOR TOP PROGNOSTIC GENES
# =============================================================================

print("Generating Kaplan-Meier plots...")

# Select top prognostic genes
top_prognostic <- head(cox_results[cox_results$cox_pvalue < 0.05, ], 10)

# Create directory for plots
dir.create("survival_plots", showWarnings = FALSE)

colnames(result$data)

# Generate KM plots for top genes
for (i in 1:nrow(top_prognostic)) {
  gene <- top_prognostic$gene[i]
  result <- survival_results[[gene]]
  
  if (!is.null(result)) {
    
    # Create KM plot
    km_plot <- ggsurvplot(
      result$survfit,
      data = result$data,
      pval = TRUE,
      pval.method = TRUE,
      conf.int = TRUE,
      risk.table = TRUE,
      risk.table.height = 0.25,
      palette = c("#E7B800", "#2E9FDF"),
      title = paste("Kaplan-Meier Survival Curve:", gene),
      xlab = "Time (years)",
      ylab = "Survival probability",
      legend.title = paste(gene, "Expression"),
      legend.labs = c("Low", "High"),
      ggtheme = theme_bw()
    )
    
    # Add HR to plot
    km_plot$plot <- km_plot$plot + 
      annotate("text", x = 0.5, y = 0.1, 
               label = paste("HR =", round(result$cox_hr, 2),
                           "\n95% CI: [", round(result$cox_ci_lower, 2),
                           "-", round(result$cox_ci_upper, 2), "]"),
               hjust = 0, size = 3)
    
    # Save plot
    pdf(paste0("survival_plots/KM_", gene, ".pdf"), width = 10, height = 8)
    print(km_plot)
    dev.off()
  }
}

# =============================================================================
# 5. MULTIVARIATE COX REGRESSION WITH GENE SIGNATURES
# =============================================================================

print("Building multivariate Cox models...")

# Create gene signature scores
# Upregulated signature
up_genes <- top_survival_genes %>%
  filter(log2FoldChange > 0) %>%
  head(20) %>%
  pull(gene_symbol)

up_genes_present <- up_genes[up_genes %in% names(surv_data)]
surv_data$up_signature <- rowMeans(surv_data[, up_genes_present], na.rm = TRUE)

# Downregulated signature
down_genes <- top_survival_genes %>%
  filter(log2FoldChange < 0) %>%
  head(20) %>%
  pull(gene_symbol)

down_genes_present <- down_genes[down_genes %in% names(surv_data)]
surv_data$down_signature <- rowMeans(surv_data[, down_genes_present], na.rm = TRUE)

# Combined signature score
surv_data$combined_signature <- surv_data$up_signature - surv_data$down_signature

# Multivariate Cox model with signatures
mv_cox <- coxph(Surv(survival_time, survival_event) ~ 
                  combined_signature + age_at_diagnosis + stage_simple + subtype,
                data = surv_data)

mv_cox_summary <- summary(mv_cox)
print(mv_cox_summary)

# Forest plot for multivariate model
forest_data <- data.frame(
  Variable = rownames(mv_cox_summary$conf.int),
  HR = mv_cox_summary$conf.int[, "exp(coef)"],
  Lower = mv_cox_summary$conf.int[, "lower .95"],
  Upper = mv_cox_summary$conf.int[, "upper .95"],
  pvalue = mv_cox_summary$coefficients[, "Pr(>|z|)"]
)

# Create forest plot
forest_plot <- ggplot(forest_data, aes(x = HR, y = Variable)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  scale_x_log10() +
  xlab("Hazard Ratio (95% CI)") +
  ylab("") +
  theme_bw() +
  ggtitle("Multivariate Cox Regression Forest Plot")

ggsave("forest_plot_multivariate.png", forest_plot, width = 10, height = 6, dpi = 300)

# =============================================================================
# 6. PROGNOSTIC GENE IDENTIFICATION BY PATHWAY
# =============================================================================

print("Identifying prognostic genes in enriched pathways...")

# Load pathway results if available
if (file.exists("KEGG_upregulated.csv")) {
  kegg_up <- read.csv("KEGG_upregulated.csv", stringsAsFactors = FALSE)
  
  # Extract genes from top pathways
  if (nrow(kegg_up) > 0) {
    top_pathways <- head(kegg_up, 5)
    pathway_genes <- list()
    
    for (i in 1:nrow(top_pathways)) {
      # Parse gene list from pathway (this depends on the format)
      genes_in_pathway <- unlist(strsplit(top_pathways$geneID[i], "/"))
      pathway_name <- top_pathways$Description[i]
      
      # Convert Entrez IDs to gene symbols
      gene_symbols <- deg_results %>%
        filter(entrez_id %in% genes_in_pathway) %>%
        pull(gene_symbol)
      
      pathway_genes[[pathway_name]] <- gene_symbols
      
      # Test each gene in the pathway
      pathway_survival_results <- data.frame()
      
      for (gene in gene_symbols) {
        if (gene %in% names(surv_data)) {
          result <- analyze_gene_survival(gene, surv_data)
          if (!is.null(result)) {
            pathway_survival_results <- rbind(pathway_survival_results, data.frame(
              pathway = pathway_name,
              gene = gene,
              HR = result$cox_hr,
              pvalue = result$cox_pvalue,
              stringsAsFactors = FALSE
            ))
          }
        }
      }
      
      if (nrow(pathway_survival_results) > 0) {
        pathway_survival_results$padj <- p.adjust(pathway_survival_results$pvalue, method = "BH")
        write.csv(pathway_survival_results, 
                 paste0("pathway_survival_", gsub(" ", "_", pathway_name), ".csv"),
                 row.names = FALSE)
      }
    }
  }
}

# =============================================================================
# 7. LASSO COX REGRESSION FOR FEATURE SELECTION
# =============================================================================

print("Performing LASSO Cox regression for gene selection...")

# Prepare data for LASSO
# Select genes with some variance
gene_vars <- apply(surv_data[, top_genes_to_analyze[top_genes_to_analyze %in% names(surv_data)]], 
                   2, var, na.rm = TRUE)
genes_for_lasso <- names(gene_vars[gene_vars > 0.1])

# Create matrix for LASSO
x_matrix <- as.matrix(surv_data[, genes_for_lasso])
y_matrix <- Surv(surv_data$survival_time, surv_data$survival_event)

# Remove rows with NA
complete_rows <- complete.cases(x_matrix, y_matrix)
x_matrix <- x_matrix[complete_rows, ]
y_matrix <- y_matrix[complete_rows, ]

# Perform LASSO Cox regression with cross-validation
set.seed(123)
cv_lasso <- cv.glmnet(x_matrix, y_matrix, family = "cox", alpha = 1, nfolds = 10)

# Plot cross-validation results
pdf("lasso_cv_plot.pdf", width = 10, height = 6)
plot(cv_lasso)
dev.off()

# Get selected genes at lambda.min
lasso_coef <- coef(cv_lasso, s = "lambda.min")
selected_genes <- rownames(lasso_coef)[lasso_coef[,1] != 0]

print(paste("Genes selected by LASSO:", length(selected_genes)))
print(selected_genes)

# Build Cox model with LASSO-selected genes
if (length(selected_genes) > 0) {
  formula_str <- paste("Surv(survival_time, survival_event) ~", 
                      paste(selected_genes, collapse = " + "))
  lasso_cox <- coxph(as.formula(formula_str), data = surv_data)
  lasso_cox_summary <- summary(lasso_cox)
  
  # Save LASSO results
  lasso_results <- data.frame(
    gene = selected_genes,
    coefficient = lasso_cox_summary$coefficients[, "coef"],
    HR = lasso_cox_summary$coefficients[, "exp(coef)"],
    pvalue = lasso_cox_summary$coefficients[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
  write.csv(lasso_results, "lasso_selected_genes.csv", row.names = FALSE)
}

# =============================================================================
# 8. RISK SCORE CALCULATION AND VALIDATION
# =============================================================================

print("Calculating prognostic risk scores...")

# Calculate risk score based on top prognostic genes
prog_genes <- head(cox_results[cox_results$cox_pvalue < 0.05, "gene"], 10)

if (length(prog_genes) > 0) {
  # Calculate risk score as weighted sum of gene expression
  risk_score <- rep(0, nrow(surv_data))
  
  for (gene in prog_genes) {
    if (gene %in% names(surv_data)) {
      # Get coefficient from Cox model
      gene_hr <- cox_results[cox_results$gene == gene, "HR"][1]
      weight <- log(gene_hr)
      
      # Add weighted expression to risk score
      risk_score <- risk_score + weight * surv_data[[gene]]
    }
  }
  
  surv_data$risk_score <- risk_score
  
  # Stratify patients by risk score
  surv_data$risk_group <- ifelse(surv_data$risk_score > median(surv_data$risk_score, na.rm = TRUE),
                                 "High Risk", "Low Risk")
  
  # Survival analysis for risk groups
  risk_survfit <- survfit(Surv(survival_time, survival_event) ~ risk_group, data = surv_data)
  risk_survdiff <- survdiff(Surv(survival_time, survival_event) ~ risk_group, data = surv_data)
  
  # Plot risk group survival
  risk_plot <- ggsurvplot(
    risk_survfit,
    data = surv_data,
    pval = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    palette = c("#2E9FDF", "#E7B800"),
    title = "Risk Score Stratification",
    xlab = "Time (years)",
    ylab = "Survival probability",
    legend.title = "Risk Group",
    ggtheme = theme_bw()
  )
  
  pdf("risk_score_survival.pdf", width = 10, height = 8)
  print(risk_plot)
  dev.off()
  
  # ROC analysis for risk score
  library(survivalROC)
  
  # Calculate time-dependent ROC
  roc_1year <- survivalROC(
    Stime = surv_data$survival_time,
    status = surv_data$survival_event,
    marker = surv_data$risk_score,
    predict.time = 1,
    method = "KM"
  )
  
  roc_3year <- survivalROC(
    Stime = surv_data$survival_time,
    status = surv_data$survival_event,
    marker = surv_data$risk_score,
    predict.time = 3,
    method = "KM"
  )
  
  roc_5year <- survivalROC(
    Stime = surv_data$survival_time,
    status = surv_data$survival_event,
    marker = surv_data$risk_score,
    predict.time = 5,
    method = "KM"
  )
  
  # Plot ROC curves
  pdf("roc_curves.pdf", width = 8, height = 8)
  plot(roc_1year$FP, roc_1year$TP, type = "l", col = "red", lwd = 2,
       xlim = c(0, 1), ylim = c(0, 1),
       xlab = "1 - Specificity", ylab = "Sensitivity",
       main = "Time-Dependent ROC Curves")
  lines(roc_3year$FP, roc_3year$TP, col = "blue", lwd = 2)
  lines(roc_5year$FP, roc_5year$TP, col = "green", lwd = 2)
  abline(0, 1, lty = 2)
  legend("bottomright", c(paste("1-year AUC =", round(roc_1year$AUC, 3)),
                          paste("3-year AUC =", round(roc_3year$AUC, 3)),
                          paste("5-year AUC =", round(roc_5year$AUC, 3))),
         col = c("red", "blue", "green"), lty = 1, lwd = 2)
  dev.off()
}

# =============================================================================
# 9. SUBGROUP ANALYSIS
# =============================================================================

print("Performing subgroup survival analysis...")

# Analyze top genes by molecular subtype
subtypes <- unique(surv_data$subtype[surv_data$subtype != "Unknown"])
subtype_results <- list()

for (subtype in subtypes) {
  subtype_data <- surv_data[surv_data$subtype == subtype, ]
  
  if (nrow(subtype_data) > 20) {  # Minimum sample size
    subtype_cox <- data.frame()
    
    for (gene in head(top_genes_to_analyze, 20)) {
      if (gene %in% names(subtype_data)) {
        result <- analyze_gene_survival(gene, subtype_data)
        if (!is.null(result)) {
          subtype_cox <- rbind(subtype_cox, data.frame(
            subtype = subtype,
            gene = gene,
            HR = result$cox_hr,
            pvalue = result$cox_pvalue,
            n_patients = nrow(subtype_data),
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    
    subtype_results[[subtype]] <- subtype_cox
  }
}

# Combine and save subtype results
all_subtype_results <- do.call(rbind, subtype_results)
write.csv(all_subtype_results, "subtype_survival_analysis.csv", row.names = FALSE)

# =============================================================================
# 10. COMPREHENSIVE HEATMAP OF SURVIVAL ASSOCIATIONS
# =============================================================================

print("Creating comprehensive survival heatmap...")

# Create matrix of HR values for visualization
hr_matrix <- matrix(NA, nrow = length(top_genes_to_analyze), ncol = 1)
rownames(hr_matrix) <- top_genes_to_analyze
colnames(hr_matrix) <- "Hazard Ratio"

for (i in 1:length(top_genes_to_analyze)) {
  gene <- top_genes_to_analyze[i]
  if (gene %in% cox_results$gene) {
    hr_matrix[i, 1] <- log2(cox_results[cox_results$gene == gene, "HR"][1])
  }
}

# Remove NA rows
hr_matrix <- hr_matrix[!is.na(hr_matrix[,1]), , drop = FALSE]

# Add significance annotations
sig_genes <- rownames(hr_matrix) %in% cox_results[cox_results$cox_pvalue < 0.05, "gene"]
sig_annotation <- ifelse(sig_genes, "*", "")

# Create heatmap
pdf("survival_heatmap.pdf", width = 6, height = 12)
Heatmap(hr_matrix,
        name = "log2(HR)",
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        cell_fun = function(j, i, x, y, width, height, fill) {
          if (sig_genes[i]) {
            grid.text("*", x, y, gp = gpar(fontsize = 10))
          }
        },
        column_title = "Survival Association Heatmap",
        row_title = "Genes",
        clustering_distance_rows = "euclidean",
        show_row_names = TRUE,
        row_names_gp = gpar(fontsize = 8))
dev.off()

# =============================================================================
# 11. SAVE COMPREHENSIVE RESULTS
# =============================================================================

print("Saving comprehensive survival analysis results...")

# Create summary report
summary_report <- list(
  n_patients = nrow(surv_data),
  median_followup = median(surv_data$survival_time[surv_data$survival_event == 0]),
  n_events = sum(surv_data$survival_event),
  n_genes_tested = length(top_genes_to_analyze),
  n_significant_genes = sum(cox_results$cox_pvalue < 0.05),
  top_prognostic_genes = head(cox_results[cox_results$cox_pvalue < 0.05, c("gene", "HR", "cox_pvalue")], 10)
)

# Save all results
write.csv(cox_results, "cox_regression_results_all_genes.csv", row.names = FALSE)
saveRDS(survival_results, "survival_results_object.rds")
saveRDS(summary_report, "survival_summary_report.rds")

# Create markdown report
report_text <- paste0(
  "# Survival Analysis Report\n\n",
  "## Dataset: TCGA-BRCA\n",
  "Date: ", Sys.Date(), "\n\n",
  "### Patient Statistics\n",
  "- Total patients analyzed: ", summary_report$n_patients, "\n",
  "- Median follow-up (years): ", round(summary_report$median_followup, 2), "\n",
  "- Number of events: ", summary_report$n_events, "\n\n",
  "### Gene Analysis\n",
  "- Genes tested: ", summary_report$n_genes_tested, "\n",
  "- Significant prognostic genes: ", summary_report$n_significant_genes, "\n\n",
  "### Top Prognostic Genes\n"
)

# Add top genes to report
for (i in 1:min(10, nrow(summary_report$top_prognostic_genes))) {
  gene_info <- summary_report$top_prognostic_genes[i, ]
  report_text <- paste0(report_text,
    i, ". ", gene_info$gene, 
    " (HR: ", round(gene_info$HR, 2), 
    ", p-value: ", format(gene_info$cox_pvalue, scientific = TRUE, digits = 2), ")\n")
}

if (length(selected_genes) > 0) {
  report_text <- paste0(report_text,
    "\n### LASSO Selected Genes\n",
    "Number of genes selected: ", length(selected_genes), "\n",
    "Genes: ", paste(selected_genes, collapse = ", "), "\n")
}

writeLines(report_text, "survival_analysis_report.md")

# =============================================================================
# 12. FINAL SUMMARY TABLE
# =============================================================================

# Create comprehensive results table
final_results <- merge(
  cox_results,
  deg_results[, c("gene_symbol", "log2FoldChange", "padj")],
  by.x = "gene",
  by.y = "gene_symbol",
  all.x = TRUE
)

names(final_results)[names(final_results) == "padj"] <- "deg_padj"

# Add categories
final_results$prognostic_significance <- case_when(
  final_results$cox_pvalue < 0.001 ~ "Highly Significant",
  final_results$cox_pvalue < 0.01 ~ "Significant",
  final_results$cox_pvalue < 0.05 ~ "Moderately Significant",
  TRUE ~ "Not Significant"
)

final_results$expression_pattern <- ifelse(final_results$log2FoldChange > 0, 
                                          "Upregulated in Tumor", 
                                          "Downregulated in Tumor")

final_results$risk_association <- ifelse(final_results$HR > 1, 
                                        "Poor Prognosis", 
                                        "Good Prognosis")

# Order by significance
final_results <- final_results[order(final_results$cox_pvalue), ]

# Save final comprehensive table
write.csv(final_results, "comprehensive_survival_results.csv", row.names = FALSE)

# Create publication-ready table
publication_table <- final_results %>%
  select(gene, log2FoldChange, deg_padj, HR, CI_lower, CI_upper, 
         cox_pvalue, prognostic_significance, risk_association) %>%
  head(20) %>%
  mutate(
    log2FoldChange = round(log2FoldChange, 2),
    deg_padj = format(deg_padj, scientific = TRUE, digits = 2),
    HR_CI = paste0(round(HR, 2), " (", round(CI_lower, 2), "-", round(CI_upper, 2), ")"),
    cox_pvalue = format(cox_pvalue, scientific = TRUE, digits = 2)
  ) %>%
  select(gene, log2FoldChange, deg_padj, HR_CI, cox_pvalue, 
         prognostic_significance, risk_association)

write.csv(publication_table, "publication_ready_table.csv", row.names = FALSE)

print("===========================================")
print("SURVIVAL ANALYSIS COMPLETED SUCCESSFULLY")
print("===========================================")
print(paste("Total genes analyzed:", length(top_genes_to_analyze)))
print(paste("Significant prognostic genes:", sum(cox_results$cox_pvalue < 0.05)))
print(paste("Output files generated:", length(list.files(pattern = "*.csv|*.pdf|*.png"))))

