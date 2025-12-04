###### Immune Infiltration Analysis Section - GSE78097 Dataset ######
# Set working directory to target folder
setwd("D:/dualdisease/IMMUNEreview")
cat("Working directory set to:", getwd(), "\n")

# Load required packages
library(e1071)
library(preprocessCore)
library(parallel)
library(usethis)
library(devtools)
library(bseqsc)
library(csSAM)
library(Rcpp)
library(CIBERSORT)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(IOBR)
library(ggsci)
library(tidyr)
library(ggpubr)
library(rstatix)
library(xCell)
library(viridis)
library(scico)
library(Hmisc)
library(tidyverse)
library(corrplot)
library(linkET)
library(dplyr)

# Create function to save PDF
save_pdf <- function(plot_obj, filename, width = 10, height = 8) {
  pdf_file <- paste0(filename, ".pdf")
  ggsave(pdf_file, plot = plot_obj, width = width, height = height)
  cat("PDF saved to:", pdf_file, "\n")
}

# Heatmap saving function
save_pheatmap_pdf <- function(heatmap_obj, filename, width = 10, height = 8) {
  pdf_file <- paste0(filename, ".pdf")
  pdf(pdf_file, width = width, height = height)
  print(heatmap_obj)
  dev.off()
  cat("Heatmap PDF saved to:", pdf_file, "\n")
}

# Focus only on these 7 key genes
key_genes <- c("IFI6", "MX1", "NMI", "OAS3", "OASL", "SAMD9", "UBE2L6")

###### Data Processing Section - GSE78097 Dataset ######
cat("Starting GSE78097 data processing...\n")

# 1. Load raw data
load("D:/dualdisease/WGCNA/Preal.data/P_mexp_clin.RData")

# 2. Extract sample IDs of GSE78097 from clin
# Format of SAMPLE_ID for GSE78097 samples in clin: GSE78097_GSM2066635
gse78097_sample_ids <- clin$SAMPLE_ID[grepl("GSE78097", clin$SAMPLE_ID)]
cat("Number of GSE78097 samples:", length(gse78097_sample_ids), "\n")

# 3. Extract GSE78097 data from mexp
# Get all sample IDs from mexp (from Sample column)
all_mexp_samples <- mexp$Sample
cat("Total samples in mexp:", length(all_mexp_samples), "\n")

# Match GSE78097 samples - need to clean sample ID format
clean_gse78097_samples <- gsub("GSE78097_", "", gse78097_sample_ids)
matched_samples <- intersect(all_mexp_samples, clean_gse78097_samples)
cat("Number of matched GSE78097 samples:", length(matched_samples), "\n")

if (length(matched_samples) == 0) {
  stop("Error: No matching GSE78097 samples found!")
}

# Extract rows for GSE78097 samples
mexp_gse78097 <- mexp[mexp$Sample %in% matched_samples, ]
cat("Dimensions of GSE78097 expression matrix:", dim(mexp_gse78097), "\n")

# 4. Prepare expression matrix data
# Convert tibble to data.frame
dataExpr <- as.data.frame(mexp_gse78097[, -1])  # Remove first column "Sample"
rownames(dataExpr) <- mexp_gse78097$Sample

cat("Dimensions of processed expression matrix:", dim(dataExpr), "\n")
cat("Number of samples:", nrow(dataExpr), "Number of genes:", ncol(dataExpr), "\n")

# 5. Extract accurate clinical information of GSE78097 from clin
# Create cleaned sample IDs for matching
clin_clean <- clin
clin_clean$sampleid_clean <- gsub("GSE[0-9]+_", "", clin_clean$SAMPLE_ID)

# Match GSE78097 samples - using cleaned sample IDs
clin_gse78097 <- clin_clean[clin_clean$sampleid_clean %in% rownames(dataExpr), ]
cat("Number of GSE78097 samples found in clin:", nrow(clin_gse78097), "\n")

# 6. Check if key genes exist
available_genes <- intersect(key_genes, colnames(dataExpr))
missing_genes <- setdiff(key_genes, colnames(dataExpr))

cat("Available key genes:", length(available_genes), "\n")
print(available_genes)
cat("Missing key genes:", length(missing_genes), "\n")
print(missing_genes)

if (length(available_genes) == 0) {
  stop("Error: No key genes found!")
}

# 7. Data normalization
cat("Starting data normalization...\n")

# Check and remove zero variance genes
zero_var_genes <- apply(dataExpr, 2, function(x) var(x, na.rm = TRUE) == 0)
if (any(zero_var_genes)) {
  cat("Removing", sum(zero_var_genes), "zero variance genes\n")
  dataExpr <- dataExpr[, !zero_var_genes]
}

# Transpose and normalize (rows are genes, columns are samples)
dataExpr_t <- t(dataExpr)
data_scaled <- t(scale(t(dataExpr_t)))
dataExpr_scaled <- t(data_scaled)

cat("Dimensions of normalized data:", dim(dataExpr_scaled), "\n")

# Save normalized data
write.csv(dataExpr_scaled, "GSE78097_Scaled_Gene_Expression.csv")
cat("Normalized data saved to: GSE78097_Scaled_Gene_Expression.csv\n")

# 8. Save clinical information
write.csv(clin_gse78097, "GSE78097_clinical_info.csv", row.names = FALSE)
cat("Clinical information saved to: GSE78097_clinical_info.csv\n")

cat("GSE78097 data processing completed!\n")

###### Immune Infiltration Analysis Section ######
cat("Starting GSE78097 immune infiltration analysis...\n")

# 9. Perform immune infiltration analysis using multiple methods
cat("Starting immune infiltration analysis...\n")

# Ensure correct data format: rows are genes, columns are samples
# Current dataExpr_scaled has rows as samples, columns as genes, need transposition
eset_for_immune <- t(dataExpr_scaled)
cat("Dimensions of input data for immune analysis:", dim(eset_for_immune), "\n")
cat("Number of rows (genes):", nrow(eset_for_immune), "Number of columns (samples):", ncol(eset_for_immune), "\n")

# Convert gene names to uppercase to match reference database
rownames(eset_for_immune) <- toupper(rownames(eset_for_immune))
cat("Gene names converted to uppercase\n")

# Check if key genes still exist after conversion
available_genes_upper <- intersect(toupper(available_genes), rownames(eset_for_immune))
cat("Available key genes after conversion:", length(available_genes_upper), "\n")

# Use tryCatch to handle potentially failed analysis methods
cibersort_gse78097 <- NULL
epic_gse78097 <- NULL
mcp_gse78097 <- NULL
xcell_df_gse78097 <- NULL
estimate_gse78097 <- NULL
timer_gse78097 <- NULL
quantiseq_gse78097 <- NULL
ips_res_gse78097 <- NULL

# CIBERSORT analysis - use tryCatch to handle potential errors
cat("Performing CIBERSORT analysis...\n")
tryCatch({
  cibersort_gse78097 <- deconvo_tme(eset = eset_for_immune, method = "cibersort", arrays = FALSE, perm = 100)
  cat("CIBERSORT analysis succeeded\n")
}, error = function(e) {
  cat("CIBERSORT analysis failed:", e$message, "\n")
})

# EPIC analysis
cat("Performing EPIC analysis...\n")
tryCatch({
  epic_gse78097 <- deconvo_tme(eset = eset_for_immune, method = "epic", arrays = FALSE)
  cat("EPIC analysis succeeded\n")
}, error = function(e) {
  cat("EPIC analysis failed:", e$message, "\n")
})

# MCP-counter analysis
cat("Performing MCP-counter analysis...\n")
tryCatch({
  mcp_gse78097 <- deconvo_tme(eset = eset_for_immune, method = "mcpcounter")
  cat("MCP-counter analysis succeeded\n")
}, error = function(e) {
  cat("MCP-counter analysis failed:", e$message, "\n")
})

# xCELL analysis
cat("Performing xCELL analysis...\n")
tryCatch({
  xcell_result_gse78097 <- xCellAnalysis(eset_for_immune)
  xcell_df_gse78097 <- as.data.frame(t(xcell_result_gse78097))
  xcell_df_gse78097$ID <- rownames(xcell_df_gse78097)
  rownames(xcell_df_gse78097) <- NULL
  cat("xCELL analysis succeeded\n")
}, error = function(e) {
  cat("xCELL analysis failed:", e$message, "\n")
})

# ESTIMATE analysis
cat("Performing ESTIMATE analysis...\n")
tryCatch({
  estimate_gse78097 <- deconvo_tme(eset = eset_for_immune, method = "estimate")
  cat("ESTIMATE analysis succeeded\n")
}, error = function(e) {
  cat("ESTIMATE analysis failed:", e$message, "\n")
})

# TIMER analysis
cat("Performing TIMER analysis...\n")
tryCatch({
  timer_gse78097 <- deconvo_tme(eset = eset_for_immune, method = "timer", group_list = rep("other", ncol(eset_for_immune)))
  cat("TIMER analysis succeeded\n")
}, error = function(e) {
  cat("TIMER analysis failed:", e$message, "\n")
})

# quanTIseq analysis
cat("Performing quanTIseq analysis...\n")
tryCatch({
  quantiseq_gse78097 <- deconvo_tme(eset = eset_for_immune, tumor = TRUE, arrays = FALSE, scale_mrna = TRUE, method = "quantiseq")
  cat("quanTIseq analysis succeeded\n")
}, error = function(e) {
  cat("quanTIseq analysis failed:", e$message, "\n")
})

# IPS analysis
cat("Performing IPS analysis...\n")
tryCatch({
  ips_res_gse78097 <- deconvo_tme(eset = eset_for_immune, method = "ips", plot = F)
  cat("IPS analysis succeeded\n")
}, error = function(e) {
  cat("IPS analysis failed:", e$message, "\n")
})

# 10. Merge all immune infiltration results
cat("Merging immune infiltration results...\n")

# Create results list
results_list <- list()

# Add available results
if (!is.null(mcp_gse78097)) results_list[["mcp"]] <- mcp_gse78097
if (!is.null(xcell_df_gse78097)) results_list[["xcell"]] <- xcell_df_gse78097
if (!is.null(epic_gse78097)) results_list[["epic"]] <- epic_gse78097
if (!is.null(estimate_gse78097)) results_list[["estimate"]] <- estimate_gse78097
if (!is.null(timer_gse78097)) results_list[["timer"]] <- timer_gse78097
if (!is.null(quantiseq_gse78097)) results_list[["quantiseq"]] <- quantiseq_gse78097
if (!is.null(ips_res_gse78097)) results_list[["ips"]] <- ips_res_gse78097

# If CIBERSORT results exist, process separately (special column name handling required)
if (!is.null(cibersort_gse78097)) {
  # CIBERSORT results processing
  colnames(cibersort_gse78097) <- gsub("_CIBERSORT", "", colnames(cibersort_gse78097))
  cibersort_gse78097 <- cibersort_gse78097[,-c(24:26)]
  names(cibersort_gse78097) <- c("ID","B_cells_naive","B_cells_memory","Plasma_cells",
                                 "T_cells_CD8","T_cells_CD4_naive","T_cells_CD4_memory_resting",
                                 "T_cells_CD4_memory_activated","T_cells_follicular_helper",
                                 "T_cells_regulatory_(Tregs)","T_cells_gamma_delta","NK_cells_resting",
                                 "NK_cells_activated","Monocytes","Macrophages_M0","Macrophages_M1",
                                 "Macrophages_M2","Dendritic_cells_resting","Dendritic_cells_activated",
                                 "Mast_cells_resting","Mast_cells_activated","Eosinophils","Neutrophils")
  
  cibersort_gse78097$ID <- sub(".*_", "", cibersort_gse78097$ID)
  results_list[["cibersort"]] <- cibersort_gse78097
}

# Merge available results step by step
if (length(results_list) > 0) {
  tme_combine_gse78097 <- results_list[[1]]
  
  if (length(results_list) > 1) {
    for (i in 2:length(results_list)) {
      tme_combine_gse78097 <- tme_combine_gse78097 %>%
        inner_join(results_list[[i]], by = "ID")
    }
  }
  
  save(tme_combine_gse78097, file = "GSE78097_tme_combine_scale.RData")
  cat("Immune infiltration analysis completed! Successfully merged results from", length(results_list), "methods\n")
  cat("Results saved to: GSE78097_tme_combine_scale.RData\n")
} else {
  stop("Error: All immune infiltration analysis methods failed!")
}

# 11. Prepare clinical group information
conditions_gse78097 <- clin_gse78097
conditions_gse78097$ID <- conditions_gse78097$sampleid_clean

# Create groups based on DISEASE column
conditions_gse78097$Group <- ifelse(conditions_gse78097$DISEASE %in% c("Control", "healthy", "normal"), "healthy", "disease")
cat("GSE78097 group distribution:\n")
print(table(conditions_gse78097$Group))

# 12. If CIBERSORT results exist, perform subsequent processing
if (!is.null(cibersort_gse78097)) {
  # Sample matching
  cat("Number of samples in cibersort_gse78097:", nrow(cibersort_gse78097), "\n")
  cat("Number of samples in conditions_gse78097:", nrow(conditions_gse78097), "\n")
  
  # Use merge for exact matching
  cibersort_gse78097_with_group <- merge(cibersort_gse78097, 
                                         conditions_gse78097[, c("ID", "Group")], 
                                         by = "ID", all.x = TRUE)
  
  cibersort_gse78097_with_group$Group <- factor(cibersort_gse78097_with_group$Group, levels = c("healthy", "disease"))
  
  # Use merged data
  res_cibersort_gse78097 <- cibersort_gse78097_with_group
} else {
  cat("Warning: No CIBERSORT results, skipping CIBERSORT-related visualization\n")
  res_cibersort_gse78097 <- NULL
}

# 13. Plot stacked bar chart for GSE78097
cat("Plotting stacked bar chart...\n")
violin_dat_gse78097 <- res_cibersort_gse78097 %>% 
  pivot_longer(
    cols = -c("ID", "Group"),
    names_to = "ImmuneCell",
    values_to = "Score"
  )

mypalette <- c("#FF0000","#D8D8BF","#8E236B","#EAADEA","#BC8F8F","#5959AB","#0000FF","#2F4F4F","#3232CD","#FF00FF","#EAEAAE","#CC3299","#9370DB","#FF6EC7","#545454","#E47833","#856363","#E6E8FA","#3299CC","#9F5F9F","#8E2323","#007FFF","#B5A642","#DB7093","#FF1CAE","#D9D919","#CD7F32","#A68064","#A67D3D","#2F2F4F","#236B8E","#8C7853","#5F9F9F")

p_gse78097 <- ggplot(violin_dat_gse78097, aes(x=ID, y=Score, fill=ImmuneCell)) + 
  geom_bar(position = 'stack', stat = 'identity') + 
  scale_y_continuous(expand = c(0,0)) + 
  theme_bw() + 
  labs(x = '', y = 'Relative Percent', fill = '') + 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        legend.position = 'top') + 
  scale_fill_manual(values = mypalette) + 
  facet_grid(~Group, scales = "free", space = "free")

ggsave("GSE78097_cibersort_stacked_barplot.pdf", p_gse78097, width = 12, height = 8)
write.csv(res_cibersort_gse78097, "GSE78097_immune_cell_cibersort.csv")
cat("Stacked bar chart saved to: GSE78097_cibersort_stacked_barplot.pdf\n")

# 14. Plot boxplot for GSE78097
cat("Plotting boxplot...\n")
b_gse78097 <- gather(res_cibersort_gse78097, key = "Cibersort", value = "Score", -c(Group, ID))

box_plot_gse78097 <- ggboxplot(b_gse78097, x = "Cibersort", y = "Score", 
                               fill = "Group", palette = c("#8680AF","#eaa9c4")) + 
  stat_compare_means(aes(group = Group), 
                     method = "wilcox.test", 
                     label = "p.signif", 
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("***", "**", "*", "ns"))) +
  labs(x = '', y = 'Relative Percent', fill = '') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave("GSE78097_cibersort_boxplot.pdf", box_plot_gse78097, width = 14, height = 8)
cat("Boxplot saved to: GSE78097_cibersort_boxplot.pdf\n")

# 15. Correlation heatmap analysis between key genes and immune infiltration
cat("Performing correlation analysis between key genes and immune infiltration...\n")
load("GSE78097_tme_combine_scale.RData")

# Extract key gene expression data
eset_genes_gse78097 <- dataExpr_scaled[, available_genes, drop = FALSE]

# Preprocessing
gene_sds <- apply(eset_genes_gse78097, 2, sd, na.rm = TRUE)
valid_genes_gse78097 <- names(gene_sds)[gene_sds > 0]
eset_genes_gse78097 <- eset_genes_gse78097[, valid_genes_gse78097, drop = FALSE]

cat("Number of valid key genes in GSE78097:", length(valid_genes_gse78097), "\n")
cat("Invalid key genes in GSE78097:", setdiff(available_genes, valid_genes_gse78097), "\n")

expr_t_gse78097 <- as.data.frame(eset_genes_gse78097)
expr_t_gse78097$ID <- rownames(expr_t_gse78097)

merged_data_gse78097 <- merge(expr_t_gse78097, tme_combine_gse78097, by = "ID")

cat("Dimensions of merged data for GSE78097:", dim(merged_data_gse78097), "\n")

# Calculate correlation matrix
full_cor_matrix_gse78097 <- cor(
  merged_data_gse78097[, valid_genes_gse78097],
  merged_data_gse78097[, grep("_CIBERSORT|_MCPcounter|_EPIC|_estimate|_TIMER|_quantiseq|_IPS|_xCell", 
                              colnames(tme_combine_gse78097))],
  use = "pairwise.complete.obs"
)

cat("Dimensions of correlation matrix for GSE78097:", dim(full_cor_matrix_gse78097), "\n")

# Create annotations
annotation_col_gse78097 <- data.frame(
  Method = sapply(colnames(full_cor_matrix_gse78097), function(x) {
    if (grepl("_CIBERSORT$", x)) "CIBERSORT"
    else if (grepl("_MCPcounter$", x)) "MCPcounter"
    else if (grepl("_EPIC$", x)) "EPIC"
    else if (grepl("_estimate$", x)) "ESTIMATE"
    else if (grepl("_TIMER$", x)) "TIMER"
    else if (grepl("_quantiseq$", x)) "quanTIseq"
    else if (grepl("_IPS$", x)) "IPS"
    else "xCell"
  })
)
rownames(annotation_col_gse78097) <- colnames(full_cor_matrix_gse78097)

# Color settings
col_scale_mako <- viridis::mako(13, direction = -1)[1:10]
col_scale_acton <- scico::scico(10, direction = -1, palette = "acton")
my_color <- c(rev(col_scale_mako[1:8]), "white", col_scale_acton[1:8])

abs_max_val <- max(abs(full_cor_matrix_gse78097), na.rm = TRUE)
ph_breaks <- seq(-abs_max_val, abs_max_val, length.out = length(my_color))

method_colors <- setNames(
  c("#66c2a5", "#fc8d62", "#a1c9f4", "#e78ac3", "#a6d854", "#ffd92f", "#e5c494", "#a999c9"),
  c("CIBERSORT", "MCPcounter", "EPIC", "ESTIMATE", "TIMER", "quanTIseq", "IPS", "xCell")
)

annotation_colors <- list(Method = method_colors)
mat_t_gse78097 <- t(full_cor_matrix_gse78097)
annotation_row_gse78097 <- annotation_col_gse78097

# Create and save heatmap
heatmap_obj_gse78097 <- pheatmap(
  mat = mat_t_gse78097,
  color = my_color,
  breaks = ph_breaks,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  cellwidth = 15,
  cellheight = 8,
  annotation_row = annotation_row_gse78097,
  annotation_colors = annotation_colors,
  fontsize_row = 10,
  fontsize_col = 10,
  angle_col = "45",
  na_col = "grey80",
  main = "GSE78097 - 7 Key Genes Immune Infiltration Correlation"
)

pdf("GSE78097_7key_genes_immune_correlation_heatmap.pdf", width = 12, height = 12)
print(heatmap_obj_gse78097)
dev.off()
cat("Correlation heatmap saved to: GSE78097_7key_genes_immune_correlation_heatmap.pdf\n")

# 16. Additional correlation analysis
# Spearman correlation analysis
cat("Performing Spearman correlation analysis...\n")
df_gse78097 <- merged_data_gse78097[, c(valid_genes_gse78097, setdiff(colnames(tme_combine_gse78097), "ID"))]
data_spearman <- rcorr(as.matrix(df_gse78097), type = "spearman")
r_spearman <- data_spearman$r
p_spearman <- data_spearman$P

r_spearman_sub <- r_spearman[1:length(valid_genes_gse78097), (length(valid_genes_gse78097)+1):ncol(df_gse78097)]
p_spearman_sub <- p_spearman[1:length(valid_genes_gse78097), (length(valid_genes_gse78097)+1):ncol(df_gse78097)]

significant_pairs <- data.frame(
  Gene = rep(rownames(r_spearman_sub), each = ncol(r_spearman_sub)),
  Immune_Method = rep(colnames(r_spearman_sub), times = nrow(r_spearman_sub)),
  r = as.vector(r_spearman_sub),
  p = as.vector(p_spearman_sub)
) %>% filter(abs(r) > 0.6 & p < 0.05)

write.csv(significant_pairs, "GSE78097_7key_genes_significant_correlations.csv", row.names = FALSE)

# Check if there are significant correlation pairs
if(nrow(significant_pairs) > 0) {
  genes_keep <- unique(significant_pairs$Gene)
  immune_keep <- unique(significant_pairs$Immune_Method)
  
  r_spearman_filtered <- r_spearman_sub[genes_keep, immune_keep, drop = FALSE]
  p_spearman_filtered <- p_spearman_sub[genes_keep, immune_keep, drop = FALSE]
  
  my_col <- colorRampPalette(c("#9BCD9B","#B4EEB4", "gray98", "#F2B379","#DD5F60"))(100)
  
  # Save Spearman correlation heatmap as PDF
  pdf("GSE78097_7key_genes_spearman_correlation_heatmap.pdf", width = 10, height = 10)
  corrplot(r_spearman_filtered, 
           p.mat = p_spearman_filtered, 
           insig = "label_sig",
           sig.level = c(0.001, 0.01, 0.05), 
           pch.cex = 0.8, 
           pch.col = "black",
           tl.col = "black", 
           tl.cex = 0.8,
           tl.srt = 45,
           col = my_col,
           method = "square",
           mar = c(4,4,4,2),
           main = "GSE78097 - 7 Key Genes Spearman Correlation (r > 0.6 & p < 0.05)")
  dev.off()
  cat("Spearman correlation heatmap PDF saved to: GSE78097_7key_genes_spearman_correlation_heatmap.pdf\n")
} else {
  cat("Warning: No significant correlation pairs found (r > 0.6 & p < 0.05)\n")
}

# 17. Mantel test and butterfly plot analysis - GSE78097
cat("Starting Mantel test and butterfly plot analysis...\n")

# Data preparation
riskscore <- merged_data_gse78097[, valid_genes_gse78097, drop = FALSE]
rownames(riskscore) <- merged_data_gse78097$ID

# Extract immune cell data from merged data
data <- tme_combine_gse78097
data <- data[data$ID %in% rownames(riskscore), ]
data <- as.data.frame(data)
rownames(data) <- data$ID
data <- data[, -which(colnames(data) == "ID"), drop = FALSE]
data <- data[rownames(riskscore), , drop = FALSE]

# Create spec_select list containing only valid key genes
spec_select_list <- list()
for(i in seq_along(valid_genes_gse78097)) {
  spec_select_list[[valid_genes_gse78097[i]]] <- i
}

# Mantel test
cat("Performing Mantel test, please wait...\n")
mantel <- mantel_test(
  riskscore, 
  data,
  method = "spearman",
  spec_select = spec_select_list,
  permutations = 100
) %>% 
  mutate(
    rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
             labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
    pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
             labels = c("< 0.01", "0.01 - 0.05", ">= 0.05"))
  )

cat("Mantel test completed!\n")

# Filter significant results
mantel_filtered <- mantel %>% filter(r >= 0.4 | p < 0.05)
write.csv(mantel_filtered, "GSE78097_7key_genes_mantel_filtered.csv", row.names = FALSE)

cat("Found", nrow(mantel_filtered), "significant associations\n")

# Create butterfly plot (if there are significant results)
if(nrow(mantel_filtered) > 0) {
  immune_cells <- unique(mantel_filtered$env)
  data_filtered <- data[, colnames(data) %in% immune_cells, drop = FALSE]
  
  # Calculate correlation matrix
  corr_matrix <- correlate(data_filtered)
  
  # Color settings
  mycolor1 <- c("#eb7e60", "#6ca3d4", "grey90")
  
  # Select all associations (including non-significant) related to these immune cells from complete mantel data
  mantel_subset <- mantel %>% filter(env %in% immune_cells)
  
  # Fix butterfly plot
  p4 <- qcorrplot(corr_matrix, 
                  type = "upper", 
                  diag = FALSE,
                  show.diag = FALSE) +
    geom_tile() +
    geom_couple(
      aes(colour = pd, size = rd), 
      data = mantel_subset,
      curvature = 0.1,
      nudge_x = 0.1
    ) +
    scale_fill_gradient2(
      low = "#2166AC", 
      mid = "white", 
      high = "#B2182B",
      midpoint = 0,
      name = "Pearson's r"
    ) +  
    scale_colour_manual(values = mycolor1, name = "Mantel's p") +
    scale_size_manual(
      values = c(0.3, 0.6, 1.2), 
      name = "Mantel's r"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(angle = 0, hjust = 1),
      legend.position = "right"
    ) +
    labs(title = "GSE78097 - 7 Key Genes Immune Cell Correlation Network")
  
  # Save plot
  tryCatch({
    ggsave("GSE78097_7key_genes_mantel_network.pdf", p4, width = 16, height = 12)
    cat("Butterfly plot saved to: GSE78097_7key_genes_mantel_network.pdf\n")
  }, error = function(e) {
    # Alternative saving method
    pdf("GSE78097_7key_genes_mantel_network.pdf", width = 16, height = 12)
    print(p4)
    dev.off()
    cat("Butterfly plot saved to: GSE78097_7key_genes_mantel_network.pdf (using alternative method)\n")
  })
  
} else {
  cat("No significant Mantel associations found, skipping butterfly plot generation\n")
}

# 18. Summary output
cat("\n=== GSE78097 Analysis Completion Summary ===\n")
cat("All results saved to:", getwd(), "\n")
cat("Output files:\n")
cat("1. GSE78097_Scaled_Gene_Expression.csv - Normalized expression data\n")
cat("2. GSE78097_clinical_info.csv - Clinical information\n")
cat("3. GSE78097_tme_combine_scale.RData - Immune infiltration results\n")
cat("4. GSE78097_immune_cell_cibersort.csv - Detailed CIBERSORT results\n")
cat("5. GSE78097_cibersort_stacked_barplot.pdf - Stacked bar chart\n")
cat("6. GSE78097_cibersort_boxplot.pdf - Boxplot\n")
cat("7. GSE78097_7key_genes_immune_correlation_heatmap.pdf - Correlation heatmap\n")
cat("8. GSE78097_7key_genes_significant_correlations.csv - Significant correlation pairs\n")
cat("9. GSE78097_7key_genes_spearman_correlation_heatmap.pdf - Spearman correlation heatmap\n")
cat("10. GSE78097_7key_genes_mantel_filtered.csv - Significant Mantel test results\n")
cat("11. GSE78097_7key_genes_mantel_network.pdf - Butterfly plot\n")
cat("\nGSE78097 analysis fully completed!\n")