# Set working directory
setwd("D:/dualdisease/IMMUNEreview")
# Load basic packages first
library(dplyr)
library(tidyr)
# Then load ggplot2
library(ggplot2)
# Load other visualization packages
library(pheatmap)
library(ComplexHeatmap)
library(corrplot)
# Finally load potentially conflicting packages
library(linkET)
library(IOBR)
library(ggsci)
library(ggpubr)
library(rstatix)
library(viridis)
library(scico)
library(Hmisc)

# Other necessary packages
library(e1071)
library(preprocessCore)
library(parallel)
library(usethis)
library(devtools)
library(bseqsc)
library(csSAM)
library(Rcpp)
library(CIBERSORT)
library(xCell)

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

# Import data
data <- read_table(
  file = "D:/dualdisease/WGCNA/SLEreal.data/GSE50772_GSE61635_ComBat_expr.txt",
  col_names = TRUE
)

# Convert to numeric matrix
data <- as.data.frame(data)
rownames(data) <- data$Gene
data <- data[, -1]
corrected_data <- as.data.frame(t(data))
data <- t(corrected_data)
colnames(data) <- make.unique(colnames(data))
data <- matrix(as.numeric(as.matrix(data)), 
               nrow = nrow(data),
               dimnames = list(rownames(data), colnames(data)))

# Normalize data
zero_var_genes <- apply(data, 1, function(x) var(x) == 0)
if (any(zero_var_genes)) {
  cat("Removing", sum(zero_var_genes), "zero variance genes:\n")
  data <- data[!zero_var_genes, ]
}

data_scaled <- t(scale(t(data)))
data <- data_scaled

# Save normalized data
write.csv(data, "Scaled_Gene_Expression.csv")

# Immune infiltration analysis
cibersort <- deconvo_tme(eset = data, method = "cibersort", arrays = FALSE, perm = 200)
epic <- deconvo_tme(eset = data, method = "epic", arrays = FALSE)
mcp <- deconvo_tme(eset = data, method = "mcpcounter")

# xCELL analysis
data_matrix <- as.matrix(data)
xcell_result <- xCellAnalysis(data_matrix)
xcell_df <- as.data.frame(t(xcell_result))
xcell_df$ID <- rownames(xcell_df)
rownames(xcell_df) <- NULL

# Other methods
estimate <- deconvo_tme(eset = data, method = "estimate")
timer <- deconvo_tme(eset = data, method = "timer", group_list = rep("stad",dim(data)[2]))
quantiseq <- deconvo_tme(eset = data, tumor = TRUE, arrays = FALSE, scale_mrna = TRUE, method = "quantiseq")
ips_res <- deconvo_tme(eset = data, method = "ips", plot= F)

# Merge results
tme_combine <- cibersort %>%
  inner_join(mcp, by = "ID") %>%
  inner_join(xcell_df, by = "ID") %>%
  inner_join(epic, by = "ID") %>%
  inner_join(estimate, by = "ID") %>%
  inner_join(timer, by = "ID") %>%
  inner_join(quantiseq, by = "ID") %>%
  inner_join(ips_res, by = "ID")

save(tme_combine, file = "tme_combine_scale.RData")

# Clear environment and reload data
rm(list = ls())
eset <- read.table("D:/dualdisease/WGCNA/Scaled_Gene_Expression.csv", 
                   header = TRUE, row.names = 1, sep = ",", check.names = FALSE,
                   stringsAsFactors = FALSE)

file_path <- "D:/dualdisease/SLE/SLE.combined/combined_phenotype.txt"
conditions <- read.delim(file_path, header = TRUE, stringsAsFactors = FALSE)

conditions$ID <- sub("^[^_]*_", "", conditions$SAMPLE_ID)
conditions$DISEASE <- ifelse(conditions$DISEASE %in% c("Control", "healthy"), "healthy", "disease")

# CIBERSORT analysis
eset <- eset[, !colnames(eset) %in% "symbol"]
res_cibersort <- deconvo_tme(eset, method = "cibersort", arrays = FALSE, perm = 100)

colnames(res_cibersort) <- gsub("_CIBERSORT", "", colnames(res_cibersort))
res_cibersort <- res_cibersort[,-c(24:26)]
names(res_cibersort) <- c("ID","B_cells_naive","B_cells_memory","Plasma_cells",
                          "T_cells_CD8","T_cells_CD4_naive","T_cells_CD4_memory_resting",
                          "T_cells_CD4_memory_activated","T_cells_follicular_helper",
                          "T_cells_regulatory_(Tregs)","T_cells_gamma_delta","NK_cells_resting",
                          "NK_cells_activated","Monocytes","Macrophages_M0","Macrophages_M1",
                          "Macrophages_M2","Dendritic_cells_resting","Dendritic_cells_activated",
                          "Mast_cells_resting","Mast_cells_activated","Eosinophils","Neutrophils")

res_cibersort$ID <- sub(".*_", "", res_cibersort$ID)
res_cibersort <- res_cibersort[-1, ]

# Sample matching
cat("Number of samples in res_cibersort:", nrow(res_cibersort), "\n")
cat("Number of samples in conditions:", nrow(conditions), "\n")

missing_in_cibersort <- setdiff(conditions$ID, res_cibersort$ID)
missing_in_conditions <- setdiff(res_cibersort$ID, conditions$ID)

cat("Samples missing in res_cibersort:", length(missing_in_cibersort), "\n")
if(length(missing_in_cibersort) > 0) print(missing_in_cibersort)

cat("Samples missing in conditions:", length(missing_in_conditions), "\n")
if(length(missing_in_conditions) > 0) print(missing_in_conditions)

# Use merge for exact matching
res_cibersort_with_group <- merge(res_cibersort, 
                                  conditions[, c("ID", "DISEASE")], 
                                  by = "ID", all.x = TRUE)

colnames(res_cibersort_with_group)[colnames(res_cibersort_with_group) == "DISEASE"] <- "Group"
res_cibersort_with_group$Group <- factor(res_cibersort_with_group$Group, levels = c("healthy", "disease"))

# Use merged data
res_cibersort <- res_cibersort_with_group

# Stacked bar chart
violin_dat <- res_cibersort %>% 
  pivot_longer(
    cols = -c("ID", "Group"),
    names_to = "ImmuneCell",
    values_to = "Score"
  )

mypalette <- c("#FF0000","#D8D8BF","#8E236B","#EAADEA","#BC8F8F","#5959AB","#0000FF","#2F4F4F","#3232CD","#FF00FF","#EAEAAE","#CC3299","#9370DB","#FF6EC7","#545454","#E47833","#856363","#E6E8FA","#3299CC","#9F5F9F","#8E2323","#007FFF","#B5A642","#DB7093","#FF1CAE","#D9D919","#CD7F32","#A68064","#A67D3D","#2F2F4F","#236B8E","#8C7853","#5F9F9F")

violin.cibersort <- violin_dat
p <- ggplot(violin.cibersort, aes(x=ID, y=Score, fill=ImmuneCell)) + 
  geom_bar(position = 'stack', stat = 'identity') + 
  scale_y_continuous(expand = c(0,0)) + 
  theme_bw() + 
  labs(x = '', y = 'Relative Percent', fill = '') + 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        legend.position = 'top') + 
  scale_fill_manual(values = mypalette) + 
  facet_grid(~Group, scales = "free", space = "free")

# Replace last line
ggsave("SLE_cibersort_stacked_barplot.pdf", p, width = 12, height = 8)
write.csv(res_cibersort, "SLE_immune_cell_cibersort.csv")

# Boxplot
b <- gather(res_cibersort, key = Cibersort, value = Score, -c(Group, ID))

box_plot <- ggboxplot(b, x = "Cibersort", y = "Score", 
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
# Use ggsave directly
ggsave("SLE_cibersort_boxplot.pdf", box_plot, width = 14, height = 8)

# Heatmap analysis - focus only on 7 key genes
load("tme_combine_scale.RData")
eset <- read.csv("Scaled_Gene_Expression.csv", row.names = 1, check.names = FALSE)

# Use 7 key genes
key_genes <- c("IFI6", "MX1", "NMI", "OAS3", "OASL", "SAMD9", "UBE2L6")
eset_genes <- eset[key_genes, ]

# Preprocessing
gene_sds <- apply(eset_genes, 1, sd, na.rm = TRUE)
valid_genes <- names(gene_sds)[gene_sds > 0]
eset_genes <- eset_genes[valid_genes, ]

cat("Number of valid key genes:", length(valid_genes), "\n")
cat("Invalid key genes:", setdiff(key_genes, valid_genes), "\n")

colnames(eset_genes) <- sub(".*_", "", colnames(eset_genes))
expr_t <- as.data.frame(t(eset_genes))
expr_t$ID <- rownames(expr_t)

merged_data <- merge(expr_t, tme_combine, by = "ID")

cat("Dimensions of merged data:", dim(merged_data), "\n")

full_cor_matrix <- cor(
  merged_data[, valid_genes],
  merged_data[, grep("_CIBERSORT|_MCPcounter|_EPIC|_estimate|_TIMER|_quantiseq|_IPS|_xCell", 
                     colnames(tme_combine))],
  use = "pairwise.complete.obs"
)

cat("Dimensions of correlation matrix:", dim(full_cor_matrix), "\n")

annotation_col <- data.frame(
  Method = sapply(colnames(full_cor_matrix), function(x) {
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
rownames(annotation_col) <- colnames(full_cor_matrix)

col_scale_mako <- viridis::mako(13, direction = -1)[1:10]
col_scale_acton <- scico::scico(10, direction = -1, palette = "acton")
my_color <- c(rev(col_scale_mako[1:8]), "white", col_scale_acton[1:8])

abs_max_val <- max(abs(full_cor_matrix), na.rm = TRUE)
ph_breaks <- seq(-abs_max_val, abs_max_val, length.out = length(my_color))

method_colors <- setNames(
  c("#66c2a5", "#fc8d62", "#a1c9f4", "#e78ac3", "#a6d854", "#ffd92f", "#e5c494", "#a999c9"),
  c("CIBERSORT", "MCPcounter", "EPIC", "ESTIMATE", "TIMER", "quanTIseq", "IPS", "xCell")
)

annotation_colors <- list(Method = method_colors)
mat_t <- t(full_cor_matrix)
annotation_row <- annotation_col

# Create and save heatmap
heatmap_obj <- pheatmap(
  mat = mat_t,
  color = my_color,
  breaks = ph_breaks,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  cellwidth = 15,
  cellheight = 8,
  annotation_row = annotation_row,
  annotation_colors = annotation_colors,
  fontsize_row = 10,
  fontsize_col = 10,
  angle_col = "45",
  na_col = "grey80",
  main = "7 Key Genes - Immune Infiltration Correlation"
)

pdf("7key_genes_immune_correlation_heatmap.pdf", width = 12, height = 12)
draw(heatmap_obj)
dev.off()

# Spearman correlation analysis - focus only on 7 key genes
df <- merged_data[, c(valid_genes, setdiff(colnames(tme_combine), "ID"))]
data_spearman <- rcorr(as.matrix(df), type = "spearman")
r_spearman <- data_spearman$r
p_spearman <- data_spearman$P

r_spearman_sub <- r_spearman[1:length(valid_genes), (length(valid_genes)+1):ncol(df)]
p_spearman_sub <- p_spearman[1:length(valid_genes), (length(valid_genes)+1):ncol(df)]

significant_pairs <- data.frame(
  Gene = rep(rownames(r_spearman_sub), each = ncol(r_spearman_sub)),
  Immune_Method = rep(colnames(r_spearman_sub), times = nrow(r_spearman_sub)),
  r = as.vector(r_spearman_sub),
  p = as.vector(p_spearman_sub)
) %>% filter(abs(r) > 0.6 & p < 0.05)

write.csv(significant_pairs, "7key_genes_significant_correlations.csv", row.names = FALSE)

# Check if there are significant correlation pairs
if(nrow(significant_pairs) > 0) {
  genes_keep <- unique(significant_pairs$Gene)
  immune_keep <- unique(significant_pairs$Immune_Method)
  
  r_spearman_filtered <- r_spearman_sub[genes_keep, immune_keep, drop = FALSE]
  p_spearman_filtered <- p_spearman_sub[genes_keep, immune_keep, drop = FALSE]
  
  my_col <- colorRampPalette(c("#9BCD9B","#B4EEB4", "gray98", "#F2B379","#DD5F60"))(100)
  
  # Save Spearman correlation heatmap as PDF
  pdf("7key_genes_spearman_correlation_heatmap.pdf", width = 10, height = 10)
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
           main = "7 Key Genes - Spearman Correlation (r > 0.6 & p < 0.05)")
  dev.off()
  cat("Spearman correlation heatmap PDF saved to: 7key_genes_spearman_correlation_heatmap.pdf\n")
} else {
  cat("Warning: No significant correlation pairs found (r > 0.6 & p < 0.05)\n")
}

# Mantel test and butterfly plot - focus only on 7 key genes
# Data preparation
riskscore <- merged_data[, valid_genes, drop = FALSE]
rownames(riskscore) <- merged_data$ID

data <- tme_combine
data <- data[data$ID %in% rownames(riskscore), ]
data <- as.data.frame(data)
rownames(data) <- data$ID
data <- data[, -which(colnames(data) == "ID"), drop = FALSE]
data <- data[rownames(riskscore), , drop = FALSE]

# Create spec_select list containing only 7 key genes
spec_select_list <- list()
for(i in seq_along(valid_genes)) {
  spec_select_list[[valid_genes[i]]] <- i
}

# Mantel test and butterfly plot - focus only on 7 key genes
# Data preparation
riskscore <- merged_data[, valid_genes, drop = FALSE]
rownames(riskscore) <- merged_data$ID

data <- tme_combine
data <- data[data$ID %in% rownames(riskscore), ]
data <- as.data.frame(data)
rownames(data) <- data$ID
data <- data[, -which(colnames(data) == "ID"), drop = FALSE]
data <- data[rownames(riskscore), , drop = FALSE]

# Create spec_select list containing only 7 key genes
spec_select_list <- list()
for(i in seq_along(valid_genes)) {
  spec_select_list[[valid_genes[i]]] <- i
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
write.csv(mantel_filtered, "7key_genes_mantel_filtered.csv", row.names = FALSE)

cat("Found", nrow(mantel_filtered), "significant associations\n")

# Create butterfly plot (if there are significant results)
if(nrow(mantel_filtered) > 0) {
  immune_cells <- unique(mantel_filtered$env)
  data_filtered <- data[, colnames(data) %in% immune_cells, drop = FALSE]
  
  # Calculate correlation matrix
  corr_matrix <- correlate(data_filtered)
  
  # Color settings
  mycolor1 <- c("#eb7e60", "#6ca3d4", "grey90")
  
  # Fix butterfly plot - add correct geometric layer to display correlation squares
  p4 <- qcorrplot(corr_matrix, 
                  type = "upper", 
                  diag = FALSE,
                  show.diag = FALSE) +
    # Add geometric layer to display correlation squares
    geom_tile() +  # Or use geom_square() if available
    geom_couple(
      aes(colour = pd, size = rd), 
      data = mantel_filtered,
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
    labs(title = "7 Key Genes - Immune Cell Correlation Network")
  
  # Save plot
  tryCatch({
    ggsave("7key_genes_mantel_network123.pdf", p4, width = 16, height = 12)
    cat("Butterfly plot saved to: 7key_genes_mantel_network123.pdf\n")
  }, error = function(e) {
    # Alternative saving method
    pdf("7key_genes_mantel_network.pdf", width = 16, height = 12)
    print(p4)
    dev.off()
    cat("Butterfly plot saved to: 7key_genes_mantel_network123.pdf (using alternative method)\n")
  })
  
} else {
  cat("No significant Mantel associations found, skipping butterfly plot generation\n")
}