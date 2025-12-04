# --------------------------
# 0. Environment Setup
# --------------------------
# Load required libraries with error handling
required_packages <- c('GEOquery', 'data.table', 'tidyr', 'tibble', 
                       'dplyr', 'limma', 'ggplot2', 'ggrepel', 
                       'pheatmap', 'plot3D', 'readr')

# Install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages) > 0) {
  install.packages(new_packages)
  # For Bioconductor packages
  if("GEOquery" %in% new_packages) {
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("GEOquery")
  }
}

# Load all libraries
invisible(lapply(required_packages, library, character.only = TRUE))

# Set working directory (modify to your path)
setwd("d:/dualdisease/T1D/P.1/")

# Set random seed for reproducibility
set.seed(1234)

# --------------------------
# 1. Data Loading & Preprocessing
# --------------------------
cat("Step 1: Loading and preprocessing GSE13355 data...\n")

# Load GEO series matrix file (no platform info)
gse_file <- 'd:/dualdisease/T1D/P.1/GEO.data3/GSE13355_series_matrix.txt'
if (file.exists(gse_file)) {
  gse <- getGEO(filename = gse_file, getGPL = FALSE)
  dat <- exprs(gse)  # Extract expression matrix
  cat("✓ Expression matrix loaded successfully (dimensions:", dim(dat), ")\n")
} else {
  stop(paste("Error: GSE file not found at", gse_file))
}

# Convert to numeric matrix and check for log2 transformation need
ex <- as.matrix(dat)
mode(ex) <- "numeric"

# Check if log2 transformation is needed using quantile criteria
qx <- quantile(ex, c(0.0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE)
LogC <- (qx[5] > 100) || 
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

# Perform log2 transformation if needed
if(LogC) {
  ex[ex <= 0] <- NA  # Replace non-positive values to avoid log errors
  dat <- log2(ex)
  cat("✓ log2 transformation completed\n")
} else {
  cat("✓ log2 transformation not needed (data already normalized)\n")
}

# Convert to data frame for further processing
expr <- as.data.frame(dat)

# --------------------------
# 2. Annotation Data Processing
# --------------------------
cat("Step 2: Processing gene annotation data...\n")

# Load annotation file (top.table from limma analysis)
anno_file <- "d:/dualdisease/T1D/P.1/GEO.data3/GSE13355.top.table.tsv"
if (file.exists(anno_file)) {
  gset <- read_tsv(anno_file, show_col_types = FALSE)
  # Select key annotation columns
  ann <- gset[, c("ID", "Gene.symbol", "logFC", "adj.P.Val")]
  cat("✓ Annotation file loaded successfully\n")
} else {
  stop(paste("Error: Annotation file not found at", anno_file))
}

# Process probe IDs for matching
ann$ID <- as.character(ann$ID)
expr$ID <- rownames(expr)
expr$ID <- sub("^X", "", expr$ID)  # Remove leading X from probe IDs

# Merge expression data with annotations (keep all records)
expr2 <- merge(ann, expr, by = "ID", all = TRUE)

# Clean gene symbols (remove // separators and whitespace)
expr2$Gene.symbol <- trimws(sapply(strsplit(as.character(expr2$Gene.symbol), "//"), `[`, 1))

# Filter out invalid gene symbols
expr3 <- expr2[!is.na(expr2$Gene.symbol) & 
                 expr2$Gene.symbol != "" & 
                 expr2$Gene.symbol != "---", ]
cat("✓ Gene symbol cleaning completed (remaining genes:", nrow(expr3), ")\n")

# --------------------------
# 3. Differential Expression Analysis
# --------------------------
cat("Step 3: Identifying differentially expressed genes (DEGs)...\n")

# Filter DEGs using standard thresholds
# |logFC| > 0.585 (~1.5-fold change) and adjusted p-value < 0.05
deg <- expr3[abs(expr3$logFC) > 0.585 & expr3$adj.P.Val < 0.05, ]

# Remove duplicate gene symbols (keep first occurrence)
deg <- deg[!duplicated(deg$Gene.symbol), ]
rownames(deg) <- deg$Gene.symbol

# Extract expression matrix (remove annotation columns)
deg_expr <- deg[, !colnames(deg) %in% c("ID", "Gene.symbol", "logFC", "adj.P.Val")]

# Transpose to sample × gene format
deg_expr <- t(deg_expr) %>% as.data.frame()
deg_expr$ID <- rownames(deg_expr)

cat("✓ DEG identification completed (total DEGs:", nrow(deg), ")\n")

# --------------------------
# 4. Clinical Data Processing
# --------------------------
cat("Step 4: Processing clinical metadata...\n")

# Extract phenotype data from GSE object
clin <- pData(gse)
clin$ID <- rownames(clin)

# Extract sample IDs and group information (title column contains PP/NN labels)
mc <- clin[, c("ID", "title")]
colnames(mc) <- c("ID", "group")

# Filter relevant samples (PP = Disease, NN = Healthy)
mc <- mc[grepl("PP|NN", mc$group), ]
mc$group <- ifelse(grepl("PP", mc$group), "Disease", "Healthy")

# Check group distribution
group_counts <- table(mc$group)
cat("✓ Clinical data processing completed\n")
cat("  - Disease samples:", group_counts["Disease"], "\n")
cat("  - Healthy samples:", group_counts["Healthy"], "\n")

# --------------------------
# 5. Merge Expression and Clinical Data
# --------------------------
cat("Step 5: Merging expression and clinical data...\n")

# Merge by sample ID
df <- merge(mc, deg_expr, by = "ID")
rownames(df) <- df$ID
df <- df[, -1]  # Remove redundant ID column

# Save merged dataset
output_file <- "GSE13355_expr&group.txt"
write.table(df, output_file, sep = "\t", quote = FALSE, row.names = TRUE)
cat("✓ Merged data saved to", output_file, "\n")

# --------------------------
# 6. PCA Analysis (3D Visualization)
# --------------------------
cat("Step 6: Performing PCA analysis...\n")

# Prepare group information for plotting
dfGroup <- data.frame(Group = df$group)
rownames(dfGroup) <- rownames(df)

# Perform PCA (scale variables to unit variance)
pca_result <- prcomp(df[, -1], scale. = TRUE)
pca_df <- as.data.frame(pca_result$x[, 1:3])  # Keep first 3 principal components

# Calculate variance explained by each PC
pVar <- round(pca_result$sdev^2 / sum(pca_result$sdev^2), 3)[1:3]
xName <- paste0("PC1 (", pVar[1]*100, "%)")
yName <- paste0("PC2 (", pVar[2]*100, "%)")
zName <- paste0("PC3 (", pVar[3]*100, "%)")

cat("  - PCA variance explained: PC1 =", pVar[1]*100, "%, PC2 =", pVar[2]*100, "%, PC3 =", pVar[3]*100, "%\n")

# Set visualization parameters
colors <- c("red", "blue")  # Disease = red, Healthy = blue
myColors <- colors[as.integer(as.factor(dfGroup$Group))]
my_pch <- 21:22  # Different shapes for groups
pchs <- my_pch[as.integer(as.factor(dfGroup$Group))]

# Generate 3D PCA plot
pdf("PCA3D_DEG_GSE13355.pdf", width = 7, height = 6)
scatter3D(
  x = pca_df$PC1, y = pca_df$PC2, z = pca_df$PC3,
  pch = pchs, cex = 2, col = "black", bg = myColors,
  xlab = xName, ylab = yName, zlab = zName,
  ticktype = "detailed", bty = "b2", box = TRUE,
  theta = 40, phi = 15, d = 2, colkey = FALSE,
  main = "3D PCA - GSE13355 (T1D vs Healthy)"
)
# Add legend
legend("bottom", legend = c("Disease (PP)", "Healthy (NN)"),
       pch = my_pch, pt.bg = colors, pt.cex = 2,
       bg = "white", bty = "n", horiz = TRUE, cex = 1.2)
dev.off()

cat("✓ 3D PCA plot saved as PCA3D_DEG_GSE13355.pdf\n")

# --------------------------
# 7. Volcano Plot Generation
# --------------------------
cat("Step 7: Generating volcano plot...\n")

# Prepare volcano plot data
vol_data <- expr3 %>%
  mutate(
    group = case_when(
      adj.P.Val < 0.05 & logFC > 0.585 ~ "up",
      adj.P.Val < 0.05 & logFC < -0.585 ~ "down",
      TRUE ~ "ns"
    ),
    log10FDR = -log10(adj.P.Val)
  ) %>%
  distinct(Gene.symbol, .keep_all = TRUE)  # Remove duplicate genes

# Count DEG categories
vol_counts <- table(vol_data$group)
cat("  - DEG counts: Up =", vol_counts["up"], ", Down =", vol_counts["down"], ", Not significant =", vol_counts["ns"], "\n")

# Identify top 10 significant genes (by absolute logFC)
top10 <- vol_data %>%
  filter(group != "ns") %>%
  arrange(desc(abs(logFC))) %>%
  slice_head(n = 10)

# Generate enhanced volcano plot
pdf("volcano_GSE13355.pdf", width = 10, height = 8)
volcano_plot <- ggplot(vol_data, aes(x = logFC, y = log10FDR, color = group)) +
  # Non-significant points
  geom_point(data = subset(vol_data, group == "ns"), 
             alpha = 0.4, size = 2, color = "gray70") +
  # Up-regulated points
  geom_point(data = subset(vol_data, group == "up"), 
             color = "#e74c3c", size = 3, alpha = 0.8) +
  # Down-regulated points
  geom_point(data = subset(vol_data, group == "down"), 
             color = "#3498db", size = 3, alpha = 0.8) +
  # Top 10 gene labels
  geom_text_repel(data = top10, aes(label = Gene.symbol), 
                  box.padding = 0.8, point.padding = 0.5, 
                  size = 4, fontface = "bold", color = "black") +
  # Threshold lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", 
             color = "black", linewidth = 1) +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", 
             color = "black", linewidth = 1) +
  # Custom color scale
  scale_color_manual(values = c("down" = "#3498db", "up" = "#e74c3c", "ns" = "gray70")) +
  # Labels and theme
  labs(x = "log2 Fold Change", 
       y = "-log10(Adjusted P-value)", 
       title = "Volcano Plot - GSE13355 (T1D vs Healthy)",
       color = "Expression") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  # Axis limits
  xlim(c(-max(abs(vol_data$logFC))*1.1, max(abs(vol_data$logFC))*1.1)) +
  ylim(c(0, max(vol_data$log10FDR)*1.1))

print(volcano_plot)
dev.off()

cat("✓ Volcano plot saved as volcano_GSE13355.pdf\n")

# --------------------------
# 8. Heatmap Analysis
# --------------------------
cat("Step 8: Generating heatmap of DEGs...\n")

# Prepare heatmap data
heat_data <- df[, -1]  # Remove group column
heat_data <- t(heat_data)  # Transpose to gene × sample

# Extract DEG expression values
DEGs <- rownames(deg)
heat_data <- heat_data[DEGs, ]

# Row-wise normalization (z-score)
heat_scaled <- t(apply(heat_data, 1, function(x) (x - mean(x)) / sd(x)))
colnames(heat_scaled) <- colnames(heat_data)

# Cap extreme values to reduce outlier impact
heat_scaled[heat_scaled > 2] <- 2
heat_scaled[heat_scaled < -2] <- -2

# Create sample annotation
sample_groups <- df$group
ann_df <- data.frame(Group = sample_groups)
rownames(ann_df) <- colnames(heat_scaled)

# Define custom color schemes
ann_colors <- list(Group = c(Disease = "#e74c3c", Healthy = "#3498db"))
heat_palette <- colorRampPalette(c("#3498db", "#ffffff", "#e74c3c"))(100)

# Generate enhanced heatmap
pdf("pheatmap_GSE13355.pdf", width = 10, height = 8)
pheatmap(
  heat_scaled,
  color = heat_palette,
  show_rownames = FALSE,
  show_colnames = FALSE,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  annotation_col = ann_df,
  annotation_colors = ann_colors,
  border_color = NA,
  fontsize = 12,
  main = "Differentially Expressed Genes - GSE13355 (T1D vs Healthy)",
  treeheight_row = 20,
  treeheight_col = 20,
  annotation_names_col = TRUE,
  annotation_legend = TRUE
)
dev.off()

cat("✓ Heatmap saved as pheatmap_GSE13355.pdf\n")

# --------------------------
# 9. Save Key Results
# --------------------------
cat("Step 9: Saving final results...\n")

# Save DEG list
write.table(deg, "GSE13355_DEGs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Save volcano plot data
write.table(vol_data, "GSE13355_volcano_data.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Save normalized heatmap data
write.table(heat_scaled, "GSE13355_heatmap_data.txt", sep = "\t", quote = FALSE)

cat("\n=============================================\n")
cat("Analysis completed successfully!\n")
cat("Output files:\n")
cat("  - GSE13355_expr&group.txt (merged expression and clinical data)\n")
cat("  - GSE13355_DEGs.txt (differentially expressed genes)\n")
cat("  - PCA3D_DEG_GSE13355.pdf (3D PCA visualization)\n")
cat("  - volcano_GSE13355.pdf (volcano plot)\n")
cat("  - pheatmap_GSE13355.pdf (DEG heatmap)\n")
cat("=============================================\n")