library(WGCNA)
library(reshape2)
library(stringr)
library(data.table)
library(impute)
library(ggplot2)
library(ggpubr)
library(dplyr)

output_dir <- "D:/dualdisease/WGCNAreview"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Output directory created:", output_dir, "\n")
}
setwd(output_dir)
cat("Working directory set to:", getwd(), "\n")

load("D:/dualdisease/WGCNA/Preal.data/P_mexp_clin.RData")

gse13355_ids <- read.table("D:/dualdisease/T1D/P.1/GSE13355_expression.txt", 
                           header = TRUE, stringsAsFactors = FALSE)

gse13355_sample_ids <- colnames(gse13355_ids)[4:ncol(gse13355_ids)]
cat("Number of GSE13355 samples:", length(gse13355_sample_ids), "\n")

all_mexp_samples <- mexp$Sample
cat("Total samples in mexp:", length(all_mexp_samples), "\n")

matched_samples <- intersect(all_mexp_samples, gse13355_sample_ids)
cat("Number of matched GSE13355 samples:", length(matched_samples), "\n")

if (length(matched_samples) == 0) {
  stop("Error: No matching GSE13355 samples found!")
}

mexp_gse13355 <- mexp[mexp$Sample %in% matched_samples, ]
cat("GSE13355 expression matrix dimensions:", dim(mexp_gse13355), "\n")

dataExpr <- as.data.frame(mexp_gse13355[, -1])
rownames(dataExpr) <- mexp_gse13355$Sample

cat("Processed expression matrix dimensions:", dim(dataExpr), "\n")
cat("Number of samples:", nrow(dataExpr), "Number of genes:", ncol(dataExpr), "\n")

clin_clean <- clin
clin_clean$sampleid_clean <- gsub("GSE[0-9]+_", "", clin_clean$SAMPLE_ID)

clin_gse13355 <- clin_clean[clin_clean$sampleid_clean %in% rownames(dataExpr), ]
cat("Number of GSE13355 samples found in clin:", nrow(clin_gse13355), "\n")

missing_clin <- setdiff(rownames(dataExpr), clin_gse13355$sampleid_clean)
if (length(missing_clin) > 0) {
  cat("Warning: The following samples are missing clinical information:\n")
  print(missing_clin)
}

clin_processed <- data.frame(
  sampleid = rownames(dataExpr),
  stringsAsFactors = FALSE
)

clin_processed <- merge(clin_processed, 
                        clin_gse13355[, c("sampleid_clean", "DISEASE")], 
                        by.x = "sampleid", by.y = "sampleid_clean", 
                        all.x = TRUE)

clin_processed$PSO <- ifelse(clin_processed$DISEASE == "disease", 1, 0)
clin_processed$Control <- ifelse(clin_processed$DISEASE == "healthy", 1, 0)

rownames(clin_processed) <- clin_processed$sampleid
datTraits <- clin_processed[rownames(dataExpr), c("PSO", "Control")]

if (any(is.na(datTraits$PSO)) | any(is.na(datTraits$Control))) {
  cat("Warning: NA values present in clinical information\n")
  cat("Number of NA in PSO column:", sum(is.na(datTraits$PSO)), "\n")
  cat("Number of NA in Control column:", sum(is.na(datTraits$Control)), "\n")
}

cat("Final group distribution:\n")
cat("Disease samples (PSO=1):", sum(datTraits$PSO, na.rm = TRUE), "\n")
cat("Healthy controls (Control=1):", sum(datTraits$Control, na.rm = TRUE), "\n")
cat("Missing group information:", sum(is.na(datTraits$PSO)), "\n")

m.mad <- apply(dataExpr, 2, mad)
dataExprVar <- dataExpr[, which(m.mad > max(quantile(m.mad, probs = seq(0, 1, 0.25))[2], 0.01))]
cat("Number of genes retained after MAD filtering:", ncol(dataExprVar), "\n")

dataExpr_final <- dataExprVar
nGenes = ncol(dataExpr_final)
nSamples = nrow(dataExpr_final)
cat("Final expression matrix dimensions - Samples:", nSamples, "Genes:", nGenes, "\n")

cat("Performing data standardization...\n")
dataExpr_final <- as.data.frame(scale(dataExpr_final))
cat("Standardization completed. Data range:", 
    round(range(dataExpr_final, na.rm = TRUE), 3), "\n")

cat("Standardized data statistics:\n")
print(summary(as.vector(as.matrix(dataExpr_final))))

cat("Calculating immune scores...\n")
library(estimate)
library(tibble)

exp_matrix <- t(dataExpr_final)
write.table(exp_matrix, 
            file = "GSE13355_expression_for_estimate.txt", 
            sep = "\t", quote = FALSE, col.names = TRUE)

filterCommonGenes(input.f = "GSE13355_expression_for_estimate.txt", 
                  output.f = "GSE13355_common_genes.gct", 
                  id = "GeneSymbol")

estimateScore(input.ds = "GSE13355_common_genes.gct", 
              output.ds = "GSE13355_estimate_score.gct")

estimate_scores <- read.table("GSE13355_estimate_score.gct", 
                              skip = 2, header = TRUE)
result_clean <- estimate_scores[, -2]
result_transposed <- as.data.frame(t(result_clean[-1]))
colnames(result_transposed) <- result_clean[, 1]
result_transposed <- rownames_to_column(result_transposed, var = "sampleid")
colnames(result_transposed) <- c("sampleid", "StromalScore", "ImmuneScore", 
                                 "ESTIMATEScore", "TumorPurity")

clin_processed$immune_score <- result_transposed$ImmuneScore[match(clin_processed$sampleid, 
                                                                   result_transposed$sampleid)]

datTraits <- clin_processed[rownames(dataExpr), c("PSO", "Control", "immune_score")]

cat("Immune Score statistics:\n")
print(summary(datTraits$immune_score))

cat("Checking for batch effects...\n")
batch_info <- substr(rownames(dataExpr_final), 1, 6)
cat("Detected batch groups:", unique(batch_info), "\n")
cat("Number of batch groups:", length(unique(batch_info)), "\n")

if(length(unique(batch_info)) > 1) {
  cat("Multiple batches detected, performing PCA visualization...\n")
  
  pca_result <- prcomp(dataExpr_final, scale. = FALSE)
  pca_df <- as.data.frame(pca_result$x[, 1:2])
  pca_df$Batch <- batch_info
  pca_df$Group <- datTraits[rownames(dataExpr_final), "PSO"]
  pca_df$Group <- ifelse(pca_df$Group == 1, "Disease", "Control")
  
  pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch, shape = Group)) +
    geom_point(size = 3, alpha = 0.7) +
    theme_bw() +
    labs(title = "GSE13355 - PCA Batch Effect Check",
         subtitle = "Color by Batch, Shape by Disease Group") +
    stat_ellipse(aes(group = Batch), type = "norm", linetype = 2, alpha = 0.5)
  
  ggsave("GSE13355_batch_effect_check.pdf", pca_plot, width = 10, height = 8)
  cat("Batch effect check plot saved: GSE13355_batch_effect_check.pdf\n")
  
  cat("Recommendation: Multiple batches detected, consider using ComBat for batch effect correction\n")
  cat("Code to add: library(sva); corrected_data <- ComBat(dat = as.matrix(dataExpr_final), batch = batch_info)\n")
} else {
  cat("No significant batch effects detected\n")
}

cat("Performing data quality check...\n")
gsg <- goodSamplesGenes(dataExpr_final, verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) 
    printFlush(paste("Removing genes:", paste(names(dataExpr_final)[!gsg$goodGenes], collapse = ",")))
  if (sum(!gsg$goodSamples) > 0) 
    printFlush(paste("Removing samples:", paste(rownames(dataExpr_final)[!gsg$goodSamples], collapse = ",")))
  dataExpr_final = dataExpr_final[gsg$goodSamples, gsg$goodGenes]
}

nGenes = ncol(dataExpr_final)
nSamples = nrow(dataExpr_final)
cat("Data dimensions after quality check:", dim(dataExpr_final), "\n")

datTraits_final <- datTraits[rownames(dataExpr_final), ]

cat("Final number of samples for analysis:", nSamples, "\n")
cat("Number of disease samples:", sum(datTraits_final$PSO, na.rm = TRUE), "\n")
cat("Number of healthy controls:", sum(datTraits_final$Control, na.rm = TRUE), "\n")

sampleTree = hclust(dist(dataExpr_final), method = "average")

pdf("GSE13355_sample_clustering.pdf", width = 14, height = 7)
par(mfrow = c(1, 2))
plot(sampleTree, main = "GSE13355 - Sample Clustering", 
     sub = "", xlab = "", labels = FALSE)
abline(h = quantile(sampleTree$height, 0.95), col = "red")

if(length(unique(batch_info)) > 1) {
  plot(pca_plot)
}
dev.off()

cat("Generating Disease vs Healthy boxplot...\n")

boxplot_data <- datTraits_final
boxplot_data$Group <- ifelse(boxplot_data$PSO == 1, "Disease", "Healthy")
boxplot_data$Group <- factor(boxplot_data$Group, levels = c("Healthy", "Disease"))

cat("Immune Score statistics:\n")
print(summary(boxplot_data$immune_score))
cat("Group statistics:\n")
print(table(boxplot_data$Group))

p_boxplot <- ggplot(boxplot_data, aes(x = Group, y = immune_score, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("Healthy" = "#a3cb38", "Disease" = "#f79f1f")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  labs(
    x = "Group",
    y = "Immune Score",
    title = "GSE13355 - Immune Score Distribution",
    subtitle = "Disease vs Healthy Groups"
  ) +
  stat_compare_means(method = "t.test", 
                     label = "p.format",
                     label.x = 1.5, 
                     label.y = max(boxplot_data$immune_score, na.rm = TRUE) * 1.05)

ggsave("GSE13355_immune_score_boxplot.pdf", p_boxplot, width = 6, height = 6)
cat("Immune Score boxplot saved: GSE13355_immune_score_boxplot.pdf\n")

cat("\n=== Group Statistics ===\n")
cat("Healthy Group:", sum(boxplot_data$Group == "Healthy"), "samples\n")
cat("Disease Group:", sum(boxplot_data$Group == "Disease"), "samples\n")
cat("Immune Score comparison: t-test p-value displayed in plot\n")

# Outlier removal (optional) - kept as per original choice, not removed
# clust <- cutreeStatic(sampleTree, cutHeight = quantile(sampleTree$height, 0.95), minSize = 5)
# table(clust)

cat("Generating sample dendrogram and clinical traits heatmap...\n")

traitData <- datTraits_final[, c("PSO", "Control", "immune_score")]

traitData$immune_score <- scale(traitData$immune_score)

traitColors <- numbers2colors(traitData, signed = FALSE)

pdf("GSE13355_sample_dendrogram_trait_heatmap.pdf", width = 12, height = 8)
plotDendroAndColors(sampleTree, traitColors,
                    groupLabels = colnames(traitData),
                    main = "GSE13355 - Sample Dendrogram and Clinical Traits",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
cat("Sample dendrogram and clinical traits heatmap saved: GSE13355_sample_dendrogram_trait_heatmap.pdf\n")

cat("Clinical traits statistics:\n")
cat("PSO - Disease samples:", sum(traitData$PSO), "Healthy controls:", sum(traitData$Control), "\n")
cat("Immune Score - Mean:", round(mean(traitData$immune_score), 3), 
    "Standard deviation:", round(sd(traitData$immune_score), 3), "\n")

save(dataExpr_final, datTraits_final, nGenes, nSamples, 
     file = "GSE13355_WGCNA_processed.Rdata")

sink("GSE13355_preprocessing_summary.txt")
cat("=== GSE13355 Data Preprocessing Summary ===\n")
cat("Processing date:", date(), "\n")
cat("Final number of samples:", nrow(dataExpr_final), "\n")
cat("Final number of genes:", ncol(dataExpr_final), "\n")
cat("Disease samples:", sum(datTraits_final$PSO), "\n")
cat("Healthy controls:", sum(datTraits_final$Control), "\n")
cat("Data standardization: Completed (Z-score standardization)\n")
cat("Gene filtering: Completed (based on MAD)\n")
cat("Quality check: Completed\n")
cat("Number of batch groups:", length(unique(batch_info)), "\n")
cat("Batch groups:", paste(unique(batch_info), collapse = ", "), "\n")
cat("Standardized data range:", round(range(dataExpr_final, na.rm = TRUE), 3), "\n")
sink()

cat("\n=== GSE13355 Data Preprocessing Completed ===\n")
cat("Final number of samples:", nrow(dataExpr_final), "\n")
cat("Final number of genes:", ncol(dataExpr_final), "\n")
cat("Disease samples:", sum(datTraits_final$PSO), "\n")
cat("Healthy controls:", sum(datTraits_final$Control), "\n")
cat("Output file: GSE13355_WGCNA_processed.Rdata\n")
cat("Preprocessing summary: GSE13355_preprocessing_summary.txt\n")

enableWGCNAThreads()
options(stringsAsFactors = FALSE)

type = "unsigned"
corType = "pearson"
corFnc = ifelse(corType == "pearson", cor, bicor)
maxPOutliers = ifelse(corType == "pearson", 1, 0.05)
robustY = ifelse(corType == "pearson", TRUE, FALSE)

powers = c(c(1:20), seq(from = 12, to = 20, by = 2))

sft = pickSoftThreshold(dataExpr_final, powerVector = powers, 
                        networkType = type, verbose = 3)

pdf("GSE13355_soft_threshold.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
cex1 = 0.9

plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", cex.lab = 0.8,
     ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
     main = "GSE13355 - Scale Independence", cex.lab = 1.1, cex.main = 1.8)

text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.80, col = "red")

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = "GSE13355 - Mean Connectivity", cex.lab = 1.1, cex.main = 1.8)
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, 
     cex = cex1, col = "red")
dev.off()

power = sft$powerEstimate
if (is.na(power)) {
  r2_values <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
  power <- which(r2_values > 0.8)[1]
  if (is.na(power)) power <- 4
}
cat("Selected soft threshold power:", power, "\n")   ##7
library(rlang)

net = blockwiseModules(dataExpr_final,
                       power = power,
                       maxBlockSize = 5000,
                       TOMType = type,
                       minModuleSize = 30,
                       mergeCutHeight = 0.25,
                       numericLabels = TRUE,
                       pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "GSE13355_TOM",
                       verbose = 3)

moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
table(moduleColors)

save(net, file = "GSE13355_WGCNA_network.Rdata")

pdf("GSE13355_module_dendrogram.pdf", width = 12, height = 9)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "GSE13355 - Gene Dendrogram and Module Colors")
dev.off()

cat("Generating TOM heatmap...\n")

load(net$TOMFiles[1], verbose = TRUE)
TOM <- as.matrix(TOM)
dissTOM <- 1 - TOM

actual_nGenes <- nrow(dissTOM)
cat("TOM matrix dimensions:", dim(dissTOM), "\n")
cat("Current nGenes:", nGenes, "\n")

set.seed(123)
nSelect <- min(400, actual_nGenes)
cat("Number of genes selected for visualization:", nSelect, "\n")

select <- sample(actual_nGenes, size = nSelect)
selectTOM <- dissTOM[select, select]

selectColors <- moduleColors[select]

plotDiss <- selectTOM^7
diag(plotDiss) <- NA

pdf("GSE13355_TOM_heatmap.pdf", width = 9, height = 9)
TOMplot(plotDiss, 
        hclust(as.dist(selectTOM), method = "average"), 
        selectColors,
        main = "GSE13355 - TOM Heatmap",
        col = gplots::colorpanel(250, 'red', "orange", 'lemonchiffon'))
dev.off()
cat("TOM heatmap saved: GSE13355_TOM_heatmap.pdf\n")

cat("Performing module eigengene analysis...\n")

MEs = net$MEs
module_numbers <- as.numeric(gsub("ME", "", colnames(MEs)))
module_color_names <- labels2colors(module_numbers)

cat("Module numbers:", module_numbers, "\n")
cat("Corresponding colors:", module_color_names, "\n")

MEs_col <- MEs
colnames(MEs_col) <- paste0("ME", module_color_names)
MEs_col <- orderMEs(MEs_col)

cat("Renamed MEs column names:", colnames(MEs_col), "\n")

pdf("GSE13355_eigengene_network.pdf", width = 10, height = 8)
plotEigengeneNetworks(MEs_col, "GSE13355 - Eigengene Adjacency Heatmap", 
                      marDendro = c(3, 3, 2, 4), 
                      marHeatmap = c(3, 4, 2, 2), 
                      plotDendrograms = TRUE, xLabelsAngle = 90)
dev.off()
cat("Eigengene network plot saved\n")

cat("Performing module-trait relationship analysis...\n")

rownames(MEs_col) <- rownames(dataExpr_final)

moduleTraitCor <- cor(MEs_col, datTraits_final[, c("PSO", "immune_score")], 
                      use = "pairwise.complete.obs")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datTraits_final))

cat("moduleTraitCor dimensions:", dim(moduleTraitCor), "\n")
cat("moduleTraitCor row names:", rownames(moduleTraitCor), "\n")

format_text <- function(cor, p) {
  p_fmt <- ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
                  sprintf("%.3f", p))
  stars <- ifelse(p < 0.001, "***",
                  ifelse(p < 0.01, "**",
                         ifelse(p < 0.05, "*", "")))
  paste0(sprintf("%.2f", cor), "\n(", p_fmt, ")", stars)
}

textMatrix <- matrix(
  mapply(format_text, moduleTraitCor, moduleTraitPvalue),
  nrow = nrow(moduleTraitCor),
  ncol = ncol(moduleTraitCor)
)

pdf("GSE13355_module_trait_relationships.pdf", width = 8, height = 10)
par(mar = c(6, 8, 3, 3))
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = c("PSO", "Immune Score"),
  yLabels = rownames(moduleTraitCor),
  ySymbols = rownames(moduleTraitCor),
  colorLabels = FALSE,
  colors = colorRampPalette(c("blue", "white", "red"))(50),
  textMatrix = textMatrix,
  cex.text = 0.7,
  zlim = c(-1, 1),
  main = "GSE13355 - Module-Trait Relationships\n(PSO vs Immune Score)",
  textAdj = c(0.5, 0.5)
)
dev.off()
cat("Module-trait relationship heatmap saved\n")

cat("Filtering significantly correlated modules...\n")

pso_cor <- moduleTraitCor[, "PSO", drop = FALSE]
immune_cor <- moduleTraitCor[, "immune_score", drop = FALSE]

pso_filtered <- pso_cor[abs(pso_cor[, "PSO"]) > 0.2, , drop = FALSE]
cat("\nModules significantly correlated with PSO (|cor| > 0.2):\n")
print(pso_filtered)

immune_filtered <- immune_cor[abs(immune_cor[, "immune_score"]) > 0.2, , drop = FALSE]
cat("\nModules significantly correlated with Immune Score (|cor| > 0.2):\n")
print(immune_filtered)

combined_filtered <- moduleTraitCor[
  abs(moduleTraitCor[, "PSO"]) > 0.2 & 
    abs(moduleTraitCor[, "immune_score"]) > 0.2,
  , drop = FALSE
]
cat("\nModules meeting both PSO and Immune Score correlation thresholds (|cor| > 0.2):\n")
print(combined_filtered)

write.csv(pso_filtered, "GSE13355_PSO_filtered_modules.csv", quote = FALSE)
write.csv(immune_filtered, "GSE13355_ImmuneScore_filtered_modules.csv", quote = FALSE)
write.csv(combined_filtered, "GSE13355_Combined_filtered_modules.csv", quote = FALSE)

cat("Filtering results saved\n")

pso_cor_sorted <- sort(abs(pso_cor[, "PSO"]), decreasing = TRUE)
top_modules <- names(pso_cor_sorted)[1:min(3, length(pso_cor_sorted))]

cat("Selected key modules:\n")
print(top_modules)

for (module in top_modules) {
  module_color <- gsub("ME", "", module)
  cat("Extracting genes for module", module_color, "...\n")
  
  module_genes <- colnames(dataExpr_final)[moduleColors == module_color]
  
  if (length(module_genes) > 0) {
    write.table(module_genes,
                file = paste0("GSE13355_", module_color, "_module_genes.txt"),
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    cat("Module", module_color, "contains", length(module_genes), "genes\n")
  } else {
    cat("Warning: No genes found for module", module_color, "\n")
  }
}

target_modules <- c("brown", "turquoise", "purple")
pheno <- "PSO"

cat("Starting MM-GS analysis for key modules:\n")
print(target_modules)

common_samples <- intersect(rownames(dataExpr_final), rownames(MEs_col))
dataExpr_aligned <- dataExpr_final[common_samples, ]
MEs_aligned <- MEs_col[common_samples, ]
datTraits_aligned <- datTraits_final[common_samples, ]

cat("Calculating module membership...\n")
module_membership <- cor(dataExpr_aligned, MEs_aligned, use = "pairwise.complete.obs")
colnames(module_membership) <- paste0("MM_", colnames(MEs_aligned))

cat("Calculating gene significance...\n")
gene_significance <- cor(dataExpr_aligned, datTraits_aligned[, pheno, drop = FALSE], 
                         use = "pairwise.complete.obs")
colnames(gene_significance) <- paste0("GS_", pheno)

all_data <- list()
modNames <- gsub("^ME", "", colnames(MEs_aligned))
gene_names <- colnames(dataExpr_aligned)

cat("Available module names:", modNames, "\n")
cat("Target modules:", target_modules, "\n")

for (module in target_modules) {
  if (!(module %in% modNames)) {
    cat("Warning: Module", module, "not found in MEs\n")
    cat("Available modules:", modNames, "\n")
    next
  }
  
  gene_idx <- which(moduleColors == module)
  module_gene_names <- gene_names[gene_idx]
  
  cat("Processing module", module, ", number of genes:", length(module_gene_names), "\n")
  
  mm_column_name <- paste0("MM_ME", module)
  if (!mm_column_name %in% colnames(module_membership)) {
    cat("Warning: MM column", mm_column_name, "not found\n")
    cat("Available MM columns:", colnames(module_membership), "\n")
    next
  }
  
  mm <- module_membership[module_gene_names, mm_column_name]
  
  gs_column_name <- paste0("GS_", pheno)
  if (!gs_column_name %in% colnames(gene_significance)) {
    cat("Warning: GS column", gs_column_name, "not found\n")
    next
  }
  
  gs <- gene_significance[module_gene_names, gs_column_name]
  
  df <- data.frame(
    Gene = module_gene_names,
    MM = as.numeric(mm),
    GS = as.numeric(gs),
    Module = module,
    stringsAsFactors = FALSE
  )
  
  df <- df[!is.na(df$MM) & !is.na(df$GS), ]
  
  cat("Module", module, "valid gene count:", nrow(df), "\n")
  cat("MM range:", round(range(df$MM), 3), "\n")
  cat("GS range:", round(range(df$GS), 3), "\n")
  
  all_data[[module]] <- df
}

combined_data <- do.call(rbind, all_data)

if (nrow(combined_data) == 0) {
  stop("Error: No valid data available for MM-GS analysis")
}

cat("Total valid gene count:", nrow(combined_data), "\n")
cat("MM summary:\n")
print(summary(combined_data$MM))
cat("GS summary:\n")
print(summary(combined_data$GS))

combined_data$Significant <- ifelse(
  abs(combined_data$MM) > 0.8 & abs(combined_data$GS) > 0.2, 
  "Yes", "No"
)

significant_count <- sum(combined_data$Significant == "Yes")
cat("Number of significant genes (MM > 0.8 and |GS| > 0.2):", significant_count, "\n")

if (significant_count > 0) {
  cat("Number of significant genes per module:\n")
  print(table(combined_data$Module[combined_data$Significant == "Yes"]))
}

library(tidyr)
plot_data <- combined_data %>%
  mutate(
    MM_abs = abs(MM),
    GS_abs = abs(GS),
    color_group = ifelse(Significant == "Yes", as.character(Module), "grey")
  )

corr_stats <- plot_data %>%
  group_by(Module) %>%
  summarise(
    cor_val = cor(MM_abs, GS_abs, use = "complete.obs"),
    p_val = cor.test(MM_abs, GS_abs)$p.value,
    .groups = "drop"
  ) %>%
  complete(Module = target_modules, fill = list(cor_val = 0, p_val = 1)) %>%
  mutate(
    module_order = match(Module, target_modules),
    x_pos = 0.7,
    y_pos = 0.9 - (module_order - 1) * 0.15,
    label = ifelse(
      is.na(cor_val), 
      sprintf("%s: No data", Module),
      sprintf("%s: cor=%.2f, p=%.2e", Module, cor_val, p_val)
    )
  )

module_colors <- c(
  "brown" = "brown",
  "turquoise" = "turquoise", 
  "purple" = "purple",
  "grey" = "grey"
)

p <- ggplot(plot_data, aes(x = MM_abs, y = GS_abs, color = color_group)) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_vline(xintercept = 0.8, color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = 0.2, color = "black", linetype = "dashed", linewidth = 0.5) +
  annotate("text", x = 0.85, y = max(plot_data$GS_abs), 
           label = "MM = 0.8", hjust = 0, vjust = 1, color = "black", size = 3) +
  annotate("text", x = 0.1, y = 0.25, 
           label = "GS = 0.2", hjust = 0, vjust = 0, color = "black", size = 3) +
  geom_text(data = corr_stats, aes(x = x_pos, y = y_pos, label = label),
            hjust = 1, vjust = 1, size = 3, color = "black", fontface = "bold") +
  scale_color_manual(values = module_colors, name = "Module") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10, colour = "black"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  ) +
  labs(
    x = "Module Membership (|MM|)",
    y = paste0("Gene Significance (|GS|) for ", pheno),
    title = "GSE13355 - Module Membership vs. Gene Significance\n(Using Absolute Values)"
  )

ggsave("GSE13355_MM_GS_scatterplot_key_modules.pdf", plot = p, width = 10, height = 8, dpi = 300)
cat("MM-GS scatter plot saved: GSE13355_MM_GS_scatterplot_key_modules.pdf\n")

if (significant_count > 0) {
  significant_genes <- combined_data[combined_data$Significant == "Yes", ]
  write.csv(significant_genes, "GSE13355_significant_MM_GS_genes_key_modules.csv", row.names = FALSE)
  cat("Significant genes saved to: GSE13355_significant_MM_GS_genes_key_modules.csv\n")
  
  for (module in target_modules) {
    module_genes <- significant_genes[significant_genes$Module == module, "Gene"]
    if (length(module_genes) > 0) {
      write.table(module_genes,
                  file = paste0("GSE13355_", module, "_significant_genes.txt"),
                  row.names = FALSE, col.names = FALSE, quote = FALSE)
      cat("Significant genes for module", module, "saved. Count:", length(module_genes), "\n")
    }
  }
} else {
  cat("Warning: No significant genes found meeting the filtering criteria\n")
}

write.csv(combined_data, "GSE13355_MM_GS_key_modules_complete.csv", row.names = FALSE)
cat("Complete MM-GS results saved: GSE13355_MM_GS_key_modules_complete.csv\n")

cat("\n=== GSE13355 Key Modules MM-GS Analysis Completed ===\n")
cat("Analyzed modules:", paste(target_modules, collapse = ", "), "\n") ##brown,turquoise,purple
cat("Total gene count:", nrow(combined_data), "\n") #6855
cat("Significant gene count:", significant_count, "\n") #1767
cat("Output files:\n")
cat("- GSE13355_MM_GS_scatterplot_key_modules.pdf: MM-GS scatter plot\n")
cat("- GSE13355_significant_MM_GS_genes_key_modules.csv: Significant gene list\n")
cat("- GSE13355_MM_GS_key_modules_complete.csv: Complete MM-GS results\n")
cat("- GSE13355_*_significant_genes.txt: Module-specific significant genes\n")