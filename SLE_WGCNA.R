library(WGCNA)
library(reshape2)
library(stringr)
library(data.table)
library(impute)
library(ggplot2)
library(ggpubr)

setwd("D:/dualdisease/WGCNA/SLEreal.data/")
load("SLE_mexp_clin613.RData")

rownames(mexp) <- mexp$Sample
mexp <- mexp[,-1]
mexp[1:4,1:4]
dataExpr <- as.data.frame(t(mexp))
dataExpr[1:4,1:4]

dataExpr <- dataExpr[, clin$sampleid]

dataExpr[1:4, 1:4]
dataExpr <- as.data.frame(t(dataExpr))
dataExpr[1:4, 1:4]
rownames(dataExpr) <- gsub("GSE[0-9]+_", "", rownames(dataExpr))
dataExpr[1:4, 1:4]
head(clin$sampleid)
clin$sampleid <- gsub("GSE[0-9]+_", "", clin$sampleid)
all(rownames(dataExpr) == clin$sampleid)
dim(dataExpr)
dataExpr[1:4, 1:4]

dataExpr <- as.data.frame(t(dataExpr))
dataExpr <- as.data.frame(scale(dataExpr, center = TRUE, scale = TRUE))

m.mad <- apply(dataExpr,1,mad)

dataExprVar <- dataExpr[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]

dataExpr=as.data.frame(t(dataExprVar));
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
dim(dataExpr)

gsg <- goodSamplesGenes(dataExpr,verbose = 3)
gsg$allOK
if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
dim(dataExpr)

sampleTree = hclust(dist(dataExpr), method = "average")
summary(sampleTree$height)
pdf("SLE_WGCNA_mt.pdf",width = 14,height = 7)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",labels = F)
abline(h=100,col="red")
dev.off()

clust <- cutreeStatic(sampleTree, cutHeight = 100)
table(clust)

outsamples <- (clust==2)
outdataExpr <- dataExpr[outsamples,]
outdataExpr[,1:4]

keepsamples <- (clust==1)
dataExpr <- dataExpr[keepsamples,]
dim(dataExpr)

load("D:/dualdisease/WGCNA/SLEreal.data/SLE_mexp_clin613.RData")

exists("clin")
names(clin)
length(clin$sampleid)
length(clin$group)
extracted_data <- data.frame(
  V1 = clin$sampleid,
  sampleid = clin$sampleid,
  group = clin$group
)
write.csv(
  extracted_data,
  file = "D:/dualdisease/WGCNA/SLEreal.data/extracted_data.csv",
  row.names = FALSE
)

os_c <- fread("D:/dualdisease/WGCNA/SLEreal.data/extracted_data.csv",header = T)

colnames(os_c)

os_c$SLE <- ifelse(os_c$group == "SLE", 1, 0)
os_c[1:10,]

os_c$Control <- ifelse(os_c$group == "Control", 1, 0)
os_c[1:10,]

head(os_c)

os_c <- os_c[, -c(1, 3)]
os_c$sampleid <- sapply(strsplit(os_c$sampleid, "_"), "[", 2)

head(os_c)

library(estimate)
library(tidyverse)
rownames(mexp)
dim(mexp)
mexp$Sample
mexp[1:4,1:4]
rownames(mexp) <- mexp$Sample
exp <- mexp[, -1]
exp[1:4,1:4]
dim(exp)
colnames(exp) <- gsub("^[^_]+_", "", colnames(exp))
texp <- t(exp)
write.table(texp, 
            file = "D:/dualdisease/WGCNA/SLEreal.data/exp_cleaned613_WYC.txt", 
            sep = "\t", 
            quote = FALSE, 
            row.names = TRUE)

exp <- read.table("D:/dualdisease/WGCNA/SLEreal.data/exp_cleaned613_WYC.txt", 
                  sep = '\t', 
                  row.names = 1,
                  check.names = FALSE, 
                  stringsAsFactors = FALSE, 
                  header = TRUE)

head(exp)

nrow(exp)
exp_wyc <- exp[,-c(211:212)]

expression_data <- as.matrix(exp_wyc)
expression_data[1:4,1:4]
write.table(expression_data, 
            file ="wyc_expression_matrix.tsv", 
            sep ="\t", quote =FALSE, 
            col.names =T)
filterCommonGenes(input.f ="wyc_expression_matrix.tsv", 
                  output.f ="common_genes.gct", id ="GeneSymbol")

estimateScore(input.ds ="common_genes.gct", 
              output.ds ="estimate_score.gct")
estimate_scores <- read.table("estimate_score.gct", 
                              skip =2, header =TRUE)

head(estimate_scores)
write.csv(estimate_scores, "exp_estimate_score_clean613.csv", quote = FALSE)
result_clean <- read.csv("exp_estimate_score_clean613.csv", row.names = 1)

head(result_clean)

colnames(result_clean)
result_clean <- result_clean[, -2]
result_clean <- result_clean[, -1]

result_transposed <- as.data.frame(t(result_clean))
head(result_transposed)
library(tibble)

result_transposed <- rownames_to_column(result_transposed, var = "ID")

colnames(result_transposed) <- c("ID", "StromalScore", "ImmuneScore", "ESTIMATEScore", "TumorPurity")

head(result_transposed)
output_dir <- "D:/dualdisease/WGCNA/SLEreal.data/"
output_file <- file.path(output_dir, "result_transposed613.rdata")
save(result_transposed, file = output_file)

output_dir <- "D:/dualdisease/WGCNA/SLEreal.data/"
output_file <- file.path(output_dir, "result_transposed613.txt")

write.table(result_transposed, 
            file = output_file,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

load("D:/dualdisease/WGCNA/SLEreal.data/result_transposed613.rdata")

str(result_transposed)
summary(result_transposed$ImmuneScore)
os_c$immune_score <- result_transposed$ImmuneScore
os_c[1:4,]
colnames(os_c)
view(os_c)
cl <- os_c
colnames(cl)

cl_unique <- cl[!duplicated(cl), ]
dim(cl_unique)
cl <- cl_unique
names(cl)[1] <- "id"

bracsam <- rownames(dataExpr)
traitRows <- intersect(bracsam, cl$id)

dataExpr <- dataExpr[traitRows,]
datTraits <- cl[id %in% traitRows,-1]
datTraits <- as.data.frame(datTraits)
rownames(datTraits) <- cl[id %in% traitRows,id]
datTraits[1:20,]
os_c[1:20,]
nrow(os_c)

im_nc <- datTraits
im_nc$SLE <- as.factor(im_nc$SLE)
im_nc$immune_score <- scale(im_nc$immune_score)

pdf("SLE_boxplot_mt613.pdf",width = 14,height = 7)
p <-ggplot(im_nc,aes(x=SLE,y=immune_score))+
  geom_boxplot(aes(fill=SLE),alpha=0.7)+
  geom_jitter(aes(color =SLE))+
  scale_fill_manual(values =c("#f79f1f","#a3cb38"))+
  scale_color_manual(values =c("#f79f1f","#a3cb38"))+
  theme_bw()+
  theme(panel.grid =element_blank())
p <- p+stat_compare_means(method ="wilcox.test",
                          label ="p.format",
                          label.x= 2,
                          label.y = 4,
                          show.legend =F
)
p
dev.off()

sampleTree2 = hclust(dist(dataExpr), method = "average")
traitColors <- numbers2colors(datTraits, signed = F)
pdf("SLE_mt_WGCNA_cluster_clinical613.pdf",width = 7,height = 4.5)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and clincal traits heatmap",
                    sub="", xlab="", dendroLabels = F,cex.lab = 1.6,cex.main = 1.8)
dev.off()

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
cl <- datTraits
save(dataExpr,cl,nGenes,nSamples,file = "SLE_mt_WGCNA613.Rdata")

rm(list = ls())
gc()

load("D:/dualdisease/WGCNA/SLEreal.data/613/SLE_mt_WGCNA613.Rdata")

options(stringsAsFactors = FALSE)
library(WGCNA)
allowWGCNAThreads()
type = "unsigned"
exprMat <- "SLE_mt_WGCNA613.txt"
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
maxPOutliers = ifelse(corType=="pearson",1,0.05)
robustY = ifelse(corType=="pearson",T,F)

powers = c(c(1:20), seq(from = 12, to=20, by=2))
library(WGCNA)
sft = pickSoftThreshold(dataExpr, powerVector=powers, 
                        networkType=type, verbose=3)

par(mfrow = c(1,2))
cex1 = 0.9
pdf("SLE_MT_WGCNA_scale_mean613.pdf", width = 7,height = 7)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",cex.lab = 0.8,
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"),cex.lab = 1.1,cex.main = 1.8)  

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red") 
abline(h=0.80,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"),cex.lab = 1.1,cex.main = 1.8)
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")
while (!is.null(dev.list()))  dev.off()
power = sft$powerEstimate
power
dr <- as.data.frame(dataExpr)

corType = "pearson"
net = blockwiseModules(dataExpr,
                       power = 4,
                       maxBlockSize = nGenes,
                       TOMType = type,
                       minModuleSize = 15,
                       mergeCutHeight = 0.30,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType, 
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       saveTOMFileBase =  paste0(exprMat,".tom"),
                       verbose = 3)
table(net$colors)
save(net,file = "D:/dualdisease/WGCNA/SLEreal.data/SLE_wgcna_net613.Rdata")
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)

pdf("CD_MT_WGCNA_CLUSTER_DEND613.pdf", width = 7, height = 4.5)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE ,guideHang = 0.05,setLayout = T,cex.main = 2,cex.lab = 2)
dev.off()

library(stringr)
library(WGCNA)
library(data.table)

MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
setwd("D:/dualdisease/WGCNA/SLEreal.data")
pdf("SLE_MT_WGCNA_Eigengene_text613.pdf",width = 7,height = 7)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4), 
                      marHeatmap = c(3,4,2,2), 
                      plotDendrograms = T, xLabelsAngle = 90)
dev.off()

load(net$TOMFiles[1], verbose=T)
TOM <- as.matrix(TOM)
dissTOM = 1-TOM
nSelect = 400   

set.seed(10)
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]

pdf(file="MT_TOM_plot613.pdf")
install.packages("gplots")
library(gplots)

plotDiss = selectTOM^7
diag(plotDiss) = NA
TOMplot(plotDiss, selectTree, selectColors, 
        main = "Network heatmap plot, selected genes",
        col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'))
dev.off()

library(WGCNA)

load("D:/dualdisease/WGCNA/SLEreal.data/613/SLE_mt_WGCNA613.Rdata")
load("D:/dualdisease/WGCNA/SLEreal.data/613/SLE_mexp_clin613.RData")
load("D:/dualdisease/WGCNA/SLEreal.data/613/SLE_wgcna_net613.Rdata")

MEs <- net$MEs
datTraits <- cl[, c("SLE", "immune_score")]
rownames(MEs) <- rownames(dataExpr)

module_labels <- net$colors
color_labels <- labels2colors(module_labels)
unique_modules <- unique(module_labels[module_labels != 0])
module_color_map <- data.frame(
  module_label = unique_modules,
  color = labels2colors(unique_modules)
)
print("model_and_color")
print(module_color_map)

color_mapping <- setNames(
  object = paste0("ME", module_color_map$color),
  nm = paste0("ME", module_color_map$module_label)
)

colnames(MEs) <- sapply(colnames(MEs), function(x) {
  ifelse(x %in% names(color_mapping), color_mapping[x], x)
})

target_order <- c(
  "MEblack", "MEpink", "MEgreenyellow", "MEturquoise",
  "MEbrown", "MEyellow", "MEpurple", "MEblue",
  "MEmagenta", "MEgreen", "MEred", "MEgrey"
)
valid_order <- target_order[target_order %in% colnames(MEs)]
MEs <- MEs[, valid_order, drop = FALSE]

moduleTraitCor <- cor(MEs, datTraits, use = "pairwise.complete.obs")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datTraits))

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

setwd("D:/dualdisease/WGCNA/SLEreal.data")
pdf("Module_Trait_Relationships613.pdf", width = 10, height = 6)
par(mar = c(8, 10, 3, 3))
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(datTraits),
  yLabels = colnames(MEs),
  ySymbols = colnames(MEs),
  colorLabels = FALSE,
  colors = colorRampPalette(c("blue", "white", "red"))(50),
  textMatrix = textMatrix,
  cex.text = 0.7,
  zlim = c(-1, 1),
  main = "Module-Trait Relationships (SLE vs Immune Score)",
  textAdj = c(0.5, 0.5)
)
dev.off()

write.csv(moduleTraitCor, "Module_Trait_Correlation613.csv", quote = FALSE)
write.csv(moduleTraitPvalue, "Module_Trait_Pvalue613.csv", quote = FALSE)

cor_matrix <- read.csv("Module_Trait_Correlation613.csv", row.names = 1, check.names = FALSE)

pvalue_matrix <- read.csv("Module_Trait_Pvalue613.csv", row.names = 1, check.names = FALSE)

sle_cor <- moduleTraitCor[, "SLE", drop = FALSE]
immune_cor <- moduleTraitCor[, "immune_score", drop = FALSE]

sle_filtered <- sle_cor[abs(sle_cor[, "SLE"]) > 0.2, , drop = FALSE]
cat("\nSLE > 0.2 model：\n")
print(sle_filtered)

immune_filtered <- immune_cor[abs(immune_cor[, "immune_score"]) > 0.2, , drop = FALSE]
cat("\nImmune Score > 0.2 model：\n")
print(immune_filtered)

combined_filtered <- moduleTraitCor[
  abs(moduleTraitCor[, "SLE"]) > 0.2 & 
    abs(moduleTraitCor[, "immune_score"]) > 0.2,
  , drop = FALSE
]
cat("\n同时满足 SLE|>0.2| 和 Immune Score|>0.2| 的模块：\n")
print(combined_filtered)

write.csv(sle_filtered, "SLE_filtered_modules613.csv", quote = FALSE)
write.csv(immune_filtered, "ImmuneScore_filtered_modules613.csv", quote = FALSE)
write.csv(combined_filtered, "Combined_filtered_modules613.csv", quote = FALSE)

print(module_color_map)

modules_to_extract <- list(
  MEblue = 2,
  MEmagenta = 9
)

gene_names <- colnames(dataExpr)

for (module_name in names(modules_to_extract)) {
  module_label <- modules_to_extract[[module_name]]
  module_genes <- gene_names[net$colors == module_label]
  write.table(module_genes, 
              file = paste0(module_name, "_genes613.txt"),
              row.names = FALSE, 
              col.names = FALSE, 
              quote = FALSE)
}

blue_genes <- sum(net$colors == modules_to_extract$MEblue)
magenta_genes <- sum(net$colors == modules_to_extract$MEmagenta)

total <-  blue_genes+magenta_genes
total

cat(sprintf("MEblue (2): %d \nMEmagenta (9): %d ",
            blue_genes, 
            magenta_genes, 
            total))

setwd("D:/dualdisease/WGCNA/SLEreal.data")

if (ncol(dataExpr) < nrow(dataExpr)) {
  dataExpr <- t(dataExpr)
}

common_samples <- intersect(rownames(dataExpr), rownames(MEs))
dataExpr <- dataExpr[common_samples, ]
MEs <- MEs[common_samples, ]
datTraits <- datTraits[common_samples, ]

module_membership <- cor(dataExpr, MEs, use = "pairwise.complete.obs")
colnames(module_membership) <- paste0("MM_", colnames(MEs))
rownames(module_membership) <- colnames(dataExpr)
module_membership_pvalue <- corPvalueStudent(module_membership, nSamples = nrow(dataExpr))
rownames(module_membership_pvalue) <- colnames(dataExpr)

gene_significance <- cor(dataExpr, datTraits, use = "pairwise.complete.obs")
colnames(gene_significance) <- paste0("GS_", colnames(datTraits))
rownames(gene_significance) <- colnames(dataExpr)
gene_significance_pvalue <- corPvalueStudent(gene_significance, nSamples = nrow(dataExpr))
rownames(gene_significance_pvalue) <- colnames(dataExpr)

gene_module_colors <- labels2colors(net$colors)
stopifnot(length(gene_module_colors) == ncol(dataExpr))

results_df <- data.frame(
  Gene = colnames(dataExpr),
  Module = gene_module_colors,
  module_membership,
  gene_significance,
  module_membership_pvalue,
  gene_significance_pvalue,
  check.names = FALSE
)
write.csv(results_df, "MM_GS_Results613.csv", row.names = FALSE)

target_modules <- c("blue", "magenta")
pheno <- "immune_score"
modNames <- gsub("^ME", "", colnames(MEs))
moduleColors <- labels2colors(net$colors)
gene_names <- colnames(dataExpr)
stopifnot(length(moduleColors) == length(gene_names))

all_data <- list()

for (module in target_modules) {
  if (!(module %in% modNames)) next
  gene_idx <- which(moduleColors == module)
  module_gene_names <- gene_names[gene_idx]
  module_index <- which(modNames == module)
  mm <- abs(module_membership[module_gene_names, module_index])
  gs_vec <- cor(dataExpr[, module_gene_names, drop = FALSE], datTraits[, pheno], use = "pairwise.complete.obs")
  gs <- abs(gs_vec)
  df <- data.frame(
    Gene = module_gene_names,
    MM = mm,
    GS = gs,
    Module = module
  )
  all_data[[module]] <- df
}

combined_data <- do.call(rbind, all_data)
combined_data$Significant <- ifelse(combined_data$MM > 0.8 & combined_data$GS > 0.2, "Yes", "No")
significant_genes <- subset(combined_data, Significant == "Yes")

write.csv(significant_genes, "Significant_Genes_MM0_GS06132.csv", row.names = FALSE)
gene_list <- significant_genes$Gene
writeLines(gene_list, "Key_Genes_List6132.txt")