options(stringsAsFactors = FALSE)
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

library(COSG)
library(harmony)
library(ggsci)
library(dplyr)
library(Seurat)
library(clustree)
library(cowplot)
library(data.table)
library(ggplot2)
library(patchwork)
library(stringr)
library(future)

plan(multisession, workers = 12)

setwd("path/to/project")

untar("GSE135779_RAW.tar", exdir = "GSE135779")

dir <- "GSE135779/"
fs <- list.files("./GSE135779/", pattern = "^GSM")
samples <- str_split(fs, "_", simplify = TRUE)[, 1]
samples
fs

lapply(unique(samples), function(x) {
  y <- fs[grepl(x, fs)]
  folder <- paste0("GSE135779/", str_split(y[1], "_", simplify = TRUE)[, 1])
  dir.create(folder, recursive = TRUE)
  file.rename(paste0("GSE135779/", y[1]), file.path(folder, "barcodes.tsv.gz"))
  file.rename(paste0("GSE135779/", y[3]), file.path(folder, "matrix.mtx.gz"))
})

tmp <- list.dirs("GSE135779")[-1]
tmp
names(tmp) <- list.dirs("GSE135779", recursive = FALSE)

ct <- Read10X(tmp)
ct

sce.all <- CreateSeuratObject(
  counts = ct,
  min.cells = 5,
  min.features = 300
)
sce.all

as.data.frame(sce.all@assays$RNA$counts[1:10, 1:2])
head(sce.all@meta.data, 10)
table(sce.all@meta.data$orig.ident)

sce.all$group <- ifelse(
  grepl("GSM4029907|4029914|4029915|4029922|4029923|4029930|4029931|4029936|4029937|4029938|4029939|4029940|4029942|4029947|4029948|4029949", sce.all$orig.ident),
  "HD",
  "SLE"
)

table(sce.all$group)
getwd()

dir.create("./1-QC")
setwd("./1-QC")
source("../GSE135779/qc.R")
sce.all.filt <- basic_qc(sce.all)
print(dim(sce.all))
print(dim(sce.all.filt))
setwd("../")
table(sce.all.filt$group)

sp <- "human"

dir.create("2-harmony")
getwd()
setwd("2-harmony")
source("../GSE135779/harmony.R")
sce.all.int <- run_harmony(sce.all.filt)
setwd("../")

homotypic.prop <- modelHomotypic(sc$seurat_clusters)
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
sc <- doubletFinder(
  sc,
  PCs = 1:35,
  pN = 0.25,
  pK = pK_bcmvn,
  nExp = nExp_poi.adj,
  reuse.pANN = FALSE,
  sct = FALSE
)

classification_col <- grep("^DF\\.classifications_", colnames(sc@meta.data), value = TRUE)
p <- DimPlot(sc, reduction = "umap", group.by = classification_col, raster = FALSE)
pdf("DoubletFinder_Result.pdf")
print(p)
dev.off()

table(Idents(sce.all.int))
table(sce.all.int$seurat_clusters)
table(sce.all.int$RNA_snn_res.0.1)
table(sce.all.int$RNA_snn_res.0.3)

p_umap <- DimPlot(
  sce.all.int,
  reduction = "umap",
  group.by = "group",
  label = TRUE,
  label.box = TRUE
)
p_umap
getwd()

dir.create("check-by-0.2")
setwd("check-by-0.2")
sel.clust <- "RNA_snn_res.0.3"
sce.all.int <- SetIdent(sce.all.int, value = sel.clust)
table(sce.all.int@active.ident)
source("../GSE135779/check-all-markers.R")
setwd("../")
getwd()

p_umap <- DimPlot(sce.all.int, reduction = "umap", label = TRUE, repel = TRUE)
p_umap

sce.all.int

celltype <- data.frame(
  ClusterID = 0:11,
  celltype = 0:11
)

celltype[celltype$ClusterID %in% c(0, 1, 3, 8), 2] <- "T"
celltype[celltype$ClusterID %in% c(4), 2] <- "NK"
celltype[celltype$ClusterID %in% c(5), 2] <- "B"
celltype[celltype$ClusterID %in% c(2, 6, 11), 2] <- "monocytes"
celltype[celltype$ClusterID %in% c(7), 2] <- "Platelets"
celltype[celltype$ClusterID %in% c(9), 2] <- "DC"
celltype[celltype$ClusterID %in% c(10, 11), 2] <- "Other"

dotplotmarker <- c(
  "CD3D", "CD3E", "CD3G", "CD2", "TRBC2",
  "MS4A1", "CD79A", "CD79B", "CD19",
  "KLRD1", "NKG7", "NCAM1",
  "IFITM3", "FCN1", "CD14", "VCAN", "FCGR3A",
  "TUBB1", "PF4", "PPBP",
  "IL3RA", "CLEC4C", "LILRB4"
)

head(celltype)
celltype
table(celltype$celltype)

sce.all.int@meta.data$celltype <- "NA"

for (i in 1:nrow(celltype)) {
  sce.all.int@meta.data[
    which(sce.all.int@meta.data$RNA_snn_res.0.3 == celltype$ClusterID[i]),
    "celltype"
  ] <- celltype$celltype[i]
}

Idents(sce.all.int) <- sce.all.int$celltype
table(Idents(sce.all.int))

sel.clust <- "celltype"
sce.all.int <- SetIdent(sce.all.int, value = sel.clust)
table(sce.all.int@active.ident)

p_umap <- DimPlot(
  sce.all.int,
  reduction = "umap",
  group.by = "celltype",
  label = TRUE,
  label.box = TRUE
)
p_umap

desired_order <- c("T", "B", "NK", "monocytes", "Platelets", "DC", "Other")
sce.all.int$celltype <- factor(sce.all.int$celltype, levels = desired_order)
Idents(sce.all.int) <- "celltype"

pdf("results/GSE135779_celltype_dotplot.pdf", width = 7.5, height = 8.8)
DotPlot(
  sce.all.int,
  features = dotplotmarker,
  assay = "RNA",
  cols = c("#FFE5B4", "#FF8C00")
) +
  coord_flip() +
  FontSize(y.text = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

sp <- "human"
dir.create("check-by-celltype")
setwd("check-by-celltype")
source("../GSE135779/check-all-markers.R")
setwd("../")
getwd()

p_umap <- DimPlot(
  sce.all.int,
  reduction = "umap",
  group.by = "celltype",
  label = TRUE,
  label.box = TRUE
)
p_umap

saveRDS(sce.all.int, "GSE135779/sce.all.int.rds")

library(plot1cell)
library(tidyverse)
library(RColorBrewer)

sce <- sce.all.int
Idents(sce) <- sce$celltype
colnames(sce@meta.data)
table(sce$group)

circ_data <- prepare_circlize_data(sce, scale = 0.8)
set.seed(1234)

mycolors <- c(
  "#E64A35",
  "#4DBBD4",
  "#01A187",
  "#FFC20A",
  "#3C5588",
  "#F29F80",
  "#7F2268"
)
cluster_colors <- mycolors
group_colors <- rand_color(length(names(table(sce$group))))
rep_colors <- rand_color(length(names(table(sce$orig.ident))))

pdf("GSE135779/circlize.pdf", height = 10, width = 10)
plot_circlize(
  circ_data,
  do.label = TRUE,
  pt.size = 0.15,
  col.use = cluster_colors,
  bg.color = FALSE,
  kde2d.n = 200,
  repel = FALSE,
  label.cex = 1
)

add_track(circ_data, group = "group", colors = group_colors, track_num = 2)
add_track(circ_data, group = "orig.ident", colors = rep_colors, track_num = 3)
dev.off()

cell_counts <- as.data.frame(table(sce@meta.data$celltype, sce$group))

ggplot(
  data = cell_counts,
  aes(x = forcats::fct_rev(Var2), y = Freq, fill = Var1)
) +
  geom_bar(position = "fill", stat = "identity") +
  coord_flip() +
  labs(x = NULL, y = "Cell Lineage Frequency") +
  scale_x_discrete(labels = c("HC", "Pso")) +
  scale_fill_manual(
    name = "Cell lineage",
    values = cluster_colors
  ) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(size = 16, color = "black"))

setwd("../")
dir.create("subT")
setwd("subT")

my_sub <- "T"
seu.obj <- sce.all.int
sub.cells <- subset(seu.obj, idents = my_sub)
f <- "obj.Rdata"

if (!file.exists(f)) {
  sub.cells <- sub.cells %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData(features = rownames(.)) %>%
    RunPCA(features = VariableFeatures(.)) %>%
    FindNeighbors(dims = 1:15) %>%
    FindClusters(resolution = 0.5) %>%
    RunUMAP(dims = 1:15)
  save(sub.cells, file = f)
}

load(f)
DimPlot(sub.cells, reduction = "umap", label = TRUE) + NoLegend()

library(stringr)

CD4_markers_list <- list(
  Tc = c("CD3D", "CD3E"),
  CD4 = c("CD4"),
  Treg = c("TNFRSF4", "BATF", "TNFRSF18", "FOXP3", "IL2RA", "IKZF2"),
  naive = c("CCR7", "SELL", "CD5"),
  Tfh = c("CXCR5", "BCL6", "ICA1", "TOX", "TOX2", "IL6ST"),
  ILC = c("TNFRSF25", "KRT81", "LST1", "AREG", "LTB", "CD69")
)

CD8_markers_list1 <- list(
  CD8 = c("CD8A", "CD8B"),
  TN_TCM = c("CCR7", "SELL", "TCF7", "LEF1"),
  TEM = c("GZMK"),
  TEFF = c("TBX21", "FCGR3A", "FGFBP2"),
  TRM = c("XCL1", "XCL2", "ITGAE", "CD69"),
  IEL_T = c("TMIGD2"),
  yT1c = c("GNLY", "PTGDS", "GZMB", "TRDC"),
  yT2c = c("TMN1", "HMGB2", "TYMS"),
  MAIT_T = c("SLC4A10")
)

CD8_markers_list2 <- list(
  CD8T = c("CD8A", "CD8B"),
  MAIT = c("ZBTB16", "NCR3", "RORA"),
  ExhaustedCD8T = c("LAG3", "TIGIT", "PDCD1", "HAVCR2", "CTLA4"),
  EffMemoryCD8 = c("EOMES", "ITM2C"),
  Resting_NK = c("XCL1", "XCL2", "KLRC1"),
  Cytotoxic_NK = c("CX3CR1", "FGFBP2", "FCGR3A", "KLRD1"),
  Pre_exhausted = c("IFNG", "PRF1", "GNLY", "GZMA", "NKG7", "GZMK")
)

cd4_and_cd8T_markers_list <- list(
  naive = c("CCR7", "SELL", "TCF7", "IL7R", "CD27", "CD28", "LEF1", "S1PR1"),
  CD8Trm = c("XCL1", "XCL2", "MYADM"),
  NKTc = c("GNLY", "GZMA"),
  Tfh = c("CXCR5", "BCL6", "ICA1", "TOX", "TOX2", "IL6ST"),
  th17 = c(
    "IL17A", "KLRB1", "CCL20", "ANKRD28",
    "IL23R", "RORC", "FURIN", "CCR6",
    "CAPG", "IL22"
  ),
  CD8Tem = c("CXCR4", "GZMH", "CD44", "GZMK"),
  Treg = c("FOXP3", "IL2RA", "TNFRSF18", "IKZF2"),
  naive2 = c("CCR7", "SELL", "TCF7", "IL7R", "CD27", "CD28"),
  CD8Trm2 = c("XCL1", "XCL2", "MYADM"),
  MAIT = c("KLRB1", "ZBTB16", "NCR3", "RORA"),
  yT1c = c("GNLY", "PTGDS", "GZMB", "TRDC"),
  yT2c = c("TMN1", "HMGB2", "TYMS"),
  yt = c("TRGV9", "TRDV2")
)

markers_list <- c(
  "CD4_markers_list",
  "CD8_markers_list1",
  "CD8_markers_list2",
  "cd4_and_cd8T_markers_list"
)

Idents(sce.all.int) <- "celltype"

.clean_feature_list <- function(feat_list, gene_universe) {
  feat_list <- lapply(feat_list, stringr::str_to_upper)
  allv <- unlist(feat_list, use.names = FALSE)
  dup <- names(table(allv))[table(allv) > 1]
  feat_list <- lapply(feat_list, function(v) setdiff(v, dup))
  feat_list <- lapply(feat_list, function(v) intersect(v, gene_universe))
  feat_list[lengths(feat_list) > 0]
}

invisible(lapply(markers_list, function(x) {
  genes_to_check <- get(x, inherits = TRUE)
  genes_to_check <- .clean_feature_list(genes_to_check, rownames(sub.cells))
  
  if (length(genes_to_check) == 0L) {
    message("skip ", x, ": no available genes")
    return(NULL)
  }
  
  p <- DotPlot(sub.cells, features = genes_to_check, assay = "RNA") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0("Marker check: ", x))
  
  w <- length(unique(unlist(genes_to_check))) / 5 + 6
  ggsave(
    filename = paste0("check_for_", x, ".pdf"),
    plot = p,
    width = w,
    height = 6
  )
}))

library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(DESeq2)
library(ggrepel)
library(ggsci)
library(RColorBrewer)
library(qs)
library(ggpubr)

sce <- sce.all.int
colnames(sce.all@meta.data)

dir.create("5-pseudobulks")
setwd("5-pseudobulks")

sce.all <- sce
av <- AggregateExpression(
  sce.all,
  group.by = c("orig.ident", "celltype"),
  assays = "RNA",
  return.seurat = FALSE
)
av <- as.data.frame(av[[1]])
head(av)[1:3, 1:3]

cg <- names(tail(sort(apply(log(av + 1), 1, sd)), 1000))
df <- cor(as.matrix(log(av[cg, ] + 1)))
colnames(df)

ac <- as.data.frame(str_split(colnames(df), "_", simplify = TRUE))
rownames(ac) <- colnames(df)
ac$V1

ac$group <- ifelse(
  grepl("GSM4029940|4029942|4029947|4029948|4029949", ac$V1),
  "aHD",
  "aSLE"
)

table(ac$group)
head(ac)

pheatmap::pheatmap(
  df,
  show_colnames = FALSE,
  show_rownames = FALSE,
  annotation_col = ac,
  filename = "cor_celltype-vs-orig.ident.pdf"
)

save(av, file = "av_for_pseudobulks.Rdata")

av <- AggregateExpression(
  sce.all,
  group.by = c("orig.ident", "celltype"),
  assays = "RNA"
)
av <- as.data.frame(av[[1]])
df <- log(av + 1)
head(ac)
celltp <- unique(ac$V2)
celltp

library(FactoMineR)
library(factoextra)
library(ggstatsplot)
library(pheatmap)

pca_list <- lapply(celltp, function(x) {
  exp <- df[, rownames(ac[ac$V2 == x, ])]
  cg <- names(tail(sort(apply(exp, 1, sd)), 1000))
  exp <- exp[cg, ]
  dat.pca <- PCA(as.data.frame(t(exp)), graph = FALSE)
  group_list <- ac[ac$V2 == x, "group"]
  this_title <- paste0(x, "_PCA")
  
  p2 <- fviz_pca_ind(
    dat.pca,
    geom.ind = "point",
    col.ind = group_list,
    palette = "Dark2",
    addEllipses = TRUE,
    legend.title = "Groups"
  ) +
    ggtitle(this_title) +
    theme_ggstatsplot() +
    theme(plot.title = element_text(size = 12, hjust = 0.5))
  
  p2
})

wrap_plots(pca_list, byrow = TRUE, nrow = 2)
ggsave("all_pca.pdf", width = 12, height = 6)
ggsave("all_pca.png", dpi = 300, width = 12, height = 6)

av <- AggregateExpression(
  sce.all,
  group.by = c("orig.ident", "group"),
  assays = "RNA",
  slot = "counts",
  return.seurat = FALSE
)
av <- as.data.frame(av[[1]])
head(av)[1:3, 1:3]

library(tibble)

counts_res <- av
head(counts_res)

colData <- data.frame(samples = colnames(counts_res))
colData <- colData %>%
  mutate(
    condition = ifelse(
      grepl("aSLE", samples),
      "aSLE",
      "aHD"
    )
  ) %>%
  column_to_rownames(var = "samples")

dds <- DESeqDataSetFromMatrix(
  countData = counts_res,
  colData = colData,
  design = ~condition
)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, name = "condition_Psoriasis_vs_Control")
res

gene <- c(
  "OASL","SCO2","IFI6","NMI",
  "SAMD9","MX1","UBE2L6","OAS3"
)

res[gene, ]
res["SAMD9", ]
res["UBE2L6", ]

resOrdered <- res[order(res$padj), ]
head(resOrdered)

DEG <- as.data.frame(resOrdered)
DEG_deseq2 <- na.omit(DEG)

DEG_deseq2 <- DEG_deseq2 %>%
  mutate(
    Type = if_else(
      padj > 0.05, "stable",
      if_else(
        abs(log2FoldChange) < 1, "stable",
        if_else(log2FoldChange >= 1, "up", "down")
      )
    )
  ) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  rownames_to_column("Symbol")

ggplot(DEG_deseq2, aes(log2FoldChange, -log10(padj))) +
  geom_point(
    size = 3.5,
    alpha = 0.8,
    aes(color = Type),
    show.legend = TRUE
  ) +
  scale_color_manual(values = c("#00468B", "gray", "#E64B35")) +
  ylim(0, 15) +
  xlim(-10, 10) +
  labs(x = "Log2(fold change)", y = "-log10(padj)") +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = 2,
    color = "black",
    lwd = 0.8
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = 2,
    color = "black",
    lwd = 0.8
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("volcano.pdf", width = 9, height = 7)

library(tidyverse)
library(ggstatsplot)

gene <- c(
  "RTP4", "SCO2", "IFI6", "LGALS3BP", "NMI", "IRF7",
  "SAMD9", "IFIT3", "MX1", "UBE2L6", "ISG15", "TAP1", "OAS3"
)

bdata <- as.data.frame(t(counts_res[gene, ]))
td <- bdata
class(td)

td$group <- "ss"
td[rownames(td) %like% "aSLE", ]$group <- "aSLE"
td[rownames(td) %like% "aHD", ]$group <- "aHD"
table(td$group)

datasets <- list(td)
dataset_names <- c("GSE135779")

library(tidyr)
library(stringr)

for (i in seq_along(datasets)) {
  df <- datasets[[i]]
  name <- dataset_names[i]
  
  df_long <- df %>%
    pivot_longer(
      cols = where(is.numeric),
      names_to = "gene",
      values_to = "Value"
    )
  
  df_long <- df_long %>%
    dplyr::filter(group %in% c("aSLE", "aHD")) %>%
    mutate(group = factor(group, levels = c("aHD", "aSLE")))
  
  genes <- unique(df_long$gene)
  
  outdir <- file.path(getwd(), paste0("plots_", name))
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  for (g in genes) {
    subdata <- df_long %>% dplyr::filter(gene == g)
    
    if (n_distinct(na.omit(subdata$group)) < 2 || nrow(na.omit(subdata)) < 3) next
    
    safe_g <- str_replace_all(g, "[/\\\\:*?\"<>| ]", "_")
    
    p <- try(
      ggbetweenstats(
        data = subdata,
        x = group,
        y = Value,
        title = paste0("Expression of ", g, " in ", name),
        messages = FALSE
      ),
      silent = TRUE
    )
    if (inherits(p, "try-error")) next
    
    ggsave(
      filename = file.path(outdir, paste0(name, "_", safe_g, ".pdf")),
      plot = p,
      width = 6,
      height = 5
    )
  }
}

library(gghalves)
library(cols4all)
library(rstatix)
library(tibble)

col_info <- data.frame(sample_full = colnames(df)) %>%
  mutate(
    gsm = str_extract(sample_full, "GSM\\d+"),
    cluster = sub(".*_", "", sample_full),
    Group = ifelse(
      grepl("GSM4029940|4029942|4029947|4029948|4029949", gsm),
      "aHD",
      "aSLE"
    )
  )

dt <- df %>%
  rownames_to_column("genes") %>%
  pivot_longer(-genes, names_to = "sample_full", values_to = "expressions") %>%
  left_join(col_info, by = "sample_full") %>%
  transmute(cluster, genes, expressions, Group)

head(dt)
tail(dt)

dt <- fread("5-pseudobulks/dt_cell.txt")

gene <- c(
  "OASL","SCO2","IFI6","NMI",
  "SAMD9","MX1","UBE2L6","OAS3"
)

dt2 <- dt[genes %in% gene, ]
dt2$genes <- factor(dt2$genes, levels = unique(dt2$genes))
dt <- dt2

p <- ggplot() +
  geom_half_violin(
    data = dt[dt$Group == "Control", ],
    aes(x = genes, y = expressions, fill = Group),
    color = "black",
    linewidth = 0.4,
    draw_quantiles = c(0.5),
    scale = "width"
  ) +
  facet_grid(rows = vars(cluster), scales = "free_y")
p

mytheme <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12, angle = 60, hjust = 1),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    strip.text.y = element_text(size = 12, angle = 0),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

fill_col <- c(Control = "#80B1D3", Psoriasis = "#FDB462")
point_col <- c(Control = "#BF5B17", Psoriasis = "#F0027F")

p <- ggplot() +
  geom_half_violin(
    data = dt %>% filter(Group == "Control"),
    aes(x = genes, y = expressions, fill = Group),
    color = "black",
    linewidth = 0.4,
    draw_quantiles = c(0.5),
    scale = "width",
    side = "l"
  ) +
  facet_grid(rows = vars(cluster), scales = "free_y")

p1 <- p +
  geom_half_violin(
    data = dt %>% filter(Group == "Psoriasis"),
    aes(x = genes, y = expressions, fill = Group),
    color = "black",
    linewidth = 0.4,
    draw_quantiles = c(0.5),
    scale = "width",
    side = "r"
  ) +
  facet_grid(rows = vars(cluster), scales = "free_y")

p2 <- p1 +
  mytheme +
  scale_fill_manual(values = fill_col, breaks = c("Control", "Psoriasis")) +
  scale_y_continuous(breaks = seq(0, 11, by = 5)) +
  labs(y = "Log Normalized Expression")

valid_panels <- dt %>%
  group_by(cluster, genes) %>%
  summarise(n_group = n_distinct(Group), .groups = "drop") %>%
  filter(n_group == 2)

dt_valid <- dt %>%
  inner_join(valid_panels %>% select(cluster, genes), by = c("cluster", "genes"))

pvals <- dt_valid %>%
  group_by(cluster, genes) %>%
  t_test(expressions ~ Group) %>%
  add_significance(p.col = "p") %>%
  ungroup()

pvals %>% select(cluster, genes, p, p.signif)

ypos <- dt_valid %>%
  group_by(cluster, genes) %>%
  summarise(
    y.position = max(expressions, na.rm = TRUE) * 1.05,
    .groups = "drop"
  )

p_anno <- pvals %>%
  left_join(ypos, by = c("cluster", "genes")) %>%
  mutate(
    xmin = genes,
    xmax = genes,
    label = p.signif
  )

p_final <- p2 +
  ggpubr::stat_pvalue_manual(
    p_anno,
    label = "label",
    xmin = "xmin",
    xmax = "xmax",
    y.position = "y.position",
    tip.length = 0,
    bracket.size = 0
  )

pdf("results/pea.pdf", width = 12, height = 8)
p_final
dev.off()
