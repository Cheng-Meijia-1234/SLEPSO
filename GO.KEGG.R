library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# Prepare gene list
gl <- read.csv("D:/WEIXIN/WeChat Files/yichen548527/FileStorage/File/2025-07/Psoriasis_SLE_ALLOverlapGenes(1)(1).txt",header = F)

gene_list <- gl$V1

# Convert gene IDs from Symbol to Entrez ID (This step is crucial, KEGG typically requires Entrez ID)
gene_entrez <- bitr(gene_list, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Hs.eg.db)

# Handle potential ID conversion losses
cat("Total gene ID conversions:", nrow(gene_entrez),"/", length(gene_list),"\n")

# GO enrichment analysis - all genes
ego_all <- enrichGO(gene = gene_entrez$ENTREZID,
                    OrgDb = org.Hs.eg.db,
                    keyType ="ENTREZID",
                    ont ="ALL", # "BP", "MF", "CC" or "ALL"
                    pAdjustMethod ="BH",
                    pvalueCutoff =0.05,
                    qvalueCutoff =0.05)

# KEGG enrichment analysis - all genes
ekegg_all <- enrichKEGG(gene = gene_entrez$ENTREZID,
                        organism ="hsa", 
                        pvalueCutoff =0.05,
                        pAdjustMethod ="BH")

# Save enrichment results to files
if(!is.null(ego_all)) write.csv(ego_all@result,"GO_enrichment_all.csv")
if(!is.null(ekegg_all)) write.csv(ekegg_all@result,"KEGG_enrichment_all.csv")

# Enrichment result visualization - GO bar plot
if(!is.null(ego_all)) {
  pdf("GO_barplot_all.pdf", width =12, height =8)
  print(barplot(ego_all, showCategory =20))
  dev.off()
}

# Enrichment result visualization - KEGG bar plot
if(!is.null(ekegg_all)) {
  pdf("KEGG_barplot_all.pdf", width =12, height =8)
  print(barplot(ekegg_all, showCategory =20))
  dev.off()
}

# Display summary of enrichment results
cat("\nGO enrichment summary for all genes:\n")
if(!is.null(ego_all)) print(head(ego_all))

cat("\nKEGG enrichment summary for all genes:\n")
if(!is.null(ekegg_all)) print(head(ekegg_all))

# # Load required packages
# BiocManager::install("ggprism")
# if (!require("devtools", quietly = TRUE))
#  install.packages("devtools")
# devtools::install_github("dxsbiocc/gground")
library(gground)
library(ggprism)
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)

# Define color scheme - use similar but more harmonious colors as original code
pal <- c('#7bc4e2','#acd372','#fbb05b','#ed6ca4')
pal <- c('#eaa052','#b74147','#90ad5b','#23929c')
pal <- c('#c3e1e6','#f3dfb7','#dcc6dc','#96c38e')

# Use setReadable to convert gene IDs in GO enrichment results to gene symbols
ego_readable <- setReadable(ego_all, OrgDb = org.Hs.eg.db, keyType ="ENTREZID")
GO <- as.data.frame(ego_readable)

gl14 <- read.csv("D:/WEIXIN/WeChat Files/yichen548527/FileStorage/File/2025-07/14.txt",header = F)
my_genes <- gl14$V1

# Find rows in GO results that contain each key gene
matched_rows_GO <- lapply(my_genes, function(gene) 
  which(sapply(GO$geneID, function(x) gene %in% unlist(strsplit(x, "/"))))
)
names(matched_rows_GO) <- my_genes

# Organize into data.frame format for output
matched_df_GO <- data.frame(
  Gene = names(matched_rows_GO),
  Row = sapply(matched_rows_GO, function(x) if(length(x) > 0) paste(x, collapse = ",") else NA)
)
matched_df_GO

# Count how many key genes are present in each row's geneID
match_counts_GO <- sapply(GO$geneID, function(x) {
  genes_in_row <- unlist(strsplit(x, "/"))
  sum(genes_in_row %in% my_genes)
})

# Add to GO data frame (if you want to keep the original structure)
GO$MatchedGeneCount <- match_counts_GO
GO1 <- GO[order(GO$MatchedGeneCount,decreasing = T),]
GO1[, c("Description", "geneID", "MatchedGeneCount")]

# Note: The following line was in original code but KEGG was not defined yet, kept as is
# Add to KEGG data frame (if you want to keep the original structure)
# KEGG$MatchedGeneCount <- match_counts_KEGG

# Use setReadable to convert gene IDs in KEGG enrichment results to gene symbols
ekegg_readable <- setReadable(ekegg_all, OrgDb = org.Hs.eg.db, keyType ="ENTREZID")
KEGG <- as.data.frame(ekegg_readable)

# Find rows in KEGG results that contain each key gene
matched_rows_KEGG <- lapply(my_genes, function(gene) 
  which(sapply(KEGG$geneID, function(x) gene %in% unlist(strsplit(x, "/"))))
)
names(matched_rows_KEGG) <- my_genes

# Organize into data.frame format for output
matched_df_KEGG <- data.frame(
  Gene = names(matched_rows_KEGG),
  Row = sapply(matched_rows_KEGG, function(x) if(length(x) > 0) paste(x, collapse = ",") else NA)
)
matched_df_KEGG

# Count how many key genes are present in each row's geneID
match_counts_KEGG <- sapply(KEGG$geneID, function(x) {
  genes_in_row <- unlist(strsplit(x, "/"))
  sum(genes_in_row %in% my_genes)
})

# Add to KEGG data frame (if you want to keep the original structure)
KEGG$MatchedGeneCount <- match_counts_KEGG
KEGG1 <- KEGG[order(KEGG$MatchedGeneCount,decreasing = T),]
KEGG1[, c("Description", "geneID", "MatchedGeneCount")]

# Filter TOP gene pathways
# Select top5 for each GO category, sorted by p.adjust, for same p.adjust values select pathway with most enriched genes
use_pathway <- group_by(GO, ONTOLOGY) %>%
  top_n(10, wt = -p.adjust) %>%
  group_by(p.adjust) %>%
  top_n(1, wt = Count) %>%
  rbind(
    top_n(KEGG,10, -p.adjust) %>%
      group_by(p.adjust) %>%
      top_n(1, wt = Count) %>%
      mutate(ONTOLOGY ='KEGG')
  ) %>%
  ungroup() %>%
  mutate(ONTOLOGY = factor(ONTOLOGY,
                           levels = rev(c('BP','CC','MF','KEGG')))) %>%
  dplyr::arrange(ONTOLOGY, p.adjust) %>%
  mutate(Description = factor(Description, levels = Description)) %>%
  tibble::rowid_to_column('index')

# Construct left side marker data
# Width for left category labels and gene count dot plot
width <-0.5

# X-axis length
xaxis_max <- max(-log10(use_pathway$p.adjust)) +1

# Left category label data
rect.data <- group_by(use_pathway, ONTOLOGY) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  mutate(
    xmin = -3* width,
    xmax = -2* width,
    ymax = cumsum(n),
    ymin = lag(ymax, default =0) +0.6,
    ymax = ymax +0.4
  )

# Plot enrichment pathway graph
plot_enrichment <-function() {
  p <- use_pathway %>%
    ggplot(aes(-log10(p.adjust), y = index, fill = ONTOLOGY)) +
    geom_round_col(
      aes(y = Description), width =0.6, alpha =0.8
    ) +
    geom_text(
      aes(x =0.05, label = Description),
      hjust =0, size =5
    ) +
    geom_text(
      aes(x =0.1, label = geneID, colour = ONTOLOGY),
      hjust =0, vjust =2.6, size =3.5, fontface ='italic',
      show.legend =FALSE
    ) +
    # Gene count
    geom_point(
      aes(x = -width, size = Count),
      shape =21
    ) +
    geom_text(
      aes(x = -width, label = Count)
    ) +
    scale_size_continuous(name ='Count', range = c(5,12)) +
    # Category labels
    geom_round_rect(
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
          fill = ONTOLOGY),
      data = rect.data,
      radius = unit(2,'mm'),
      inherit.aes =FALSE
    ) +
    geom_text(
      aes(x = (xmin + xmax) /2, y = (ymin + ymax) /2, label = ONTOLOGY),
      data = rect.data,
      inherit.aes =FALSE
    ) +
    geom_segment(
      aes(x =0, y =0, xend = xaxis_max, yend =0),
      linewidth =1.5,
      inherit.aes =FALSE
    ) +
    labs(y =NULL) +
    scale_fill_manual(name ='Category', values = pal) +
    scale_colour_manual(values = pal) +
    scale_x_continuous(
      breaks = seq(0, xaxis_max,2),
      expand = expansion(c(0,0))
    ) +
    theme_prism() +
    theme(
      axis.text.y = element_blank(),
      axis.line = element_blank(),
      axis.ticks.y = element_blank(),
      legend.title = element_text()
    )
  return(p)
}

# Generate plot
enrichment_plot <- plot_enrichment()

# Display plot
print(enrichment_plot)

#### Simplified version
rect.data <- group_by(use_pathway, ONTOLOGY) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  mutate(
    xmin = -1.5* width, 
    xmax = -0.5* width, 
    ymax = cumsum(n),
    ymin = lag(ymax, default =0) +0.6,
    ymax = ymax +0.4
  )

use_pathway %>%
  ggplot(aes(-log10(p.adjust), y = index, fill = ONTOLOGY)) +
  geom_round_col(
    aes(y = Description), width =0.6, alpha =0.8
  ) +
  geom_text(
    aes(x =0.05, label = Description),
    hjust =0, size =5
  ) +
  # Gene count
  scale_size_continuous(name ='Count', range = c(5,16)) +
  # Category labels
  geom_round_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        fill = ONTOLOGY),
    data = rect.data,
    radius = unit(2,'mm'),
    inherit.aes =FALSE
  ) +
  geom_text(
    aes(x = (xmin + xmax) /2, y = (ymin + ymax) /2, label = ONTOLOGY),
    data = rect.data,
    inherit.aes =FALSE
  ) +
  geom_segment(
    aes(x =0, y =0, xend = xaxis_max, yend =0),
    linewidth =1.5,
    inherit.aes =FALSE
  ) +
  labs(y =NULL) +
  scale_fill_manual(name ='Category', values = pal) +
  scale_colour_manual(values = pal) +
  scale_x_continuous(
    breaks = seq(0, xaxis_max,2),
    expand = expansion(c(0,0))
  ) +
  theme_prism() +
  theme(
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text()
  )