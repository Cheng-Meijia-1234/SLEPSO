rm(list=ls())
library('GEOquery')
library('data.table')
library('tidyr')
library('tibble')
library('dplyr')
setwd("d:/dualdisease/SLE/SLE.2/")

# Retrieve data from GEO database using specified filename without obtaining platform information
gse <- getGEO(filename='d:/dualdisease/SLE/SLE.2/GEO.data/GSE61635_series_matrix.txt',getGPL=F)
dat <- exprs(gse)          #### Get expression data

# If data has already been log-transformed, display "log2 transform not need"
# If data has not been log-transformed, perform log2 transformation as per code and display "log2 transform finished"
ex <- dat
qx <- as.numeric(quantile(ex,c(0.,0.25,0.5,0.75,0.99,1.0),na.rm=T))
LogC <-(qx[5] > 100)||
  (qx[6]-qx[1] > 50 && qx[2] > 0)||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if(LogC) {
  ex[which(ex <= 0)] <- NA
  dat <-log2(ex)
  print("log2 transform finished")
} else {
  print("log2 transform not need")
}

expr <- as.data.frame(dat) 
class(expr)
expr[1:4,1:4]
expr$ID <- rownames(expr)

library(tidyverse)
gset <- read_tsv("d:/dualdisease/SLE/SLE.2/GEO.data/GSE61635.top.table.tsv")
gset[1:4,1:6]
dim(gset)
sum(is.na(gset$Gene.symbol))
table(duplicated(gset$Gene.symbol))
colnames(gset)

# Gene annotation - select ID and Gene.symbol columns
ann <- gset[,c("ID","Gene.symbol","logFC","adj.P.Val")]       
ann[1:4,1:4]

# expr$ID <- sub("_at$","",expr$ID)
ann$ID <- as.character(ann$ID)
expr$ID <- as.character(expr$ID)
# expr[expr$ID %like% ann$ID,] # Use left_join instead of direct comparison
dim(expr)
dim(ann)

# Merge "ann" and "expr" data frames by ID using merge function
# expr2 <- ann %>%
#   left_join(expr, by = "ID")
expr2 <- merge(ann, expr, by = "ID")  # Keep all rows from both data frames
expr2[1:4,1:6]

# Remove rows with NA in Gene.Symbol column
expr3 <- expr2[!is.na(expr2$Gene.symbol), ]

# Remove probes without annotated symbols
expr3=expr3[expr3$Gene.symbol !='',]
expr3=expr3[expr3$Gene.symbol !='---',]
expr3=expr3[expr3$Gene.symbol !='--- ',]

# Filter differentially expressed genes (DEGs) with |logFC| > 0.585 and adj.P.Val < 0.05
deg <- expr3[abs(expr3$logFC) > 0.585 & expr3$adj.P.Val < 0.05, ]
nrow(deg)  #2391 

library(stringr)
colnames(deg)
names(deg)[2]<-"Gene.symbol" 
deg$Gene.symbol=trimws(str_split(deg$Gene.symbol,'//',simplify = T)[,1])

# Remove rows with duplicate Gene.Symbol values from deg data frame
deg <- deg[!duplicated(deg$Gene.symbol),]
rownames(deg) <- deg$Gene.symbol
deg <- deg[,-c(1:2)]
deg[1:4,1:4]

# Transpose rows and columns
de <- deg[,c(1:2)]
deg <- deg[,-c(1:2)]
colnames(deg)
deg <- t(deg)
class(deg)
deg[1:4,1:4]
deg <- as.data.frame(deg)
deg$ID <- rownames(deg)

# Get clinical data from GSE object using pData function
clin <- pData(gse)  # pData(gse) extracts phenotype data stored in gse object and stores in clin 
colnames(clin)
clin[1:4,1:35]
clin
clin$ID <- rownames(clin)
clin$`diagnosis:ch1`
table(clin$`diagnosis:ch1`)

# Extract ID and diagnosis information
mc <- clin[,c("ID","diagnosis:ch1")]
mc[1:4,]
names(mc)[2] <- "group"
table(mc$group)
dim(mc)
dim(deg)

# Merge clinical data with expression data
df <- merge(mc,deg,by= "ID")
df[1:4,1:4]
table(df$group)
rownames(df) <- df$ID
df <- df[,-1]
df[1:4,1:4]
class(df$PRSS33)#"numeric"

# Save merged data
write.table(data.frame(df),file = "GSE61635_expr&group.txt", sep = "\t",quote = F,row.names = T)

# Install and load required packages (install if not already installed)
install.packages("plot3D")
library(plot3D)
#BiocManager::install("plot3D")

# Read PCA data file
df[1:4,1:4]
class(df)

# Extract sample group information
dfGroup = as.data.frame(df[,1])
rownames(dfGroup) <- rownames(df)
names(dfGroup)[1] <- "Group"
head(dfGroup)

# Perform PCA analysis
pca_result <- prcomp(df[,-1],  # Perform principal component analysis on all columns except first column of df
                     scale=T   # Logical value indicating whether variables should be scaled to unit variance before analysis
)
pca_result$x<-data.frame(pca_result$x)

# Set colors (one color per group)
colors <- c("red","blue")

# Convert factor values in first column of dfGroup to numeric factors and get corresponding colors from colors
myColors <- colors[as.numeric(as.factor(dfGroup[,1]))]

# Set shape numbers (22:24 represent different point shapes)
my_pch = 23:24

# Select corresponding point shapes based on Species column of iris dataset
pchs = my_pch[as.numeric(as.factor(dfGroup[,1]))]

# Calculate PC variance percentages for axis labels
pVar <- pca_result$sdev^2/sum(pca_result$sdev^2)
pVar = round(pVar,digits = 3)
xName = paste0("PC1 (",as.character(pVar[1] * 100 ),"%)") # Convert first PC variance to percentage for PCA1 label
yName = paste0("PC2 (",as.character(pVar[2] * 100 ),"%)")
zName = paste0("PC3 (",as.character(pVar[3] * 100 ),"%)")

# Plot 3D PCA
pdf("PCA3D_DEG_GSE61635_NEW.pdf",width = 5,height = 5)
with(data.frame(pca_result$x[,1:3]), plot3D::scatter3D(
  x = PC1,       # Set X-axis mapping based on column names in data
  y = PC2,
  z = PC3,
  pch = 21,      # Set scatter point shape to solid circle #my_pch
  cex = 1,       # Set scatter point size
  col=NA,        
  bg=myColors,   # Set scatter point background color
  xlab = xName,  # Set x-axis label
  ylab = yName,
  zlab = zName,
  ticktype = "simple",  # Set tick type (alternative: detailed)
  bty = "b2",    # Set plot theme b
  box = T,       # Show axis box
  theta = 5,     # Set rotation angle
  phi = 15,      # Set viewing angle
  d=10,
  colkey = F     # Hide color key
))

# Add legend
legend("bottom",     # Specify legend position at bottom
       title = "",   # Set legend title to empty
       legend = unique(dfGroup$Group),  # Use groups from dfGroup for labels
       pch = 21,     # Set legend marker shape to solid circle
       pt.cex = 1,   # Set legend marker size
       cex = 0.8,    # Set legend text size
       pt.bg = colors,  # Set legend marker background color
       bg = "white",  # Set legend background color
       bty = "n",     # Remove legend border
       horiz = TRUE   # Set legend to horizontal display
)
dev.off()

# Volcano plot analysis
# Create new column based on conditional judgment: 
# "up" if logFC > 0.585 and adj.P.Val < 0.05;
# "down" if logFC < -0.585 and adj.P.Val < 0.05;
# "ns" (not significant) otherwise
expr2$group <- ifelse(expr2$adj.P.Val < 0.05 & expr2$logFC > 0.585, "up", 
                      ifelse(expr2$adj.P.Val < 0.05 & expr2$logFC < -0.585, "down", "ns"))

# Calculate frequency of each value in new group column using table()
table(expr2$group)

# Clean gene symbols
expr2$Gene.symbol=trimws(str_split(expr2$Gene.symbol,'//',simplify = T)[,1])

# Remove duplicate genes
# Separate non-significant and significant genes
expr2_ns <- expr2[expr2$group == "ns", ]
expr2_up_down <- expr2[expr2$group %in% c("up", "down"), ]

# Remove duplicate genes in up/down groups (based on Gene.symbol column)
expr2_up_down <- expr2_up_down[!duplicated(expr2_up_down$Gene.symbol), ]

# Remove genes from ns group that already exist in up/down groups
expr2_ns <- expr2_ns[!expr2_ns$Gene.symbol %in% expr2_up_down$Gene.symbol, ]

# Merge datasets
expr2_combined <- rbind(expr2_up_down, expr2_ns)

# Keep first occurrence of each gene (prioritize genes in expr2_up_down)
expr2_combined <- expr2_combined[!duplicated(expr2_combined$Gene.symbol), ]
table(expr2_combined$group)

# Remove rows with NA in Gene.Symbol column
expr2_combined <- expr2_combined[!is.na(expr2_combined$Gene.symbol), ]

# Remove probes without annotated symbols
expr2_combined=expr2_combined[expr2_combined$Gene.symbol !='',]
expr2_combined=expr2_combined[expr2_combined$Gene.symbol !='---',]
expr2_combined=expr2_combined[expr2_combined$Gene.symbol !='--- ',]

# Set Gene.symbol column as row names
rownames(expr2_combined) <- expr2_combined$Gene.symbol
expr2_combined[1:4,]
colnames(expr2_combined)

# Select key columns
expr2_com <- expr2_combined[,c(2:4,134)]
colnames(expr2_com)

# Match sample IDs
df_row <- rownames(df)
ix <- expr2_combined[, match(df_row, colnames(expr2_combined), nomatch = 0)]
ix$Gene.symbol <- rownames(ix)
ig <- merge(expr2_com,ix,by = "Gene.symbol")
rownames(ig) <- ig$Gene.symbol
ig <- ig[,-1]

# Move group column to first position
ig <- ig[, c("group", setdiff(names(ig), "group"))]
colnames(ig)

# Extract only up/down regulated genes
updown <- ig[ig$group != "ns", ]
table(updown$group)
# down   up 
# 1487 2260

# Output files
write.table(ig,file = "GSE61635_expression.txt",sep = "\t",row.names=T,col.names = T)
write.table(updown,file = "GSE61635_expression_updown.txt",sep = "\t",row.names=T,col.names = T)

# Volcano plot generation
library(AMR)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(scales)

ig$Gene.symbol <- rownames(ig)

# Convert to tibble for easier processing (Tibble is an enhanced data frame type from dplyr with improved features)
dt <- as_tibble(ig)

# View first few rows of tibble
head(dt)

# Add new column log10FDR with negative log10 of adj.P.Val
dt$log10FDR <- -log10(dt$adj.P.Val)  # Note: using adj.P.Val!!

# View first few rows of dt
head(dt)
table(dt$group)

# Get top 10 most significantly differentially expressed genes
# Filter rows where group is not "None",
# remove duplicates by Gene.symbol (keeping first instance with distinct function),
# and select top 10 rows by absolute logFC value
top10sig <- filter(dt,group!="ns") %>% distinct(Gene.symbol,.keep_all = T) %>% top_n(10,abs(logFC))

# View top10sig gene symbols
top10sig$Gene.symbol

# Extract up-regulated genes from top10sig
up <- filter(top10sig,group=="up")
# Select 5th column of up data frame
up[5]

# Extract down-regulated genes from top10sig
down <- filter(top10sig,group=="down")
# Select 5th column of down data frame
down[5]

# Add new column to mark top 10 DEGs with size 2 and others with size 1
# If dt$ID not in top10sig$ID, size = 1;
# If dt$ID in top10sig$ID, size = 2
dt$size <- case_when(!(dt$Gene.symbol %in% top10sig$Gene.symbol)~ 1,
                     dt$Gene.symbol %in% top10sig$Gene.symbol ~ 2)

# View first few rows of dt
head(dt)
table(dt$size)

# Extract non-top10 genes
dt_f <- filter(dt,size==1)

# View first few rows of df
head(dt_f)

# Specify plotting order by converting group column to factor type
# Convert group column in df to ordered factor
# Specifies factor levels and sets factor as ordered
# Levels are specified as "Up", "Down", and "None"
dt_f$group <- factor(dt_f$group,
                     levels = c("up","down","ns"),
                     ordered = T)

#df$log10FDR <- -log10(df$P.Value)
# Start plotting - create mapping
# Create scatter plot object p1 using ggplot2 package where:
# data = df specifies data source as df data frame
# aes(logFC, log10FDR, color = group) sets x-axis to logFC, y-axis to log10FDR, and colors by group column
# geom_point(size = 1.6) adds scatter point layer with point size 1.6
p1 <- ggplot(data=dt_f,aes(logFC,log10FDR,color=group))+
  geom_point(size=1.6)

# scale_colour_manual() function for manual color setting
# name parameter empty means no legend title
# values parameter specifies color vector, alpha() function sets color transparency
# p21 adds color scale to scatter plot p1 to get new plot object
mycolor <- c("#f0d0f0","#a8e09f","gray80")

p21 <- p1 + scale_colour_manual(name="",values=alpha(mycolor,0.7))
p21

# Add points for top 10 genes
# geom_point() to plot data points for "Up" and "Down" groups separately
p2 <- p21+geom_point(data=up,aes(logFC,log10FDR,),
                     color="#d8b0d8",size=3,alpha=0.9)+
  geom_point(data=down,aes(logFC,log10FDR),
             color="#88d078",size=3,alpha=0.9)
p2

# expansion function sets size of blank areas at both ends of axis range; mult is "multiplier" mode, add is "additive" mode
# labs(y = "-log10FDR") sets y-axis label to "-log10FDR"
# scale_y_continuous() function adjusts y-axis ticks
# limits parameter specifies y-axis range as 0 to 12
# breaks parameter specifies tick values as 0, 4, 8, and 12
# label parameter specifies corresponding labels as "0", "4", "8", and "12"
summary(dt$logFC) # Determine X-axis range
summary(dt$log10FDR) # Determine Y-axis range

p3<-p2+labs(y="-log10FDR")+
  scale_y_continuous(#expand=expansion(add = c(2, 0)),
    limits = c(0, 39),
    breaks = c(0,13,26,39),
    label = c("0","13","26","39"))+
  scale_x_continuous(limits = c(-6, 6),#### To include top genes
                     breaks = c(-6,-3,0,3,6),
                     label = c("-6","-3","0","3","6"))
p3

# Add arrows and gene labels
# geom_text_repel() function adds text labels
set.seed(007)
p4 <- 
  p3+geom_text_repel(data=top10sig,aes(logFC,log10FDR,label=Gene.symbol),
                     force=80,color="grey20",size=3,
                     point.padding = 0.5,hjust = 0.5,   # Control distance between text label and data point # Control horizontal alignment of text label (0=left, 1=right, 0.5=center)
                     arrow = arrow(length = unit(0.01, "npc"),  # Add arrow pointing to text label with length as percentage of data point
                                   type = "open", ends = "last"),  # Type is "open", arrow ends at "last"
                     segment.color="grey20", # Arrow color
                     segment.size=0.2,       # Arrow size
                     segment.alpha=0.8,      # Arrow transparency
                     nudge_x=0, # Small adjustment of text label in x and y directions
                     nudge_y=0.5
  )
p4+theme_bw()  # Set plot theme to white background

# Custom plot theme with fine adjustments
# Control size of top, right, bottom, and left margins in inches
top.mar=0.2
right.mar=0.2
bottom.mar=0.2
left.mar=0.2

# Hide y-axis and set font style, axis line thickness, color, and tick length
mytheme<-theme_classic()+
  theme(text=element_text(family = "sans",colour ="gray30",size = 12),  # Set font for text in plot
        axis.line = element_line(size = 0.6,colour = "gray30"),         # Control appearance of axis lines
        axis.ticks = element_line(size = 0.6,colour = "gray30"),        # Control appearance of tick marks
        axis.ticks.length = unit(1.5,units = "mm"),                     # Control length of tick marks in millimeters
        plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),    # Control plot margin size in inches
                         units="inches"))

# Apply custom theme and save plot
pdf("volcano_GSE61635.pdf",width = 6,height = 5)
p4+mytheme  +geom_hline(yintercept = c(-log10(0.05)),  # Draw horizontal dashed line with width 0.7 to mark significance threshold (0.05)
                        linewidth = 0.7,
                        color = "black",
                        lty = "dashed")+
  geom_vline(xintercept = c(-0.585,0.585),  # Draw two vertical dashed lines to mark fold change thresholds (-0.585 and 0.585)
             linewidth = 0.7,
             color = "black",
             lty = "dashed")
dev.off()

# View top 10 significant genes
top10sig$Gene.symbol
# [1] "CHI3L1"  "CXCL2"   "CXCL8"   "EGR1"    "IFI27"   "IL1R2"   "OLFM4"   "PGLYRP1" "RPS4Y1" 
# [10] "SIGLEC1"

# Heatmap generation
library('GEOquery')
library('data.table')
library('tidyr')
library('tibble')
library('dplyr')
library(limma)

# Preview data
ig[1:4,1:4]
colnames(ig)

# Remove non-expression columns
ex <- ig[,-c(1:3,133)]
colnames(ex)

# Convert data frame to matrix format
ex <- as.matrix(ex)
ex[1:4,1:4]
class(ex)

# Filter rows using apply function
# Keep rows where number of zeros is less than 80% of total columns
pr3 <- ex[apply(ex,1,function(x)sum(x==0)<0.8*ncol(ex)),]

# Normalize expression matrix (scale each row)
geo.expr <- t(scale(t(pr3)))%>%as.data.frame()    #### Qualified normalized expression matrix

# Get sample IDs for SLE group
table(dfGroup$Group)
t_id <- row.names(filter(dfGroup,`Group`=='systemic lupus erythematosus (SLE)'))

# Get sample IDs for control group
n_id <- row.names(filter(dfGroup,`Group`=='healthy'))

# Select columns corresponding to SLE and control samples from geo.expr
geo.expr1 <- geo.expr[,c(t_id,n_id)]   #### control--8 samples  UC---13 samples  Columns ordered as UC then control

# Generate heatmap
library(pheatmap)

# Get row names of DEGs
DEGs <- rownames(updown)
#rownames(sigdiff)

# Extract DEG expression data and limit value range
pre_heatdata <- geo.expr1[DEGs,]
pre_heatdata[pre_heatdata >1] <- 1
pre_heatdata[pre_heatdata < -1] <- -1

# Create annotation data frame with SLE and Control groups
ann <- data.frame(Type=factor(rep(c('SLE', 'Control'), c(99, 30))),
                  row.names = colnames(geo.expr1))

# Ensure color mapping names match factor levels in data frame exactly
ann_colors <- list(Type = c("SLE" = "#b83db8", "Control" = "#75af6c"))

# Save heatmap to PDF
pdf("pheatmap_GSE61635.pdf",width = 6,height = 5)
pheatmap(pre_heatdata,
         color = colorRampPalette(c("#88d078","#a8e09f","#FFFFF0", "#f0d0f0","#d8b0d8"))(256),
         show_rownames = F,
         show_colnames = F,
         cluster_cols =F,
         scale = "none",
         treeheight_row = 0,
         treeheight_col = F,
         #cluster_rows = T,
         annotation_col = ann,
         annotation_colors = ann_colors)  # Set annotation colors)    
#c("darkgreen","green","black", "red","firebrick3")
dev.off()

# Number of rows in heatmap data
nrow(pre_heatdata)#[1]2944

# Save up/down regulated gene lists
deg_SLE2 <- as.data.frame(dt_f) 
deg_SLE2[1:4,]
deg_SLE2 <- deg_SLE2[,c("group","Gene.symbol")]
table(deg_SLE2$group)

# Extract up and down regulated gene symbols
deg_SLE2_up <- deg_SLE2[deg_SLE2$group == "up",]$Gene.symbol
deg_SLE2_down <- deg_SLE2[deg_SLE2$group == "down",]$Gene.symbol

# Save gene lists to RData file
save(deg_SLE2_up, deg_SLE2_down, file = "DEG_SLE2_up_down.RData")
load("DEG_SLE2_up_down.RData")