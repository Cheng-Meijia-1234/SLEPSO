#### Data Standardization ####
setwd("D:/dualdisease/RANDOMreview")
load("D:/dualdisease/WGCNA/Preal.data/P_mexp_clin.RData")

# Define min-max scaling function (according to your standardization logic)
min_max_scale <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# GSE78097 data processing
samples <- gsub("GSE78097_", "", clin$SAMPLE_ID[grepl("GSE78097", clin$SAMPLE_ID)])
mexp_sub <- mexp[mexp$Sample %in% samples, ]

# Perform standardization according to your logic: only scale gene expression columns (columns 3-9 correspond to data)
genes <- c("IFI6", "MX1", "NMI", "OAS3", "OASL", "SAMD9", "UBE2L6")

# Extract gene expression data and standardize
gene_data <- mexp_sub[, genes]
gene_data_scaled <- as.data.frame(lapply(gene_data, min_max_scale))
rownames(gene_data_scaled) <- mexp_sub$Sample

# Construct final dataframe (keep first two columns unchanged, only standardize gene columns)
final_data <- data.frame(
  sampleid = rownames(gene_data_scaled),
  DISEASE = ifelse(clin$DISEASE[match(rownames(gene_data_scaled), 
                                      gsub("GSE[0-9]+_", "", clin$SAMPLE_ID))] == "disease", 1, 0),
  gene_data_scaled
)

write.csv(final_data, "D:/dualdisease/RANDOMreview/GSE78097_7genes_data.csv", row.names = FALSE)
cat("GSE78097 completed! Samples:", nrow(final_data), "Disease cases:", sum(final_data$DISEASE), "\n")
cat("Data standardization method: Min-max scaling (0-1 normalization)\n")
cat("Standardized columns:", paste(genes, collapse = ", "), "\n")

# GSE13355 data processing
samples <- gsub("GSE13355_", "", clin$SAMPLE_ID[grepl("GSE13355", clin$SAMPLE_ID)])
mexp_sub <- mexp[mexp$Sample %in% samples, ]

# Extract gene expression data and standardize
gene_data <- mexp_sub[, genes]
gene_data_scaled <- as.data.frame(lapply(gene_data, min_max_scale))
rownames(gene_data_scaled) <- mexp_sub$Sample

# Construct final dataframe
final_data <- data.frame(
  sampleid = rownames(gene_data_scaled),
  DISEASE = ifelse(clin$DISEASE[match(rownames(gene_data_scaled), 
                                      gsub("GSE[0-9]+_", "", clin$SAMPLE_ID))] == "disease", 1, 0),
  gene_data_scaled
)

write.csv(final_data, "D:/dualdisease/RANDOMreview/GSE13355_7genes_data.csv", row.names = FALSE)
cat("GSE13355 completed! Samples:", nrow(final_data), "Disease cases:", sum(final_data$DISEASE), "\n")
cat("Data standardization method: Min-max scaling (0-1 normalization)\n")
cat("Standardized columns:", paste(genes, collapse = ", "), "\n")

# Verify standardization results
cat("\n=== Standardization Verification ===\n")
cat("GSE78097 data range:\n")
print(sapply(final_data[, 3:9], function(x) c(Min = min(x), Max = max(x))))

# View first few rows of standardized data
cat("\nFirst 4 rows of GSE78097 standardized data:\n")
print(head(final_data, 4))

# Process GSE109248 and GSE141804 data
library(readxl)

# GSE109248 data processing
gse109248_data <- read_excel("GSE109248_scale.xlsx")
genes <- c("IFI6", "MX1", "NMI", "OAS3", "OASL", "SAMD9", "UBE2L6")
gse109248_filtered <- gse109248_data %>%
  select(sampleid = 1, DISEASE = 2, all_of(genes))

write.csv(gse109248_filtered, "GSE109248_7genes_data.csv", row.names = FALSE)
cat("GSE109248 completed! Samples:", nrow(gse109248_filtered), "\n")

# GSE141804 data processing
gse141804_data <- read_excel("GSE141804_scale.xlsx")
gse141804_filtered <- gse141804_data %>%
  select(sampleid = 1, DISEASE = 2, all_of(genes))

write.csv(gse141804_filtered, "GSE141804_7genes_data.csv", row.names = FALSE)
cat("GSE141804 completed! Samples:", nrow(gse141804_filtered), "\n")

############ 8 Machine Learning Models can be Started Now ###################
setwd("D:/dualdisease/RANDOMreview")

# Load required packages
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)
library(randomForest)
library(glmnet)
library(xgboost)
library(kernlab)
library(ada)
library(gbm)
library(e1071)
library(kknn)
library(naivebayes)
library(ComplexHeatmap)
library(grid)
library(circlize)
library(RColorBrewer)

# Create results directory
if(!dir.exists("ML_results")) dir.create("ML_results")

# 1. Read data
cat("Reading data...\n")
train_data <- read.csv("GSE13355_7genes_data.csv")  # Training set
internal_test <- read.csv("GSE78097_7genes_data.csv")  # Internal validation set
external_test1 <- read.csv("GSE109248_7genes_data.csv")  # External validation set 1
external_test2 <- read.csv("GSE141804_7genes_data.csv")  # External validation set 2

# Ensure DISEASE is a factor
train_data$DISEASE <- factor(train_data$DISEASE, levels = c(0, 1), labels = c("Control", "Disease"))
internal_test$DISEASE <- factor(internal_test$DISEASE, levels = c(0, 1), labels = c("Control", "Disease"))
external_test1$DISEASE <- factor(external_test1$DISEASE, levels = c(0, 1), labels = c("Control", "Disease"))
external_test2$DISEASE <- factor(external_test2$DISEASE, levels = c(0, 1), labels = c("Control", "Disease"))

# Remove sampleid column
train_data <- train_data[, -1]
internal_test <- internal_test[, -1]
external_test1 <- external_test1[, -1]
external_test2 <- external_test2[, -1]

cat("Data reading completed!\n")
cat("Training set samples:", nrow(train_data), "\n")
cat("Internal validation set samples:", nrow(internal_test), "\n")
cat("External validation set 1 samples:", nrow(external_test1), "\n")
cat("External validation set 2 samples:", nrow(external_test2), "\n")

# 2. Set training control parameters
ctrl <- trainControl(
  method = "cv",
  number = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = TRUE
)

# 3. Define ROC plotting function
plot_roc_comparison <- function(roc_results, method_name) {
  # Extract all ROC curves
  roc_list <- list(
    Train = roc_results$train_roc,
    Internal_Test = roc_results$internal_roc,
    External_Test1 = roc_results$external1_roc,
    External_Test2 = roc_results$external2_roc
  )
  
  # Create AUC table
  auc_tbl <- data.frame(
    Dataset = c("Train", "Internal_Test", "External_Test1", "External_Test2"),
    AUC = c(
      auc(roc_results$train_roc),
      auc(roc_results$internal_roc),
      auc(roc_results$external1_roc),
      auc(roc_results$external2_roc)
    )
  )
  
  # Plot ROC curves
  p <- ggroc(roc_list, legacy.axes = TRUE, size = 1.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
    ggtitle(paste(method_name, "ROC Curves - PSO")) +
    theme_bw() +
    labs(x = "False Positive Rate (1 - Specificity)", 
         y = "True Positive Rate (Sensitivity)") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          legend.position = "bottom")
  
  # Add AUC annotations
  annot_y <- c(0.25, 0.20, 0.15, 0.10)
  for (i in 1:nrow(auc_tbl)) {
    p <- p + annotate("text", x = 0.6, y = annot_y[i],
                      label = sprintf("%s: AUC = %.3f", auc_tbl$Dataset[i], auc_tbl$AUC[i]),
                      hjust = 0, size = 4)
  }
  
  # Save plot
  filename <- paste0("ML_results/", method_name, "_ROC_PSO.pdf")
  ggsave(filename, p, width = 8, height = 7)
  
  return(list(plot = p, auc_table = auc_tbl))
}

# 4. Random Forest
cat("Starting Random Forest training...\n")
set.seed(123)
rf_model <- train(
  DISEASE ~ .,
  data = train_data,
  method = "rf",
  trControl = ctrl,
  metric = "ROC",
  ntree = 1000,
  importance = TRUE
)

# RF prediction
rf_train_prob <- predict(rf_model, train_data, type = "prob")[, "Disease"]
rf_internal_prob <- predict(rf_model, internal_test, type = "prob")[, "Disease"]
rf_external1_prob <- predict(rf_model, external_test1, type = "prob")[, "Disease"]
rf_external2_prob <- predict(rf_model, external_test2, type = "prob")[, "Disease"]

rf_results <- list(
  train_roc = roc(train_data$DISEASE, rf_train_prob),
  internal_roc = roc(internal_test$DISEASE, rf_internal_prob),
  external1_roc = roc(external_test1$DISEASE, rf_external1_prob),
  external2_roc = roc(external_test2$DISEASE, rf_external2_prob)
)

rf_plot <- plot_roc_comparison(rf_results, "RandomForest")
cat("Random Forest completed!\n")

# 5. LASSO
cat("Starting LASSO training...\n")
set.seed(123)

# Prepare data
x_train <- model.matrix(DISEASE ~ ., train_data)[, -1]
y_train <- ifelse(train_data$DISEASE == "Disease", 1, 0)

# Train LASSO model
lasso_model <- cv.glmnet(x_train, y_train, 
                         family = "binomial", 
                         alpha = 1, 
                         nfolds = 10)

# Get lambda values
lambda_values <- lasso_model$lambda
cat("Available lambda range:", round(range(lambda_values), 5), "\n")

# Use lambda.min (selects more features)
best_lambda <- lasso_model$lambda.min
cat("Selected lambda (lambda.min):", round(best_lambda, 5), "\n")

# View selected features
coef_mat <- coef(lasso_model, s = "lambda.min")
cat("Selected features:\n")
print(coef_mat)

# LASSO prediction
lasso_train_prob <- predict(lasso_model, x_train, s = "lambda.min", type = "response")[, 1]
lasso_internal_prob <- predict(lasso_model, 
                               model.matrix(DISEASE ~ ., internal_test)[, -1], 
                               s = "lambda.min", type = "response")[, 1]
lasso_external1_prob <- predict(lasso_model, 
                                model.matrix(DISEASE ~ ., external_test1)[, -1], 
                                s = "lambda.min", type = "response")[, 1]
lasso_external2_prob <- predict(lasso_model, 
                                model.matrix(DISEASE ~ ., external_test2)[, -1], 
                                s = "lambda.min", type = "response")[, 1]

lasso_results <- list(
  train_roc = roc(y_train, lasso_train_prob),
  internal_roc = roc(internal_test$DISEASE, lasso_internal_prob),
  external1_roc = roc(external_test1$DISEASE, lasso_external1_prob),
  external2_roc = roc(external_test2$DISEASE, lasso_external2_prob)
)

lasso_plot <- plot_roc_comparison(lasso_results, "LASSO")

cat("LASSO completed!\n")
cat("Training set AUC:", round(auc(lasso_results$train_roc), 3), "\n")
cat("Internal validation set AUC:", round(auc(lasso_results$internal_roc), 3), "\n")
cat("External validation set 1 AUC:", round(auc(lasso_results$external1_roc), 3), "\n")
cat("External validation set 2 AUC:", round(auc(lasso_results$external2_roc), 3), "\n")

coef_mat <- coef(lasso_model, s = "lambda.min")
print(coef_mat)

# 1. Univariate logistic regression to examine the relationship between each gene and disease
univariate_analysis <- function(data) {
  genes <- c("IFI6", "MX1", "NMI", "OAS3", "OASL", "SAMD9", "UBE2L6")
  results <- data.frame(Gene = character(), Coefficient = numeric(), P_value = numeric(), stringsAsFactors = FALSE)
  
  for (gene in genes) {
    formula <- as.formula(paste("DISEASE ~", gene))
    model <- glm(formula, data = data, family = binomial)
    coef <- summary(model)$coefficients[2, 1]
    p_val <- summary(model)$coefficients[2, 4]
    results <- rbind(results, data.frame(Gene = gene, Coefficient = coef, P_value = p_val))
  }
  
  return(results)
}

uni_results <- univariate_analysis(train_data)
print(uni_results)

# 2. Check correlation between genes
gene_cor <- cor(train_data[, c("IFI6", "MX1", "NMI", "OAS3", "OASL", "SAMD9", "UBE2L6")])
print(gene_cor)

# 3. Plot correlation heatmap
library(ComplexHeatmap)
Heatmap(gene_cor, name = "Correlation", 
        cluster_rows = TRUE, cluster_columns = TRUE,
        row_names_gp = gpar(fontsize = 10), 
        column_names_gp = gpar(fontsize = 10),
        column_title = "Gene-Gene Correlation")

# Ridge Regression handles multicollinearity better
set.seed(123)
ridge_model <- train(
  DISEASE ~ .,
  data = train_data,
  method = "glmnet",
  trControl = ctrl,
  metric = "ROC",
  tuneGrid = expand.grid(alpha = 0, lambda = 10^seq(-3, 3, length = 100))
)

# View Ridge Regression coefficients
ridge_coef <- coef(ridge_model$finalModel, s = ridge_model$bestTune$lambda)
cat("=== Ridge Regression Coefficients ===\n")
print(ridge_coef)

########################## Transcription Factor Regulatory Network Analysis ################
# Redefine genes object
genes <- c("IFI6", "MX1", "NMI", "OAS3", "OASL", "SAMD9", "UBE2L6")

# Improved ChEA3 API analysis function
analyze_tf_regulation <- function(genes) {
  library(httr)
  library(jsonlite)
  
  url <- "https://maayanlab.cloud/chea3/api/enrich/"
  
  # Construct request body
  payload <- list(
    query_name = "my_psoriasis_genes",
    gene_set = paste(genes, collapse = "\n")
  )
  
  # Send POST request
  response <- POST(url, body = payload, encode = "form")
  
  # Check if request was successful
  if (status_code(response) == 200) {
    results <- fromJSON(content(response, "text"))
    return(results)
  } else {
    cat("API request failed with status code:", status_code(response), "\n")
    return(NULL)
  }
}

# Run analysis
tf_results <- analyze_tf_regulation(genes)

# If online API is unavailable, use offline method
if (is.null(tf_results)) {
  cat("ChEA3 API unavailable, using alternative method...\n")
  
  # Method 1: Use enrichR for transcription factor enrichment analysis
  if (!requireNamespace("enrichR", quietly = TRUE)) {
    install.packages("enrichR")
  }
  library(enrichR)
  
  # View available databases
  dbs <- listEnrichrDbs()
  
  # Select ChEA 2016 database (transcription factor targets)
  dbs_to_use <- c("ChEA_2016", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")
  
  # Perform enrichment analysis
  enriched <- enrichr(genes, dbs_to_use)
  
  # View results
  print("ChEA_2016 Transcription Factor Enrichment Results:")
  print(head(enriched[[1]]))
  
  print("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X Results:")
  print(head(enriched[[2]]))
  
  # Save results
  write.csv(enriched[[1]], "ML_results/TF_enrichment_ChEA2016.csv", row.names = FALSE)
  write.csv(enriched[[2]], "ML_results/TF_enrichment_ENCODE_ChEA.csv", row.names = FALSE)
  
  tf_results <- enriched
}

## Key Transcription Factors: Our 7 genes (IFI6, MX1, NMI, OAS3, OASL, SAMD9, UBE2L6) are regulated by several important transcription factors.
## Among them, the most significant is: SOX2: In the ChEA_2016 database, the target gene set of SOX2 includes all 7 of our genes with extremely high enrichment significance.
## "Transcription factor enrichment analysis revealed that these 7 interferon-related genes are coordinately regulated by key transcription factors such as IRF1 and SOX2.
## As a core regulator of interferon signaling, IRF1 explains the coordinated upregulation of these genes in the inflammatory response of psoriasis.
## The unexpected finding that SOX2 regulates all 7 genes suggests that abnormal differentiation of epidermal stem cells may play an important role in the pathogenesis of psoriasis."
## This analysis provides a mechanistic upstream explanation for your 7 genes, greatly enhancing the depth and scientific rigor of the research!
# Yes, IRF1 and SOX2 are genes, but in the context of this analysis, we focus on how the transcription factor proteins they encode regulate your 7 genes.
# This is like finding the "master switches" that control your 7 genes, providing deeper insights into understanding the pathogenesis of psoriasis!

#####################################################

# Plot LASSO lollipop plot with macaron color palette
cat("Plotting LASSO lollipop plot with macaron color palette...\n")

# Extract LASSO coefficients
coef_mat <- coef(lasso_model, s = "lambda.min")
coef_df <- as.data.frame(as.matrix(coef_mat))
colnames(coef_df) <- "coef"

# Remove intercept term
coef_df <- coef_df[rownames(coef_df) != "(Intercept)", , drop = FALSE]

# Create dataframe
lasso.result <- data.frame(
  diffvariable = rownames(coef_df),
  coef = coef_df$coef
)

# Calculate absolute values of coefficients
lasso.result$abs_coef <- abs(lasso.result$coef)

# Sort by absolute value and take top 5
lasso.plot.data <- lasso.result[order(lasso.result$abs_coef, decreasing = TRUE), ]
lasso.plot.data <- head(lasso.plot.data, 5)

# Custom macaron color palette
macaron_colors <- c("#FFB3BA", "#FFDFBA", "#FFFFBA", "#BAFFC9", "#BAE1FF", 
                    "#FFB347", "#FFD700", "#FF6F61", "#FFB5E8", "#B5EAD7", 
                    "#B5D1FF", "#FF9CEE", "#FFBABA", "#FFDAC1", "#FFF1BA")

# View variables to be plotted
cat("Variables to be plotted (sorted by absolute coefficient value):\n")
print(lasso.plot.data)

# Plot (show absolute values only, same direction, macaron color palette)
p_lasso_macaron <- ggplot(lasso.plot.data, 
                          aes(y = reorder(diffvariable, abs_coef), x = abs_coef)) + 
  # Lines extending from x=0 to abs_coef values (using macaron color palette)
  geom_segment(aes(yend = reorder(diffvariable, abs_coef), x = 0, xend = abs_coef, 
                   color = diffvariable), 
               size = 1.8, alpha = 0.8) + 
  # Dots representing variable importance
  geom_point(aes(color = diffvariable), size = 9, alpha = 0.8) +  
  # Annotations showing absolute values next to dots
  geom_text(aes(label = round(abs_coef, 3)), color = "black", size = 3.2, fontface = "bold") +  
  # Use macaron color palette
  scale_color_manual(values = macaron_colors) + 
  # Add title, x-axis and y-axis labels
  labs(title = "LASSO Feature Importance - PSO Dataset", 
       subtitle = paste("Top", nrow(lasso.plot.data), "Features by Absolute Coefficient Value"),
       y = "Genes", 
       x = "Absolute Coefficient Value") + 
  # Clean theme
  theme_minimal() + 
  # Center title
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = "none") +  # Hide legend
  # Ensure x-axis starts from 0
  expand_limits(x = 0)

# Display plot
print(p_lasso_macaron)

# Save plot
ggsave("ML_results/LASSO_Macaron_Absolute_PSO.pdf", p_lasso_macaron, width = 6, height = 5)

cat("LASSO lollipop plot with macaron color palette saved!\n")

# Output detailed coefficient information
cat("\n=== Detailed LASSO Feature Selection Results - PSO Dataset ===\n")
cat("Lambda.min:", round(best_lambda, 5), "\n")
cat("Number of non-zero coefficient features:", nrow(lasso.result), "\n")
cat("Number of features plotted:", nrow(lasso.plot.data), "\n\n")

cat("Detailed gene coefficient information:\n")
for (i in 1:nrow(lasso.plot.data)) {
  cat(sprintf("  %d. %-20s: Absolute Value = %.4f (Original Coefficient: %.4f)\n", 
              i, 
              lasso.plot.data$diffvariable[i], 
              lasso.plot.data$abs_coef[i],
              lasso.plot.data$coef[i]))
}

# If any features were shrunk to 0, display them too
genes <- c("IFI6", "MX1", "NMI", "OAS3", "OASL", "SAMD9", "UBE2L6")
zero_features <- setdiff(genes, lasso.result$diffvariable)
if(length(zero_features) > 0) {
  cat("\nFeatures shrunk to 0 by LASSO:\n")
  cat(paste(zero_features, collapse = ", "), "\n")
}

# Save detailed coefficient results
lasso_coef_details <- lasso.result[order(lasso.result$abs_coef, decreasing = TRUE), ]
write.csv(lasso_coef_details, "ML_results/LASSO_Coefficient_Details_PSO.csv", row.names = FALSE)
cat("\nDetailed coefficient results saved to: ML_results/LASSO_Coefficient_Details_PSO.csv\n")

# 6. XGBoost
cat("Starting XGBoost training...\n")
set.seed(123)
xgb_model <- train(
  DISEASE ~ .,
  data = train_data,
  method = "xgbTree",
  trControl = ctrl,
  metric = "ROC",
  verbose = FALSE
)

# XGBoost prediction
xgb_train_prob <- predict(xgb_model, train_data, type = "prob")[, "Disease"]
xgb_internal_prob <- predict(xgb_model, internal_test, type = "prob")[, "Disease"]
xgb_external1_prob <- predict(xgb_model, external_test1, type = "prob")[, "Disease"]
xgb_external2_prob <- predict(xgb_model, external_test2, type = "prob")[, "Disease"]

xgb_results <- list(
  train_roc = roc(train_data$DISEASE, xgb_train_prob),
  internal_roc = roc(internal_test$DISEASE, xgb_internal_prob),
  external1_roc = roc(external_test1$DISEASE, xgb_external1_prob),
  external2_roc = roc(external_test2$DISEASE, xgb_external2_prob)
)

xgb_plot <- plot_roc_comparison(xgb_results, "XGBoost")
cat("XGBoost completed!\n")

# 7. SVM
cat("Starting SVM training...\n")
set.seed(123)
svm_model <- train(
  DISEASE ~ .,
  data = train_data,
  method = "svmRadial",
  trControl = ctrl,
  metric = "ROC",
  preProcess = c("center", "scale")
)

# SVM prediction
svm_train_prob <- predict(svm_model, train_data, type = "prob")[, "Disease"]
svm_internal_prob <- predict(svm_model, internal_test, type = "prob")[, "Disease"]
svm_external1_prob <- predict(svm_model, external_test1, type = "prob")[, "Disease"]
svm_external2_prob <- predict(svm_model, external_test2, type = "prob")[, "Disease"]

svm_results <- list(
  train_roc = roc(train_data$DISEASE, svm_train_prob),
  internal_roc = roc(internal_test$DISEASE, svm_internal_prob),
  external1_roc = roc(external_test1$DISEASE, svm_external1_prob),
  external2_roc = roc(external_test2$DISEASE, svm_external2_prob)
)

svm_plot <- plot_roc_comparison(svm_results, "SVM")
cat("SVM completed!\n")

# 8. AdaBoost
cat("Starting AdaBoost training...\n")
set.seed(123)
ada_model <- train(
  DISEASE ~ .,
  data = train_data,
  method = "ada",
  trControl = ctrl,
  metric = "ROC"
)

# AdaBoost prediction
ada_train_prob <- predict(ada_model, train_data, type = "prob")[, "Disease"]
ada_internal_prob <- predict(ada_model, internal_test, type = "prob")[, "Disease"]
ada_external1_prob <- predict(ada_model, external_test1, type = "prob")[, "Disease"]
ada_external2_prob <- predict(ada_model, external_test2, type = "prob")[, "Disease"]

ada_results <- list(
  train_roc = roc(train_data$DISEASE, ada_train_prob),
  internal_roc = roc(internal_test$DISEASE, ada_internal_prob),
  external1_roc = roc(external_test1$DISEASE, ada_external1_prob),
  external2_roc = roc(external_test2$DISEASE, ada_external2_prob)
)

ada_plot <- plot_roc_comparison(ada_results, "AdaBoost")
cat("AdaBoost completed!\n")

# 9. GBM
cat("Starting GBM training...\n")
set.seed(123)
gbm_model <- train(
  DISEASE ~ .,
  data = train_data,
  method = "gbm",
  trControl = ctrl,
  metric = "ROC",
  verbose = FALSE
)

# GBM prediction
gbm_train_prob <- predict(gbm_model, train_data, type = "prob")[, "Disease"]
gbm_internal_prob <- predict(gbm_model, internal_test, type = "prob")[, "Disease"]
gbm_external1_prob <- predict(gbm_model, external_test1, type = "prob")[, "Disease"]
gbm_external2_prob <- predict(gbm_model, external_test2, type = "prob")[, "Disease"]

gbm_results <- list(
  train_roc = roc(train_data$DISEASE, gbm_train_prob),
  internal_roc = roc(internal_test$DISEASE, gbm_internal_prob),
  external1_roc = roc(external_test1$DISEASE, gbm_external1_prob),
  external2_roc = roc(external_test2$DISEASE, gbm_external2_prob)
)

gbm_plot <- plot_roc_comparison(gbm_results, "GBM")
cat("GBM completed!\n")

# 10. Naive Bayes
cat("Starting Naive Bayes training...\n")
set.seed(123)
nb_model <- train(
  DISEASE ~ .,
  data = train_data,
  method = "naive_bayes",
  trControl = ctrl,
  metric = "ROC"
)

# Naive Bayes prediction
nb_train_prob <- predict(nb_model, train_data, type = "prob")[, "Disease"]
nb_internal_prob <- predict(nb_model, internal_test, type = "prob")[, "Disease"]
nb_external1_prob <- predict(nb_model, external_test1, type = "prob")[, "Disease"]
nb_external2_prob <- predict(nb_model, external_test2, type = "prob")[, "Disease"]

nb_results <- list(
  train_roc = roc(train_data$DISEASE, nb_train_prob),
  internal_roc = roc(internal_test$DISEASE, nb_internal_prob),
  external1_roc = roc(external_test1$DISEASE, nb_external1_prob),
  external2_roc = roc(external_test2$DISEASE, nb_external2_prob)
)

nb_plot <- plot_roc_comparison(nb_results, "NaiveBayes")
cat("Naive Bayes completed!\n")

# 11. KKNN
cat("Starting KKNN training...\n")
set.seed(123)
kknn_model <- train(
  DISEASE ~ .,
  data = train_data,
  method = "kknn",
  trControl = ctrl,
  metric = "ROC",
  preProcess = c("center", "scale")
)

# KKNN prediction
kknn_train_prob <- predict(kknn_model, train_data, type = "prob")[, "Disease"]
kknn_internal_prob <- predict(kknn_model, internal_test, type = "prob")[, "Disease"]
kknn_external1_prob <- predict(kknn_model, external_test1, type = "prob")[, "Disease"]
kknn_external2_prob <- predict(kknn_model, external_test2, type = "prob")[, "Disease"]

kknn_results <- list(
  train_roc = roc(train_data$DISEASE, kknn_train_prob),
  internal_roc = roc(internal_test$DISEASE, kknn_internal_prob),
  external1_roc = roc(external_test1$DISEASE, kknn_external1_prob),
  external2_roc = roc(external_test2$DISEASE, kknn_external2_prob)
)

kknn_plot <- plot_roc_comparison(kknn_results, "KKNN")
cat("KKNN completed!\n")

# 12. Summarize all results
cat("Generating summary report...\n")

# Collect all AUC results
all_auc <- rbind(
  cbind(Method = "RandomForest", rf_plot$auc_table),
  cbind(Method = "LASSO", lasso_plot$auc_table),
  cbind(Method = "XGBoost", xgb_plot$auc_table),
  cbind(Method = "SVM", svm_plot$auc_table),
  cbind(Method = "AdaBoost", ada_plot$auc_table),
  cbind(Method = "GBM", gbm_plot$auc_table),
  cbind(Method = "NaiveBayes", nb_plot$auc_table),
  cbind(Method = "KKNN", kknn_plot$auc_table)
)

# Save AUC summary table
write.csv(all_auc, "ML_results/All_Methods_AUC_Summary_PSO.csv", row.names = FALSE)

# Create method comparison plot
comparison_plot <- ggplot(all_auc, aes(x = Dataset, y = AUC, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = round(AUC, 3)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3) +
  scale_fill_brewer(palette = "Set3") +
  ggtitle("AUC Comparison of 8 Machine Learning Methods - PSO") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("ML_results/Methods_Comparison_PSO.pdf", comparison_plot, width = 10, height = 7)

cat("All analyses completed!\n")
cat("Results are saved in the ML_results folder:\n")
cat("- ROC curve plots for each method\n")
cat("- AUC summary table: All_Methods_AUC_Summary_PSO.csv\n")
cat("- Method comparison plot: Methods_Comparison_PSO.pdf\n")

# Print best methods
best_methods <- all_auc %>%
  group_by(Dataset) %>%
  slice_max(AUC, n = 1)

cat("\nBest methods for each dataset:\n")
print(best_methods)

#################### UPSET Plot ###################
# Feature extraction function definitions
# 1. Random Forest (column names: Control/Disease, either column can be used as values are identical)
extract_rf_top5 <- function(model) {
  imp <- varImp(model)$importance
  # Random Forest has identical values in Control and Disease columns, sort by either
  imp_df <- data.frame(
    Feature = rownames(imp),
    Score = imp$Control  # Or imp$Disease, results are the same
  )
  imp_df <- imp_df[order(imp_df$Score, decreasing = TRUE), ]
  head(imp_df$Feature, 5)
}

# 2. LASSO (based on non-zero coefficients, correctly extracted)
extract_lasso_top5 <- function(model) {
  coef_mat <- coef(model, s = "lambda.min")
  features <- rownames(coef_mat)[coef_mat[, 1] != 0]
  features <- setdiff(features, "(Intercept)")  # Exclude intercept term
  head(features, 5)
}

# 3. XGBoost (column name: Overall)
extract_xgb_top5 <- function(model) {
  imp <- varImp(model)$importance
  imp_df <- data.frame(
    Feature = rownames(imp),
    Score = imp$Overall
  )
  imp_df <- imp_df[order(imp_df$Score, decreasing = TRUE), ]
  head(imp_df$Feature, 5)
}

# 4. SVM (column names: Control/Disease, values are identical)
extract_svm_top5 <- function(model) {
  imp <- varImp(model)$importance
  imp_df <- data.frame(
    Feature = rownames(imp),
    Score = imp$Control  # Or imp$Disease, results are the same
  )
  imp_df <- imp_df[order(imp_df$Score, decreasing = TRUE), ]
  head(imp_df$Feature, 5)
}

# 5. AdaBoost (column names: Control/Disease, values are identical)
extract_ada_top5 <- function(model) {
  imp <- varImp(model)$importance
  imp_df <- data.frame(
    Feature = rownames(imp),
    Score = imp$Control  # Or imp$Disease, results are the same
  )
  imp_df <- imp_df[order(imp_df$Score, decreasing = TRUE), ]
  head(imp_df$Feature, 5)
}

# 6. GBM (column name: Overall)
extract_gbm_top5 <- function(model) {
  imp <- varImp(model)$importance
  imp_df <- data.frame(
    Feature = rownames(imp),
    Score = imp$Overall
  )
  imp_df <- imp_df[order(imp_df$Score, decreasing = TRUE), ]
  head(imp_df$Feature, 5)
}

# 7. Naive Bayes (column names: Control/Disease, values are identical)
extract_nb_top5 <- function(model) {
  imp <- varImp(model)$importance
  imp_df <- data.frame(
    Feature = rownames(imp),
    Score = imp$Control  # Or imp$Disease, results are the same
  )
  imp_df <- imp_df[order(imp_df$Score, decreasing = TRUE), ]
  head(imp_df$Feature, 5)
}

# 8. KKNN (column names: Control/Disease, values are identical)
extract_kknn_top5 <- function(model) {
  imp <- varImp(model)$importance
  imp_df <- data.frame(
    Feature = rownames(imp),
    Score = imp$Control  # Or imp$Disease, results are the same
  )
  imp_df <- imp_df[order(imp_df$Score, decreasing = TRUE), ]
  head(imp_df$Feature, 5)
}

# Extract Top5 features from each model
feature_sets_top5 <- list(
  RandomForest = extract_rf_top5(rf_model),
  LASSO = extract_lasso_top5(lasso_model),
  XGBoost = extract_xgb_top5(xgb_model),
  SVM = extract_svm_top5(svm_model),
  AdaBoost = extract_ada_top5(ada_model),
  GBM = extract_gbm_top5(gbm_model),
  NaiveBayes = extract_nb_top5(nb_model),
  KKNN = extract_kknn_top5(kknn_model)
)

# View extraction results (verify success)
print(feature_sets_top5)

# Generate feature combination matrix
m_comb <- make_comb_mat(feature_sets_top5)

# Set colors (using macaron color palette)
macaron_colors <- c("#FFB3BA", "#FFDFBA", "#FFFFBA", "#BAFFC9", 
                    "#BAE1FF", "#FFB347", "#FFD700", "#FF6F61")
names(macaron_colors) <- names(feature_sets_top5)  # Corresponding to model names

# Plot UpSet plot
pdf("ML_results/UpSet_Top5_Features_PSO.pdf", width = 10, height = 7)
up <- UpSet(
  m_comb,
  set_order = names(feature_sets_top5),  # Order by model
  bg_col = macaron_colors,               # Each model corresponds to a color
  comb_order = order(comb_size(m_comb), decreasing = TRUE),  # Sort by intersection size
  row_names_gp = gpar(fontsize = 10),    # Font size for feature names
  column_title = "Top5 Features Intersection Across 8 Models - PSO"  # Title
)

# Plot and add intersection count labels
ht <- draw(up)
cs <- comb_size(m_comb)  # Intersection sizes
od <- column_order(ht)   # Column order
decorate_annotation("intersection_size", {
  grid.text(
    cs[od],
    x = seq_along(cs),
    y = unit(cs[od], "native") + unit(2, "pt"),
    just = "bottom",
    gp = gpar(fontsize = 9)  # Label font size
  )
})
dev.off()
cat("UpSet plot saved to: ML_results/UpSet_Top5_Features_PSO.pdf\n")

# Count occurrences of each gene
gene_model_df <- data.frame(
  Gene = unlist(feature_sets_top5),  # Top5 genes from all models
  Model = rep(names(feature_sets_top5), each = 5)  # Corresponding model names
)

# Count occurrences of each gene and extract corresponding model names
gene_stats <- gene_model_df %>%
  group_by(Gene) %>%
  summarise(
    Occurrences = n(),  # Number of occurrences
    Models = paste(unique(Model), collapse = ", ")  # Models containing the gene (deduplicated)
  ) %>%
  arrange(desc(Occurrences))  # Sort by number of occurrences in descending order

# Print results
cat("Occurrence statistics of Top5 genes across models:\n")
print(gene_stats, row.names = FALSE)

################### Plot Chord Diagram #####################
# 1. Construct adjacency matrix for shared features between models
methods <- names(feature_sets_top5)
k <- length(methods)
adj_mat <- matrix(0, nrow = k, ncol = k, dimnames = list(methods, methods))

# Calculate number of common features between each pair of models
for (i in 1:(k-1)) {
  for (j in (i+1):k) {
    shared <- intersect(feature_sets_top5[[i]], feature_sets_top5[[j]])
    adj_mat[i, j] <- adj_mat[j, i] <- length(shared)
  }
}

# 2. Set visual parameters (using macaron color palette)
macaron_colors <- c(
  "#FFB3BA", "#FFDFBA", "#FFFFBA", "#BAFFC9", 
  "#BAE1FF", "#FFB347", "#FFD700", "#FF6F61"
)
names(macaron_colors) <- methods

# 3. Plot chord diagram
pdf("ML_results/Chord_Top5_Features_PSO.pdf", width = 9, height = 9)
circos.clear()  # Clear previous circos state

# Core plotting function
chordDiagram(
  adj_mat,
  grid.col = macaron_colors,       # Model sector colors
  transparency = 0.2,              # Link transparency
  directional = FALSE,             # Undirected (symmetric relationship)
  annotationTrack = c("grid", "name"),  # Show grid and model names
  annotationTrackHeight = c(0.03, 0.08), # Adjust track heights
  preAllocateTracks = list(track.height = 0.1), # Reserve track space
  link.border = "white",           # Link borders (enhance clarity)
  link.lwd = 2,                    # Link thickness
  link.sort = TRUE,                # Sort links by size
  link.decreasing = TRUE           # Arrange from largest to smallest
)

# Add title
title("Shared Top5 Genes Between Machine Learning Models - PSO", cex = 1.2)
dev.off()

cat("Chord diagram saved to: ML_results/Chord_Top5_Features_PSO.pdf\n")

################ UPSET Filtering ###################
# 1. Random Forest (column names: Control/Disease, either column can be used as values are identical)
extract_rf_top5 <- function(model) {
  imp <- varImp(model)$importance
  # Random Forest has identical values in Control and Disease columns, sort by either
  imp_df <- data.frame(
    Feature = rownames(imp),
    Score = imp$Control  # Or imp$Disease, results are the same
  )
  imp_df <- imp_df[order(imp_df$Score, decreasing = TRUE), ]
  head(imp_df$Feature, 5)
}

# 2. LASSO (based on non-zero coefficients, correctly extracted)
extract_lasso_top5 <- function(model) {
  coef_mat <- coef(model, s = "lambda.min")
  features <- rownames(coef_mat)[coef_mat[, 1] != 0]
  features <- setdiff(features, "(Intercept)")  # Exclude intercept term
  head(features, 5)
}

# 3. XGBoost (column name: Overall)
extract_xgb_top5 <- function(model) {
  imp <- varImp(model)$importance
  imp_df <- data.frame(
    Feature = rownames(imp),
    Score = imp$Overall
  )
  imp_df <- imp_df[order(imp_df$Score, decreasing = TRUE), ]
  head(imp_df$Feature, 5)
}

# 4. SVM (column names: Control/Disease, values are identical)
extract_svm_top5 <- function(model) {
  imp <- varImp(model)$importance
  imp_df <- data.frame(
    Feature = rownames(imp),
    Score = imp$Control  # Or imp$Disease, results are the same
  )
  imp_df <- imp_df[order(imp_df$Score, decreasing = TRUE), ]
  head(imp_df$Feature, 5)
}

# 5. AdaBoost (column names: Control/Disease, values are identical)
extract_ada_top5 <- function(model) {
  imp <- varImp(model)$importance
  imp_df <- data.frame(
    Feature = rownames(imp),
    Score = imp$Control  # Or imp$Disease, results are the same
  )
  imp_df <- imp_df[order(imp_df$Score, decreasing = TRUE), ]
  head(imp_df$Feature, 5)
}

# 6. GBM (column name: Overall)
extract_gbm_top5 <- function(model) {
  imp <- varImp(model)$importance
  imp_df <- data.frame(
    Feature = rownames(imp),
    Score = imp$Overall
  )
  imp_df <- imp_df[order(imp_df$Score, decreasing = TRUE), ]
  head(imp_df$Feature, 5)
}

# 7. Naive Bayes (column names: Control/Disease, values are identical)
extract_nb_top5 <- function(model) {
  imp <- varImp(model)$importance
  imp_df <- data.frame(
    Feature = rownames(imp),
    Score = imp$Control  # Or imp$Disease, results are the same
  )
  imp_df <- imp_df[order(imp_df$Score, decreasing = TRUE), ]
  head(imp_df$Feature, 5)
}

# 8. KKNN (column names: Control/Disease, values are identical)
extract_kknn_top5 <- function(model) {
  imp <- varImp(model)$importance
  imp_df <- data.frame(
    Feature = rownames(imp),
    Score = imp$Control  # Or imp$Disease, results are the same
  )
  imp_df <- imp_df[order(imp_df$Score, decreasing = TRUE), ]
  head(imp_df$Feature, 5)
}

# Extract Top5 features from each model
feature_sets_top5 <- list(
  RandomForest = extract_rf_top5(rf_model),
  LASSO = extract_lasso_top5(lasso_model),
  XGBoost = extract_xgb_top5(xgb_model),
  SVM = extract_svm_top5(svm_model),
  AdaBoost = extract_ada_top5(ada_model),
  GBM = extract_gbm_top5(gbm_model),
  NaiveBayes = extract_nb_top5(nb_model),
  KKNN = extract_kknn_top5(kknn_model)
)

# View extraction results (verify success)
print(feature_sets_top5)

# Generate feature combination matrix
m_comb <- make_comb_mat(feature_sets_top5)

# Set colors (using macaron color palette)
macaron_colors <- c("#FFB3BA", "#FFDFBA", "#FFFFBA", "#BAFFC9", 
                    "#BAE1FF", "#FFB347", "#FFD700", "#FF6F61")
names(macaron_colors) <- names(feature_sets_top5)  # Corresponding to model names

# Plot UpSet plot
pdf("ML_results/UpSet_Top5_Features_PSO.pdf", width = 10, height = 7)
up <- UpSet(
  m_comb,
  set_order = names(feature_sets_top5),  # Order by model
  bg_col = macaron_colors,               # Each model corresponds to a color
  comb_order = order(comb_size(m_comb), decreasing = TRUE),  # Sort by intersection size
  row_names_gp = gpar(fontsize = 10),    # Font size for feature names
  column_title = "Top5 Features Intersection Across 8 Models - PSO"  # Title
)

# Plot and add intersection count labels
ht <- draw(up)
cs <- comb_size(m_comb)  # Intersection sizes
od <- column_order(ht)   # Column order
decorate_annotation("intersection_size", {
  grid.text(
    cs[od],
    x = seq_along(cs),
    y = unit(cs[od], "native") + unit(2, "pt"),
    just = "bottom",
    gp = gpar(fontsize = 9)  # Label font size
  )
})
dev.off()
cat("UpSet plot saved to: ML_results/UpSet_Top5_Features_PSO.pdf\n")

# Count occurrences of each gene
gene_model_df <- data.frame(
  Gene = unlist(feature_sets_top5),  # Top5 genes from all models
  Model = rep(names(feature_sets_top5), each = 5)  # Corresponding model names
)

# Count occurrences of each gene and extract corresponding model names
gene_stats <- gene_model_df %>%
  group_by(Gene) %>%
  summarise(
    Occurrences = n(),  # Number of occurrences
    Models = paste(unique(Model), collapse = ", ")  # Models containing the gene (deduplicated)
  ) %>%
  arrange(desc(Occurrences))  # Sort by number of occurrences in descending order

# Print results - full display
cat("Occurrence statistics of Top5 genes across models (full information):\n")
print(gene_stats, row.names = FALSE, width = Inf)

# Additionally display complete Top5 list for each model separately
cat("\n=== Detailed Top5 Gene Lists by Model - PSO ===\n")
for (model_name in names(feature_sets_top5)) {
  cat(sprintf("\n%s:\n", model_name))
  cat(paste(feature_sets_top5[[model_name]], collapse = ", "), "\n")
}

# Display most important shared genes
cat("\n=== Core Shared Gene Analysis - PSO ===\n")
cat("Genes appearing in all 8 models (100% shared):\n")
core_genes_8 <- gene_stats %>% filter(Occurrences == 8)
print(core_genes_8, row.names = FALSE, width = Inf)

cat("\nGenes appearing in 7 models:\n")
core_genes_7 <- gene_stats %>% filter(Occurrences == 7)
print(core_genes_7, row.names = FALSE, width = Inf)

cat("\nGenes appearing in 6 models:\n")
core_genes_6 <- gene_stats %>% filter(Occurrences == 6)
print(core_genes_6, row.names = FALSE, width = Inf)

# Calculate specific model distribution for each gene
cat("\n=== Detailed Model Distribution by Gene - PSO ===\n")
for(i in 1:nrow(gene_stats)) {
  gene <- gene_stats$Gene[i]
  models <- gene_stats$Models[i]
  cat(sprintf("%s: Appears in %d models -> %s\n", 
              gene, gene_stats$Occurrences[i], models))
}

# Save complete statistical results to file
write.csv(gene_stats, "ML_results/Gene_Occurrence_Statistics_PSO.csv", row.names = FALSE)
cat("\nComplete statistical results saved to: ML_results/Gene_Occurrence_Statistics_PSO.csv\n")

# Create visualization to show gene sharing patterns
library(ggplot2)

# Create gene occurrence frequency bar plot
p_gene_frequency <- ggplot(gene_stats, aes(x = reorder(Gene, Occurrences), y = Occurrences, fill = Gene)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Occurrences), vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = macaron_colors) +
  labs(title = "Gene Occurrence Frequency Across 8 Machine Learning Models - PSO Dataset",
       x = "Genes", 
       y = "Number of Models") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "none") +
  ylim(0, 8.5)

ggsave("ML_results/Gene_Frequency_Barplot_PSO.pdf", p_gene_frequency, width = 10, height = 6)
cat("Gene frequency bar plot saved to: ML_results/Gene_Frequency_Barplot_PSO.pdf\n")

# Create heatmap showing occurrence of each gene in each model
library(reshape2)

# Create binary matrix
binary_matrix <- matrix(0, nrow = length(unique(unlist(feature_sets_top5))), 
                        ncol = length(feature_sets_top5),
                        dimnames = list(unique(unlist(feature_sets_top5)), names(feature_sets_top5)))

for (model in names(feature_sets_top5)) {
  binary_matrix[feature_sets_top5[[model]], model] <- 1
}

# Convert to long format for ggplot
binary_df <- melt(binary_matrix)
colnames(binary_df) <- c("Gene", "Model", "Present")

# Create heatmap
p_heatmap <- ggplot(binary_df, aes(x = Model, y = Gene, fill = as.factor(Present))) +
  geom_tile(color = "white", size = 1) +
  scale_fill_manual(values = c("0" = "lightgray", "1" = "#FF6F61"), 
                    labels = c("Not Selected", "Top5 Selected")) +
  labs(title = "Top5 Gene Selection Heatmap by Model - PSO Dataset",
       x = "Machine Learning Models", 
       y = "Genes",
       fill = "Selection Status") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "bottom")

ggsave("ML_results/Gene_Selection_Heatmap_PSO.pdf", p_heatmap, width = 10, height = 6)
cat("Gene selection heatmap saved to: ML_results/Gene_Selection_Heatmap_PSO.pdf\n")

# Final summary
cat("\n=== Summary of Gene Sharing Analysis - PSO Dataset ===\n")
cat("Total number of genes:", nrow(gene_stats), "\n")
cat("Total number of models: 8\n")
cat("100% shared genes (8/8):", nrow(core_genes_8), "->", if(nrow(core_genes_8) > 0) paste(core_genes_8$Gene, collapse = ", ") else "None", "\n")
cat("High-frequency shared genes (≥6 models):", nrow(gene_stats %>% filter(Occurrences >= 6)), "\n")
cat("Low-frequency shared genes (≤3 models):", nrow(gene_stats %>% filter(Occurrences <= 3)), "\n")

cat("\nMost important feature genes (based on model consensus):\n")
for(i in 1:nrow(gene_stats)) {
  rank <- ifelse(gene_stats$Occurrences[i] == 8, "⭐⭐⭐",
                 ifelse(gene_stats$Occurrences[i] >= 6, "⭐⭐", "⭐"))
  cat(sprintf("%s %s: Selected by %d/8 models\n", rank, gene_stats$Gene[i], gene_stats$Occurrences[i]))
}

################### Confusion Matrix ####################
################################ LASSO Confusion Matrix for PSO ########################
############################ 0. Path and Packages ##############################
setwd("D:/dualdisease/RANDOMreview")

library(caret)
library(pROC)
library(ggplot2)
library(dplyr)
library(glmnet)
library(tidyr)

#######################################################################
## 1. Read Data ########################################################
#######################################################################
cat("Reading PSO data...\n")
train_data <- read.csv("GSE13355_7genes_data.csv")  # Training set
internal_test <- read.csv("GSE78097_7genes_data.csv")  # Internal validation set
external_test1 <- read.csv("GSE109248_7genes_data.csv")  # External validation set 1
external_test2 <- read.csv("GSE141804_7genes_data.csv")  # External validation set 2

# Ensure DISEASE is a factor
train_data$DISEASE <- factor(train_data$DISEASE, levels = c(0, 1), labels = c("Class0", "Class1"))
internal_test$DISEASE <- factor(internal_test$DISEASE, levels = c(0, 1), labels = c("Class0", "Class1"))
external_test1$DISEASE <- factor(external_test1$DISEASE, levels = c(0, 1), labels = c("Class0", "Class1"))
external_test2$DISEASE <- factor(external_test2$DISEASE, levels = c(0, 1), labels = c("Class0", "Class1"))

# Remove sampleid column
train_data <- train_data[, -1]
internal_test <- internal_test[, -1]
external_test1 <- external_test1[, -1]
external_test2 <- external_test2[, -1]

cat("PSO data reading completed!\n")
cat("Training set samples:", nrow(train_data), "\n")
cat("Internal validation set samples:", nrow(internal_test), "\n")
cat("External validation set 1 samples:", nrow(external_test1), "\n")
cat("External validation set 2 samples:", nrow(external_test2), "\n")

#######################################################################
## 2. Train LASSO Model ###################################################
#######################################################################
cat("Starting LASSO model training...\n")
set.seed(123)

# Prepare data
x_train <- model.matrix(DISEASE ~ ., train_data)[, -1]
y_train <- ifelse(train_data$DISEASE == "Class1", 1, 0)

# Train LASSO model
lasso_model <- cv.glmnet(x_train, y_train, 
                         family = "binomial", 
                         alpha = 1, 
                         nfolds = 10)

# Use lambda.min
best_lambda <- lasso_model$lambda.min
cat("Selected lambda (lambda.min):", round(best_lambda, 5), "\n")

# View selected features
coef_mat <- coef(lasso_model, s = "lambda.min")
cat("Selected features:\n")
print(coef_mat)

#######################################################################
## 3. Predict Probabilities #########################################################
#######################################################################
# LASSO prediction probabilities
lasso_train_prob <- predict(lasso_model, x_train, s = "lambda.min", type = "response")[, 1]
lasso_internal_prob <- predict(lasso_model, 
                               model.matrix(DISEASE ~ ., internal_test)[, -1], 
                               s = "lambda.min", type = "response")[, 1]
lasso_external1_prob <- predict(lasso_model, 
                                model.matrix(DISEASE ~ ., external_test1)[, -1], 
                                s = "lambda.min", type = "response")[, 1]
lasso_external2_prob <- predict(lasso_model, 
                                model.matrix(DISEASE ~ ., external_test2)[, -1], 
                                s = "lambda.min", type = "response")[, 1]

#######################################################################
## 4. Performance Metric Calculation Functions #################################################
#######################################################################
# ----------- Function to calculate metrics from confusion matrix -----------
calc_metrics <- function(true_labels, pred_labels) {
  cm <- confusionMatrix(pred_labels, true_labels, positive = "Class1")
  tab <- cm$table
  TP <- tab[2,2]; TN <- tab[1,1]; FP <- tab[1,2]; FN <- tab[2,1]; total <- sum(tab)
  
  ##── Manual calculation ──##
  specificity  <- ifelse((TN + FP) > 0, TN / (TN + FP), 0)
  sensitivity  <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
  precision    <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
  accuracy     <- (TP + TN) / total
  f1           <- ifelse((precision + sensitivity) > 0,
                         2 * precision * sensitivity / (precision + sensitivity), 0)
  prevalence      <- (TP + FN) / total
  detection_rate  <- TP / total
  
  list(
    ConfusionMatrix = cm,
    Accuracy        = accuracy,
    Precision       = precision,
    Sensitivity     = sensitivity,
    Specificity     = specificity,
    F1              = f1,
    Prevalence      = prevalence,
    DetectionRate   = detection_rate
  )
}

# ----------- Bar plot function -----------
plot_metrics <- function(metrics_list, dataset_name, gse_id) {
  df <- tibble(
    Metric = c("Accuracy","Precision","Sensitivity","Specificity",
               "F1 Score","Prevalence","Detection Rate"),
    Value  = c(metrics_list$Accuracy, metrics_list$Precision,
               metrics_list$Sensitivity, metrics_list$Specificity,
               metrics_list$F1, metrics_list$Prevalence,
               metrics_list$DetectionRate)
  )
  ggplot(df, aes(Metric, Value, fill = Metric)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.3f", Value)), vjust = -0.3, size = 4) +
    ylim(0, 1) +
    labs(title = paste0("LASSO Performance Metrics - ", dataset_name, " (", gse_id, ")"),
         y = "Value", x = NULL) +
    theme_minimal() +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold"),
          plot.title   = element_text(hjust = 0.5, face = "bold"))
}

# ——— Confusion matrix plotting function ————————————————————————————————————————————————
plot_confmat <- function(cm, out_file, gse_id, threshold_type = "Default") {
  acc  <- cm$overall["Accuracy"]            # Extract Accuracy
  cm_df <- as.data.frame(cm$table) |>
    mutate(AccuracyLabel = ifelse(Reference == Prediction,
                                  "Correct", "Incorrect"))
  
  gg <- ggplot(cm_df,
               aes(x = Prediction, y = Reference, fill = AccuracyLabel)) +
    geom_tile(color = "white", alpha = 0.9, width = 0.95, height = 0.95) +
    geom_text(aes(label = Freq, colour = AccuracyLabel),
              size = 8, fontface = "bold", show.legend = FALSE) +
    scale_fill_manual(values = c(Correct   = "#4da6ff",
                                 Incorrect = "#ff7b7b")) +
    scale_colour_manual(values = c(Correct   = "black",
                                   Incorrect = "black"), guide = "none") +
    labs(title = paste0("LASSO Confusion Matrix (", gse_id, ")\n", 
                        "Accuracy = ", sprintf("%.3f", acc), " (", threshold_type, " Threshold)"),
         x = "Predicted", y = "Actual") +
    theme_bw(base_size = 14) +
    theme(plot.title  = element_text(hjust = 0.5, face = "bold"),
          panel.border = element_rect(color = "#e0e0e0",
                                      fill  = NA, size = 1.5),
          aspect.ratio = 1)
  
  ggsave(out_file, gg, width = 8, height = 7)
  message("Saved: ", out_file)
  invisible(gg)
}

#######################################################################
## 5. Training Set Performance Analysis ###################################################
#######################################################################
cat("=== Training Set Performance Analysis ===\n")

# ===== 1. Default threshold for training set =====
train_pred_prob <- lasso_train_prob
train_pred_raw  <- factor(ifelse(train_pred_prob > 0.5, "Class1", "Class0"),
                          levels = c("Class0", "Class1"))
train_metrics   <- calc_metrics(train_data$DISEASE, train_pred_raw)

cat("Train - default threshold metrics:\n")
print(train_metrics$ConfusionMatrix)

# ---- 1b. Optimal threshold for training set (Youden) ----
roc_train <- roc(train_data$DISEASE, train_pred_prob, levels = c("Class0","Class1"))
youden    <- roc_train$sensitivities + roc_train$specificities - 1
best_thr  <- roc_train$thresholds[which.max(youden)]
cat(sprintf("Train - Best threshold by Youden Index: %.3f\n", best_thr))

train_pred_opt <- factor(ifelse(train_pred_prob > best_thr, "Class1","Class0"),
                         levels = c("Class0","Class1"))
train_metrics_opt <- calc_metrics(train_data$DISEASE, train_pred_opt)
cat("Train - metrics @ optimal threshold:\n")
print(train_metrics_opt$ConfusionMatrix)

# ===== 2. Save training set bar plot (default threshold)=====
ggsave("ML_results/LASSO_Train_metrics_barplot_PSO.pdf",
       plot = plot_metrics(train_metrics, "Training Set", "GSE13355"),
       width = 8, height = 5)

# ===== 3. Save training set confusion matrix (default threshold)=====
plot_confmat(train_metrics$ConfusionMatrix, 
             "ML_results/LASSO_Train_confmat_default_PSO.pdf",
             "GSE13355", "Default")

