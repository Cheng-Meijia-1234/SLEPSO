setwd("D:/dualdisease/RANDOMreview")
library(readxl)
library(dplyr)
# GSE81622 data processing
gse81622_data <- read_excel("GSE81622_scale.xlsx")
genes <- c("IFI6", "MX1", "NMI", "OAS3", "OASL", "SAMD9", "UBE2L6")
gse81622_filtered <- gse81622_data %>%
  select(sampleid = 1, DISEASE = 2, all_of(genes))

write.csv(gse81622_filtered, "GSE81622_7genes_data.csv", row.names = FALSE)
cat("GSE81622 completed! Samples:", nrow(gse81622_filtered), "\n")

# GSE110174 data processing
gse110174_data <- read_excel("GSE110174_scale.xlsx")
gse110174_filtered <- gse110174_data %>%
  select(sampleid = 1, DISEASE = 2, all_of(genes))
write.csv(gse110174_filtered, "GSE110174_7genes_data.csv", row.names = FALSE)
cat("GSE110174 completed! Samples:", nrow(gse110174_filtered), "\n")

# Train data processing
train_data <- read_excel("SLEtrain_scale.xlsx")
train_filtered <- train_data %>%
  select(sampleid = 1, DISEASE = 2, all_of(genes))
write.csv(train_filtered, "SLEtrain_7genes_data.csv", row.names = FALSE)
cat("Train dataset completed! Samples:", nrow(train_filtered), "\n")

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

# Create results directory
if(!dir.exists("ML_results")) dir.create("ML_results")

# 1. Read data - Using your SLE dataset
cat("Reading SLE data...\n")
train_data <- read.csv("SLEtrain_7genes_data.csv")  # Training set
validation1 <- read.csv("GSE81622_7genes_data.csv")  # Validation set 1
validation2 <- read.csv("GSE110174_7genes_data.csv")  # Validation set 2

# Ensure DISEASE is a factor
train_data$DISEASE <- factor(train_data$DISEASE, levels = c(0, 1), labels = c("Control", "Disease"))
validation1$DISEASE <- factor(validation1$DISEASE, levels = c(0, 1), labels = c("Control", "Disease"))
validation2$DISEASE <- factor(validation2$DISEASE, levels = c(0, 1), labels = c("Control", "Disease"))

# Remove sampleid column
train_data <- train_data[, -1]
validation1 <- validation1[, -1]
validation2 <- validation2[, -1]

cat("SLE data reading completed!\n")
cat("Training set samples:", nrow(train_data), "\n")
cat("Validation set 1 samples:", nrow(validation1), "\n")
cat("Validation set 2 samples:", nrow(validation2), "\n")

# 2. Set training control parameters
ctrl <- trainControl(
  method = "cv",
  number = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = TRUE
)

# 3. Define ROC plotting function - Modified for 2 validation sets
plot_roc_comparison <- function(roc_results, method_name) {
  # Extract all ROC curves
  roc_list <- list(
    Train = roc_results$train_roc,
    Validation1 = roc_results$validation1_roc,
    Validation2 = roc_results$validation2_roc
  )
  
  # Create AUC table
  auc_tbl <- data.frame(
    Dataset = c("Train", "Validation1", "Validation2"),
    AUC = c(
      auc(roc_results$train_roc),
      auc(roc_results$validation1_roc),
      auc(roc_results$validation2_roc)
    )
  )
  
  # Plot ROC curves
  p <- ggroc(roc_list, legacy.axes = TRUE, size = 1.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
    ggtitle(paste(method_name, "ROC Curves - SLE Dataset")) +
    theme_bw() +
    labs(x = "False Positive Rate (1 - Specificity)", 
         y = "True Positive Rate (Sensitivity)") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          legend.position = "bottom")
  
  # Add AUC annotations
  annot_y <- c(0.25, 0.20, 0.15)
  for (i in 1:nrow(auc_tbl)) {
    p <- p + annotate("text", x = 0.6, y = annot_y[i],
                      label = sprintf("%s: AUC = %.3f", auc_tbl$Dataset[i], auc_tbl$AUC[i]),
                      hjust = 0, size = 4)
  }
  
  # Save plot
  filename <- paste0("ML_results/", method_name, "_ROC_SLE.pdf")
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

# RF predictions
rf_train_prob <- predict(rf_model, train_data, type = "prob")[, "Disease"]
rf_validation1_prob <- predict(rf_model, validation1, type = "prob")[, "Disease"]
rf_validation2_prob <- predict(rf_model, validation2, type = "prob")[, "Disease"]

rf_results <- list(
  train_roc = roc(train_data$DISEASE, rf_train_prob),
  validation1_roc = roc(validation1$DISEASE, rf_validation1_prob),
  validation2_roc = roc(validation2$DISEASE, rf_validation2_prob)
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

# Use lambda.min
best_lambda <- lasso_model$lambda.min
cat("Selected lambda (lambda.min):", round(best_lambda, 5), "\n")

# LASSO predictions
lasso_train_prob <- predict(lasso_model, x_train, s = "lambda.min", type = "response")[, 1]
lasso_validation1_prob <- predict(lasso_model, 
                                  model.matrix(DISEASE ~ ., validation1)[, -1], 
                                  s = "lambda.min", type = "response")[, 1]
lasso_validation2_prob <- predict(lasso_model, 
                                  model.matrix(DISEASE ~ ., validation2)[, -1], 
                                  s = "lambda.min", type = "response")[, 1]

lasso_results <- list(
  train_roc = roc(y_train, lasso_train_prob),
  validation1_roc = roc(validation1$DISEASE, lasso_validation1_prob),
  validation2_roc = roc(validation2$DISEASE, lasso_validation2_prob)
)

lasso_plot <- plot_roc_comparison(lasso_results, "LASSO")

cat("LASSO completed!\n")

# View detailed coefficient information of LASSO model
cat("=== LASSO Model Coefficient Analysis ===\n")

# Extract coefficients
coef_mat <- coef(lasso_model, s = "lambda.min")
cat("Complete coefficient matrix:\n")
print(coef_mat)

# Convert to data frame for easier viewing
coef_df <- as.data.frame(as.matrix(coef_mat))
colnames(coef_df) <- "Coefficient"
coef_df$Feature <- rownames(coef_df)
coef_df <- coef_df[, c("Feature", "Coefficient")]

cat("\nCoefficient details:\n")
print(coef_df)

# Filter non-zero coefficients (selected features)
non_zero_coef <- coef_df[coef_df$Coefficient != 0, ]
cat("\n=== Selected Features (Non-zero Coefficients) ===\n")
print(non_zero_coef)

# Sort by absolute coefficient values
sorted_coef <- non_zero_coef[order(abs(non_zero_coef$Coefficient), decreasing = TRUE), ]
cat("\n=== Features Sorted by Importance ===\n")
print(sorted_coef)

# Visualize coefficients
library(ggplot2)

# Remove intercept term, only show gene coefficients
gene_coef <- sorted_coef[sorted_coef$Feature != "(Intercept)", ]

if(nrow(gene_coef) > 0) {
  p_coef <- ggplot(gene_coef, aes(x = reorder(Feature, Coefficient), y = Coefficient, fill = Coefficient > 0)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = "LASSO Model Gene Coefficients",
         x = "Genes",
         y = "Coefficient Value") +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                      name = "Direction",
                      labels = c("Negative", "Positive")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  print(p_coef)
  ggsave("ML_results/LASSO_Gene_Coefficients.pdf", p_coef, width = 10, height = 6)
}

# Model performance summary
cat("\n=== LASSO Model Performance Summary ===\n")
cat(sprintf("Selected lambda value: %.5f\n", best_lambda))
cat(sprintf("Total number of features: %d\n", ncol(x_train)))
cat(sprintf("Number of selected features: %d\n", nrow(non_zero_coef) - 1))  # Subtract intercept term
cat(sprintf("Training set AUC: %.3f\n", auc(lasso_results$train_roc)))
cat(sprintf("Validation set 1 AUC: %.3f\n", auc(lasso_results$validation1_roc)))
cat(sprintf("Validation set 2 AUC: %.3f\n", auc(lasso_results$validation2_roc)))

# Explain each gene in detail
cat("\n=== Gene Coefficient Interpretation ===\n")
cat("Positive coefficient: Increases disease risk\n")
cat("Negative coefficient: Decreases disease risk\n")
cat("Larger absolute coefficient value: Greater impact on prediction\n")

# If there are selected genes, show detailed information
if(nrow(gene_coef) > 0) {
  cat("\nSpecific gene effects:\n")
  for(i in 1:nrow(gene_coef)) {
    gene <- gene_coef$Feature[i]
    coef_val <- gene_coef$Coefficient[i]
    direction <- ifelse(coef_val > 0, "increases disease risk", "decreases disease risk")
    cat(sprintf("%s: Coefficient=%.4f (%s)\n", gene, coef_val, direction))
  }
}

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

# XGBoost predictions
xgb_train_prob <- predict(xgb_model, train_data, type = "prob")[, "Disease"]
xgb_validation1_prob <- predict(xgb_model, validation1, type = "prob")[, "Disease"]
xgb_validation2_prob <- predict(xgb_model, validation2, type = "prob")[, "Disease"]

xgb_results <- list(
  train_roc = roc(train_data$DISEASE, xgb_train_prob),
  validation1_roc = roc(validation1$DISEASE, xgb_validation1_prob),
  validation2_roc = roc(validation2$DISEASE, xgb_validation2_prob)
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

# SVM predictions
svm_train_prob <- predict(svm_model, train_data, type = "prob")[, "Disease"]
svm_validation1_prob <- predict(svm_model, validation1, type = "prob")[, "Disease"]
svm_validation2_prob <- predict(svm_model, validation2, type = "prob")[, "Disease"]

svm_results <- list(
  train_roc = roc(train_data$DISEASE, svm_train_prob),
  validation1_roc = roc(validation1$DISEASE, svm_validation1_prob),
  validation2_roc = roc(validation2$DISEASE, svm_validation2_prob)
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

# AdaBoost predictions
ada_train_prob <- predict(ada_model, train_data, type = "prob")[, "Disease"]
ada_validation1_prob <- predict(ada_model, validation1, type = "prob")[, "Disease"]
ada_validation2_prob <- predict(ada_model, validation2, type = "prob")[, "Disease"]

ada_results <- list(
  train_roc = roc(train_data$DISEASE, ada_train_prob),
  validation1_roc = roc(validation1$DISEASE, ada_validation1_prob),
  validation2_roc = roc(validation2$DISEASE, ada_validation2_prob)
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

# GBM predictions
gbm_train_prob <- predict(gbm_model, train_data, type = "prob")[, "Disease"]
gbm_validation1_prob <- predict(gbm_model, validation1, type = "prob")[, "Disease"]
gbm_validation2_prob <- predict(gbm_model, validation2, type = "prob")[, "Disease"]

gbm_results <- list(
  train_roc = roc(train_data$DISEASE, gbm_train_prob),
  validation1_roc = roc(validation1$DISEASE, gbm_validation1_prob),
  validation2_roc = roc(validation2$DISEASE, gbm_validation2_prob)
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

# Naive Bayes predictions
nb_train_prob <- predict(nb_model, train_data, type = "prob")[, "Disease"]
nb_validation1_prob <- predict(nb_model, validation1, type = "prob")[, "Disease"]
nb_validation2_prob <- predict(nb_model, validation2, type = "prob")[, "Disease"]

nb_results <- list(
  train_roc = roc(train_data$DISEASE, nb_train_prob),
  validation1_roc = roc(validation1$DISEASE, nb_validation1_prob),
  validation2_roc = roc(validation2$DISEASE, nb_validation2_prob)
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

# KKNN predictions
kknn_train_prob <- predict(kknn_model, train_data, type = "prob")[, "Disease"]
kknn_validation1_prob <- predict(kknn_model, validation1, type = "prob")[, "Disease"]
kknn_validation2_prob <- predict(kknn_model, validation2, type = "prob")[, "Disease"]

kknn_results <- list(
  train_roc = roc(train_data$DISEASE, kknn_train_prob),
  validation1_roc = roc(validation1$DISEASE, kknn_validation1_prob),
  validation2_roc = roc(validation2$DISEASE, kknn_validation2_prob)
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
write.csv(all_auc, "ML_results/All_Methods_AUC_Summary_SLE.csv", row.names = FALSE)

# Create method comparison plot
comparison_plot <- ggplot(all_auc, aes(x = Dataset, y = AUC, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = round(AUC, 3)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3) +
  scale_fill_brewer(palette = "Set3") +
  ggtitle("AUC Comparison of 8 Machine Learning Methods - SLE Dataset") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("ML_results/Methods_Comparison_SLE.pdf", comparison_plot, width = 10, height = 7)

cat("All SLE data analysis completed!\n")
cat("Results are saved in the ML_results folder:\n")
cat("- ROC curves for each method\n")
cat("- AUC summary table: All_Methods_AUC_Summary_SLE.csv\n")
cat("- Method comparison plot: Methods_Comparison_SLE.pdf\n")

# Print best methods
best_methods <- all_auc %>%
  group_by(Dataset) %>%
  slice_max(AUC, n = 1)

cat("\nBest methods for each dataset:\n")
print(best_methods)

# 13. LASSO Feature Importance Lollipop Plot
cat("Generating LASSO feature importance lollipop plot...\n")

# Extract LASSO coefficients
coef_mat <- coef(lasso_model, s = "lambda.min")
coef_df <- as.data.frame(as.matrix(coef_mat))
colnames(coef_df) <- "coef"

# Remove intercept term
coef_df <- coef_df[rownames(coef_df) != "(Intercept)", , drop = FALSE]

# Create data frame
lasso.result <- data.frame(
  diffvariable = rownames(coef_df),
  coef = coef_df$coef
)

# Calculate absolute coefficients
lasso.result$abs_coef <- abs(lasso.result$coef)

# Sort by absolute value and take top 5
lasso.plot.data <- lasso.result[order(lasso.result$abs_coef, decreasing = TRUE), ]
lasso.plot.data <- head(lasso.plot.data, 5)

# Macaron color palette
macaron_colors <- c("#FFB3BA", "#FFDFBA", "#FFFFBA", "#BAFFC9", "#BAE1FF", 
                    "#FFB347", "#FFD700", "#FF6F61", "#FFB5E8", "#B5EAD7")

# Plot LASSO lollipop plot
p_lasso_macaron <- ggplot(lasso.plot.data, 
                          aes(y = reorder(diffvariable, abs_coef), x = abs_coef)) + 
  geom_segment(aes(yend = reorder(diffvariable, abs_coef), x = 0, xend = abs_coef, 
                   color = diffvariable), 
               size = 1.8, alpha = 0.8) + 
  geom_point(aes(color = diffvariable), size = 9, alpha = 0.8) +  
  geom_text(aes(label = round(abs_coef, 3)), color = "black", size = 3.2, fontface = "bold") +  
  scale_color_manual(values = macaron_colors) + 
  labs(title = "LASSO Feature Importance - SLE Dataset", 
       subtitle = paste("Top", nrow(lasso.plot.data), "Features by Absolute Coefficient Value"),
       y = "Genes", 
       x = "Absolute Coefficient Value") + 
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = "none") +
  expand_limits(x = 0)

# Save plot
ggsave("ML_results/LASSO_Macaron_Absolute_SLE.pdf", p_lasso_macaron, width = 6, height = 5)
cat("LASSO lollipop plot saved!\n")

# 14. UpSet Plot Analysis
cat("Starting UpSet plot analysis...\n")

# Define feature extraction functions
extract_rf_top5 <- function(model) {
  imp <- varImp(model)$importance
  imp_df <- data.frame(
    Feature = rownames(imp),
    Score = imp$Control
  )
  imp_df <- imp_df[order(imp_df$Score, decreasing = TRUE), ]
  head(imp_df$Feature, 5)
}

extract_lasso_top5 <- function(model) {
  coef_mat <- coef(model, s = "lambda.min")
  features <- rownames(coef_mat)[coef_mat[, 1] != 0]
  features <- setdiff(features, "(Intercept)")
  head(features, 5)
}

extract_xgb_top5 <- function(model) {
  imp <- varImp(model)$importance
  imp_df <- data.frame(
    Feature = rownames(imp),
    Score = imp$Overall
  )
  imp_df <- imp_df[order(imp_df$Score, decreasing = TRUE), ]
  head(imp_df$Feature, 5)
}

extract_svm_top5 <- function(model) {
  imp <- varImp(model)$importance
  imp_df <- data.frame(
    Feature = rownames(imp),
    Score = imp$Control
  )
  imp_df <- imp_df[order(imp_df$Score, decreasing = TRUE), ]
  head(imp_df$Feature, 5)
}

extract_ada_top5 <- function(model) {
  imp <- varImp(model)$importance
  imp_df <- data.frame(
    Feature = rownames(imp),
    Score = imp$Control
  )
  imp_df <- imp_df[order(imp_df$Score, decreasing = TRUE), ]
  head(imp_df$Feature, 5)
}

extract_gbm_top5 <- function(model) {
  imp <- varImp(model)$importance
  imp_df <- data.frame(
    Feature = rownames(imp),
    Score = imp$Overall
  )
  imp_df <- imp_df[order(imp_df$Score, decreasing = TRUE), ]
  head(imp_df$Feature, 5)
}

extract_nb_top5 <- function(model) {
  imp <- varImp(model)$importance
  imp_df <- data.frame(
    Feature = rownames(imp),
    Score = imp$Control
  )
  imp_df <- imp_df[order(imp_df$Score, decreasing = TRUE), ]
  head(imp_df$Feature, 5)
}

extract_kknn_top5 <- function(model) {
  imp <- varImp(model)$importance
  imp_df <- data.frame(
    Feature = rownames(imp),
    Score = imp$Control
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

# Generate UpSet plot
library(ComplexHeatmap)
library(grid)

m_comb <- make_comb_mat(feature_sets_top5)

# Macaron color palette
macaron_colors <- c("#FFB3BA", "#FFDFBA", "#FFFFBA", "#BAFFC9", 
                    "#BAE1FF", "#FFB347", "#FFD700", "#FF6F61")
names(macaron_colors) <- names(feature_sets_top5)

# Plot UpSet plot
pdf("ML_results/UpSet_Top5_Features_SLE.pdf", width = 10, height = 7)
up <- UpSet(
  m_comb,
  set_order = names(feature_sets_top5),
  bg_col = macaron_colors,
  comb_order = order(comb_size(m_comb), decreasing = TRUE),
  row_names_gp = gpar(fontsize = 10),
  column_title = "Top5 Features Intersection Across 8 Models - SLE"
)

ht <- draw(up)
cs <- comb_size(m_comb)
od <- column_order(ht)
decorate_annotation("intersection_size", {
  grid.text(
    cs[od],
    x = seq_along(cs),
    y = unit(cs[od], "native") + unit(2, "pt"),
    just = "bottom",
    gp = gpar(fontsize = 9)
  )
})
dev.off()
cat("UpSet plot saved to: ML_results/UpSet_Top5_Features_SLE.pdf\n")

# Count gene occurrence frequency
gene_model_df <- data.frame(
  Gene = unlist(feature_sets_top5),
  Model = rep(names(feature_sets_top5), each = 5)
)

gene_stats <- gene_model_df %>%
  group_by(Gene) %>%
  summarise(
    Occurrences = n(),
    Models = paste(unique(Model), collapse = ", ")
  ) %>%
  arrange(desc(Occurrences))

cat("Occurrence statistics of Top5 genes across models:\n")
print(gene_stats, row.names = FALSE)

# 15. Chord Diagram Analysis
cat("Starting chord diagram analysis...\n")
library(circlize)

# Build adjacency matrix for shared features between models
methods <- names(feature_sets_top5)
k <- length(methods)
adj_mat <- matrix(0, nrow = k, ncol = k, dimnames = list(methods, methods))

# Calculate number of shared features for each pair of models
for (i in 1:(k-1)) {
  for (j in (i+1):k) {
    shared <- intersect(feature_sets_top5[[i]], feature_sets_top5[[j]])
    adj_mat[i, j] <- adj_mat[j, i] <- length(shared)
  }
}

# Plot chord diagram
pdf("ML_results/Chord_Top5_Features_SLE.pdf", width = 9, height = 9)
circos.clear()

chordDiagram(
  adj_mat,
  grid.col = macaron_colors,
  transparency = 0.2,
  directional = FALSE,
  annotationTrack = c("grid", "name"),
  annotationTrackHeight = c(0.03, 0.08),
  preAllocateTracks = list(track.height = 0.1),
  link.border = "white",
  link.lwd = 2,
  link.sort = TRUE,
  link.decreasing = TRUE
)

title("Shared Top5 Genes Between Machine Learning Models - SLE", cex = 1.2)
dev.off()

cat("Chord diagram saved to: ML_results/Chord_Top5_Features_SLE.pdf\n")

# 16. Save detailed LASSO coefficient results
lasso_coef_details <- lasso.result[order(lasso.result$abs_coef, decreasing = TRUE), ]
write.csv(lasso_coef_details, "ML_results/LASSO_Coefficient_Details_SLE.csv", row.names = FALSE)
cat("Detailed LASSO coefficient results saved to: ML_results/LASSO_Coefficient_Details_SLE.csv\n")

# 17. Print important statistics
cat("\n=== SLE Dataset Important Statistics ===\n")
cat("Training set samples:", nrow(train_data), "\n")
cat("Validation set 1 samples:", nrow(validation1), "\n") 
cat("Validation set 2 samples:", nrow(validation2), "\n")
cat("Number of features selected by LASSO:", nrow(lasso.result), "\n")
cat("Core genes appearing in all models:\n")
print(filter(gene_stats, Occurrences >= 4))

cat("\nSLE data analysis fully completed!\n")
cat("Generated files include:\n")
cat("- ROC curves for each method: *_ROC_SLE.pdf\n")
cat("- Method comparison plot: Methods_Comparison_SLE.pdf\n")
cat("- AUC summary table: All_Methods_AUC_Summary_SLE.csv\n")
cat("- LASSO lollipop plot: LASSO_Macaron_Absolute_SLE.pdf\n")
cat("- LASSO coefficient details: LASSO_Coefficient_Details_SLE.csv\n")
cat("- UpSet plot: UpSet_Top5_Features_SLE.pdf\n")
cat("- Chord diagram: Chord_Top5_Features_SLE.pdf\n")

############UPSET Plot Interpretation###############
# Modify statistics code to ensure full model names are displayed
cat("Occurrence statistics of Top5 genes across models (Complete Information):\n")

# Recalculate and display complete information
gene_stats_detailed <- gene_model_df %>%
  group_by(Gene) %>%
  summarise(
    Occurrences = n(),
    Models = paste(Model, collapse = ", ")  # Do not deduplicate, show all models
  ) %>%
  arrange(desc(Occurrences))

# Print complete information without truncation
print(gene_stats_detailed, row.names = FALSE, width = Inf)

# Additionally display complete Top5 list for each model separately
cat("\n=== Detailed Top5 Genes List for Each Model ===\n")
for (model_name in names(feature_sets_top5)) {
  cat(sprintf("\n%s:\n", model_name))
  cat(paste(feature_sets_top5[[model_name]], collapse = ", "), "\n")
}

# Display most important shared genes
cat("\n=== Core Shared Genes Analysis ===\n")
cat("Genes appearing in all 8 models (100% shared):\n")
core_genes_8 <- gene_stats_detailed %>% filter(Occurrences == 8)
print(core_genes_8, row.names = FALSE, width = Inf)

cat("\nGenes appearing in 7 models:\n")
core_genes_7 <- gene_stats_detailed %>% filter(Occurrences == 7)
print(core_genes_7, row.names = FALSE, width = Inf)

cat("\nGenes appearing in 6 models:\n")
core_genes_6 <- gene_stats_detailed %>% filter(Occurrences == 6)
print(core_genes_6, row.names = FALSE, width = Inf)

# Calculate specific model distribution for each gene
cat("\n=== Detailed Model Distribution for Each Gene ===\n")
for(i in 1:nrow(gene_stats_detailed)) {
  gene <- gene_stats_detailed$Gene[i]
  models <- gene_stats_detailed$Models[i]
  cat(sprintf("%s: Appears in %d models -> %s\n", 
              gene, gene_stats_detailed$Occurrences[i], models))
}

# Save complete statistics to file
write.csv(gene_stats_detailed, "ML_results/Gene_Occurrence_Statistics_SLE.csv", row.names = FALSE)
cat("\nComplete statistics saved to: ML_results/Gene_Occurrence_Statistics_SLE.csv\n")

# Create visualization to show gene sharing
library(ggplot2)

# Create gene sharing frequency bar plot
p_gene_frequency <- ggplot(gene_stats_detailed, aes(x = reorder(Gene, Occurrences), y = Occurrences, fill = Gene)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Occurrences), vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = macaron_colors) +
  labs(title = "Gene Occurrence Frequency Across 8 Machine Learning Models - SLE Dataset",
       x = "Genes", 
       y = "Number of Models") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "none") +
  ylim(0, 8.5)

ggsave("ML_results/Gene_Frequency_Barplot_SLE.pdf", p_gene_frequency, width = 10, height = 6)
cat("Gene frequency bar plot saved to: ML_results/Gene_Frequency_Barplot_SLE.pdf\n")

# Create heatmap showing gene occurrence in each model
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
  labs(title = "Top5 Gene Selection Heatmap Across Models - SLE Dataset",
       x = "Machine Learning Models", 
       y = "Genes",
       fill = "Selection Status") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "bottom")

ggsave("ML_results/Gene_Selection_Heatmap_SLE.pdf", p_heatmap, width = 10, height = 6)
cat("Gene selection heatmap saved to: ML_results/Gene_Selection_Heatmap_SLE.pdf\n")

# Final summary
cat("\n=== SLE Dataset Gene Sharing Analysis Summary ===\n")
cat("Total number of genes: 7\n")
cat("Total number of models: 8\n")
cat("100% shared genes (8/8):", nrow(core_genes_8), "->", paste(core_genes_8$Gene, collapse = ", "), "\n")
cat("High-frequency shared genes (≥6 models):", nrow(gene_stats_detailed %>% filter(Occurrences >= 6)), "\n")
cat("Low-frequency shared genes (≤3 models):", nrow(gene_stats_detailed %>% filter(Occurrences <= 3)), "\n")

cat("\nMost important feature genes (based on model consensus):\n")
for(i in 1:nrow(gene_stats_detailed)) {
  rank <- ifelse(gene_stats_detailed$Occurrences[i] == 8, "⭐⭐⭐",
                 ifelse(gene_stats_detailed$Occurrences[i] >= 6, "⭐⭐", "⭐"))
  cat(sprintf("%s %s: %d/8 models selected\n", rank, gene_stats_detailed$Gene[i], gene_stats_detailed$Occurrences[i]))
}

###################Confusion Matrix - SLE Dataset LASSO Model####################
################################ LASSO Confusion Matrix for SLE ########################
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
cat("Reading SLE data...\n")
train_data <- read.csv("SLEtrain_7genes_data.csv")  # Training set
validation1 <- read.csv("GSE81622_7genes_data.csv")  # Validation set 1
validation2 <- read.csv("GSE110174_7genes_data.csv")  # Validation set 2

# Ensure DISEASE is a factor (using Class0 and Class1 labels)
train_data$DISEASE <- factor(train_data$DISEASE, levels = c(0, 1), labels = c("Class0", "Class1"))
validation1$DISEASE <- factor(validation1$DISEASE, levels = c(0, 1), labels = c("Class0", "Class1"))
validation2$DISEASE <- factor(validation2$DISEASE, levels = c(0, 1), labels = c("Class0", "Class1"))

# Remove sampleid column
train_data <- train_data[, -1]
validation1 <- validation1[, -1]
validation2 <- validation2[, -1]

cat("SLE data reading completed!\n")
cat("Training set samples:", nrow(train_data), "\n")
cat("Validation set 1 samples:", nrow(validation1), "\n")
cat("Validation set 2 samples:", nrow(validation2), "\n")

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
lasso_validation1_prob <- predict(lasso_model, 
                                  model.matrix(DISEASE ~ ., validation1)[, -1], 
                                  s = "lambda.min", type = "response")[, 1]
lasso_validation2_prob <- predict(lasso_model, 
                                  model.matrix(DISEASE ~ ., validation2)[, -1], 
                                  s = "lambda.min", type = "response")[, 1]

#######################################################################
## 4. Performance Metrics Calculation Function #################################################
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
ggsave("ML_results/LASSO_Train_metrics_barplot_SLE.pdf",
       plot = plot_metrics(train_metrics, "Training Set", "SLEtrain"),
       width = 8, height = 5)

# ===== 3. Save training set confusion matrix (default threshold)=====
plot_confmat(train_metrics$ConfusionMatrix, 
             "ML_results/LASSO_Train_confmat_default_SLE.pdf",
             "SLEtrain", "Default")

#######################################################################
## 6. Validation Set Performance Analysis ###################################################
#######################################################################
# Define validation set list (including GSE IDs)
test_datasets <- list(
  Validation1 = list(data = validation1, prob = lasso_validation1_prob, name = "Validation Set 1", gse_id = "GSE81622"),
  Validation2 = list(data = validation2, prob = lasso_validation2_prob, name = "Validation Set 2", gse_id = "GSE110174")
)

for (test_set in test_datasets) {
  nm <- test_set$name
  gse_id <- test_set$gse_id
  testData <- test_set$data
  prob <- test_set$prob
  
  cat(sprintf("\n=== %s Performance Analysis (%s) ===\n", nm, gse_id))
  
  # ---- Default threshold ----
  pred_raw <- factor(ifelse(prob > 0.5, "Class1", "Class0"),
                     levels = c("Class0", "Class1"))
  met_default <- calc_metrics(testData$DISEASE, pred_raw)
  cat(sprintf("%s - default threshold metrics:\n", nm))
  print(met_default$ConfusionMatrix)
  
  # ---- Youden optimal threshold ----
  roc_te  <- roc(testData$DISEASE, prob, levels = c("Class0","Class1"))
  youden  <- roc_te$sensitivities + roc_te$specificities - 1
  thr_opt <- roc_te$thresholds[which.max(youden)]
  cat(sprintf("%s - Best threshold by Youden Index: %.3f\n", nm, thr_opt))
  
  pred_opt <- factor(ifelse(prob > thr_opt, "Class1","Class0"),
                     levels = c("Class0","Class1"))
  met_opt  <- calc_metrics(testData$DISEASE, pred_opt)
  cat(sprintf("%s - metrics @ optimal threshold:\n", nm))
  print(met_opt$ConfusionMatrix)
  
  # ---- Save default threshold bar plot ----
  ggsave(paste0("ML_results/LASSO_", gsub(" ", "_", nm), "_metrics_barplot_SLE.pdf"),
         plot = plot_metrics(met_default, nm, gse_id),
         width = 8, height = 5)
  
  # ---- Save default threshold confusion matrix ----
  plot_confmat(met_default$ConfusionMatrix,
               paste0("ML_results/LASSO_", gsub(" ", "_", nm), "_confmat_default_SLE.pdf"),
               gse_id, "Default")
  
  # ---- Save optimal threshold confusion matrix ----
  plot_confmat(met_opt$ConfusionMatrix,
               paste0("ML_results/LASSO_", gsub(" ", "_", nm), "_confmat_optimal_SLE.pdf"),
               gse_id, "Optimal")
}

#######################################################################
## 7. Summarize Performance Metrics #####################################################
#######################################################################
cat("\n=== LASSO Model Performance Summary ===\n")

# Collect performance metrics for all datasets
all_metrics <- list()

# Training set
all_metrics$SLEtrain <- calc_metrics(train_data$DISEASE, 
                                     factor(ifelse(lasso_train_prob > 0.5, "Class1", "Class0"),
                                            levels = c("Class0", "Class1")))

# Validation sets
all_metrics$GSE81622 <- calc_metrics(validation1$DISEASE,
                                     factor(ifelse(lasso_validation1_prob > 0.5, "Class1", "Class0"),
                                            levels = c("Class0", "Class1")))

all_metrics$GSE110174 <- calc_metrics(validation2$DISEASE,
                                      factor(ifelse(lasso_validation2_prob > 0.5, "Class1", "Class0"),
                                             levels = c("Class0", "Class1")))

# Create performance summary table
performance_summary <- data.frame(
  Dataset = names(all_metrics),
  Accuracy = sapply(all_metrics, function(x) x$Accuracy),
  Sensitivity = sapply(all_metrics, function(x) x$Sensitivity),
  Specificity = sapply(all_metrics, function(x) x$Specificity),
  Precision = sapply(all_metrics, function(x) x$Precision),
  F1_Score = sapply(all_metrics, function(x) x$F1),
  AUC = c(auc(roc(train_data$DISEASE, lasso_train_prob)),
          auc(roc(validation1$DISEASE, lasso_validation1_prob)),
          auc(roc(validation2$DISEASE, lasso_validation2_prob)))
)

# Print performance summary
cat("\nLASSO model performance summary across datasets (default threshold):\n")
print(performance_summary)

# Save performance summary table
write.csv(performance_summary, "ML_results/LASSO_Performance_Summary_SLE.csv", row.names = FALSE)
cat("\nPerformance summary table saved to: ML_results/LASSO_Performance_Summary_SLE.csv\n")

# Create performance comparison plot
performance_long <- performance_summary %>%
  pivot_longer(cols = c(Accuracy, Sensitivity, Specificity, Precision, F1_Score, AUC),
               names_to = "Metric", values_to = "Value")

p_perf_comparison <- ggplot(performance_long, aes(x = Dataset, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = sprintf("%.3f", Value)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3) +
  scale_fill_brewer(palette = "Set3") +
  ggtitle("LASSO Model Performance Metrics Comparison - SLE Dataset\n(Based on 7 Gene Features)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 1)

ggsave("ML_results/LASSO_Performance_Comparison_SLE.pdf", p_perf_comparison, width = 12, height = 7)

cat("LASSO model performance comparison plot saved to: ML_results/LASSO_Performance_Comparison_SLE.pdf\n")

#######################################################################
## 8. Final Summary #########################################################
#######################################################################
cat("\n=== LASSO Confusion Matrix Analysis Completed - SLE Dataset ===\n")
cat("Generated files include:\n")
cat("- Training set performance metrics bar plot: LASSO_Train_metrics_barplot_SLE.pdf (SLEtrain)\n")
cat("- Training set confusion matrix plot: LASSO_Train_confmat_default_SLE.pdf (SLEtrain)\n")
cat("- Validation set 1 performance metrics bar plot: LASSO_Validation_Set_1_metrics_barplot_SLE.pdf (GSE81622)\n")
cat("- Validation set 1 confusion matrix plots: LASSO_Validation_Set_1_confmat_*_SLE.pdf (GSE81622)\n")
cat("- Validation set 2 performance metrics bar plot: LASSO_Validation_Set_2_metrics_barplot_SLE.pdf (GSE110174)\n")
cat("- Validation set 2 confusion matrix plots: LASSO_Validation_Set_2_confmat_*_SLE.pdf (GSE110174)\n")
cat("- Performance summary table: LASSO_Performance_Summary_SLE.csv\n")
cat("- Performance comparison plot: LASSO_Performance_Comparison_SLE.pdf\n")

cat("\nLASSO model performance summary on SLE dataset:\n")
cat(sprintf("Training set (SLEtrain) Accuracy: %.3f\n", performance_summary$Accuracy[1]))
cat(sprintf("Validation set 1 (GSE81622) Accuracy: %.3f\n", performance_summary$Accuracy[2]))
cat(sprintf("Validation set 2 (GSE110174) Accuracy: %.3f\n", performance_summary$Accuracy[3]))

cat("\nAll analysis results have been saved to the ML_results folder!\n")

#################### LASSO Model SHAP Analysis - SLE Dataset ##################
######## 0. Load Packages & Data Preparation ########
setwd("D:/dualdisease/RANDOMreview")

library(fastshap)
library(glmnet)
library(shapviz)
library(dplyr)
library(ggplot2)

# Create output directory
if (!dir.exists("ML_results")) {
  dir.create("ML_results")
}

# Read data
cat("Reading SLE data...\n")
train_data <- read.csv("SLEtrain_7genes_data.csv")

# Save sampleid and process DISEASE column
sample_ids <- train_data$sampleid
train_data$DISEASE <- factor(train_data$DISEASE, levels = c(0, 1), labels = c("Class0", "Class1"))

# Remove sampleid column for modeling
train_data_model <- train_data[, -1]
cat("SLE training set samples:", nrow(train_data_model), "\n")

#######################################################################
## 1. Train LASSO Model ##################################################
#######################################################################
cat("Starting LASSO model training...\n")
set.seed(123)

x_train <- model.matrix(DISEASE ~ ., train_data_model)[, -1]
y_train <- ifelse(train_data_model$DISEASE == "Class1", 1, 0)

lasso_model <- cv.glmnet(x_train, y_train, 
                         family = "binomial", 
                         alpha = 1, 
                         nfolds = 10)

best_lambda <- lasso_model$lambda.min
cat("Selected lambda (lambda.min):", round(best_lambda, 5), "\n")

coef_mat <- coef(lasso_model, s = "lambda.min")
cat("Selected features:\n")
print(coef_mat)

################ 2. LASSO Prediction Wrapper #######################
pfun_lasso <- function(object, newdata) {
  if (!is.matrix(newdata)) {
    newdata_matrix <- model.matrix(~ . - 1, data = newdata)
  } else {
    newdata_matrix <- newdata
  }
  
  colnames(newdata_matrix) <- colnames(x_train)
  preds <- predict(object, newx = newdata_matrix, s = "lambda.min", type = "response")
  return(as.numeric(preds))
}

################ 3. Prepare SHAP Analysis Data ###################
feature_cols <- setdiff(colnames(train_data_model), "DISEASE")
X_train <- as.data.frame(x_train)
colnames(X_train) <- feature_cols

baseline <- mean(pfun_lasso(lasso_model, X_train))
cat("Baseline prediction probability:", round(baseline, 4), "\n")

################ 4. Unified SHAP Plotting Function ###################
# Unified SHAP plotting function to ensure sample consistency
plot_shap_analysis <- function(sample_id, sample_data, nsim = 500, 
                               save_prefix = "Sample", description = "") {
  
  # Find row index by sample_id
  sample_idx <- which(sample_ids == sample_id)
  if (length(sample_idx) == 0) {
    stop(paste("Sample not found:", sample_id))
  }
  
  # Get sample data
  sample_row <- sample_data[sample_idx, , drop = FALSE]
  pred_prob <- pfun_lasso(lasso_model, sample_row)
  true_label <- as.character(train_data_model$DISEASE[sample_idx])
  
  cat(sprintf("\nAnalyzing sample: %s\n", sample_id))
  cat("True label:", true_label, "\n")
  cat("Predicted probability:", round(pred_prob, 4), "\n")
  
  # Calculate SHAP values
  set.seed(123)
  shap_values <- explain(
    lasso_model, X = X_train, pred_wrapper = pfun_lasso,
    newdata = sample_row, nsim = nsim, adjust = TRUE
  )
  
  # Create shapviz object
  sv_obj <- shapviz(shap_values, X = sample_row, baseline = baseline)
  
  # Custom plotting function
  custom_plot <- function(plot_func, plot_type) {
    plot_func(sv_obj) +
      labs(title = paste0("LASSO SHAP ", plot_type, " Plot\n", 
                          "Sample: ", sample_id, 
                          " | Pred: ", round(pred_prob, 3),
                          " | True: ", true_label,
                          ifelse(description != "", paste0("\n", description), ""))) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
  }
  
  # Plot waterfall plot
  p_waterfall <- custom_plot(sv_waterfall, "Waterfall")
  print(p_waterfall)
  ggsave(paste0("ML_results/SHAP_Waterfall_", save_prefix, "_", sample_id, "_SLE.pdf"), 
         p_waterfall, width = 8, height = 6)
  
  # Plot force plot
  p_force <- custom_plot(sv_force, "Force")
  print(p_force)
  ggsave(paste0("ML_results/SHAP_Force_", save_prefix, "_", sample_id, "_SLE.pdf"), 
         p_force, width = 10, height = 6)
  
  # Return sample information
  return(list(
    sample_id = sample_id,
    sample_idx = sample_idx,
    pred_prob = pred_prob,
    true_label = true_label,
    shap_values = shap_values,
    shap_range = range(shap_values)
  ))
}

################ 5. Single Sample Analysis (Using Sample ID instead of Index) ###################
# Select first sample for analysis
first_sample_id <- sample_ids[1]
sample1_info <- plot_shap_analysis(
  sample_id = first_sample_id,
  sample_data = X_train,
  nsim = 500,
  save_prefix = "Sample1",
  description = "First Sample Analysis"
)

################ 6. Global SHAP Analysis ########################
cat("\nStarting global SHAP analysis...\n")
set.seed(123)

n_samples_global <- min(100, nrow(X_train))
global_indices <- sample(1:nrow(X_train), n_samples_global)

global_expl <- explain(
  lasso_model, X = X_train[global_indices, ], pred_wrapper = pfun_lasso,
  nsim = 100, adjust = TRUE, shap_only = FALSE
)

sg <- shapviz(global_expl, X = X_train[global_indices, ])

# Feature importance plot
p_importance <- sv_importance(sg) +
  labs(title = "LASSO Global Feature Importance - SLE Dataset") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
print(p_importance)
ggsave("ML_results/SHAP_Global_Importance_SLE.pdf", p_importance, width = 8, height = 6)

# Beeswarm plot
p_beeswarm <- sv_importance(sg, kind = "beeswarm") +
  labs(title = "LASSO SHAP Beeswarm Plot - SLE Dataset") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
print(p_beeswarm)
ggsave("ML_results/SHAP_Beeswarm_SLE.pdf", p_beeswarm, width = 10, height = 7)

################ 7. Find "Bidirectional" Samples ###################
cat("\nStarting search for 'bidirectional' samples...\n")

find_bidirectional_sample <- function(pos_thr = 0.02, neg_thr = -0.02, nsim_search = 200) {
  
  # Calculate prediction probabilities for all samples
  all_probs <- pfun_lasso(lasso_model, X_train)
  prob_tbl <- data.frame(
    row_id = 1:nrow(X_train),
    sample_id = sample_ids,
    pred_prob = all_probs,
    true_label = train_data_model$DISEASE
  )
  
  # Prioritize samples near decision boundary
  borderline_samples <- prob_tbl %>% 
    filter(pred_prob >= 0.3 & pred_prob <= 0.7) %>%
    arrange(desc(abs(pred_prob - 0.5)))
  
  cat(sprintf("Checking %d borderline samples...\n", nrow(borderline_samples)))
  
  # Check borderline samples
  for (i in 1:min(20, nrow(borderline_samples))) {
    sample_id_candidate <- borderline_samples$sample_id[i]
    sample_row <- X_train[borderline_samples$row_id[i], , drop = FALSE]
    
    shap_test <- explain(lasso_model, X = X_train, pred_wrapper = pfun_lasso,
                         newdata = sample_row, nsim = nsim_search, adjust = TRUE)
    
    if (any(shap_test > pos_thr) && any(shap_test < neg_thr)) {
      return(sample_id_candidate)
    }
  }
  
  # Expand search range
  cat("No bidirectional samples found among borderline samples, expanding search...\n")
  random_indices <- sample(1:nrow(X_train), min(50, nrow(X_train)))
  
  for (i in random_indices) {
    sample_id_candidate <- sample_ids[i]
    sample_row <- X_train[i, , drop = FALSE]
    
    shap_test <- explain(lasso_model, X = X_train, pred_wrapper = pfun_lasso,
                         newdata = sample_row, nsim = nsim_search, adjust = TRUE)
    
    if (any(shap_test > pos_thr) && any(shap_test < neg_thr)) {
      return(sample_id_candidate)
    }
  }
  
  return(NA)
}

# Find bidirectional sample
bidir_sample_id <- find_bidirectional_sample()

if (!is.na(bidir_sample_id)) {
  cat(sprintf("Found 'bidirectional' sample: %s\n", bidir_sample_id))
  
  # Calculate with high precision and plot
  bidir_info <- plot_shap_analysis(
    sample_id = bidir_sample_id,
    sample_data = X_train,
    nsim = 1000,
    save_prefix = "Bidirectional",
    description = "Bidirectional SHAP Values"
  )
  
  cat("\n=== Detailed Information for 'Bidirectional' Sample ===\n")
  cat("Sample ID:", bidir_info$sample_id, "\n")
  cat("Predicted probability:", round(bidir_info$pred_prob, 4), "\n")
  cat("True label:", bidir_info$true_label, "\n")
  cat("SHAP value range:", round(bidir_info$shap_range, 4), "\n")
  
} else {
  warning("No samples meeting 'bidirectional' criteria found, try adjusting thresholds.")
}

################ 8. Batch Analysis of Risk Samples ###################
analyze_risk_samples <- function(risk_type = "high", threshold = 0.8, max_samples = 3) {
  
  all_probs <- pfun_lasso(lasso_model, X_train)
  
  if (risk_type == "high") {
    risk_samples <- data.frame(
      sample_id = sample_ids,
      pred_prob = all_probs,
      true_label = train_data_model$DISEASE
    ) %>% 
      filter(pred_prob >= threshold) %>%
      arrange(desc(pred_prob))
    
    description <- "High Risk"
    save_prefix <- "HighRisk"
  } else {
    risk_samples <- data.frame(
      sample_id = sample_ids,
      pred_prob = all_probs,
      true_label = train_data_model$DISEASE