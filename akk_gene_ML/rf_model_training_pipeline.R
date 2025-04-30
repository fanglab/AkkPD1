# rf_model_training_pipeline.R with parallel support
# Multi-mode RF Model: ML1, ML2, ML3

library(caret)
library(readr)
library(dplyr)
library(MLeval)
library(optparse)
library(pROC)
library(doParallel)

# ---------------------------
# 1. CLI options
# ---------------------------
option_list <- list(
  make_option(c("--input"), type = "character", help = "Path to full gene table CSV"),
  make_option(c("--mode"), type = "character", help = "Model type: ML1, ML2, or ML3"),
  make_option(c("--threads"), type = "integer", default = 4, help = "Number of parallel threads [default: %default]"),
  make_option(c("--test"), action = "store_true", default = FALSE, help = "Run in quick test mode [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---------------------------
# 2. Load and preprocess data
# ---------------------------
raw_data <- read_csv(opt$input)

metadata_row_ids <- grep("^(Cohort|Response|Gender|Antibiotics|LIGNE|Clean_Reads_Number|Age|Phylo|High_Akk|ECOGPS|BMI|Smoking)", raw_data[[1]])

metadata <- as.data.frame(raw_data[metadata_row_ids, -1])
rownames(metadata) <- raw_data[[1]][metadata_row_ids]
gene_data <- as.data.frame(raw_data[-metadata_row_ids, ])
rownames(gene_data) <- gene_data[[1]]
gene_data <- gene_data[, -1]

metadata_t <- as.data.frame(t(metadata))
gene_data_t <- as.data.frame(t(gene_data))
metadata_t$SampleID <- rownames(metadata_t)
gene_data_t$SampleID <- rownames(gene_data_t)

data <- merge(metadata_t, gene_data_t, by = "SampleID")
data$SampleID <- NULL

# Convert metadata
factor_cols <- c("Response", "Antibiotics", "Gender", "Cohort", "Phylo")
data[factor_cols] <- lapply(data[factor_cols], factor)
data$Clean_Reads_Number <- as.numeric(data$Clean_Reads_Number)
data$Age <- as.numeric(data$Age)

# ---------------------------
# 3. Training/Testing Split
# ---------------------------
gene_start_col <- which(colnames(data) == "Gene_1")
data[, gene_start_col:ncol(data)] <- lapply(data[, gene_start_col:ncol(data)], factor)

if ("High_Akk_in_C3 (Yes or No)" %in% colnames(data)) {
  data$HighAkk <- as.logical(data$`High_Akk_in_C3 (Yes or No)` == "Yes")
} else {
  stop("Missing 'High_Akk_in_C3 (Yes or No)' column.")
}

if (opt$mode == "ML1") {
  Training <- data %>% filter((Cohort %in% c("C1", "C2")) | (Cohort == "C3" & !HighAkk))
  Testing  <- data %>% filter(Cohort == "C4" & `LIGNE (NA, 1, or >1)` != 1)
} else if (opt$mode == "ML2") {
  Training <- data %>% filter((Cohort %in% c("C1", "C2")) | (Cohort == "C3" & !HighAkk & `LIGNE (NA, 1, or >1)` != 1))
  Testing  <- data %>% filter(Cohort == "C4" & `LIGNE (NA, 1, or >1)` != 1)
} else if (opt$mode == "ML3") {
  Training <- data
  Testing <- NULL
} else {
  stop("Invalid --mode argument: must be ML1, ML2, or ML3")
}

# ---------------------------
# 4. Quick Test Mode
# ---------------------------
if (opt$test) {
  cat("Running in quick test mode: using pre-made importance table...\n")

  # Only one tuning length in test
  tune_lengths <- c(50)

  # Load the pre-made importance table
  importance_table <- read_csv("/data/Importance_score.csv")  # Adjust path if needed

  top_50_genes <- importance_table %>%
    arrange(desc(importance)) %>%
    slice(1:50) %>%
    pull(gene)

  cat("Top 50 precomputed important genes:\n")
  print(top_50_genes)

  # Subset Training to top important genes
  gene_cols <- intersect(colnames(Training), top_50_genes)
  if (length(gene_cols) == 0) {
    stop("No matching top genes found in input data.")
  }

  # Parallel backend
  cl <- makeCluster(opt$threads)
  registerDoParallel(cl)

  ctrl <- trainControl(method = "LOOCV", classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = "final")

  model <- train(Response ~ ., data = Training[, c("Response", gene_cols)],
                 method = "rf", trControl = ctrl, tuneLength = 50)

  stopCluster(cl)

  # Save model
  saveRDS(model, paste0("fit_", opt$mode, "_TL50_test.rds"))

  # Generate ROC curve from LOOCV predictions
  # Get best model tuning parameters
  best_params <- model$bestTune
  
  # Subset predictions matching the bestTune
  filtered_pred <- model$pred
  
  for (param_name in names(best_params)) {
    filtered_pred <- filtered_pred[filtered_pred[[param_name]] == best_params[[param_name]], ]
  }
  
  # Now get predictions and true labels
  predicted_probs <- filtered_pred$R
  true_labels <- filtered_pred$obs
  
  # Now make ROC
  roc_obj <- pROC::roc(response = true_labels, predictor = predicted_probs)
  auc_value <- pROC::auc(roc_obj)
  
  # Save ROC plot
  png(paste0("demo_ROC_curve_", opt$mode, ".png"), width = 800, height = 600)
  plot.roc(roc_obj, main = paste("Demo ROC Curve (AUC =", round(auc_value, 3), ")"),
           col = "#1c61b6", lwd = 3)
  dev.off()


  cat("\n==== Quick Test Run Finished Successfully ====\n")
  quit(save = "no")
}

# ---------------------------
# 5. Full Model Training (no --test)
# ---------------------------
tune_lengths <- c(50, 100, 150, 200, 500)

gene_cols <- colnames(Training)[gene_start_col:ncol(Training)]

cl <- makeCluster(opt$threads)
registerDoParallel(cl)

ctrl <- trainControl(method = "LOOCV", classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = "final")

models <- list()
results <- list()

for (tl in tune_lengths) {
  cat(paste0("Training with tuneLength = ", tl, "...\n"))

  model <- train(Response ~ ., data = Training[, c("Response", gene_cols)],
                 method = "rf", trControl = ctrl, tuneLength = tl)

  models[[paste0("TL", tl)]] <- model
  saveRDS(model, paste0("fit_", opt$mode, "_TL", tl, ".rds"))

  if (opt$mode == "ML3") {
    results[[paste0("TL", tl)]] <- evalm(model)
  }
}

stopCluster(cl)

# ---------------------------
# 6. Full Model Evaluation
# ---------------------------
if (opt$mode == "ML3") {
  roc_values <- sapply(results, function(res) max(res$roc$AUC))
  best_tl <- names(which.max(roc_values))
  best_model <- models[[best_tl]]
  best_auc <- roc_values[best_tl]
} else {
  best_tl <- names(models)[1]
  best_model <- models[[best_tl]]

  if (!is.null(Testing)) {
    preds <- predict(best_model, Testing[, gene_cols])
    conf_matrix <- confusionMatrix(preds, Testing$Response)
    print(conf_matrix)

    probs <- predict(best_model, Testing[, gene_cols], type = "prob")[, 2]
    roc <- pROC::roc(response = Testing$Response, predictor = probs)
    best_auc <- as.numeric(pROC::auc(roc))
  } else {
    stop("Testing set is not available for ML1/ML2 evaluation.")
  }
}

# ---------------------------
# 7. Save Full Results
# ---------------------------
importance <- varImp(best_model)$importance
write.csv(importance, paste0("Importance_score_", opt$mode, "_", best_tl, ".csv"))

summary_df <- data.frame(
  Mode = opt$mode,
  Best_TuneLength = best_tl,
  Best_AUC = best_auc
)
write.csv(summary_df, paste0("model_summary_", opt$mode, ".csv"), row.names = FALSE)
print(summary_df)

cat("\n==== Full Training Run Finished Successfully ====\n")
