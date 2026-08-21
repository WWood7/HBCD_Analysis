source("Scripts/00_config.r")
source("Scripts/Helpers/03_birthweight_helper.r")

suppressPackageStartupMessages({
  library(sl3)
  library(dplyr)
  library(readr)
})

# Number of participant-level outer cross-fitting folds.
num_outer_folds <- 2L
runtime_preset <- "full"
set.seed(123)

pipeline_start <- Sys.time()
message("Pipeline started at: ", format(pipeline_start, "%Y-%m-%d %H:%M:%S"))

# Load the participant-level HBCD analysis table.
input_path <- file.path(preprocessed_dir, "preprocessed_df.csv")
if (!file.exists(input_path)) {
  stop("Preprocessed HBCD data not found: ", input_path)
}
df_preprocessed <- readr::read_csv(input_path, show_col_types = FALSE)

# Define classification and regression Super Learners.
if (runtime_preset == "fast") {
  cls_learners <- list(
    Lrnr_glmnet$new(nfolds = 5, nlambda = 50),
    Lrnr_ranger$new(probability = TRUE, num.trees = 300, min.node.size = 20)
  )
  reg_learners <- list(
    Lrnr_glmnet$new(nfolds = 5, nlambda = 50),
    Lrnr_ranger$new(num.trees = 300, min.node.size = 20)
  )
} else {
  cls_learners <- list(
    Lrnr_glm$new(),
    Lrnr_ranger$new(probability = TRUE),
    Lrnr_earth$new(),
    Lrnr_mean$new(),
    Lrnr_gbm$new()
  )
  reg_learners <- list(
    Lrnr_glm$new(),
    Lrnr_ranger$new(),
    Lrnr_earth$new(),
    Lrnr_mean$new(),
    Lrnr_gbm$new()
  )
}

cls_stack <- do.call(Stack$new, cls_learners)
reg_stack <- do.call(Stack$new, reg_learners)
learners <- list(
  cls_lrnr = Lrnr_sl$new(
    learners = cls_stack,
    metalearner = Lrnr_nnls$new(eval_function = loss_loglik_binomial)
  ),
  reg_lrnr = Lrnr_sl$new(
    learners = reg_stack,
    metalearner = Lrnr_nnls$new(eval_function = loss_squared_error)
  )
)

message("\n========== HBCD birthweight mediation: prenatal cannabis 0 vs 1 ==========")
results <- run_birthweight_pipeline(
  df_input = df_preprocessed,
  a_var = A_var,
  x_vars = X_vars,
  m_vars = M_vars,
  y_vars = Y_vars,
  learners = learners,
  num_outer_folds = num_outer_folds
)

models_output_dir <- file.path(birthweight_output_dir, "Models")
dir.create(models_output_dir, recursive = TRUE, showWarnings = FALSE)
out_rds <- file.path(models_output_dir, "birthweight_mediation_results.rds")
out_csv <- file.path(models_output_dir, "birthweight_mediation_results.csv")
saveRDS(results, file = out_rds)
readr::write_csv(results, file = out_csv)

pipeline_end <- Sys.time()
message(
  "\nSaved results to ", out_rds, " and ", out_csv,
  "\nPipeline finished at: ", format(pipeline_end, "%Y-%m-%d %H:%M:%S"),
  " | Total runtime: ",
  round(as.numeric(difftime(pipeline_end, pipeline_start, units = "mins")), 2),
  " mins"
)