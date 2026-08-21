# Variables for the HBCD birthweight mediation analysis
Y_vars <- c("birth_weight_lbs")
A_var <- "prenatal_cannabis"
M_vars <- c("gestational_age")
X_vars <- c("prenatal_alcohol", "prenatal_nicotine", "prenatal_opioid", "child_race", "child_ethnicity",
            "mother_education", "household_income", "site")


##### preprocess the dataset to enable sl3 compatibility:
## drop rows with NA in X/A/M/Y
## add Y_sq to the dataset
## one-hot encode X and M
preprocess_data <- function(
    data,
    y_vars = Y_vars,
    a_var = A_var,
    x_vars = X_vars,
    m_vars = M_vars
) {
  
  required_vars <- unique(c(a_var, x_vars, m_vars, y_vars))
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars)) {
    stop("Missing required variables: ", paste(missing_vars, collapse = ", "))
  }

  # This analysis uses the complete-case population selected by the user.
  out <- data[stats::complete.cases(data[, required_vars, drop = FALSE]), , drop = FALSE]
  if (!nrow(out)) {
    stop("No complete cases remain after filtering the required variables.")
  }

  sanitize_sl3_col <- function(x) {
    if (inherits(x, "labelled")) x <- as.numeric(x)
    else if (is.ordered(x)) x <- factor(x)
    else if (is.factor(x)) x <- factor(x)
    else if (is.character(x)) x <- as.character(x)
    if (length(class(x)) > 1) class(x) <- class(x)[1]
    x
  }

  for (v in required_vars) {
    out[[v]] <- sanitize_sl3_col(out[[v]])
  }
  for (v in c(m_vars, y_vars)) {
    if (!is.numeric(out[[v]])) {
      stop(v, " must be numeric.")
    }
    # sl3 treats integer-valued outcomes as categorical unless converted.
    out[[v]] <- as.numeric(out[[v]])
  }
  
  for (yv in y_vars) {
    out[[paste0(yv, "_sq")]] <- out[[yv]]^2
  }

  # one-hot encode X and M separately so their updated name vectors are explicit
  X_mat <- stats::model.matrix(
    ~ . - 1,
    data = out[, x_vars, drop = FALSE]
  )
  M_mat <- stats::model.matrix(
    ~ . - 1,
    data = out[, m_vars, drop = FALSE]
  )
  colnames(X_mat) <- make.names(colnames(X_mat), unique = TRUE)
  colnames(M_mat) <- make.names(colnames(M_mat), unique = TRUE)
  x_vars_num <- colnames(X_mat)
  m_vars_num <- colnames(M_mat)

  # Keep non-covariate columns as-is, replace X/M with encoded columns
  out_encoded <- cbind(
    out[, setdiff(names(out), c(x_vars, m_vars)), drop = FALSE],
    as.data.frame(X_mat),
    as.data.frame(M_mat)
  )

  list(
    data = out_encoded,
    X_vars = x_vars_num,
    M_vars = m_vars_num
  )
}


# Run preprocessing and participant-level causal estimation.
run_birthweight_pipeline <- function(
    df_input,
    a_var = A_var,
    x_vars = X_vars,
    learners,
    y_vars = Y_vars,
    m_vars = M_vars,
    num_outer_folds = 2L
) {
  num_outer_folds <- as.integer(num_outer_folds)
  prep <- preprocess_data(
    df_input,
    a_var = a_var,
    x_vars = x_vars,
    y_vars = y_vars,
    m_vars = m_vars
  )
  df_model <- prep$data
  x_vars_num <- prep$X_vars
  m_vars_num <- prep$M_vars

  treatment_values <- sort(unique(df_model[[a_var]]))
  if (!identical(as.numeric(treatment_values), c(0, 1))) {
    stop(a_var, " must contain exactly the values 0 and 1 after preprocessing.")
  }

  folds <- create_participant_folds(
    data = df_model,
    num_folds = num_outer_folds,
    a_var = a_var
  )

  results <- estimate_causal_estimands(
    data = df_model,
    X_vars = x_vars_num,
    M_vars = m_vars_num,
    Y_vars = y_vars,
    A_var = a_var,
    learners = learners,
    folds = folds
  )
  dplyr::mutate(
    results,
    treatment = a_var,
    mediator = paste(m_vars, collapse = ", "),
    n_complete = nrow(df_model),
    n_excluded = nrow(df_input) - nrow(df_model),
    num_outer_folds = num_outer_folds
  )
}


##### Create treatment-balanced participant-level cross-fitting folds
create_participant_folds <- function(
    data,
    num_folds = 2L,
    a_var = A_var,
    seed = 123L
) {
  if (!a_var %in% names(data)) {
    stop("Treatment variable not found: ", a_var)
  }
  num_folds <- as.integer(num_folds)
  if (length(num_folds) != 1L || is.na(num_folds) || num_folds < 2L) {
    stop("num_folds must be one integer greater than or equal to 2.")
  }
  treatment_values <- sort(unique(data[[a_var]]))
  if (!identical(as.numeric(treatment_values), c(0, 1))) {
    stop(a_var, " must contain exactly the values 0 and 1.")
  }
  treatment_counts <- table(data[[a_var]])
  if (any(treatment_counts < num_folds)) {
    stop("Each treatment group needs at least num_folds observations.")
  }

  set.seed(seed)
  row_fold <- integer(nrow(data))
  for (a_level in treatment_values) {
    indices <- which(data[[a_var]] == a_level)
    indices <- sample(indices, length(indices), replace = FALSE)
    row_fold[indices] <- rep(seq_len(num_folds), length.out = length(indices))
  }

  outer_folds <- lapply(seq_len(num_folds), function(fold_id) {
    validation_set <- which(row_fold == fold_id)
    training_set <- which(row_fold != fold_id)
    if (!length(training_set) || !length(validation_set)) {
      stop("An outer cross-fitting fold has an empty training or validation set.")
    }
    if (!all(c(0, 1) %in% unique(data[[a_var]][training_set]))) {
      stop("An outer training fold does not contain both treatment levels.")
    }
    origami::make_fold(
      v = fold_id,
      training_set = training_set,
      validation_set = validation_set
    )
  })

  validation_rows <- sort(unlist(lapply(outer_folds, `[[`, "validation_set")))
  if (!identical(validation_rows, seq_len(nrow(data)))) {
    stop("Each row must occur in exactly one outer validation fold.")
  }
  outer_folds
}


# IID variance of an asymptotically linear estimator
calc_iid_var <- function(eif_vector) {
  ok <- is.finite(eif_vector)
  eif <- eif_vector[ok]
  if (length(eif) < 2L) {
    return(NA_real_)
  }
  stats::var(eif) / length(eif)
}

# define a function to estimate the causal estimands
estimate_causal_estimands <- function(data, X_vars, M_vars, Y_vars, A_var, learners, folds) {
  reg_lrnr <- learners$reg_lrnr
  cls_lrnr <- learners$cls_lrnr
  results <- data.frame()

  make_inner_folds <- function(fold_data, max_folds = 5L) {
    treatment_counts <- table(fold_data[[A_var]])
    if (length(treatment_counts) != 2L || min(treatment_counts) < 2L) {
      stop("Each outer training sample needs at least two observations per treatment level.")
    }
    create_participant_folds(
      data = fold_data,
      num_folds = min(max_folds, min(treatment_counts)),
      a_var = A_var,
      seed = 123L
    )
  }

  assert_binary_support <- function(x, label, fold_id) {
    if (!all(c(0, 1) %in% unique(x))) {
      stop("Outer training fold ", fold_id, " does not contain both levels of ", label, ".")
    }
  }

  assert_oof_predictions <- function(column_names) {
    incomplete <- column_names[vapply(column_names, function(column_name) {
      anyNA(data[[column_name]]) || any(!is.finite(data[[column_name]]))
    }, logical(1))]
    if (length(incomplete)) {
      stop("Missing or non-finite out-of-fold predictions: ", paste(incomplete, collapse = ", "))
    }
  }

  # Shared nuisance functions are fit on each outer training fold and evaluated
  # only on participants in the corresponding validation fold.
  shared_nuisance <- c("ps", "ps_m")
  for (column_name in shared_nuisance) data[[column_name]] <- NA_real_

  for (fold_id in seq_along(folds)) {
    train_idx <- folds[[fold_id]]$training_set
    valid_idx <- folds[[fold_id]]$validation_set
    train_data <- data[train_idx, , drop = FALSE]
    valid_data <- data[valid_idx, , drop = FALSE]
    valid_data_0 <- valid_data
    valid_data_0[[A_var]] <- 0
    valid_data_1 <- valid_data
    valid_data_1[[A_var]] <- 1

    assert_binary_support(train_data[[A_var]], A_var, fold_id)
    inner_folds <- make_inner_folds(train_data)

    task_ps_train <- sl3::make_sl3_Task(
      data = train_data,
      outcome = A_var,
      outcome_type = "binomial",
      covariates = X_vars,
      folds = inner_folds
    )
    ps_fit <- cls_lrnr$train(task_ps_train)
    task_ps_valid <- sl3::make_sl3_Task(data = valid_data, covariates = X_vars)
    data$ps[valid_idx] <- ps_fit$predict(task_ps_valid)

    task_ps_m_train <- sl3::make_sl3_Task(
      data = train_data,
      outcome = A_var,
      outcome_type = "binomial",
      covariates = c(X_vars, M_vars),
      folds = inner_folds
    )
    ps_m_fit <- cls_lrnr$train(task_ps_m_train)
    task_ps_m_valid <- sl3::make_sl3_Task(
      data = valid_data,
      covariates = c(X_vars, M_vars)
    )
    data$ps_m[valid_idx] <- ps_m_fit$predict(task_ps_m_valid)
  }
  assert_oof_predictions(shared_nuisance)

  # truncate the ps's to avoid numerical instability
  bound <- 0.025
  data$ps <- pmax(pmin(data$ps, 1 - bound), bound)
  data$ps_m <- pmax(pmin(data$ps_m, 1 - bound), bound)

  data$dens_ratio <- (1 - data$ps_m) / data$ps_m * data$ps / (1 - data$ps)



  ####### for the parts that depend on Y_vars
  for (v in Y_vars) {
    sq_flag <- 0
    for (yv in c(v, paste0(v, "_sq"))) {
      prediction_columns <- paste0(
        yv,
        c("_or_0", "_or_1", "_sr_1_0", "_sr_0_0", "_sr_1_1")
      )
      for (column_name in prediction_columns) data[[column_name]] <- NA_real_

      for (fold_id in seq_along(folds)) {
        train_idx <- folds[[fold_id]]$training_set
        valid_idx <- folds[[fold_id]]$validation_set
        train_data <- data[train_idx, , drop = FALSE]
        valid_data <- data[valid_idx, , drop = FALSE]
        assert_binary_support(train_data[[A_var]], A_var, fold_id)

        train_data_0 <- train_data
        train_data_0[[A_var]] <- 0
        train_data_1 <- train_data
        train_data_1[[A_var]] <- 1
        valid_data_0 <- valid_data
        valid_data_0[[A_var]] <- 0
        valid_data_1 <- valid_data
        valid_data_1[[A_var]] <- 1

        outcome_inner_folds <- make_inner_folds(train_data)
        task_outcome_train <- sl3::make_sl3_Task(
          data = train_data,
          outcome = yv,
          outcome_type = "continuous",
          covariates = c(X_vars, M_vars, A_var),
          folds = outcome_inner_folds
        )
        outcome_fit <- reg_lrnr$train(task_outcome_train)

        task_outcome_train_0 <- sl3::make_sl3_Task(
          data = train_data_0,
          covariates = c(X_vars, M_vars, A_var)
        )
        task_outcome_train_1 <- sl3::make_sl3_Task(
          data = train_data_1,
          covariates = c(X_vars, M_vars, A_var)
        )
        task_outcome_valid_0 <- sl3::make_sl3_Task(
          data = valid_data_0,
          covariates = c(X_vars, M_vars, A_var)
        )
        task_outcome_valid_1 <- sl3::make_sl3_Task(
          data = valid_data_1,
          covariates = c(X_vars, M_vars, A_var)
        )

        train_data[[paste0(yv, "_or_0")]] <- outcome_fit$predict(task_outcome_train_0)
        train_data[[paste0(yv, "_or_1")]] <- outcome_fit$predict(task_outcome_train_1)
        data[[paste0(yv, "_or_0")]][valid_idx] <- outcome_fit$predict(task_outcome_valid_0)
        data[[paste0(yv, "_or_1")]][valid_idx] <- outcome_fit$predict(task_outcome_valid_1)

        sequential_inner_folds <- make_inner_folds(train_data)
        task_sequential_1_0 <- sl3::make_sl3_Task(
          data = train_data,
          outcome = paste0(yv, "_or_1"),
          outcome_type = "continuous",
          covariates = c(X_vars, A_var),
          folds = sequential_inner_folds
        )
        task_sequential_0_0 <- sl3::make_sl3_Task(
          data = train_data,
          outcome = paste0(yv, "_or_0"),
          outcome_type = "continuous",
          covariates = c(X_vars, A_var),
          folds = sequential_inner_folds
        )
        task_sequential_1_1 <- sl3::make_sl3_Task(
          data = train_data,
          outcome = paste0(yv, "_or_1"),
          outcome_type = "continuous",
          covariates = c(X_vars, A_var),
          folds = sequential_inner_folds
        )

        sequential_fit_1_0 <- reg_lrnr$train(task_sequential_1_0)
        sequential_fit_0_0 <- reg_lrnr$train(task_sequential_0_0)
        sequential_fit_1_1 <- reg_lrnr$train(task_sequential_1_1)
        task_sequential_valid_0 <- sl3::make_sl3_Task(
          data = valid_data_0,
          covariates = c(X_vars, A_var)
        )
        task_sequential_valid_1 <- sl3::make_sl3_Task(
          data = valid_data_1,
          covariates = c(X_vars, A_var)
        )

        data[[paste0(yv, "_sr_1_0")]][valid_idx] <-
          sequential_fit_1_0$predict(task_sequential_valid_0)
        data[[paste0(yv, "_sr_0_0")]][valid_idx] <-
          sequential_fit_0_0$predict(task_sequential_valid_0)
        data[[paste0(yv, "_sr_1_1")]][valid_idx] <-
          sequential_fit_1_1$predict(task_sequential_valid_1)
      }
      assert_oof_predictions(prediction_columns)

      ## calculate the efficient influence function
      theta_0_0 <- mean(
        (data[[A_var]] == 0) / (1 - data$ps) * data[[paste0(yv, "_or_0")]]
      )
      theta_1_1 <- mean(
        (data[[A_var]] == 1) / data$ps * data[[paste0(yv, "_or_1")]]
      )
      theta_1_0 <- mean(
        (data[[A_var]] == 0) / (1 - data$ps) * data[[paste0(yv, "_or_1")]]
      )
      data[[paste0(yv, "_D_0_0")]] <- (data[[A_var]] == 0) / (1 - data$ps) *
        (data[[yv]] - data[[paste0(yv, "_or_0")]]) +
        (data[[A_var]] == 0) / (1 - data$ps) * (data[[paste0(yv, "_or_0")]] - data[[paste0(yv, "_sr_0_0")]]) +
        data[[paste0(yv, "_sr_0_0")]] - theta_0_0

      data[[paste0(yv, "_D_1_1")]] <- (data[[A_var]] == 1) / data$ps *
        (data[[yv]] - data[[paste0(yv, "_or_1")]]) +
        (data[[A_var]] == 1) / data$ps * (data[[paste0(yv, "_or_1")]] - data[[paste0(yv, "_sr_1_1")]]) +
        data[[paste0(yv, "_sr_1_1")]] - theta_1_1

      data[[paste0(yv, "_D_1_0")]] <- (data[[A_var]] == 1) / data$ps *
        data$dens_ratio * (data[[yv]] - data[[paste0(yv, "_or_1")]]) +
        (data[[A_var]] == 0) / (1 - data$ps) * (data[[paste0(yv, "_or_1")]] - data[[paste0(yv, "_sr_1_0")]]) +
        data[[paste0(yv, "_sr_1_0")]] - theta_1_0
      
      if (sq_flag == 0) {
        t_0_0 <- theta_0_0 + mean(data[[paste0(yv, "_D_0_0")]])
        t_1_1 <- theta_1_1 + mean(data[[paste0(yv, "_D_1_1")]])
        t_1_0 <- theta_1_0 + mean(data[[paste0(yv, "_D_1_0")]])
      } else if (sq_flag == 1) {
        t_0_0_sq <- theta_0_0 + mean(data[[paste0(yv, "_D_0_0")]])
        t_1_1_sq <- theta_1_1 + mean(data[[paste0(yv, "_D_1_1")]])
        t_1_0_sq <- theta_1_0 + mean(data[[paste0(yv, "_D_1_0")]])
      }
      sq_flag <- sq_flag + 1
    }
    
    # the actual theta's regarding to Y, not squared Y
    ate <- t_1_1 - t_0_0 
    nde <- t_1_0 - t_0_0
    nie <- t_1_1 - t_1_0
    var_floor <- 1e-6
    # Use the ATE counterfactual outcome variance as the universal
    # standardizer for all three causal effect sizes.
    v_ate <- 0.5 * pmax(t_1_1_sq - t_1_1^2, var_floor) + 0.5 * pmax(t_0_0_sq - t_0_0^2, var_floor)
    ate_es <- ate / sqrt(v_ate)
    nde_es <- nde / sqrt(v_ate)
    nie_es <- nie / sqrt(v_ate)
    # IID influence-function variances for participant-level data.
    ate_var <- calc_iid_var(data[[paste0(v, "_D_1_1")]] - data[[paste0(v, "_D_0_0")]])
    nde_var <- calc_iid_var(data[[paste0(v, "_D_1_0")]] - data[[paste0(v, "_D_0_0")]])
    nie_var <- calc_iid_var(data[[paste0(v, "_D_1_1")]] - data[[paste0(v, "_D_1_0")]])
    # calculate the EIF for the effect sizes
    # note the v vs yv distinction here, when this part is run, it is expected that
    # v is the variable name for the outcome variable, yv is the variable name for the outcome variable squared
    data[[paste0(v, "_D_ate")]] <- (data[[paste0(v, "_D_1_1")]] - data[[paste0(v, "_D_0_0")]]) / sqrt(v_ate) -
      ate_es / (2 * v_ate) *
      (0.5 * (data[[paste0(yv, "_D_1_1")]] - 2 * t_1_1 * data[[paste0(v, "_D_1_1")]]) +
      0.5 * (data[[paste0(yv, "_D_0_0")]] - 2 * t_0_0 * data[[paste0(v, "_D_0_0")]]))
    
    data[[paste0(v, "_D_nde")]] <- (data[[paste0(v, "_D_1_0")]] - data[[paste0(v, "_D_0_0")]]) / sqrt(v_ate) -
      nde_es / (2 * v_ate) *
      (0.5 * (data[[paste0(yv, "_D_1_1")]] - 2 * t_1_1 * data[[paste0(v, "_D_1_1")]]) +
      0.5 * (data[[paste0(yv, "_D_0_0")]] - 2 * t_0_0 * data[[paste0(v, "_D_0_0")]]))
    
    data[[paste0(v, "_D_nie")]] <- (data[[paste0(v, "_D_1_1")]] - data[[paste0(v, "_D_1_0")]]) / sqrt(v_ate) -
      nie_es / (2 * v_ate) *
      (0.5 * (data[[paste0(yv, "_D_1_1")]] - 2 * t_1_1 * data[[paste0(v, "_D_1_1")]]) +
      0.5 * (data[[paste0(yv, "_D_0_0")]] - 2 * t_0_0 * data[[paste0(v, "_D_0_0")]]))
    
    ate_es_var <- calc_iid_var(data[[paste0(v, "_D_ate")]])
    nde_es_var <- calc_iid_var(data[[paste0(v, "_D_nde")]])
    nie_es_var <- calc_iid_var(data[[paste0(v, "_D_nie")]])

    # build a dataframe to store the results
    temp_results <- data.frame(
      y_var = v,
      ate = ate,
      ate_var = ate_var,
      nde = nde,
      nde_var = nde_var,
      nie = nie,
      nie_var = nie_var,
      ate_es = ate_es,
      ate_es_var = ate_es_var,
      nde_es = nde_es,
      nde_es_var = nde_es_var,
      nie_es = nie_es,
      nie_es_var = nie_es_var
    )
    results <- dplyr::bind_rows(results, temp_results)
  }
  return(results)
}