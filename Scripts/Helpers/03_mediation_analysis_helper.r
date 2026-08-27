##### preprocess the dataset to enable sl3 compatibility:
## drop rows with NA in X/A/M/Delta_Y; retain missing Y only when requested
## add Y_sq to the dataset
## one-hot encode X and M
preprocess_data <- function(
    data,
    y_vars,
    a_var,
    x_vars,
    m_vars,
    delta_y_var,
    flag_missingness = FALSE
) {
  if (length(flag_missingness) != 1L || is.na(flag_missingness) ||
      !is.logical(flag_missingness)) {
    stop("flag_missingness must be TRUE or FALSE.")
  }
  
  # check if all required variables are present in the data
  required_observed <- unique(c(a_var, x_vars, m_vars, delta_y_var))
  if (!flag_missingness) {
    required_observed <- unique(c(required_observed, y_vars))
  }
  keep_vars <- unique(c(required_observed, y_vars))
  missing_req <- setdiff(keep_vars, names(data))
  if (length(missing_req)) stop("Missing columns: ", paste(missing_req, collapse = ", "))
  
  # Keep binary-treatment rows with observed A/X/M/missingness indicator.
  # Missing outcomes are retained only for missingness-adjusted estimation.
  out <- data %>%
    dplyr::filter(.data[[a_var]] %in% c(0, 1)) %>%
    dplyr::select(dplyr::all_of(keep_vars))
  if (!flag_missingness) {
    out <- dplyr::filter(out, .data[[delta_y_var]] == 1)
  }
  ok <- stats::complete.cases(out[, required_observed, drop = FALSE])
  out <- out[ok, , drop = FALSE]
  if (!nrow(out)) stop("No eligible rows remain after preprocessing.")

  sanitize_sl3_col <- function(x) {
    if (inherits(x, "labelled")) x <- as.numeric(x)
    else if (is.ordered(x)) x <- factor(x)
    else if (is.factor(x)) x <- factor(x)
    else if (is.character(x)) x <- as.character(x)
    if (length(class(x)) > 1) class(x) <- class(x)[1]
    x
  }

  for (v in keep_vars) {
    out[[v]] <- sanitize_sl3_col(out[[v]])
  }

  if (!identical(sort(unique(as.numeric(out[[a_var]]))), c(0, 1))) {
    stop(a_var, " must contain both 0 and 1 after preprocessing.")
  }
  if (flag_missingness &&
      !identical(sort(unique(as.numeric(out[[delta_y_var]]))), c(0, 1))) {
    stop(delta_y_var, " must contain both 0 and 1 after preprocessing.")
  }
  if (!flag_missingness && any(out[[delta_y_var]] != 1)) {
    stop(delta_y_var, " must equal 1 in the complete-case analysis.")
  }

  for (yv in y_vars) {
    if (!is.numeric(out[[yv]])) stop(yv, " must be numeric.")
    out[[yv]] <- as.numeric(out[[yv]])
    if (any(out[[delta_y_var]] == 1 & is.na(out[[yv]]))) {
      stop(yv, " is missing where ", delta_y_var, " equals 1.")
    }
    out[[paste0(yv, "_sq")]] <- ifelse(
      out[[delta_y_var]] == 1,
      out[[yv]]^2,
      NA_real_
    )
  }

  bad_x <- x_vars[!vapply(x_vars, function(v) {
    x <- out[[v]]
    !all(is.na(x)) && length(unique(x[!is.na(x)])) > 1
  }, logical(1))]
  if (length(bad_x) > 0) {
    stop("These X variables are all-NA or constant after filtering: ", paste(bad_x, collapse = ", "))
  }

  bad_m <- m_vars[!vapply(m_vars, function(v) {
    x <- out[[v]]
    !all(is.na(x)) && length(unique(x[!is.na(x)])) > 1
  }, logical(1))]
  if (length(bad_m) > 0) {
    stop("These M variables are all-NA or constant after filtering: ", paste(bad_m, collapse = ", "))
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
run_mediation_pipeline <- function(
    df_input,
    a_var,
    x_vars,
    m_vars,
    y_vars,
    delta_y_var,
    learners,
    num_outer_folds = 2L,
    flag_missingness = FALSE
) {
  num_outer_folds <- as.integer(num_outer_folds)
  if (length(num_outer_folds) != 1L || is.na(num_outer_folds) || num_outer_folds < 2L) {
    stop("num_outer_folds must be one integer greater than or equal to 2.")
  }
  prep <- preprocess_data(
    df_input,
    a_var = a_var,
    x_vars = x_vars,
    m_vars = m_vars,
    y_vars = y_vars,
    delta_y_var = delta_y_var,
    flag_missingness = flag_missingness
  )
  df_model <- prep$data
  x_vars_num <- prep$X_vars
  m_vars_num <- prep$M_vars

  folds <- create_participant_folds(
    data = df_model,
    num_folds = num_outer_folds,
    a_var = a_var,
    delta_y_var = delta_y_var,
    flag_missingness = flag_missingness
  )

  results <- estimate_causal_estimands(
    data = df_model,
    X_vars = x_vars_num,
    M_vars = m_vars_num,
    Y_vars = y_vars,
    A_var = a_var,
    Delta_Y_var = delta_y_var,
    learners = learners,
    folds = folds,
    flag_missingness = flag_missingness
  )

  dplyr::mutate(
    results,
    treatment = a_var,
    mediator = paste(m_vars, collapse = ", "),
    outcome_observed_indicator = delta_y_var,
    n_analysis = nrow(df_model),
    n_outcome_observed = sum(df_model[[delta_y_var]] == 1),
    n_outcome_missing = sum(df_model[[delta_y_var]] == 0),
    num_outer_folds = num_outer_folds,
    missingness_adjusted = flag_missingness
  )
}


##### Create participant-level folds balanced by treatment and outcome observation.
create_participant_folds <- function(
    data,
    num_folds = 2L,
    a_var,
    delta_y_var,
    flag_missingness = FALSE,
    seed = 123L
) {
  stopifnot(a_var %in% names(data), delta_y_var %in% names(data))
  num_folds <- as.integer(num_folds)
  if (length(num_folds) != 1L || is.na(num_folds) || num_folds < 2L) {
    stop("num_folds must be one integer greater than or equal to 2.")
  }
  insufficient_support <- any(table(data[[a_var]]) < num_folds)
  if (flag_missingness) {
    insufficient_support <- insufficient_support ||
      any(table(data[[delta_y_var]]) < num_folds) ||
      any(table(data[[a_var]][data[[delta_y_var]] == 1]) < num_folds)
  }
  if (insufficient_support) {
    stop("Insufficient treatment or observed-outcome support for the requested folds.")
  }

  set.seed(seed)
  strata <- if (flag_missingness) {
    interaction(data[[a_var]], data[[delta_y_var]], drop = TRUE)
  } else {
    factor(data[[a_var]])
  }
  row_fold <- integer(nrow(data))
  for (stratum in levels(strata)) {
    indices <- which(strata == stratum)
    indices <- sample(indices, length(indices), replace = FALSE)
    row_fold[indices] <- rep(seq_len(num_folds), length.out = length(indices))
  }

  outer_folds <- lapply(seq_len(num_folds), function(fold_id) {
    validation_set <- which(row_fold == fold_id)
    training_set <- which(row_fold != fold_id)
    if (!length(training_set) || !length(validation_set)) {
      stop("An outer cross-fitting fold has an empty training or validation set.")
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


# IID variance for the participant-level HBCD sample.
calc_iid_var <- function(eif_vector) {
  ok <- is.finite(eif_vector)
  eif <- eif_vector[ok]
  if (length(eif) < 2L) return(NA_real_)
  stats::var(eif) / length(eif)
}

# define a function to estimate the causal estimands
estimate_causal_estimands <- function(data, X_vars, M_vars, Y_vars, A_var, Delta_Y_var, learners, folds, flag_missingness = FALSE) {
  if (length(flag_missingness) != 1L || is.na(flag_missingness) ||
      !is.logical(flag_missingness)) {
    stop("flag_missingness must be TRUE or FALSE.")
  }
  reg_lrnr <- learners$reg_lrnr
  cls_lrnr <- learners$cls_lrnr
  missing_lrnr <- learners$missing_lrnr
  if (is.null(missing_lrnr)) missing_lrnr <- sl3::Lrnr_glm$new()
  results <- data.frame()

  if (!flag_missingness) {
    incomplete_y <- vapply(Y_vars, function(y_var) anyNA(data[[y_var]]), logical(1))
    if (any(data[[Delta_Y_var]] != 1) || any(incomplete_y)) {
      stop(
        "When flag_missingness is FALSE, data must contain only participants ",
        "with observed outcomes. Use run_mediation_pipeline() to preprocess it."
      )
    }
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
  if (flag_missingness) {
    shared_nuisance <- c(shared_nuisance, "missing_ps_0", "missing_ps_1")
  } else {
    data$missing_ps_0 <- 1
    data$missing_ps_1 <- 1
  }
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
    if (flag_missingness) {
      assert_binary_support(train_data[[Delta_Y_var]], Delta_Y_var, fold_id)
    }

    task_ps_train <- sl3::make_sl3_Task(
      data = train_data,
      outcome = A_var,
      outcome_type = "binomial",
      covariates = X_vars
    )
    ps_fit <- cls_lrnr$train(task_ps_train)
    task_ps_valid <- sl3::make_sl3_Task(data = valid_data, covariates = X_vars)
    data$ps[valid_idx] <- ps_fit$predict(task_ps_valid)

    if (flag_missingness) {
      task_missing_train <- sl3::make_sl3_Task(
        data = train_data,
        outcome = Delta_Y_var,
        outcome_type = "binomial",
        covariates = c(X_vars, M_vars, A_var)
      )
      missing_fit <- missing_lrnr$train(task_missing_train)
      task_missing_valid_0 <- sl3::make_sl3_Task(
        data = valid_data_0,
        covariates = c(X_vars, M_vars, A_var)
      )
      task_missing_valid_1 <- sl3::make_sl3_Task(
        data = valid_data_1,
        covariates = c(X_vars, M_vars, A_var)
      )
      data$missing_ps_0[valid_idx] <- missing_fit$predict(task_missing_valid_0)
      data$missing_ps_1[valid_idx] <- missing_fit$predict(task_missing_valid_1)
    }

    task_ps_m_train <- sl3::make_sl3_Task(
      data = train_data,
      outcome = A_var,
      outcome_type = "binomial",
      covariates = c(X_vars, M_vars)
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
  bound <- 0.01
  data$ps <- pmax(pmin(data$ps, 1 - bound), bound)
  data$ps_m <- pmax(pmin(data$ps_m, 1 - bound), bound)
  if (flag_missingness) {
    data$missing_ps_0 <- pmax(pmin(data$missing_ps_0, 1 - bound), bound)
    data$missing_ps_1 <- pmax(pmin(data$missing_ps_1, 1 - bound), bound)
  }

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
        train_observed <- train_data[train_data[[Delta_Y_var]] == 1, , drop = FALSE]

        if (!nrow(train_observed)) {
          stop("Outer training fold ", fold_id, " has no observed outcomes for ", yv, ".")
        }
        assert_binary_support(train_observed[[A_var]], paste0(A_var, " among observed outcomes"), fold_id)

        train_data_0 <- train_data
        train_data_0[[A_var]] <- 0
        train_data_1 <- train_data
        train_data_1[[A_var]] <- 1
        valid_data_0 <- valid_data
        valid_data_0[[A_var]] <- 0
        valid_data_1 <- valid_data
        valid_data_1[[A_var]] <- 1

        task_outcome_train <- sl3::make_sl3_Task(
          data = train_observed,
          outcome = yv,
          outcome_type = "continuous",
          covariates = c(X_vars, M_vars, A_var)
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

        task_sequential_1_0 <- sl3::make_sl3_Task(
          data = train_data,
          outcome = paste0(yv, "_or_1"),
          outcome_type = "continuous",
          covariates = c(X_vars, A_var)
        )
        task_sequential_0_0 <- sl3::make_sl3_Task(
          data = train_data,
          outcome = paste0(yv, "_or_0"),
          outcome_type = "continuous",
          covariates = c(X_vars, A_var)
        )
        task_sequential_1_1 <- task_sequential_1_0

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
      observed_residual_0 <- ifelse(
        data[[Delta_Y_var]] == 1,
        data[[yv]] - data[[paste0(yv, "_or_0")]],
        0
      )
      observed_residual_1 <- ifelse(
        data[[Delta_Y_var]] == 1,
        data[[yv]] - data[[paste0(yv, "_or_1")]],
        0
      )
      data[[paste0(yv, "_D_0_0")]] <- (data[[A_var]] == 0) * data[[Delta_Y_var]] / (1 - data$ps) / data$missing_ps_0 *
        observed_residual_0 +
        (data[[A_var]] == 0) / (1 - data$ps) * (data[[paste0(yv, "_or_0")]] - data[[paste0(yv, "_sr_0_0")]]) +
        data[[paste0(yv, "_sr_0_0")]] - theta_0_0

      data[[paste0(yv, "_D_1_1")]] <- (data[[A_var]] == 1) * data[[Delta_Y_var]] / data$ps / data$missing_ps_1 *
        observed_residual_1 +
        (data[[A_var]] == 1) / data$ps * (data[[paste0(yv, "_or_1")]] - data[[paste0(yv, "_sr_1_1")]]) +
        data[[paste0(yv, "_sr_1_1")]] - theta_1_1

      data[[paste0(yv, "_D_1_0")]] <- (data[[A_var]] == 1) * data[[Delta_Y_var]] / data$ps / data$missing_ps_1 *
        data$dens_ratio * observed_residual_1 +
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