#' Calculate BC-IRR using GEE with Automatic Centering and Delta Method SE
#'
#' @param formula A formula object (e.g., Outcome ~ Phase * Session).
#' @param id Character string indicating the subject ID variable column name.
#' @param data A data frame containing the variables.
#' @param se_method Method for standard error calculation: "kc"
#'   (Kauermann & Carroll) or "md" (Mancl & DeRouen).
#'   Defaults to "kc". For proportion or percentage outcomes,
#'   specify se_method = "md".
#' @param delay Integer. The delay parameter. If NULL (default), it is automatically calculated based on the shortest treatment duration.
#' @param family The error distribution and link function. Defaults to poisson.
#' @param corstr Character string specifying the correlation structure. Defaults to "independence".
#' @param total_var Column name containing total trials or intervals.
#'   Required for proportion and percentage outcomes.
#'   Ignored when outcome_type = "count".
#' @param outcome_type Type of outcome: "count" for frequencies or rates,
#'   "percentage" for percentages (0–100), or "proportion" for
#'   proportions (0–1). Defaults to "count".
#'   Percentage and proportion outcomes require total_var.
#'
#' @return A `bcirr` object containing the GEE model, log BC-IRR table,
#'   BC-IRR table, delay used, and the processed data. The processed data
#'   may retain rows with missing outcomes that are omitted from model fitting.
#' @importFrom stats as.formula gaussian model.frame model.matrix model.response na.omit pnorm poisson terms complete.cases na.pass
#' @importFrom dplyr group_by mutate summarise pull ungroup filter arrange select n rename
#' @importFrom rlang sym !! :=
#' @importFrom msm deltamethod
#' @importFrom tibble tibble
#' @import gee
#' @export
bc_irr_gee <- function(formula, id, data,
                       se_method = c("kc", "md"),
                       delay = NULL,
                       family = poisson,
                       corstr = "independence",
                       total_var = NULL,
                       outcome_type = "count") {

  se_method <- match.arg(se_method)
  outcome_type <- match.arg(
    outcome_type,
    choices = c("count", "proportion", "percentage")
  )

  if (outcome_type == "count") {
    total_var <- NULL
  }

  if (outcome_type != "count" && is.null(total_var)) {
    stop("total_var is required for proportion or percentage outcomes.")
  }
  dat <- as.data.frame(data)

  if (!id %in% names(dat)) stop(paste("ID variable", id, "not found."))

  # 1. Parse Formula ------------------------------
  vars <- all.vars(formula)
  needed_vars <- unique(c(vars, id, total_var))

  if (!all(needed_vars %in% names(dat))) {
    stop("Some variables in the formula or the ID variable were not found in the data.")
  }
  # Keep all rows for determing treatment onset and centering

  dat <- dat[, needed_vars, drop = FALSE]

  outcome_var <- vars[1]
  predictors <- setdiff(vars[-1], total_var)
  phase_var <- NULL
  time_var <- NULL

  for (v in predictors) {
    if (v %in% names(dat)) {
      unique_vals <- unique(dat[[v]][!is.na(dat[[v]])])
      if (length(unique_vals) == 2 && all(unique_vals %in% c(0, 1))) {
        phase_var <- v
      } else if (is.numeric(dat[[v]]) && length(unique_vals) > 2) {
        time_var <- v
      }
    }
  }

  if (!is.null(total_var)) {
    if (total_var %in% c(vars, id)) {
      stop("Total must be a separate column, supplied through total_var.")
    }

    proportion <- dat[[outcome_var]]
    total <- dat[[total_var]]

    # Convert numeric text, such as "0.6" or "60", to numbers.
    if (is.character(proportion) || is.factor(proportion)) {
      text_values <- trimws(as.character(proportion))
      text_values[text_values == ""] <- NA_character_

      if (any(grepl("%", text_values), na.rm = TRUE)) {
        stop(
          "Enter numbers without a % sign: use 0.6 in Proportion mode ",
          "or 60 in Percentage mode."
        )
      }

      numeric_values <- suppressWarnings(as.numeric(text_values))

      if (any(!is.na(text_values) & is.na(numeric_values))) {
        stop("Outcome values must be numeric. Leave missing outcomes blank.")
      }

      proportion <- numeric_values
    }

    # Check the range according to the selected format.
    upper_limit <- if (outcome_type == "percentage") 100 else 1

    if (!is.numeric(proportion) ||
        any(!is.na(proportion) &
            (!is.finite(proportion) |
             proportion < 0 |
             proportion > upper_limit))) {
      stop(
        if (outcome_type == "percentage") {
          "Percentage outcome requires values from 0 to 100, such as 60."
        } else {
          "Proportion outcome requires values from 0 to 1, such as 0.6."
        }
      )
    }

    # Convert percentages to proportions.
    if (outcome_type == "percentage") {
      proportion <- proportion / 100
    }

    if (!is.numeric(total) ||
        any(!is.na(total) &
            (!is.finite(total) | total <= 0 | total != floor(total)))) {
      stop("Total must contain positive whole numbers or missing values.")
    }

    count <- proportion * total

    if (any(!is.na(count) & !is.finite(count))) {
      stop("The converted counts must be finite.")
    }

    if (any(abs(count - round(count)) > 1e-7, na.rm = TRUE)) {
      warning(
        "Some reconstructed counts were rounded to the nearest integer. ",
        "Results use approximate counts; check the outcome values and totals.",
        call. = FALSE
      )
    }

    # Convert only the function's local analysis copy.
    # Rounding here removes negligible floating-point differences.
    dat[[outcome_var]] <- round(count)
  }

  terms_obj <- terms(formula)
  term_labels <- attr(terms_obj, "term.labels")
  has_interaction <- any(grepl(":", term_labels))
  has_time_main <- !is.null(time_var) && (time_var %in% term_labels)

  # 2. Centering Logic & Data Preparation -----------------------------------
  dat_processed <- dat
  final_formula <- NULL

  calc_delay <- delay

  # Target Variable Names for Output (Updated per request)
  name_trend <- "Baseline_Trend"
  name_phase_new <- "Immediate_Treatment_Effect"
  name_slope_shift <- "Trend_Treatment_Effect"

  if (has_interaction || has_time_main) {
    if (is.null(time_var)) stop("Time variable needed but not found.")

    if (has_interaction) {
      if (is.null(delay)) {
        delay_vars <- c(id, phase_var, time_var)
        trt_lengths <- dat[
          complete.cases(dat[, delay_vars, drop = FALSE]),
          , drop = FALSE
        ] %>%
          group_by(!!sym(id)) %>%
          filter(!!sym(phase_var) == 1) %>%
          summarise(n = dplyr::n(), .groups = 'drop') %>%
          pull(n)

        if (length(trt_lengths) > 0) {
          min_trt <- min(trt_lengths)
          calc_delay <- max(0, min_trt - 1)
        } else {
          calc_delay <- 0
        }
      }
    } else {
      if (!is.null(delay)) {
        warning("A 'delay' parameter was provided, but no interaction effect was detected. Delay is ignored.")
      }
      calc_delay <- NULL
    }

    dat_processed <- dat_processed %>%
      group_by(!!sym(id)) %>%
      mutate(
        First_Trt = min(ifelse(!!sym(phase_var) == 1, !!sym(time_var), NA), na.rm = TRUE),
        !!name_trend := !!sym(time_var) - First_Trt + 1,
        raw_int = !!sym(name_trend) * !!sym(phase_var),
        !!name_slope_shift := ifelse(raw_int == 0, 0, raw_int - 1)
      ) %>%
      ungroup() %>%
      as.data.frame()
    if (any(!is.finite(dat_processed$First_Trt))) {
      stop(
        "Cannot determine treatment start for every participant. Check phase and session values."
      )
    }

    dat_processed[[name_phase_new]] <- dat_processed[[phase_var]]

    rhs <- paste(name_trend, "+", name_phase_new)
    if (has_interaction) rhs <- paste(rhs, "+", name_slope_shift)
    final_formula <- as.formula(paste(outcome_var, "~", rhs))

  } else {
    if (!is.null(delay)) {
      warning("A 'delay' parameter was provided, but no time variable or interaction effect was detected in the formula. Delay is ignored.")
    }
    calc_delay <- NULL

    if (!is.null(phase_var)) {
      dat_processed[[name_phase_new]] <- dat_processed[[phase_var]]
      final_formula <- as.formula(paste(outcome_var, "~", name_phase_new))
    } else {
      final_formula <- formula
    }
  }
  if (outcome_type != "count") {
    offset_term <- call("offset", call("log", as.name(total_var)))
    final_formula[[3]] <- call("+", final_formula[[3]], offset_term)
  }

  dat_processed <- dat_processed[order(dat_processed[[id]]), ]
  id_vec <- dat_processed[[id]]
  environment(final_formula) <- environment()

  # 3. Fit Model (Silenced) -------------------------------------------------
  invisible(utils::capture.output({
    suppressMessages({
      gee_fit <- gee::gee(final_formula, id = id_vec, data = dat_processed,
                          family = family, corstr = corstr, na.action = na.omit, silent = TRUE)
    })
  }))

  names(gee_fit$coefficients)[names(gee_fit$coefficients) == "(Intercept)"] <- "Baseline_Level"

  # 4. Apply Correction (SE) ------------------------------------------------
  if (se_method == "kc") {
    var_res <- GEE.var.kc.r(final_formula, id = id, family = family, data = dat_processed, corstr = corstr)
    se_label <- "Kauermann & Carroll (KC)"
  } else {
    var_res <- GEE.var.md.r(final_formula, id = id, family = family, data = dat_processed, corstr = corstr)
    se_label <- "Mancl & DeRouen (MD)"
  }

  cov_mat <- var_res$cov.beta
  nms <- names(gee_fit$coefficients)
  rownames(cov_mat) <- nms
  colnames(cov_mat) <- nms

  # 5. BC-IRR Calculation ---------------------------------------------------
  target_phase <- name_phase_new
  target_int <- name_slope_shift

  if (target_phase %in% nms) {
    idx_p <- which(nms == target_phase)

    if (has_interaction && target_int %in% nms) {
      idx_i <- which(nms == target_int)

      beta_p <- gee_fit$coefficients[idx_p]
      beta_i <- gee_fit$coefficients[idx_i]

      log_est <- beta_p + beta_i * calc_delay

      dm_formula_string <- paste("~ x1 + x2 *", calc_delay)
      se_log <- msm::deltamethod(
        as.formula(dm_formula_string),
        mean = c(beta_p, beta_i),
        cov = cov_mat[c(idx_p, idx_i), c(idx_p, idx_i)]
      )
    } else {
      beta_p <- gee_fit$coefficients[idx_p]
      log_est <- beta_p
      se_log <- sqrt(cov_mat[idx_p, idx_p])
      calc_delay <- NULL
    }

    est_irr <- exp(log_est)
    se_irr <- msm::deltamethod(~exp(x1), log_est, se_log^2)

    z_stat <- log_est / se_log

    p_raw <- 2 * (1 - pnorm(abs(z_stat)))
    p_val_fmt <- ifelse(p_raw < 0.001, "<0.001", as.character(round(p_raw, 3)))

    delay_label <- if (is.null(calc_delay)) "" else paste0(" (Delay=", calc_delay, ")")

    log_table <- tibble::tibble(
      Estimator = paste0("Log BC-IRR", delay_label),
      SE_Method = se_label,
      Estimate = log_est,
      Std_Error = se_log,
      `95% CI Lower` = log_est - 1.96 * se_log,
      `95% CI Upper` = log_est + 1.96 * se_log,
      Test_Statistic = z_stat,
      P_Value = p_val_fmt
    )

    irr_table <- tibble::tibble(
      Estimator = paste0("BC-IRR", delay_label),
      SE_Method = se_label,
      Estimate = est_irr,
      Std_Error = se_irr,
      `95% CI Lower` = exp(log_est - 1.96 * se_log),
      `95% CI Upper` = exp(log_est + 1.96 * se_log),
      P_Value = p_val_fmt
    )

    # Return trimmed list, assign S3 class
    res <- list(
      model = gee_fit,
      log_bc_irr = log_table,
      bc_irr = irr_table,
      delay_used = calc_delay,
      centered_data = dat_processed
    )
    class(res) <- "bcirr"
    return(res)

  } else {
    res <- list(
      model = gee_fit,
      log_bc_irr = NULL,
      bc_irr = NULL,
      delay_used = NULL,
      centered_data = dat_processed
    )
    class(res) <- "bcirr"
    return(res)
  }
}

# S3 Methods for beautiful printing ----------------------------------------

#' Print method for bcirr objects
#' @param x A bcirr object
#' @param ... Additional arguments
#' @export
print.bcirr <- function(x, ...) {
  if (is.null(x$bc_irr)) {
    cat("BC-IRR could not be calculated. Model output:\n")
    print(x$model)
  } else {
    print(x$bc_irr)
  }
  invisible(x)
}

#' Summary method for bcirr objects
#' @param object A bcirr object
#' @param ... Additional arguments
#' @export
summary.bcirr <- function(object, ...) {
  cat("\n=== BC-IRR Model Summary ===\n\n")

  cat("--- Log-Scale Estimates ---\n")
  print(object$log_bc_irr)

  cat("\n--- Exponentiated Estimates (BC-IRR) ---\n")
  print(object$bc_irr)

  if (!is.null(object$delay_used)) {
    cat("\nNote: Calculated using Delay =", object$delay_used, "\n")
  }
  invisible(object)
}
