# ============================================================
# Main simulation
# Scenario-specific treatment, outcome, and censoring parameters
# Estimand-specific true values are read from truth_values_by_scenario.csv
# Parallel processing with progress reporting
# Progress is also recorded in progress_log.txt
# ============================================================
library(dplyr)
library(survival)
library(future.apply)
library(writexl)
library(progressr)
set.seed(94652)
# =========================
# Parallel processing
# =========================
n_workers <- 12
plan(multisession, workers = n_workers)
# Progress reporting
handlers(global = TRUE)
handlers("txtprogressbar")
# =========================
# General settings
# =========================
filepath <- "senario1/"
dir.create(filepath, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(filepath, "progress_log.txt")
writeLines(
  paste0("Simulation started: ", Sys.time()),
  con = log_file
)
ipwvec <- c("ATE", "ATT", "ATO")
# Sample sizes
nvec <- c(300,1000,10000)
# Number of simulation replicates per condition
N <- 10000
# =========================
# Read estimand-specific true values
# =========================
truth_df <- read.csv(
  "truth_output_senario1or2/truth_values_by_scenario.csv",
  check.names = FALSE
)
required_truth_cols <- c("scenario", "estimand", "log_beta_truth")
if (!all(required_truth_cols %in% names(truth_df))) {
  stop(
    "The required columns are missing from truth_values_by_scenario.csv: ",
    paste(required_truth_cols, collapse = ", ")
  )
}
# =========================
# Scenario-specific parameters
# =========================


scenario_setting <- data.frame(
  scenario = 1:2,
  
  # Coefficients of L2 and L3 in the treatment assignment model
  trt_L2 = c(-0.1,0.5),
  trt_L3 = c(-0.1, 0.5),
  
  # Hazard ratios associated with L2 and L3
  # Incorporated as -log(out_L2) * L2 and -log(out_L3) * L3 in the event-time model
  out_L2 = c(0.4, 0.95),
  out_L3 = c(0.4, 0.95),
  
  # Rate parameter of the exponential censoring distribution
  rateC = c(
    0.001, 0.001
  )
)
# =========================
# Indicator of 95% confidence interval coverage
# =========================
in_ci95_hr <- function(loghr_hat, se_beta, loghr_true) {
  lower <- loghr_hat - 1.96 * se_beta
  upper <- loghr_hat + 1.96 * se_beta
  as.integer(loghr_true >= lower & loghr_true <= upper)
}
# =========================
# Extract the standard error from a Cox model summary
# =========================
get_cox_se <- function(fit2) {
  cn <- colnames(fit2$coefficients)
  
  if ("robust se" %in% cn) {
    return(as.numeric(fit2$coefficients[1, "robust se"]))
  }
  
  if ("se(coef)" %in% cn) {
    return(as.numeric(fit2$coefficients[1, "se(coef)"]))
  }
  
  return(NA_real_)
}
# =========================
# Return missing values when model estimation fails
# =========================
empty_result <- function(scenario, Estimand, n, m, loghr_true = NA_real_) {
  data.frame(
    scenario = scenario,
    estimand = Estimand,
    n = n,
    m = m,
    truth = loghr_true,
    theta = NA_real_,
    std_errors_FALSE = NA_real_,
    std_errors_TRUE = NA_real_,
    coverage_FALSE = NA_real_,
    coverage_TRUE = NA_real_,
    CI_lowhigh_FALSE = NA_real_,
    CI_lowhigh_TRUE = NA_real_
  )
}
# =========================
# Summarize simulation results
# =========================
summarize_sim_simple <- function(est_csv,
                                 cov_csv,
                                 truth,
                                 se_TRUE_csv,
                                 se_FALSE_csv,
                                 est_col = NULL,
                                 cov_col = NULL,
                                 out_xlsx = "summary.xlsx",
                                 cov2_csv = NULL,
                                 cov2_col = NULL,
                                 se_TRUE_col = NULL,
                                 se_FALSE_col = NULL) {
  
  first_numlog_col <- function(df) {
    idx <- which(vapply(df, function(x) {
      is.numeric(x) || is.logical(x)
    }, logical(1)))
    
    if (length(idx) == 0) return(NULL)
    names(df)[idx[1]]
  }
  
  est_df <- read.csv(est_csv, check.names = FALSE)
  cov_df <- read.csv(cov_csv, check.names = FALSE)
  se_TRUE_df <- read.csv(se_TRUE_csv, check.names = FALSE)
  se_FALSE_df <- read.csv(se_FALSE_csv, check.names = FALSE)
  
  if (is.null(est_col)) est_col <- first_numlog_col(est_df)
  if (is.null(cov_col)) cov_col <- first_numlog_col(cov_df)
  if (is.null(se_TRUE_col)) se_TRUE_col <- first_numlog_col(se_TRUE_df)
  if (is.null(se_FALSE_col)) se_FALSE_col <- first_numlog_col(se_FALSE_df)
  
  if (is.null(est_col)) stop("No numeric estimate column was found in: ", est_csv)
  if (is.null(cov_col)) stop("No binary coverage column was found in: ", cov_csv)
  if (is.null(se_TRUE_col)) stop("No numeric corrected-SE column was found in: ", se_TRUE_csv)
  if (is.null(se_FALSE_col)) stop("No numeric uncorrected-SE column was found in: ", se_FALSE_csv)
  
  est <- est_df[[est_col]]
  cov1 <- cov_df[[cov_col]]
  se_TRUE <- se_TRUE_df[[se_TRUE_col]]
  se_FALSE <- se_FALSE_df[[se_FALSE_col]]
  
  est <- est[is.finite(est)]
  cov1 <- cov1[is.finite(as.numeric(cov1))]
  se_TRUE <- se_TRUE[is.finite(se_TRUE)]
  se_FALSE <- se_FALSE[is.finite(se_FALSE)]
  
  cov1_rate <- mean(as.integer(cov1) == 1, na.rm = TRUE)
  
  cov2_rate <- NA_real_
  R_cov2 <- NA_integer_
  cov2_file_used <- NA_character_
  cov2_col_used <- NA_character_
  
  if (!is.null(cov2_csv)) {
    cov2_df <- read.csv(cov2_csv, check.names = FALSE)
    
    if (is.null(cov2_col)) cov2_col <- first_numlog_col(cov2_df)
    if (is.null(cov2_col)) {
      stop("No binary coverage column for the alternative method was found in: ", cov2_csv)
    }
    
    cov2 <- cov2_df[[cov2_col]]
    cov2 <- cov2[is.finite(as.numeric(cov2))]
    
    cov2_rate <- mean(as.integer(cov2) == 1, na.rm = TRUE)
    R_cov2 <- length(cov2)
    cov2_file_used <- cov2_csv
    cov2_col_used <- cov2_col
  }
  
  mean_est <- mean(est, na.rm = TRUE)
  bias <- mean_est - truth
  
  empirical_se <- if (length(est) > 1) {
    sd(est, na.rm = TRUE)
  } else {
    NA_real_
  }
  
  mean_se_TRUE <- mean(se_TRUE, na.rm = TRUE)
  mean_se_FALSE <- mean(se_FALSE, na.rm = TRUE)
  
  res <- data.frame(
    truth = truth,
    R_est = length(est),
    R_cov1 = length(cov1),
    R_cov2 = R_cov2,
    mean_estimate = mean_est,
    Bias = bias,
    Empirical_SE = empirical_se,
    mean_SE_TRUE = mean_se_TRUE,
    mean_SE_FALSE = mean_se_FALSE,
    Coverage_95 = cov1_rate,
    Coverage_FALSE = cov2_rate,
    estimate_file = est_csv,
    cover1_file = cov_csv,
    cover2_file = cov2_file_used,
    se_TRUE_file = se_TRUE_csv,
    se_FALSE_file = se_FALSE_csv,
    estimate_col = est_col,
    cover1_col = cov_col,
    cover2_col = cov2_col_used,
    se_TRUE_col = se_TRUE_col,
    se_FALSE_col = se_FALSE_col,
    stringsAsFactors = FALSE
  )
  
  writexl::write_xlsx(res, out_xlsx)
  return(res)
}
# =========================
# Conduct one simulation replicate
# =========================
one_simulation <- function(m,
                           scenario,
                           ipw,
                           x2,
                           nvec,
                           ipwvec,
                           scenario_setting,
                           truth_df) {
  
  sc <- scenario_setting[scenario_setting$scenario == scenario, ]
  
  if (nrow(sc) != 1) {
    stop("The specified scenario must correspond to exactly one row in scenario_setting.")
  }
  
  n <- nvec[x2]
  
  T0 <- rexp(n, rate = 0.01)
  
  # Generate covariates
  L  <- rbinom(n, 1, 0.5)
  L2 <- rnorm(n)
  L3 <- rnorm(n)
  
  # Generate treatment assignment
  logit_p <- -2 - 0.5 * L + sc$trt_L2 * L2 + sc$trt_L3 * L3
  
  pA <- 1 / (1 + exp(-logit_p))
  A <- rbinom(n, 1, pA)
  
  # Generate event times
  Ttime <- T0 * exp(
    -log(0.8) * A -
      log(0.4) * L -
      log(5) * L * A -
      log(sc$out_L2) * L2 -
      log(sc$out_L3) * L3
  )
  
  # Generate censoring times and observed outcomes
  rateC <- sc$rateC
  C <- rexp(n, rate = rateC)
  
  time <- pmin(Ttime, C)
  status <- as.numeric(Ttime <= C)
  
  da <- data.frame(
    id = 1:n,
    time = time,
    status = status,
    A = A,
    L = L,
    L2 = L2,
    L3 = L3
  )
  
  data <- da
  indA <- "A"
  indX <- c("L", "L2", "L3")
  indStatus <- "status"
  indTime <- "time"
  ties <- "breslow"
  confidence <- 0.95
  
  Estimand <- ipwvec[ipw]
  
  # =========================
  # Retrieve the estimand-specific true value
  # =========================
  
  truth_now <- truth_df[
    truth_df$scenario == scenario &
      truth_df$estimand == Estimand,
  ]
  
  if (nrow(truth_now) != 1) {
    stop(
      "The true value could not be uniquely identified in truth_df: scenario = ",
      scenario,
      ", estimand = ",
      Estimand
    )
  }
  
  loghr_true <- truth_now$log_beta_truth
  
  # =========================
  # Fit the propensity score-weighted Cox model and estimate its variance
  # =========================
  
  dat <- data
  n <- nrow(dat)
  dat$id <- 1:n
  dat$A <- dat[, indA]
  dat$time <- dat[, indTime]
  dat$status <- dat[, indStatus]
  
  # Check whether the generated data are sufficient for model estimation
  if (length(unique(dat$A)) < 2) {
    return(empty_result(scenario, Estimand, n, m, loghr_true))
  }
  
  if (sum(dat$status == 1) < 5) {
    return(empty_result(scenario, Estimand, n, m, loghr_true))
  }
  
  if (sum(dat$status == 1 & dat$A == 1) < 2 ||
      sum(dat$status == 1 & dat$A == 0) < 2) {
    return(empty_result(scenario, Estimand, n, m, loghr_true))
  }
  
  nX <- length(indX) + 1
  covX0 <- dat[, indX]
  A <- dat$A
  
  psmd <- glm(
    A ~ .,
    family = "binomial",
    data = as.data.frame(cbind(A, covX0))
  )
  
  psfit <- predict(psmd, type = "response")
  
  # Calculate propensity score weights
  if (Estimand == "ATE") {
    dat$wt <- dat$A / psfit + (1 - dat$A) / (1 - psfit)
  } else if (Estimand == "ATT") {
    dat$wt <- dat$A + (1 - dat$A) * (psfit / (1 - psfit))
  } else if (Estimand == "ATO") {
    dat$wt <- dat$A * (1 - psfit) + (1 - dat$A) * psfit
  } else {
    stop("Estimand must be one of 'ATE', 'ATT', or 'ATO'")
  }
  
  fit <- tryCatch(
    survival::coxph(
      survival::Surv(time, status) ~ A + cluster(id),
      weights = dat$wt,
      data = dat,
      ties = ties
    ),
    error = function(e) {
      return(NULL)
    }
  )
  
  if (is.null(fit)) {
    return(empty_result(scenario, Estimand, n, m, loghr_true))
  }
  
  logHR <- fit$coefficients
  
  if (!is.finite(logHR)) {
    return(empty_result(scenario, Estimand, n, m, loghr_true))
  }
  
  eventid <- which(dat$status == 1)
  
  covX <- as.matrix(cbind(rep(1, n), covX0))
  
  if (Estimand == "ATE") {
    dgvec <- -dat$A * (1 - psfit) / psfit +
      (1 - dat$A) * psfit / (1 - psfit)
  } else if (Estimand == "ATT") {
    dgvec <- (1 - dat$A) * psfit / (1 - psfit)
  } else if (Estimand == "ATO") {
    dgvec <- (1 - 2 * dat$A) * psfit * (1 - psfit)
  } else {
    stop("Estimand must be one of 'ATE', 'ATT', or 'ATO'")
  }
  
  WdevR <- t(diag(dgvec) %*% covX)
  
  A11x <- rep(0, n)
  A12x <- matrix(0, nX, n)
  
  for (x in eventid) {
    idrs <- which(dat$time >= dat$time[x])
    
    s0x <- sum(dat$wt[idrs] * exp(dat$A[idrs] * logHR))
    s1x <- sum(dat$wt[idrs] * exp(dat$A[idrs] * logHR) * dat$A[idrs])
    
    A11x[x] <- dat$wt[x] * (s1x / s0x - s1x^2 / (s0x^2))
    
    A12a <- (dat$A[x] - s1x / s0x) * WdevR[, x]
    
    if (length(idrs) == 1) {
      A12c <- (WdevR[, idrs] * exp(dat$A[idrs] * logHR)) * s1x
      A12d <- WdevR[, idrs] * (exp(dat$A[idrs] * logHR) * dat$A[idrs])
    } else {
      A12c <- (WdevR[, idrs] %*% exp(dat$A[idrs] * logHR)) * s1x
      A12d <- WdevR[, idrs] %*% (exp(dat$A[idrs] * logHR) * dat$A[idrs])
    }
    
    A12b <- dat$wt[x] * (A12d / s0x - A12c / (s0x^2))
    A12x[, x] <- -(A12a - A12b)
  }
  
  A11A12 <- c(sum(A11x), apply(A12x, 1, sum))
  
  sumsquare <- function(u) {
    u %*% t(u)
  }
  
  A22mat <- apply(covX, 1, sumsquare) %*% (psfit * (1 - psfit))
  A22 <- matrix(apply(A22mat, 1, sum), nX, nX)
  
  AA <- as.matrix(rbind(A11A12, cbind(0, A22)))
  
  invAA <- tryCatch(
    solve(AA),
    error = function(e) {
      return(NULL)
    }
  )
  
  if (is.null(invAA)) {
    return(empty_result(scenario, Estimand, n, m, loghr_true))
  }
  
  # =========================
  # Estimate the empirical component of the sandwich variance
  # =========================
  
  eventPa <- subset(dat, dat$status == 1)
  eventimes <- unique(eventPa$time)
  
  RScol5 <- eventimes
  RScol3 <- rep(0, length(eventimes))
  RScol4 <- RScol3
  RScol1 <- RScol3
  RScol2 <- RScol3
  
  for (ii in seq_along(eventimes)) {
    idc <- which(dat$time == eventimes[ii] & dat$status == 1)
    RD1 <- dat[idc, ]
    
    RScol1[ii] <- sum(RD1$wt * RD1$A)
    RScol2[ii] <- sum(RD1$wt)
    
    RSind <- which(dat$time >= eventimes[ii])
    RS1 <- dat[RSind, ]
    
    RScol3[ii] <- sum(RS1$wt * RS1$A)
    RScol4[ii] <- sum(RS1$wt * (1 - RS1$A))
  }
  
  risksetFull <- data.frame(
    sumAiWi = RScol1,
    sumWi = RScol2,
    sumW1 = RScol3,
    sumW0 = RScol4,
    time = RScol5
  )
  
  indevent <- which(dat$status == 1)
  gg <- rep(0, nrow(dat))
  
  HRest <- exp(logHR)
  
  for (ii in indevent) {
    ind1 <- which(risksetFull$time == dat$time[ii])
    
    vv <- (risksetFull$sumW1 * HRest) /
      (risksetFull$sumW1 * HRest + risksetFull$sumW0)
    
    gg[ii] <- dat$wt[ii] * (dat$A[ii] - vv[ind1])
  }
  
  rs11 <- rep(0, nrow(dat))
  rs12 <- rs11
  
  for (ii in 1:nrow(dat)) {
    ind2 <- which(risksetFull$time <= dat$time[ii])
    rsred <- risksetFull[ind2, ]
    
    rs11[ii] <- dat$wt[ii] *
      dat$A[ii] *
      exp(log(HRest) * dat$A[ii]) *
      sum(rsred$sumWi / (rsred$sumW1 * HRest + rsred$sumW0))
    
    rs12[ii] <- dat$wt[ii] *
      exp(log(HRest) * dat$A[ii]) *
      sum(
        rsred$sumWi *
          (rsred$sumW1 * HRest) /
          ((rsred$sumW1 * HRest + rsred$sumW0)^2)
      )
  }
  
  eta <- gg - rs11 + rs12
  
  covX <- as.matrix(cbind(rep(1, n), covX0))
  
  matpi <- diag(dat$A - psfit) %*% covX
  
  bbmat <- cbind(eta, matpi)
  
  oot <- apply(bbmat, 1, sumsquare)
  BB <- matrix(apply(oot, 1, sum), nX + 1, nX + 1)
  
  propVar <- invAA %*% BB %*% t(invAA)
  proposeStdErr <- (diag(propVar)^0.5)[1]
  
  if (!is.finite(proposeStdErr)) {
    return(empty_result(scenario, Estimand, n, m, loghr_true))
  }
  
  lowProp <- logHR - qnorm(1 - (1 - confidence) / 2) * proposeStdErr
  upProp <- logHR + qnorm(1 - (1 - confidence) / 2) * proposeStdErr
  
  est <- c(logHR)
  se <- c(proposeStdErr)
  hrest <- exp(est)
  low <- exp(c(lowProp))
  up <- exp(c(upProp))
  
  output <- cbind(est, se, hrest, low, up)
  
  colnames(output) <- c(
    "log HR Estimate",
    "Standard Error",
    "HR Estimate",
    paste("HR ", confidence * 100, "% CI", "-low", sep = ""),
    paste("HR ", confidence * 100, "% CI", "-up", sep = "")
  )
  
  rownames(output) <- c("conventional weights")
  
  fit2 <- summary(fit)
  robust_se <- get_cox_se(fit2)
  
  theta_m <- output[1]
  std_errors_FALSE_m <- robust_se
  std_errors_TRUE_m <- output[2]
  
  coverage_FALSE_m <- in_ci95_hr(
    loghr_hat = output[1],
    se_beta = robust_se,
    loghr_true = loghr_true
  )
  
  coverage_TRUE_m <- in_ci95_hr(
    loghr_hat = output[1],
    se_beta = output[2],
    loghr_true = loghr_true
  )
  
  CI_lowhigh_FALSE_m <- 1.96 * 2 * robust_se
  CI_lowhigh_TRUE_m <- 1.96 * 2 * output[2]
  
  data.frame(
    scenario = scenario,
    estimand = Estimand,
    n = n,
    m = m,
    truth = loghr_true,
    theta = theta_m,
    std_errors_FALSE = std_errors_FALSE_m,
    std_errors_TRUE = std_errors_TRUE_m,
    coverage_FALSE = coverage_FALSE_m,
    coverage_TRUE = coverage_TRUE_m,
    CI_lowhigh_FALSE = CI_lowhigh_FALSE_m,
    CI_lowhigh_TRUE = CI_lowhigh_TRUE_m
  )
}
# =========================
# Run simulations across scenarios, estimands, and sample sizes
# =========================
all_summary_list <- list()
summary_index <- 1
for (scenario in scenario_setting$scenario) {
  
  scenario_dir <- file.path(filepath, paste0("scenario", scenario))
  dir.create(scenario_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (ipww in seq_along(ipwvec)) {
    
    ipw <- ipww
    
    for (x2 in seq_along(nvec)) {
      
      msg_start <- paste0(
        Sys.time(),
        " | START | scenario=", scenario,
        " | estimand=", ipwvec[ipw],
        " | n=", nvec[x2],
        " | N=", N,
        " | workers=", n_workers
      )
      
      cat("\n", msg_start, "\n", sep = "")
      write(msg_start, file = log_file, append = TRUE)
      
      with_progress({
        
        p <- progressor(steps = N)
        
        result_list <- future_lapply(
          1:N,
          function(m) {
            
            res <- one_simulation(
              m = m,
              scenario = scenario,
              ipw = ipw,
              x2 = x2,
              nvec = nvec,
              ipwvec = ipwvec,
              scenario_setting = scenario_setting,
              truth_df = truth_df
            )
            
            p(sprintf(
              "scenario %s, %s, n=%s, m=%s/%s",
              scenario, ipwvec[ipw], nvec[x2], m, N
            ))
            
            res
          },
          future.seed = TRUE
        )
      })
      
      result <- do.call(rbind, result_list)
      result <- result[order(result$m), ]
      
      prefix <- file.path(
        scenario_dir,
        paste0(ipwvec[ipw], nvec[x2], "_")
      )
      
      # =========================
      # Save replicate-level results
      # =========================
      
      theta_file <- paste0(prefix, "theta.csv")
      se_false_file <- paste0(prefix, "std_errors_FALSE.csv")
      se_true_file <- paste0(prefix, "std_errors_TRUE.csv")
      cov_false_file <- paste0(prefix, "coverage_FALSE.csv")
      cov_true_file <- paste0(prefix, "coverage_TRUE.csv")
      ci_false_file <- paste0(prefix, "CI_lowhigh_FALSE.csv")
      ci_true_file <- paste0(prefix, "CI_lowhigh_TRUE.csv")
      all_file <- paste0(prefix, "all_results.csv")
      summary_file <- paste0(prefix, "summary.xlsx")
      
      write.csv(
        data.frame(theta = result$theta),
        file = theta_file,
        row.names = FALSE
      )
      
      write.csv(
        data.frame(std_errors_FALSE = result$std_errors_FALSE),
        file = se_false_file,
        row.names = FALSE
      )
      
      write.csv(
        data.frame(std_errors_TRUE = result$std_errors_TRUE),
        file = se_true_file,
        row.names = FALSE
      )
      
      write.csv(
        data.frame(coverage_FALSE = result$coverage_FALSE),
        file = cov_false_file,
        row.names = FALSE
      )
      
      write.csv(
        data.frame(coverage_TRUE = result$coverage_TRUE),
        file = cov_true_file,
        row.names = FALSE
      )
      
      write.csv(
        data.frame(CI_lowhigh_FALSE = result$CI_lowhigh_FALSE),
        file = ci_false_file,
        row.names = FALSE
      )
      
      write.csv(
        data.frame(CI_lowhigh_TRUE = result$CI_lowhigh_TRUE),
        file = ci_true_file,
        row.names = FALSE
      )
      
      write.csv(
        result,
        file = all_file,
        row.names = FALSE
      )
      
      # =========================
      # Summarize the simulation results
      # =========================
      
      truth_row <- truth_df[
        truth_df$scenario == scenario &
          truth_df$estimand == ipwvec[ipw],
      ]
      
      if (nrow(truth_row) != 1) {
        stop(
          "The true value for the summary could not be uniquely identified: scenario = ",
          scenario,
          ", estimand = ",
          ipwvec[ipw]
        )
      }
      
      truth_now <- truth_row$log_beta_truth
      
      summary_res <- summarize_sim_simple(
        est_csv = theta_file,
        cov_csv = cov_true_file,
        cov2_csv = cov_false_file,
        se_TRUE_csv = se_true_file,
        se_FALSE_csv = se_false_file,
        truth = truth_now,
        out_xlsx = summary_file
      )
      
      summary_res$scenario <- scenario
      summary_res$estimand <- ipwvec[ipw]
      summary_res$n <- nvec[x2]
      summary_res$truth_from_csv <- truth_now
      
      all_summary_list[[summary_index]] <- summary_res
      summary_index <- summary_index + 1
      
      msg_done <- paste0(
        Sys.time(),
        " | DONE | scenario=", scenario,
        " | estimand=", ipwvec[ipw],
        " | n=", nvec[x2],
        " | truth=", truth_now
      )
      
      cat(msg_done, "\n")
      write(msg_done, file = log_file, append = TRUE)
    }
  }
}
# =========================
# Combine results across all simulation conditions
# =========================
all_summary <- do.call(rbind, all_summary_list)
write.csv(
  all_summary,
  file = file.path(filepath, "all_summary.csv"),
  row.names = FALSE
)
writexl::write_xlsx(
  all_summary,
  path = file.path(filepath, "all_summary.xlsx")
)
write(
  paste0("Simulation finished: ", Sys.time()),
  file = log_file,
  append = TRUE
)
plan(sequential)
