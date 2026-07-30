# ============================================================
# Estimation of estimand-specific true log hazard ratios
# for two simulation scenarios
# ============================================================

library(dplyr)
library(survival)

set.seed(987612456)

# ============================================================
# General settings
# ============================================================

filepath <- "truth_output_senario1or2/"
dir.create(filepath, showWarnings = FALSE, recursive = TRUE)

ipwvec <- c("ATE", "ATT", "ATO")

# Monte Carlo population size
n_truth <- 50000000

# ============================================================
# Scenario-specific parameters
# ============================================================

scenario_setting <- data.frame(
  scenario = 1:2,
  
  # Coefficients of L2 and L3 in the treatment assignment model
  trt_L2 = c(-0.1,0.5),
  trt_L3 = c(-0.1, 0.5),
  
  # Hazard ratios associated with L2 and L3
  # These are incorporated as -log(out_L2) * L2
  # and -log(out_L3) * L3 in the event-time model
  out_L2 = c(0.4, 0.95),
  out_L3 = c(0.4, 0.95),
  
  # Rate parameter of the exponential censoring distribution
  rateC = c(
    0.001, 0.001
  )
)

write.csv(
  scenario_setting,
  file = file.path(filepath, "scenario_setting.csv"),
  row.names = FALSE
)

# ============================================================
# Function for estimating the true value
# for each scenario and estimand
# ============================================================

calc_truth_one <- function(scenario,
                           estimand,
                           scenario_setting,
                           n_truth = 50000,
                           seed = 987612456,
                           ties = "breslow",
                           use_cluster = FALSE) {
  
  set.seed(seed + scenario * 100 + match(estimand, c("ATE", "ATT", "ATO")))
  
  sc <- scenario_setting[scenario_setting$scenario == scenario, ]
  
  if (nrow(sc) != 1) {
    stop("Exactly one row corresponding to the specified scenario must exist in scenario_setting.")
  }
  
  n <- n_truth
  
  # Generate covariates and baseline event times
  T0 <- rexp(n, rate = 0.01)
  
  L  <- rbinom(n, 1, 0.5)
  L2 <- rnorm(n)
  L3 <- rnorm(n)
  
  # Generate treatment assignments
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
  C <- rexp(n, rate = sc$rateC)
  
  time <- pmin(Ttime, C)
  status <- as.numeric(Ttime <= C)
  
  dat <- data.frame(
    id = 1:n,
    time = time,
    status = status,
    A = A,
    L = L,
    L2 = L2,
    L3 = L3
  )
  
  # Estimate propensity scores
  psmd <- glm(
    A ~ L + L2 + L3,
    family = binomial(),
    data = dat
  )
  
  psfit <- predict(psmd, type = "response")
  
  dat$ps <- psfit
  
  # Calculate propensity score weights
  if (estimand == "ATE") {
    dat$wt <- dat$A / psfit + (1 - dat$A) / (1 - psfit)
  } else if (estimand == "ATT") {
    dat$wt <- dat$A + (1 - dat$A) * (psfit / (1 - psfit))
  } else if (estimand == "ATO") {
    dat$wt <- dat$A * (1 - psfit) + (1 - dat$A) * psfit
  } else {
    stop("estimand must be one of 'ATE', 'ATT', or 'ATO'")
  }
  
  # Fit the propensity score-weighted Cox model
  if (use_cluster) {
    fit <- survival::coxph(
      survival::Surv(time, status) ~ A + cluster(id),
      weights = dat$wt,
      data = dat,
      ties = ties
    )
  } else {
    fit <- survival::coxph(
      survival::Surv(time, status) ~ A,
      weights = dat$wt,
      data = dat,
      ties = ties
    )
  }
  
  fit_summary <- summary(fit)
  
  log_beta <- as.numeric(coef(fit)[1])
  HR <- exp(log_beta)
  
  coef_table <- fit_summary$coefficients
  
  se_coef <- if ("se(coef)" %in% colnames(coef_table)) {
    as.numeric(coef_table[1, "se(coef)"])
  } else {
    NA_real_
  }
  
  robust_se <- if ("robust se" %in% colnames(coef_table)) {
    as.numeric(coef_table[1, "robust se"])
  } else {
    NA_real_
  }
  
  result <- data.frame(
    scenario = scenario,
    estimand = estimand,
    n_truth = n_truth,
    
    log_beta_truth = log_beta,
    HR_truth = HR,
    
    se_coef = se_coef,
    robust_se = robust_se,
    
    censoring_rate = mean(status == 0),
    event_rate = mean(status == 1),
    prop_A1 = mean(A == 1),
    prop_A0 = mean(A == 0),
    
    n_event = sum(status == 1),
    n_censored = sum(status == 0),
    n_A1 = sum(A == 1),
    n_A0 = sum(A == 0),
    
    mean_ps = mean(psfit),
    min_ps = min(psfit),
    max_ps = max(psfit),
    mean_wt = mean(dat$wt),
    max_wt = max(dat$wt),
    
    trt_L2 = sc$trt_L2,
    trt_L3 = sc$trt_L3,
    out_L2 = sc$out_L2,
    out_L3 = sc$out_L3,
    rateC = sc$rateC,
    
    stringsAsFactors = FALSE
  )
  
  return(result)
}

# ============================================================
# Estimate the true values for all scenarios and estimands
# ============================================================

truth_list <- list()
k <- 1

for (scenario in scenario_setting$scenario) {
  
  cat("Start truth calculation: scenario =", scenario, "\n")
  
  for (estimand in ipwvec) {
    
    cat("  estimand =", estimand, "\n")
    
    truth_list[[k]] <- calc_truth_one(
      scenario = scenario,
      estimand = estimand,
      scenario_setting = scenario_setting,
      n_truth = n_truth,
      seed = 987612456,
      use_cluster = FALSE
    )
    
    k <- k + 1
  }
}

truth_df <- do.call(rbind, truth_list)

# ============================================================
# Save the results
# ============================================================

write.csv(
  truth_df,
  file = file.path(filepath, "truth_values_by_scenario.csv"),
  row.names = FALSE
)

# Save the results separately for each scenario
for (scenario in scenario_setting$scenario) {
  
  scenario_dir <- file.path(filepath, paste0("scenario", scenario))
  dir.create(scenario_dir, showWarnings = FALSE, recursive = TRUE)
  
  tmp <- truth_df[truth_df$scenario == scenario, ]
  
  write.csv(
    tmp,
    file = file.path(scenario_dir, "truth_values.csv"),
    row.names = FALSE
  )
}

print(truth_df)

cat("\nFinished.\n")
cat("Saved to: ", file.path(filepath, "truth_values_by_scenario.csv"), "\n", sep = "")