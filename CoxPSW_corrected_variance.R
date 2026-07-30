





# Propensity score-weighted Cox regression with corrected variance estimation

library(survival)

ipw_cox_corrected_variance <- function(
    data,
    treatment,
    covariates,
    time,
    status,
    estimand = c("ATE", "ATT", "ATO"),
    ties = "breslow",
    confidence = 0.95
) {
  
  estimand <- match.arg(estimand)
  
  required_variables <- c(treatment, covariates, time, status)
  
  if (!all(required_variables %in% names(data))) {
    missing_variables <- setdiff(required_variables, names(data))
    stop(
      "The following variables are missing from data: ",
      paste(missing_variables, collapse = ", ")
    )
  }
  
  dat <- data[, required_variables, drop = FALSE]
  
  if (anyNA(dat)) {
    stop("Missing values are not allowed in the variables used in the analysis.")
  }
  
  dat$id <- seq_len(nrow(dat))
  dat$A <- dat[[treatment]]
  dat$time <- dat[[time]]
  dat$status <- dat[[status]]
  
  if (!all(dat$A %in% c(0, 1))) {
    stop("The treatment variable must be coded as 0 or 1.")
  }
  
  if (!all(dat$status %in% c(0, 1))) {
    stop("The status variable must be coded as 0 or 1.")
  }
  
  if (length(unique(dat$A)) < 2) {
    stop("Both treatment groups must be present.")
  }
  
  if (sum(dat$status == 1) == 0) {
    stop("At least one event is required.")
  }
  
  n <- nrow(dat)
  nX <- length(covariates) + 1
  
  covX0 <- dat[, covariates, drop = FALSE]
  A <- dat$A
  
  ps_model <- glm(
    A ~ .,
    family = binomial(),
    data = data.frame(A = A, covX0, check.names = FALSE)
  )
  
  psfit <- predict(ps_model, type = "response")
  
  if (any(!is.finite(psfit)) || any(psfit <= 0 | psfit >= 1)) {
    stop("The estimated propensity scores must lie strictly between 0 and 1.")
  }
  
  if (estimand == "ATE") {
    dat$wt <- dat$A / psfit +
      (1 - dat$A) / (1 - psfit)
  } else if (estimand == "ATT") {
    dat$wt <- dat$A +
      (1 - dat$A) * psfit / (1 - psfit)
  } else {
    dat$wt <- dat$A * (1 - psfit) +
      (1 - dat$A) * psfit
  }
  
  cox_model <- coxph(
    Surv(time, status) ~ A + cluster(id),
    weights = dat$wt,
    data = dat,
    ties = ties
  )
  
  logHR <- as.numeric(coef(cox_model)[1])
  
  if (!is.finite(logHR)) {
    stop("The treatment coefficient from the weighted Cox model is not finite.")
  }
  
  cox_summary <- summary(cox_model)
  coefficient_table <- cox_summary$coefficients
  
  if ("robust se" %in% colnames(coefficient_table)) {
    uncorrected_se <- as.numeric(coefficient_table[1, "robust se"])
  } else {
    uncorrected_se <- as.numeric(coefficient_table[1, "se(coef)"])
  }
  
  eventid <- which(dat$status == 1)
  covX <- as.matrix(cbind(`(Intercept)` = 1, covX0))
  
  if (estimand == "ATE") {
    dgvec <- -dat$A * (1 - psfit) / psfit +
      (1 - dat$A) * psfit / (1 - psfit)
  } else if (estimand == "ATT") {
    dgvec <- (1 - dat$A) * psfit / (1 - psfit)
  } else {
    dgvec <- (1 - 2 * dat$A) * psfit * (1 - psfit)
  }
  
  WdevR <- t(covX * dgvec)
  
  A11x <- numeric(n)
  A12x <- matrix(0, nrow = nX, ncol = n)
  
  for (x in eventid) {
    idrs <- which(dat$time >= dat$time[x])
    
    exp_term <- exp(dat$A[idrs] * logHR)
    s0x <- sum(dat$wt[idrs] * exp_term)
    s1x <- sum(dat$wt[idrs] * exp_term * dat$A[idrs])
    
    A11x[x] <- dat$wt[x] *
      (s1x / s0x - s1x^2 / s0x^2)
    
    A12a <- (dat$A[x] - s1x / s0x) * WdevR[, x]
    
    if (length(idrs) == 1) {
      A12c <- WdevR[, idrs] * exp_term * s1x
      A12d <- WdevR[, idrs] * exp_term * dat$A[idrs]
    } else {
      A12c <- as.numeric(WdevR[, idrs, drop = FALSE] %*% exp_term) * s1x
      A12d <- as.numeric(
        WdevR[, idrs, drop = FALSE] %*%
          (exp_term * dat$A[idrs])
      )
    }
    
    A12b <- dat$wt[x] *
      (A12d / s0x - A12c / s0x^2)
    
    A12x[, x] <- -(A12a - A12b)
  }
  
  A11A12 <- c(sum(A11x), rowSums(A12x))
  
  A22 <- crossprod(
    covX,
    covX * (psfit * (1 - psfit))
  )
  
  AA <- rbind(
    A11A12,
    cbind(0, A22)
  )
  
  invAA <- tryCatch(
    solve(AA),
    error = function(e) {
      stop("The sensitivity matrix is singular and cannot be inverted.")
    }
  )
  
  event_data <- dat[dat$status == 1, , drop = FALSE]
  event_times <- unique(event_data$time)
  
  risksetFull <- data.frame(
    sumAiWi = numeric(length(event_times)),
    sumWi = numeric(length(event_times)),
    sumW1 = numeric(length(event_times)),
    sumW0 = numeric(length(event_times)),
    time = event_times
  )
  
  for (i in seq_along(event_times)) {
    event_index <- which(
      dat$time == event_times[i] &
        dat$status == 1
    )
    
    risk_index <- which(dat$time >= event_times[i])
    
    risksetFull$sumAiWi[i] <- sum(
      dat$wt[event_index] * dat$A[event_index]
    )
    
    risksetFull$sumWi[i] <- sum(dat$wt[event_index])
    
    risksetFull$sumW1[i] <- sum(
      dat$wt[risk_index] * dat$A[risk_index]
    )
    
    risksetFull$sumW0[i] <- sum(
      dat$wt[risk_index] * (1 - dat$A[risk_index])
    )
  }
  
  HRest <- exp(logHR)
  gg <- numeric(n)
  
  for (i in eventid) {
    event_time_index <- which(
      risksetFull$time == dat$time[i]
    )
    
    vv <- (
      risksetFull$sumW1[event_time_index] * HRest
    ) / (
      risksetFull$sumW1[event_time_index] * HRest +
        risksetFull$sumW0[event_time_index]
    )
    
    gg[i] <- dat$wt[i] * (dat$A[i] - vv)
  }
  
  rs11 <- numeric(n)
  rs12 <- numeric(n)
  
  for (i in seq_len(n)) {
    previous_event_index <- which(
      risksetFull$time <= dat$time[i]
    )
    
    if (length(previous_event_index) == 0) {
      next
    }
    
    reduced_riskset <- risksetFull[
      previous_event_index,
      ,
      drop = FALSE
    ]
    
    denominator <- reduced_riskset$sumW1 * HRest +
      reduced_riskset$sumW0
    
    relative_risk <- exp(logHR * dat$A[i])
    
    rs11[i] <- dat$wt[i] *
      dat$A[i] *
      relative_risk *
      sum(reduced_riskset$sumWi / denominator)
    
    rs12[i] <- dat$wt[i] *
      relative_risk *
      sum(
        reduced_riskset$sumWi *
          (reduced_riskset$sumW1 * HRest) /
          denominator^2
      )
  }
  
  eta <- gg - rs11 + rs12
  
  matpi <- covX * (dat$A - psfit)
  bbmat <- cbind(eta, matpi)
  BB <- crossprod(bbmat)
  
  corrected_variance_matrix <- invAA %*% BB %*% t(invAA)
  corrected_variance <- corrected_variance_matrix[1, 1]
  
  if (!is.finite(corrected_variance) || corrected_variance < 0) {
    stop("The corrected variance is negative or non-finite.")
  }
  
  corrected_se <- sqrt(corrected_variance)
  critical_value <- qnorm(1 - (1 - confidence) / 2)
  
  corrected_log_ci <- logHR +
    c(-1, 1) * critical_value * corrected_se
  
  uncorrected_log_ci <- logHR +
    c(-1, 1) * critical_value * uncorrected_se
  
  results <- data.frame(
    estimand = estimand,
    log_HR = logHR,
    HR = exp(logHR),
    corrected_SE = corrected_se,
    uncorrected_SE = uncorrected_se,
    corrected_HR_CI_lower = exp(corrected_log_ci[1]),
    corrected_HR_CI_upper = exp(corrected_log_ci[2]),
    uncorrected_HR_CI_lower = exp(uncorrected_log_ci[1]),
    uncorrected_HR_CI_upper = exp(uncorrected_log_ci[2])
  )
  
  output <- list(
    results = results,
    propensity_scores = psfit,
    weights = dat$wt,
    propensity_score_model = ps_model,
    weighted_cox_model = cox_model,
    corrected_variance_matrix = corrected_variance_matrix,
    sensitivity_matrix = AA,
    variability_matrix = BB,
    analysis_data = dat
  )
  
  class(output) <- "ipw_cox_corrected_variance"
  return(output)
}







# Example

set.seed(9465476)

n <- 1000

T0 <- rexp(n, rate = 0.01)

L <- rbinom(n, 1, 0.5)
L2 <- rnorm(n)
L3 <- rnorm(n)

logit_p <- -2 - 0.5 * L + 0.1 * L2 + 0.1 * L3
pA <- plogis(logit_p)
A <- rbinom(n, 1, pA)

Ttime <- T0 * exp(
  -log(0.8) * A -
    log(0.4) * L -
    log(5) * L * A -
    log(0.4) * L2 -
    log(0.4) * L3
)

C <- rexp(n, rate = 0.001)

time <- pmin(Ttime, C)
status <- as.integer(Ttime <= C)

example_data <- data.frame(
  time = time,
  status = status,
  A = A,
  L = L,
  L2 = L2,
  L3 = L3
)

fit_ate <- ipw_cox_corrected_variance(
  data = example_data,
  treatment = "A",
  covariates = c("L", "L2", "L3"),
  time = "time",
  status = "status",
  estimand = "ATT"
)

fit_ate$results
summary(fit_ate$propensity_score_model)
summary(fit_ate$weighted_cox_model)











# Comparison with the ipwCoxCSV package proposed by Shu et al.

# install.packages(
#   "ipwCoxCSV",
#   repos = c(
#     "https://shu-d.r-universe.dev",
#     "https://cloud.r-project.org"
#   )
# )

library(ipwCoxCSV)

shu_ate <- ipwCoxInd(
  data = example_data,
  indA = "A",
  indX = c("L", "L2", "L3"),
  indStatus = "status",
  indTime = "time",
  ties = "breslow",
  confidence = 0.95
)

shu_ate
