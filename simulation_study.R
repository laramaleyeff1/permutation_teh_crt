start_time = Sys.time()

library(foreach)
library(dplyr)
library(geepack)
library(mgcv)
library(stats)
library(lme4)
library(lmerTest)
library(splines)
# ====================================================
# Simulation Parameters and Setup (Table 4 Rows 1-8)
# ====================================================

# Capture command-line arguments
params <- commandArgs(trailingOnly = TRUE)

# Define simulation parameters from command-line input
k_max <- as.numeric(params[1])  # Maximum number of clusters
n_per <- as.numeric(params[2])  # Number of individuals per cluster
R <- as.numeric(params[3])  # Number of permutations for permutation test
Delta.sd.scale <- as.numeric(params[4])  # Scaling factor for SD of treatment effect heterogeneity
index <- as.numeric(params[5])  # Index identifier for batch processing

# Default values for additional simulation parameters
Delta.mean <- 1  # Mean treatment effect
p_trt <- 0.5  # Probability of receiving treatment
icc_x <- 0.01  # Intra-cluster correlation for covariates
icc_y <- 0.1  # Intra-cluster correlation for outcome
sigma_eps_sq <- 0.1  # Variance of residual errors
mu_x <- 0  # Mean of covariates
sigma_x_sq <- 1  # Variance of covariates
sd_n_per <- 0  # Standard deviation for per-cluster sample size
gamma <- 0.0001  # Regularization parameter
grid.size <- 51  # Grid size for numerical integration

# Load necessary helper functions
source("permutation_test_crt.R")

# ====================================================
# Function: set_up_trial
# Generates a simulated cluster randomized trial dataset
# ====================================================

set_up_trial <- function(k_max, n_per, sd_n_per, p_trt, icc_y, icc_x, mu_x, sigma_x_sq, sigma_eps_sq) {

  # Generate cluster sizes with normally distributed variability
  count <- rnorm(k_max, n_per, sd_n_per)

  # Create treatment assignment for clusters
  freq_table <- data.frame(
    k = rep(1:k_max),  # Cluster IDs
    W = sample(rep(0:1, each = p_trt * k_max))  # Assign treatment/control
  )

  # Assign the generated cluster sizes to the frequency table
  freq_table$count <- count

  # Expand the dataset to include individuals within each cluster
  data <- freq_table[rep(1:nrow(freq_table), freq_table[["count"]]), -3]

  # ====================================================
  # Generate Individual-Level Covariates with Cluster Effects
  # ====================================================

  # Generate cluster-level random effects for the first covariate (X_1)
  nu_k <- rnorm(k_max, 0, sqrt(icc_x * sigma_x_sq))

  # Assign individual-level values for X_1, incorporating cluster effects
  data$X_1 <- mu_x + apply(data, 1, function(x) {
    nu_k[as.numeric(x["k"])]  # Use cluster random effect
  }) + rnorm(nrow(data), 0, sqrt((1 - icc_x) * sigma_x_sq))  # Add individual noise

  # Generate cluster-level random effects for the second covariate (X_2)
  nu_k_2 <- rnorm(k_max, 0, sqrt(icc_x * sigma_x_sq))

  # Assign individual-level values for X_2, incorporating cluster effects
  data$X_2 <- mu_x + apply(data, 1, function(x) {
    nu_k_2[as.numeric(x["k"])]  # Use second cluster random effect
  }) + rnorm(nrow(data), 0, sqrt((1 - icc_x) * sigma_x_sq))  # Add individual noise

  return(data)
}


# ====================================================
# Function: one_iter
# Simulates a single iteration of the trial under different scenarios
# ====================================================

one_iter <- function(k_max, n_per, sd_n_per, p_trt, scenario, Delta.mean,
                     Delta.sd.scale, icc_y, icc_x, mu_x, sigma_x_sq, sigma_eps_sq) {

  # Compute cluster-level variance based on intra-cluster correlation
  sigma_cl_sq <- sigma_eps_sq * (icc_y / (1 - icc_y))

  # Generate cluster-level random effects
  group_ranef <- rnorm(k_max, 0, sqrt(sigma_cl_sq))

  # Generate trial data
  data <- set_up_trial(k_max, n_per, sd_n_per, p_trt, icc_y, icc_x, mu_x, sigma_x_sq, sigma_eps_sq)

  # Rename variables for clarity
  names(data) <- c("k", "W", "X_1", "X_2")

  # Generate observed outcomes with cluster-level random effects
  data$Y_obs <- 2 * data$X_1 + 1.5 * data$X_2 +
    group_ranef[data$k] +
    rnorm(nrow(data), 0, sqrt(sigma_eps_sq))

  # ====================================================
  # Apply Different Scenarios to Modify Treatment Effect
  # ====================================================

  if (scenario == 1) {
    # Scenario 1: Null hypothesis holds (Type I error test)
    data$Y_obs <- data$Y_obs + Delta.mean * data$W
  } else if (scenario == 2) {
    # Scenario 2: Linear treatment effect modification
    Deltas = Delta.sd.scale * (0.28 * data$X_1 + 0.3 * data$X_2)
    data$Y_obs <- data$Y_obs + Delta.mean * data$W + Deltas * data$W
    Delta.mean <- mean(Delta.mean + Deltas[data$W==1])
   } else if (scenario == 3) {
    # Scenario 3: Non-linear (cosine transformation) treatment effect modification
    Deltas <- Delta.sd.scale * (0.6 * cos(6 * data$X_1 / pi))
    data$Y_obs <- data$Y_obs + Delta.mean * data$W + Deltas * data$W
    Delta.mean <- mean(Delta.mean + Deltas[data$W==1])
  } else if (scenario == 4) {
    # Scenario 4: Quadratic treatment effect modification
    Deltas <- (0.1 * data$X_1 + 2 * data$X_1^2 + data$X_2^2) * Delta.sd.scale * 0.2 / sqrt(0.1 + 2)
    data$Y_obs <- data$Y_obs + Delta.mean * data$W + Deltas * data$W
    Delta.mean <- mean(Delta.mean + Deltas[data$W == 1])
  }

  # ====================================================
  # Perform Permutation Tests
  # ====================================================

  # Unadjusted permutation test
  perm.test.unadj <- runPermutationTest(data$Y_obs, data$W, data$k, R = R, grid.size = grid.size, truth = Delta.mean)

  # Adjusted permutation test (controlling for covariates)
  perm.test.adj <- runPermutationTest(data$Y_obs, data$W, data$k,
                                      data[, c("X_1", "X_2")], X.bin = NA,
                                      adj = TRUE, R = R, truth = Delta.mean, grid.size = grid.size, test.stat=getSKSAdj)

  # ====================================================
  # Likelihood Ratio Tests for Interaction Effects
  # ====================================================

  # Standard linear models
  no_interaction_model <- lmer(Y_obs ~ W + X_1 + X_2 + (1 | k), data = data, REML = FALSE)
  interaction_model <- lmer(Y_obs ~ W * (X_1 + X_2) + (1 | k), data = data, REML = FALSE)

  # Compute p-value from likelihood ratio test
  p_value_model <- anova(interaction_model, no_interaction_model)[2, 8]

  # Spline-based interaction models
  spline_no_interaction_model <- lmer(Y_obs ~ W + ns(X_1, df = 4) + ns(X_2, df = 4) + (1 | k), data = data, REML = FALSE)
  spline_interaction_model <- lmer(Y_obs ~ W * ns(X_1, df = 4) + W * ns(X_2, df = 4) + (1 | k), data = data, REML = FALSE)

  # Compute p-value from spline-based likelihood ratio test
  p_value_spline_model <- anova(spline_interaction_model, spline_no_interaction_model)[2, 8]

  # ====================================================
  # Return Results
  # ====================================================

  returned <- data.frame(
    k_max, n_per, sd_n_per, p_trt, scenario, Delta.mean, Delta.sd.scale, gamma,
    icc_y, icc_x, mu_x, sigma_x_sq, sigma_eps_sq,
    p_value_ci = perm.test.unadj$p.value.ci,
    p_value_pi = perm.test.unadj$p.value.pi,
    p_value_truedelta = perm.test.unadj$p.value.truth,
    p_value_ci_residuals = perm.test.adj$p.value.ci,
    p_value_pi_residuals = perm.test.adj$p.value.pi,
    p_value_truedelta_residuals = perm.test.adj$p.value.truth,
    p_value_lr_test = p_value_model,
    p_value_spline_model
  )

  return(returned)
}

# ====================================================
# Function: one_run
# Runs multiple scenarios in a single execution
# ====================================================

one_run <- function() {
  ret <- foreach(scenario = 1:4, .combine = rbind) %do%
    one_iter(k_max, n_per, sd_n_per, p_trt, scenario, Delta.mean,
             Delta.sd.scale, icc_y, icc_x, mu_x, sigma_x_sq, sigma_eps_sq)
  return(ret)
}

# Run simulations
out <- foreach(iter = 1:10, .combine = rbind) %do% one_run()

# Store execution time
end_time <- Sys.time()
out$time <- difftime(end_time, start_time, units = "mins")
out$R <- R

# ====================================================
# Save Output Results
# ====================================================

setwd("results")

out.name <- paste0("power_sims_", n_per, "_", k_max, "_", R, "_",
                   Delta.sd.scale, "_", index, ".csv")

write.csv(out, out.name, row.names = FALSE)
