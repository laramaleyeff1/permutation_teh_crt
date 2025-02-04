library(haven)
library(dplyr)
library(lme4)
library(xtable)
library(hettx)
library(tidyr)
library(ggplot2)

# ====================================================
# Load required data and functions
# ====================================================

R = 2000  # Number of permutations

# Set working directory to source file location
source("./permutation_test_crt.R")

# Load the dataset
data <- read_sas("../ppact_public_bpi.sas7bdat")
# The dataset "ppact_public_bpi.sas7bdat" can be downloaded from:
# https://github.com/PainResearch/PPACT
# Example usage:
# data <- read_sas("ppact_public_bpi.sas7bdat")


# ====================================================
# Data Preprocessing: Filtering and Selection
# ====================================================

# Filter data for timepoint 12 and select relevant variables
data_12_cov = data %>%
  filter(TIMEPOINT == 12) %>%  # Keep only rows where TIMEPOINT is 12
  dplyr::select(PEGS, INTERVENTION, CLUST,  # Key variables
                AGE, FEMALE, disable, Current_Smoke, BMI,  # Covariates
                Alcohol_Abuse, Drug_Abuse, Diabetes, CVD,
                Hypertension, Pulmonary, Dep_OR_Anx,
                pain1:pain11,  # Pain variables
                BL_avg_daily,  # Baseline daily average
                BL_benzo_flag   # Baseline benzodiazepine use flag
  ) %>%
  na.omit()  # Remove rows with missing values

# Exclude individuals with extreme baseline average daily values (greater than 200)
data_12_cov <- data_12_cov[data_12_cov$BL_avg_daily < 200,]

# Aggregate data at the cluster level by computing the mean of all variables except cluster ID
data_12_cov_cl_level = data_12_cov %>%
  group_by(CLUST) %>%
  summarise_at(vars(-group_cols()), mean) %>%  # Compute mean for each cluster
  ungroup()  # Remove grouping


# ====================================================
# Permutation Test for Treatment Effect Heterogeneity (Ding et al.)
# ====================================================

# Set seed for reproducibility
set.seed(123)

# Conduct unadjusted cluster-level test
ding_unadj = detect_idiosyncratic(PEGS ~ INTERVENTION, data_12_cov_cl_level, B=2000)

# Output:
# Statistic           P-Value (Sweep)     P-Value (Plug-In)
# 1 0.09433962        0.9276              0.9161


# ====================================================
# Adjusted Permutation Test for Treatment Effect Heterogeneity
# ====================================================

# Construct regression formula excluding PEGS, INTERVENTION, CLUST
formula = paste0("PEGS ~ ", paste(names(data_12_cov_cl_level)[-c(1:3)], collapse = "+"))

# Fit a linear model using only the control group (INTERVENTION == 0)
model <- lm(as.formula(formula), data = data_12_cov_cl_level[data_12_cov_cl_level$INTERVENTION == 0,])

# Construct model matrix
mod_mat = model.matrix(as.formula(formula), data = data_12_cov_cl_level)

# Compute residuals by subtracting predicted values from observed PEGS values
data_12_cov_cl_level$resids = c(data_12_cov_cl_level$PEGS - mod_mat %*% coef(model))

# Conduct adjusted test
set.seed(123)
ding_adj = detect_idiosyncratic(resids ~ INTERVENTION, data_12_cov_cl_level, B=2000)

# Output:
# Statistic           P-Value (Sweep)     P-Value (Plug-In)
# 1 0.2075472         0.0376              0.0296


# ====================================================
# Permutation Tests Based on GK-S and SK-S
# ====================================================

set.seed(1234)
perm.test.adj = runPermutationTest(data_12_cov$PEGS,
                                   data_12_cov$INTERVENTION,
                                   data_12_cov$CLUST,
                                   data_12_cov[,c("AGE", "BL_avg_daily")],
                                   data_12_cov[,c("disable",
                                                  "Diabetes",
                                                  "Hypertension",
                                                  paste0("pain",c(3,5,11))
                                   )], adj=TRUE, R=R)

set.seed(5678)
perm.test.unadj = runPermutationTest(data_12_cov$PEGS,
                                     data_12_cov$INTERVENTION,
                                     data_12_cov$CLUST, R=R)

# Outputs:
perm.test.adj
# $p.value.ci
# [1] 0.03258376
# $p.value.pi
# [1] 0.02758626

perm.test.unadj
# $p.value.ci
# [1] 0.1265368
# $p.value.pi
# [1] 0.03358326


# ====================================================
# Likelihood Ratio Test Analysis
# ====================================================

# Fit null model
model_null <- lmer(PEGS ~ INTERVENTION + AGE + FEMALE + disable +
                     Current_Smoke + BMI + Alcohol_Abuse +
                     Drug_Abuse  + Diabetes +
                     CVD + Hypertension +
                     Pulmonary + Dep_OR_Anx +
                     pain1 + pain2 + pain3 + pain4 +
                     pain5 + pain6 + pain7 + pain8 +
                     pain9 + pain10 + pain11 + BL_avg_daily + BL_benzo_flag +
                     (1|CLUST), data=data_12_cov)

# Fit full interaction model
model_full_combined <- lmer(PEGS ~ INTERVENTION * (AGE + FEMALE + disable +
                                                     Current_Smoke + BMI + Alcohol_Abuse +
                                                     Drug_Abuse  + Diabetes +
                                                     CVD + Hypertension +
                                                     Pulmonary + Dep_OR_Anx +
                                                     pain1 + pain2 + pain3 + pain4 +
                                                     pain5 + pain6 + pain7 + pain8 +
                                                     pain9 + pain10 + pain11 + BL_avg_daily +
                                                     BL_benzo_flag) +
                              (1|CLUST), data=data_12_cov)

# Perform likelihood ratio test
anova(model_full_combined, model_null)


# ====================================================
# Spline-Based Likelihood Ratio Test Analysis
# ====================================================

# Fit full model with splines
model_full_spline <- lmer(
  PEGS ~ INTERVENTION * (ns(AGE, df = 4) + ns(BMI, df = 4) + ns(BL_avg_daily, df = 4) +
                           FEMALE + disable + Current_Smoke +
                           Alcohol_Abuse + Drug_Abuse + Diabetes +
                           CVD + Hypertension + Pulmonary + Dep_OR_Anx +
                           pain1 + pain2 + pain3 + pain4 + pain5 + pain6 +
                           pain7 + pain8 + pain9 + pain10 + pain11 +
                           BL_benzo_flag) + (1 | CLUST),
  data = data_12_cov
)

# Fit null model with splines
model_null_spline <- lmer(
  PEGS ~ INTERVENTION + ns(AGE, df = 4) + ns(BMI, df = 4) + ns(BL_avg_daily, df = 4) +
    FEMALE + disable + Current_Smoke +
    Alcohol_Abuse + Drug_Abuse + Diabetes +
    CVD + Hypertension + Pulmonary + Dep_OR_Anx +
    pain1 + pain2 + pain3 + pain4 + pain5 + pain6 +
    pain7 + pain8 + pain9 + pain10 + pain11 + BL_benzo_flag +
    (1 | CLUST),
  data = data_12_cov
)

# Perform likelihood ratio test
anova(model_full_spline, model_null_spline)

# Stepwise selection to identify key predictors for adjustment in permutation test
step(model_full_spline, direction = "backward", reduce.random = FALSE)



# ====================================================
# Ad-hoc analysis to detect source of heterogeneity
# ====================================================

# Scale continuous covariates
data_12_cov$AGE_scale <- scale(data_12_cov$AGE)
data_12_cov$BMI_scale <- scale(data_12_cov$BMI)
data_12_cov$BL_avg_daily_scale <- scale(data_12_cov$BL_avg_daily)

# ====================================================
# Fit full model including all covariates with a random effect for clusters
# ====================================================
model_full_spline_scaled <- lmer(
  PEGS ~ INTERVENTION * (ns(AGE_scale, df = 4) + ns(BMI_scale, df = 4) + ns(BL_avg_daily_scale, df = 4) +
                           FEMALE + disable + Current_Smoke +
                           Alcohol_Abuse + Drug_Abuse + Diabetes +
                           CVD + Hypertension + Pulmonary + Dep_OR_Anx +
                           pain1 + pain2 + pain3 + pain4 + pain5 + pain6 +
                           pain7 + pain8 + pain9 + pain10 + pain11 +
                           BL_benzo_flag) + (1 | CLUST),
  data = data_12_cov
)

# Extract fixed effects and treatment effect estimates
fixed_effects <- fixef(model_full_spline_scaled)
intervention_effects <- fixed_effects[grep("INTERVENTION", names(fixed_effects))]

# Construct model matrix for covariates
X_mat <- model.matrix(~ ns(AGE_scale, df = 4) + ns(BMI_scale, df = 4) + ns(BL_avg_daily_scale, df = 4) +
                        FEMALE + disable + Current_Smoke +
                        Alcohol_Abuse + Drug_Abuse + Diabetes +
                        CVD + Hypertension + Pulmonary + Dep_OR_Anx +
                        pain1 + pain2 + pain3 + pain4 + pain5 + pain6 +
                        pain7 + pain8 + pain9 + pain10 + pain11 + BL_benzo_flag,
                      data = data_12_cov)

# Compute the full model's estimated treatment effect
full_treatment_effect <- mean(X_mat %*% intervention_effects)

# ====================================================
# Function to construct a reduced model formula excluding a covariate
# ====================================================
build_reduced_formula <- function(excluded_var, continuous_vars, categorical_vars, model_matrix = FALSE) {
  continuous_terms <- sapply(
    continuous_vars,
    function(var) if (var == excluded_var) NULL else paste0("ns(", var, ", df = 4)")
  )
  categorical_terms <- sapply(
    categorical_vars,
    function(var) if (var == excluded_var) NULL else var
  )
  terms <- paste(c(continuous_terms, categorical_terms)[!sapply(c(continuous_terms, categorical_terms), is.null)], collapse = " + ")

  if (!model_matrix) {
    return(as.formula(paste("PEGS ~ INTERVENTION * (", terms, ") + (1 | CLUST)")))
  } else {
    return(as.formula(paste(" ~ ", terms)))
  }
}

# ====================================================
# Define covariates and separate into continuous and categorical groups
# ====================================================
covariates <- c("AGE_scale", "FEMALE", "disable", "Current_Smoke", "BMI_scale",
                "Alcohol_Abuse", "Drug_Abuse", "Diabetes", "CVD",
                "Hypertension", "Pulmonary", "Dep_OR_Anx",
                paste0("pain", 1:11), "BL_avg_daily_scale", "BL_benzo_flag")

continuous_vars <- paste0(c("AGE", "BMI", "BL_avg_daily"), "_scale")
categorical_vars <- setdiff(covariates, continuous_vars)

# Store treatment effects for reduced models
results <- data.frame(Covariate = covariates, Delta = NA)

# ====================================================
# Loop through each covariate to fit reduced models
# ====================================================
for (covariate in covariates) {
  # Build the reduced formula
  reduced_formula <- build_reduced_formula(
    excluded_var = covariate,
    continuous_vars = continuous_vars,
    categorical_vars = categorical_vars
  )

  # Build the reduced model matrix formula
  model_matrix_reduced_formula <- build_reduced_formula(
    excluded_var = covariate,
    continuous_vars = continuous_vars,
    categorical_vars = categorical_vars,
    model_matrix = TRUE
  )

  # Fit the reduced model
  reduced_model <- lmer(reduced_formula, data = data_12_cov)

  # Construct model matrix for the reduced model
  X_mat = model.matrix(model_matrix_reduced_formula, data = data_12_cov)

  # Extract fixed effects and treatment effect estimates
  fixed_effects <- fixef(reduced_model)
  intervention_effects <- fixed_effects[grep("INTERVENTION", names(fixed_effects))]

  # Store results
  results$Delta[results$Covariate == covariate] <- mean(X_mat %*% intervention_effects)
}

# ====================================================
# Compute percent difference in treatment effect estimation
# ====================================================
results$Original_Delta <- full_treatment_effect
results$Diff <- results$Original_Delta - results$Delta
results$Percent_Difference <- (results$Original_Delta - results$Delta) / results$Original_Delta * 100

# Sort results by percent difference in descending order
results_sorted <- results %>%
  arrange(desc(abs(Percent_Difference)))

# Transform data for visualization
results_sorted_long <- results_sorted %>%
  pivot_longer(cols = c(Delta), names_to = "Type", values_to = "Delta")

results_sorted_long$Covariate <- factor(results_sorted_long$Covariate,
                                        levels = covariates,
                                        labels = cov_names)

results_sorted$Covariate <- factor(results_sorted$Covariate,
                                   levels = covariates,
                                   labels = cov_names)

# ====================================================
# Create bar plot to visualize impact of variable exclusion
# ====================================================
ggplot(results_sorted_long, aes(x = reorder(Covariate, Delta), y = Delta)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_hline(yintercept = full_treatment_effect, linetype = "dashed", color = "gray", size = 0.8) +
  geom_text(
    aes(label = sprintf("%.1f%%", results_sorted$Percent_Difference[match(Covariate, results_sorted$Covariate)]),
        y = 0.02),  # Offset text above bars
    position = position_dodge(width = 0.7),
    vjust = 0.5,
    size = 3
  ) +
  labs(
    title = "Impact of Variable Exclusion on Estimated Treatment Effect",
    subtitle = sprintf("Estimated Average Treatment Effect Using All Variables: %.2f (Dashed Gray Line)", full_treatment_effect),
    y = "Estimated Leave-One-Out Average Treatment Effect",
    x = ""
  ) +
  scale_fill_manual(values = c("blue")) +
  theme_minimal() +
  coord_flip() +
  theme(legend.position = "none")


