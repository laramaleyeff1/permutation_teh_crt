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
                                   data_12_cov[,c("AGE","BMI","BL_avg_daily")],
                                   data_12_cov[,c("FEMALE", "disable", "Current_Smoke",
                                           "Alcohol_Abuse", "Drug_Abuse", "Diabetes",
                                           "CVD","Hypertension", "Pulmonary",
                                           "Dep_OR_Anx", paste0("pain",1:11),
                                           "BL_benzo_flag"
                                   )],
                                   adj = TRUE,
                                   R = R,
                                   grid.size = 51,
                                   test.stat = getSKSAdj)

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


# ====================================================
# Ad-hoc analysis to detect source of heterogeneity
# ====================================================

cov_names <- c(
  "Age, y",
  "Sex",
  "Receives disability benefits",
  "Current smoker",
  "BMI",
  "Alcohol misuse",
  "Drug misuse",
  "Diabetes",
  "Cardiovascular disorder",
  "Hypertension",
  "Chronic pulmonary disease",
  "Anxiety or depression diagnosis",
  "Back and/or neck pain diagnosis",
  "General pain diagnosis",
  "Limb/extremity pain, joint pain and/or arthritic disorders diagnosis",
  "Neuropathy diagnosis",
  "Abdominal and/or bowel pain diagnosis",
  "Musculoskeletal chest pain diagnosis",
  "Urogenital, pelvic and menstrual pain diagnosis",
  "Headache diagnosis",
  "Other painful condition diagnosis",
  "Orofacial, ear, and/or temporomandibular disorder pain diagnosis",
  "Fibromyalgia diagnosis",
  "Average morphine milligram equivalents (MME) dose per day",
  "Benzodiazepine dispensed in 6 months prior to randomization"
)


# ====================================================
# Create bar plot to visualize impact of variable exclusion
# ====================================================

data_12_cov$AGE_scale <- scale(data_12_cov$AGE)
data_12_cov$BMI_scale <- scale(data_12_cov$BMI)
data_12_cov$BL_avg_daily_scale <- scale(data_12_cov$BL_avg_daily)

# ====================================================
# Fit full model including all covariates with a random effect for clusters
# ====================================================
full_model <- lmer(
  PEGS ~ INTERVENTION * (ns(AGE_scale, df = 4) + ns(BMI_scale, df = 4) + ns(BL_avg_daily_scale, df = 4) +
                           FEMALE + disable + Current_Smoke +
                           Alcohol_Abuse + Drug_Abuse + Diabetes +
                           CVD + Hypertension + Pulmonary + Dep_OR_Anx +
                           pain1 + pain2 + pain3 + pain4 + pain5 + pain6 +
                           pain7 + pain8 + pain9 + pain10 + pain11 +
                           BL_benzo_flag) + (1 | CLUST),
  data = data_12_cov
)



# Compute full CATE estimates
data_12_cov$tau_hat <- predict(full_model, newdata = transform(data_12_cov, INTERVENTION = 1)) -
  predict(full_model, newdata = transform(data_12_cov, INTERVENTION = 0))

var_full <- var(data_12_cov$tau_hat)

covariates <- c("AGE_scale", "FEMALE", "disable", "Current_Smoke", "BMI_scale",
                "Alcohol_Abuse", "Drug_Abuse", "Diabetes", "CVD",
                "Hypertension", "Pulmonary", "Dep_OR_Anx",
                paste0("pain", 1:11), "BL_avg_daily_scale", "BL_benzo_flag")

continuous_vars <- paste0(c("AGE", "BMI", "BL_avg_daily"), "_scale")
categorical_vars <- setdiff(covariates, continuous_vars)



# Function to compute individual TE-VIM values
compute_individual_TE_VIM <- function(covariate, data) {
  continuous_terms <- sapply(
    continuous_vars,
    function(var) if (var == covariate) NULL else paste0("ns(", var, ", df = 4)")
  )
  categorical_terms <- sapply(
    categorical_vars,
    function(var) if (var == covariate) NULL else var
  )
  
  reduced_covariates <- paste(c(continuous_terms, categorical_terms)[!sapply(c(continuous_terms, categorical_terms), is.null)], collapse = " + ")
  
  reduced_formula <- as.formula(
    paste("tau_hat ~", paste(reduced_covariates, collapse = " + "), "+ (1|CLUST)")
  )                   
  
  reduced_model <- lmer(reduced_formula, data = data)
  
  # Estimate tau_s(X) using reduced model
  data$tau_s_hat <- predict(reduced_model, newdata = data)
  
  # Compute individual TE-VIM values
  te_vim_individual <- (data$tau_hat - data$tau_s_hat)^2
  
  return(te_vim_individual)
}

# Compute individual TE-VIM values for each covariate
individual_te_vim_list <- lapply(covariates, compute_individual_TE_VIM, data = data_12_cov)

# Combine into a data frame
individual_te_vim_df <- do.call(cbind, individual_te_vim_list)
colnames(individual_te_vim_df) <- cov_names

# Convert to long format for ggplot
individual_te_vim_long <- pivot_longer(as.data.frame(individual_te_vim_df), cols = everything(), names_to = "Covariate", values_to = "TE_VIM")

# Compute median TE-VIM for sorting
te_vim_median <- individual_te_vim_long %>%
  group_by(Covariate) %>%
  dplyr::summarise(Median_TE_VIM = median(TE_VIM, na.rm = TRUE)) %>%
  arrange(Median_TE_VIM)

individual_te_vim_long$Covariate <- factor(individual_te_vim_long$Covariate, levels = te_vim_median$Covariate)

# Create box plot
ggplot(individual_te_vim_long, aes(y = Covariate, x = TE_VIM)) +
  geom_boxplot(outlier.shape=NA, fill = "white", color = "black") +  # Boxplot without outliers
  # geom_jitter(width = 0.1, alpha = 0.3, color = "black") +  # Jittered points for better visualization
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  # Reference line at zero
  theme_minimal(base_size = 14) +
  coord_cartesian(xlim = c(0, 0.36)) +  # Display only x-axis range [0, 0.36] while keeping full data
  labs(title = "LOO TE-VIM Distribution Across Covariates",
       x = "LOO TE-VIM Estimate", y = "Covariate") +
  theme(panel.grid.major.y = element_blank(),  # Remove horizontal gridlines
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(face = "bold"))

