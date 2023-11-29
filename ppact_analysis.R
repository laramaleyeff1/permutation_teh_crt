library(haven)
library(dplyr)
library(lme4)
library(xtable)

R = 2000

# Set working directory to source file location
source("./permutation_test_crt.R")

# Download "ppact_public_bpi.sas7bdat" from https://github.com/PainResearch/PPACT
# and load here:
#     > data <- read_sas("ppact_public_bpi.sas7bdat")

data_12_cov = data %>%
  filter(TIMEPOINT == 12) %>%
  dplyr::select(PEGS,
                INTERVENTION,
                CLUST,
                AGE,
                FEMALE,
                disable,
                Current_Smoke, 
                BMI, 
                Alcohol_Abuse,
                Drug_Abuse,
                Diabetes,
                CVD, 
                Hypertension, 
                Pulmonary,
                Dep_OR_Anx,
                pain1:pain11,
                BL_avg_daily,
                BL_benzo_flag
                ) %>%
  na.omit()

X = data_12_cov[,-c(1,2,3)]
set.seed(1234)
perm.test.adj = runPermutationTest(data_12_cov$PEGS,data_12_cov$INTERVENTION, data_12_cov$CLUST, 
                   X[,c(1,5,24)],X[,-c(1,5,24)],adj=TRUE,R=R)
set.seed(5678)
perm.test.unadj = runPermutationTest(data_12_cov$PEGS,data_12_cov$INTERVENTION, data_12_cov$CLUST,
                   R=R)

perm.test.adj
# $p.value.ci
# [1] 0.00159925

# $p.value.pi
# [1] 0.0005997501

perm.test.unadj
# $p.value.ci
# [1] 0.09505252
# 
# $p.value.pi
# [1] 0.0170915


model_null <- lmer(PEGS ~ INTERVENTION + AGE + FEMALE + disable +
                     Current_Smoke + BMI + Alcohol_Abuse +
                     Drug_Abuse  + Diabetes +
                     CVD + Hypertension + 
                     Pulmonary + Dep_OR_Anx +
                     pain1 + pain2 + pain3 + pain4 +
                     pain5 + pain6 + pain7 + pain8 + 
                     pain9 + pain10 + pain11 + BL_avg_daily + BL_benzo_flag + 
                     (1|CLUST),data=data_12_cov)
model_full_combined <- lmer(PEGS ~ INTERVENTION*(AGE + FEMALE + disable +
                                         Current_Smoke + BMI + Alcohol_Abuse +
                                         Drug_Abuse  + Diabetes +
                                         CVD + Hypertension + 
                                         Pulmonary + Dep_OR_Anx +
                                         pain1 + pain2 + pain3 + pain4 +
                                         pain5 + pain6 + pain7 + pain8 + 
                                         pain9 + pain10 + pain11 + BL_avg_daily + 
                                         BL_benzo_flag ) + 
                              (1|CLUST),data=data_12_cov)

# Combined individual p-values
p_vals_combined = car::Anova(model_full_combined, type=3)[["Pr(>Chisq)"]][28:52]
# Combined test p-values
anova(model_full_combined, model_null)

# Individual model p-values
cov_list = names(data_12_cov)[-c(1,2,3)]
p_vals_ind = rep(NA, length(cov_list))
for (i in 1:length(cov_list)) {
  em = data_12_cov[[cov_list[i]]] 
  model_null = lmer(PEGS ~ INTERVENTION + em + (1|CLUST),data=data_12_cov,REML=FALSE)
  model_full = lmer(PEGS ~ INTERVENTION*em + (1|CLUST),data=data_12_cov,REML=FALSE)
  p_vals_ind[i] = anova(model_full, model_null)[["Pr(>Chisq)"]][2]
}
cov_names = c("Age, y",
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
#print(xtable(data.frame(cov_names, p_vals_ind, p_vals_combined), digits=c(0,0,2,2)),include.rownames = FALSE)
