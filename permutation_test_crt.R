library(mgcv)
library(geepack)
library(splines)
library(dplyr)

######################################################################################################
#           Set of functions to run permutation test described in Maleyeff et al. (2024).            #
#     Author: Lara Maleyeff, with credit to Ding et al. (2016) for the first two helper functions    #
#                                   Contact: laramaleyeff@gmail.com                                  #
#                                          Date: November 2023                                       # 
######################################################################################################

#
# function getTeVec
#
# Author:   Peng Ding
# Date:     2016
#
# Function:             getTeVec
# Description:          Creates a vector that spans the confidence interval for
#                       the average treatment effect (te.hat), used in the CI method
# Parameters:           te.hat    Estimated average treatment effect, hat(Delta)
#                                 in paper
#                       te.se     Estimated standard error of te.hat
#                       gamma     Confidence interval (CI) precision, defaults 
#                                 to 0.0001
#                       grid.size numeric value indicating grid size for CI 
#                                 method, defaults to 151
# Returns:              A vector that spans the confidence interval for te.hat
# Examples of usage:
#           > getTeVec(te.hat = 1, te.se = 1, grid.size = 2)
#             [1] -2.890592  4.890592
#
getTeVec <- function(te.hat, te.se, gamma = 0.0001, grid.size = 151) {
  grid.gamma = 100*gamma
  te.MOE <- qnorm(1 - gamma/2)*te.se
  te.vec <- te.hat + (te.MOE/qnorm(1-grid.gamma)) * qnorm( seq( grid.gamma, 1-grid.gamma, length.out=grid.size ) )
  return(te.vec)
}


# function getSKS
#
# Author:   Peng Ding
# Date:     2016
#
# Function:             getSKS
# Description:          Computed the shifted Kolmogorov-Smirnoff test statistic
#                       for given outcomes (Y) and treatment indicators (W)
# Parameters:           Y         vector of continuous outcome values
#                       W         vector of treatment indicators
# Returns:              Numeric value of shifted K-S test statistic
# Examples of usage:
#           > data = data.frame(Y=...,W=...)
#           > getSKS(data$Y, data$W)
#             [1] 0.09
#
getSKS <- function(Y, W, k, X.cont, X.bin) {
  Y1 = Y[W==1]
  Y0 = Y[W==0]
  
  Y1.star   = Y1 - mean(Y1)
  Y0.star   = Y0 - mean(Y0)
  
  unique.points = c(Y1.star, Y0.star)
  
  Fn1 = ecdf(Y1.star)
  Fn0 = ecdf(Y0.star)
  
  difference = Fn1(unique.points) - Fn0(unique.points)
  
  return(max(abs(difference)))
}


getSKSAdj <- function(Y, W, k, X.cont, X.bin) {
  Y1 = Y[W==1]
  Y0 = Y[W==0]
  
  # safe_gamm <- function(...) {
  #   tryCatch({
  #     mgcv::gamm(...)
  #   }, error = function(e) {
  #     print("error")
  #     warning("GAMM failed to converge:", conditionMessage(e))
  #     return(NA)
  #   })
  # }
  
  formula_terms <- c()  # Always include W
  
  # Add smooth terms if X.cont is not NULL
  if (!is.null(X.cont)) {
    smooth_terms <- paste0("s(", names(X.cont), ")")
    formula_terms <- c(formula_terms, smooth_terms)
  }
  
  # Add binary terms if X.bin is not NULL
  if (!is.null(X.bin)) {
    binary_terms <- names(X.bin)
    formula_terms <- c(formula_terms, binary_terms)
  }
  
  # Combine into formula string
  formula_str <- paste("Y ~", paste(formula_terms, collapse = " + "))
  
  # model.gamm = mgcv::gam(as.formula(formula_str), list(k=~1),
  #                         data=cbind(X.cont, X.bin, data.frame(W, k, Y)))
  
  model.gamm = mgcv::gam(as.formula(formula_str),
                         data=cbind(X.cont, X.bin, data.frame(k, Y))[W==0,])

  Y.pred = predict.gam(model.gamm, cbind(X.cont,X.bin))
  Y1.resids = Y1 - Y.pred[W==1]
  Y0.resids = Y0 - Y.pred[W==0]
  Y1.star   = Y1.resids - mean(Y1.resids)
  Y0.star   = Y0.resids - mean(Y0.resids)
  
  unique.points = c(Y1.star, Y0.star)
  
  Fn1 = ecdf(Y1.star)
  Fn0 = ecdf(Y0.star)
  
  difference = Fn1(unique.points) - Fn0(unique.points)
  
  return(max(abs(difference)))
}

#
# function generatePermutations
#
# Author:   Lara Maleyeff (laramaleyeff@gmail.com)
# Date:     November 2023
#
# Function:             generatePermutations
# Description:          Generates permutation distribution after pre-processing
#                       in runPermutationTest(...)
# Parameters:           Y         vector of continuous outcome values
#                       W         vector of treatment indicators
#                       k         vector with numeric cluster membership
#                       te.vec    vector that spans the confidence interval for
#                                 the average treatment effect
#                       R         number permutation samples, defaults to 2000
#                       test.stat test statistic function, defaults to getSKS
# Returns:              A list with the p-values for each
#                       entry of te.vec, p.value.ci.all
#
generatePermutations <- function(Y, W, k, X.cont, X.bin, te.vec, R = 2000, test.stat = getSKS) {
  Y0.mat = sapply(te.vec, function(te) Y-W*te)
  Y1.mat = sapply(te.vec, function(te) Y+(1-W)*te)
  
  data.byk = data.frame(W=W,k=k) %>%
    group_by(k) %>%
    dplyr::summarize(W = mean(W), n=n())
  
  n.te.vec <- ncol(Y0.mat)
  
  teststat.mat <- matrix(NA, nrow = R, ncol = n.te.vec)
  teststat.obs <- test.stat(Y, W, k, X.cont, X.bin)
  
  for(r in 1:R){
    W.perm = rep(sample(data.byk$W), data.byk$n)
    ci.out <- sapply(1:n.te.vec, 
                     function(i){ 
                       Y.perm = W.perm*Y1.mat[,i] + (1-W.perm)*Y0.mat[,i]
                       teststat.star <- test.stat(Y.perm, W.perm, k, X.cont, X.bin)
                       teststat.star
                     })  
    
    teststat.mat[r,] <- as.numeric(ci.out)
  }
  
  p.value.ci.all <- sapply(1:n.te.vec, function(i){
    (1+sum(teststat.mat[,i] >= teststat.obs))/(1+R)
  })
  
  return(list(p.value.ci.all=p.value.ci.all))
  
}

#
# function runPermutationTest
#
# Author:   Lara Maleyeff (laramaleyeff@gmail.com)
# Date:     November 2023
#
# Function:             runPermutationTest
# Description:          Executes the permutation test described in Maleyeff et.
#                       al. (2024) to detect effect measure modification in 
#                       cluster randomized trials 
# Parameters:           Y         vector of continuous outcome values
#                       W         vector of treatment indicators
#                       k         vector with numeric cluster membership
#                       X.cont    named matrix of baseline continuous covariates, if 
#                                 unadjusted can set X=NA (default)
#                       X.bin     named matrix of baseline binary covariates, if 
#                                 unadjusted can set X=NA (default)
#                       adj       boolean indicating if test is GAMM-adjusted,
#                                 defaults to FALSE
#                       R         number permutation samples, defaults to 2000
#                       gamma     Confidence interval (CI) precision, defaults 
#                                 to 0.0001
#                       grid.size numeric value indicating grid size for CI 
#                                 method, defaults to 151
#                       test.stat test statistic function, defaults to getSKS
# Returns:              A list with the p-values for the plug in method (p.value.pi)  
#                       and confidence interval method (p.value.ci)
# Examples of usage:
#           > data = data.frame(Y=...,W=...,k=...)
#           > runPermutationTest(data$Y, data$W, data$k)
#             $p.value.ci
#             [1] 0.8182818
# 
#             $p.value.pi
#             [1] 0.4546455
#
runPermutationTest <- function(Y, W, k, X.cont = NA, X.bin = NA, adj = FALSE,
                               R = 2000, gamma = 0.0001, grid.size = 151, 
                               test.stat = getSKS, truth = NA) {
  
  if (adj) {
    formula_terms <- c("W")  # Always include W
    
    # Add smooth terms if X.cont is not NULL
    if (!is.null(X.cont)) {
      smooth_terms <- paste0("bs(", names(X.cont), ")")
      formula_terms <- c(formula_terms, smooth_terms)
    }
    
    # Add binary terms if X.bin is not NULL
    if (!is.null(X.bin)) {
      binary_terms <- names(X.bin)
      formula_terms <- c(formula_terms, binary_terms)
    }
    
    # Combine into formula string
    formula_str <- paste("Y ~", paste(formula_terms, collapse = " + "))
    ate.model <- geese(as.formula(formula_str), id=k,data=cbind(X.cont,X.bin,data.frame(k,Y)))
    
    te.vec <- getTeVec(ate.model[["beta"]][["W"]],
                       as.numeric(sqrt(ate.model$vbeta[2,2])),
                       gamma, grid.size)
    
    permutation.dist <- generatePermutations(Y, W, k, X.cont, X.bin, te.vec, R, test.stat)
  } else {
    ate.model <- geese(Y ~ W, id=k)
    te.vec <- getTeVec(ate.model[["beta"]][["W"]],
                         as.numeric(sqrt(ate.model$vbeta[2,2])),
                         gamma, grid.size)
    if (!is.na(truth)) {
      te.vec = c(truth,te.vec)
    }
    permutation.dist <- generatePermutations(Y, W, k, X.cont, X.bin, te.vec, R, test.stat)
  }
 
  if (!is.na(truth)) {
    p.value.ci.all = permutation.dist$p.value.ci.all[-1] + gamma
    p.value = max(p.value.ci.all)
    
    return(list(p.value.ci = p.value,
                p.value.pi = p.value.ci.all[(grid.size+1)/2],
                p.value.truth = permutation.dist$p.value.ci.all[1]))
    
  } else {
    p.value.ci.all = permutation.dist$p.value.ci.all + gamma
    p.value = max(p.value.ci.all)
    
    return(list(p.value.ci = p.value,
                p.value.pi = p.value.ci.all[(grid.size+1)/2]))
  }
 
}

