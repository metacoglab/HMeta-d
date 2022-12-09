# HMeta-d for between-subjects regression on meta-d'/d'
#
#Adaptation in R of matlab function 'fit_metad_mcmc_regression.m'
#by Steve Fleming (2017) 
#
#
# you need to install the following packing before using the function:
# coda
# rjags
# magrittr
# dplyr
# tidyr
# tibble
# ggmcmc
#
# nR_S1 and nR_S2 should be two vectors
# cov is a n x s matrix of covariates, where s=number of subjects, n=number of covarient
# model output is a large mcmc list and two vectors for d1 and c1

#########################

#Packages
library(coda)
library(rjags)
library(magrittr)
library(dplyr)
library(tidyr)
library(tibble)
library(ggmcmc)

fit_meta_d_regression <- function (nR_S1, nR_S2) {
  
  #get type 1 SDT parameter values
  Nsubj <- ncol(nR_S1)
  nRatings = nrow(nR_S1)/2
  nTot <- sum(nR_S1[,1], nR_S2[,1]) #nb of trials
  
  nR_S1 <- list(nR_S1)
  nR_S2 <- list(nR_S2)

  d1 <- data.frame()
  c1 <- data.frame()
  
  for(n in 1:(Nsubj)) {
  # Adjust to ensure non-zero counts for type 1 d' point estimate 
  # (not necessary if estimating d' inside JAGS)
  adj_f <- 1/((nRatings)*2)
  nR_S1_adj = nR_S1[[1]][,n] + adj_f
  nR_S2_adj = nR_S2[[1]][,n] + adj_f
  
  ratingHR <- matrix()
  ratingFAR <- matrix()
  
  for (c in 2:(nRatings*2)) {
    ratingHR[c-1] <- sum(nR_S2_adj[c:length(nR_S2_adj)]) / sum(nR_S2_adj)
    ratingFAR[c-1] <- sum(nR_S1_adj[c:length(nR_S1_adj)]) / sum(nR_S1_adj)
    
  }
  
  t1_index <- nRatings

  a <- qnorm(ratingHR[(t1_index)]) - qnorm(ratingFAR[(t1_index)])
  d1 %<>%
    rbind(a)
  a <- -0.5 * (qnorm(ratingHR[(t1_index)]) + qnorm(ratingFAR[(t1_index)]))
  c1 %<>%
    rbind(a)
  } 
  
  d1 <- c(d1[1:Nsubj,1])
  c1 <- c(c1[1:Nsubj,1])
  cov <- c(cov[1:Nsubj,1])
  
  # Data preparation for model
  counts <- t(nR_S1[[1]]) %>%
    cbind(t(nR_S2[[1]]))

  Tol <- 1e-05

  data <- list(
    d1 = d1,
    c1 = c1,
    nsubj = Nsubj,
    counts = counts,
    cov = as.vector(cov),
    nratings = nRatings,
    Tol = Tol)
}

## Model using JAGS
# Create and update model
regression_model <- jags.model(file = 'Bayes_metad_group_regress_nodp.txt', data = data,
                               n.chains = 3, quiet=FALSE)
update(regression_model, n.iter=1000)

# Sampling
output <- coda.samples( 
  model          = regression_model,
  variable.names = c("mu_logMratio", "sigma_logMratio", "mu_c2", "mu_beta1", "Mratio"),
  n.iter         = 10000,
  thin           = 1 )

## Model output ------------------------------------------------------------

# Convergence diagnostic
Value <- gelman.diag(output, confidence = 0.95)
Rhat <- data.frame(conv = Value$psrf)

# Values (mean and CI)
Value <- summary(output)
stat <- data.frame(mean = Value$statistics[,"Mean"])
stat %<>%
  rownames_to_column(var = "name") %>% 
  cbind(CILow = Value$quantiles[,"2.5%"]) %>% 
  cbind(CIUp = Value$quantiles[,"97.5%"])

# HDI function
HDI <- data.frame(HPDinterval(output, prob = 0.95))
HDI %<>%
  rownames_to_column(var = "name")

## Plot trace mcmc ---------------------------------------------------------
traceplot(output)

## Plot posterior distributions --------------------------------------------
mcmc.sample <- ggs(output)

##Save model
Fit <- stat %>%
  cbind(lower = HDI$lower,
        upper = HDI$upper,
        Rhat = Rhat[,1])
