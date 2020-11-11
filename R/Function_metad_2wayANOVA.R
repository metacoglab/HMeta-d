#####################################

# Estimate metacognitive efficiency (Mratio) at the group level
#
# Adaptation in R of matlab function 'fit_meta_d_mcmc_groupCorr.m'
# by Steve Fleming 
# for more details see Fleming (2017). HMeta-d: hierarchical Bayesian 
# estimation of metacognitive efficiency from confidence ratings. 
#
# you need to install the following packing before using the function:
# tidyverse
# magrittr
# reshape2
# rjags
# coda
# lattice
# broom
# ggpubr
# ggmcmc
#
# nR_S1 and nR_S2 should be two lists of each nR_S1 or nR_S2 per task
# model output is a large mcmc list and two vectors for d1 and c1
#
# Author: Nicolas Legrand nicolas.legrand@cfin.au.dk

#####################################

## Packages
library(tidyverse)
library(magrittr)
library(reshape2)
library(rjags)
library(coda)
library(lattice)
library(broom)
library(ggpubr)
library(ggmcmc)

# Model -------------------------------------------------------------------

metad_2wayANOVA <- function (nR_S1_tot, nR_S2_tot) {
  # nR_S1_tot and nR_S2_tot should be lists of size N subjects.
  # Each list contain 4 responses vectors.
  # Design matrix: Condition 1: (1, 1, 0, 0) - Condition 2: (1, 0, 1, 0)

  # Type 1 parameters
  nratings <- length(nR_S1_tot[[1]][[1]])/2
  nsubj <- length((nR_S1_tot))
  
  d1 <- matrix(ncol = 4, nrow = nsubj)
  c1 <- matrix(ncol = 4, nrow = nsubj)
  counts_total = array(dim = c(nsubj, nratings*4, 4))
  
  for (n in 1:(nsubj)) {
    for (i in 1:4) {
      
      nR_S1 = nR_S1_tot[[n]][i]
      nR_S2 = nR_S2_tot[[n]][i]
      
      # Adjust to ensure non-zero counts for type 1 d' point estimate
      adj_f <- 1/((nratings)*2)
      nR_S1_adj = unlist(nR_S1) + adj_f
      nR_S2_adj = unlist(nR_S2) + adj_f
      
      ratingHR <- matrix()
      ratingFAR <- matrix()
      
      for (c in 2:(nratings*2)) {
        ratingHR[c-1] <- sum(nR_S2_adj[c:length(nR_S2_adj)]) / sum(nR_S2_adj)
        ratingFAR[c-1] <- sum(nR_S1_adj[c:length(nR_S1_adj)]) / sum(nR_S1_adj)
        
      }
      
      d1[n, i] = qnorm(ratingHR[(nratings)]) - qnorm(ratingFAR[(nratings)])
      c1[n, i] = -0.5 * (qnorm(ratingHR[(nratings)]) + qnorm(ratingFAR[(nratings)]))
      
      counts_total[n, , i] <- t(nR_S1[[1]]) %>% 
        cbind(t(nR_S2[[1]]))
    }
  }
  
  Tol <- 1e-05
  
  data <- list(
    d1 = d1,  # [Subjects * Condition]
    c1 = c1,  # [Subjects * Condition]
    nsubj = nsubj,
    counts = counts_total,  # [Subjects * Counts * Condition]
    nratings = nratings,
    Tol = Tol,
    Condition1 = c(1, 1, 0, 0),
    Condition2 = c(1, 0, 1, 0),
    Interaction = c(1, 0, 0, 0)
  )
  
  ## Model using JAGS
  # Create and update model
  aov_model <- jags.model("Bayes_metad_2wayANOVA.txt", data = data,
                          n.chains = 3, n.adapt= 2000, quiet=FALSE)
  update(aov_model, n.iter=5000)
  
  # Sampling
  output <- coda.samples( 
    model          = aov_model,
    variable.names = c("muBd_Condition1", "lamBd_Condition1", "sigD_Condition1", "muBd_Condition2", "lamBd_Condition2", "sigD_Condition2", 
                       "muBd_interaction", "lamBd_interaction", "sigD_interaction", "Mratio"),
    n.iter         = 10000,
    thin           = 1 )
  
  return(output)
}