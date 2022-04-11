#####################################

# Estimate metacognitive efficiency (Mratio) at the group level
#
# Adaptation in R of matlab function 'fit_meta_d_mcmc_group.m'
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
# nR_S1 and nR_S2 should be two vectors
# model output is a large mcmc list and two vectors for d1 and c1
#
# Audrey Mazancieux 2018

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

fit_metad_group <- function (nR_S1, nR_S2) {
  
  # Type 1 parameters
  nTot <- sum(nR_S1[[1]]$V1, nR_S2[[1]]$V1)
  nratings <- nrow(nR_S1[[1]])/2
  nsubj <- ncol(nR_S1[[1]])
  nTask <- length(nR_S1)

    
  # Adjust to ensure non-zero counts for type 1 d' point estimate
  d1 <- data.frame()
  c1 <- data.frame()
    
    for (n in 1:(nsubj)) {
        
      adj_f <- 1/((nratings)*2)
      nR_S1_adj = nR_S1[[1]][,n] + adj_f
      nR_S2_adj = nR_S2[[1]][,n] + adj_f
        
      ratingHR <- matrix()
      ratingFAR <- matrix()
        
      for (c in 2:(nratings*2)) {
        ratingHR[c-1] <- sum(nR_S2_adj[c:length(nR_S2_adj)]) / sum(nR_S2_adj)
        ratingFAR[c-1] <- sum(nR_S1_adj[c:length(nR_S1_adj)]) / sum(nR_S1_adj)
          
      }
        
      t1_index <- nratings
      a <- qnorm(ratingHR[(t1_index)]) - qnorm(ratingFAR[(t1_index)])
      d1 %<>%
        rbind(a)
      a <- -0.5 * (qnorm(ratingHR[(t1_index)]) + qnorm(ratingFAR[(t1_index)]))
      c1 %<>%
        rbind(a)
    }

  d1 <- c(d1[1:nsubj,1])
  c1 <- c(c1[1:nsubj,1])
    
  # Data preparation for model
  counts <- t(nR_S1[[1]]) %>% 
    cbind(t(nR_S2[[1]]))
  
  d1 <<- as.matrix(d1)
  c1 <<- as.matrix(c1)

    Tol <- 1e-05
    
    data <- list(
      d1 = d1,
      c1 = c1,
      nsubj = nsubj,
      counts = counts,
      nratings = nratings,
      Tol = Tol
    )
    
    ## Model using JAGS
    # Create and update model
    model <- jags.model(file = 'Bayes_metad_group_R.txt', data = data,
                            n.chains = 3, quiet=FALSE)
    update(model, n.iter=1000)
    
    # Sampling
    output <- coda.samples( 
      model          = model,
      variable.names = c("mu_logMratio", "sigma_logMratio", "Mratio", "mu_c2"),
      n.iter         = 10000,
      thin           = 1 )
    
  return(output)
}
