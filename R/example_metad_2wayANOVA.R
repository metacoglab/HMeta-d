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
  
  
  #########################################################################
  # Create data for 5 participants ----------------------------------------
  # Simulate responses from 5 participants in a 2x2 repeated measures design
  # Condition 1: (1, 1, 0, 0) - Condition 2: (1, 0, 1, 0)
  
  nR_S1 = list(
  list(c(23.,  9.,  7.,  5.,  3.,  3.,  0.,  0.), c(10., 10.,  8., 12.,  6.,  2.,  2.,  0.),
       c(11.,  8., 11.,  8.,  7.,  3.,  1.,  1.), c(20.,  7.,  7.,  6.,  6.,  2.,  1.,  1.)),
  list(c(16., 13., 10.,  4.,  7.,  0.,  0.,  0.), c(11.,  9., 12., 13.,  3.,  2.,  0.,  0.),
       c(16., 14.,  6.,  3.,  7.,  3.,  1.,  0.), c(22.,  8.,  5.,  6.,  5.,  3.,  1.,  0.)),
  list(c(16.,  7.,  8.,  9.,  7.,  2.,  1.,  0.), c(14.,  8., 10.,  9.,  6.,  2.,  1.,  0.),
       c(17., 12.,  6.,  7.,  3.,  4.,  0.,  1.), c(21.,  6.,  8.,  9.,  4.,  1.,  1.,  0.)),
  list(c(21.,  8.,  9.,  4.,  4.,  3.,  1.,  0.), c(11., 16.,  6.,  9.,  4.,  3.,  0.,  1.),
       c(10.,  8.,  9., 11.,  9.,  2.,  1.,  0.), c(15.,  7., 13.,  7.,  4.,  4.,  0.,  0.)),
  list(c(15., 12.,  7.,  8.,  4.,  2.,  2.,  0.), c(5., 16.,  8., 11.,  4.,  4.,  2.,  0.),
       c(9., 10., 12.,  9.,  8.,  1.,  0.,  1.), c(25.,  7.,  5.,  6.,  6.,  1.,  0.,  0.))
  )
  
  nR_S2 = list(
    list(c(0.,  1.,  4.,  6.,  9.,  6., 10., 14.), c(0.,  1.,  2.,  3.,  8., 10.,  9., 17.),
         c(0.,  2.,  3.,  2.,  7., 10.,  6., 20.), c(0.,  1.,  2.,  3.,  6.,  8.,  9., 21.)),
    list(c(0.,  1.,  1.,  6.,  9.,  8., 10., 15.), c(0.,  2.,  3.,  5., 11., 15.,  9.,  5.),
         c(0.,  1.,  2.,  4.,  3.,  8.,  8., 24.), c(0.,  0.,  1.,  5.,  2., 10.,  8., 24.)),
    list(c(0.,  1.,  3.,  2., 12.,  4.,  9., 19.), c(1.,  1.,  0.,  5., 12., 14.,  6., 11.),
         c(0.,  1.,  4.,  4.,  8.,  6., 11., 16.), c(1.,  2.,  1.,  6.,  5.,  8., 10., 17.)),
    list(c(0.,  0.,  3.,  5.,  3., 11.,  6., 22.), c(1.,  0.,  3.,  4., 13.,  7., 13.,  9.),
         c(0.,  0.,  4.,  1.,  4., 12., 12., 17.), c(1.,  2.,  0.,  5.,  5., 14.,  8., 15.)),
    list(c(0.,  1.,  3.,  4.,  6.,  7., 16., 13.), c(0.,  0.,  1.,  5.,  8., 10., 12., 14.),
         c(0.,  0.,  3.,  2.,  2.,  7., 15., 21.), c(0.,  2.,  1.,  6., 10., 10., 10., 11.))
  )
  
  # Model -------------------------------------------------------------------
  
  # Fit all data at once
  source("fit_metad_2wayANOVA.R")
  output <- metad_2wayANOVA(nR_S1 = nR_S1, nR_S2 = nR_S2)
  
  ## Summary stats --------------------------------------------------
  
  # Mean values 
  Value <- summary(output)
  stat <- data.frame(mean = Value[["statistics"]][, "Mean"])
  stat %<>%
    rownames_to_column(var = "name")
  
  # Rhat 
  Value <- gelman.diag(output, confidence = 0.95)
  Rhat <- data.frame(conv = Value$psrf)
  
  # HDI 
  HDI <- data.frame(HPDinterval(output, prob = 0.95))
  HDI %<>%
    rownames_to_column(var = "name")
  
  # All values in the same dataframe
  Fit <- stat %>%
    cbind(lower = HDI$lower,
          upper = HDI$upper)
  
  ## Plots ---------------------------------------------------------
  
  # Plot trace mcmc
  traceplot(output)
  
  # mcmc values in df for plot posterior distributions
  mcmc.sample <- ggs(output)
  
  # Plot posterior distribution for rho value
  PlotCondition1 <- mcmc.sample %>%
    filter(Parameter == "muBd_Condition1") %>% 
    ggplot(aes(value)) +
    geom_histogram(fill = "blue", colour = "grey", alpha = 0.5, bins = 100) +
    geom_vline(xintercept = stat$mean[stat$name == "muBd_Condition1"],linetype="dashed", linewidth = 1.) +
    geom_vline(xintercept = 0, linewidth = 1.5, color='red') +
    annotate("segment", x= HDI$lower[HDI$name == "muBd_Condition1"], 
             y = 50, xend = HDI$upper[HDI$name == "muBd_Condition1"],
             yend = 50, colour = "white", linewidth = 2.5) +
    ylab("Sample count") +
    xlab(expression(paste("muBd_Condition1")))
  PlotCondition1
  
  # Plot posterior distribution for rho value
  PlotCondition2 <- mcmc.sample %>%
    filter(Parameter == "muBd_Condition2") %>% 
    ggplot(aes(value)) +
    geom_histogram(fill = "blue", colour = "grey", alpha = 0.5, bins = 100) +
    geom_vline(xintercept = stat$mean[stat$name == "muBd_Condition2"],linetype="dashed", linewidth = 1.) +
    geom_vline(xintercept = 0, linewidth = 1.5, color='red') +
    annotate("segment", x = HDI$lower[HDI$name == "muBd_Condition2"],
             y = 50, xend = HDI$upper[HDI$name == "muBd_Condition2"],
             yend = 50, colour = "white", linewidth = 2.5) +
    ylab("Sample count") +
    xlab(expression(paste("muBd_Condition2")))
  PlotCondition2
  
  # Plot posterior distribution for rho value
  PlotInteraction <- mcmc.sample %>%
    filter(Parameter == "muBd_interaction") %>% 
    ggplot(aes(value)) +
    geom_histogram(fill = "blue", colour = "grey", alpha = 0.5, bins = 100) +
    geom_vline(xintercept = stat$mean[stat$name == "muBd_interaction"],linetype="dashed", linewidth = 1.) +
    geom_vline(xintercept = 0, linewidth = 1.5, color='red') +
    annotate("segment", x = HDI$lower[HDI$name == "muBd_interaction"], 
             y = 50, xend = HDI$upper[HDI$name == "muBd_interaction"],
             yend = 50, colour = "white", linewidth = 2.5) +
    ylab("Sample count") +
    xlab(expression(paste("muBd_interaction")))
  PlotInteraction
  
  # Save samples ------------------------------------------------------------
  
  df = rbind(data.frame(output[[1]]), data.frame(output[[2]]), data.frame(output[[3]]))
  
  write.table(df, file = 'posterior.txt', append = FALSE, sep = "\t", dec = ".",
              row.names = TRUE, col.names = TRUE)
