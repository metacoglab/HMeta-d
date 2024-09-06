#####################################

# Example of hierarchical metacognitive efficiency (Mratio) 
# at the group level
# exemple of trace plots and posterior distribution plots
# using the Function_metad_group.R
#
# AM 2018

#####################################

## Packages ----------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(reshape2)
library(rjags)
library(coda)
library(lattice)
library(broom)
library(ggpubr)
library(ggmcmc)

## Create data for 3 participants -------------------------------------------------------------

nR_S1 <- data.frame(
  p1 = c(52,32,35,37,26,12,4,2),
  p2 = c(27,39,46,52,14,10,9,3),
  p3 = c(112,30,15,17,17,9,0,0))
nR_S2 <- data.frame(
  p1 = c(2,5,15,22,33,38,40,45),
  p2 = c(2,4,9,21,57,48,34,25),
  p3 = c(0,1,7,18,12,17,27,118))

## Hierarchical meta_d group function ------------------------------------------------------

# List creation for model inputs
nR_S1 <- list(nR_S1)
nR_S2 <- list(nR_S2)

# Fit all data at once
source("fit_metad_group.R")
output <- fit_metad_group(nR_S1 = nR_S1, nR_S2 = nR_S2)

## Model output ------------------------------------------------------------

# Values 
Value <- summary(output)
stat <- data.frame(mean = Value$statistics[,"Mean"])
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
        upper = HDI$upper,
        Rhat = Rhat[,1])

## Plots ---------------------------------------------------------

# Plot trace mcmc
traceplot(output)

# mcmc values in df for plot posterior distributions
mcmc.sample <- ggs(output)

# Plot posterior distribution for mu Mratio value
Mratio_plot <- mcmc.sample %>%
  filter(Parameter == "mu_logMratio") %>% 
  ggplot(aes(exp(value))) +
  geom_histogram(binwidth = 0.03, fill = "blue", colour = "grey", alpha = 0.5) +
  geom_vline(xintercept = exp(stat$mean[stat$name == "mu_logMratio"]), linetype = "dashed", linewidth = 1.5) +
  annotate("segment", x = exp(HDI$lower[HDI$name == "mu_logMratio"]), y = 50, 
                   xend = exp(HDI$upper[HDI$name == "mu_logMratio"]), yend = 50, 
               colour = "white", linewidth = 2.5) +
  ylab("Sample count") +
  xlab(expression(paste(mu, " Mratio")))


Mratio_plot
