#####################################

# Example of hierarchical metacognitive efficiency (Mratio) calculation 
# for two domains and correlation coefficient 
# exemple of trace plots and posterior distribution plots
# using the Function_metad_groupcorr.R
# The same function allows also the calculation for 3 and 4 domains
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

# Task 1
nR_S1_1 <- data.frame(
  p1 = c(52,32,35,37,26,12,4,2),
  p2 = c(27,39,46,52,14,10,9,3),
  p3 = c(112,30,15,17,17,9,0,0))
nR_S2_1 <- data.frame(
  p1 = c(2,5,15,22,33,38,40,45),
  p2 = c(2,4,9,21,57,48,34,25),
  p3 = c(0,1,7,18,12,17,27,118))

# Task 2
nR_S1_2 <- data.frame(
  p1 = c(97,49,13,9,20,11,1,0),
  p2 = c(37,41,49,44,17,11,0,1),
  p3 = c(61,45,34,28,21,9,1,1))
nR_S2_2 <- data.frame(
  p1 = c(0,1,8,23,17,33,22,96),
  p2 = c(0,2,9,18,44,46,43,38),
  p3 = c(2,5,3,22,32,38,27,71))

## Hierarchical meta_d correlation function ------------------------------------------------------

# List creation for model inputs
nR_S1 <- list(nR_S1_1,
              nR_S1_2)

nR_S2 <- list(nR_S2_1,
              nR_S2_2)

# Fit all data at once
source("fit_metad_groupcorr.R")
output <- fit_metad_groupcorr(nR_S1 = nR_S1, nR_S2 = nR_S2)

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

# Plot posterior distribution for rho value
Rho_plot <- mcmc.sample %>%
  filter(Parameter == "rho") %>% 
  ggplot(aes(value)) +
  geom_histogram(binwidth = 0.03, fill = "blue", colour = "grey", alpha = 0.5) +
  geom_vline(xintercept = stat$mean[stat$name == "rho"],linetype="dashed", linewidth = 1.5) +
  annotate("segment", x = HDI$lower[HDI$name == "rho"], y = 50, 
           xend = HDI$upper[HDI$name == "rho"], yend = 50, 
           colour = "white", linewidth = 2.5) +
  xlim(c(-1, 1)) +
  ylab("Sample count") +
  xlab(expression(paste(rho, " value")))

Rho_plot
