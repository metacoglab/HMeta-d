#####################################

# Example of meta d calculation for individual subject and
# exemple of trace plots and posterior distribution plots
# using the Function_metad_indiv.R
# AM 2018

#####################################

## Packages ----------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(rjags)
library(coda)
library(ggmcmc)

## Create data for 1 participant -------------------------------------------------------------
nR_S1 <- c(52,32,35,37,26,12,4,2)
nR_S2 <- c(2,5,15,22,33,38,40,45)

## Individual meta_d function ------------------------------------------------------
source("fit_metad_indiv.R")
fit <- fit_metad_indiv(nR_S1 = nR_S1, nR_S2 = nR_S2)

## Model output ------------------------------------------------------------
output = fit[[1]]
d1 = fit[[2]]$d1

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

# Plot posterior distribution for meta-d value
Plot <- mcmc.sample %>%
  filter(Parameter == "meta_d") %>% 
  ggplot(aes(value)) +
  geom_histogram(binwidth = 0.03, fill = "blue", colour = "grey", alpha = 0.5) +
  geom_vline(xintercept = stat$mean[stat$name == "meta_d"],linetype="dashed", linewidth= 1.5) +
  annotate("segment", x = HDI$lower[HDI$name == "meta_d"], y = 50, 
           xend = HDI$upper[HDI$name == "meta_d"],
           yend = 50, colour = "white", linewidth = 2.5) +
  ylab("Sample count") +
  xlab(expression(paste("Meta d'")))

Plot
