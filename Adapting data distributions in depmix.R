# Code to try adapting data distributions fed to depmixS4
# Run in code depmixS4.R first to obtain data

library(depmixS4)
library(EnvStats)
library(ggplot2)
library(ddply)
library(dplyr)

dat <- HMMdat2
  
# Find data distributions
# Step length
gamma <- egamma(dat$step)
shape <- as.numeric(gamma$parameters[1])
scale <- as.numeric(gamma$parameters[2])

# Straightness Index
beta <- ebeta(dat$SI)
a <- as.numeric(beta$parameters[1])
b <- as.numeric(beta$parameters[2])

# Initialise model
HMMmod <- depmix(list(step ~ 1,SI ~ 1), data = dat, 
                 family = list(Gamma(), gaussian()), respstart = c(shape, scale, a, b, 1, 1), 
                 nstates = 2, ntimes = tracklengths$Freq)

# Fit model
depmod <- fit(HMMmod)
depmod

dat$state <- as.factor(depmixS4::viterbi(fitHMMmod)$state)
dat$prob1 <- depmixS4::viterbi(fitHMMmod)$S1
dat$prob2 <- depmixS4::viterbi(fitHMMmod)$S2

