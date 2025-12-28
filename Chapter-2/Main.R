#Main code for the paper: An age-period modeling of population gaps in mortality: the cases of cancer and cardiovascular diseases

#Necessary Packages:
library(Rmpfr)
library(Bessel) 
library(skellam)
library(tidyverse)
library(plotly)
library(brms)
library(Metrics)
library(HMDHFDplus)
library(beepr)
library(readxl)
library(ggpubr)

#Import all functions
source("Functions.R")

#Import Data 

#Dati Italia WHO
load("Data/ready_to_use_Italy.RData")

#Dati USA HMD,
#load("Data/ready_to_use_USA.RData")

dt <- ITA %>%
  filter(Year %in% 1961:2015) %>%
  grouping_01ages() 

plot_initial(dt)

SIM <- prepare_data(dt)

##OUT-OF-SAMPLE EXERCISE
h <- 15 #Number of years for forecasting
baseline <- 1961:2000 
forec_year <- (max(baseline)+1):(max(baseline)+h)
nages <- nlevels(SIM$ages)

dt_base <- dt %>% 
  filter(Year %in% baseline) 

SIM_base <- prepare_data(dt_base)

##RUN ALL THE MODELS
#Skellam regression (WARNING: it might take a while (~ 60 min))
sk_reg_base <- Skellam_regression(SIM)

#Double Poisson regression
dp_reg_base <- DoublePoisson_regression(SIM)

#Bivariate Poisson regression 
bp_reg_base <- BivariatePoisson_regression(SIM)

##GET AGE-PERIOD EFFECTS
sk_eff_base <- get_effects(sk_reg_base, SIM_base)
dp_eff_base <- get_effects(dp_reg_base, SIM_base)
bp_eff_base <- get_effects(bp_reg_base, SIM_base)

All_eff_base <- join_eff(sk_eff_base, dp_eff_base, bp_eff_base, SIM_base)

#Age-Period Effects estimated in the baseline period
plot_eff(All_eff_base)

#Forecasting Period Effects using Multivariate Random Walk with Drift
year_eff_forec <- forecasting_yearEffects(All_eff_base, h)

SIM_for <- data.frame("years" = as.factor(rep(forec_year, each=nages)),
                      "ages" = as.factor(rep(unique(SIM_base$ages), h)),
                      "diff" = 0) 

SIM_tot <- rbind(SIM_base[,c(1,2,5)], SIM_for)

#Gap estimation (in-and-out-of sample)
diff_inout <- forecasting_differences(SIM_tot, 
                                      All_eff_base, 
                                      year_eff_forec, 
                                      h, nages)

#BIC, AIC and AICc
InfoCriterion <- rbind("Skellam"=data.frame(sk_reg_base$Metrics),
                       "Double Poisson"=data.frame(dp_reg_base$Metrics),
                       "Bivariate Poisson"=data.frame(bp_reg_base$Metrics))

#Check the estimation accuracy in the in-sample part 
diff_in <- diff_inout$in_AND_out_of_Sample %>% 
  filter(Year %in% baseline)

dt_check_in <- dt %>% 
  filter(Year %in% baseline) %>% 
  select(Year, Age, "ObservedGap"=Diff) %>% 
  ungroup()

check_in <- dt_check_in %>% 
  left_join(diff_in)

accuracymetrics_in <- metrics(check_in)

#Check the estimation accuracy in the out-of-sample part 
diff_forec <- diff_inout$out_of_sample

dt_check_for <- dt %>% 
  filter(Year %in% forec_year) %>% 
  select(Year, Age, "ObservedGap"=Diff) %>% 
  ungroup()

check_for <- dt_check_for %>% 
  left_join(diff_forec)

accuracymetrics_for <- metrics(check_for)

accuracy_heatmap(check_for, "RMSE")

#In-sample and out-of-sample fitting for specific ages
check_tot <- dt %>% 
  select(Year, Age, "ObservedGap"=Diff) %>% 
  ungroup() %>% 
  left_join(diff_inout$in_AND_out_of_Sample)

#Diebold-Mariano Test Skellam vs Double and Bivariate Poisson

DM_result <- DM_test(check_for)

