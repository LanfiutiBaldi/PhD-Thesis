#########
###GBR###
#########

## load useful packages
library(HMDHFDplus)
library(tidyverse)
library(viridis)
library(rstan)
library(BAPC)
library(Metrics)
library(brms)
library(ggpubr)

## load all the functions:
source("Functions.R")

# StanOptions
stanOptionsfit <- list(nchains = 4, niter = 2000, nwarmup = 1000, ncores = 4)

# Data
load("Data/GBR_data.RData")
GB_apc <- apc; rm(apc, data)
GB_bapc_data <- bapc_data; rm(bapc_data, bapc_deaths, bapc_pop)

# Fitting
GB_fit <- fit_APC_SkN(GB_apc, stanOptionsfit)
GB_effects <- get_effects(GB_fit)

# Parameter
GB_table_par <- table_par(GB_effects, "GBR") #Tab 2

## -- OUT-OF-SAMPLE TGBT PHASE -------------

## General settings

baseline1 <- 1960:1990
baseline2 <- 1970:2000
baseline3 <- 1980:2010

h=10 

stanOptions <-  list(nchains = 3, niter = 1000, nwarmup = 500, ncores = 3)

#FIRST WINDOW:

#Our Model
GB_fit_baseline1 <- fit_APC_SkN_forecast(GB_apc, baseline1, stanOptions)

GB_effects_baseline1 <- get_effects(GB_fit_baseline1)

GB_table_par1 <- table_par(GB_effects_baseline1, "GBR") #Tab 2

#ARIMA
GB_SR_forecast1 <- SR_forecasting(GB_apc, GB_effects_baseline1, baseline1, h,
                                  arimaPlot=FALSE, IC=TRUE)

#ARMA
GB_SR_forecast1_arma <- SR_forecasting_arma(GB_apc, 
                                            GB_effects_baseline1,
                                            baseline = baseline1,
                                            h = 10, IC = T, nsim = 1000)

## Poisson Approach
GB_bapc_forecast1 <- bapc_forecasting(GB_bapc_data, baseline1, h)

#Compare Approaches:
GB_compare1 <- GB_apc %>% 
  filter(P %in% max(baseline1+1):(max(baseline1)+h)) %>% 
  mutate("P"=as.factor(P),"A"=as.factor(A)) %>% 
  select(P,A,SR) %>% 
  left_join(GB_bapc_forecast1 %>% 
              mutate(P = as.factor(P), A = as.factor(A)))

GB_metric1 <- metrics(GB_SR_forecast1, GB_bapc_forecast1)
GB_metric1_arma <- metrics_arma(GB_SR_forecast1_arma, GB_bapc_forecast1)

#SECOND WINDOW:

#Our Model
GB_fit_baseline2 <- fit_APC_SkN_forecast(GB_apc, baseline2, stanOptions)

GB_effects_baseline2 <- get_effects(GB_fit_baseline2)

GB_table_par2 <- table_par(GB_effects_baseline2, "GBR") #Tab 2

#ARIMA
GB_SR_forecast2 <- SR_forecasting(GB_apc, GB_effects_baseline2, baseline2, h,
                                  arimaPlot=FALSE, IC=TRUE)

#ARMA
GB_SR_forecast2_arma <- SR_forecasting_arma(GB_apc, 
                                            GB_effects_baseline2,
                                            baseline = baseline2,
                                            h = 10, IC = T, nsim = 1000)

## Poisson Approach
GB_bapc_forecast2 <- bapc_forecasting(GB_bapc_data, baseline2, h)

#Compare Approaches:
GB_compare2 <- GB_apc %>% 
  filter(P %in% max(baseline2+1):(max(baseline2)+h)) %>% 
  mutate("P"=as.factor(P),"A"=as.factor(A)) %>% 
  select(P,A,SR) %>% 
  left_join(GB_bapc_forecast2 %>% 
              mutate(P = as.factor(P), A = as.factor(A)))

GB_metric2 <- metrics(GB_SR_forecast2, GB_bapc_forecast2)
GB_metric2_arma <- metrics_arma(GB_SR_forecast2_arma, GB_bapc_forecast2)

#THIRD WINDOW:

#Our Model
GB_fit_baseline3 <- fit_APC_SkN_forecast(GB_apc, baseline3, stanOptions)

GB_effects_baseline3 <- get_effects(GB_fit_baseline3)

GB_table_par3 <- table_par(GB_effects_baseline3, "GBR") #Tab 2

#ARIMA
GB_SR_forecast3 <- SR_forecasting(GB_apc, GB_effects_baseline3, baseline3, h,
                                  arimaPlot=FALSE, IC=TRUE)

#ARMA
GB_SR_forecast3_arma <- SR_forecasting_arma(GB_apc, 
                                            GB_effects_baseline3,
                                            baseline = baseline3,
                                            h = 10, IC = T, nsim = 1000)

## Poisson Approach
GB_bapc_forecast3 <- bapc_forecasting(GB_bapc_data, baseline3, h)

#Compare Approaches:
GB_compare3 <- GB_apc %>% 
  filter(P %in% max(baseline3+1):(max(baseline3)+h)) %>% 
  mutate("P"=as.factor(P),"A"=as.factor(A)) %>% 
  select(P,A,SR) %>% 
  left_join(GB_bapc_forecast3 %>% 
              mutate(P = as.factor(P), A = as.factor(A)))

GB_metric3 <- metrics(GB_SR_forecast3, GB_bapc_forecast3)
GB_metric3_arma <- metrics_arma(GB_SR_forecast3_arma, GB_bapc_forecast3)


#SUMMARY
GB_parameters <- rbind(data.frame(GB_table_par, Window="Full Period"),
                       data.frame(GB_table_par1, Window="1960-1990"),
                       data.frame(GB_table_par2, Window="1970-2000"),
                       data.frame(GB_table_par3, Window="1980-2010"))

GB_metrics_tab <-  rbind(data.frame(rbind(GB_metric1, GB_metric1_arma[-2,]), 
                                    Window="1960-1990"),
                         data.frame(rbind(GB_metric2, GB_metric2_arma[-2,]), 
                                    Window="1970-2000"),
                         data.frame(rbind(GB_metric3, GB_metric3_arma[-2,]), 
                                    Window="1980-2010"))