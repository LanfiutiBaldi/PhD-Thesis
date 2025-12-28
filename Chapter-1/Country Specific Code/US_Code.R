#########
###USA###
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
load("Data/USA_data.RData")
US_apc <- apc; rm(apc, data)
US_bapc_data <- bapc_data; rm(bapc_data, bapc_deaths, bapc_pop)

# Fitting
US_fit <- fit_APC_SkN(US_apc, stanOptionsfit)
US_effects <- get_effects(US_fit)

# Parameter
US_table_par <- table_par(US_effects, "USA") #Tab 2

## -- OUT-OF-SAMPLE TEST PHASE -------------

## General settings

baseline1 <- 1960:1990
baseline2 <- 1970:2000
baseline3 <- 1980:2010

h=10 

stanOptions <-  list(nchains = 3, niter = 1000, nwarmup = 500, ncores = 3)

#FIRST WINDOW:

#Our Model
US_fit_baseline1 <- fit_APC_SkN_forecast(US_apc, baseline1, stanOptions)

US_effects_baseline1 <- get_effects(US_fit_baseline1)

US_table_par1 <- table_par(US_effects_baseline1, "US") #Tab 2

#ARIMA
US_SR_forecast1 <- SR_forecasting(US_apc, US_effects_baseline1, baseline1, h,
                                  arimaPlot=FALSE, IC=TRUE)

#ARMA
US_SR_forecast1_arma <- SR_forecasting_arma(US_apc, 
                                            US_effects_baseline1,
                                            baseline = baseline1,
                                            h = 10, IC = T, nsim = 1000)

## Poisson Approach
US_bapc_forecast1 <- bapc_forecasting(US_bapc_data, baseline1, h)

#Compare Approaches:
US_compare1 <- US_apc %>% 
  filter(P %in% max(baseline1+1):(max(baseline1)+h)) %>% 
  mutate("P"=as.factor(P),"A"=as.factor(A)) %>% 
  select(P,A,SR) %>% 
  left_join(US_bapc_forecast1 %>% 
              mutate(P = as.factor(P), A = as.factor(A)))

US_metric1 <- metrics(US_SR_forecast1, US_bapc_forecast1)
US_metric1_arma <- metrics_arma(US_SR_forecast1_arma, US_bapc_forecast1)

#SECOND WINDOW:

#Our Model
US_fit_baseline2 <- fit_APC_SkN_forecast(US_apc, baseline2, stanOptions)

US_effects_baseline2 <- get_effects(US_fit_baseline2)

US_table_par2 <- table_par(US_effects_baseline2, "US") #Tab 2

#ARIMA
US_SR_forecast2 <- SR_forecasting(US_apc, US_effects_baseline2, baseline2, h,
                                  arimaPlot=FALSE, IC=TRUE)

#ARMA
US_SR_forecast2_arma <- SR_forecasting_arma(US_apc, 
                                            US_effects_baseline2,
                                            baseline = baseline2,
                                            h = 10, IC = T, nsim = 1000)

## Poisson Approach
US_bapc_forecast2 <- bapc_forecasting(US_bapc_data, baseline2, h)

#Compare Approaches:
US_compare2 <- US_apc %>% 
  filter(P %in% max(baseline2+1):(max(baseline2)+h)) %>% 
  mutate("P"=as.factor(P),"A"=as.factor(A)) %>% 
  select(P,A,SR) %>% 
  left_join(US_bapc_forecast2 %>% 
              mutate(P = as.factor(P), A = as.factor(A)))

US_metric2 <- metrics(US_SR_forecast2, US_bapc_forecast2)
US_metric2_arma <- metrics_arma(US_SR_forecast2_arma, US_bapc_forecast2)

#THIRD WINDOW:

#Our Model
US_fit_baseline3 <- fit_APC_SkN_forecast(US_apc, baseline3, stanOptions)

US_effects_baseline3 <- get_effects(US_fit_baseline3)

US_table_par3 <- table_par(US_effects_baseline3, "US") #Tab 2

#ARIMA
US_SR_forecast3 <- SR_forecasting(US_apc, US_effects_baseline3, baseline3, h,
                                  arimaPlot=FALSE, IC=TRUE)

#ARMA
US_SR_forecast3_arma <- SR_forecasting_arma(US_apc, 
                                            US_effects_baseline3,
                                            baseline = baseline3,
                                            h = 10, IC = T, nsim = 1000)

## Poisson Approach
US_bapc_forecast3 <- bapc_forecasting(US_bapc_data, baseline3, h)

#Compare Approaches:
US_compare3 <- US_apc %>% 
  filter(P %in% max(baseline3+1):(max(baseline3)+h)) %>% 
  mutate("P"=as.factor(P),"A"=as.factor(A)) %>% 
  select(P,A,SR) %>% 
  left_join(US_bapc_forecast3 %>% 
              mutate(P = as.factor(P), A = as.factor(A)))

US_metric3 <- metrics(US_SR_forecast3, US_bapc_forecast3)
US_metric3_arma <- metrics_arma(US_SR_forecast3_arma, US_bapc_forecast3)


#SUMMARY
US_parameters <- rbind(data.frame(US_table_par, Window="Full Period"),
                   data.frame(US_table_par1, Window="1960-1990"),
                   data.frame(US_table_par2, Window="1970-2000"),
                   data.frame(US_table_par3, Window="1980-2010"))

US_metrics_tab <-  rbind(data.frame(rbind(US_metric1, US_metric1_arma[-2,]), 
                                   Window="1960-1990"),
                        data.frame(rbind(US_metric2, US_metric2_arma[-2,]), 
                                   Window="1970-2000"),
                        data.frame(rbind(US_metric3, US_metric3_arma[-2,]), 
                                   Window="1980-2010"))
