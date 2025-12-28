#########
###ITA###
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
load("Data/Italy_data.RData")
IT_apc <- apc; rm(apc, data)
IT_bapc_data <- bapc_data; rm(bapc_data, bapc_deaths, bapc_pop)

# Fitting
IT_fit <- fit_APC_SkN(IT_apc, stanOptionsfit)
IT_effects <- get_effects(IT_fit)

# Parameter
IT_table_par <- table_par(IT_effects, "ITA") #Tab 2

## -- OUT-OF-SAMPLE TEST PHASE -------------

## General settings

baseline1 <- 1960:1990
baseline2 <- 1970:2000
baseline3 <- 1980:2010

h=10 

stanOptions <-  list(nchains = 3, niter = 1000, nwarmup = 500, ncores = 3)

#FIRST WINDOW:

#Our Model
IT_fit_baseline1 <- fit_APC_SkN_forecast(IT_apc, baseline1, stanOptions)

IT_effects_baseline1 <- get_effects(IT_fit_baseline1)

IT_table_par1 <- table_par(IT_effects_baseline1, "IT") #Tab 2

#ARIMA
IT_SR_forecast1 <- SR_forecasting(IT_apc, IT_effects_baseline1, baseline1, h,
                                  arimaPlot=FALSE, IC=TRUE)


#ARMA
IT_SR_forecast1_arma <- SR_forecasting_arma(IT_apc, 
                                            IT_effects_baseline1,
                                            baseline = baseline1,
                                            h = 10, IC = T, nsim = 1000)

## Poisson Approach
IT_bapc_forecast1 <- bapc_forecasting(IT_bapc_data, baseline1, h)

#Compare Approaches:
IT_compare1 <- IT_apc %>% 
  filter(P %in% max(baseline1+1):(max(baseline1)+h)) %>% 
  mutate("P"=as.factor(P),"A"=as.factor(A)) %>% 
  select(P,A,SR) %>% 
  left_join(IT_bapc_forecast1 %>% 
              mutate(P = as.factor(P), A = as.factor(A)))

IT_metric1 <- metrics(IT_SR_forecast1, IT_bapc_forecast1)
IT_metric1_arma <- metrics(IT_SR_forecast1_arma, IT_bapc_forecast1)

#SECOND WINDOW:

#Our Model
IT_fit_baseline2 <- fit_APC_SkN_forecast(IT_apc, baseline2, stanOptions)

IT_effects_baseline2 <- get_effects(IT_fit_baseline2)

IT_table_par2 <- table_par(IT_effects_baseline2, "IT") #Tab 2

#AIRMA
IT_SR_forecast2 <- SR_forecasting(IT_apc, IT_effects_baseline2, baseline2, h,
                                  arimaPlot=FALSE, IC=TRUE)

#ARMA
IT_SR_forecast2_arma <- SR_forecasting_arma(IT_apc, 
                                            IT_effects_baseline2,
                                            baseline = baseline2,
                                            h = 10, IC = T, nsim = 1000)

## Poisson Approach
IT_bapc_forecast2 <- bapc_forecasting(IT_bapc_data, baseline2, h)

#Compare Approaches:
IT_compare2 <- IT_apc %>% 
  filter(P %in% max(baseline2+1):(max(baseline2)+h)) %>% 
  mutate("P"=as.factor(P),"A"=as.factor(A)) %>% 
  select(P,A,SR) %>% 
  left_join(IT_bapc_forecast2 %>% 
              mutate(P = as.factor(P), A = as.factor(A)))

IT_metric2 <- metrics(IT_SR_forecast2, IT_bapc_forecast2)
IT_metric2_arma <- metrics(IT_SR_forecast2_arma, IT_bapc_forecast2)

#THIRD WINDOW:

#Our Model
IT_fit_baseline3 <- fit_APC_SkN_forecast(IT_apc, baseline3, stanOptions)

IT_effects_baseline3 <- get_effects(IT_fit_baseline3)

IT_table_par3 <- table_par(IT_effects_baseline3, "IT") #Tab 2

#ARIMA
IT_SR_forecast3 <- SR_forecasting(IT_apc, IT_effects_baseline3, baseline3, h,
                                  arimaPlot=FALSE, IC=TRUE)

#ARMA 
IT_SR_forecast3_arma <- SR_forecasting_arma(IT_apc, 
                                            IT_effects_baseline3,
                                            baseline = baseline3,
                                            h = 10, IC = T, nsim = 1000)

## Poisson Approach
IT_bapc_forecast3 <- bapc_forecasting(IT_bapc_data, baseline3, h)

#Compare Approaches:
IT_compare3 <- IT_apc %>% 
  filter(P %in% max(baseline3+1):(max(baseline3)+h)) %>% 
  mutate("P"=as.factor(P),"A"=as.factor(A)) %>% 
  select(P,A,SR) %>% 
  left_join(IT_bapc_forecast3 %>% 
              mutate(P = as.factor(P), A = as.factor(A)))

IT_metric3 <- metrics(IT_SR_forecast3, IT_bapc_forecast3)
IT_metric3_arma <- metrics(IT_SR_forecast3_arma, IT_bapc_forecast3)



#SUMMARY
IT_parameters <- rbind(data.frame(IT_table_par, Window="Full Period"),
                   data.frame(IT_table_par1, Window="1960-1990"),
                   data.frame(IT_table_par2, Window="1970-2000"),
                   data.frame(IT_table_par3, Window="1980-2010"))


IT_metrics_tab <-  rbind(data.frame(rbind(IT_metric1, IT_metric1_arma[-2,]), 
                                   Window="1960-1990"),
                        data.frame(rbind(IT_metric2, IT_metric2_arma[-2,]), 
                                   Window="1970-2000"),
                        data.frame(rbind(IT_metric3, IT_metric3_arma[-2,]), 
                                   Window="1980-2010"))

