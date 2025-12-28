#########
###FRA###
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
load("Data/FRA_data.RData")
FR_apc <- apc; rm(apc, data)
FR_bapc_data <- bapc_data; rm(bapc_data, bapc_deaths, bapc_pop)

# Fitting
FR_fit <- fit_APC_SkN(FR_apc, stanOptionsfit)
FR_effects <- get_effects(FR_fit)

# Parameter
FR_table_par <- table_par(FR_effects, "FRA") #Tab 2

## -- OUT-OF-SAMPLE TFRT PHASE -------------

## General settings

baseline1 <- 1960:1990
baseline2 <- 1970:2000
baseline3 <- 1980:2010

h=10 

stanOptions <-  list(nchains = 3, niter = 1000, nwarmup = 500, ncores = 3)

#FIRST WINDOW:

#Our Model
FR_fit_baseline1 <- fit_APC_SkN_forecast(FR_apc, baseline1, stanOptions)

FR_effects_baseline1 <- get_effects(FR_fit_baseline1)

FR_table_par1 <- table_par(FR_effects_baseline1, "FRA") #Tab 2

#ARIMA
FR_SR_forecast1 <- SR_forecasting(FR_apc, FR_effects_baseline1, baseline1, h,
                                  arimaPlot=FALSE, IC=TRUE)

#ARMA
FR_SR_forecast1_arma <- SR_forecasting_arma(FR_apc, 
                                            FR_effects_baseline1,
                                            baseline = baseline1,
                                            h = 10, IC = T, nsim = 1000)

## Poisson Approach
FR_bapc_forecast1 <- bapc_forecasting(FR_bapc_data, baseline1, h)

#Compare Approaches:
FR_compare1 <- FR_apc %>% 
  filter(P %in% max(baseline1+1):(max(baseline1)+h)) %>% 
  mutate("P"=as.factor(P),"A"=as.factor(A)) %>% 
  select(P,A,SR) %>% 
  left_join(FR_bapc_forecast1 %>% 
              mutate(P = as.factor(P), A = as.factor(A)))

FR_metric1 <- metrics(FR_SR_forecast1, FR_bapc_forecast1)
FR_metric1_arma <- metrics_arma(FR_SR_forecast1_arma, FR_bapc_forecast1)

#SECOND WINDOW:

#Our Model
FR_fit_baseline2 <- fit_APC_SkN_forecast(FR_apc, baseline2, stanOptions)

FR_effects_baseline2 <- get_effects(FR_fit_baseline2)

FR_table_par2 <- table_par(FR_effects_baseline2, "FRA") #Tab 2

#ARIMA
FR_SR_forecast2 <- SR_forecasting(FR_apc, FR_effects_baseline2, baseline2, h,
                                  arimaPlot=FALSE, IC=TRUE)

#ARMA
FR_SR_forecast2_arma <- SR_forecasting_arma(FR_apc, 
                                            FR_effects_baseline2,
                                            baseline = baseline2,
                                            h = 10, IC = T, nsim = 1000)

## Poisson Approach
FR_bapc_forecast2 <- bapc_forecasting(FR_bapc_data, baseline2, h)

#Compare Approaches:
FR_compare2 <- FR_apc %>% 
  filter(P %in% max(baseline2+1):(max(baseline2)+h)) %>% 
  mutate("P"=as.factor(P),"A"=as.factor(A)) %>% 
  select(P,A,SR) %>% 
  left_join(FR_bapc_forecast2 %>% 
              mutate(P = as.factor(P), A = as.factor(A)))

FR_metric2 <- metrics(FR_SR_forecast2, FR_bapc_forecast2)
FR_metric2_arma <- metrics_arma(FR_SR_forecast2_arma, FR_bapc_forecast2)

#THIRD WINDOW:

#Our Model
FR_fit_baseline3 <- fit_APC_SkN_forecast(FR_apc, baseline3, stanOptions)

FR_effects_baseline3 <- get_effects(FR_fit_baseline3)

FR_table_par3 <- table_par(FR_effects_baseline3, "FRA") #Tab 2

#ARIMA
FR_SR_forecast3 <- SR_forecasting(FR_apc, FR_effects_baseline3, baseline3, h,
                                  arimaPlot=FALSE, IC=TRUE)

#ARMA
FR_SR_forecast3_arma <- SR_forecasting_arma(FR_apc, 
                                            FR_effects_baseline3,
                                            baseline = baseline3,
                                            h = 10, IC = T, nsim = 1000)

## Poisson Approach
FR_bapc_forecast3 <- bapc_forecasting(FR_bapc_data, baseline3, h)

#Compare Approaches:
FR_compare3 <- FR_apc %>% 
  filter(P %in% max(baseline3+1):(max(baseline3)+h)) %>% 
  mutate("P"=as.factor(P),"A"=as.factor(A)) %>% 
  select(P,A,SR) %>% 
  left_join(FR_bapc_forecast3 %>% 
              mutate(P = as.factor(P), A = as.factor(A)))

FR_metric3 <- metrics(FR_SR_forecast3, FR_bapc_forecast3)
FR_metric3_arma <- metrics_arma(FR_SR_forecast3_arma, FR_bapc_forecast3)


#SUMMARY
FR_parameters <- rbind(data.frame(FR_table_par, Window="Full Period"),
                       data.frame(FR_table_par1, Window="1960-1990"),
                       data.frame(FR_table_par2, Window="1970-2000"),
                       data.frame(FR_table_par3, Window="1980-2010"))

FR_metrics_tab <-  rbind(data.frame(rbind(FR_metric1, FR_metric1_arma[-2,]), 
                                    Window="1960-1990"),
                         data.frame(rbind(FR_metric2, FR_metric2_arma[-2,]), 
                                    Window="1970-2000"),
                         data.frame(rbind(FR_metric3, FR_metric3_arma[-2,]), 
                                    Window="1980-2010"))