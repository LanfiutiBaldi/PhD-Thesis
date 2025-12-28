#########
###RUS###
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
load("Data/RUS_data.RData")
RU_apc <- apc; rm(apc, data)
RU_bapc_data <- bapc_data; rm(bapc_data, bapc_deaths, bapc_pop)

# Fitting
RU_fit <- fit_APC_SkN(RU_apc, stanOptionsfit)
RU_effects <- get_effects(RU_fit)

# Parameter
RU_table_par <- table_par(RU_effects, "RUS") #Tab 2

## -- OUT-OF-SAMPLE TRUT PHASE -------------

## General settings

baseline1 <- 1960:1990
baseline2 <- 1970:2000
baseline3 <- 1980:2010

h=10 

stanOptions <-  list(nchains = 3, niter = 1000, nwarmup = 500, ncores = 3)

#FIRST WINDOW:

#Our Model
RU_fit_baseline1 <- fit_APC_SkN_forecast(RU_apc, baseline1, stanOptions)

RU_effects_baseline1 <- get_effects(RU_fit_baseline1)

RU_table_par1 <- table_par(RU_effects_baseline1, "RUS") #Tab 2

#ARIMA
RU_SR_forecast1 <- SR_forecasting(RU_apc, RU_effects_baseline1, baseline1, h,
                                  arimaPlot=FALSE, IC=TRUE)

#ARMA
RU_SR_forecast1_arma <- SR_forecasting_arma(RU_apc, 
                                            RU_effects_baseline1,
                                            baseline = baseline1,
                                            h = 10, IC = T, nsim = 1000)

## Poisson Approach
RU_bapc_forecast1 <- bapc_forecasting(RU_bapc_data, baseline1, h)

#Compare Approaches:
RU_compare1 <- RU_apc %>% 
  filter(P %in% max(baseline1+1):(max(baseline1)+h)) %>% 
  mutate("P"=as.factor(P),"A"=as.factor(A)) %>% 
  select(P,A,SR) %>% 
  left_join(RU_bapc_forecast1 %>% 
              mutate(P = as.factor(P), A = as.factor(A)))

RU_metric1 <- metrics(RU_SR_forecast1, RU_bapc_forecast1)
RU_metric1_arma <- metrics_arma(RU_SR_forecast1_arma, RU_bapc_forecast1)

#SECOND WINDOW:

#Our Model
RU_fit_baseline2 <- fit_APC_SkN_forecast(RU_apc, baseline2, stanOptions)

RU_effects_baseline2 <- get_effects(RU_fit_baseline2)

RU_table_par2 <- table_par(RU_effects_baseline2, "RUS") #Tab 2

#ARIMA
RU_SR_forecast2 <- SR_forecasting(RU_apc, RU_effects_baseline2, baseline2, h,
                                  arimaPlot=FALSE, IC=TRUE)

#ARMA
RU_SR_forecast2_arma <- SR_forecasting_arma(RU_apc, 
                                            RU_effects_baseline2,
                                            baseline = baseline2,
                                            h = 10, IC = T, nsim = 1000)

## Poisson Approach
RU_bapc_forecast2 <- bapc_forecasting(RU_bapc_data, baseline2, h)

#Compare Approaches:
RU_compare2 <- RU_apc %>% 
  filter(P %in% max(baseline2+1):(max(baseline2)+h)) %>% 
  mutate("P"=as.factor(P),"A"=as.factor(A)) %>% 
  select(P,A,SR) %>% 
  left_join(RU_bapc_forecast2 %>% 
              mutate(P = as.factor(P), A = as.factor(A)))

RU_metric2 <- metrics(RU_SR_forecast2, RU_bapc_forecast2)
RU_metric2_arma <- metrics_arma(RU_SR_forecast2_arma, RU_bapc_forecast2)

#SUMMARY
RU_parameters <- rbind(data.frame(RU_table_par, Window="Full Period"),
                       data.frame(RU_table_par1, Window="1960-1990"),
                       data.frame(RU_table_par2, Window="1970-2000"))

RU_metrics_tab <-  rbind(data.frame(rbind(RU_metric1, RU_metric1_arma[-2,]), 
                                    Window="1960-1990"),
                         data.frame(rbind(RU_metric2, RU_metric2_arma[-2,]), 
                                    Window="1970-2000"))