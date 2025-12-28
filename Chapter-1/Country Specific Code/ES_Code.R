#########
###ESP###
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
load("Data/ESP_data.RData")
ES_apc <- apc; rm(apc, data)
ES_bapc_data <- bapc_data; rm(bapc_data, bapc_deaths, bapc_pop)

# Fitting
ES_fit <- fit_APC_SkN(ES_apc, stanOptionsfit)
ES_effects <- get_effects(ES_fit)

# Parameter
ES_table_par <- table_par(ES_effects, "ESP") #Tab 2

## -- OUT-OF-SAMPLE TEST PHASE -------------

## General settings

baseline1 <- 1960:1990
baseline2 <- 1970:2000
baseline3 <- 1980:2010

h=10 

stanOptions <-  list(nchains = 3, niter = 1000, nwarmup = 500, ncores = 3)

#FIRST WINDOW:

#Our Model
ES_fit_baseline1 <- fit_APC_SkN_forecast(ES_apc, baseline1, stanOptions)

ES_effects_baseline1 <- get_effects(ES_fit_baseline1)

ES_table_par1 <- table_par(ES_effects_baseline1, "ESP") #Tab 2

#ARIMA
ES_SR_forecast1 <- SR_forecasting(ES_apc, ES_effects_baseline1, baseline1, h,
                                  arimaPlot=FALSE, IC=TRUE)

#ARMA
ES_SR_forecast1_arma <- SR_forecasting_arma(ES_apc, 
                                            ES_effects_baseline1,
                                            baseline = baseline1,
                                            h = 10, IC = T, nsim = 1000)

## Poisson Approach
ES_bapc_forecast1 <- bapc_forecasting(ES_bapc_data, baseline1, h)

#Compare Approaches:
ES_compare1 <- ES_apc %>% 
  filter(P %in% max(baseline1+1):(max(baseline1)+h)) %>% 
  mutate("P"=as.factor(P),"A"=as.factor(A)) %>% 
  select(P,A,SR) %>% 
  left_join(ES_bapc_forecast1 %>% 
              mutate(P = as.factor(P), A = as.factor(A)))

ES_metric1 <- metrics(ES_SR_forecast1, ES_bapc_forecast1)
ES_metric1_arma <- metrics_arma(ES_SR_forecast1_arma, ES_bapc_forecast1)

#SECOND WINDOW:

#Our Model
ES_fit_baseline2 <- fit_APC_SkN_forecast(ES_apc, baseline2, stanOptions)

ES_effects_baseline2 <- get_effects(ES_fit_baseline2)

ES_table_par2 <- table_par(ES_effects_baseline2, "ESP") #Tab 2

#ARIMA
ES_SR_forecast2 <- SR_forecasting(ES_apc, ES_effects_baseline2, baseline2, h,
                                  arimaPlot=FALSE, IC=TRUE)

#ARMA
ES_SR_forecast2_arma <- SR_forecasting_arma(ES_apc, 
                                            ES_effects_baseline2,
                                            baseline = baseline2,
                                            h = 10, IC = T, nsim = 1000)

## Poisson Approach
ES_bapc_forecast2 <- bapc_forecasting(ES_bapc_data, baseline2, h)

#Compare Approaches:
ES_compare2 <- ES_apc %>% 
  filter(P %in% max(baseline2+1):(max(baseline2)+h)) %>% 
  mutate("P"=as.factor(P),"A"=as.factor(A)) %>% 
  select(P,A,SR) %>% 
  left_join(ES_bapc_forecast2 %>% 
              mutate(P = as.factor(P), A = as.factor(A)))

ES_metric2 <- metrics(ES_SR_forecast2, ES_bapc_forecast2)
ES_metric2_arma <- metrics_arma(ES_SR_forecast2_arma, ES_bapc_forecast2)

#THIRD WINDOW:

#Our Model
ES_fit_baseline3 <- fit_APC_SkN_forecast(ES_apc, baseline3, stanOptions)

ES_effects_baseline3 <- get_effects(ES_fit_baseline3)

ES_table_par3 <- table_par(ES_effects_baseline3, "ESP") #Tab 2

#ARIMA
ES_SR_forecast3 <- SR_forecasting(ES_apc, ES_effects_baseline3, baseline3, h,
                                  arimaPlot=FALSE, IC=TRUE)

#ARMA
ES_SR_forecast3_arma <- SR_forecasting_arma(ES_apc, 
                                            ES_effects_baseline3,
                                            baseline = baseline3,
                                            h = 10, IC = T, nsim = 1000)

## Poisson Approach
ES_bapc_forecast3 <- bapc_forecasting(ES_bapc_data, baseline3, h)

#Compare Approaches:
ES_compare3 <- ES_apc %>% 
  filter(P %in% max(baseline3+1):(max(baseline3)+h)) %>% 
  mutate("P"=as.factor(P),"A"=as.factor(A)) %>% 
  select(P,A,SR) %>% 
  left_join(ES_bapc_forecast3 %>% 
              mutate(P = as.factor(P), A = as.factor(A)))

ES_metric3 <- metrics(ES_SR_forecast3, ES_bapc_forecast3)
ES_metric3_arma <- metrics_arma(ES_SR_forecast3_arma, ES_bapc_forecast3)


#SUMMARY
ES_parameters <- rbind(data.frame(ES_table_par, Window="Full Period"),
                    data.frame(ES_table_par1, Window="1960-1990"),
                    data.frame(ES_table_par2, Window="1970-2000"),
                    data.frame(ES_table_par3, Window="1980-2010"))

ES_metrics_tab <-  rbind(data.frame(rbind(ES_metric1, ES_metric1_arma[-2,]), 
                                 Window="1960-1990"),
                      data.frame(rbind(ES_metric2, ES_metric2_arma[-2,]), 
                                 Window="1970-2000"),
                      data.frame(rbind(ES_metric3, ES_metric3_arma[-2,]), 
                                 Window="1980-2010"))