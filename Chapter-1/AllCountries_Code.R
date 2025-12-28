# -- NOTE ---------
#This code is prepared to directly obtain the results shown in the paper by just pressing the "source" key. It needs a lot of time to run: about 2 hours for the estimation phase and about 6 hours for the forecasting phase. We recommend using the "Code.R" file for faster and more agile execution. Moreover, this script has few comments, while in Code.R file all the steps are described.

## -- PRELIMINARY STEPS -------------

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

## -- ESTIMATION PHASE -------------

#  First country: UNITED STATES -------------------------------------------

# Data
load("Data/USA_data.RData")

US_apc <- apc; rm(apc, data)
US_bapc_data <- bapc_data; rm(bapc_data, bapc_deaths, bapc_pop)

#Visualization
US_lexis <- US_lexis <- lexis_surface(US_apc, "USA") #Fig 1

# Fitting
US_fit <- fit_APC_SkN(US_apc, stanOptionsfit)
US_effects <- get_effects(US_fit)

#APC effects - visualization
US_plot_effects <- plot_effects(US_effects, "USA") #Fig 2

#Parameters
US_table_par <- table_par(US_effects, "USA") #Tab 2

#  Second country: SPAIN -------------------------------------------

# Data
load("Data/ESP_data.RData")

ES_apc <- apc; rm(apc, data)
ES_bapc_data <- bapc_data; rm(bapc_data, bapc_deaths, bapc_pop)

#Visualization
ES_lexis <- ES_lexis <- lexis_surface(ES_apc, "ESP") #Fig 1

# Fitting
ES_fit <- fit_APC_SkN(ES_apc, stanOptionsfit)
ES_effects <- get_effects(ES_fit)

#APC effects - visualization
ES_plot_effects <- plot_effects(ES_effects, "ESP") #Fig 2

#Parameters
ES_table_par <- table_par(ES_effects, "ESP") #Tab 2

#  Third country: RUSSIA -------------------------------------------


# Data
load("Data/RUS_data.RData")

RU_apc <- apc; rm(apc, data)
RU_bapc_data <- bapc_data; rm(bapc_data, bapc_deaths, bapc_pop)

#Visualization
RU_lexis <- RU_lexis <- lexis_surface(RU_apc, "RUS") #Fig 1

# Fitting
RU_fit <- fit_APC_SkN(RU_apc, stanOptionsfit)
RU_effects <- get_effects(RU_fit)

#APC effects - visualization
RU_plot_effects <- plot_effects(RU_effects, "RUS") #Fig 2

#Parameters
RU_table_par <- table_par(RU_effects, "RUS") #Tab 2


## - VISUALIZATION -----------------------------------------------------------



#FIGURE 1:
ggarrange(US_lexis,
          ES_lexis,
          RU_lexis,
          ncol=1, 
          common.legend = TRUE)

#FIGURE 2:
ggarrange(US_plot_effects,
          ES_plot_effects,
          RU_plot_effects,
          ncol=1, nrow=4, 
          common.legend = T, legend = "bottom")


#TABLE 2 (Only Entire Period):
Tab2_ep <- data.frame("Window"="Full Period",
                      rbind(US_table_par,
                            ES_table_par,
                            RU_table_par))


## -- OUT-OF-SAMPLE TEST PHASE -------------

## General settings

baseline1 <- 1960:1990
baseline2 <- 1970:2000
baseline3 <- 1980:2010

h=10 

stanOptions <-  list(nchains = 3, niter = 1000, nwarmup = 500, ncores = 3)


#  First country: UNITED STATES -------------------------------------------


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

#  Second country: SPAIN -------------------------------------------

#FIRST WINDOW:

#Our Model
ES_fit_baseline1 <- fit_APC_SkN_forecast(ES_apc, baseline1, stanOptions)

ES_effects_baseline1 <- get_effects(ES_fit_baseline1)

ES_table_par1 <- table_par(ES_effects_baseline1, "ES") #Tab 2

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

ES_table_par2 <- table_par(ES_effects_baseline2, "ES") #Tab 2

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

ES_table_par3 <- table_par(ES_effects_baseline3, "ES") #Tab 2

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

#  Third country: RUSSIA -------------------------------------------

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










## -- VISUALIZATION ---------

#TABLE 2
Tab2 <- rbind(data.frame("Country" = "USA", US_parameters),
              data.frame("Country" = "ESP", ES_parameters),
              data.frame("Country" = "RUS", RU_parameters))

#TABLE 2
Tab3 <- rbind(data.frame("Country" = "USA", US_metrics_tab),
              data.frame("Country" = "ESP", ES_metrics_tab),
              data.frame("Country" = "RUS", RU_metrics_tab))


