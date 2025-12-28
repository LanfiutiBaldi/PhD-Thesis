## ------------------------------------------------------------------- ##
##  R code to estimate and forecast sex ratio of age-specific detah rates of male over female using the an Age-Period-Cohort model. 

# load useful packages
library(HMDHFDplus)
library(tidyverse)
library(viridis)
library(rstan)

# load functions useful to plot, estimate and forecasting the Sex-Ratio
source("Functions.R")

## -- READING DATA AND VISUALISATION -------------

#Read mortality data   
cou <- "USA"   # choose any country from HMD (to reproduce our results please use: USA, ESP, RUS)
username <- ""   # set your HMD credentials
password <- ""   # set your HMD credentials

# Load the data for the selected country from the Data folder. 
data <- NULL
load("Data/USA_data.RData")

# Or download the mortality rates for the two sexes in the chosen country: 
if(is.null(data)){
  data <- readHMDweb(CNTRY=cou, item="Mx_1x1",
                  username=username,
                  password=password,
                  fixup=T) 
}

# Transforming data into a apc format:
apc <- data %>% 
  select(Year, Age, Female, Male) %>% 
  mutate(Cohort=Year-Age) %>% 
  filter(Age<101 & (Year>1959 & Year<2021)) %>%
  mutate(ratio = log(Male/Female)) %>% 
  relocate(Cohort, .after = Age) %>% 
  filter(Age<46) %>%
  select("SR"=ratio,"P"=Year,"A"=Age, "K"=Cohort)

# The apc dataset contains: 
# SR: is the collection of the Sex-Ratio (logarithm of the ratio of the age-specific mortality rates of the two sexes).
# P: contains the years between 1960-2020, repeated as many time as the number of ages.
# A: contains the ages between 0:45, repeated as many time as the number of years.
# K: contains the cohorts between 1915:2020, obtained as P-A.

## -- LEXIS SURFACE 
# Figure 1: Lexis Surface of the Sex Ratio of the Age-Specific Mortality Rates

lexis_surface(apc, cou)

#To save it in Figures folder.
ggsave(paste("Figures/Figure1_",cou,".png", sep=""), 
       plot = lexis_surface(apc, cou), 
       device = "png")

## -- FITTING -------------

# Provide options for Stan estimation of the model. The options provided here are those have been used in the paper:

stanOptions <- list(nchains = 4,
                    niter = 2000,
                    nwarmup = 1000,
                    ncores = 4)

# The function fit_APC_SkN takes in input the APC dataset and the desired options for Stan. It gives as output a dataframe containg all the estimated coefficients and parameters with relative Monte Carlo standard error (se_mean), the effective sample size (n_eff), and the R-hat statistic (Rhat). 

# WARNING: it might take a while (~ 2.30 mins to read the model + ~30 mins for the running. It depends on the machine you are using.)
fit <- fit_APC_SkN(apc, stanOptions)

#The get_effects function outputs a suitably ordered list with the estimated effects of age, period, cohort, the intercept and the values of the regression coefficients associated with shape and scale parameters and all the relative values of the credibility interval bounds.
effects <- get_effects(fit)

## -- PLOT AGE-PEIROD-COHORT EFFECTS 

# Figure 2: Age-Period-Cohort Effects
plot_effects(effects, cou)

#To save it in Figures folder.
ggsave(paste("Figures/Figure2_",cou,".png", sep=""), 
       plot = plot_effects(effects, cou), 
       device = "png")

## -- TABLE OF INTERCEPT AND SHAPE AND SCALE REGRESSION COEFFICIENTS:

#Table 2: Intercept and regression coefficients for the Skew-Normal parameters and relative credibility interval at level 95%.
table_par(effects, cou)

## -- IN-SAMPLE ESTIMATION AND OUT-OF-SAMPLE FORECASTING -------

baseline <- 1960:1990 #Set the baseline period. Ensure it is included in apc dataset.
#to reproduce our results please use: 1960:1990, 1970:2000, 1980:2010.

h=10 #Choose for how many years you want forecasting.

stanOptions <- list(nchains = 3,
                    niter = 1000,
                    nwarmup = 500,
                    ncores = 3)

#WARNING: it might take a while (~ 2.30 mins to read the model + ~30 mins for the running.)
fit_baseline <- fit_APC_SkN_forecast(apc, baseline, stanOptions)

effects_baseline <- get_effects(fit_baseline)

table_par(effects_baseline, cou)

#SR_forecasting function provide the forecasting of the Sex-Ratio for h periods ahead. Within it, predictions for period and cohort effects are obtained. Declare arimaPlot=TRUE if you want to show and save (in Figures folder) the plot of the forecasting of the two effects. Figure 2 of the Supplementary Material.

SR_forecast <- SR_forecasting(apc, effects_baseline, baseline, h, cou,
                            arimaPlot=FALSE)

#SR_forecasting_ARMA function provide the forecasting of the Sex-Ratio for h periods ahead, using ARMA model. Within it, predictions for period and cohort effects are obtained. 
SR_forecast_arma <- SR_forecasting_arma(apc, effects_baseline,
                                         baseline, h , IC = T, nsim = 1000)

## -- POISSON APPROACH by BAPC PACKAGE -------------
library(BAPC) 

#To use BAPC we need the series of Deaths and Exposures, not the rates. So we obtain the necessary data again from the Human Mortality Database.

bapc_deaths <- readHMDweb(CNTRY=cou, item="Deaths_1x1",
                          username=username,
                          password=password,
                          fixup=T) %>% 
  mutate(Male=round(Male), Female=round(Female))%>% 
  select(Year, Age, Female, Male) %>% 
  mutate(Cohort=Year-Age) %>% 
  filter(Age<46 & (Year>1959 & Year<2021))

bapc_pop <- readHMDweb(CNTRY=cou, item="Exposures_1x1",
                       username=username,
                       password=password,
                       fixup=T) %>% 
  select(Year, Age, Female, Male) %>% 
  mutate(Cohort=Year-Age) %>% 
  filter(Age<46 & (Year>1959 & Year<2021))

# We combine Deaths and Exposures into a single dataset. This dataset has a similar structure to the apc dataset plus the information about the sex.

bapc_data <- bapc_deaths %>% 
  pivot_longer(cols = c("Male", "Female"), 
               names_to = "Sex", values_to = "Deaths") %>% 
  left_join(bapc_pop %>% 
              pivot_longer(cols = c("Male", "Female"), 
                           names_to = "Sex", values_to = "Pop"))%>%
  mutate("A"=as.factor(Age), "P"=as.factor(Year), "K"=as.factor(Cohort)) %>%
  select("Deaths"=Deaths, "Pop"=Pop, "P"=Year,"A"=Age, "K"=Cohort, Sex)

## -- FORECASTING WITH BAPC -------------#

#Forecasting of the same period as above and based on the same baseline 
#We provide a function to let you to skip the data preparation for the forecasting obtained with BAPC function itself.
bapc_forecast <- bapc_forecasting(bapc_data, baseline, h)

## -- METRICS -------------

#To asses the best model between "Our Model" and "BAPC" we combine the forecasted (with both the models) and the observed values of the Sex-Ratio in a single dataset. 

forecast_window <- max(baseline+1):(max(baseline)+h)

compare <- apc %>% 
  filter(P %in% forecast_window) %>% 
     mutate("P"=as.factor(P),"A"=as.factor(A)) %>% 
     select(P,A,SR) %>% 
     left_join(bapc_forecast %>% 
     mutate(P = as.factor(P), A = as.factor(A)))

#Table 3: RMSE and MAE
metric_arima <- metrics(SR_forecast, bapc_forecast)
metric_arma <- metric(SR_forecast_arma, bapc_forecast)

metric <- data.frame(rbind(metric_arima, metric_arma[-2,]))

#Figure 2 of the Supplementary Material: heat map of the differences.
heat_map(SR_forecast, bapc_forecast) 

#To save it in Figures folder.
ggsave(paste("Figures/SM_Figure2_",cou," in period ",
             min(forecast_window),"-",max(forecast_window),".png", sep=""), 
       plot = heat_map(SR_forecast, bapc_forecast), 
       device = "png")
