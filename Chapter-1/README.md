# An Age-Period-Cohort model for gender gap in youth and early adult mortality
This repository accompanies the paper "An Age-Period-Cohort model for gender gap in youth and early adult mortality” written by Giacomo Lanfiuti Baldi and Andrea Nigri. 
This repository contains all the material required to replicate the results presented in the main text and in the Supplementary Material.

The repository includes:

1. *Ready-to-use scripts* to directly reproduce the results for the countries analysed in the paper.

2. *Fully generalisable code*, allowing users to replicate the entire workflow for any country in the Human Mortality Database and any baseline length.

3. *Functions*, *Stan scripts*, and *data files* needed to run the APC model, perform forecasting, and reproduce all tables, figures, and metrics.

Detailed descriptions of each component are provided below.

# Ready-to-use code files
This repository includes ready-to-use code files:

**AllCountries_Code.R**
A fully executable script that reproduces exactly the results reported in the paper for the three main countries analysed: United States, Spain, and Russia.

**Country Specific Code (folder)**
Contains standalone scripts for each country included in the main paper and the Supplementary Material.
In addition to USA, Spain, and Russia, this folder provides ready-made code for:
France, Italy, and Great Britain.
These scripts can be run directly to reproduce all tables, figures, and forecasts for the selected country.

# Fully generalisable code
The main code provided in the **Code.R** file is used to perform all the analyses provided in the paper for each country.

To avoid repeating the code six times (United States, Spain, Russia, Italy, France, Great Britain) and for each forecasting interval used in the paper, we provide the general code for replication of the analyses performed for **any country** and **any forecasting period**.

In addition, the code is generalised to work for every country included in the **Human Mortality Database**, and allows obtaining forecasting results based on:
- the chosen baseline period  
- the chosen forecasting length  

**Please note:** some estimates may vary slightly depending on the machine used, without significantly affecting the results.

---

## Code execution process

1. Set the current working directory.  
2. Load all needed functions from the file `Functions.R`.  
3. Download data from the HMD or import the data for the chosen country from the `data` folder.  
4. Manual alterations:
   - **Line 16:** choose any country from HMD (to reproduce our main results use: `USA`, `ESP`, `RUS`).  
   - **Line 17–18:** set your HMD credentials.  
   - **Line 27:** if you do not download data directly from HMD, load it from the `data` folder.  
   - **Line 92:** set the baseline period. Ensure it is included in the APC dataset (recommended options: `1960:1990`, `1970:2000`, `1980:2010`).  

---

## **Code.R**
Contains the main code. It consists of 5 parts plus an initial section loading packages and functions.

1. **READING DATA AND VISUALISATION**  
   Uploading the data:  
   - directly from HMD using credentials, or  
   - from the `data` folder.  
   Reproduces **Figure 1 – Lexis Surface**.

2. **FITTING**  
   Estimates the proposed model for the selected country.  
   Reproduces **Figure 2 – Age-Period-Cohort Effects** and **Table 2 – coefficient estimates**.

3. **IN-SAMPLE ESTIMATION AND OUT-OF-SAMPLE FORECASTING**  
   Define the baseline period and the forecasting horizon `h`.  
   The model is estimated over the baseline without identifiability constraints.  
   Produces sex-ratio forecasts using both **ARIMA** and **ARMA**.

4. **POISSON APPROACH via BAPC PACKAGE**  
   Downloads Deaths and Exposed Population from HMD (skip this step if using data from the folder).  
   Uses BAPC to obtain the Poisson-based sex-ratio forecast.

5. **METRICS**  
   Compares accuracy of:  
   - *Our Model* (ARIMA and ARMA variants)  
   - *BAPC*  
   by computing **RMSE** and **MAE**.  
   Reproduces heat maps reported in the Supplementary Material.

---

## Functions.R

Below is the full list of functions, **exactly as in your original description**, with purpose, inputs, and outputs.

### **lexis_surface(apc, cou)**  
**Purpose:** outputs the Lexis surface of the sex-ratio for the selected country.  
**INPUT:**  
- `apc`: dataset containing sex-ratio (SR) values for each age (A) and period (P)  
- `cou`: name of the country (string), default `""`  
**OUTPUT:**  
- `ggplot` object containing the Lexis surface  

` ` ` ` ` ` ` ` ` `

### **fit_APC_SkN(apc, stanOptions)**  
**Purpose:** estimates the proposed model via Stan.  
**INPUT:**  
- `apc`: APC dataset  
- `stanOptions`: list specifying number of chains, iterations, warm-up, and cores  
**OUTPUT:**  
- dataframe of estimated coefficients and parameters, including  
  `se_mean`, `n_eff`, `Rhat`  

` ` ` ` ` ` ` ` ` `

### **get_effects(fit)**  
**Purpose:** extracts key results from the fit.  
**INPUT:**  
- `fit`: output of `fit_APC_SkN` or `fit_APC_SkN_forecast`  
**OUTPUT:**  
- ordered list containing:  
  - age, period, cohort effects  
  - intercept  
  - regression coefficients for shape and scale  
  - 95% credibility interval bounds  

` ` ` ` ` ` ` ` ` `

### **plot_effects(effects, cou)**  
**Purpose:** graphical display of APC effects.  
**INPUT:**  
- `effects`: output of `get_effects`  
- `cou`: country name (optional)  
**OUTPUT:**  
- `grid.arrange` with age, period, and cohort plots  

` ` ` ` ` ` ` ` ` `

### **table_par(effects, cou)**  
**Purpose:** table of estimated intercept and regression coefficients.  
**INPUT:**  
- `effects`: output of `get_effects`  
- `cou`: optional country name  
**OUTPUT:**  
- dataframe with estimates and 95% CI  

` ` ` ` ` ` ` ` ` `

### **fit_APC_SkN_forecast(apc, baseline, stanOptions)**  
**Purpose:** same as `fit_APC_SkN` but applied only to baseline years and without identifiability constraints.  
**INPUT:**  
- `apc`: dataset  
- `baseline`: vector of years included in estimation  
- `stanOptions`: list of stan settings  
**OUTPUT:**  
- dataframe of estimated parameters  

` ` ` ` ` ` ` ` ` `

### **SR_forecasting(apc, effects, baseline, h, arimaPlot)**  
**Purpose:** forecasts sex-ratio `h` years ahead via ARIMA.  
**INPUT:**  
- `apc`: dataset  
- `effects`: output of `get_effects`  
- `baseline`: baseline vector  
- `h`: forecasting horizon  
- `arimaPlot`: whether to show ARIMA projections (default: `FALSE`)  
**OUTPUT:**  
- dataframe with forecasted sex-ratio values  

` ` ` ` ` ` ` ` ` `

### **SR_forecasting_ARMA(apc, effects, baseline, h, arimaPlot)**  
**Purpose:** like `SR_forecasting` but using ARMA projections.  

` ` ` ` ` ` ` ` ` `

### **bapc_forecasting(bapc_data, baseline, h)**  
**Purpose:** forecasting via the Poisson BAPC approach.  
**INPUT:**  
- `bapc_data`: dataframe with Deaths and Population (Pop) for Age (A), Period (P), Sex  
- `baseline`: baseline years  
- `h`: forecasting horizon  
**OUTPUT:**  
- dataframe with forecasted sex-ratio values  

` ` ` ` ` ` ` ` ` `

### **metrics(skn_forecast, bapc_forecast)**  
**Purpose:** assess which model performs best.  
**INPUT:**  
- `skn_forecast`: output of SR forecasting  
- `bapc_forecast`: output of BAPC forecasting  
**OUTPUT:**  
- dataframe with **RMSE** and **MAE** values  

` ` ` ` ` ` ` ` ` `

### **heat_map(skn_forecast, bapc_forecast)**  
**Purpose:** produce heat maps of absolute differences.  
**OUTPUT:**  
- two `ggplot` heat maps  

---

## Stan Scripts

- **fit_apc_lss_generalized.stan**: Stan script for estimating the proposed model  
- **fit_apc_lss_FORECAST.stan**: Stan script for estimating the model on baseline without identifiability constraints  

---

## Data Folder

Contains all data taken from the Human Mortality Database required for the analysis in the 4 countries described in the paper.  
Data are organised as `Country_data.Rdata`.

---

## Figures Folder

Will contain images saved by `Code.R`.

---

