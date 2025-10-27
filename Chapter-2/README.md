# SkellamProject
Code, result and Data of Skellam Project

This repository contains materials to partially reproduce the results of the working paper:

*An age-period modeling of population gaps in mortality: the cases of cancer and cardiovascular diseases* **()**

## General information
In this repository you can find all the codes, functions and data used to obtain the results of the paper. 
The main code allows the replication of modelling in each country, period and age of your interest. 
By default it performs the analysis for the Italian WHO data in the baseline period 1961-2000 on all ages. 


## Structure
**Main.R**: main code (with some general comment).

**Functions.R**: all functions written for the project.  

**Data folder**: it contains .Rdata files ready to use with data from WHO for Itay, and HMD for US:
- *ready_to_use_Italy.Rdata*
- *ready_to_use_Usa.Rdata*
  
**Other Packages folder**: it contains some packages not anymore available in the CRAN:
- BivPois Package by Michail Tsagris: to run Bivariate Poisson Regession;
- MRW Package by Marius D. Pascariu: to run the forecasting with Multivariate Random Walk.
  
**Output folder**: 
- *Figures.R*: Code to obtain all figures as in the paper;
- *Appendix Figures.R*: Code to obtain all figures as in the appendix;
All the file .Rdata with the result of the analysis run in the paper:
- *All-Ages 1961-2000.Rdata*;
- *All-Ages 1971-2000.Rdata*;
- *All-Ages 1981-2000.Rdata*;
- *Adult-Ages 1961-2000.Rdata*;
- *Adult-Ages 1971-2000.Rdata*;
- *Adult-Ages 1981-2000.Rdata*;

**Other Result folder**: it contains useful HMTL files to visualize several results
(if you want to see it as HTML, you need to download the files. Alternative, open the PDF version of them):
- Age-Period Effects.HTML
- All-Ages Results.HTML
- Adult-Ages Results.HTML


