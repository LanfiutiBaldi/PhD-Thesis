#Code to obtain all figures as in the paper

library(tidyverse)

#FIG 1
load("Output/All-Ages 1961-2000.RData")
plot_initial(dt)

#FIG 2
load("Output/All-Ages 1961-2000.RData")

accuracy_heatmap(check_for, "RMSE")

#FIG 3
load("Output/All-Ages 1961-2000.RData")

plot_agesfitting(check_tot)

#FIG 4
load("Output/Adult-Ages 1961-2000.RData")

accuracy_heatmap(check_for, "RMSE")

#FIG 5
load("Output/Adult-Ages 1961-2000.RData")

plot_agesfitting(check_tot)



