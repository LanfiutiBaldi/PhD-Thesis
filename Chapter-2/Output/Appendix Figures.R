#Code to obtain all figures as in the appendix

library(tidyverse)

load("Output/All-Ages 1971-2000.RData")

accuracy_heatmap(check_for, "RMSE")
plot_agesfitting(check_tot%>% 
                   filter(Year>1970))

load("Output/All-Ages 1981-2000.RData")

accuracy_heatmap(check_for, "RMSE")
plot_agesfitting(check_tot%>% 
                   filter(Year>1980))

load("Output/Adult-Ages 1971-2000.RData")

accuracy_heatmap(check_for, "RMSE")
plot_agesfitting(check_tot %>% 
                   filter(Year>1970))

load("Output/Adult-Ages 1981-2000.RData")

accuracy_heatmap(check_for, "RMSE")
plot_agesfitting(check_tot%>% 
                   filter(Year>1980))

