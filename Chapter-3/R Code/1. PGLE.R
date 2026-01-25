library(tidyverse)
library(readxl)
library(zoo)
library(ggpubr)
library(DemoTools)
library(scales)

load("dt_fv.Rdata")

source("0. Functions.R")

dt_Italia <- dt %>% 
  select(-PopStd) %>% 
  group_by(Year, Sex, Age, Cause) %>% 
  summarise(Deaths = sum(Deaths),
            Pop = sum(Pop)) %>%
  group_by(Year, Sex, Age) %>% 
  summarise(Deaths = sum(Deaths),
            Pop = mean(Pop)) %>% 
  mutate(Mx = Deaths/Pop)

Italia_mx <-  dt %>% 
  filter(Year < 2019) %>% 
  select(-PopStd) %>% 
  group_by(Year, Sex, Age, Cause) %>% 
  summarise(Deaths = sum(Deaths),
            Pop = sum(Pop)) %>% 
  pivot_wider(names_from = Cause, values_from = Deaths) %>% 
  mutate(All = Other+Drugs+Alcohol+Suicide,
         mxAll = All/Pop,
         mxNoDrugs = (All-Drugs)/Pop,
         mxNoAlcohol = (All-Alcohol)/Pop,
         mxNoSuicide = (All-Suicide)/Pop,
         mxNoDoD = Other/Pop) %>% 
  select(Year, Sex, Age, 
         mxAll, mxNoDrugs, mxNoAlcohol, mxNoSuicide, mxNoDoD) 

Italia_PGLE <- Italia_mx %>% 
  group_by(Year, Sex) %>% 
  summarise(exAll = ex0_demotools(mxAll),
            exNoDrugs = ex0_demotools(mxNoDrugs),
            exNoAlcohol = ex0_demotools(mxNoAlcohol),
            exNoSuicide = ex0_demotools(mxNoSuicide),
            exNoDoD = ex0_demotools(mxNoDoD)) %>%
  mutate(Drugs = exNoDrugs - exAll,
         Alcohol = exNoAlcohol - exAll,
         Suicide = exNoSuicide - exAll,
         Dod= exNoDoD - exAll) %>% 
  select(Year, Sex, Drugs, Alcohol, Suicide, Dod) %>% 
  pivot_longer(cols=c("Drugs", "Alcohol", "Suicide", "Dod"),
               names_to = "RemovedCause", values_to = "PGLE")%>%
  sex_factor() %>% 
  removedcause_factor()

Italia_PGLE %>% 
  mutate("Cause" = RemovedCause) %>% 
  ggplot(aes(x=Year, y=PGLE, col=Cause, group=Cause))+
  #alpha= ifelse(Cause %in% c("Drugs", "Alcohol", "Suicide"), 0.2, 1))) +
  geom_line(size=1)+
  facet_grid("Sex")+
  theme_bw()+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("Drugs" = alpha("#D32F2F", 0.4),   
                                "Alcohol" = alpha("#388E3C",0.4),  
                                "Suicide" = alpha("#1976D2",0.4), 
                                "Dod" = "#7B1FA2")) 


# NUTS 1

dt_Dod <- dt %>% 
  from_code_to_region()

dt_Nuts1 <- dt_Dod %>% 
  make_nuts1() %>% 
  group_by(Year, Nuts1, Sex, Age, Cause) %>% 
  summarise(Deaths = sum(Deaths),
            Pop = sum(Pop)) 

mx_Nuts1 <- dt_Nuts1 %>% 
  filter(Year < 2019) %>% 
  pivot_wider(names_from = Cause, values_from = Deaths) %>% 
  mutate(All = Other+Drugs+Alcohol+Suicide,
         mxAll = All/Pop,
         mxNoDrugs = (All-Drugs)/Pop,
         mxNoAlcohol = (All-Alcohol)/Pop,
         mxNoSuicide = (All-Suicide)/Pop,
         mxNoDoD = Other/Pop) %>% 
  select(Nuts1, Year, Sex, Age, 
         mxAll, mxNoDrugs, mxNoAlcohol, mxNoSuicide, mxNoDoD)

ex0_Nuts1 <- mx_Nuts1 %>% 
  group_by(Nuts1, Year, Sex) %>% 
  summarise(exAll = ex0_demotools(mxAll),
            exNoDrugs = ex0_demotools(mxNoDrugs),
            exNoAlcohol = ex0_demotools(mxNoAlcohol),
            exNoSuicide = ex0_demotools(mxNoSuicide),
            exNoDoD = ex0_demotools(mxNoDoD))

PGLE_Nuts1 <- ex0_Nuts1 %>% 
  mutate(Drugs = exNoDrugs - exAll,
         Alcohol = exNoAlcohol - exAll,
         Suicide = exNoSuicide - exAll,
         Dod= exNoDoD - exAll) %>% 
  select(Nuts1, Year, Sex, Drugs, Alcohol, Suicide, Dod) %>% 
  pivot_longer(cols=c("Drugs", "Alcohol", "Suicide", "Dod"),
               names_to = "RemovedCause", values_to = "PGLE") %>% 
  removedcause_factor()

PGLE_Nuts1 %>% 
  sex_factor() %>% 
  mutate(
    "Nuts1"=fct_rev(Nuts1),
    "Cause"=RemovedCause) %>% 
  ggplot(aes(x=Year, y=PGLE, col=Cause, group=Cause))+
  #geom_smooth(alpha=0.5)+
  geom_line(size=1)+
  facet_grid(c("Sex","Nuts1"))+
  theme_bw()+
  #ggtitle("")+
  scale_color_manual(values = c("Drugs" = alpha("#D32F2F", 0.7),   
                                "Alcohol" = alpha("#388E3C",0.7),  
                                "Suicide" = alpha("#1976D2",0.7), 
                                "Dod" = "#7B1FA2"))+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_x_continuous(breaks = seq(min(Italia_PGLE$Year), 
                                  max(Italia_PGLE$Year), by = 7))


PGLE_Nuts1 %>% 
  sex_factor() %>% 
  ggplot(aes(x=Year, y=Nuts1, z=PGLE))+
  geom_tile(aes(fill=PGLE))+
  facet_grid(c("Sex","RemovedCause"))+
  scale_fill_gradient2(low = "white", mid = "red4", high = "black",
                       midpoint=0.5, 
                       limits = c(0.000, 1))+
  theme_bw()+
  #ggtitle("")+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, hjust = 1))+
  scale_x_continuous(breaks = seq(min(Italia_PGLE$Year), 
                                  max(Italia_PGLE$Year), by = 7))


