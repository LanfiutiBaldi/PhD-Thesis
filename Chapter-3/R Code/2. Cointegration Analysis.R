rm(list=ls(all=TRUE))
library(egcm)
library(latticeExtra)
library(tseries)
library(urca)
library(vars)
library(reshape)
library(FitAR)
library(data.table)
library(tidyverse)
library(knitr)
library(ggpubr)

load("C:/Users/990840/Desktop/DOC/Progetti/Deaths of despair/Post-Caronte/PGLE_Nuts1.RData")
load("C:/Users/990840/Desktop/DOC/Progetti/Deaths of despair/Post-Caronte/Italia_PGLE.RData")

PGLE_Italia <- Italia_PGLE %>% 
  mutate("Nuts1"="Italia", .before=Year) %>% 
  mutate(Sex= case_when(
    Sex == "Male" ~ 1,
    Sex == "Female" ~ 2,))

PGLE_Nuts1 <- PGLE_Nuts1 %>% 
  mutate(Nuts1 = as.character(Nuts1))

PGLE <- rbind(PGLE_Italia, PGLE_Nuts1) %>%
  pivot_wider(names_from = RemovedCause, values_from = PGLE) %>%
  ungroup() 

p.value <- 0.05

output <- NULL

for (causename in c("Drugs", "Alcohol", "Suicide", "Dod")) {
  for(s in 1:2) {
    
    training1 <- PGLE %>% 
      filter(Sex==s) %>% 
      select(Nuts1, Year, "Cause"=all_of(causename)) %>% 
      pivot_wider(names_from = Nuts1, values_from = Cause) %>% 
      column_to_rownames("Year") %>% 
      as.matrix()
    
    diff_t <- diff(training1)
    
    #UNIT ROOT TEST
    URtest1 <- matrix(0, nrow=ncol(training1), ncol=3)
    colnames(URtest1) <- c("PP", "ADF", "KPSS")
    rownames(URtest1) <- colnames(training1)
    
    URdiff1 <- URtest1
    NonStationary <- c()
    I1 <- c()
    
    for (i in 1:ncol(training1)){
      URtest1[i,1] <- pp.test(training1[,i])$p.value
      URtest1[i,2] <- adf.test(training1[,i])$p.value
      #URtest1[i,3] <- kpss.test(training1[,i], null="Trend")$p.value
      URtest1[i,3] <- kpss.test(training1[,i], null="Level")$p.value
      
      NonStationary[i] <- ifelse(sum(URtest1[i,1] > p.value,
                                     URtest1[i,2] > p.value,
                                     URtest1[i,3] < p.value)>=2,
                                 TRUE, FALSE)
      
      URdiff1[i,1] <- tseries::pp.test(diff_t[,i])$p.value
      URdiff1[i,2] <- tseries::adf.test(diff_t[,i])$p.value
      URdiff1[i,3] <- tseries::kpss.test(diff_t[,i], null="Level")$p.value
      
      I1[i] <- ifelse(sum(URdiff1[i,1] < p.value,
                          URdiff1[i,2] < p.value,
                          URdiff1[i,3] > p.value)>=2,
                      TRUE, FALSE)
    }
    
    output <- rbind(output, 
                    data.frame("Cause"=causename, 
                               "Sex"=ifelse(s==1, "Male", "Female"),
                               "Nuts1"=rownames(URdiff1),
                               round(URtest1,4),
                               "NonStationary"=NonStationary,
                               "I1"=I1))
    rownames(output) <- 1:nrow(output)
    
  }
}

round(sum(output$NonStationary)/nrow(output),3)*100
round(sum(output$I1)/nrow(output),3)*100

#ANALISI BETWEEN NUTS1

Coint.Rel_between <- NULL
plot_between <- list()

i<-0

for (causename in c("Drugs", "Alcohol", "Suicide", "Dod")) {
  for(s in 1:2) {
    i<- i+1
    
    print(paste(causename, ifelse(s==1, "Male", "Female")))
    
    cs_plot <- PGLE %>% 
      filter(Nuts1!= "Italia") %>% 
      filter(Sex==s) %>% 
      select(Nuts1, Year, "Cause"=all_of(causename)) %>% 
      ggplot(aes(x=Year, y=Cause, col=Nuts1, group=Nuts1))+
      geom_line(linetype=2)+
      geom_smooth(method="lm", alpha=0.1, size=0.3)+
      theme_bw()+
      ylab(" ")+xlab(" ")+
      ggtitle(paste(causename, ifelse(s==1, "Male", "Female")))+
      theme(plot.title = element_text(size=10))
    #scale_color_manual(values=c("#00BA38",  "#F8766D", "#619CFF", "darkviolet"))
    
    plot_between[[i]] <- cs_plot
    
    training1 <- PGLE %>% 
      filter(Nuts1!= "Italia") %>% 
      filter(Sex==s) %>% 
      select(Nuts1, Year, "Cause"=all_of(causename)) %>% 
      pivot_wider(names_from = Nuts1, values_from = Cause) %>% 
      column_to_rownames("Year") %>% 
      as.matrix()
    
    
    H1.trace  <- urca::ca.jo(data.frame(training1), type="trace", spec="longrun", dumvar = NULL)
    H1.trace.c<- urca::ca.jo(data.frame(training1), type="trace", spec="longrun", dumvar = NULL,
                             ecdet = "const")
    
    H1.eigen  <- urca::ca.jo(data.frame(training1), type="eigen", spec="longrun", dumvar = NULL)
    H1.eigen.c<- urca::ca.jo(data.frame(training1), type="eigen", spec="longrun", dumvar = NULL,
                             ecdet = "const")
    
    trace.r <- data.frame(H1.trace@cval, H1.trace@teststat) %>%
      mutate(r = as.character(rownames(.))) %>% 
      mutate(r = gsub("\\|","", r)) %>% 
      select(r, X5pct, "Stat"=H1.trace.teststat) %>% 
      mutate(n = as.numeric(gsub("r <= ", "", gsub("r = ", "r <= ", r)))) %>% 
      mutate(NumberRelation = ifelse(X5pct < Stat, TRUE, FALSE))
    
    eigen.r <- data.frame(H1.eigen@cval, H1.eigen@teststat) %>%
      mutate(r = as.character(rownames(.))) %>% 
      mutate(r = gsub("\\|","", r)) %>% 
      select(r, X5pct, "Stat"=H1.eigen.teststat)%>% 
      mutate(n = as.numeric(gsub("r <= ", "", gsub("r = ", "r <= ", r)))) %>% 
      mutate(NumberRelation = ifelse(X5pct < Stat, TRUE, FALSE))
    
    #data.frame(H1.eigen.c@cval, H1.eigen.c@teststat)
    
    r = sort(intersect(trace.r %>% filter(NumberRelation) %>% pull(n),
                       eigen.r %>% filter(NumberRelation) %>% pull(n)),
             decreasing = T)
    
    r <- r[r!=0]
    
    attempts <- NULL
    
    for(rr in r) {
      t_adf <- adf.test(c(t(urca::cajorls(H1.trace,  r=rr)$beta[,1])%*%t(H1.trace@ZK)))$p.value
      t_pp <- pp.test(c(t(urca::cajorls(H1.trace,   r=rr)$beta[,1])%*%t(H1.trace@ZK)))$p.value
      t_kpss <- kpss.test(c(t(urca::cajorls(H1.trace, r=rr)$beta[,1])%*%t(H1.trace@ZK)))$p.value
      
      #Controllare se questa relazion è giusta
      RelStationary <- ifelse(sum(t_adf < p.value,
                                  t_pp < p.value,
                                  t_kpss > p.value)>=2, TRUE, FALSE)
      
      attempts <- as.data.frame(rbind(attempts,
                                      cbind(rr, RelStationary)))
      
      if(RelStationary){
        break
      }
      
      if(!RelStationary & rr==last(r)){
        attempts <- as.data.frame(rbind(attempts,
                                        data.frame("rr"=0, "RelStationary"=99)))
      }
      
    }
    
    if(is_empty(r)){
      attempts <- data.frame("rr"=0, "RelStationary"=99)
    }
    
    colnames(attempts)=c("r", "RelStationary")
    
    Coint.Rel_between <- rbind(Coint.Rel_between, 
                               data.frame("Cause"=causename, 
                                          "Sex"=ifelse(s==1, "Male", "Female"),
                                          NumRel=attempts$r,
                                          "RelationStationary"=attempts$RelStationary)) 
  }
}

Coint.Rel_between <- Coint.Rel_between %>% 
  mutate(RelationStationary = case_when(
    RelationStationary==1 ~ "Stationary", 
    RelationStationary==0 ~ "Non Stationary",
    RelationStationary==99 ~ "No relations",))

#plot_between <- lapply(plot_between, ggplotGrob)
ggarrange(plot_between[[1]], plot_between[[3]],
          plot_between[[5]], plot_between[[7]],
          plot_between[[2]], plot_between[[6]],
          plot_between[[4]], plot_between[[8]],
          nrow=2, ncol=4, common.legend = TRUE)

#ANALISI WITHIN NUTS1
Coint.Rel_within <- NULL
plot_within <- list()

i<-0

for (nuts1 in c("Nord-Ovest", "Nord-Est", "Centro", "Sud", "Isole", "Italia")) {
  for(s in 1:2) {
    i <- i+1
    
    print(paste(nuts1, ifelse(s==1, "Male", "Female")))
    
    ns_plot <- PGLE %>% 
      select(-Dod) %>% 
      pivot_longer(cols=c(4:6), names_to = "Cause", values_to = "PGLE") %>% 
      pivot_wider(names_from = "Nuts1", values_from = "PGLE") %>% 
      filter(Sex==s) %>% 
      select(Cause, Year, "Nuts1"=all_of(nuts1)) %>% 
      ggplot(aes(x=Year, y=Nuts1, col=Cause, group=Cause))+
      geom_line(linetype=2)+
      geom_smooth(method="lm", alpha=0.1, size=0.3)+
      theme_bw()+
      ylab("")+xlab("")+
      ggtitle(paste(nuts1, ifelse(s==1, "Male", "Female")))+
      theme(plot.title = element_text(size=10))
    
    plot_within[[i]] <- ns_plot
    
    
    training1 <- PGLE %>% 
      select(-Dod) %>% 
      pivot_longer(cols=c(4:6), names_to = "Cause", values_to = "PGLE") %>% 
      pivot_wider(names_from = "Nuts1", values_from = "PGLE") %>% 
      filter(Sex==s) %>% 
      select(Cause, Year, "Nuts1"=all_of(nuts1)) %>% 
      pivot_wider(names_from = Cause, values_from = Nuts1) %>% 
      column_to_rownames("Year") %>% 
      as.matrix()
    
    
    H1.trace  <- urca::ca.jo(data.frame(training1), type="trace", spec="longrun", dumvar = NULL)
    H1.trace.c<- urca::ca.jo(data.frame(training1), type="trace", spec="longrun", dumvar = NULL,
                             ecdet = "const")
    
    H1.eigen  <- urca::ca.jo(data.frame(training1), type="eigen", spec="longrun", dumvar = NULL)
    H1.eigen.c<- urca::ca.jo(data.frame(training1), type="eigen", spec="longrun", dumvar = NULL,
                             ecdet = "const")
    
    trace.r <- data.frame(H1.trace@cval, H1.trace@teststat) %>%
      mutate(r = as.character(rownames(.))) %>% 
      mutate(r = gsub("\\|","", r)) %>% 
      select(r, X5pct, "Stat"=H1.trace.teststat) %>% 
      mutate(n = as.numeric(gsub("r <= ", "", gsub("r = ", "r <= ", r)))) %>% 
      mutate(NumberRelation = ifelse(X5pct < Stat, TRUE, FALSE))
    
    eigen.r <- data.frame(H1.eigen@cval, H1.eigen@teststat) %>%
      mutate(r = as.character(rownames(.))) %>% 
      mutate(r = gsub("\\|","", r)) %>% 
      select(r, X5pct, "Stat"=H1.eigen.teststat)%>% 
      mutate(n = as.numeric(gsub("r <= ", "", gsub("r = ", "r <= ", r)))) %>% 
      mutate(NumberRelation = ifelse(X5pct < Stat, TRUE, FALSE))
    
    #data.frame(H1.eigen.c@cval, H1.eigen.c@teststat)
    
    r = sort(intersect(trace.r %>% filter(NumberRelation) %>% pull(n),
                       eigen.r %>% filter(NumberRelation) %>% pull(n)),
             decreasing = T)
    
    r <- r[r!=0]
    
    attempts <- NULL
    
    for(rr in r) {
      t_adf <- adf.test(c(t(urca::cajorls(H1.trace,  r=rr)$beta[,1])%*%t(H1.trace@ZK)))$p.value
      t_pp <- pp.test(c(t(urca::cajorls(H1.trace,   r=rr)$beta[,1])%*%t(H1.trace@ZK)))$p.value
      t_kpss <- kpss.test(c(t(urca::cajorls(H1.trace, r=rr)$beta[,1])%*%t(H1.trace@ZK)))$p.value
      
      #Controllare se questa relazion è giusta
      RelStationary <- ifelse(sum(t_adf < p.value,
                                  t_pp < p.value,
                                  t_kpss > p.value)>=2, TRUE, FALSE)
      
      attempts <- as.data.frame(rbind(attempts,
                                      cbind(rr, RelStationary)))
      
      if(RelStationary){
        break
      }
      
      if(!RelStationary & rr==last(r)){
        attempts <- as.data.frame(rbind(attempts,
                                        data.frame("rr"=0, "RelStationary"=99)))
      }
      
    }
    
    if(is.null(attempts)){
      attempts <- data.frame("rr"=0, "RelStationary"=99)
    }
    colnames(attempts)=c("r", "RelStationary")
    
    Coint.Rel_within <- rbind(Coint.Rel_within, 
                              data.frame("Nuts1"=nuts1, 
                                         "Sex"=ifelse(s==1, "Male", "Female"),
                                         NumRel=attempts$r,
                                         "RelationStationary"=attempts$RelStationary)) 
    
  }
}

Coint.Rel_within <- Coint.Rel_within %>% 
  mutate(RelationStationary = case_when(
    RelationStationary==1 ~ "Stationary", 
    RelationStationary==0 ~ "Non Stationary",
    RelationStationary==99 ~ "No relations",))

ggarrange(plot_within[[1]], plot_within[[3]],
          plot_within[[5]], plot_within[[7]],
          plot_within[[9]], plot_within[[11]],
          plot_within[[2]], plot_within[[4]],
          plot_within[[6]], plot_within[[8]],
          plot_within[[10]], plot_within[[12]],
          nrow=2, ncol=6, common.legend = TRUE)

#Cleaning Global Enviroment
rm(list = setdiff(ls(), c("PGLE", 
                          "Coint.Rel_within", "Coint.Rel_between",
                          "output",
                          "plot_between", "plot_within")))

Coint.Rel_between$Type = "Between"
Coint.Rel_within$Type = "Within"

colnames(Coint.Rel_between)=c("Variable", "Sex", "NumRel", "RelationStationary", "Type")
colnames(Coint.Rel_within)=c("Variable", "Sex", "NumRel", "RelationStationary", "Type")
Coint.rel <- rbind(Coint.Rel_between, Coint.Rel_within)

