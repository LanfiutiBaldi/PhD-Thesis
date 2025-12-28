#Functions associated with the Main code for the paper: 
#An age-period modeling of population gaps in mortality: the cases of cancer and cardiovascular diseases

plot_initial <- function(dt){
  dt %>% 
    pivot_longer(3:5, names_to = "What", values_to = "Deaths") %>% 
    ggplot(aes(x=Age, y = Deaths, col=Year, group=Year))+
    geom_line()+
    labs(color="Year") +
    facet_wrap("What")+
    theme_bw()+
    theme(legend.position = "bottom",
          legend.key.width = unit(0.9, "cm"),
          legend.key.height = unit(0.5, "cm"))
}

grouping_01ages <- function(dt){
  dt <- dt %>% 
    mutate(Age = ifelse(Age<5, 0, Age)) %>% 
    group_by(Year, Age) %>% 
    summarise(Cancer = sum(Cancer),
            Circulatory = sum(Circulatory),
            Diff = sum(Diff))
}

prepare_data <- function(dt){
  SIM <- data.frame("years"=as.factor(dt$Year),
                    "ages"=as_factor(dt$Age), 
                    "Deaths1"=dt$Cancer,
                    "Deaths2"=dt$Circulatory, 
                    "diff"=dt$Diff)
}

Skellam_regression <- function(SIM, init.values=NULL){
  set.seed(4444)

  X <- make_standata(bf(diff ~ years + ages), data=SIM)$X #Design matrix
  
  mpf <- 2^5
  
  y <- SIM$diff
  x <- SIM[, 1:2]
  
  n <- length(y)
  x <- stats::model.matrix( ~., data.frame(x) )
  p <- dim(x)[2]
  
  skelreg <- function(pa) {
    b1 <- pa[1:p]   ;   b2 <- pa[ -c(1:p) ]
    a1 <- x %*% b1      ;     a2 <- x %*% b2
    lam1 <- exp(a1)
    lam2 <- exp(a2)
    a <- 2 * sqrt(lam1 * lam2)
    
    
    sum(lam1 + lam2)  + 0.5 * sum(y*(a1 - a2)) - sum(log(besselI.nuAsym(a, abs(y), k.max=4)))
  }
  
  skelreg_opt <- function(pa) {
    b1 <- pa[1:p]   ;   b2 <- pa[ -c(1:p) ]
    a1 <- x %*% b1      ;     a2 <- x %*% b2
    lam1 <- exp(a1)
    lam2 <- exp(a2)
    a <- 2 * sqrt(lam1 * lam2)
    
    
    sum(lam1+lam2) + 0.5 * sum(y * (a1-a2))-
      sum(as.numeric(log(besselI.nuAsym(mpfr(a, mpf), abs(y), k.max=4))))
  }
  
  options(warn = -1)
  mod <- stats::nlm(skelreg, stats::rnorm(2 * p), 
                    iterlim = 5000) 
  
  mod <- stats::nlm(skelreg, mod$estimate, 
                    iterlim = 5000)
  
  if(is.null(init.values)){init.values <-  mod$estimate}
  
  mod <- stats::optim(init.values, skelreg_opt, 
                      hessian = FALSE, 
                      control = list(maxit = 5000, trace = 3),
                      method = "BFGS")
  
  b1 <- mod$par[1:p]    ;    b2 <- mod$par[ -c(1:p)]
  
  par1_sk <- b1
  par2_sk <- b2
  
  dim(X)
  length(b1)
  
  est1_sk <- exp(X %*% par1_sk)[,1]
  est2_sk <- exp(X %*% par2_sk)[,1]
  
  est_diff_sk <- est1_sk-est2_sk
  
  #AIC e BIC
  logLik_sk <- -mod$value
  k <- length(mod$par) # Number of parameters
  n <- length(y) # Number of observations
  
  AIC_sk <- 2 * k - 2 * logLik_sk
  BIC_sk <- k * log(n) - 2 * logLik_sk
  AICc_sk <- AIC_sk + (2 * k * (k + 1)) / (n - k - 1)
  
  sk_reg <- list(
    "Parameters" = list("Par1"=as.numeric(par1_sk),
                        "Par2"=as.numeric(par2_sk)),
    "Estimates" = data.frame("Pop1"=est1_sk,
                              "Pop2"=est2_sk,
                              "Diff"=est_diff_sk),
    "Metrics" = data.frame("BIC"=BIC_sk,
                           "AIC"=AIC_sk,
                           "AICc"=AICc_sk),
    "logLik"= logLik_sk)
  
  
  return(sk_reg)
  
}

DoublePoisson_regression <- function(SIM){

  X <- make_standata(bf(diff ~ years + ages), data=SIM)$X #Design matrix
  
  fit1_pois <- glm(Deaths1 ~ years + ages, family="poisson", data=SIM)
  fit2_pois <- glm(Deaths2 ~ years + ages, family="poisson", data=SIM)
  
  par1_pois <- fit1_pois$coefficients
  par2_pois <- fit2_pois$coefficients
  
  est1_pois <- exp(X %*% par1_pois)[,1]
  est2_pois <- exp(X %*% par2_pois)[,1]
  
  est_diff_poi <- est1_pois-est2_pois
  
  logLik_fit1 <- as.numeric(logLik(fit1_pois)) 
  logLik_fit2 <- as.numeric(logLik(fit2_pois)) 
  
  logLik_total <- logLik_fit1 + logLik_fit2  
  
  k_fit1 <- length(par1_pois)  
  k_fit2 <- length(par2_pois) 
  k_total <- k_fit1 + k_fit2  
  
  n <- nrow(SIM)  
  
  # Calcolo di AIC e BIC
  AIC_pois <- 2 * k_total - 2 * logLik_total
  BIC_pois <- k_total * log(n) - 2 * logLik_total
  AICc_pois <- AIC_pois + (2 * k_total * (k_total + 1)) / (n - k_total - 1)
  
  dp_reg <- list(
    "Parameters" = list("Par1"=par1_pois,
                        "Par2"=par2_pois),
    "Estimates" = data.frame("Pop1"=est1_pois,
                              "Pop2"=est2_pois,
                              "Diff"=est_diff_poi),
    "Metrics" = data.frame("BIC"=BIC_pois,
                           "AIC"=AIC_pois,
                           "AICc"=AICc_pois),
    "logLik"= logLik_total)

  return(dp_reg)
}

BivariatePoisson_regression <- function(SIM){
  source("Other Packages/Bivariate Poisson Package/bivpois pack.R")
  
  X <- make_standata(bf(diff ~ years + ages), data=SIM)$X #Design matrix
  
  bpreg <- lm.bp(Deaths1 ~ years + ages, 
                 Deaths2 ~ years + ages,
                 ~ years + ages, data=SIM, verbose=FALSE)
  
  par1_bp <- data.frame("coef"=bpreg$coefficients, 
                        "par"=rownames(as.data.frame(bpreg$coefficients))) %>% 
    filter(str_detect(par, "l1"))
  
  par1_bp$par <- str_replace(par1_bp$par, "\\(l1\\)\\:ages", "")
  par1_bp$par <- str_replace(par1_bp$par, "\\(l1\\)\\:years", "0.")
  par1_bp$par <- str_replace(par1_bp$par, "\\(l1\\)\\:\\(Intercept\\)", "0")
  par1_bp$par <- as.numeric(par1_bp$par)
  par1_bp <- par1_bp[order(par1_bp$par), ]$coef
  
  par2_bp <- data.frame("coef"=bpreg$coefficients, 
                        "par"=rownames(as.data.frame(bpreg$coefficients))) %>% 
    filter(str_detect(par, "l2")) 
  
  par3_bp <- data.frame("coef"=bpreg$coefficients, 
                        "par"=rownames(as.data.frame(bpreg$coefficients))) %>% 
    filter(str_detect(par, "l3")) %>% pull(coef)
  
  par2_bp$par <- str_replace(par2_bp$par, "\\(l2\\)\\:ages", "")
  par2_bp$par <- str_replace(par2_bp$par, "\\(l2\\)\\:years", "0.")
  par2_bp$par <- str_replace(par2_bp$par, "\\(l2\\)\\:\\(Intercept\\)", "0")
  par2_bp$par <- as.numeric(par2_bp$par)
  par2_bp <- par2_bp[order(par2_bp$par),]$coef
  
  est_diff_bp <- bpreg$fitted.values[,1]-bpreg$fitted.values[,2]
  
  BIC_bp <- bpreg$BIC[2]
  AIC_bp <- bpreg$AIC[2]
  k <- length(bpreg$coefficients) 
  n <- nrow(X) 
  AICc_bp <- AIC_bp + (2 * k * (k + 1)) / (n - k - 1)

  bp_reg <- list(
    "Parameters" = list("Par1"=par1_bp,
                        "Par2"=par2_bp,
                        "Par3"=par3_bp),
    "Estimates" = data.frame("Pop1"=bpreg$fitted.values[,1],
                              "Pop2"=bpreg$fitted.values[,2],
                              "Diff"=est_diff_bp),
    "Metrics" = data.frame("BIC"=BIC_bp,
                           "AIC"=AIC_bp,
                           "AICc"=AICc_bp),
    "logLik"= max(bpreg$loglikelihood))
  return(bp_reg)
  
  
  
}

get_effects <- function(reg_output, SIM){

  nages <- nlevels(SIM$ages)
  nyears <- nlevels(SIM$years)
  ages <- unique(SIM$ages)
  years <- unique(SIM$years)
  
  #Effetti 
  par1 <- reg_output$Parameters$Par1
  par2 <- reg_output$Parameters$Par2
  
  intercept1 <- par1[1]
  year_eff1 <- c(0,par1[2:(nyears)])
  age_eff1 <- c(0,par1[(nyears+1):(nages+nyears-1)])
  
  intercept2 <- par2[1]
  year_eff2 <- c(0,par2[2:(nyears)])
  age_eff2 <- c(0,par2[(nyears+1):(nages+nyears-1)])
  
  eff <- list(
    "intercepts" = list(
        "intercept1"=intercept1,
        "intercept2"=intercept2),
    "age_effects" = list(
        "age_effect1" = age_eff1,
        "age_effect2" = age_eff2),
    "year_effects" = list(
        "year_effect1" = year_eff1,
        "year_effect2" = year_eff2))

  return(eff)
}

join_eff <- function(sk_eff, dp_eff, bp_eff, SIM){
  
  nages <- nlevels(SIM$ages)
  nyears <- nlevels(SIM$years)
  ages <- unique(SIM$ages)
  years <- unique(SIM$years)
  
  intercepts1 <- data.frame("Skellam"= sk_eff$intercepts$intercept1,
                            "Double Poisson"=dp_eff$intercepts$intercept1,
                            "Bivariate Poisson"=bp_eff$intercepts$intercept1,
                            "Cause"="Cancer")
  
  intercepts2 <- data.frame("Skellam"= sk_eff$intercepts$intercept2,
                            "Double Poisson"=dp_eff$intercepts$intercept2,
                            "Bivariate Poisson"=bp_eff$intercepts$intercept2,
                            "Cause"="Circ")
  
  year_eff1 <- data_frame("Skellam"=sk_eff$year_effects$year_effect1,
                          "Double Poisson"=dp_eff$year_effects$year_effect1, 
                        "Bivariate Poisson"=bp_eff$year_effects$year_effect1,
                          "years"=years) %>% 
    pivot_longer(cols=1:3, names_to = "Model") %>% 
    mutate("eff"="Year", "Cause"="Cancer")
  
  year_eff2 <- data_frame("Skellam"=sk_eff$year_effects$year_effect2,
                          "Double Poisson"=dp_eff$year_effects$year_effect2, 
                        "Bivariate Poisson"=bp_eff$year_effects$year_effect2,
                          "years"=years) %>% 
    pivot_longer(cols=1:3, names_to = "Model") %>% 
    mutate("eff"="Year", "Cause"="Circ")

  age_eff1 <- data_frame("Skellam"=sk_eff$age_effects$age_effect1,
                          "Double Poisson"=dp_eff$age_effects$age_effect1, 
                          "Bivariate Poisson"=bp_eff$age_effects$age_effect1,
                          "ages"=ages) %>% 
    pivot_longer(cols=1:3, names_to = "Model") %>% 
    mutate("eff"="Age", "Cause"="Cancer")
    #mutate("eff"="Age Effect - Cancer")
  
  age_eff2 <- data_frame("Skellam"=sk_eff$age_effects$age_effect2,
                          "Double Poisson"=dp_eff$age_effects$age_effect2, 
                          "Bivariate Poisson"=bp_eff$age_effects$age_effect2,
                          "ages"=ages) %>% 
    pivot_longer(cols=1:3, names_to = "Model") %>% 
    mutate("eff"="Age", "Cause"="Circ")
    #mutate("eff"="Age Effect - Circ")
  
  intercepts <- as.data.frame(rbind(intercepts1, intercepts2))
  year_eff <- as.data.frame(rbind(year_eff1, year_eff2))
  age_eff <- as.data.frame(rbind(age_eff1, age_eff2))
  
  All_eff <- list(
    "intercepts" = intercepts,
    "age_eff" = age_eff,
    "year_eff" = year_eff)
  
  return(All_eff)
}

plot_eff <- function(All_eff){
  age_eff <- All_eff$age_eff
  year_eff <- All_eff$year_eff
  
  age_eff_plot <- age_eff %>% 
    mutate(eff = paste(eff, Cause, sep=" Effects - ")) %>% 
    ggplot(aes(x=ages, y=value, col=Model, group=Model))+
    geom_line()+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    geom_hline(yintercept=0, lty=3, alpha=1)+
    facet_wrap(vars(eff))
  
  year_eff_plot <- year_eff %>%
    mutate(eff = paste(eff, Cause, sep=" Effects - ")) %>% 
    ggplot(aes(x=years, y=value, col=Model, group=Model))+
    geom_line()+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    geom_hline(yintercept=0, lty=3, alpha=1)+
    facet_wrap(vars(eff))
  
  print(ggarrange(age_eff_plot, year_eff_plot, 
            common.legend = TRUE, legend = "bottom",
            ncol = 1))
  
}

metrics <- function(check){
  rmse_sk <- rmse(check$ObservedGap, check$Skellam)
  rmse_dp <- rmse(check$ObservedGap, check$DoublePoisson)
  rmse_bp <- rmse(check$ObservedGap, check$BivariatePoisson)
  
  mae_sk <- mae(check$ObservedGap, check$Skellam)
  mae_dp <- mae(check$ObservedGap, check$DoublePoisson)
  mae_bp <- mae(check$ObservedGap, check$BivariatePoisson)
  
  mape_sk <- mape(check$ObservedGap, check$Skellam)
  mape_dp <- mape(check$ObservedGap, check$DoublePoisson)
  mape_bp <- mape(check$ObservedGap, check$BivariatePoisson)
  
  ev <- data.frame("RMSE"=c(rmse_sk, rmse_dp, rmse_bp), 
                   "MAE"= c(mae_sk, mae_dp, mae_bp),
                   "MAPE"=c(mape_sk ,mape_dp, mape_bp),
                   row.names = c("Skellam", 
                                 "Double Poisson",
                                 "Bivariate Poisson"))
  
  return(ev)
}

forecasting_yearEffects <- function(All_eff, h){
  
  year_eff <- All_eff$year_eff %>% 
    mutate(years = as.factor(as.numeric(as.character(years))))
  
  source("Other Packages/MRW Package/mrw.R")
  
  years <- as.numeric(levels(year_eff$years))
  
  forecast <- NULL
  
  for (model in c("Skellam", "Double Poisson", "Bivariate Poisson")) {
    
    training1 <- year_eff %>% 
      filter(Model == model) %>%
      select(years, value, Cause) %>% 
      pivot_wider(names_from = Cause, values_from = value) %>% 
      column_to_rownames("years") %>% 
      as.matrix() %>% 
      t()
    
    M <- model.MRW(data = training1, x = c("Cancer", "Circ"), y = , 
                   include.drift = TRUE)
    P <- predict(M, h = h)
    
    forec <- as.matrix(P$predicted.values) %>%
      t() %>% 
      as.data.frame() %>% 
      mutate(years = as.numeric(rownames(.)), .before="Cancer") %>% 
      mutate(Model = model, .after="years")
    
    forecast <- rbind(forecast, forec)
  }
  
  return(forecast)
  
}

forecasting_differences <- function(SIM_tot, All_eff_base, year_eff_forec, 
                                    h, nages){
  
  X_tot <- make_standata(bf(diff ~ years + ages), data=SIM_tot)$X 
  
  #Skellam
  par1_sk_for <- c(All_eff_base$intercepts$Skellam[1],
                    c((All_eff_base$year_eff %>% 
                      filter(Model=="Skellam" & 
                               Cause == "Cancer") %>% 
                      pull(value))[-1], 
                      year_eff_forec %>% 
                        filter(Model=="Skellam") %>% 
                        pull(Cancer)),
                    (All_eff_base$age_eff %>% 
                      filter(Model=="Skellam" & 
                               Cause == "Cancer") %>% 
                      pull(value))[-1])
  
  par2_sk_for <- c( All_eff_base$intercepts$Skellam[2],
                c((All_eff_base$year_eff %>% 
                    filter(Model=="Skellam" & 
                             Cause == "Circ") %>% 
                    pull(value))[-1], year_eff_forec %>% 
                    filter(Model=="Skellam") %>% 
                    pull(Circ)),
                (All_eff_base$age_eff %>% 
                  filter(Model=="Skellam" & 
                           Cause == "Circ") %>% 
                  pull(value))[-1])
  
  Skellam_estfor <- exp(X_tot %*% par1_sk_for) -  exp(X_tot %*% par2_sk_for)
  
  Skellam_for <- tail(Skellam_estfor, h*nages)
  
  #Double Poisson
  
  par1_dp_for <- c(All_eff_base$intercepts$Double.Poisson[1],
                      c((All_eff_base$year_eff %>% 
                      filter(Model=="Double Poisson" & 
                               Cause=="Cancer") %>% 
                          pull(value))[-1], 
                      year_eff_forec %>% 
                        filter(Model=="Double Poisson") %>% 
                        pull(Cancer)),
                      (All_eff_base$age_eff %>% 
                        filter(Model=="Double Poisson" & 
                                 Cause=="Cancer") %>%
                        pull(value))[-1])
  
  par2_dp_for <- c(All_eff_base$intercepts$Double.Poisson[2],
                   c((All_eff_base$year_eff %>% 
                       filter(Model=="Double Poisson" & 
                                Cause == "Circ") %>% 
                       pull(value))[-1], year_eff_forec %>% 
                       filter(Model=="Double Poisson") %>% 
                       pull(Circ)),
                   (All_eff_base$age_eff %>% 
                     filter(Model=="Double Poisson" & 
                              Cause == "Circ") %>% 
                     pull(value))[-1])
  
  DoublePoisson_estfor <- exp(X_tot %*% par1_dp_for) -  
    exp(X_tot %*% par2_dp_for)
  
  DoublePoisson_for <- tail(DoublePoisson_estfor, h*nages)
  
  #Bivariate Poisson
  
  par1_bp_for <- c(All_eff_base$intercepts$Bivariate.Poisson[1],
                   c((All_eff_base$year_eff %>% 
                       filter(Model=="Bivariate Poisson" & 
                                Cause=="Cancer") %>% 
                       pull(value))[-1], year_eff_forec %>% 
                       filter(Model=="Bivariate Poisson") %>% 
                       pull(Cancer)),
                   (All_eff_base$age_eff %>% 
                     filter(Model=="Bivariate Poisson" & 
                              Cause=="Cancer") %>%
                     pull(value))[-1])
  
  par2_bp_for <- c(All_eff_base$intercepts$Bivariate.Poisson[2],
                   c((All_eff_base$year_eff %>% 
                       filter(Model=="Bivariate Poisson" & 
                                Cause == "Circ") %>% 
                       pull(value))[-1], year_eff_forec %>% 
                       filter(Model=="Bivariate Poisson") %>% 
                       pull(Circ)),
                   (All_eff_base$age_eff %>% 
                     filter(Model=="Bivariate Poisson" & 
                              Cause == "Circ") %>% 
                     pull(value))[-1])
  
  BivariatePoisson_estfor <- exp(X_tot %*% par1_bp_for) -  
    exp(X_tot %*% par2_bp_for)
  
  BivariatePoisson_for <- tail(BivariatePoisson_estfor, h*nages)
  
  diff_forec <- data.frame(
    "Year" = as.numeric(as.character(tail(SIM_tot$years, h*nages))),
    "Age" = as.numeric(as.character(tail(SIM_tot$ages, h*nages))),
    "Skellam" = -Skellam_for,
    "DoublePoisson" = DoublePoisson_for,
    "BivariatePoisson" = BivariatePoisson_for)
  
  
  diff_estforec <- data.frame(
    "Year" =  as.numeric(as.character(SIM_tot$years)),
    "Age" =  as.numeric(as.character(SIM_tot$ages)),
    "Skellam" = -Skellam_estfor,
    "DoublePoisson" = DoublePoisson_estfor,
    "BivariatePoisson" = BivariatePoisson_estfor)

  output <- list("out_of_sample" = diff_forec,
                 "in_AND_out_of_Sample" = diff_estforec)
  
}

plot_agesfitting <- function(check_tot){
  
 check_tot %>% 
    group_by(Year) %>% 
    pivot_longer(cols=4:6, names_to = "Model", values_to = "Est_Diff") %>% 
    mutate(Model = factor(Model, levels=c("Skellam", 
                                          "BivariatePoisson", 
                                          "DoublePoisson"))) %>% 
    filter(Age %in% seq(25, 75, 10)) %>%
    mutate(Age = paste("Age", Age)) %>% 
    ggplot(aes(x=Year, group=Model, col=Model, lty=Model))+
    geom_point(aes(y=ObservedGap, shape="Observed"),
               col="black", size=0.8)+
    geom_line(aes(y=Est_Diff), size=0.65)+
    geom_vline(xintercept=2001, lty=3, alpha=1)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust=1),
          legend.position = "bottom")+
    scale_shape_manual(values = c("Observed" = 16), name = "") +
    facet_wrap(vars(Age), scales="free")
}

heat_mape <- function(check_tot){
  require(viridis)
  
  check_tot %>% 
    mutate(Skellam = abs((Skellam - ObservedGap)/ObservedGap)*100,
           BivariatePoisson = abs((BivariatePoisson - ObservedGap)/
                                    ObservedGap)*100,
           DoublePoisson = abs((DoublePoisson - ObservedGap)/ObservedGap)*100) %>% 
    pivot_longer(cols=4:6, names_to = "Model", values_to = "Rel_Diff") %>% 
    mutate(Model = factor(Model, levels=c("Skellam", 
                                          "BivariatePoisson", 
                                          "DoublePoisson"))) %>% 
    ggplot(aes(x=Year, y=Age, z=Rel_Diff, fill=Rel_Diff))+
    geom_tile()+
    scale_fill_viridis(option = "G", discrete = F, 
                       direction=-1, na.value = "white")+
    geom_vline(xintercept=2001, lty=3, alpha=1)+
    facet_wrap(vars(Model), ncol=1)+
    theme_bw()
}

accuracy_heatmap <- function(check_for, metric){
  require(Metrics)
  require(viridis)
  
  if(metric == "RMSE"){
    check_for_by_age <- check_for %>%
    mutate(Age_group = case_when(
      Age < 25 ~ "Young",
      Age %in% 25:59 ~ "Adult",
      Age > 59 ~ "Elderly", )) %>% 
    mutate(Age_group = factor(Age_group, levels= c("Young",
                                                   "Adult", 
                                                   "Elderly"))) %>% 
    group_by(Year, Age_group) %>%
    summarise(Skellam = rmse(sum(ObservedGap), sum(Skellam)),
              DoublePoisson = rmse(sum(ObservedGap), sum(DoublePoisson)),
              BivariatePoisson = rmse(sum(ObservedGap),sum(BivariatePoisson))) %>%
    pivot_longer(cols = 3:5, names_to = "Model", values_to = "RMSE") %>% 
    mutate(Model = factor(Model, levels=c("Skellam", 
                                          "DoublePoisson", 
                                          "BivariatePoisson")))
    
    heat_map <- check_for_by_age %>% 
      mutate(Year = as.factor(Year)) %>% 
      ggplot(aes(x=Year, y=Age_group, z=RMSE, fill=RMSE))+
      geom_tile()+
      scale_fill_viridis(option = "G", discrete = F,  direction = -1)+
      facet_grid(vars(Model))+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 90, hjust=1),
            legend.position = "bottom",
            legend.key.width = unit(0.9, "cm"),
            legend.key.height = unit(0.5, "cm"))
    
    return(heat_map)
  }
  
  if(metric == "MAE"){
    check_for_by_age <- check_for %>%
      mutate(Age_group = case_when(
        Age < 25 ~ "Young",
        Age %in% 25:59 ~ "Adult",
        Age > 59 ~ "Elderly", )) %>% 
      mutate(Age_group = factor(Age_group, levels= c("Young",
                                                     "Adult", 
                                                     "Elderly"))) %>% 
      group_by(Year, Age_group) %>%
      summarise(Skellam = mae(sum(ObservedGap), sum(Skellam)),
                DoublePoisson = mae(sum(ObservedGap), sum(DoublePoisson)),
                BivariatePoisson = mae(sum(ObservedGap),sum(BivariatePoisson))) %>%
      pivot_longer(cols = 3:5, names_to = "Model", values_to = "MAE") %>% 
      mutate(Model = factor(Model, levels=c("DoublePoisson", 
                                            "BivariatePoisson",
                                            "Skellam")))
    
    heat_map <- check_for_by_age %>% 
      mutate(Year = as.factor(Year)) %>% 
      ggplot(aes(x=Year, y=Age_group, z=MAE, fill=MAE))+
      geom_tile()+
      scale_fill_viridis(option = "G", discrete = F,  direction = -1)+
      facet_grid(vars(Model))+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 90, hjust=1))+
      theme(axis.text.x = element_text(angle = 90, hjust=1),
            legend.position = "bottom",
            legend.key.width = unit(0.9, "cm"),
            legend.key.height = unit(0.5, "cm"))
    
    return(heat_map)
  }
  
  if(metric == "MAPE"){
    check_for_by_age <- check_for %>%
      mutate(Age_group = case_when(
        Age < 25 ~ "Young",
        Age %in% 25:59 ~ "Adult",
        Age > 59 ~ "Elderly", )) %>% 
      mutate(Age_group = factor(Age_group, levels= c("Young",
                                                     "Adult", 
                                                     "Elderly"))) %>% 
      group_by(Year, Age_group) %>%
      summarise(Skellam = mape(sum(ObservedGap), sum(Skellam)),
                DoublePoisson = mape(sum(ObservedGap), sum(DoublePoisson)),
                BivariatePoisson = mape(sum(ObservedGap),sum(BivariatePoisson))) %>%
      pivot_longer(cols = 3:5, names_to = "Model", values_to = "MAPE") %>% 
      mutate(Model = factor(Model, levels=c("DoublePoisson", 
                                            "BivariatePoisson",
                                            "Skellam")))
    
    heat_map <- check_for_by_age %>% 
      mutate(Year = as.factor(Year)) %>% 
      ggplot(aes(x=Year, y=Age_group, z=MAPE, fill=MAPE))+
      geom_tile()+
      scale_fill_viridis(option = "G", discrete = F,  direction = -1)+
      facet_grid(vars(Model))+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 90, hjust=1))+
      theme(axis.text.x = element_text(angle = 90, hjust=1),
            legend.position = "bottom",
            legend.key.width = unit(0.9, "cm"),
            legend.key.height = unit(0.5, "cm"))
    
    return(heat_map)
  }
  
}

DM_test <- function(check_for){
  require(forecast)
  
  diff_for <- check_for %>% 
    group_by(Year, Age) %>% 
    reframe(Skellam = ObservedGap-Skellam,
            DoublePoisson = ObservedGap-DoublePoisson,
            BivariatePoisson = ObservedGap-BivariatePoisson) %>% 
    ungroup()
  
  ages <- unique(diff_for$Age)
  
  DM_pvalue <- c()
  DM_stat <- c()
  
  alpha <- 0.05
  
  for(age in ages) {
    diff_for_age <- diff_for %>% 
      filter(Age == age)
    
    Sk_DP <- dm.test(diff_for_age$Skellam, diff_for_age$DoublePoisson)$p.value
    Sk_BP <- dm.test(diff_for_age$Skellam, diff_for_age$BivariatePoisson)$p.value
    
    Sk_DP_stat <- dm.test(diff_for_age$Skellam, diff_for_age$DoublePoisson)$statistic
    Sk_BP_stat <- dm.test(diff_for_age$Skellam, diff_for_age$BivariatePoisson)$statistic
    
    DM_pvalue_age <- data.frame("Age" = age, 
                                "DoublePoisson" = round(Sk_DP, 3),
                                "BivariatePoisson" = round(Sk_BP, 3))

    DM_stat_age <- data.frame("Age" = age, 
                              "DoublePoisson" = round(Sk_DP_stat,3), 
                              "BivariatePoisson" = round(Sk_BP_stat,3))
    
    DM_pvalue <- rbind(DM_pvalue, DM_pvalue_age)
    
    DM_stat <- rbind(DM_stat, DM_stat_age)
  }
  
  DM <- rbind(data.frame(DM_pvalue, "Type"="p.value"),
              data.frame(DM_stat, "Type"="DM.stat")) 
  
  DM_pvalue_models <- DM %>% 
    filter(Type == "p.value") %>% 
    select(-Type) %>% 
    pivot_longer(cols = 2:3, names_to = "Model", values_to = "p.value")
  
  DM_stat_models <- DM %>% 
    filter(Type == "DM.stat") %>% 
    select(-Type) %>% 
    pivot_longer(cols = 2:3, names_to = "Model", values_to = "DM_stat")
  
  DM_output <- DM_stat_models %>% 
    left_join(DM_pvalue_models, by = c("Age", "Model")) %>% 
    mutate("Significance" = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      .default = "")) %>% 
    select(Age, Model, DM_stat, p.value, Significance)
  
  return(DM_output)
  
}

