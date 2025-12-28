## ------------------------------------------------------------------- ##
##  R functions accompanying the file "Code.R".
##  Authors: Giacomo Lanfiuti Baldi & Andrea Nigri
## ------------------------------------------------------------------- ##

## Function to reproduce Figure 1: Lexis Surface of the Sex Ratio of the Age-Specific Mortality Rates.

lexis_surface <- function(apc, cou=""){
  require(viridis)
  
  brks <- quantile(c(min(apc$SR), max(apc$SR)), probs = seq(0, 1, 0.25))
  lbls <- round(brks,1)
  
  pmin=min(apc$P); pmax=max(apc$P); amin=min(apc$A); amax=max(apc$A)
  
  lexis_shape <<-  list(
    geom_vline(xintercept = seq(pmin, pmax, 10), 
               linewidth = 0.2, linetype = "dashed", 
               alpha = 0.8, color = "grey20"),
    geom_hline(yintercept = seq(amin, amax, 10), 
               linewidth = 0.2, linetype = "dashed", 
               alpha = 0.8, color = "grey20"),
    #Adding cohorts
    geom_abline(intercept = seq(-pmax, -(pmin-amax), 10), slope = 1, 
                linetype = "dashed", color = "grey20", 
                linewidth = .2, alpha = 0.8),
    coord_equal(expand = 0),
    #Adding proper labels to both axis
    scale_x_continuous(breaks = seq(pmin, pmax, 10)),
    scale_y_continuous(breaks = seq(amin, amax, 10)),
    #Adding axis titles
    labs(y = "Age", x = "Period/Cohort"),
    theme_bw())
  
  surface <- apc %>%
    ggplot(aes(x = P, y = A, z = SR))+
    geom_tile(aes(fill = SR))+
    scale_fill_viridis(option = "G", discrete = F,  direction = -1, 
                       breaks = brks, labels = lbls ,
                       name = "Sex-Ratios\n of ASDR") +
    geom_contour(bins = 10, col = "black", size = .15, alpha = 0.8)+
    labs(title=expression(paste("Lexis Surface: ", log(m[x][t]^M/m[x][t]^F))),
         subtitle=paste(cou, ": ", pmin, "-", pmax, sep=""))+
    lexis_shape
  
  print(surface)
  return(surface)
}


### Function to estimate the model. It contains the STAN model.

fit_APC_SkN <- function(apc, stanOptions){
  require(rstan)
  require(brms)
  
  #Convert the variables A, P, K into factors, so that they are treated as categorical variables.
  APC <- apc %>%
    mutate(A=as.factor(A), P=as.factor(P), K=as.factor(K))
  
  #Loading the STAN code. 
  #WARNING: it might take a while (~ 2.30 mins)
  model_Stan <- stan_model(file="fit_apc_lss_generalized.stan")
  
  print("Model correctly loaded")
  
  #Preparing data for the estimation:
  data_Stan <- c(make_standata(bf(SR ~ A + P + K, 
                                  sigma ~ 0 + as.numeric(as.character(A)), 
                                  alpha ~ 0 + as.numeric(as.character(A))),
                               data = APC, 
                               family = "skew_normal"),
                 "nage"=length(unique(apc$A)), 
                 "nperiod"=length(unique(apc$P)),
                 "ncohort"=length(unique(apc$K)))
  
  #Running the estimation
  #WARNING: it might take a while (~ 30 mins)
  fit_APC <- suppressWarnings({
    stan("fit_apc_lss_generalized.stan",
         model_code=model_Stan,
         data=data_Stan,
         chains = stanOptions$nchains, 
         iter = stanOptions$niter,
         warmup = stanOptions$nwarmup,
         cores = stanOptions$ncores,
         seed = 123)})
  
  #We save the output in a data frame that contains each parameter estimation, credibility interval, 
  output <- as.data.frame(summary(fit_APC)$summary) 
  output <- output  %>% add_column("par"=rownames(output))
  
  return(output)
}

get_effects <- function(fit){
  
  Intercept <- fit %>% 
    filter(str_detect(par, "Intercept") & !str_detect(par, "_")) %>% 
    select(mean, `2.5%`, `97.5%`)
  colnames(Intercept)<- c("Mean", "IC_L", "IC_U")
  
  b_alpha <- fit %>% 
    filter(str_detect(par, "alpha")) %>% 
    select(mean, `2.5%`, `97.5%`)
  colnames(b_alpha)<- c("Mean", "IC_L", "IC_U")
  
  b_sigma <- fit %>% 
    filter(str_detect(par, "sigma"))  %>% 
    select(mean, `2.5%`, `97.5%`)
  colnames(b_sigma)<- c("Mean", "IC_L", "IC_U")
  
  age_effect <- fit %>% 
    filter(str_detect(par, "b_1")) %>% 
    select(mean, `2.5%`, `97.5%`)
  colnames(age_effect)<- c("Mean", "IC_L", "IC_U")
  
  period_effect <- fit %>% 
    filter(str_detect(par, "b_2") & !str_detect(par, "raw")) %>% 
    select(mean, `2.5%`, `97.5%`)
  colnames(period_effect)<- c("Mean", "IC_L", "IC_U")
  
  cohort_effect <- fit %>% 
    filter(str_detect(par, "b_3")) %>% 
    select(mean, `2.5%`, `97.5%`)
  colnames(cohort_effect)<- c("Mean", "IC_L", "IC_U")
  
  effects <- list("Intercept"=Intercept, 
                  "gamma"=b_sigma,
                  "omega"=b_alpha,
                  "age_effect"=age_effect,
                  "period_effect"=period_effect,
                  "cohort_effect"=cohort_effect)
  return(effects)
}

plot_effects <- function(effects, cou=""){
  require(gridExtra)
  
  age_plot<- effects$age_effect %>%
    add_column("Age"=1:45) %>% 
    add_row("Age"=0,"Mean"=0,"IC_L"=0,"IC_U"=0,.before=1) %>% 
    pivot_longer(cols=1:3, names_to="type", values_to="value") %>% 
    ggplot()+
    geom_line(aes(x=Age, y=value, col=type, lty=type, size=type))+
    geom_hline(yintercept=0, lty=3, alpha=1)+
    scale_linetype_manual(values=c("dashed", "dashed","solid"))+
    scale_size_manual(values=c(0.5,0.5,0.8))+
    theme_bw()+
    theme(legend.position = "none")+
    labs(title=paste("Age-Period-Cohort Decomposition in ", cou))
  
  
  period_plot <- effects$period_effect %>%
    add_column("Period"=1961:(1960+length(effects$period_effect$Mean))) %>% 
    add_row("Period"=1960,"Mean"=0,"IC_L"=0,"IC_U"=0,.before=1) %>% 
    pivot_longer(cols=1:3, names_to="type", values_to="value") %>% 
    ggplot()+
    geom_line(aes(x=Period, y=value, col=type, lty=type, size=type))+
    geom_hline(yintercept=0, lty=3, alpha=1)+
    scale_linetype_manual(values=c("dashed", "dashed","solid"))+
    scale_size_manual(values=c(0.5,0.5,0.8))+
    theme_bw()
  
  cohort_plot <- effects$cohort_effect %>%
    add_column("Cohort"=1916:(1960+length(effects$period_effect$Mean))) %>% 
    add_row("Cohort"=1915,"Mean"=0,"IC_L"=0,"IC_U"=0,.before=1) %>% 
    pivot_longer(cols=1:3, names_to="type", values_to="value") %>% 
    ggplot()+
    geom_line(aes(x=Cohort, y=value,col=type, lty=type, size=type))+
    #ylim(-1,1)+
    geom_hline(yintercept=0, lty=3, alpha=1)+
    scale_linetype_manual(values=c("dashed", "dashed","solid"))+
    scale_size_manual(values=c(0.5,0.5,0.8))+
    theme_bw()+
    theme(legend.position = "none")
  
  print(grid.arrange(age_plot, period_plot, cohort_plot, ncol=1, nrow=3))
  return(grid.arrange(age_plot, period_plot, cohort_plot, ncol=1, nrow=3))
}

table_par <- function(effects, cou=""){
  
  par <- data.frame("Country" = c(cou, cou, cou),
                    "Coefficient"=c("Intercept",
                                    "gamma - Scale",
                                    "omega - Shape"),
                    round(rbind(effects$Intercept,
                                effects$gamma,
                                effects$omega),3))
  rownames(par) <- 1:nrow(par)
  
  print(par)
}

fit_APC_SkN_forecast <- function(apc, baseline, stanOptions){
  require(rstan)
  require(brms)
  
  APC_base <- apc %>% 
    filter(P %in% baseline)%>%
    mutate(A=as.factor(A), P=as.factor(P), K=as.factor(K))
  
  #Loading the STAN code. 
  #WARNING: it might take a while (~ 2.30 mins)
  model_Stan <- stan_model(file="fit_apc_lss_FORECAST.stan")
  
  print("Model correctly loaded")
  
  #Preparing data for the estimation:
  data_Stan <- c(make_standata(bf(SR ~ A + P + K, 
                                  sigma ~ 0 + as.numeric(as.character(A)), 
                                  alpha ~ 0 + as.numeric(as.character(A))),
                               data = APC_base, 
                               family = "skew_normal"),
                 "nage"=length(unique(APC_base$A)), 
                 "nperiod"=length(unique(APC_base$P)),
                 "ncohort"=length(unique(APC_base$K)))
  
  #Running the estimation
  #WARNING: it might take a while (~ 30 mins)
  suppressWarnings({
    fit_APC <- stan("fit_apc_lss_FORECAST.stan",
                    model_code=model_Stan,
                    data=data_Stan,
                    chains = stanOptions$nchains,
                    iter = stanOptions$niter,
                    warmup = stanOptions$nwarmup,
                    cores = stanOptions$ncores,
                    seed = 123)})
  
  output <- as.data.frame(summary(fit_APC)$summary) 
  output <- output  %>% add_column("par"=rownames(output))
  
  return(output)
}


SR_forecasting <- function(apc, effects, baseline, h, cou, arimaPlot=FALSE, IC=TRUE){
  require(forecast)
  require(brms)
  
  total_period <- c(baseline, (max(baseline)+1):max(baseline+h))
  
  period_pred <- as_tibble(
    forecast(auto.arima(effects$period_effect$Mean), h=h)) %>%
    select("mean"=`Point Forecast`,"L95"=`Lo 95`,"H95"=`Hi 95`)
  
  cohort_pred <- as_tibble(
    forecast(auto.arima(effects$cohort_effect$Mean), h=h)) %>%
    select("mean"=`Point Forecast`,"L95"=`Lo 95`,"H95"=`Hi 95`)
  
  age_baseline <- effects$age_effect$Mean
  period_tot <- c(effects$period_effect$Mean,period_pred$mean)
  cohort_tot <- c(effects$cohort_effect$Mean,cohort_pred$mean)
  
  apc_extended <- apc %>% 
    filter(P %in% total_period)%>%
    mutate(A=as.factor(A), P=as.factor(P), K=as.factor(K))
  
  data_forecast <- make_standata(bf(SR ~ A + P + K), data = apc_extended)               
  X <- data_forecast[["X"]]
  N <- data_forecast[["N"]]
  K <- data_forecast[["K"]]
  
  Kc = K - 1
  Xc <- matrix(nrow=N, ncol=Kc)
  means_X <- c()
  
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i])
    Xc[,i-1] = X[,i] - means_X[i - 1]
  }
  
  b_forecast <- c(age_baseline, period_tot, cohort_tot)
  SR_pred <- effects$Intercept$Mean + Xc%*%b_forecast
  
  SR_forecast <- apc_extended %>%
    mutate("SR_forecast" = SR_pred[,1]) %>% 
    mutate(P = as.integer(as.character(P)), 
           A = as.integer(as.character(A))) %>% 
    filter(P > max(baseline)) %>% 
    select(A,P,K,SR,SR_forecast)
  
  if(arimaPlot){
    arimaPlot(baseline, effects, h, cou)
  }
  
  if(IC){
    period_tot_L <- c(effects$period_effect$Mean,period_pred$L95)
    cohort_tot_L <- c(effects$cohort_effect$Mean,cohort_pred$L95)
    
    period_tot_H <- c(effects$period_effect$Mean,period_pred$H95)
    cohort_tot_H <- c(effects$cohort_effect$Mean,cohort_pred$H95)
    
    b_forecast_L <- c(age_baseline, period_tot_L, cohort_tot_L)
    SR_pred_L <- effects$Intercept$Mean + Xc%*%b_forecast_L
    
    b_forecast_H <- c(age_baseline, period_tot_H, cohort_tot_H)
    SR_pred_H <- effects$Intercept$Mean + Xc%*%b_forecast_H
    
    SR_forecast <- apc_extended %>%
      mutate("SR_forecast" = SR_pred[,1],
             "L95"=SR_pred_L[,1],
             "H95"=SR_pred_H[,1])%>% 
      mutate(P = as.integer(as.character(P)), 
             A = as.integer(as.character(A))) %>% 
      filter(P > max(baseline)) %>% 
      select(A,P,K,SR,SR_forecast, L95, H95)
  }
  
  return(SR_forecast)
}

SR_forecasting_arma <- function(apc, effects, baseline, h, IC = TRUE, nsim = 500, stepwise = FALSE, approximation = FALSE) {
  # --- helpers ---
  estimate_drift <- function(betaP, betaC, t_idx = NULL, c_idx = NULL) {
    if (is.null(t_idx)) t_idx <- seq_along(betaP)
    if (is.null(c_idx)) c_idx <- seq_along(betaC)
    bP <- coef(lm(betaP ~ t_idx))[2]
    bC <- coef(lm(betaC ~ c_idx))[2]
    list(b = as.numeric(bP + bC), bP = as.numeric(bP), bC = as.numeric(bC))
  }
  
  fit_and_forecast_d2 <- function(series_curv, h, stepwise = FALSE, approximation = FALSE) {
    d2 <- diff(series_curv, differences = 2)
    fit <- forecast::auto.arima(d2, d = 0, D = 0, seasonal = FALSE,
                                allowmean = FALSE, allowdrift = FALSE,
                                stepwise = stepwise, approximation = approximation)
    f  <- forecast::forecast(fit, h = h)
    list(fit = fit, mean = as.numeric(f$mean))
  }
  
  reconstruct_from_d2 <- function(last_two, d2_future) {
    h <- length(d2_future)
    k <- numeric(h + 2); k[1:2] <- last_two
    for (i in 3:(h + 2)) k[i] <- 2*k[i - 1] - k[i - 2] + d2_future[i - 2]
    k[3:(h + 2)]
  }
  
  simulate_d2_paths <- function(fit, h, nsim = 500, seed = 123) {
    set.seed(seed)
    replicate(nsim, as.numeric(stats::simulate(fit, nsim = h)))
  }
  
  ## 0) observed effects on baseline (posterior means)
  betaP <- as.numeric(effects$period_effect$Mean)
  betaC <- as.numeric(effects$cohort_effect$Mean)
  betaA <- as.numeric(effects$age_effect$Mean)
  intercept <- as.numeric(effects$Intercept$Mean)
  
  t_idx <- seq_along(betaP)      # 1..T
  c_idx <- seq_along(betaC)      # 1..C
  
  ## 1) total drift b=bP+bC from levels (used deterministically; allocation: bP=0, bC=b)
  dr <- estimate_drift(betaP, betaC, t_idx, c_idx)
  b  <- dr$b
  
  ## 2) curvature components (remove estimated linear part from each series)
  # NOTE: curvature is defined relative to the estimated solution; forecast law acts only on curvature.
  kappa  <- betaP - dr$bP * t_idx   # period curvature
  gamma0 <- betaC - dr$bC * c_idx   # cohort curvature
  
  ## 3) ARMA(mean=0) 
  fP <- fit_and_forecast_d2(kappa,  h, stepwise, approximation)
  fC <- fit_and_forecast_d2(gamma0, h, stepwise, approximation)
  
  ## 4) reconstruct future curvature by double summation (anchor = last two observed values)
  kappa_future <- reconstruct_from_d2(tail(kappa, 2),  fP$mean)
  gamma_future <- reconstruct_from_d2(tail(gamma0, 2), fC$mean)
  
  ## 5) rebuild future LEVELS consistent with paper:
  ##    - PERIOD: NO drift (bP=0). Keep only curvature + constant level anchored at last t.
  linP_at_last <- tail(betaP, 1) - tail(kappa, 1)        # a_P := betaP_T - kappa_T
  betaP_future <- kappa_future + as.numeric(linP_at_last) # no + b*h
  
  ##    - COHORT: ALL drift (bC=b). Anchor level at last cohort index and add b per step.
  Clast <- tail(c_idx, 1)
  linC_at_last <- tail(betaC, 1) - tail(gamma0, 1)       # a_C + b*c_T
  betaC_future <- gamma_future + as.numeric(linC_at_last) + b*(1:h)
  
  ## 6) concatenate observed + future
  period_tot <- c(betaP, betaP_future)
  cohort_tot <- c(betaC, betaC_future)
  
  ## 7) extended grid and design matrix 
  total_period <- c(baseline, (max(baseline) + 1):max(baseline + h))
  apc_extended <- apc %>%
    filter(P %in% total_period) %>%
    mutate(A = as.factor(A), P = as.factor(P), K = as.factor(K))
  
  data_forecast <- brms::make_standata(brms::bf(SR ~ A + P + K), data = apc_extended)
  X <- data_forecast[["X"]]; N <- data_forecast[["N"]]; K <- data_forecast[["K"]]
  
  # Column-centering
  Kc <- K - 1
  Xc <- matrix(nrow = N, ncol = Kc); means_X <- numeric(K)
  for (i in 2:K) { means_X[i - 1] <- mean(X[, i]); Xc[, i - 1] <- X[, i] - means_X[i - 1] }
  
  ## 8) align age block length to Xc (prepend 0 for A=0 if needed)
  needed_A <- Kc - (length(period_tot) + length(cohort_tot))
  if (needed_A != length(betaA)) {
    if (needed_A > length(betaA)) betaA <- c(rep(0, needed_A - length(betaA)), betaA)
    else if (needed_A < length(betaA)) betaA <- tail(betaA, needed_A)
  }
  stopifnot((length(betaA) + length(period_tot) + length(cohort_tot)) == ncol(Xc))
  
  ## 9) point forecasts
  b_forecast <- c(betaA, period_tot, cohort_tot)
  SR_point <- as.numeric(intercept + Xc %*% b_forecast)
  
  out <- apc_extended %>%
    mutate(SR_forecast = SR_point) %>%
    mutate(P = as.integer(as.character(P)),
           A = as.integer(as.character(A))) %>%
    filter(P > max(baseline)) %>%
    select(A, P, K, SR, SR_forecast)
  
  ## 10) predictive intervals (optional) 
  if (IC) {
    d2_sims_P <- simulate_d2_paths(fP$fit, h, nsim = nsim)
    d2_sims_C <- simulate_d2_paths(fC$fit, h, nsim = nsim)
    
    kappa_sims <- apply(d2_sims_P, 2, function(col) reconstruct_from_d2(tail(kappa, 2), col))
    gamma_sims <- apply(d2_sims_C, 2, function(col) reconstruct_from_d2(tail(gamma0, 2), col))
    
    # PERIOD paths: curvature + constant level (no drift)
    betaP_paths <- sweep(kappa_sims, 2, as.numeric(linP_at_last), `+`)
    
    # COHORT paths: curvature + anchored level + b per step
    betaC_paths <- sweep(gamma_sims, 2, as.numeric(linC_at_last), `+`) +
      matrix(b*(1:h), nrow = h, ncol = nsim, byrow = TRUE)
    
    # Assemble full coefficient matrices (observed + future) for each simulation
    period_paths <- rbind(matrix(betaP, nrow = length(betaP), ncol = nsim, byrow = FALSE),
                          betaP_paths)
    cohort_paths <- rbind(matrix(betaC, nrow = length(betaC), ncol = nsim, byrow = FALSE),
                          betaC_paths)
    
    age_mat <- matrix(betaA, nrow = length(betaA), ncol = nsim, byrow = FALSE)
    Bmat <- rbind(age_mat, period_paths, cohort_paths)  # ncol(Xc) x nsim
    
    SR_all_paths <- Xc %*% Bmat + as.numeric(intercept)
    
    mask_fore <- apc_extended %>%
      mutate(Pi = as.integer(as.character(P))) %>%
      transmute(is_fore = Pi > max(baseline)) %>%
      pull(is_fore)
    SR_fcst_paths <- SR_all_paths[mask_fore, , drop = FALSE]
    
    L95 <- apply(SR_fcst_paths, 1, stats::quantile, probs = 0.025)
    H95 <- apply(SR_fcst_paths, 1, stats::quantile, probs = 0.975)
    
    out <- out %>% mutate(L95 = as.numeric(L95),
                          H95 = as.numeric(H95))
  }
  
  return(out)
}

bapc_forecasting <- function(bapc_data, baseline, h){
  require(BAPC)
  
  ID_male <- bapc_data %>% 
  filter(Sex=="Male") %>% 
  select(c(A, P, Deaths, Pop))

  ID_female <- bapc_data %>% 
  filter(Sex=="Female")%>% 
  select(c(A, P, Deaths, Pop))

  ID_male_d <- ID_male %>%
    select(P, Deaths, A) %>% 
    pivot_wider(names_from = A, values_from = Deaths) %>%  
    as.data.frame()
  rownames(ID_male_d) <- as.character(ID_male_d$P)
  ID_male_d <- ID_male_d[,-1]
  
  ID_female_d <- ID_female %>%
    select(P, Deaths, A) %>% 
    pivot_wider(names_from = A, values_from = Deaths) %>%  
    as.data.frame()
  rownames(ID_female_d) <- as.character(ID_female_d$P)
  ID_female_d <- ID_female_d[,-1]
  
  ID_male_p <- ID_male %>%
    select(P, Pop, A) %>% 
    pivot_wider(names_from = A, values_from = Pop) %>%  
    as.data.frame()
  rownames(ID_male_p) <- as.character(ID_male_p$P)
  ID_male_p <- ID_male_p[,-1]
  
  ID_female_p <- ID_female %>%
    select(P, Pop, A) %>% 
    pivot_wider(names_from = A, values_from = Pop) %>%  
    as.data.frame()
  rownames(ID_female_p) <- as.character(ID_female_p$P)
  ID_female_p <- ID_female_p[,-1]
  
  total_period <- c(baseline, (max(baseline)+1):max(baseline+h))
  
  ID_male_d_base <- ID_male_d[row.names(ID_male_d) %in% total_period,]
  ID_male_p_base <- ID_male_p[row.names(ID_male_p) %in% total_period,]
  ID_female_d_base <- ID_female_d[row.names(ID_female_d) %in% total_period,]
  ID_female_p_base <- ID_female_p[row.names(ID_female_p) %in% total_period,]
  
  ID_male_d_base[row.names(ID_male_d_base) > max(baseline), ] <- NA
  ID_female_d_base[row.names(ID_female_d_base)> max(baseline), ] <- NA
  ID_male_p_base[row.names(ID_male_p_base) > max(baseline), ] <- 0
  ID_female_p_base[row.names(ID_female_p_base)> max(baseline), ] <- 0
  
  lc_male <- APCList(ID_male_d_base, ID_male_p_base, gf=1)
  lc_female <- APCList(ID_female_d_base, ID_female_p_base, gf=1)
  
  out_male <- qapc(BAPC(lc_male), percentiles = c(0.025, 0.5, 0.975))
  out_female <- qapc(BAPC(lc_female), percentiles = c(0.025, 0.5, 0.975))
  
  n.a <- length(unique(bapc_data$A))*5
  
  rate_male <- data.frame(agespec.rate(out_male)) %>% 
    select(seq(1,n.a,5)) %>% 
    rename("mean.0"="mean") %>% 
    mutate("Year"=rownames(.)) %>% 
    pivot_longer(cols = 1:46, names_to = "Age", values_to = "rates") %>% 
    mutate("Age"=gsub('X','',Age),"Age"=gsub('mean.','',Age)) %>% 
    rename("Male"=rates) %>% left_join(data.frame(agespec.rate(out_male)) %>% 
    select(seq(3,n.a,5)) %>% 
    rename("X0.025Q.0"="X0.025Q") %>% 
    mutate("Year"=rownames(.)) %>% 
    pivot_longer(cols = 1:46, names_to = "Age", values_to = "rates") %>% 
    mutate("Age"=gsub('X0.025Q.','',Age)) %>% 
    rename("L95.male"=rates)) %>% left_join( data.frame(agespec.rate(out_male)) %>% 
    select(seq(5,n.a,5)) %>% 
    rename("X0.975Q.0"="X0.975Q") %>% 
    mutate("Year"=rownames(.)) %>% 
    pivot_longer(cols = 1:46, names_to = "Age", values_to = "rates") %>% 
    mutate("Age"=gsub('X0.975Q.','',Age)) %>% 
    rename("H95.male"=rates))
  
  rate_female <- data.frame(agespec.rate(out_female)) %>% 
    select(seq(1,n.a,5)) %>% 
    rename("mean.0"="mean") %>% 
    mutate("Year"=rownames(.)) %>% 
    pivot_longer(cols = 1:46, names_to = "Age", values_to = "rates") %>% 
    mutate("Age"=gsub('X','',Age),"Age"=gsub('mean.','',Age)) %>% 
    rename("Female"=rates) %>% left_join(data.frame(agespec.rate(out_female)) %>% 
    select(seq(3,n.a,5)) %>% 
    rename("X0.025Q.0"="X0.025Q") %>% 
    mutate("Year"=rownames(.)) %>% 
    pivot_longer(cols = 1:46, names_to = "Age", values_to = "rates") %>% 
    mutate("Age"=gsub('X0.025Q.','',Age)) %>% 
    rename("L95.female"=rates)) %>% left_join(data.frame(agespec.rate(out_female)) %>% 
    select(seq(5,n.a,5)) %>% 
    rename("X0.975Q.0"="X0.975Q") %>% 
    mutate("Year"=rownames(.)) %>% 
    pivot_longer(cols = 1:46, names_to = "Age", values_to = "rates") %>% 
    mutate("Age"=gsub('X0.975Q.','',Age)) %>% 
    rename("H95.female"=rates))
  
  bapc_forecast <- rate_male %>% left_join(rate_female) %>% 
    mutate(SR_bapc=log(Male/Female)) %>% 
    filter(Year %in% (max(baseline)+1):(max(baseline)+h)) %>% 
    mutate(P = as.integer(Year), A = as.integer(Age)) %>%
    mutate(L95_bapc = log(L95.male/L95.female),
           H95_bapc = log(H95.male/H95.female)) %>% 
    select(A,P,SR_bapc, L95_bapc, H95_bapc)
  
  return(bapc_forecast)
}

metrics <- function(skn_forecast, bapc_forecast){
  require(Metrics)

  compare <- skn_forecast %>% 
    left_join(bapc_forecast, by=c("A","P"))
  
  PICP_skn <- round((compare %>% 
                        mutate(PICP = ifelse(SR >= L95 &
                                                    SR <= H95, 1, 0)) %>% 
                        pull(PICP) %>% sum()) / (nrow(compare)-1)*100,2) 
  
  PICP_bapc <- round((compare %>% 
    mutate(PICP = ifelse(SR >= L95_bapc &
                           SR <= H95_bapc, 1, 0)) %>% 
    pull(PICP) %>% sum()) / (nrow(compare)-1)*100,2) 
  
  model <- c("Our Model", "BAPC")
  RMSE <- round(c(rmse(compare$SR, compare$SR_forecast),
                  rmse(compare$SR, compare$SR_bapc)), 3)
  MAE <- round(c(mae(compare$SR, compare$SR_forecast),
                 mae(compare$SR, compare$SR_bapc)), 3)
  PICP <- c(PICP_skn, PICP_bapc)
  
  metric <- data.frame(model, RMSE, MAE, PICP)
  return(metric)    
}

#Supplementary Material

arimaPlot <- function(baseline, effects, h, cou){
  require(gridExtra)
  
  total_period <- c(baseline, (max(baseline)+1):max(baseline+h))
  
  arima_period <- auto.arima(effects$period_effect$Mean)
  arima_cohort <- auto.arima(effects$cohort_effect$Mean)
  
  period_pred <- as_tibble(
    forecast(arima_period, h=h)) %>%
    select("mean"=`Point Forecast`,"L95"=`Lo 95`,"H95"=`Hi 95`)
  
  cohort_pred <- as_tibble(
    forecast(arima_cohort, h=h)) %>%
    select("mean"=`Point Forecast`,"L95"=`Lo 95`,"H95"=`Hi 95`)
  
  plot_data_period <- data.frame(
    year = total_period,
    value = c(0,effects$period_effect$Mean, 
              period_pred$mean),
    upper = c(rep(NA, length(effects$period_effect$Mean)+1),
              period_pred$H95),
    lower = c(rep(NA, length(effects$period_effect$Mean)+1),
              period_pred$L95),
    type = c(rep("Actual", length(effects$period_effect$Mean)+1),
             rep("Forecast", length(period_pred$mean))))
  
  total_cohort <- c((min(baseline)-max(apc$A)):max(baseline+h))
  
  plot_data_cohort <- data.frame(
    year = total_cohort,
    value = c(0,effects$cohort_effect$Mean, cohort_pred$mean),
    upper = c(rep(NA, length(effects$cohort_effect$Mean)+1), 
              cohort_pred$H95),
    lower = c(rep(NA, length(effects$cohort_effect$Mean)+1),
              cohort_pred$L95),
    type = c(rep("Actual", length(effects$cohort_effect$Mean)+1), 
             rep("Forecast", length(cohort_pred$mean))))
  
  order_period <- paste(as.character(arimaorder(arima_period)), collapse=", ")
  order_cohort <- paste(as.character(arimaorder(arima_cohort)), collapse=", ")
  
  plot_period <- ggplot(plot_data_period, 
                        aes(x = year, y = value, color = type)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = "Forecast"),
                alpha = 0.3) +
    labs(title = paste("Period effects - ARIMA(", order_period,")", sep = ""),
         x = "year", y = expression(beta^P)) +
    theme_minimal() +
    scale_fill_manual(values = c("Forecast" = "blue"))+
    theme(legend.position = "none",
          plot.title = element_text(size = 12))
  
  plot_cohort <- ggplot(plot_data_cohort, 
                        aes(x = year, y = value, color = type)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = "Forecast"), 
                alpha = 0.3) +
    labs(title = paste("Cohort effects - ARIMA(", order_cohort,")", sep = ""),
         x = "year", y = expression(beta^C)) +
    theme_minimal() +
    scale_fill_manual(values = c("Forecast" = "blue"))+
    theme(legend.position = "none",
          plot.title = element_text(size = 12))
  
  forecast_window <- max(baseline+1):(max(baseline)+h)
  
  print(grid.arrange(plot_period, plot_cohort, ncol=2, nrow=1))
                     #top = paste("Projection of the effects in ",
                     #           cou," in period ",
                     #             min(forecast_window),"-",max(forecast_window),sep#="")))
  
  ggsave(paste0("Figures/SM_Figure1_", cou," in ",
               min(forecast_window),"-",max(forecast_window),".png"), 
         plot = grid.arrange(plot_period, plot_cohort, ncol=2, nrow=1), 
         device = "png",
         width = 20, height = 6, units = "cm")
}

heat_map <- function(skn_forecast, bapc_forecast){
  compare <- skn_forecast %>%
    left_join(bapc_forecast, by=c("A","P"))
  
  heatmap <- compare %>% 
    pivot_longer(cols=c(SR_forecast, SR_bapc), 
                 values_to = "SR_est",
                 names_to = "Model") %>% 
    mutate(Model = case_when(
      Model == "SR_forecast" ~ "Our model",
      Model == "SR_bapc" ~ "BAPC",)) %>% 
    mutate(Diff = SR-SR_est, 
           P = as.factor(P),
           Model = as_factor(Model)) %>% 
    ggplot(aes(x=P, y=A, z=Diff))+
    geom_raster(aes(fill=Diff),interpolate = TRUE)+
    scale_fill_gradientn(colours = c("red", "white", "blue"), 
                         limits = c(-2, 2), 
                         breaks = c(-2, 0, 2), 
                         na.value = "grey50") +
    theme_bw()+
    ylab("Age")+xlab("Years")+
    facet_wrap(vars(Model),scales = "free")+
    ggtitle("Heat map of differences")
  
  print(heatmap)
}
