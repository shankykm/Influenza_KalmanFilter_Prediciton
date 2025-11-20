library(tidyverse)
library(rjags)
library(R2jags)

# =============================================================================
# Import historical flu data  with smoothed covariate
# import informative priors
# =============================================================================
flu_data <- read.csv("flu_data.csv", header=TRUE)
simul_kf_ensem_informativePriors_agg_df <- read.csv("F://CCF Files//CCF_Old_Laptop_shortcut//CCF_DeskopFiles//Books_Articles_litreatures//GraduateDegree//OSU//Thesis//Short term Flu Prediction Model//modelOutputData//Simulation//simul_kf_ensem_informativePriors_agg_df.csv", header=TRUE)
# simul_kf_ensem_informativePriors_agg_df <- read.csv("simul_kf_ensem_informativePriors_agg_df.csv", header=TRUE)

simul_kf_bivar_diff_informativePriors_agg_df <- read.csv("F://CCF Files//CCF_Old_Laptop_shortcut//CCF_DeskopFiles//Books_Articles_litreatures//GraduateDegree//OSU//Thesis//Short term Flu Prediction Model//modelOutputData//Simulation//simul_kf_bivar_diff_informativePriors_agg_df.csv", header=TRUE)
# simul_kf_bivar_diff_informativePriors_agg_df <- read.csv("simul_kf_bivar_diff_informativePriors_agg_df.csv", header=TRUE)

# =============================================================================
# Total Lab test forecast using KF Ensemble Model with informative priors
# =============================================================================
kf_ensemble <- "
model {
  # Model1 informative priors for kfar1guassJags_covariance_6.
  theta_1[1,1] ~  dunif(theta_1_1.hyper.a, theta_1_1.hyper.b)
  theta_1[1,2] <- 0
  theta_1[2,1] <- 0
  theta_1[2,2] ~  dunif(theta_2_2.hyper.a, theta_2_2.hyper.b)

  state.sigma[1] ~ dunif(state.sigma.1.a.hyper, state.sigma.1.b.hyper)
  state.sigma[2] ~ dunif(state.sigma.2.a.hyper, state.sigma.2.b.hyper)
  obs.sigma[1] ~ dunif(obs.sigma.1.a.hyper, obs.sigma.1.b.hyper)
  obs.sigma[2] ~ dunif(obs.sigma.2.a.hyper, obs.sigma.2.b.hyper)

  state.rho ~ dunif(state.rho.a.hyper, state.rho.b.hyper)
  obs.rho ~ dunif(obs.rho.a.hyper, obs.rho.b.hyper)

  state.cov[1,1] <- state.sigma[1] * state.sigma[1]
  state.cov[1,2] <- state.sigma[1] * state.sigma[2] * state.rho
  state.cov[2,1] <- state.sigma[1] * state.sigma[2] * state.rho
  state.cov[2,2] <- state.sigma[2] * state.sigma[2]
  state.prec[1:2,1:2] <- inverse(state.cov[,])

  obs.cov[1,1] <- obs.sigma[1] * obs.sigma[1]
  obs.cov[1,2] <- obs.sigma[1] * obs.sigma[2] * obs.rho
  obs.cov[2,1] <- obs.sigma[1] * obs.sigma[2] * obs.rho
  obs.cov[2,2] <- obs.sigma[2] * obs.sigma[2]
  obs.prec[1:2,1:2] <- inverse(obs.cov[,])

  ### state process equation
  # 1st time step (based on a known value)
  z[1,1] <- z1_init
  z[1,2] <- z2_init
  z11[1] <- z1_init
  z21[1] <- z2_init

  p[1, 1, 1] <- p1_init
  p[1, 1, 2] <- 1
  p[1, 2, 1] <- 1
  p[1, 2, 2] <- p2_init
  #define 2*2 identity matrix
  identity_matrix[1, 1] <- 1
  identity_matrix[2, 2] <- 1
  identity_matrix[1, 2] <- 0
  identity_matrix[2, 1] <- 0

  # Model2 informative Priors for kfar1covguass_Jags4 model

  theta2_1 ~  dunif(theta2_1.hyper.a, theta2_1.hyper.b)

  inv.state.variance ~ dgamma(inv.state.variance.hyper.a, inv.state.variance.hyper.b)
  inv.obs.variance.1  ~ dgamma(inv.obs.variance.1.hyper.a, inv.obs.variance.1.hyper.b)

  # Transform inv.variance to variance
  state.variance <- 1/(inv.state.variance)
  obs.variance.1 <- 1/(inv.obs.variance.1)
  p12[1] <- p1_init
  z12[1] <- z1_init

  # informative priors for ensemble weight
  w1 ~  dunif(w1.hyper.a, w1.hyper.b)
  z_ensemble[1] <- z1_init
  inv.obs.variance.2  ~ dgamma(inv.obs.variance.2.hyper.a, inv.obs.variance.2.hyper.b)
  obs.variance.2 <- 1/(inv.obs.variance.2)

  # Remaining time steps
  for(k in 2:TT){
    # compute kalman gain for the previous step because we would have observed t-1
    kg[k-1, 1:2, 1:2] <- p[k-1, , ] %*% inverse(p[k-1, , ] + obs.cov[,])

    # update the predictions made for previous step, t-1, by kalman gain and observed value for t-1
    z_update[k, 1:2] <- z[k-1, ] + c( t( kg[k-1, , ] %*% t(t(y[k-1,] - z[k-1,])) ) )
    p_update[k, 1:2, 1:2] <- ( (identity_matrix[,] - kg[k-1, , ]) %*% p[k-1, , ] %*%
                             t((identity_matrix[,] - kg[k-1, , ])) )  +
                             (kg[k-1, , ] %*% obs.cov[,] %*% t(kg[k-1, , ]))

    # predict the current step, t, value using the updated previous step predictions
    p[k, 1:2, 1:2] <- state.cov[,] + (theta_1[,] %*% p_update[k, , ] %*% t(theta_1[,]) )
    mu[k,1:2] <- c( t( theta_1[,] %*% t(t(z_update[k,1:2])) ) )
    z[k,1:2] ~ dmnorm( mu[k,1:2], inverse(p[k, 1:2, 1:2]) )
    z11[k] <- z[k,1]
    z21[k] <- z[k,2]

    # Model2 definition for kfar1covguass_Jags3
    # compute kalman gain for the previous step because we would have observed t-1
    k12[k-1] <- p12[k-1]/(p12[k-1] + obs.variance.1)

    # update the predictions made for previous step, t-1, by kalman gain and observed value for t-1
    z12_update[k] <- z12[k-1] + (k12[k-1] * (y12[k-1] - z12[k-1]))
    p12_update[k] <- (1-k12[k-1])*p12[k-1]

    # predict the current step, t, value using the updated previous step predictions
    p12[k] <- state.variance+(p12_update[k]*theta2_1^2)
    z12[k] ~ dnorm(z12_update[k]*theta2_1, 1/p12[k])

    # weighted average from two models
    z_ensemble[k] <- (w1*z11[k]) + ( (1-w1)*z12[k] )
  }

  k12[TT] <- k12[TT-1]
  kg[TT, 1:2, 1:2] <- kg[TT-1, 1:2, 1:2]

  # Observation equations
  for(j in 1:TT) {
    y[j,1:2] ~ dmnorm(z[j,1:2], inverse(obs.cov[,]))

    y12[j] ~ dnorm(z12[j], 1/obs.variance.1)

    yensemb[j] ~ dnorm(z_ensemble[j], 1/obs.variance.2)
  }
}"

kf_ensem_predict_informPrior <- function(data, informPrior_df,
                                         dependent_var1="", dependent_var2="",
                                         covar1 = "",
                                         test_Start, n, h, n_ahead,
                                         seasonalYear, last_id,
                                         modelname="kf_ensemble",
                                         z0=0, nc=3, burn=1000, infer=10000){
  
  forecasted_values_df <- data.frame()
  
  prior_df <- informPrior_df %>%
    filter(as.numeric(seasonal_Year) == seasonalYear &
             modelname == substr(model_name,1,nchar(model_name)-24))
  
  theta_1_1.hyper.a <- prior_df$theta_1_1.hyper.a
  theta_1_1.hyper.b <- prior_df$theta_1_1.hyper.b
  theta_2_2.hyper.a <- prior_df$theta_2_2.hyper.a
  theta_2_2.hyper.b <- prior_df$theta_2_2.hyper.b
  theta2_1.hyper.a <- prior_df$theta2_1.hyper.a
  theta2_1.hyper.b <- prior_df$theta2_1.hyper.b
  
  inv.state.variance.hyper.a <- prior_df$inv.state.variance.hyper.a
  inv.state.variance.hyper.b <- prior_df$inv.state.variance.hyper.b
  inv.obs.variance.1.hyper.a <- prior_df$inv.obs.variance.1.hyper.a
  inv.obs.variance.1.hyper.b <- prior_df$inv.obs.variance.1.hyper.b
  inv.obs.variance.2.hyper.a <- prior_df$inv.obs.variance.2.hyper.a
  inv.obs.variance.2.hyper.b <- prior_df$inv.obs.variance.2.hyper.b
  
  state.sigma.1.a.hyper <- prior_df$state.sigma.1.a.hyper
  state.sigma.1.b.hyper <- prior_df$state.sigma.1.b.hyper
  obs.sigma.1.a.hyper <- prior_df$obs.sigma.1.a.hyper
  obs.sigma.1.b.hyper <- prior_df$obs.sigma.1.b.hyper
  state.sigma.2.a.hyper <- prior_df$state.sigma.2.a.hyper
  state.sigma.2.b.hyper <- prior_df$state.sigma.2.b.hyper
  obs.sigma.2.a.hyper <- prior_df$obs.sigma.2.a.hyper
  obs.sigma.2.b.hyper <- prior_df$obs.sigma.2.b.hyper
  obs.rho.a.hyper <- prior_df$obs.rho.a.hyper
  obs.rho.b.hyper <- prior_df$obs.rho.b.hyper
  state.rho.a.hyper <- prior_df$state.rho.a.hyper
  state.rho.b.hyper <- prior_df$state.rho.b.hyper
  
  w1.hyper.a <- prior_df$w1.hyper.a
  w1.hyper.b <- prior_df$w1.hyper.b

    for (i in test_Start:test_Start) { # Not looping through all future weeks

    print(paste(seasonalYear, modelname, i, sep=","))
    
    last_id <- last_id + 1
    if((i+n_ahead-1) > n){
      end_index <- n-i
    }else{
      end_index <- 3
    }
    
    col1 <- ts(data[1:(i-1), dependent_var1], start=1)
    col1 <- c(col1, rep(NA, 1+end_index))
    
    col2 <- ts(data[1:(i-1), dependent_var2], start=1)
    col2 <- c(col2, rep(NA, 1+end_index))
    
    TT <- length(col1)
    model <- get(noquote(modelname))
    
    dataar1Jags <- list(TT=TT,
                        y=cbind(matrix(col1, nrow=length(col1)),
                                matrix(col2, nrow=length(col2))),
                        y12=col1, yensemb=col1,
                        z1_init=mean(col1, na.rm=TRUE), z2_init=mean(col2, na.rm=TRUE),
                        p1_init=10, p2_init=1,
                        inv.state.variance.hyper.a=inv.state.variance.hyper.a,
                        inv.state.variance.hyper.b=inv.state.variance.hyper.b,
                        inv.obs.variance.1.hyper.a=inv.obs.variance.1.hyper.a,
                        inv.obs.variance.1.hyper.b=inv.obs.variance.1.hyper.b,
                        inv.obs.variance.2.hyper.a=inv.obs.variance.2.hyper.a,
                        inv.obs.variance.2.hyper.b=inv.obs.variance.2.hyper.b,
                        theta_1_1.hyper.a=theta_1_1.hyper.a,
                        theta_1_1.hyper.b=theta_1_1.hyper.b,
                        theta_2_2.hyper.a=theta_2_2.hyper.a,
                        theta_2_2.hyper.b=theta_2_2.hyper.b,
                        theta2_1.hyper.a=theta2_1.hyper.a,
                        theta2_1.hyper.b=theta2_1.hyper.b,
                        state.sigma.1.a.hyper=state.sigma.1.a.hyper,
                        state.sigma.1.b.hyper=state.sigma.1.b.hyper,
                        state.sigma.2.a.hyper=state.sigma.2.a.hyper,
                        state.sigma.2.b.hyper=state.sigma.2.b.hyper,
                        obs.sigma.1.a.hyper=obs.sigma.1.a.hyper,
                        obs.sigma.1.b.hyper=obs.sigma.1.b.hyper,
                        obs.sigma.2.a.hyper=obs.sigma.2.a.hyper,
                        obs.sigma.2.b.hyper=obs.sigma.2.b.hyper,
                        obs.rho.a.hyper=obs.rho.a.hyper,
                        obs.rho.b.hyper=obs.rho.b.hyper,
                        state.rho.a.hyper=state.rho.a.hyper,
                        state.rho.b.hyper=state.rho.b.hyper,
                        w1.hyper.a=w1.hyper.a,
                        w1.hyper.b=w1.hyper.b
    )
    
    print(dataar1Jags)
    ncar1Jags <- nc
    burnar1Jags <- burn
    inferar1Jags <- infer
    mar1Jags <- jags.model(file=textConnection(object=model),
                           data=dataar1Jags,
                           n.chains=ncar1Jags)
    
    update(mar1Jags, n.iter=burnar1Jags)
    
    jagsar1VarM <- c("z11", "kg", "z21", "z12", "k12", "w1", "z_ensemble")
    far1Jags <-  coda.samples(mar1Jags,
                              variable.names=jagsar1VarM,
                              n.iter=inferar1Jags)
    
    fit_total_vol_ensem <- summary(far1Jags[,paste0("z_ensemble[",1:TT,"]"),])$statistics[,"Mean"][i:(i+end_index)]
    CI_total_vol_ensem <- data.frame(summary(far1Jags[,paste0("z_ensemble[",1:TT,"]"),])$quantiles[,c(1, 5)])[i:(i+end_index),]
    
    fit_total_vol_1 <- summary(far1Jags[,paste0("z11[",1:TT,"]"),])$statistics[,"Mean"][i:(i+end_index)]
    CI_total_vol_1<- data.frame(summary(far1Jags[,paste0("z11[",1:TT,"]"),])$quantiles[,c(1, 5)])[i:(i+end_index),]
    
    fit_total_vol_2 <- summary(far1Jags[,paste0("z12[",1:TT,"]"),])$statistics[,"Mean"][i:(i+end_index)]
    CI_total_vol_2<- data.frame(summary(far1Jags[,paste0("z12[",1:TT,"]"),])$quantiles[,c(1, 5)])[i:(i+end_index),]
    
    w1 <- unname(summary(far1Jags[,paste0("w1"),])$statistics["Mean"])
    forecast <- setNames(data.frame(fit_total_vol_ensem, CI_total_vol_ensem,
                                    fit_total_vol_1, CI_total_vol_1,
                                    fit_total_vol_2, CI_total_vol_2,
                                    rep(w1, end_index+1),
                                    rep(seasonalYear, end_index+1),
                                    rep(i-1, end_index+1), rep(n, end_index+1),
                                    seq(i, i+end_index),
                                    rep(n_ahead, end_index+1)
    ),
    c("fit_total_vol_ensem", "lwr_total_vol_ensem", "upr_total_vol_ensem",
      "fit_total_vol_1", "lwr_total_vol_1", "upr_total_vol_1",
      "fit_total_vol_2", "lwr_total_vol_2", "upr_total_vol_2",
      "w1",
      "seasonal_Year", "forecast_origin_id",
      "last_obs_id", "id", "n_ahead"))
    
    forecast <- replace(forecast, forecast < 0, 0)
    
    rownames(forecast) <- NULL
    forecasted_values_df <- rbind(forecasted_values_df, forecast)
    
  }
  data.frame(forecasted_values_df)
}

# =============================================================================
# Ensemble KF Rolling 4 Weeks ahead forecast using informative priors.
# =============================================================================
set.seed(100)
h=20  # set a specific starting week for your forecast. h=50 will start predicting at week 3.
n_ahead <- 4
modelname <- c("kf_ensemble")
var1 <- "total_vol"
var2 <- "Detected"
covar1 <- "smoothed_total_vol.fit"
kf_multi_step_ensem_informPrior_forecast <- data.frame()
for (m in modelname){
  kf_multi_step_ensem_inforPrior_predict <- flu_data %>%
    filter(year %in% c(2019) & !week==53) %>%
    group_by(year) %>%
    arrange(week) %>%
    mutate(n=n(), test_start = n - h + 1, h=h) %>%
    do(kf_ensem_predict_informPrior(., simul_kf_ensem_informativePriors_agg_df,
                                    var1, var2, "covar1",
                                    first(.$test_start), first(.$n),
                                    first(.$h), n_ahead,
                                    first(.$year),
                                    last(.$week)-h,
                                    modelname=m,
                                    z0=1, nc=3, burn=1000, infer=10000)
    ) %>%
    mutate(forecast_step = id - forecast_origin_id,
           model = m,
           dep_var = "total_vol",
           Prior = "informative")
  
  kf_multi_step_ensem_informPrior_forecast <- rbind(kf_multi_step_ensem_inforPrior_predict,
                                                    kf_multi_step_ensem_informPrior_forecast)
}

simul_kf_multi_step_ensem_informPrior_forecast_v2 <- inner_join(flu_data[c("year", "week",
                                                                           "total_vol", "Detected")],
                                                                kf_multi_step_ensem_informPrior_forecast,
                                                                by = c("year" = "seasonal_Year",
                                                                       "week" = "id"))

# =============================================================================
# Positive Cases - Bivariate KF Diff Model definition using informative priors
# =============================================================================
bivar_kfar1gaussJags8_2 <- "
model {
   # informative Priors
    for (i in 1:4) {
      theta_1.precision.hyper[i] <- 1/(theta_1.sd.hyper[i]*theta_1.sd.hyper[i])
      theta_1[i] ~  dnorm(theta_1.mean.hyper[i], theta_1.precision.hyper[i])
    }

   inv.state.variance.1 ~ dgamma(inv.state.variance.1.hyper.a,
                                 inv.state.variance.1.hyper.b)
   inv.obs.variance.1  ~ dgamma(inv.obs.variance.1.hyper.a,
                                inv.obs.variance.1.hyper.b)
   inv.state.variance.2 ~ dgamma(inv.state.variance.2.hyper.a,
                                 inv.state.variance.2.hyper.b)
   inv.obs.variance.2  ~ dgamma(inv.obs.variance.2.hyper.a,
                                inv.obs.variance.2.hyper.b)

   # Transform inv.precision to variance
   state.variance.1 <- 1/(inv.state.variance.1)
   obs.variance.1 <- 1/(inv.obs.variance.1)
   state.variance.2 <- 1/(inv.state.variance.2)
   obs.variance.2 <- 1/(inv.obs.variance.2)

   ### state process equation
   # 1st time step (based on known value)
   z1[1] <- z1_init
   z2[1] <- z2_init
   p1[1] <- p1_init
   p2[1] <- p2_init

  # Remaining time steps
  for(t in 2:TT){
    # compute kalman gain for the previous step because we would have observed t-1
    k1[t-1] <- p1[t-1]/(p1[t-1] + obs.variance.1)
    k2[t-1] <- p2[t-1]/(p2[t-1] + obs.variance.2)

    # update the predictions made for previous step, t-1, by kalman gain and observed value for t-1
    z1_update[t] <- z1[t-1] + (k1[t-1]*(y[t-1,1] - z1[t-1]))
    p1_update[t] <- (1-k1[t-1])*p1[t-1]

    z2_update[t] <- z2[t-1] + (k2[t-1]*(y[t-1,2] - z2[t-1]))
    p2_update[t] <- (1-k2[t-1])*p2[t-1]

    # predict the current step, t, value using the updated previous step predictions
    p1[t] <- state.variance.1 + (p1_update[t]*theta_1[1]^2) + (p2_update[t]*theta_1[2]^2)
    z1[t] ~ dnorm(z1_update[t]*theta_1[1] + z2_update[t]*theta_1[2], 1/p1[t])

    p2[t] <- state.variance.2 + (p2_update[t]*theta_1[4]^2)
    z2[t] ~ dnorm(cov1[t]*theta_1[3] + z2_update[t]*theta_1[4], 1/p2[t])

  }

  k1[TT] <- k1[TT-1]
  k2[TT] <- k2[TT-1]

  # Observation equations
  for(t in 1:TT) {
    y[t,1] ~ dnorm(z1[t], 1/obs.variance.1)
    y[t,2] ~ dnorm(z2[t], 1/obs.variance.2)
  }
}"

kf_bivar_predict_informPrior <- function(data, informPrior_df,
                                         dependent_var1="", dependent_var2="",
                                         covar1 = "",
                                         test_Start, n, h, n_ahead,
                                         seasonalYear, last_id,
                                         modelname="bivar_kfar1gaussJags8_2",
                                         z0=0, nc=3, burn=1000, infer=10000){
  
  forecasted_values_df <- data.frame()
  
  prior_df <- informPrior_df %>%
    filter(as.numeric(seasonal_Year) == seasonalYear &
             modelname == substr(model_name,1,nchar(model_name)-24))
  
  theta_1.mean.hyper <- c(prior_df$theta_1_1_fit_mean_Mean, prior_df$theta_1_2_fit_mean_Mean,
                          prior_df$theta_2_1_fit_mean_Mean, prior_df$theta_2_2_fit_mean_Mean
  )
  
  theta_1.sd.hyper <- c(prior_df$theta_1_1_fit_mean_SD, prior_df$theta_1_2_fit_mean_SD,
                        prior_df$theta_2_1_fit_mean_SD, prior_df$theta_2_2_fit_mean_SD)
  
  inv.state.variance.1.hyper.a <- prior_df$inv.state.variance.1.hyper.a
  inv.state.variance.1.hyper.b <- prior_df$inv.state.variance.1.hyper.b
  inv.obs.variance.1.hyper.a <- prior_df$inv.obs.variance.1.hyper.a
  inv.obs.variance.1.hyper.b <- prior_df$inv.obs.variance.1.hyper.b
  
  inv.state.variance.2.hyper.a <- prior_df$inv.state.variance.2.hyper.a
  inv.state.variance.2.hyper.b <- prior_df$inv.state.variance.2.hyper.b
  inv.obs.variance.2.hyper.a <- prior_df$inv.obs.variance.2.hyper.a
  inv.obs.variance.2.hyper.b <- prior_df$inv.obs.variance.2.hyper.b
  
  for (i in test_Start:test_Start) {
    
    print(paste(seasonalYear, modelname, i, sep=","))
    
    last_id <- last_id + 1
    if((i+n_ahead-1) > n){
      end_index <- n-i
    }else{
      end_index <- 3
    }
    
    col1 <- ts(data[1:(i-1), dependent_var1], start=1)
    col1 <- c(col1, rep(NA, 1+end_index))
    
    col2 <- ts(data[1:(i-1), dependent_var2], start=1)
    col2 <- c(col2, rep(NA, 1+end_index))
    
    cov1 <- c(ts(data[1:(i+end_index), covar1], start=1))
    cov2 <- c(ts(data[1:(i+end_index), dependent_var1], start=1))
    
    TT <- length(col1)
    model <- get(noquote(modelname))
    
    dataar1Jags <- list(TT=TT,
                        y=cbind(matrix(col1, nrow=length(col1)),
                                matrix(col2, nrow=length(col2))),
                        inv.state.variance.1.hyper.a=inv.state.variance.1.hyper.a,
                        inv.state.variance.1.hyper.b=inv.state.variance.1.hyper.b,
                        inv.obs.variance.1.hyper.a=inv.obs.variance.1.hyper.a,
                        inv.obs.variance.1.hyper.b=inv.obs.variance.1.hyper.b,
                        inv.state.variance.2.hyper.a=inv.state.variance.2.hyper.a,
                        inv.state.variance.2.hyper.b=inv.state.variance.2.hyper.b,
                        inv.obs.variance.2.hyper.a=inv.obs.variance.2.hyper.a,
                        inv.obs.variance.2.hyper.b=inv.obs.variance.2.hyper.b,
                        theta_1.mean.hyper=theta_1.mean.hyper,
                        theta_1.sd.hyper=theta_1.sd.hyper,
                        cov1=cov1, cov2=cov2,
                        z1_init=10, z2_init=1,
                        p1_init=10, p2_init=1
    )
    
    print(dataar1Jags)
    ncar1Jags <- nc
    burnar1Jags <- burn
    inferar1Jags <- infer
    mar1Jags <- jags.model(file=textConnection(object=model),
                           data=dataar1Jags,
                           n.chains=ncar1Jags)
    
    update(mar1Jags, n.iter=burnar1Jags)
    
    jagsar1VarM <- c("z1", "p1", "k1", "z2", "p2", "k2")
    far1Jags <-  coda.samples(mar1Jags,
                              variable.names=jagsar1VarM,
                              n.iter=inferar1Jags)
    
    fit_total_vol <- summary(far1Jags[,paste0("z1[",1:TT,"]"),])$statistics[,"Mean"][i:(i+end_index)]
    fit_Detected <- summary(far1Jags[,paste0("z2[",1:TT,"]"),])$statistics[,"Mean"][i:(i+end_index)]
    CI_total_vol <- data.frame(summary(far1Jags[,paste0("z1[",1:TT,"]"),])$quantiles[,c(1, 5)])[i:(i+end_index),]
    CI_Detected <- data.frame(summary(far1Jags[,paste0("z2[",1:TT,"]"),])$quantiles[,c(1, 5)])[i:(i+end_index),]
    
    Kalman_gain_total_vol <- unname(summary(far1Jags[,paste0("k1[",1:TT,"]"),])$statistics[,"Mean"][i:(i+end_index)])
    Kalman_gain_Detected <- unname(summary(far1Jags[,paste0("k2[",1:TT,"]"),])$statistics[,"Mean"][i:(i+end_index)])
    forecast <- setNames(data.frame(fit_total_vol, CI_total_vol, Kalman_gain_total_vol,
                                    fit_Detected, CI_Detected, Kalman_gain_Detected,
                                    rep(seasonalYear, end_index+1),
                                    rep(i-1, end_index+1), rep(n, end_index+1),
                                    seq(i, i+end_index),
                                    rep(n_ahead, end_index+1)
    ),
    c("fit_total_vol", "lwr_total_vol", "upr_total_vol", "Kalman_gain_total_vol",
      "fit_Detected", "lwr_Detected", "upr_Detected", "Kalman_gain_Detected",
      "seasonal_Year", "forecast_origin_id",
      "last_obs_id", "id", "n_ahead"))
    
    forecast <- replace(forecast, forecast < 0, 0)
    rownames(forecast) <- NULL
    forecasted_values_df <- rbind(forecasted_values_df, forecast)
    
  }
  data.frame(forecasted_values_df)
}

# =============================================================================
# KF Diff Rolling 4 Weeks ahead forecast using informative priors.
# =============================================================================
flu_data <- flu_data <- flu_data %>%
  arrange(year, week) %>%
  mutate(lag_total_vol = lag(total_vol),
         lag_Detected = lag(Detected),
         diff_total_vol = c(NA, diff(total_vol)),
         diff_Detected = c(NA, diff(Detected)))

set.seed(100)
h=20 # set a specific starting week for your forecast. h=50 will start predicting at week 3.
n_ahead <- 4 
nc <- 3 # MCMC chain
modelname <- c("bivar_kfar1gaussJags8_2")
var1 <- "diff_total_vol"
var2 <- "diff_Detected"
covar1 <- "diff_smoothed_total_vol.fit"
simul_kf_multi_step_diff_bivar_informPrior_forecast <- data.frame()
for (m in modelname){
  kf_multi_step_bivar_inforPrior_predict <- flu_data %>%
    filter(year %in% c(2019) ) %>%
    group_by(year) %>%
    arrange(week) %>%
    mutate(n=n(), test_start = n - h + 1, h=h) %>%
    do(kf_bivar_predict_informPrior(., simul_kf_bivar_diff_informativePriors_agg_df,
                                    var1, var2, covar1,
                                    first(.$test_start), first(.$n),
                                    first(.$h), n_ahead,
                                    first(.$year),
                                    last(.$week)-h,
                                    modelname=m,
                                    z0=1, nc=3, burn=1000, infer=10000)
    ) %>%
    mutate(forecast_step = id - forecast_origin_id,
           model = m,
           dep_var = "bivar",
           Prior = "informative")
  
  simul_kf_multi_step_diff_bivar_informPrior_forecast <- rbind(kf_multi_step_bivar_inforPrior_predict,
                                                               simul_kf_multi_step_diff_bivar_informPrior_forecast)
}

simul_kf_multi_step_diff_bivar_informPrior_forecast_v2 <- inner_join(flu_data[c("year", "week",
                                                                                "total_vol", "Detected",
                                                                                "diff_total_vol", "diff_Detected",
                                                                                "lag_total_vol", "lag_Detected")],
                                                                     simul_kf_multi_step_diff_bivar_informPrior_forecast,
                                                                     by = c("year" = "year",
                                                                            "week" = "id"))

simul_kf_multi_step_diff_bivar_informPrior_forecast_v3 <- simul_kf_multi_step_diff_bivar_informPrior_forecast_v2 %>%  arrange(year, week) %>%
  mutate(diff_fit_total_vol = fit_total_vol,
         diff_lwr_total_vol = lwr_total_vol,
         diff_upr_total_vol = upr_total_vol,
         diff_fit_Detected = fit_Detected,
         diff_lwr_Detected = lwr_Detected,
         diff_upr_Detected = upr_Detected,
         model = paste(model, "diff", sep="_"),
         
         fit_total_vol_new = ifelse(forecast_step==1,diff_fit_total_vol+lag_total_vol,NA),
         lwr_total_vol_new = ifelse(forecast_step==1,diff_lwr_total_vol+lag_total_vol,NA),
         upr_total_vol_new = ifelse(forecast_step==1,diff_upr_total_vol+lag_total_vol,NA),
         
         fit_Detected_new = ifelse(forecast_step==1,diff_fit_Detected+lag_Detected,NA),
         lwr_Detected_new = ifelse(forecast_step==1,diff_lwr_Detected+lag_Detected,NA),
         upr_Detected_new = ifelse(forecast_step==1,diff_upr_Detected+lag_Detected,NA)
  ) %>%
  mutate(fit_total_vol_new = ifelse(forecast_step==2,
                                    lag(fit_total_vol_new)+diff_fit_total_vol, fit_total_vol_new),
         lwr_total_vol_new = ifelse(forecast_step==2,
                                    lag(lwr_total_vol_new)+diff_lwr_total_vol, lwr_total_vol_new),
         upr_total_vol_new = ifelse(forecast_step==2,
                                    lag(upr_total_vol_new)+diff_upr_total_vol, upr_total_vol_new),
         
         fit_Detected_new = ifelse(forecast_step==2,
                                   lag(fit_Detected_new)+diff_fit_Detected, fit_Detected_new),
         lwr_Detected_new = ifelse(forecast_step==2,
                                   lag(lwr_Detected_new)+diff_lwr_Detected, lwr_Detected_new),
         upr_Detected_new = ifelse(forecast_step==2,
                                   lag(upr_Detected_new)+diff_upr_Detected, upr_Detected_new)
  )%>%
  mutate(fit_total_vol_new = ifelse(forecast_step==3,
                                    lag(fit_total_vol_new)+diff_fit_total_vol, fit_total_vol_new),
         lwr_total_vol_new = ifelse(forecast_step==3,
                                    lag(lwr_total_vol_new)+diff_lwr_total_vol, lwr_total_vol_new),
         upr_total_vol_new = ifelse(forecast_step==3,
                                    lag(upr_total_vol_new)+diff_upr_total_vol, upr_total_vol_new),
         
         fit_Detected_new = ifelse(forecast_step==3,
                                   lag(fit_Detected_new)+diff_fit_Detected, fit_Detected_new),
         lwr_Detected_new = ifelse(forecast_step==3,
                                   lag(lwr_Detected_new)+diff_lwr_Detected, lwr_Detected_new),
         upr_Detected_new = ifelse(forecast_step==3,
                                   lag(upr_Detected_new)+diff_upr_Detected, upr_Detected_new)
  )%>%
  mutate(fit_total_vol_new = ifelse(forecast_step==4,
                                    lag(fit_total_vol_new)+diff_fit_total_vol, fit_total_vol_new),
         lwr_total_vol_new = ifelse(forecast_step==4,
                                    lag(lwr_total_vol_new)+diff_lwr_total_vol, lwr_total_vol_new),
         upr_total_vol_new = ifelse(forecast_step==4,
                                    lag(upr_total_vol_new)+diff_upr_total_vol, upr_total_vol_new),
         
         fit_Detected_new = ifelse(forecast_step==4,
                                   lag(fit_Detected_new)+diff_fit_Detected, fit_Detected_new),
         lwr_Detected_new = ifelse(forecast_step==4,
                                   lag(lwr_Detected_new)+diff_lwr_Detected, lwr_Detected_new),
         upr_Detected_new = ifelse(forecast_step==4,
                                   lag(upr_Detected_new)+diff_upr_Detected, upr_Detected_new)
  )

simul_kf_multi_step_diff_bivar_informPrior_forecast_v4 <- replace(simul_kf_multi_step_diff_bivar_informPrior_forecast_v3,
                                                                  simul_kf_multi_step_diff_bivar_informPrior_forecast_v3 < 0, 0) %>%
  mutate(fit_total_vol = fit_total_vol_new,
         lwr_total_vol = lwr_total_vol_new,
         upr_total_vol = upr_total_vol_new,
         fit_Detected_vol = fit_Detected_new,
         lwr_Detected_vol = lwr_Detected_new,
         upr_Detected_vol = upr_Detected_new)

