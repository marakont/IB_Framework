## Intensity bioassay framework - example

pacman::p_load(
  rio,
  here,
  tidyverse,
  reshape2,
  rstan,
  loo,
  boot,
  egg,
  gridExtra,
  ggpubr
)

# EXAMPLE MODEL RUN (lab data) ----
df <- import(here("example_data.csv")) # load data
strains_sub <- unique(df$Strain) # create lookup for loop
maxdose <- max(df$Dose)

## Load stan model
options(mc.cores=4)
rstan::rstan_options(auto_write = TRUE)
model <- rstan::stan_model(here("example_model.stan"))

# I. Run one DR curve ---- 
df_single <- df %>% subset(Strain==strains_sub[5]) # subset data to model

## Generate model inputs
concentration_sim = seq(0,
                        sqrt(maxdose), # plot to max of df$Dose so all strains share same x axis
                        length.out = 500)
N_1 <- length(concentration_sim)
data_stan <- list(N=nrow(df_single),
                  mortality=as.integer(df_single$Responded),
                  tested=as.integer(df_single$Subjects),
                  concentration=sqrt(as.vector(df_single$Dose)),
                  N_1=as.integer(N_1),
                  concentration_sim=as.vector(concentration_sim))

## Fit the model
fit <- sampling(model, data=data_stan, iter=100, chains=4) # rhat < 1.01 for convergence

## 1/ Checking convergence ----
Rhat_df1 <- summary(fit)$summary[,10] # get Rhat values
modelcheck <- data.frame(not_converged=sum(Rhat_df1 > 1.01, na.rm = T), # check number of non-converged iterations
                         is_NA=ifelse(length(names(table(is.na(Rhat_df1))))==1,"no","yes"), # check if any iterations return NA
                         divergent_iterations=rstan::get_num_divergent(fit)) # check number of divergent iterations
### note: if not converged, need to run the model for more iterations 
###       (normal models normally ran for 5000-10000 iterations, but this takes more computing power and time)

## 2/ Extracting log likelihood & LOO for model comparison -----
logLikelihood_lab <- extract_log_lik(fit, "LogLikelihood")
LOO_lab <- loo(logLikelihood_lab)

## 3/ Generating simulated mortality curve ----
mortality <- rstan::extract(fit, "mean_mortality_sim")[[1]] # extract simulated mortality

### get mean simulated mortality for each concentration_sim
mean_mort <- apply(mortality, 2, mean) %>%
  as.data.frame() %>%
  mutate(concentration=concentration_sim^2) %>%
  rename(Test_mort_perc=1) %>%
  mutate(dat="sim")

### get lower CI of simulated mortality for each concentration_sim
mean_lower <- apply(mortality, 2, function(x) quantile(x, 0.025)) %>%
  as.data.frame() %>%
  mutate(concentration=concentration_sim^2) %>%
  rename(Test_mort_perc=1) %>%
  mutate(dat="sim")

### get upper CI of simulated mortality for each concentration_sim
mean_upper <- apply(mortality, 2, function(x) quantile(x, 0.975)) %>%
  as.data.frame() %>%
  mutate(concentration=concentration_sim^2) %>%
  rename(Test_mort_perc=1) %>%
  mutate(dat="sim")

### collate simulated mortality values together
mean_mort <- mean_mort %>% mutate(lower=mean_lower$Test_mort_perc,
                                  upper=mean_upper$Test_mort_perc)

### create tibble of model outputs
fit_mort <- tibble(Dose=mean_mort$concentration,
                   Mortality_perc=mean_mort$Test_mort_perc*100,
                   lower=mean_mort$lower*100, upper=mean_mort$upper*100,
                   dat=mean_mort$dat,
                   Strain=unique(df_single$Strain),
                   Insecticide=unique(df_single$Insecticide))

### Combine model outputs and real data
df1_s <- df_single %>% bind_rows(fit_mort)

### Plot DR curve
plot <- ggplot(df1_s, aes(x=Dose, y=Mortality_perc)) +
  geom_point(data = filter(df1_s, is.na(dat))) +                     # add real data points
  geom_ribbon(data = filter(df1_s, !is.na(dat)),
              aes(ymin=lower, ymax=upper), alpha=0.2) +              # add upper and lower bounds of simulated data
  geom_line(data = filter(df1_s, !is.na(dat)),
            aes(y=Mortality_perc)) +                                 # add line of simulated mortality
  theme_classic() +                                                  # theme customisation
  theme(panel.grid.major = element_line(colour = "grey93"),
        panel.grid.minor = element_line(colour = "grey97")) +        # theme customisation
  ylab("Mortality (%)") +                                            # rename y axis
  xlab("Insecticide concentration (%, ug/bottle, mg/m2)") +          # rename x axis
  ggtitle("Bioassay mortality (5-param logistic model, sqrt scale)", # Add title
          paste0(df_single$Insecticide,", Strain: ",df_single$Strain)) + # Add subtitle with insecticide and strain info
  scale_x_log10()                                                    # change x axis to log scale

plot

## 4/ Output different LCs ----
### Generate function to output simulated mortality at 'y' concentration
lcx <- function(y,B,C,E){
  return(exp(((log(((-1/(y-1))^(1/E))-1))/B)+C))
}

### Extract model parameter estimates
Bfit <- rstan::extract(fit)[["B"]]
Cfit <- rstan::extract(fit)[["C"]]
Efit <- rstan::extract(fit)[["E"]]

### Generate outputs for LC10, LC50, LC90 and LC99
sapply(c(10,50,90,99), function(a){
  
  message("Running LC summary ", a)
  
  temp <- do.call(rbind, sapply(1:length(Bfit), function(i){
    data.frame("LC" = lcx(a/100, Bfit[i],Cfit[i],Efit[i]))
  }, simplify = FALSE)) # Generate mortality estimate for each model iteration
  
  temp$LC_value <- a # Add column stating the LC value
  
  lc <- temp %>%
    mutate(LC=LC^2, # Transform LC value back to linear scale
           Strain=unique(df_single$Strain), # Add strain info
           Insecticide=unique(df_single$Insecticide)) # Add insecticide info
  
  lc_summ <- lc %>% # Generate summary df
    group_by(Strain,Insecticide) %>%
    summarise(LC_mean=mean(LC),
              LC_median=median(LC),
              LC_lower=quantile(LC,0.025,na.rm = T),
              LC_upper=quantile(LC,0.975,na.rm = T))
  
  ### Save outputs
  # write.csv(lc, file=here(paste0("outputs/lc", a, "_",unique(df_single$Strain),".csv")))
  # write.csv(lc_summ, file=herepaste0("outputs/lc_summ", a, "_", unique(df_single$Strain),".csv")))
 
  ### Outputs only for LC50
  if(a == 50){
    
    ### Add LC50 line to DR plot
    plot_LC <- plot + geom_segment(data=lc_summ,
                                   aes(x=LC_mean, y=-Inf,xend=LC_mean, yend=a),
                                   linetype="dashed")
    
    ### Plot density of LC50 fitted values
    ### (if tight, low uncertainty, if wide, high uncertainty)
    LC50_density_plot <- ggplot(lc, aes(x=LC)) +
      stat_density(position="identity",aes(alpha=0.9)) +
      ggtitle(paste0("Distribution of LC50 values: ",unique(df_single$Strain))) +
      theme_bw()  +
      guides(alpha=F)

    ### Save outputs
    # ggsave(plot_LC,file=here(paste0("outputs/DR-LC50_",unique(df_single$Strain),".png")),
    #        width = 4.85,height = 4.39)
    # ggsave(LC50_density_plot,file=here(paste0("outputs/LC50-density_",unique(df_single$Strain),".png")),
    #        width = 5,height = 4.39)
  }
  NULL
}, simplify = FALSE)

## 5/ Concentrations at which mosquitoes die ----

### Generate function of model logistic function 
### (outputs mortality for specified concentration x, given parameter estimates)
modelfun <- function(x,B,C,E){
  return(1 -  (1 /((1+exp(B*(log(x)-C)))^E)))
}

### Generate concentration vector within the bounds of the real input data
lower <- log10(1e-6)
upper <- log10(maxdose)
x <- seq(lower,upper, length.out=1000)
x <- 10^x
x <- sqrt(x)

### Inverse transform sampling to get the concentration at which mosquitoes die
### Generate the cumulative density function (CDF)
getCDF <- do.call(rbind, sapply(1:length(x), function(a){
  data.frame(sim_y = mean(sapply(1:length(Bfit), function(i){
    modelfun(x[a], Bfit[i],Cfit[i],Efit[i])
  })),
  sim_x=x[a])
}, simplify = F))

getCDF <- getCDF %>% 
  mutate(sim_x2=sim_x^2)

### Inverse CDF
mortx <- function(x,y){
  f <- approxfun(y,x) # interpolation
  u <- runif(20000) # random deviates from uniform distribution
  lcdplot = data.frame("Mort_C"=f(u))
}

getPDF <- data.frame(Mortx=mortx(getCDF$sim_x,getCDF$sim_y),
                  Strain=unique(df_single$Strain),
                  Insecticide=unique(df_single$Insecticide))

### Plot
mortx_plot <- ggplot(getPDF,aes(x=Mort_C)) +
  geom_density(alpha=.5, position="identity",fill="black") +
  labs(x=expression("Mortality at"~italic(x)~"concentration"), y="Density") +
  theme_classic() +
  ggtitle(paste0("Concentrations at which mosquitoes die (",unique(df_single$Strain),")")) +
  xlim(c(0,max(df_single$Dose))) + 
  scale_x_log10()

## 6/ Background mortality estimate ----
BM <- data.frame(mean_BM=mean(rstan::extract(fit)[["A"]]),
                median_BM=median(rstan::extract(fit)[["A"]]),
                lower_BM=quantile((rstan::extract(fit)[["A"]]),0.025),
                upper_BM=quantile((rstan::extract(fit)[["A"]]), 0.975),
                Strain=unique(df_single$Strain),
                Insecticide=unique(df_single$Insecticide))

## 7/ Model variability ----
## Mortality variability (along y axis)

### KEY:
# data_stan$concentration = actual x
# data_stan$mortality = actual y (NOTE: not percentages!!)
# fit, "mean_mortality" = simulated y

### STEPS:
### a) find all relevant x values (actual x values)
### b) extract actual y values for each of these x values
### c) generate simulated y value for each of these x values
### d) compute difference between actual y and simulated y (FOR EACH ITERATION)
### e) average these differences (FOR EACH ITERATION)
### f) compute mean, median and 95% CIs of variability estimate from all iterations

sim_y <- rstan::extract(fit,"mean_mortality_bm")[[1]] %>% 
  as.data.frame() # extract simulated mortality for each actual dose
act_x <- (data_stan$concentration)^2 # extract actual dose
act_y <- as.integer(df_single$Mortality_perc) # extract actual mortality
diff_it <- do.call(rbind,sapply(1:nrow(sim_y), function(i){
  abs((sim_y[i,]*100) - act_y) 
},simplify = F)) # compute absolute difference of actual vs simulated mortality (for each model iteration)
varmort_it <- apply(diff_it,1,sum)/length(act_x) # average difference across all concentrations (for each model iteration)

var_mort <- data.frame(Strain=unique(df_single$Strain), # summary table
                       Insecticide=unique(df_single$Insecticide),
                       mean_vmort=mean(varmort_it),
                       median_vmort=median(varmort_it),
                       lower_vmort=quantile(varmort_it,0.025),
                       upper_vmort=quantile(varmort_it,0.975))

## 8/ Model fitting assessment----
## 8.1/ Residuals plots ----
mean_mortality_bm <- rstan::extract(fit, "mean_mortality_bm")[[1]] # extract simulated mortality (incl A) for each actual dose

df1_residuals <- apply(mean_mortality_bm,2,median) %>% # get the median of all model iterations (for each actual dose)
  melt(value.name = "sim_y") %>% # re-format
  bind_cols(Concentration=data_stan$concentration^2,
            act_y=df_single$Mortality_perc/100) %>%
  mutate(Residuals=act_y-sim_y) # calculate residuals

resplot <- ggplot(df1_residuals,aes(x=Concentration, y=Residuals)) + # plot residuals
  geom_point() +
  geom_hline(aes(yintercept=0), linetype="dashed") +
  ggtitle(paste0("Residuals for ",unique(df_single$Strain))) +
  scale_x_log10()

# 8.2/ Actual vs simulated ----
mort_sim <- df_single %>%
  bind_cols(Mortality_sim_perc=apply(mean_mortality_bm,2,median)*100) # bind simulated mortality to actual dataframe

mod_lm <- mort_sim %>%
  group_by(Insecticide,Strain) %>%
  do(mod = lm(Mortality_perc~Mortality_sim_perc, data=.)) # Linear regression of actual vs simulated mortality
df_coeff <- mod_lm %>%
  do(data.frame(
    Strain=.$Strain,
    Insecticide=.$Insecticide,
    var=names(coef(.$mod)), # get names of coefficients
    coef(summary(.$mod)), # get summary outputs of model
    r2=summary(.$mod)$r.square, # get model R2
    RMSE=sqrt(mean(.$mod$residuals^2)))) # get root mean squared error

### plot actual vs simulated mortality, with R2 & RMSE estimates
ap_all <- ggplot(mort_sim,aes(x=Mortality_sim_perc, y=Mortality_perc)) +
  geom_point(size=4) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  ylab("Actual mortality (%)") +
  xlab("Predicted mortality (%)") +
  theme_bw() +
  ggtitle(paste0("Field model predictive accuracy (",unique(df_single$Strain),")")) +
  geom_text(data=df_coeff,
            aes(x=15,y=75,
                label=paste0('R2=',round(r2,digits=2),"\nRMSE=",round(RMSE,digits=0),"%"),
                fontface=3),
            size=5) +
  xlim(c(0,100)) + ylim(c(0,100))

## 9/ Save objects----
# saveRDS(fit, file=here(paste0("outputs/fit_",unique(df_single$Strain),".rds")))
# write.csv(modelcheck, file=here(paste0("outputs/modelcheck_",unique(df_single$Strain),".csv")))
# saveRDS(LOO_lab, file=here(paste0("outputs/LOO_",unique(df_single$Strain),".rds")))
# write.csv(df1_s, file=here(paste0("outputs/df1_s_",unique(df_single$Strain),".csv")))
# ggsave(plot,file=here(paste0("outputs/DR_",unique(df_single$Strain),".png")),
#        width = 4.85,height = 4.39)
# write.csv(getPDF, file=here(paste0("outputs/PDF_",unique(df_single$Strain),".csv")))
# ggsave(mortx_plot, file=here(paste0("outputs/mortxplot_",unique(df_single$Strain),".png")),
#        width = 4.85,height = 4.39)
# write.csv(BM, file=here(paste0("outputs/BM_",unique(df_single$Strain),".csv")))
# write.csv(var_mort, file=here(paste0("outputs/varmort_",unique(df_single$Strain),".csv")))
# write.csv(df1_residuals,file=here(paste0("outputs/residuals_",unique(df_single$Strain),".csv")))
# ggsave(resplot,file=here(paste0("outputs/residuals_",unique(df_single$Strain),".png")),
#        width = 9.5,height = 4.39)
# write.csv(mort_sim,file=here(paste0("outputs/mortsim_",unique(df_single$Strain),".csv")))
# write.csv(df_coeff,file=here(paste0("outputs/AP-coeff_",unique(df_single$Strain),".csv")))
# ggsave(ap_all,file=here(paste0("outputs/AP_",unique(df_single$Strain),".png")),
#        width = 4.85,height = 4.39)

# II. Run model for multiple DR curves ----
# Model code
bioassay_lab_sqrt <- function(df_single) {
  concentration_sim = seq(0,
                          sqrt(maxdose), # plot to max of df$Dose so all strains share same x axis
                          length.out = 500)
  N_1 <- length(concentration_sim)
  data_stan <- list(N=nrow(df_single),
                    mortality=as.integer(df_single$Responded),
                    tested=as.integer(df_single$Subjects),
                    concentration=sqrt(as.vector(df_single$Dose)),
                    N_1=as.integer(N_1),
                    concentration_sim=as.vector(concentration_sim))

  ## Fit the model
  fit <- sampling(model, data=data_stan, iter=5000, chains=4) # rhat < 1.01 for convergence
  
  ## 1/ Checking convergence ----
  Rhat_df1 <- summary(fit)$summary[,10] # get Rhat values
  modelcheck <- data.frame(not_converged=sum(Rhat_df1 > 1.01, na.rm = T), # check number of non-converged iterations
                           is_NA=ifelse(length(names(table(is.na(Rhat_df1))))==1,"no","yes"), # check if any iterations return NA
                           divergent_iterations=rstan::get_num_divergent(fit)) # check number of divergent iterations
  ### note: if not converged, need to run the model for more iterations 
  ###       (normal models normally ran for 5000-10000 iterations, but this takes more computing power and time)
  
  ## 2/ Extracting log likelihood & LOO for model comparison -----
  logLikelihood_lab <- extract_log_lik(fit, "LogLikelihood")
  LOO_lab <- loo(logLikelihood_lab)
  
  ## 3/ Generating simulated mortality curve ----
  mortality <- rstan::extract(fit, "mean_mortality_sim")[[1]] # extract simulated mortality
  
  ### get mean simulated mortality for each concentration_sim
  mean_mort <- apply(mortality, 2, mean) %>%
    as.data.frame() %>%
    mutate(concentration=concentration_sim^2) %>%
    rename(Test_mort_perc=1) %>%
    mutate(dat="sim")
  
  ### get lower CI of simulated mortality for each concentration_sim
  mean_lower <- apply(mortality, 2, function(x) quantile(x, 0.025)) %>%
    as.data.frame() %>%
    mutate(concentration=concentration_sim^2) %>%
    rename(Test_mort_perc=1) %>%
    mutate(dat="sim")
  
  ### get upper CI of simulated mortality for each concentration_sim
  mean_upper <- apply(mortality, 2, function(x) quantile(x, 0.975)) %>%
    as.data.frame() %>%
    mutate(concentration=concentration_sim^2) %>%
    rename(Test_mort_perc=1) %>%
    mutate(dat="sim")
  
  ### collate simulated mortality values together
  mean_mort <- mean_mort %>% mutate(lower=mean_lower$Test_mort_perc,
                                    upper=mean_upper$Test_mort_perc)
  
  ### create tibble of model outputs
  fit_mort <- tibble(Dose=mean_mort$concentration,
                     Mortality_perc=mean_mort$Test_mort_perc*100,
                     lower=mean_mort$lower*100, upper=mean_mort$upper*100,
                     dat=mean_mort$dat,
                     Strain=unique(df_single$Strain),
                     Insecticide=unique(df_single$Insecticide))
  
  ### Combine model outputs and real data
  df1_s <- df_single %>% bind_rows(fit_mort)
  
  ### Plot DR curve
  plot <- ggplot(df1_s, aes(x=Dose, y=Mortality_perc)) +
    geom_point(data = filter(df1_s, is.na(dat))) +                     # add real data points
    geom_ribbon(data = filter(df1_s, !is.na(dat)),
                aes(ymin=lower, ymax=upper), alpha=0.2) +              # add upper and lower bounds of simulated data
    geom_line(data = filter(df1_s, !is.na(dat)),
              aes(y=Mortality_perc)) +                                 # add line of simulated mortality
    theme_classic() +                                                  # theme customisation
    theme(panel.grid.major = element_line(colour = "grey93"),
          panel.grid.minor = element_line(colour = "grey97")) +        # theme customisation
    ylab("Mortality (%)") +                                            # rename y axis
    xlab("Insecticide concentration (%, ug/bottle, mg/m2)") +          # rename x axis
    ggtitle("Bioassay mortality (5-param logistic model, sqrt scale)", # Add title
            paste0(df_single$Insecticide,", Strain: ",df_single$Strain)) + # Add subtitle with insecticide and strain info
    scale_x_log10()                                                    # change x axis to log scale
  
  plot
  
  ## 4/ Output different LCs ----
  ### Generate function to output simulated mortality at 'y' concentration
  lcx <- function(y,B,C,E){
    return(exp(((log(((-1/(y-1))^(1/E))-1))/B)+C))
  }
  
  ### Extract model parameter estimates
  Bfit <- rstan::extract(fit)[["B"]]
  Cfit <- rstan::extract(fit)[["C"]]
  Efit <- rstan::extract(fit)[["E"]]
  
  ### Generate outputs for LC10, LC50, LC90 and LC99
  sapply(c(10,50,90,99), function(a){
    
    message("Running LC summary ", a)
    
    temp <- do.call(rbind, sapply(1:length(Bfit), function(i){
      data.frame("LC" = lcx(a/100, Bfit[i],Cfit[i],Efit[i]))
    }, simplify = FALSE)) # Generate mortality estimate for each model iteration
    
    temp$LC_value <- a # Add column stating the LC value
    
    lc <- temp %>%
      mutate(LC=LC^2, # Transform LC value back to linear scale
             Strain=unique(df_single$Strain), # Add strain info
             Insecticide=unique(df_single$Insecticide)) # Add insecticide info
    
    lc_summ <- lc %>% # Generate summary df
      group_by(Strain,Insecticide) %>%
      summarise(LC_mean=mean(LC),
                LC_median=median(LC),
                LC_lower=quantile(LC,0.025,na.rm = T),
                LC_upper=quantile(LC,0.975,na.rm = T))
    
    ### Save outputs
    write.csv(lc, file=here(paste0("outputs/lc", a, "_",unique(df_single$Strain),".csv")))
    write.csv(lc_summ, file=here(paste0("outputs/summ_lc", a, "_", unique(df_single$Strain),".csv")))
    
    ### Outputs only for LC50
    if(a == 50){
      
      ### Add LC50 line to DR plot
      plot_LC <- plot + geom_segment(data=lc_summ,
                                     aes(x=LC_mean, y=-Inf,xend=LC_mean, yend=a),
                                     linetype="dashed")
      
      ### Plot density of LC50 fitted values
      ### (if tight, low uncertainty, if wide, high uncertainty)
      LC50_density_plot <- ggplot(lc, aes(x=LC)) +
        stat_density(position="identity",aes(alpha=0.9)) +
        ggtitle(paste0("Distribution of LC50 values: ",unique(df_single$Strain))) +
        theme_bw()  +
        guides(alpha=F)
      
      ### Save outputs
      ggsave(plot_LC,file=here(paste0("outputs/DR-LC50_",unique(df_single$Strain),".png")),
             width = 4.85,height = 4.39)
      ggsave(LC50_density_plot,file=here(paste0("outputs/LC50-density_",unique(df_single$Strain),".png")),
             width = 5,height = 4.39)
    }
    NULL
  }, simplify = FALSE)
  
  ## 5/ Concentrations at which mosquitoes die ----
  
  ### Generate function of model logistic function 
  ### (outputs mortality for specified concentration x, given parameter estimates)
  modelfun <- function(x,B,C,E){
    return(1 -  (1 /((1+exp(B*(log(x)-C)))^E)))
  }
  
  ### Generate concentration vector within the bounds of the real input data
  lower <- log10(1e-6)
  upper <- log10(maxdose)
  x <- seq(lower,upper, length.out=1000)
  x <- 10^x
  x <- sqrt(x)
  
  ### Inverse transform sampling to get the concentration at which mosquitoes die
  ### Generate the cumulative density function (CDF)
  getCDF <- do.call(rbind, sapply(1:length(x), function(a){
    data.frame(sim_y = mean(sapply(1:length(Bfit), function(i){
      modelfun(x[a], Bfit[i],Cfit[i],Efit[i])
    })),
    sim_x=x[a])
  }, simplify = F))
  
  getCDF <- getCDF %>% 
    mutate(sim_x2=sim_x^2)
  
  ### Inverse CDF
  mortx <- function(x,y){
    f <- approxfun(y,x) # interpolation
    u <- runif(20000) # random deviates from uniform distribution
    lcdplot = data.frame("Mort_C"=f(u))
  }
  
  getPDF <- data.frame(Mortx=mortx(getCDF$sim_x,getCDF$sim_y),
                    Strain=unique(df_single$Strain),
                    Insecticide=unique(df_single$Insecticide))
  
  ### Plot
  mortx_plot <- ggplot(getPDF,aes(x=Mort_C)) +
    geom_density(alpha=.5, position="identity",fill="black") +
    labs(x=expression("Mortality at"~italic(x)~"concentration"), y="Density") +
    theme_classic() +
    ggtitle(paste0("Concentrations at which mosquitoes die (",unique(df_single$Strain),")")) + 
    scale_x_log10()
  
  
  ## 6/ Background mortality estimate ----
  BM <- data.frame(mean_BM=mean(rstan::extract(fit)[["A"]]),
                  median_BM=median(rstan::extract(fit)[["A"]]),
                  lower_BM=quantile((rstan::extract(fit)[["A"]]),0.025),
                  upper_BM=quantile((rstan::extract(fit)[["A"]]), 0.975),
                  Strain=unique(df_single$Strain),
                  Insecticide=unique(df_single$Insecticide))

  ## 7/ Model variability ----
  ## Mortality variability (along y axis)
  
  ### KEY:
  # data_stan$concentration = actual x
  # data_stan$mortality = actual y (NOTE: not percentages!!)
  # fit, "mean_mortality" = simulated y
  
  ### STEPS:
  ### a) find all relevant x values (actual x values)
  ### b) extract actual y values for each of these x values
  ### c) generate simulated y value for each of these x values
  ### d) compute difference between actual y and simulated y (FOR EACH ITERATION)
  ### e) average these differences (FOR EACH ITERATION)
  ### f) compute mean, median and 95% CIs of variability estimate from all iterations
  
  sim_y <- rstan::extract(fit,"mean_mortality")[[1]] %>% 
    as.data.frame() # extract simulated mortality for each actual dose
  act_x <- (data_stan$concentration)^2 # extract actual dose
  act_y <- as.integer(df_single$Mortality_perc) # extract actual mortality
  diff_it <- do.call(rbind,sapply(1:nrow(sim_y), function(i){
    abs((sim_y[i,]*100) - act_y) 
  },simplify = F)) # compute absolute difference of actual vs simulated mortality (for each model iteration)
  varmort_it <- apply(diff_it,1,sum)/length(act_x) # average difference across all concentrations (for each model iteration)
  
  var_mort <- data.frame(Strain=unique(df_single$Strain), # summary table
                         Insecticide=unique(df_single$Insecticide),
                         mean_vmort=mean(varmort_it),
                         median_vmort=median(varmort_it),
                         lower_vmort=quantile(varmort_it,0.025),
                         upper_vmort=quantile(varmort_it,0.975))

  ## 8/ Model fitting assessment----
  ## 8.1/ Residuals plots ----
  mean_mortality_bm <- rstan::extract(fit, "mean_mortality_bm")[[1]] # extract simulated mortality (incl A) for each actual dose
  
  df1_residuals <- apply(mean_mortality_bm,2,median) %>% # get the median of all model iterations (for each actual dose)
    melt(value.name = "sim_y") %>% # re-format
    bind_cols(Concentration=data_stan$concentration^2,
              act_y=df_single$Mortality_perc/100) %>%
    mutate(Residuals=act_y-sim_y) # calculate residuals
  
  resplot <- ggplot(df1_residuals,aes(x=Concentration, y=Residuals)) + # plot residuals
    geom_point() +
    geom_hline(aes(yintercept=0), linetype="dashed") +
    ggtitle(paste0("Residuals for ",unique(df_single$Strain))) + 
    scale_x_log10()

  # 8.2/ Actual vs simulated ----
  mort_sim <- df_single %>%
    bind_cols(Mortality_sim_perc=apply(mean_mortality_bm,2,median)*100) # bind simulated mortality to actual dataframe
  
  mod_lm <- mort_sim %>%
    group_by(Insecticide,Strain) %>%
    do(mod = lm(Mortality_perc~Mortality_sim_perc, data=.)) # Linear regression of actual vs simulated mortality
  df_coeff <- mod_lm %>%
    do(data.frame(
      Strain=.$Strain,
      Insecticide=.$Insecticide,
      var=names(coef(.$mod)), # get names of coefficients
      coef(summary(.$mod)), # get summary outputs of model
      r2=summary(.$mod)$r.square, # get model R2
      RMSE=sqrt(mean(.$mod$residuals^2)))) # get root mean squared error

  ### plot actual vs simulated mortality, with R2 & RMSE estimates
  ap_all <- ggplot(mort_sim,aes(x=Mortality_sim_perc, y=Mortality_perc)) +
    geom_point(size=4) +
    geom_abline(intercept = 0, slope = 1, linetype="dashed") +
    ylab("Actual mortality (%)") +
    xlab("Predicted mortality (%)") +
    theme_bw() +
    ggtitle(paste0("Field model predictive accuracy (",unique(df_single$Strain),")")) +
    geom_text(data=df_coeff,
              aes(x=15,y=75,
                  label=paste0('R2=',round(r2,digits=2),"\nRMSE=",round(RMSE,digits=0),"%"),
                  fontface=3),
              size=5) +
    xlim(c(0,100)) + ylim(c(0,100))

  ## 9/ Save objects----
  saveRDS(fit, file=here(paste0("outputs/fit_",unique(df_single$Strain),".rds")))
  write.csv(modelcheck, file=here(paste0("outputs/modelcheck_",unique(df_single$Strain),".csv")))
  saveRDS(LOO_lab, file=here(paste0("outputs/LOO_",unique(df_single$Strain),".rds")))
  write.csv(df1_s, file=here(paste0("outputs/df1_s_",unique(df_single$Strain),".csv")))
  ggsave(plot,file=here(paste0("outputs/DR_",unique(df_single$Strain),".png")),
         width = 4.85,height = 4.39)
  write.csv(getPDF, file=here(paste0("outputs/PDF_",unique(df_single$Strain),".csv")))
  ggsave(mortx_plot, file=here(paste0("outputs/mortxplot_",unique(df_single$Strain),".png")),
         width = 4.85,height = 4.39)
  write.csv(BM, file=here(paste0("outputs/BM_",unique(df_single$Strain),".csv")))
  write.csv(var_mort, file=here(paste0("outputs/varmort_",unique(df_single$Strain),".csv")))
  write.csv(df1_residuals,file=here(paste0("outputs/residuals_",unique(df_single$Strain),".csv")))
  ggsave(resplot,file=here(paste0("outputs/residuals_",unique(df_single$Strain),".png")),
         width = 9.5,height = 4.39)
  write.csv(mort_sim,file=here(paste0("outputs/mortsim_",unique(df_single$Strain),".csv")))
  write.csv(df_coeff,file=here(paste0("outputs/AP-coeff_",unique(df_single$Strain),".csv")))
  ggsave(ap_all,file=here(paste0("outputs/AP_",unique(df_single$Strain),".png")),
         width = 4.85,height = 4.39)
  
}

## 10/ Run loop ----  
ldf <- sapply(seq_along(strains_sub),
              function(x) {
                df %>% subset(Strain==strains_sub[x])
              }, simplify = F) # create list of each individually bioassay
routput_lab <- lapply(ldf, bioassay_lab_sqrt) # model each bioassay

## 11/ Convergence check ----
### Bind all strains together
conv_ll <- do.call(rbind,sapply(seq_along(strains_sub),
                                function(i) {
                                  read.csv(here(paste0("outputs/modelcheck_",strains_sub[i],".csv"))) %>% 
                                    mutate(strain=strains_sub[i])                      
                   },simplify=F))

write.csv(conv_ll, here("outputs/summary/modelcheck_all.csv"))

# III. SUMMARY PLOTS ----
## 1/ DR plots ----
### Bind all strains together
sdf1_s <- do.call(rbind,sapply(seq_along(strains_sub),
                                function(i) {
                                  read.csv(here(paste0("outputs/df1_s_",strains_sub[i],".csv")))
                                },simplify=F))

write.csv(sdf1_s, here("outputs/summary/df_all_sim.csv"))

### Clean
sdf1_s <- sdf1_s %>% 
  rename("Concentration"=Dose) %>% 
  subset(Concentration<=0.09) # Remove last concentration in Tiassale-13 assay so can compare different strains

### Plot
ggplot(sdf1_s, aes(x=Concentration, y=Mortality_perc)) +
  geom_point(data=filter(sdf1_s, is.na(dat)),
             aes(colour=Strain)) +
  geom_line(data=filter(sdf1_s, !is.na(dat)),
            aes(colour=Strain)) +
  geom_ribbon(data=filter(sdf1_s, !is.na(dat)),
              aes(ymin=lower, ymax=upper,fill=Strain), alpha=.2) +
  facet_wrap(~Strain, scales = "free_x") +
  ylab("Mortality (%)")  + ggtitle("Bioassay mortality","5PL model (sqrt), LITE dataset") +
  theme_classic() +
  theme(panel.grid.major = element_line(colour = "grey93"),
        panel.grid.minor = element_line(colour = "grey97"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_x_sqrt() +
  ggsci::scale_fill_npg() +
  ggsci::scale_colour_npg()

## Save
ggsave(here("outputs/summary/DR_all_strains_separated.png"),
       width = 10, height = 6.5, units = "in",scale=0.8)

### Plot all strains on one plot
p1 <- ggplot(sdf1_s, aes(x=Concentration, y=Mortality_perc)) +
  geom_point(data=filter(sdf1_s, is.na(dat)),
             aes(colour=Strain)) +
  geom_line(data=filter(sdf1_s, !is.na(dat)),
            aes(colour=Strain)) +
  geom_ribbon(data=filter(sdf1_s, !is.na(dat)),
              aes(ymin=lower, ymax=upper,fill=Strain), alpha=.2) +
  ylab("Mortality (%)")  + ggtitle("Bioassay mortality","5PL model (sqrt), LITE dataset") +
  xlab("Concentration (%)") +
  theme_classic() +
  theme(panel.grid.major = element_line(colour = "grey93"),
        panel.grid.minor = element_line(colour = "grey97"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_x_sqrt() +
  ggsci::scale_fill_npg() +
  ggsci::scale_colour_npg()

p1

### Save
ggsave(here("outputs/summary/DR_all_strains.png"),p1,
       width = 7, height = 5.52, units = "in",scale=0.6)

## 2/ Concentrations at which mosquitoes die ----
### Bind all strains together
mortx <- do.call(rbind,sapply(seq_along(strains_sub),
                              function(i) {
                                read.csv(file=here(paste0("outputs/PDF_",strains_sub[i],".csv"))) %>% 
                                  mutate(Strain=strains_sub[i])
                              },simplify=F))

### Plot
p2 <- ggplot(mortx,aes(x=Mort_C, fill=Strain, colour=Strain)) +
  geom_density(alpha=.7, position="identity") +
  labs(x=expression("Concentration (%)"), y="Mortality density") +
  theme_classic() +
  ggtitle("Concentrations at which mosquitoes die") +
  ggsci::scale_fill_npg() +
  ggsci::scale_colour_npg() +
  scale_x_sqrt(limits=c(0,0.4),breaks=c(0.025,0.05,0.075,0.1,0.2,0.3,0.4))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

p2

### Save
ggsave(here("outputs/summary/mortx_all.png"),p2,
       width = 8, height = 6.5, units = "in",scale=0.6)

### 3/ Mortality variability ----
### Bind all strains together
var_mort_all <- do.call(rbind, sapply(seq_along(strains_sub), function(i){
  read.csv(here(paste0("outputs/varmort_",strains_sub[i],".csv")))
},simplify = F))

### Plot
p3 <- ggplot(var_mort_all,aes(y=Strain,x=median_vmort,colour=Strain))+
  geom_pointrange(aes(xmin=lower_vmort,
                      xmax=upper_vmort),
                  size=1,shape=15,
                  position=position_dodge(width=0.5)) +
  theme_minimal() +
  xlab("Median variability estimate\n(% variability in mortality from best fit line)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggsci::scale_colour_npg()  +
  theme(axis.title.y = element_blank())

p3

### Save
ggsave(here("outputs/summary/mortvar_all.png"),p3,
       width = 4.1, height = 4, units = "in")

## COMBINE PLOTS 1-3 ----
ggpubr::ggarrange(p1,p2,p3,nrow=1,
                  common.legend = T, legend = "bottom",
                  labels = c("A","B","C"))
ggsave(here("outputs/summary/lab_run.png"),
       width=13.2, height=4.5, units="in",scale=.8)

## 4/ AP plots ----
### Bind all strains together
mortsim <- do.call(rbind,sapply(seq_along(strains_sub),
                                function(i) {
                                  read_csv(file=here(paste0("outputs/mortsim_",strains_sub[i],".csv")))
                                },simplify=F))

### Plot
p4 <- ggplot(mortsim,aes(x=Mortality_sim_perc, y=Mortality_perc)) +
  geom_point(size=3, aes(colour=Strain)) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  ylab("Actual mortality (%)") +
  xlab("Predicted mortality (%)") +
  theme_bw() +
  ggsci::scale_colour_npg()

p4

### Save
ggsave("outputs/summary/residuals.png",p4,
       width = 6.5, height = 4.5, units = "in")

## AP table
### Bind all strains together
df_coeff <- do.call(rbind,sapply(seq_along(strains_sub),
                                 function(i) {
                                   read_csv(file=here(paste0("outputs/AP-coeff_",strains_sub[i],".csv")))
                                 },simplify=F))

### Create table
t <- df_coeff %>% group_by(Strain) %>% 
  subset(var%in%unique(df_coeff$var)[2]) %>% 
  ungroup() %>%
  dplyr::select(Strain,r2,RMSE) %>% 
  mutate(r2=round(r2,2),
         RMSE=paste0(as.character(round(RMSE,0)),"%")) %>% 
  rename("R2"=r2) %>%  
  arrange(Strain)

### Arrange table
t1 <- ggtexttable(t, rows=NULL, theme=ttheme("classic"))
t1 <- table_cell_bg(t1,column = 1,row=2,
                 fill = "#E64B35FF")
t1 <- table_cell_bg(t1,column = 1,row=3,
                 fill = "#4DBBD5FF")
t1 <- table_cell_bg(t1,column = 1,row=4,
                 fill = "#00A087FF")
t1 <- table_cell_bg(t1,column = 1,row=5,
                 fill = "#3C5488FF")
t1 <- table_cell_bg(t1,column = 1,row=6,
                 fill = "#F39B7FFF")

t1

## 5/ Residuals plots ----
### Bind all strains together
res <- do.call(rbind,sapply(seq_along(strains_sub),
                            function(i) {
                              read_csv(file=here(paste0("outputs/residuals_",strains_sub[i],".csv"))) %>% 
                                mutate(Strain=strains_sub[i])
                            },simplify=F))

### Plot
p5 <- ggplot(res,aes(x=Concentration, y=Residuals)) +
  geom_point(aes(colour=Strain),size=2) +
  geom_hline(aes(yintercept=0), linetype="dashed") +
  xlab("Concentration (%)") +
  ylim(c(-.2,.2)) +
  scale_x_log10() +
  theme_bw() +
  ggsci::scale_colour_npg()

p5

### Save
ggsave(here("outputs/summary/residuals_all.png"),p5,
       width = 10, height = 6.5, units = "in")

## COMBINE PLOTS 4-5 FOR MODEL ASSESSMENT ----
ggpubr::ggarrange(p4,t1,NULL,p5,
                  common.legend = T,legend = "bottom",
                  nrow=1,
                  labels = c("A","","","B"),
                  widths = c(1,0.7,0.1,1),
                  heights = c(1,0.7,1,1))
ggsave(here("outputs/summary/lab_model_assessment.png"),
       width = 9.72, height = 3.91, units = "in",scale = 0.9)

# 6/ Summary table ----
## 6.1/ LC50 ----
### Bind all strains together
LC50 <- do.call(rbind, sapply(seq_along(strains_sub),
                              function(i){
                                read.csv(here(paste0("outputs/summ_lc50_", strains_sub[i],".csv")))
                              },simplify = F))

### Function to format numbers
number_formatting <- function(x){
  ifelse(0.001>x,
         formatC(x,format = "e",digits = 1),
         ifelse(x>10 & x<=100, round(x,1),
                ifelse(x>100,
                       formatC(x,format = "e",digits = 1),
                       round(x,3))))
}

### LC50 median and CI for all strains
LC50 <- LC50 %>% 
  mutate(LC_lower=as.character(number_formatting(LC_lower)),
         LC_upper=as.character(number_formatting(LC_upper))) %>% 
  unite("LC50 (95% CI)",c(LC_lower,LC_upper), sep = " - ") %>% 
  rename("LC50 (median)"=LC_median) %>% 
  dplyr::select(-LC_mean)

## 6.2/ Heterogeneity (LC10 and LC90) ----
LC10 <- do.call(rbind, sapply(seq_along(strains_sub),
                              function(i){
                                read.csv(paste0("outputs/summ_lc10_", strains_sub[i],".csv"))
                              },simplify = F)) %>% 
  dplyr::select(c(2,3,5)) %>% 
  rename("LC10"=LC_median)
LC90 <- do.call(rbind, sapply(seq_along(strains_sub),
                              function(i){
                                read.csv(paste0("outputs/summ_lc90_", strains_sub[i],".csv"))
                              },simplify = F)) %>% 
  dplyr::select(c(2,3,5)) %>% 
  rename("LC90"=LC_median)

## 6.3/ Mortality variability ----
var_mort <- do.call(rbind,sapply(seq_along(strains_sub),
                                 function(i){
                                   read.csv(file=here(paste0("outputs/varmort_",strains_sub[i],".csv")))
                                 },simplify = F))
var_mort <- var_mort %>% 
  mutate(lower_vmort=as.character(number_formatting(lower_vmort)),
         upper_vmort=as.character(number_formatting(upper_vmort))) %>% 
  unite("variability (95% CIs)",c(lower_vmort,upper_vmort), sep=" - ") %>% 
  rename("variability (median)"=median_vmort) %>% 
  dplyr::select(-c(X,mean_vmort))

## 6.4/ BM ----
BM <- do.call(rbind, sapply(seq_along(strains_sub),
                            function(i){
                              read.csv(here(paste0("outputs/BM_", strains_sub[i],".csv")))
                            },simplify = F))
BM[,c(2:5)] <- BM[,c(2:5)]*100
for(i in 2:5){
  BM[,i] <- round(BM[,i],2)
  BM[,i] <- paste0(BM[,i],"%")
}
BM <- BM %>% 
  unite("BM (95% CIs)",c(lower_BM,upper_BM), sep=" - ") %>% 
  rename("BM (mean)"=mean_BM) %>% 
  dplyr::select(-c(X,median_BM))

## 6.5/ COMBINE INTO TABLE ----
stable <- df %>% group_by(Insecticide,Strain) %>% summarise("Number of data points" = n(),
                                                            "Number of mosquitoes"=sum(Subjects)) %>% 
  left_join(LC50) %>% 
  left_join(LC10) %>% 
  left_join(LC90) %>% 
  left_join(BM) %>% 
  left_join(var_mort) %>% 
  dplyr::select(-X)

write.csv(stable,here("outputs/summary/lab_summ_table.csv"))
