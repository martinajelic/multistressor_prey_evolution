# ----------------------------------------------------------------------------------------------------------------------------------
# Search grids of logistic prey growth parameters for prey from each selection environment
# ----------------------------------------------------------------------------------------------------------------------------------
#
# The following R-code is based on the method of:
#
# Rosenbaum, B. & Rall, B.C. (2018). Fitting functional responses: Direct parameter estimation by
# simulating differential equations. Methods in Ecology and Evolution.
# 
# Lemmen, K. D. & Pennekamp, F. (2024). Food web context modifies predator foraging and weakens 
# trophic interaction strength
#
# Daugaard, U., Petchey, O.L. & Pennekamp, F. (2019). Warming can destabilize predator–prey 
# interactions by shifting the functional response from Type III to Type II.
# 
# This method estimates the logistic growth parameters directly by simulating differential equations
#
# For this, the package odeintr is used. For a manual about the use of the method and a detailed
# description of the functions, see the paper of Rosenbaum and Rall (2018) mentioned above.
# 
# NOTE: the package odeintr is available on github for R version 4.3.2 and is unsupported for newer 
# versions of R at the moment, install with the following code:
# devtools::install_github("thk686/odeintr")
#
# Rdata files with the output from the search grid process can be found in the folder
# "Data/FR assay/02 - LG search grid" and can be manually loaded into the workspace.
#
# Running the search grid (even with the same starting parameter values) may result in final 
# parameter estimates that differ slightly from those reported in the manuscript. This is expected 
# due to the optimization procedure.
#
# ----------------------------------------------------------------------------------------------------------------------------------

rm(list = ls())

library(doParallel)
library(tidyverse)
library(bbmle)

source("R scripts/Function_logistic_growth.R")

##### Grid Search for r and K estimates #####
pre_sg <- expand.grid(r = c(-2, -1, 0, 0.5, 1, 1.5, 2, 2.5),
                      K.log = log(c(1000, 1500, 2000, 2500, 3000, 3500)),
                      sigma = 0.2)

pre_sg$error <- NA

pre_sg$r.est <- NA
pre_sg$K.log.est <- NA
pre_sg$sigma.est <- NA

pre_sg$converge <- NA
pre_sg$iteration <- NA
pre_sg$AICc <- NA
pre_sg$logLik <- NA
pre_sg$strtpt <- seq(1:nrow(pre_sg))




#---------------------------------------------------------------------------------------------------
# Prey from selection environment C
#---------------------------------------------------------------------------------------------------
load("Data/FR assay/01 - Raw data/Dataset_C.RData")

growth_sim_grid <- pre_sg

###### Set the computer to run in Parallel ######
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores() - 1,
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered() # returns TRUE if a doPar backend has been registered
foreach::getDoParWorkers()    # nr. of execution workers in the currently registered doPar backend

###### Run Grid ######
preyC_growth_pe <- NA # name of file where parallel runs will be saved

(t0 <- Sys.time()) # timer
mcoptions <- list(preschedule=FALSE)

preyC_growth_pe <- foreach(i = 1:nrow(growth_sim_grid),
                          .combine = 'rbind',
                          .options.multicore = mcoptions,
                          .packages = c('bbmle', 'odeintr', 'tidyverse', 'R.utils'),
                          .errorhandling = "pass") %dopar%
  {
    compile_sys("LG_prey", LG.prey, pars = c("r", "K"), method = "rk54")
    
    r<-tryCatch({
      setTimeLimit(35*60)
      fit = mle2(minuslogl = nll.odeint.LG.prey,
                 method = "BFGS",
                 control = list(maxit = 5000),
                 start = list(r = growth_sim_grid$r[[i]],
                              K.log = growth_sim_grid$K.log[[i]],
                              sigma = growth_sim_grid$sigma[[i]]),
                 data = with(list(N0 = Prey_start_density,
                                  Ndiff = Prey_end_density - Prey_start_density,
                                  Tt = Incubation_time),
                             data = dplyr::filter(FR_Dataset_C, Predator_presence == 0)))
      
      
      if(fit@details$convergence==0) {
        growth_sim_grid$r.est[i] <- fit@coef[[1]]
        growth_sim_grid$K.log.est[i] <- fit@coef[[2]]
        
        growth_sim_grid$sigma.est [i] <- fit@coef[[3]]
        
        growth_sim_grid$iteration[i] <- fit@details[[3]][[2]]
        growth_sim_grid$converge[i] <- fit@details[[4]]
        growth_sim_grid$AICc[i] <- AICc(fit)
        growth_sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        growth_sim_grid$r.est[i] <- NaN
        growth_sim_grid$K.log.est[i] <- NaN
        
        growth_sim_grid$converge[i] <- NaN
        growth_sim_grid$iteration[i] <- NaN
        growth_sim_grid$AICc[i] <- NaN
        growth_sim_grid$logLik[i] <- NaN
      }
      
      #Remove fit for next iteration
      rm(fit)
    }, error=function(e) list(message = e$message)) # to print error in console use {cat("ERROR :", paste(conditionMessage(e), "Combo", (i)), "\n")}
    
    if(length(r) > 0){
      growth_sim_grid$error[i] <- r
    } else {
      growth_sim_grid$error[i] <- NaN
    }
    
    return(growth_sim_grid[i,])
    
  }

parallel::stopCluster(cl = my.cluster)
stopImplicitCluster()

tend <- Sys.time()
tend - t0
# beep(sound = 2, expr = NULL)

preyC_growth_pe %>% arrange(AICc)

# save(preyC_growth_pe, file = "Data/FR assay/02 - LG search grid/PreyC_LG_search_grid.RData")

###### Best Model ######
# Get the best fit that had no errors
min_LG_fit_C <- preyC_growth_pe %>%
  filter(error == "NaN") %>%
  filter(AICc == min(preyC_growth_pe$AICc, na.rm = T)) %>%
  dplyr::select(c(r.est:sigma.est))

# Fit the LG model to prey growth using the estimated parameters of the best fit from the search grid
preyC_LG_best_fit = mle2(minuslogl = nll.odeint.LG.prey,
                        method = "BFGS",
                        control = list(maxit = 5000),
                        start = list(r = min_LG_fit_C$r.est,
                                     K.log = min_LG_fit_C$K.log,
                                     sigma = min_LG_fit_C$sigma),
                        data = with(list(N0 = Prey_start_density,
                                         Ndiff = Prey_end_density - Prey_start_density,
                                         Tt = Incubation_time),
                                    data = dplyr::filter(FR_Dataset_C, Predator_presence == 0)))

# save(preyC_LG_best_fit, file = "FR data/LG search grid/PreyC_LG_best_fit.RData")










#---------------------------------------------------------------------------------------------------
# Prey from selection environment S
#---------------------------------------------------------------------------------------------------
load("Data/FR assay/01 - Raw data/Dataset_S.RData")

growth_sim_grid <- pre_sg

###### Set the computer to run in Parallel ######
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores() - 1,
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered() # returns TRUE if a doPar backend has been registered
foreach::getDoParWorkers()    # nr. of execution workers in the currently registered doPar backend

###### Run Grid ######
preyS_growth_pe <- NA # name of file where parallel runs will be saved

(t0 <- Sys.time()) # timer
mcoptions <- list(preschedule=FALSE)

preyS_growth_pe <- foreach(i = 1:nrow(growth_sim_grid),
                           .combine = 'rbind',
                           .options.multicore = mcoptions,
                           .packages = c('bbmle', 'odeintr', 'tidyverse', 'R.utils'),
                           .errorhandling = "pass") %dopar%
  {
    compile_sys("LG_prey", LG.prey, pars = c("r", "K"), method = "rk54")
    
    r<-tryCatch({
      setTimeLimit(35*60)
      fit = mle2(minuslogl = nll.odeint.LG.prey,
                 method = "BFGS",
                 control = list(maxit = 5000),
                 start = list(r = growth_sim_grid$r[[i]],
                              K.log = growth_sim_grid$K.log[[i]],
                              sigma = growth_sim_grid$sigma[[i]]),
                 data = with(list(N0 = Prey_start_density,
                                  Ndiff = Prey_end_density - Prey_start_density,
                                  Tt = Incubation_time),
                             data = dplyr::filter(FR_Dataset_S, Predator_presence == 0)))
      
      
      if(fit@details$convergence==0) {
        growth_sim_grid$r.est[i] <- fit@coef[[1]]
        growth_sim_grid$K.log.est[i] <- fit@coef[[2]]
        
        growth_sim_grid$sigma.est [i] <- fit@coef[[3]]
        
        growth_sim_grid$iteration[i] <- fit@details[[3]][[2]]
        growth_sim_grid$converge[i] <- fit@details[[4]]
        growth_sim_grid$AICc[i] <- AICc(fit)
        growth_sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        growth_sim_grid$r.est[i] <- NaN
        growth_sim_grid$K.log.est[i] <- NaN
        
        growth_sim_grid$converge[i] <- NaN
        growth_sim_grid$iteration[i] <- NaN
        growth_sim_grid$AICc[i] <- NaN
        growth_sim_grid$logLik[i] <- NaN
      }
      
      #Remove fit for next iteration
      rm(fit)
    }, error=function(e) list(message = e$message)) # to print error in console use {cat("ERROR :", paste(conditionMessage(e), "Combo", (i)), "\n")}
    
    if(length(r) > 0){
      growth_sim_grid$error[i] <- r
    } else {
      growth_sim_grid$error[i] <- NaN
    }
    
    return(growth_sim_grid[i,])
    
  }

parallel::stopCluster(cl = my.cluster)
stopImplicitCluster()

tend <- Sys.time()
tend - t0
# beep(sound = 2, expr = NULL)

preyS_growth_pe %>% arrange(AICc)

# save(preyS_growth_pe, file = "Data/FR assay/02 - LG search grid/PreyS_LG_search_grid.RData")

###### Best Model ######
# Get the best fit that had no errors
min_LG_fit_S <- preyS_growth_pe %>%
  filter(error == "NaN") %>%
  filter(AICc == min(preyS_growth_pe$AICc, na.rm = T)) %>%
  dplyr::select(c(r.est:sigma.est))

# Fit the LG model to prey growth using the estimated parameters of the best fit from the search grid
preyS_LG_best_fit = mle2(minuslogl = nll.odeint.LG.prey,
                         method = "BFGS",
                         control = list(maxit = 5000),
                         start = list(r = min_LG_fit_S$r.est,
                                      K.log = min_LG_fit_S$K.log,
                                      sigma = min_LG_fit_S$sigma),
                         data = with(list(N0 = Prey_start_density,
                                          Ndiff = Prey_end_density - Prey_start_density,
                                          Tt = Incubation_time),
                                     data = dplyr::filter(FR_Dataset_S, Predator_presence == 0)))

# save(preyS_LG_best_fit, file = "FR data/LG search grid/PreyS_LG_best_fit.RData")











#---------------------------------------------------------------------------------------------------
# Prey from selection environment T
#---------------------------------------------------------------------------------------------------
load("Data/FR assay/01 - Raw data/Dataset_T.RData")

growth_sim_grid <- pre_sg

###### Set the computer to run in Parallel ######
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores() - 1,
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered() # returns TRUE if a doPar backend has been registered
foreach::getDoParWorkers()    # nr. of execution workers in the currently registered doPar backend

###### Run Grid ######
preyT_growth_pe <- NA # name of file where parallel runs will be saved

(t0 <- Sys.time()) # timer
mcoptions <- list(preschedule=FALSE)

preyT_growth_pe <- foreach(i = 1:nrow(growth_sim_grid),
                           .combine = 'rbind',
                           .options.multicore = mcoptions,
                           .packages = c('bbmle', 'odeintr', 'tidyverse', 'R.utils'),
                           .errorhandling = "pass") %dopar%
  {
    compile_sys("LG_prey", LG.prey, pars = c("r", "K"), method = "rk54")
    
    r<-tryCatch({
      setTimeLimit(35*60)
      fit = mle2(minuslogl = nll.odeint.LG.prey,
                 method = "BFGS",
                 control = list(maxit = 5000),
                 start = list(r = growth_sim_grid$r[[i]],
                              K.log = growth_sim_grid$K.log[[i]],
                              sigma = growth_sim_grid$sigma[[i]]),
                 data = with(list(N0 = Prey_start_density,
                                  Ndiff = Prey_end_density - Prey_start_density,
                                  Tt = Incubation_time),
                             data = dplyr::filter(FR_Dataset_T, Predator_presence == 0)))
      
      
      if(fit@details$convergence==0) {
        growth_sim_grid$r.est[i] <- fit@coef[[1]]
        growth_sim_grid$K.log.est[i] <- fit@coef[[2]]
        
        growth_sim_grid$sigma.est [i] <- fit@coef[[3]]
        
        growth_sim_grid$iteration[i] <- fit@details[[3]][[2]]
        growth_sim_grid$converge[i] <- fit@details[[4]]
        growth_sim_grid$AICc[i] <- AICc(fit)
        growth_sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        growth_sim_grid$r.est[i] <- NaN
        growth_sim_grid$K.log.est[i] <- NaN
        
        growth_sim_grid$converge[i] <- NaN
        growth_sim_grid$iteration[i] <- NaN
        growth_sim_grid$AICc[i] <- NaN
        growth_sim_grid$logLik[i] <- NaN
      }
      
      #Remove fit for next iteration
      rm(fit)
    }, error=function(e) list(message = e$message)) # to print error in console use {cat("ERROR :", paste(conditionMessage(e), "Combo", (i)), "\n")}
    
    if(length(r) > 0){
      growth_sim_grid$error[i] <- r
    } else {
      growth_sim_grid$error[i] <- NaN
    }
    
    return(growth_sim_grid[i,])
    
  }

parallel::stopCluster(cl = my.cluster)
stopImplicitCluster()

tend <- Sys.time()
tend - t0
# beep(sound = 2, expr = NULL)

preyT_growth_pe %>% arrange(AICc)

# save(preyT_growth_pe, file = "Data/FR assay/02 - LG search grid/PreyT_LG_search_grid.RData")

###### Best Model ######
# Get the best fit that had no errors
min_LG_fit_T <- preyT_growth_pe %>%
  filter(error == "NaN") %>%
  filter(AICc == min(preyT_growth_pe$AICc, na.rm = T)) %>%
  dplyr::select(c(r.est:sigma.est))

# Fit the LG model to prey growth using the estimated parameters of the best fit from the search grid
preyT_LG_best_fit = mle2(minuslogl = nll.odeint.LG.prey,
                         method = "BFGS",
                         control = list(maxit = 5000),
                         start = list(r = min_LG_fit_T$r.est,
                                      K.log = min_LG_fit_T$K.log,
                                      sigma = min_LG_fit_T$sigma),
                         data = with(list(N0 = Prey_start_density,
                                          Ndiff = Prey_end_density - Prey_start_density,
                                          Tt = Incubation_time),
                                     data = dplyr::filter(FR_Dataset_T, Predator_presence == 0)))

# save(preyT_LG_best_fit, file = "FR data/LG search grid/PreyT_LG_best_fit.RData")








#---------------------------------------------------------------------------------------------------
# Prey from selection environment ST
#---------------------------------------------------------------------------------------------------
load("Data/FR assay/01 - Raw data/Dataset_ST.RData")

growth_sim_grid <- pre_sg

###### Set the computer to run in Parallel ######
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores() - 1,
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered() # returns TRUE if a doPar backend has been registered
foreach::getDoParWorkers()    # nr. of execution workers in the currently registered doPar backend

###### Run Grid ######
preyST_growth_pe <- NA # name of file where parallel runs will be saved

(t0 <- Sys.time()) # timer
mcoptions <- list(preschedule=FALSE)

preyST_growth_pe <- foreach(i = 1:nrow(growth_sim_grid),
                            .combine = 'rbind',
                            .options.multicore = mcoptions,
                            .packages = c('bbmle', 'odeintr', 'tidyverse', 'R.utils'),
                            .errorhandling = "pass") %dopar%
  {
    compile_sys("LG_prey", LG.prey, pars = c("r", "K"), method = "rk54")
    
    r<-tryCatch({
      setTimeLimit(35*60)
      fit = mle2(minuslogl = nll.odeint.LG.prey,
                 method = "BFGS",
                 control = list(maxit = 5000),
                 start = list(r = growth_sim_grid$r[[i]],
                              K.log = growth_sim_grid$K.log[[i]],
                              sigma = growth_sim_grid$sigma[[i]]),
                 data = with(list(N0 = Prey_start_density,
                                  Ndiff = Prey_end_density - Prey_start_density,
                                  Tt = Incubation_time),
                             data = dplyr::filter(FR_Dataset_ST, Predator_presence == 0)))
      
      
      if(fit@details$convergence==0) {
        growth_sim_grid$r.est[i] <- fit@coef[[1]]
        growth_sim_grid$K.log.est[i] <- fit@coef[[2]]
        
        growth_sim_grid$sigma.est [i] <- fit@coef[[3]]
        
        growth_sim_grid$iteration[i] <- fit@details[[3]][[2]]
        growth_sim_grid$converge[i] <- fit@details[[4]]
        growth_sim_grid$AICc[i] <- AICc(fit)
        growth_sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        growth_sim_grid$r.est[i] <- NaN
        growth_sim_grid$K.log.est[i] <- NaN
        
        growth_sim_grid$converge[i] <- NaN
        growth_sim_grid$iteration[i] <- NaN
        growth_sim_grid$AICc[i] <- NaN
        growth_sim_grid$logLik[i] <- NaN
      }
      
      #Remove fit for next iteration
      rm(fit)
    }, error=function(e) list(message = e$message)) # to print error in console use {cat("ERROR :", paste(conditionMessage(e), "Combo", (i)), "\n")}
    
    if(length(r) > 0){
      growth_sim_grid$error[i] <- r
    } else {
      growth_sim_grid$error[i] <- NaN
    }
    
    return(growth_sim_grid[i,])
    
  }

parallel::stopCluster(cl = my.cluster)
stopImplicitCluster()

tend <- Sys.time()
tend - t0
# beep(sound = 2, expr = NULL)

preyST_growth_pe %>% arrange(AICc)

# save(preyST_growth_pe, file = "Data/FR assay/02 - LG search grid/PreyST_LG_search_grid.RData")

###### Best Model ######
# Get the best fit that had no errors
min_LG_fit_ST <- preyST_growth_pe %>%
  filter(error == "NaN") %>%
  filter(AICc == min(preyST_growth_pe$AICc, na.rm = T)) %>%
  dplyr::select(c(r.est:sigma.est))

# Fit the LG model to prey growth using the estimated parameters of the best fit from the search grid
preyST_LG_best_fit = mle2(minuslogl = nll.odeint.LG.prey,
                          method = "BFGS",
                          control = list(maxit = 5000),
                          start = list(r = min_LG_fit_ST$r.est,
                                       K.log = min_LG_fit_ST$K.log,
                                       sigma = min_LG_fit_ST$sigma),
                          data = with(list(N0 = Prey_start_density,
                                           Ndiff = Prey_end_density - Prey_start_density,
                                           Tt = Incubation_time),
                                      data = dplyr::filter(FR_Dataset_ST, Predator_presence == 0)))

# save(preyST_LG_best_fit, file = "FR data/LG search grid/PreyST_LG_best_fit.RData")




