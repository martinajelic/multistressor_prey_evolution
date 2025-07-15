# ----------------------------------------------------------------------------------------------------------------------------------
# Search grid of functional response parameters for predators fed with prey from different selection
# environments
# ----------------------------------------------------------------------------------------------------------------------------------
#
# The following R-code is based on the method of:
#
# Rosenbaum, B. & Rall, B.C. (2018). Fitting functional responses: Direct parameter estimation by
# simulating differential equations. Methods in Ecology and Evolution.
# 
# Lemmen, K. D. & Pennekamp, F. (2024). Food web context modifies predator foraging and weakens 
# trophic interaction strength.
#
# Daugaard, U., Petchey, O.L. & Pennekamp, F. (2019). Warming can destabilize predator–prey 
# interactions by shifting the functional response from Type III to Type II.
# 
# This method estimates the functional response parameters as functions of selection environment
# directly by simulating differential equations. The functional response parameters are estimated
# using an ANOVA-like fitting procedure, where:
#   - C-coefficients (bC.log, hC, mC) are the coefficients of functional response of predators fed 
#     with prey C
#   - S and T coefficients (bS.log, hS, mS, bT.log, hT, mT) are differences in parameters between 
#     predators fed with prey C and predators fed with prey S and prey T, respectively
#   - ST-coefficients (bST.log, hST, mST) are the differences between observed parameters of
#     predators fed with ST prey and the expected parameters if the effects of high selection 
#     salinity and high selection temperature were additive.
#
# For this, the package odeintr is used. For a manual about the use of the method and a detailed
# description of the functions, see the paper of Rosenbaum and Rall (2018) mentioned above.
# 
# NOTE: the package odeintr is available on github (thk686/odeintr) for R version 4.3.2 and is 
# unsupported for newer versions of R at the moment.
# 
# Rdata files with the output from the search grid process can be found in the folder
# "Data/FR assay/03 - FR search grid" and can be manually loaded into the workspace.
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
library(broom)

load("Data/FR assay/01 - Raw data/Dataset_all.RData")

# Remove NAs to avoid model fitting problems
FR_all <- FR_all %>%
  filter(!is.na(Predator_end_density))

source("R scripts/Function_FR_model_ANOVA_like.R")

# ------------------------------------------------------------------------------
# All r and K estimates have been searched previously and individually for each
# prey selection environment (see R scripts/01_LG_search_grids.R). These values
# are fixed when fitting the ANOVA-like FR function
# ------------------------------------------------------------------------------
load("Data/FR assay/02 - LG search grid/PreyC_LG_best_fit.RData")
load("Data/FR assay/02 - LG search grid/PreyS_LG_best_fit.RData")
load("Data/FR assay/02 - LG search grid/PreyT_LG_best_fit.RData")
load("Data/FR assay/02 - LG search grid/PreyST_LG_best_fit.RData")

##### Search Grid ##### 
fr_pre_sg <- expand.grid(bC.log = c(seq(-5, 5, length.out = 6)),
                         bS.log = 0,
                         bT.log = 0,
                         bST.log = 0,
                         hC = c((seq(0.0000001, 1, length.out = 6))),
                         hS = 0,
                         hT = 0,
                         hST = 0,
                         mC = seq(0.000001, 3, length.out = 6),
                         mS = 0,
                         mT = 0,
                         mST = 0,
                         c = exp(-15),
                         sigma = 0.282819,
                         rC = preyC_LG_best_fit@coef[["r"]],
                         rS = preyS_LG_best_fit@coef[["r"]],
                         rT = preyT_LG_best_fit@coef[["r"]],
                         rST = preyST_LG_best_fit@coef[["r"]],
                         KC.log = preyC_LG_best_fit@coef[["K.log"]],
                         KS.log = preyS_LG_best_fit@coef[["K.log"]],
                         KT.log = preyT_LG_best_fit@coef[["K.log"]],
                         KST.log = preyST_LG_best_fit@coef[["K.log"]])

fr_pre_sg$error <- NA

fr_pre_sg$bC.log.est <- NA
fr_pre_sg$hC.est <- NA
fr_pre_sg$mC.est <- NA


fr_pre_sg$bS.log.est <- NA
fr_pre_sg$hS.est <- NA
fr_pre_sg$mS.est <- NA

fr_pre_sg$bT.log.est <- NA
fr_pre_sg$hT.est <- NA
fr_pre_sg$mT.est <- NA

fr_pre_sg$bST.log.est <- NA
fr_pre_sg$hST.est <- NA
fr_pre_sg$mST.est <- NA

fr_pre_sg$c.est <- NA
fr_pre_sg$sigma.est <- NA

fr_pre_sg$converge <- NA
fr_pre_sg$AICc <- NA
fr_pre_sg$AIC <- NA
fr_pre_sg$logLik <- NA
fr_pre_sg$strtpt <- seq(1:nrow(fr_pre_sg))

###### Run Grid Search #####
fr_sg <- fr_pre_sg

bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(), 
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

FR_search_grid <- NA # name of file where parallel runs will be saved

(t0<-Sys.time()) # timer
mcoptions <- list(preschedule = FALSE)

FR_search_grid <- foreach(i = 1:nrow(fr_sg),
                 .combine = 'rbind',
                 .options.multicore = mcoptions,
                 .packages = c('bbmle', 'odeintr', 'tidyverse', 'R.utils'),
                 .errorhandling = "pass") %dopar%
  {
    compile_sys("FRII_gen_pred_mP", FRII.gen.pred.mP,
                pars = c("b", "h", "r", "K", "c", "m"),
                method = "rk54")
    
    r<-tryCatch({ 
      setTimeLimit(35*60)
      fit = mle2(minuslogl = nll.odeint.full.model,
                 method = "BFGS",
                 browse_obj = F,
                 skip.hessian = F, 
                 trace = T,
                 control = list(trace = 1, REPORT = 1, maxit = 5000),
                 fixed = list(rC = fr_sg$rC[[i]],
                              rS = fr_sg$rS[[i]],
                              rT = fr_sg$rT[[i]],
                              rST = fr_sg$rST[[i]],
                              KC.log = fr_sg$KC.log[[i]],
                              KS.log = fr_sg$KS.log[[i]],
                              KT.log = fr_sg$KT.log[[i]],
                              KST.log = fr_sg$KST.log[[i]]),
                 start = list(bC.log = fr_sg$bC.log[[i]], 
                              bS.log = fr_sg$bS.log[[i]],
                              bT.log = fr_sg$bT.log[[i]],
                              bST.log = fr_sg$bST.log[[i]],
                              hC = fr_sg$hC[[i]],
                              hS = fr_sg$hS[[i]],
                              hT = fr_sg$hT[[i]],
                              hST = fr_sg$hST[[i]],
                              mC = fr_sg$mC[[i]],
                              mS = fr_sg$mS[[i]],
                              mT = fr_sg$mT[[i]],
                              mST = fr_sg$mST[[i]],
                              sigma = fr_sg$sigma[[i]],
                              c = fr_sg$c[[i]]),
                 data = with(list(N0 = Prey_start_density,
                                  Ndead = Prey_start_density - Prey_end_density,
                                  P = Predator_start_density,
                                  P.end = Predator_end_density * 0.66, # Need consumer counts (0.66mL of 1mL of sample were counted)
                                  Tt = Incubation_time,
                                  SEL_environment = SEL_environment,
                                  separate_environments = FALSE),
                             data = dplyr::filter(FR_all, Predator_presence == 1)))
      
      
      if(fit@details$convergence==0) {
        fr_sg$bC.log.est[i] <- fit@coef[["bC.log"]]
        fr_sg$bS.log.est[i] <- fit@coef[["bS.log"]]
        fr_sg$bT.log.est[i] <- fit@coef[["bT.log"]]
        fr_sg$bST.log.est[i] <- fit@coef[["bST.log"]]
        
        fr_sg$hC.est[i] <- fit@coef[["hC"]]
        fr_sg$hS.est[i] <- fit@coef[["hS"]]
        fr_sg$hT.est[i] <- fit@coef[["hT"]]
        fr_sg$hST.est[i] <- fit@coef[["hST"]]
        
        fr_sg$mC.est[i] <- fit@coef[["mC"]]
        fr_sg$mS.est[i] <- fit@coef[["mS"]]
        fr_sg$mT.est[i] <- fit@coef[["mT"]]
        fr_sg$mST.est[i] <- fit@coef[["mST"]]
        
        fr_sg$c.est[i] <- fit@coef[["c"]]
        fr_sg$sigma.est[i] <- fit@coef[["sigma"]]
        
        fr_sg$converge[i] <- fit@details[[4]]
        fr_sg$AICc[i] <- MuMIn::AICc(fit)
        fr_sg$AIC[i] <- AIC(fit)
        fr_sg$logLik[i] <- logLik(fit)[[1]]
      } else {
        fr_sg$bC.log.est[i] <- NaN
        fr_sg$bS.log.est[i] <- NaN
        fr_sg$bT.log.est[i] <- NaN
        fr_sg$bST.log.est[i] <- NaN
        
        fr_sg$hC.est[i] <- NaN
        fr_sg$hS.est[i] <- NaN
        fr_sg$hT.est[i] <- NaN
        fr_sg$hST.est[i] <- NaN
        
        fr_sg$mC.est[i] <- NaN
        fr_sg$mS.est[i] <- NaN
        fr_sg$mT.est[i] <- NaN
        fr_sg$mST.est[i] <- NaN
        
        fr_sg$c.est[i] <- NaN
        fr_sg$sigma.est[i] <- NaN
        
        fr_sg$converge[i] <- NaN
        fr_sg$AICc[i] <- NaN
        fr_sg$AIC[i] <- NaN
        fr_sg$logLik[i] <- NaN
      }
      
      #Remove fit for next iteration
      rm(fit)
    }, error=function(e) list(message = e$message)) # to print error in console use {cat("ERROR :", paste(conditionMessage(e), "Combo", (i)), "\n")}
    
    if(length(r) > 0){
      fr_sg$error[i] <- r 
    } else {
      fr_sg$error[i] <- NaN
    }  
    
    return(fr_sg[i,])
    
  }

parallel::stopCluster(cl = my.cluster)
stopImplicitCluster()

tend<-Sys.time()
tend-t0 
# beep(sound = 2, expr = NULL)

FR_search_grid %>% arrange(AICc)

save(FR_search_grid, file = "Data/FR assay/03 - FR search grid/FR_ANOVA_like_search_grid.RData")


###### Best Model ######
# Get the best fit that had no errors, based on AICc
min_FR_params <- FR_search_grid %>% 
  filter(error == "NaN") %>% 
  filter(AICc == min(FR_search_grid$AICc, na.rm = T)) %>%
  dplyr::select(c(bC.log.est:sigma.est))

# Fitting the FR-TypeII model using the estimated parameters of the best fit from the search grid
best_FR_fit = mle2(minuslogl = nll.odeint.full.model,
                           method = "BFGS",
                           browse_obj = F,
                           skip.hessian = F, 
                           trace = T,
                           control = list(trace=1, REPORT=1, maxit = 5000),
                           fixed = list(rC = preyC_LG_best_fit@coef[["r"]],
                                        rS = preyS_LG_best_fit@coef[["r"]],
                                        rT = preyT_LG_best_fit@coef[["r"]],
                                        rST = preyST_LG_best_fit@coef[["r"]],
                                        KC.log = preyC_LG_best_fit@coef[["K.log"]],
                                        KS.log = preyS_LG_best_fit@coef[["K.log"]],
                                        KT.log = preyT_LG_best_fit@coef[["K.log"]],
                                        KST.log = preyST_LG_best_fit@coef[["K.log"]]),
                           start = list(bC.log = min_FR_params$bC.log.est, 
                                        bS.log = min_FR_params$bS.log.est,
                                        bT.log = min_FR_params$bT.log.est,
                                        bST.log = min_FR_params$bST.log.est,
                                        hC = min_FR_params$hC.est,
                                        hS = min_FR_params$hS.est,
                                        hT = min_FR_params$hT.est,
                                        hST = min_FR_params$hST.est,
                                        mC = min_FR_params$mC.est,
                                        mS = min_FR_params$mS.est,
                                        mT = min_FR_params$mT.est,
                                        mST = min_FR_params$mST.est,
                                        sigma = min_FR_params$sigma.est,
                                        c = min_FR_params$c.est),
                           data = with(list(N0 = Prey_start_density,
                                            Ndead = Prey_start_density - Prey_end_density,
                                            P = Predator_start_density,
                                            P.end = Predator_end_density * 0.66, # Need consumer counts (0.66mL of 1mL of sample were counted)
                                            Tt = Incubation_time,
                                            SEL_environment = SEL_environment,
                                            separate_environments = FALSE),
                                       data = dplyr::filter(FR_all, Predator_presence == 1)))

summary(best_FR_fit)

save(best_FR_fit, file = "Data/FR assay/03 - FR search grid/FR_ANOVA_like_best_fit.RData")

