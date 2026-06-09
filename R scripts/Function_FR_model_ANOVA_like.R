# ----------------------------------------------------------------------------------------------------------------------------------
# Logistic growth model for estimating prey r and K from the controls
# ----------------------------------------------------------------------------------------------------------------------------------
#
# The following R-code of the source functions is based on the method of:
#
# Rosenbaum, B. & Rall, B.C. (2018). Fitting functional responses: Direct parameter estimation by simulating differential equations. 
# Methods in Ecology and Evolution.
#
# Daugaard, U., Petchey, O.L. & Pennekamp, F. (2019). Warming can destabilize predator–prey 
# interactions by shifting the functional response from Type III to Type II.
#
# NOTE: the "odeintr" package is available on github for R version 4.3.2 and is unsupported for newer 
# versions of R at the moment. A working version of the package can be installed with the following
# command:
# install.packages("R versions/R version for search grids and ODE fitting/odeintr_1.7.1_patched_clean.tar.gz", repos = NULL, type = "source")
# 
# ----------------------------------------------------------------------------------------------------------------------------------

Sys.setenv("PKG_CXXFLAGS"="-std=gnu++14")

library(odeintr)

# Define ODEs for prey and predator (type II FR)
# x[1] = predator
# x[0] = prey
# b = space clearance rate
# h = handling time
# r = prey intrinsic growth rate
# K = prey carrying capacity
# c = conversion efficiency
# m = predator mortality rate

FRII.gen.pred.mP = '
dxdt[0] = -b * x[0] / (1 + b * h * x[0]) * x[1] + r * x[0] * (1.0 - x[0]/K);
dxdt[1] = x[1] * c * b * x[0] / (1 + b * h * x[0]) - m * x[1];
'

# Create FR_gen_pred (integrated function), using the defined C++ expressions FR.gen.pred
compile_sys("FRII_gen_pred_mP", FRII.gen.pred.mP, pars = c("b","h","r","K","c","m"), method = "rk54")

# Create function to simulate predator and prey densities in time, given the ODEs
eq.odeint.general.pred = function(start.v, b, h, r, K, c, m, Tt, timesteplength){
  FRII_gen_pred_mP_set_params(b = b, h = h, r=r, K=K, c=c, m=m)
  FRII_gen_pred_mP(start.v,Tt,timesteplength)
}

# Create function to get estimated prey consumption and predator growth:
#   Neaten.est: a vector of estimated consumed prey individuals for each sample (prey death)
#   predrep.est: a vector of estimated new predators for each sample (predator growth)
#   returns a list with Neaten.est and predrep.est
eaten.odeint.general.pred = function(N0, b, h, r, K, c, m, Tt, P, steps = 100){
  Neaten.est = vector()
  predrep.est = vector()
  
  for(i.eaten in 1:length(N0)){
    # Simulate predator and prey densities over 100 timesteps
    temp = eq.odeint.general.pred(start.v = c(N0[i.eaten], P[i.eaten]),
                                        b = b,
                                        h = h,
                                        r = r,
                                        K = K,
                                        c = c,
                                        m = m,
                                        Tt = Tt[i.eaten],
                                        timesteplength = Tt[i.eaten]/steps)
    
    # Save estimated prey consumption and predator growth
    Neaten.est[i.eaten] = N0[i.eaten] - temp[steps+1,2] # how many at the beginning - how many at the end = how many eaten
    predrep.est[i.eaten] = temp[steps+1,3] - P[i.eaten] # how many at the end - how many at the beginning = how many new
  }
  
  ret = list(Neaten.est, predrep.est)
  return(ret)
}

# The following function computes the negative log-likelihood for a given set of parameters
# There is an option to choose the reference level for different prey selection environments. Note 
# that the meaning of the coefficients changes based on reference environment (see comments).
nll.odeint.full.model = function(Ndead, N0,
                                           bC.log, bS.log, bT.log, bST.log,
                                           hC, hS, hT, hST, 
                                           rC, rS, rT, rST, 
                                           KC.log, KS.log, KT.log, KST.log,
                                           c,
                                           mC, mS, mT, mST, Tt, 
                                           P, P.end, steps = 100, sigma,
                                           SEL_environment,
                                           reference_environment = "C"){ # Default reference is environment C
  if(sigma <= 0) return(Inf)
  
  if(reference_environment == "C") {
    # If reference_environment == "C", the S and T parameters are the effects of increasing 
    # selection temperature or selection salinity
    # The multi-stressor environment is selection environment ST (high-temperature, high-salinity)
    # 
    # C-parameters = reference
    # S-parameters = effect of adding salinity
    # T-parameters = effect of adding temperature
    # ST-parameters = interactive effect of adding temperature and salinity
    
    ### Environment C
    # Get estimated prey consumed and predator growth given vital rates b, h, r, K, and c
    # Parameters for S, T and ST are all zero
    tempC = eaten.odeint.general.pred(N0 = N0[SEL_environment == "C"],        # Prey start dens
                                      b = exp(bC.log),
                                      h = hC,
                                      r = rC,
                                      K = exp(KC.log),
                                      c = c,
                                      m = mC,
                                      Tt = Tt[SEL_environment == "C"],        # Incubation time
                                      P = P[SEL_environment == "C"],          # Predator start dens
                                      steps = steps)  # Number of timesteps to simulate in Tt
    
    ### Environment S
    # Get estimated prey consumed and predator growth given vital rates b, h, r, K, and c
    # Parameters for T and ST are all zero
    tempS = eaten.odeint.general.pred(N0 = N0[SEL_environment == "S"],        # Prey start dens
                                      b = exp(bC.log + bS.log),
                                      h = hC + hS,
                                      r = rS,              # independent r's for SEL environments
                                      K = exp(KS.log),     # independent K's for SEL environments
                                      c = c,
                                      m = mC + mS,
                                      Tt = Tt[SEL_environment == "S"],        # Incubation time
                                      P = P[SEL_environment == "S"],          # Predator start dens
                                      steps = steps)  # Number of timesteps to simulate in Tt
    
    ### Environment T
    # Get estimated prey consumed and predator growth given vital rates b, h, r, K, and c
    # Parameters for S and ST are all zero
    tempT = eaten.odeint.general.pred(N0 = N0[SEL_environment == "T"],        # Prey start dens
                                      b = exp(bC.log + bT.log),
                                      h = hC + hT,
                                      r = rT,              # independent r's for SEL environments
                                      K = exp(KT.log),     # independent K's for SEL environments
                                      c = c,
                                      m = mC + mT,
                                      Tt = Tt[SEL_environment == "T"],        # Incubation time
                                      P = P[SEL_environment == "T"],          # Predator start dens
                                      steps = steps)  # Number of timesteps to simulate in Tt
    
    ### Environment ST
    # Get estimated prey consumed and predator growth given vital rates b, h, r, K, and c
    tempST = eaten.odeint.general.pred(N0 = N0[SEL_environment == "ST"],        # Prey start dens
                                       b = exp(bC.log + bS.log + bT.log + bST.log),
                                       h = hC + hS + hT + hST,
                                       r = rST,              # independent r's for SEL environments
                                       K = exp(KST.log),     # independent K's for SEL environments
                                       c = c,
                                       m = mC + mS + mT + mST,
                                       Tt = Tt[SEL_environment == "ST"],        # Incubation time
                                       P = P[SEL_environment == "ST"],          # Predator start dens
                                       steps = steps)  # Number of timesteps to simulate in Tt
    
    # Compute nll for each environment and sum them for general nll
    
    ### Environment C
    yC = tempC[[1]] # N.eaten estimate
    zC = tempC[[2]] # New predators in 1mL estimate
    
    nllC = -1*sum(dnorm(x = log(N0[SEL_environment == "C"] - Ndead[SEL_environment == "C"]), # Measured prey end density
                        mean = log(N0[SEL_environment == "C"] - yC),  # Estimated prey end density (Start prey dens + estimated change)
                        sd = sigma,
                        log = T))
    
    nll2C = -1*sum(dpois(x = P.end[SEL_environment == "C"],           # Measured predator end counts in 0.66 mL
                         lambda = (P[SEL_environment == "C"] + zC) * 0.66, # Estimated predator end counts in 0.66 mL (Start pred counts + estimated new preds)
                         log = T))
    
    ### Environment S
    yS = tempS[[1]] # N.eaten estimate
    zS = tempS[[2]] # New predators in 1mL estimate
    
    nllS = -1*sum(dnorm(x = log(N0[SEL_environment == "S"] - Ndead[SEL_environment == "S"]), # Measured prey end density
                        mean = log(N0[SEL_environment == "S"] - yS),  # Estimated prey end density (Start prey dens + estimated change)
                        sd = sigma,
                        log = T))
    
    nll2S = -1*sum(dpois(x = P.end[SEL_environment == "S"],           # Measured predator end counts in 0.66 mL
                         lambda = (P[SEL_environment == "S"] + zS) * 0.66, # Estimated predator end counts in 0.66 mL (Start pred counts + estimated new preds)
                         log = T))
    
    ### Environment T
    yT = tempT[[1]] # N.eaten estimate
    zT = tempT[[2]] # New predators in 1mL estimate
    
    nllT = -1*sum(dnorm(x = log(N0[SEL_environment == "T"] - Ndead[SEL_environment == "T"]), # Measured prey end density
                        mean = log(N0[SEL_environment == "T"] - yT),  # Estimated prey end density (Start prey dens + estimated change)
                        sd = sigma,
                        log = T))
    
    nll2T = -1*sum(dpois(x = P.end[SEL_environment == "T"],           # Measured predator end counts in 0.66 mL
                         lambda = (P[SEL_environment == "T"] + zT) * 0.66, # Estimated predator end counts in 0.66 mL (Start pred counts + estimated new preds)
                         log = T))
    
    
    ### Environment ST
    yST = tempST[[1]] # N.eaten estimate
    zST = tempST[[2]] # New predators in 1mL estimate
    
    nllST = -1*sum(dnorm(x = log(N0[SEL_environment == "ST"] - Ndead[SEL_environment == "ST"]), # Measured prey end density
                         mean = log(N0[SEL_environment == "ST"] - yST),  # Estimated prey end density (Start prey dens + estimated change)
                         sd = sigma,
                         log = T))
    
    nll2ST = -1*sum(dpois(x = P.end[SEL_environment == "ST"],           # Measured predator end counts in 0.66 mL
                          lambda = (P[SEL_environment == "ST"] + zST) * 0.66, # Estimated predator end counts in 0.66 mL (Start pred counts + estimated new preds)
                          log = T))
    
    # General nll value to minimize
    nll = nllC + nll2C + 
      nllS + nll2S + 
      nllT + nll2T + 
      nllST + nll2ST
    
    return(nll)
  }
  
  if (reference_environment == "S") {
    # If reference_environment == "S", the C parameters are the effect of decreasing salinity, 
    # while the ST-parameters are the effect of increasing temperature
    # The multi-stressor environment is now selection environment T (low-salinity, high-temperature)
    # 
    # C-parameters = effect of removing salinity
    # S-parameters = reference
    # T-parameters = interactive effect of adding temperature and removing salinity
    # ST-parameters = effect of adding temperature
    
    ### Environment C
    # Get estimated prey consumed and predator growth given vital rates b, h, r, K, and c
    # We only need to remove the effect of salinity (represented by C-parameters)
    tempC = eaten.odeint.general.pred(N0 = N0[SEL_environment == "C"],        # Prey start dens
                                      b = exp(bS.log + bC.log),
                                      h = hS + hC,
                                      r = rC,
                                      K = exp(KC.log),
                                      c = c,
                                      m = mS + mC,
                                      Tt = Tt[SEL_environment == "C"],        # Incubation time
                                      P = P[SEL_environment == "C"],          # Predator start dens
                                      steps = steps)  # Number of timesteps to simulate in Tt
    
    ### Environment S
    # Get estimated prey consumed and predator growth given vital rates b, h, r, K, and c
    tempS = eaten.odeint.general.pred(N0 = N0[SEL_environment == "S"],        # Prey start dens
                                      b = exp(bS.log),
                                      h = hS,
                                      r = rS,              # independent r's for SEL environments
                                      K = exp(KS.log),     # independent K's for SEL environments
                                      c = c,
                                      m = mS,
                                      Tt = Tt[SEL_environment == "S"],        # Incubation time
                                      P = P[SEL_environment == "S"],          # Predator start dens
                                      steps = steps)  # Number of timesteps to simulate in Tt
    
    ### Environment T
    # Get estimated prey consumed and predator growth given vital rates b, h, r, K, and c
    # We need to remove the effect of salinity and add the effect of temperature, plus interaction
    tempT = eaten.odeint.general.pred(N0 = N0[SEL_environment == "T"],        # Prey start dens
                                      b = exp(bS.log + bT.log + bC.log + bST.log),
                                      h = hS + hT + hC + hST,
                                      r = rT,              # independent r's for SEL environments
                                      K = exp(KT.log),     # independent K's for SEL environments
                                      c = c,
                                      m = mS + mT + mC + mST,
                                      Tt = Tt[SEL_environment == "T"],        # Incubation time
                                      P = P[SEL_environment == "T"],          # Predator start dens
                                      steps = steps)  # Number of timesteps to simulate in Tt
    
    ### Environment ST
    # Get estimated prey consumed and predator growth given vital rates b, h, r, K, and c
    # We only need to add the effect of temperature
    tempST = eaten.odeint.general.pred(N0 = N0[SEL_environment == "ST"],        # Prey start dens
                                       b = exp(bS.log + bT.log),
                                       h = hS + hT,
                                       r = rST,              # independent r's for SEL environments
                                       K = exp(KST.log),     # independent K's for SEL environments
                                       c = c,
                                       m = mS + mT,
                                       Tt = Tt[SEL_environment == "ST"],        # Incubation time
                                       P = P[SEL_environment == "ST"],          # Predator start dens
                                       steps = steps)  # Number of timesteps to simulate in Tt
    
    # Compute nll for each environment and sum them for general nll
    
    ### Environment C
    yC = tempC[[1]] # N.eaten estimate
    zC = tempC[[2]] # New predators in 1mL estimate
    
    nllC = -1*sum(dnorm(x = log(N0[SEL_environment == "C"] - Ndead[SEL_environment == "C"]), # Measured prey end density
                        mean = log(N0[SEL_environment == "C"] - yC),  # Estimated prey end density (Start prey dens + estimated change)
                        sd = sigma,
                        log = T))
    
    nll2C = -1*sum(dpois(x = P.end[SEL_environment == "C"],           # Measured predator end counts in 0.66 mL
                         lambda = (P[SEL_environment == "C"] + zC) * 0.66, # Estimated predator end counts in 0.66 mL (Start pred counts + estimated new preds)
                         log = T))
    
    ### Environment S
    yS = tempS[[1]] # N.eaten estimate
    zS = tempS[[2]] # New predators in 1mL estimate
    
    nllS = -1*sum(dnorm(x = log(N0[SEL_environment == "S"] - Ndead[SEL_environment == "S"]), # Measured prey end density
                        mean = log(N0[SEL_environment == "S"] - yS),  # Estimated prey end density (Start prey dens + estimated change)
                        sd = sigma,
                        log = T))
    
    nll2S = -1*sum(dpois(x = P.end[SEL_environment == "S"],           # Measured predator end counts in 0.66 mL
                         lambda = (P[SEL_environment == "S"] + zS) * 0.66, # Estimated predator end counts in 0.66 mL (Start pred counts + estimated new preds)
                         log = T))
    
    ### Environment T
    yT = tempT[[1]] # N.eaten estimate
    zT = tempT[[2]] # New predators in 1mL estimate
    
    nllT = -1*sum(dnorm(x = log(N0[SEL_environment == "T"] - Ndead[SEL_environment == "T"]), # Measured prey end density
                        mean = log(N0[SEL_environment == "T"] - yT),  # Estimated prey end density (Start prey dens + estimated change)
                        sd = sigma,
                        log = T))
    
    nll2T = -1*sum(dpois(x = P.end[SEL_environment == "T"],           # Measured predator end counts in 0.66 mL
                         lambda = (P[SEL_environment == "T"] + zT) * 0.66, # Estimated predator end counts in 0.66 mL (Start pred counts + estimated new preds)
                         log = T))
    
    
    ### Environment ST
    yST = tempST[[1]] # N.eaten estimate
    zST = tempST[[2]] # New predators in 1mL estimate
    
    nllST = -1*sum(dnorm(x = log(N0[SEL_environment == "ST"] - Ndead[SEL_environment == "ST"]), # Measured prey end density
                         mean = log(N0[SEL_environment == "ST"] - yST),  # Estimated prey end density (Start prey dens + estimated change)
                         sd = sigma,
                         log = T))
    
    nll2ST = -1*sum(dpois(x = P.end[SEL_environment == "ST"],           # Measured predator end counts in 0.66 mL
                          lambda = (P[SEL_environment == "ST"] + zST) * 0.66, # Estimated predator end counts in 0.66 mL (Start pred counts + estimated new preds)
                          log = T))
    
    # General nll value to minimize
    nll = nllC + nll2C + 
      nllS + nll2S + 
      nllT + nll2T + 
      nllST + nll2ST
    
    return(nll)
  }
  
  if (reference_environment == "T") {
    # If reference_environment == "T", the C parameters are the effect of decreasing temperature, 
    # while the ST-parameters are the effect of increasing salinity
    # The multi-stressor environment is now selection environment S (low-temperature, high-salinity)
    # 
    # C-parameters = effect of removing temperature
    # S-parameters = interactive effect of adding salinity and removing temperature
    # T-parameters = reference
    # ST-parameters = effect of adding salinity
    
    ### Environment C
    # Get estimated prey consumed and predator growth given vital rates b, h, r, K, and c
    # We only need to remove the effect of temperature (represented by C-parameters)
    tempC = eaten.odeint.general.pred(N0 = N0[SEL_environment == "C"],        # Prey start dens
                                      b = exp(bT.log + bC.log),
                                      h = hT + hC,
                                      r = rC,
                                      K = exp(KC.log),
                                      c = c,
                                      m = mT + mC,
                                      Tt = Tt[SEL_environment == "C"],        # Incubation time
                                      P = P[SEL_environment == "C"],          # Predator start dens
                                      steps = steps)  # Number of timesteps to simulate in Tt
    
    ### Environment S
    # Get estimated prey consumed and predator growth given vital rates b, h, r, K, and c
    # We need to remove the effect of temperature, and add the effect of salinity, plus interaction
    tempS = eaten.odeint.general.pred(N0 = N0[SEL_environment == "S"],        # Prey start dens
                                      b = exp(bT.log + bC.log + bS.log + bST.log),
                                      h = hT + hC + hS + hST,
                                      r = rS,              # independent r's for SEL environments
                                      K = exp(KS.log),     # independent K's for SEL environments
                                      c = c,
                                      m = mT +mC + mS + mST,
                                      Tt = Tt[SEL_environment == "S"],        # Incubation time
                                      P = P[SEL_environment == "S"],          # Predator start dens
                                      steps = steps)  # Number of timesteps to simulate in Tt
    
    ### Environment T
    # Get estimated prey consumed and predator growth given vital rates b, h, r, K, and c
    tempT = eaten.odeint.general.pred(N0 = N0[SEL_environment == "T"],        # Prey start dens
                                      b = exp(bT.log),
                                      h = hT,
                                      r = rT,              # independent r's for SEL environments
                                      K = exp(KT.log),     # independent K's for SEL environments
                                      c = c,
                                      m = mT,
                                      Tt = Tt[SEL_environment == "T"],        # Incubation time
                                      P = P[SEL_environment == "T"],          # Predator start dens
                                      steps = steps)  # Number of timesteps to simulate in Tt
    
    ### Environment ST
    # Get estimated prey consumed and predator growth given vital rates b, h, r, K, and c
    # We only need to add the effect of salinity
    tempST = eaten.odeint.general.pred(N0 = N0[SEL_environment == "ST"],        # Prey start dens
                                       b = exp(bS.log + bT.log),
                                       h = hS + hT,
                                       r = rST,              # independent r's for SEL environments
                                       K = exp(KST.log),     # independent K's for SEL environments
                                       c = c,
                                       m = mS + mT,
                                       Tt = Tt[SEL_environment == "ST"],        # Incubation time
                                       P = P[SEL_environment == "ST"],          # Predator start dens
                                       steps = steps)  # Number of timesteps to simulate in Tt
    
    # Compute nll for each environment and sum them for general nll
    
    ### Environment C
    yC = tempC[[1]] # N.eaten estimate
    zC = tempC[[2]] # New predators in 1mL estimate
    
    nllC = -1*sum(dnorm(x = log(N0[SEL_environment == "C"] - Ndead[SEL_environment == "C"]), # Measured prey end density
                        mean = log(N0[SEL_environment == "C"] - yC),  # Estimated prey end density (Start prey dens + estimated change)
                        sd = sigma,
                        log = T))
    
    nll2C = -1*sum(dpois(x = P.end[SEL_environment == "C"],           # Measured predator end counts in 0.66 mL
                         lambda = (P[SEL_environment == "C"] + zC) * 0.66, # Estimated predator end counts in 0.66 mL (Start pred counts + estimated new preds)
                         log = T))
    
    ### Environment S
    yS = tempS[[1]] # N.eaten estimate
    zS = tempS[[2]] # New predators in 1mL estimate
    
    nllS = -1*sum(dnorm(x = log(N0[SEL_environment == "S"] - Ndead[SEL_environment == "S"]), # Measured prey end density
                        mean = log(N0[SEL_environment == "S"] - yS),  # Estimated prey end density (Start prey dens + estimated change)
                        sd = sigma,
                        log = T))
    
    nll2S = -1*sum(dpois(x = P.end[SEL_environment == "S"],           # Measured predator end counts in 0.66 mL
                         lambda = (P[SEL_environment == "S"] + zS) * 0.66, # Estimated predator end counts in 0.66 mL (Start pred counts + estimated new preds)
                         log = T))
    
    ### Environment T
    yT = tempT[[1]] # N.eaten estimate
    zT = tempT[[2]] # New predators in 1mL estimate
    
    nllT = -1*sum(dnorm(x = log(N0[SEL_environment == "T"] - Ndead[SEL_environment == "T"]), # Measured prey end density
                        mean = log(N0[SEL_environment == "T"] - yT),  # Estimated prey end density (Start prey dens + estimated change)
                        sd = sigma,
                        log = T))
    
    nll2T = -1*sum(dpois(x = P.end[SEL_environment == "T"],           # Measured predator end counts in 0.66 mL
                         lambda = (P[SEL_environment == "T"] + zT) * 0.66, # Estimated predator end counts in 0.66 mL (Start pred counts + estimated new preds)
                         log = T))
    
    
    ### Environment ST
    yST = tempST[[1]] # N.eaten estimate
    zST = tempST[[2]] # New predators in 1mL estimate
    
    nllST = -1*sum(dnorm(x = log(N0[SEL_environment == "ST"] - Ndead[SEL_environment == "ST"]), # Measured prey end density
                         mean = log(N0[SEL_environment == "ST"] - yST),  # Estimated prey end density (Start prey dens + estimated change)
                         sd = sigma,
                         log = T))
    
    nll2ST = -1*sum(dpois(x = P.end[SEL_environment == "ST"],           # Measured predator end counts in 0.66 mL
                          lambda = (P[SEL_environment == "ST"] + zST) * 0.66, # Estimated predator end counts in 0.66 mL (Start pred counts + estimated new preds)
                          log = T))
    
    # General nll value to minimize
    nll = nllC + nll2C + 
      nllS + nll2S + 
      nllT + nll2T + 
      nllST + nll2ST
    
    return(nll)
  }
  
  if(reference_environment == "ST") {
    # If reference_environment == "ST", the S and T parameters are the effects of decreasing 
    # selection temperature or selection salinity
    # The multi-stressor environment is now selection environment C (low-temperature, low-salinity)
    # 
    # C-parameters = interactive effect of removing both stressors
    # S-parameters = effect of removing salinity
    # T-parameters = effect of removing temperature
    # ST-parameters = reference
    
    ### Environment C
    # Get estimated prey consumed and predator growth given vital rates b, h, r, K, and c
    # We need to remove the effect of high temperature and the effect of high salinity, plus
    # interaction
    tempC = eaten.odeint.general.pred(N0 = N0[SEL_environment == "C"],        # Prey start dens
                                      b = exp(bC.log + bS.log + bT.log + bST.log),
                                      h = hC + hS + hT + hST,
                                      r = rC,              # independent r's for SEL environments
                                      K = exp(KC.log),     # independent K's for SEL environments
                                      c = c,
                                      m = mC + mS + mT + mST,
                                      Tt = Tt[SEL_environment == "C"],        # Incubation time
                                      P = P[SEL_environment == "C"],          # Predator start dens
                                      steps = steps)  # Number of timesteps to simulate in Tt
    
    ### Environment S
    # Get estimated prey consumed and predator growth given vital rates b, h, r, K, and c
    # We need to only remove the effect of high temperature
    tempS = eaten.odeint.general.pred(N0 = N0[SEL_environment == "S"],        # Prey start dens
                                      b = exp(bST.log + bT.log),
                                      h = hST + hT,
                                      r = rS,              # independent r's for SEL environments
                                      K = exp(KS.log),     # independent K's for SEL environments
                                      c = c,
                                      m = mST + mT,
                                      Tt = Tt[SEL_environment == "S"],        # Incubation time
                                      P = P[SEL_environment == "S"],          # Predator start dens
                                      steps = steps)  # Number of timesteps to simulate in Tt
    
    ### Environment T
    # Get estimated prey consumed and predator growth given vital rates b, h, r, K, and c
    # We need to only remove the effect of high salinity
    tempT = eaten.odeint.general.pred(N0 = N0[SEL_environment == "T"],        # Prey start dens
                                      b = exp(bST.log + bS.log),
                                      h = hST + hS,
                                      r = rT,              # independent r's for SEL environments
                                      K = exp(KT.log),     # independent K's for SEL environments
                                      c = c,
                                      m = mST + mS,
                                      Tt = Tt[SEL_environment == "T"],        # Incubation time
                                      P = P[SEL_environment == "T"],          # Predator start dens
                                      steps = steps)  # Number of timesteps to simulate in Tt
    
    ### Environment ST
    # Get estimated prey consumed and predator growth given vital rates b, h, r, K, and c
    tempST = eaten.odeint.general.pred(N0 = N0[SEL_environment == "ST"],        # Prey start dens
                                       b = exp(bST.log),
                                       h = hST,
                                       r = rST,              # independent r's for SEL environments
                                       K = exp(KST.log),     # independent K's for SEL environments
                                       c = c,
                                       m = mST,
                                       Tt = Tt[SEL_environment == "ST"],        # Incubation time
                                       P = P[SEL_environment == "ST"],          # Predator start dens
                                       steps = steps)  # Number of timesteps to simulate in Tt
    
    # Compute nll for each environment and sum them for general nll
    
    ### Environment C
    yC = tempC[[1]] # N.eaten estimate
    zC = tempC[[2]] # New predators in 1mL estimate
    
    nllC = -1*sum(dnorm(x = log(N0[SEL_environment == "C"] - Ndead[SEL_environment == "C"]), # Measured prey end density
                        mean = log(N0[SEL_environment == "C"] - yC),  # Estimated prey end density (Start prey dens + estimated change)
                        sd = sigma,
                        log = T))
    
    nll2C = -1*sum(dpois(x = P.end[SEL_environment == "C"],           # Measured predator end counts in 0.66 mL
                         lambda = (P[SEL_environment == "C"] + zC) * 0.66, # Estimated predator end counts in 0.66 mL (Start pred counts + estimated new preds)
                         log = T))
    
    ### Environment S
    yS = tempS[[1]] # N.eaten estimate
    zS = tempS[[2]] # New predators in 1mL estimate
    
    nllS = -1*sum(dnorm(x = log(N0[SEL_environment == "S"] - Ndead[SEL_environment == "S"]), # Measured prey end density
                        mean = log(N0[SEL_environment == "S"] - yS),  # Estimated prey end density (Start prey dens + estimated change)
                        sd = sigma,
                        log = T))
    
    nll2S = -1*sum(dpois(x = P.end[SEL_environment == "S"],           # Measured predator end counts in 0.66 mL
                         lambda = (P[SEL_environment == "S"] + zS) * 0.66, # Estimated predator end counts in 0.66 mL (Start pred counts + estimated new preds)
                         log = T))
    
    ### Environment T
    yT = tempT[[1]] # N.eaten estimate
    zT = tempT[[2]] # New predators in 1mL estimate
    
    nllT = -1*sum(dnorm(x = log(N0[SEL_environment == "T"] - Ndead[SEL_environment == "T"]), # Measured prey end density
                        mean = log(N0[SEL_environment == "T"] - yT),  # Estimated prey end density (Start prey dens + estimated change)
                        sd = sigma,
                        log = T))
    
    nll2T = -1*sum(dpois(x = P.end[SEL_environment == "T"],           # Measured predator end counts in 0.66 mL
                         lambda = (P[SEL_environment == "T"] + zT) * 0.66, # Estimated predator end counts in 0.66 mL (Start pred counts + estimated new preds)
                         log = T))
    
    
    ### Environment ST
    yST = tempST[[1]] # N.eaten estimate
    zST = tempST[[2]] # New predators in 1mL estimate
    
    nllST = -1*sum(dnorm(x = log(N0[SEL_environment == "ST"] - Ndead[SEL_environment == "ST"]), # Measured prey end density
                         mean = log(N0[SEL_environment == "ST"] - yST),  # Estimated prey end density (Start prey dens + estimated change)
                         sd = sigma,
                         log = T))
    
    nll2ST = -1*sum(dpois(x = P.end[SEL_environment == "ST"],           # Measured predator end counts in 0.66 mL
                          lambda = (P[SEL_environment == "ST"] + zST) * 0.66, # Estimated predator end counts in 0.66 mL (Start pred counts + estimated new preds)
                          log = T))
    
    # General nll value to minimize
    nll = nllC + nll2C + 
      nllS + nll2S + 
      nllT + nll2T + 
      nllST + nll2ST
    
    return(nll)
  }
}
