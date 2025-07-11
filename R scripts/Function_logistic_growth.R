# ----------------------------------------------------------------------------------------------------------------------------------
# Logistic growth model for estimating prey r and K from the controls
# ----------------------------------------------------------------------------------------------------------------------------------
#
# The following R-code of the source functions is based on the method of:
#
# Rosenbaum, B. & Rall, B.C. (2018). Fitting functional responses: Direct parameter estimation by simulating differential equations. 
# Methods in Ecology and Evolution.
# 
# Lemmen, K. D. & Pennekamp, F. (2024). Food web context modifies predator foraging and weakens 
# trophic interaction strength.
#
# ----------------------------------------------------------------------------------------------------------------------------------

library(odeintr)

# Define ODEs for prey
# x[0] = prey density

LG.prey = '
dxdt[0] = r * x[0] * (1.0 - x[0]/K);
'

# Create LG_prey (integrated function), using the defined C++ expressions LG.prey
compile_sys("LG_prey", LG.prey, pars = c("r", "K"), method = "rk54")

# Create function to simulate prey densities in time, given the LG model ODEs
eq.odeint.LG.prey = function(start.v, r, K, Tt, timesteplength){
  LG_prey_set_params(r = r, K = K)
  LG_prey(start.v, Tt, timesteplength)
}

# Create function to get estimated prey growth:
#   Ngrowth.est: a vector of estimated number of new prey (prey growth)
#   returns a list with Ngrowth.est and the temp file where prey densities over time are saved

growth.odeint.LG.prey = function(N0, r, K, Tt, steps=100){
  
  Ngrowth.est = vector()
  
  for(i.growth in 1:length(N0)){
    # Simulate prey densities over 100 timesteps
    temp = eq.odeint.LG.prey(start.v = c(N0[i.growth]),
                             r = r,
                             K = K, 
                             Tt = Tt[i.growth],
                             timesteplength = Tt[i.growth]/steps)
    
    Ngrowth.est[i.growth] =  temp[steps+1,2] - N0[i.growth] # how many at the end - how many at the beginning = how many new
  }
  
  # Save estimated prey growth
  ret = list(Ngrowth.est, temp)
  
  return(ret)
  
}

nll.odeint.LG.prey = function(Ndiff, N0, r, K.log, steps = 100, sigma){
  if(sigma <= 0) return(Inf)
  
  # Get estimated new prey given vital rates r and K
  temp2 = growth.odeint.LG.prey(N0 = N0,        # Prey start dens
                                r = r,
                                K = exp(K.log),
                                Tt = Tt,        # Incubation time
                                steps = steps)  # Number of timesteps to simulate in Tt
  
  y = temp2[[1]] # N.growth estimate (how many new individuals)
  
  nll = -1*sum(dnorm(x = log(N0 + Ndiff), # Measured prey end density (prey start + how many new individuals)
                     mean = log(N0 + y),  # Estimated prey end density (Start prey dens + estimated change)
                     sd = sigma,
                     log = T))
  
  return(nll)
}
