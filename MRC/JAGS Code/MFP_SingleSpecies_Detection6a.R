####################################################################################
#
#                       MFP Playback - Camera Trap Data
#                       Detection Rate Bayesian Models
#                         Simplified Process Model
#
####################################################################################
library(R2jags); library(plyr); library(dplyr); library(ggplot2); library(abind)

load(file = "MFP_MultiSpecies_Detection_Data.rda", verbose = T) 

# SINGLE-SPECIES VISITATION WITH SPATIAL GAUSSIAN PROCESS
sink("MFP_SingleSpecies_SimpProc.jags")
cat("model{    
    # Prior for gamma distribution from which rho is drawn
    theta <- exp(logtheta)
    logtheta ~ dunif(-5,5)
    
    ############
    # Spatial Parameters for Gaussian Process
    # Range parameters (controls spatial decay)
    rho.psi ~ dunif(0.1, 3)     # Constrain to avoid very small values
    rho.lambda ~ dunif(0.1, 3)  # Constrain to avoid very small values
    
    # Marginal variance parameters
    sigma2.psi ~ dunif(0.01, 3)     # Avoid zero variance
    sigma2.lambda ~ dunif(0.01, 3)  # Avoid zero variance
    
    ############
    # Priors for occupancy parameters
    a0 ~ dnorm(0, 0.001)
    a1 ~ dnorm(0, 0.001)
    a2 ~ dnorm(0, 0.001)
    a3 ~ dnorm(0, 0.001)
    
    ############
    # Priors for detection parameters
    b0 ~ dnorm(0, 0.001)
    b1 ~ dnorm(0, 0.001)
    b2 ~ dnorm(0, 0.001)
    b3 ~ dnorm(0, 0.001)
    
    ############
    # Create covariance matrices for spatial effects
    for(j in 1:M){
        for(k in 1:M){
            Sigma.psi[j,k] <- sigma2.psi *
                              exp(-pow(D[j,k],2)/(2*pow(rho.psi,2))) +
                              1.0E-6 * equals(j,k)
            Sigma.lambda[j,k] <- sigma2.lambda *
                                exp(-pow(D[j,k],2)/(2*pow(rho.lambda,2))) +
                                1.0E-6 * equals(j,k)
        }
    }
    Omega.psi[1:M,1:M] <- inverse(Sigma.psi[,])
    Omega.lambda[1:M,1:M] <- inverse(Sigma.lambda[,])
    
    ############
    # Spatial random effects
    phi.psi[1:M] ~ dmnorm(zeros[], Omega.psi[,])
    phi.lambda[1:M] ~ dmnorm(zeros[], Omega.lambda[,])
    
    # Vector of zeros for MVN mean
    for(i in 1:M){
        zeros[i] <- 0
    }
    
    ############
    t.track<-c(1,1,1,1,1,2,2,2,2,2)
    
    ############
    # Logistic regression submodel for site use
    for(j in 1:M){
        for(t in 1:2){
            z[j,t] ~ dbern(psi[j,t])
            logit(psi[j,t]) <- b0 + b1*Trt[j,t] + b2*S.psi[j,t] + b3*Trt[j,t]*S.psi[j,t] + phi.psi[j]
            
            loglik.psi[j,t] <- z[j,t] * log(psi[j,t] + 1.0E-10) +
                                (1 - z[j,t]) * log(1 - psi[j,t] + 1.0E-10)
        }
    }
    
    ############
    # Negative binomial submodel for detection frequency
    for(j in 1:M){
        for(k in 1:K){
            log(lambda[j,k]) <- a0 + a1*Trt[j,k] + a2*S[j,k] + a3*Trt[j,k]*S[j,k] + phi.lambda[j]
            rho[j,k] ~ dgamma(theta, theta)
            mu[j,k] <- rho[j,k]*lambda[j,k]
            y[j,k] ~ dpois(z[j,t.track[k]] * mu[j,k])
            loglik.lambda[j,k] <- logdensity.pois(y[j,k], z[j,t.track[k]] * mu[j,k])
            
            # Create new data for calculating Bayesian p-values
            y_new[j,k] ~ dpois(z[j,t.track[k]] * mu[j,k])
            
            # Calculate statistics for Bayesian p-values
            eval[j,k] <- z[j,t.track[k]] * mu[j,k]
            
            # Freeman-Tukey Residual
            Terr[j,k] <- pow(pow(y[j,k],.5) - pow(eval[j,k],.5),2)
            Terrnew[j,k] <- pow(pow(y_new[j,k],.5) - pow(eval[j,k],.5),2)
            
            # Chi-squared stat
            ch.err[j,k] <- pow((y[j,k] - eval[j,k]),2)/ (eval[j,k] + 0.5)
            ch.errnew[j,k] <- pow((y_new[j,k] - eval[j,k]),2)/ (eval[j,k] + 0.5)
        }
    }
    
    # Posterior predictive checks
    Tobs <- sum(Terr[,])
    Tnew <- sum(Terrnew[,])
    Chisq.obs <- sum(ch.err[,])
    Chisq.new <- sum(ch.errnew[,])
    
    # WAIC calculations
    # Total log-likelihood (sum of lambda and psi components)
    total.loglik <- sum(loglik.psi[,]) + sum(loglik.lambda[,])
    
    ############
    # Derived parameters
    
    # Sum within experimental sites for both treatments
    for(t in 1:2){
        z.sum.S2[t] <- sum(z[1:12,t])
        z.sum.N2[t] <- sum(z[13:24,t])
    }
    
    # Calculate average across sites for each treatment
    z.sum.C <- mean(c(z.sum.S2[2], z.sum.N2[1]))
    z.sum.H <- mean(c(z.sum.S2[1], z.sum.N2[2]))
    # Calculate average across sites for each season
    z.sum.W <- mean(c(z.sum.S2[1], z.sum.N2[1]))
    z.sum.D <- mean(c(z.sum.S2[2], z.sum.N2[2]))
    
    # Get average value of Rho
    rho.ave <- mean(rho[,])
    
    # Estimate mean detection rate for the species
    det.N2.c <- mean(lambda[13:24,1:5])
    det.N2.h <- mean(lambda[13:24,6:10])
    det.S2.h <- mean(lambda[1:12,1:5])
    det.S2.c <- mean(lambda[1:12,6:10])
    det.C <- mean(c(det.N2.c, det.S2.c))
    det.H <- mean(c(det.N2.h, det.S2.h))
    det.W <- mean(c(det.N2.c, det.S2.h))
    det.D <- mean(c(det.N2.h, det.S2.c))
    
    # Derived parameters for spatial correlation assessment
    eff.range.psi <- rho.psi * sqrt(-2 * log(0.05))
    eff.range.lambda <- rho.lambda * sqrt(-2 * log(0.05))
}",fill = TRUE)
sink()

# Define species codes in the order they appear in the array
species_codes <- c("ele", "gir", "gvy", "bbj", "gdd", "imp", "shy", "hna", 
                   "lep", "smo", "wtm", "gen", "pze", "hip", "zor", "war", 
                   "lio", "par", "gsq")

# Get array dimensions
M <- dim(species.det.array)[1]  # Number of camera sites (24)
K <- dim(species.det.array)[2]  # Number of survey periods (10)

# Parameters to monitor
parameters <- c("a0", "a1", "a2", "a3", "b0", "b1", "b2", "b3",
                "Tobs", "Tnew", "Chisq.obs", "Chisq.new",
                "det.N2.c", "det.N2.h", "det.S2.h", "det.S2.c",
                "z.sum.C","z.sum.H","det.C","det.H",
                "z.sum.W","z.sum.D","det.W","det.D",
                "rho.psi", "rho.lambda", "sigma2.psi", "sigma2.lambda",
                "eff.range.psi", "eff.range.lambda",
                "loglik.psi", "loglik.lambda", "total.loglik")

# MCMC settings
ni <- 300000  # iterations
nb <- 250000  # burn-in
nthin <- 50   # thinning
nc <- 3       # number of chains

# Loop through all 19 species
for(species_index in 1:19) {
  
  cat("Processing species", species_index, "of 19:", species_codes[species_index], "\n")
  
  # Extract data for current species
  dat <- species.det.array[,,species_index]
  
  # Data list for JAGS
  data <- list(
    y = dat,
    K = K,
    M = M,
    Trt = T.hc,
    S = S.hc,
    S.psi = S.psi,
    D = D  # Scaled Distance matrix
  )
  
  # Initial values for z (occupancy states)
  zst <- array(dim = c(24, 2))
  z1 <- dat[,1:5]; z1 <- apply(z1,1,sum); z1 <- ifelse(z1>0,1,0)
  z2 <- dat[,6:10]; z2 <- apply(z2,1,sum, na.rm = T); z2 <- ifelse(z2>0,1,0)
  zst <- cbind(z1,z2)
  
  # Initial values for current species
  inits <- list(
    list(
      a0 = runif(1,-1,1), a1 = runif(1,-1,1), a2 = runif(1,-1,1), a3 = runif(1,-1,1),
      b0 = runif(1,-0.5,0.5), b1 = runif(1,-0.5,0.5), b2 = runif(1,-0.5,0.5), b3 = runif(1,-0.5,0.5),
      rho.psi = runif(1, 0.5, 1.2), rho.lambda = runif(1, 0.5, 1.2),
      sigma2.psi = runif(1, 0.4, 1.5), sigma2.lambda = runif(1, 0.4, 1.5),
      logtheta = runif(1, -2, 2),
      z = zst
    ),
    list(
      a0 = runif(1,-1,1), a1 = runif(1,-1,1), a2 = runif(1,-1,1), a3 = runif(1,-1,1),
      b0 = runif(1,-0.5,0.5), b1 = runif(1,-0.5,0.5), b2 = runif(1,-0.5,0.5), b3 = runif(1,-0.5,0.5),
      rho.psi = runif(1, 0.5, 1.2), rho.lambda = runif(1, 0.5, 1.2),
      sigma2.psi = runif(1, 0.4, 1.5), sigma2.lambda = runif(1, 0.4, 1.5),
      logtheta = runif(1, -2, 2),
      z = zst
    ),
    list(
      a0 = runif(1,-1,1), a1 = runif(1,-1,1), a2 = runif(1,-1,1), a3 = runif(1,-1,1),
      b0 = runif(1,-0.5,0.5), b1 = runif(1,-0.5,0.5), b2 = runif(1,-0.5,0.5), b3 = runif(1,-0.5,0.5),
      rho.psi = runif(1, 0.5, 1.2), rho.lambda = runif(1, 0.5, 1.2),
      sigma2.psi = runif(1, 0.4, 1.5), sigma2.lambda = runif(1, 0.4, 1.5),
      logtheta = runif(1, -2, 2),
      z = zst
    )
  )
  
  # Run JAGS model
  model_result <- jags(data = data, inits = inits, parameters.to.save = parameters,
                       model.file = "MFP_SingleSpecies_SimpProc.jags",
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nthin,
                       progress.bar = "text")
  
  # Create variable name for current species
  var_name <- paste0("det.out.wZT6.", species_codes[species_index])
  
  # Assign result to variable with species-specific name
  assign(var_name, model_result)
  
  # Save model with species-specific filename
  save_filename <- paste0('MFP_DetRate_Model6', species_codes[species_index], '.rda')
  save(list = var_name, file = save_filename)
  
  cat("Completed and saved:", save_filename, "\n\n")
}

cat("All 19 species models completed!\n")

#ele.det, gir.det, gvy.det, bbj.det, gdd.det, imp.det,shy.det, hna.det, 
#lep.det, smo.det, wtm.det, gen.det,pze.det, hip.det, zor.det, war.det,
#lio.det, par.det, gsq.det

###################################################################################################
###################################################################################################
# Posterior Sampling of Bayesian P-Values and Percent Change
#African Elephant
samples <- det.out.wZT6.ele.test$BUGSoutput$sims.list
p_value_T <- mean(samples$Tnew >= samples$Tobs); p_value_T
p_value_Chisq <- mean(samples$Chisq.new >= samples$Chisq.obs); p_value_Chisq
ele_percent_change_abun <- 100 * (samples$det.H - samples$det.C) / samples$det.C
mean(ele_percent_change_abun)
quantile(ele_percent_change_abun, c(0.025, 0.5, 0.975))
ele_percent_change_occ <- 100 * (samples$z.sum.H - samples$z.sum.C) / samples$z.sum.C
mean(ele_percent_change_occ)
quantile(ele_percent_change_occ, c(0.025, 0.5, 0.975))
waic_results <- calculate_waic(det.out.wZT6.ele)
print(paste("WAIC:", round(waic_results$waic, 2)))
print(paste("Standard Error:", round(waic_results$se_waic, 2)))
print(paste("Effective number of parameters:", round(waic_results$p_waic, 2)))

#Reticulated Giraffe
samples <- det.out.wZT6.gir$BUGSoutput$sims.list
p_value_T <- mean(samples$Tnew >= samples$Tobs); p_value_T
p_value_Chisq <- mean(samples$Chisq.new >= samples$Chisq.obs); p_value_Chisq
gir_percent_change_abun <- 100 * (samples$det.H - samples$det.C) / samples$det.C
mean(gir_percent_change_abun)
quantile(gir_percent_change_abun, c(0.025, 0.5, 0.975))
gir_percent_change_occ <- 100 * (samples$z.sum.H - samples$z.sum.C) / samples$z.sum.C
mean(gir_percent_change_occ)
quantile(gir_percent_change_occ, c(0.025, 0.5, 0.975))
waic_results <- calculate_waic(det.out.wZT6.gir)
print(paste("WAIC:", round(waic_results$waic, 2)))
print(paste("Standard Error:", round(waic_results$se_waic, 2)))
print(paste("Effective number of parameters:", round(waic_results$p_waic, 2)))


