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
    # Prior for gamma distribution (rho)
    theta <- exp(logtheta)
    logtheta ~ dunif(-5,5)
    
    ############
    # Spatial Parameters for Gaussian Process
    # Range parameters
    rho.psi ~ dunif(0.1, 3)     
    rho.lambda ~ dunif(0.1, 3)  
            
    # Marginal variance parameters
    sigma2.psi ~ dunif(0.01, 3)     
    sigma2.lambda ~ dunif(0.01, 3)  
    
    ############
    # Priors for probability of occurrence parameters (psi)
    a0 ~ dnorm(0, 0.001)      #intercept
    a1 ~ dnorm(0, 0.001)      #site effect
    a2 ~ dnorm(0, 0.001)      #treatment effect
    a3 ~ dnorm(0, 0.001)      #site x treatment interaction
    
    ############
    # Priors for intensity of site use parameters (lambda)
    b0 ~ dnorm(0, 0.001)      #intercept
    b1 ~ dnorm(0, 0.001)      #site effect
    b2 ~ dnorm(0, 0.001)      #treatment effect
    b3 ~ dnorm(0, 0.001)      #site x treatment interaction
    
    ############
    # Covariance matrices for spatial effects
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
    #couples each k survey period to the appropriate treatment period
    t.track<-c(1,1,1,1,1,2,2,2,2,2)  
    
    ############
    # Logistic regression submodel for probability of occurence (psi)
    for(j in 1:M){                     #loop through camera traps (1 to 24)         
        for(t in 1:2){                 #loop through treatment periods (1 to 2)
            z[j,t] ~ dbern(psi[j,t])
            
            logit(psi[j,t]) <- b0 + b1*Trt[j,t] + b2*S.psi[j,t] + 
            b3*Trt[j,t]*S.psi[j,t] + phi.psi[j]
                
            loglik.psi[j,t] <- z[j,t] * log(psi[j,t] + 1.0E-10) + 
                               (1 - z[j,t]) * log(1 - psi[j,t] + 1.0E-10)
        }
    }
    
    ############
    # Negative binomial submodel for intensity of site use (lambda)
    for(j in 1:M){                     #loop through camera traps (1 to 24)
        for(k in 1:K){                 #loop through survey periods (1 to 10)
            log(lambda[j,k]) <- a0 + a1*Trt[j,k] + a2*S[j,k] + 
            a3*Trt[j,k]*S[j,k] + phi.lambda[j]
            
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
       
    # Total log-likelihood (sum of lambda and psi log-likelihoods)
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

    # Calculate average across sites for each treatment period
    z.sum.1 <- mean(c(z.sum.S2[1], z.sum.N2[1]))
    z.sum.2 <- mean(c(z.sum.S2[2], z.sum.N2[2]))
    
    # Get average value of Rho
    rho.ave <- mean(rho[,])
    
    # Estimate mean detection rate
    det.N2.c <- mean(lambda[13:24,1:5])
    det.N2.h <- mean(lambda[13:24,6:10])
    det.S2.h <- mean(lambda[1:12,1:5])
    det.S2.c <- mean(lambda[1:12,6:10])
    det.C <- mean(c(det.N2.c, det.S2.c))
    det.H <- mean(c(det.N2.h, det.S2.h))
    det.1 <- mean(c(det.N2.c, det.S2.h))
    det.2 <- mean(c(det.N2.h, det.S2.c))
    
    # Derived parameters for spatial correlation assessment
    eff.range.psi <- rho.psi * sqrt(-2 * log(0.05))
    eff.range.lambda <- rho.lambda * sqrt(-2 * log(0.05))
}",fill = TRUE)
sink()

# Prepare data for single species (select one species from your array)
# Assuming you want the first species (elephant) - change index as needed
species_index <- 1  # Change this to select different species
dat <- species.det.S.array[,,species_index]

M <- dim(dat)[1]  # Number of camera sites (24)
K <- dim(dat)[2]  # Number of survey periods (10)

# Data list for JAGS (no longer need Sp parameter)
data <- list(
  y = dat,
  K = K,
  M = M,
  Trt = T.hc,
  S = S.hc,
  S.psi = S.psi,
  D = D  #Scaled Distance matrix
)

# Initial values for z (occupancy states) - simplified for single species
zst <- array(dim = c(24, 2))
z1 <- dat[,1:5]; z1 <- apply(z1,1,sum); z1 <- ifelse(z1>0,1,0)
z2 <- dat[,6:10]; z2 <- apply(z2,1,sum, na.rm = T); z2 <- ifelse(z2>0,1,0)
zst <- cbind(z1,z2)

# Updated initial values for single species
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

# Parameters to monitor (simplified for single species)
parameters <- c("a0", "a1", "a2", "a3", "b0", "b1", "b2", "b3",
                "Tobs", "Tnew", "Chisq.obs", "Chisq.new",
                "det.N2.c", "det.N2.h", "det.S2.h", "det.S2.c",
                "z.sum.C","z.sum.H","det.C","det.H",
                "z.sum.1","z.sum.2","det.1","det.2",
                "rho.psi", "rho.lambda", "sigma2.psi", "sigma2.lambda",
                "eff.range.psi", "eff.range.lambda", 
                "loglik.psi", "loglik.lambda", "total.loglik")

ni <- 300000  # iterations
nb <- 250000  # burn-in
nthin <- 50   # thinning
nc <- 3       # number of chains

det.out.wZT6.ele.test <- jags(data = data, inits = inits, parameters.to.save = parameters, 
                       model.file = "MFP_SingleSpecies_SimpProc.jags", 
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nthin, 
                       progress.bar = "text")

###################################################################################################
###################################################################################################
# Save models

save(list = c("det.out.wZT6.ele.test"), file = 'MFP_DetRate_Model6ele.test.rda')

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

#Data Table Summaries
#African Elephant
# R2jags summary table (mean, sd, 2.5%, 97.5%, etc.)
sumtab <- det.out.wZT6.ele$BUGSoutput$summary

# convenience helper: grab mean + CrI for a parameter
get_mci <- function(par, sumtab){
  out <- sumtab[par, c("mean", "2.5%", "97.5%"), drop = FALSE]
  colnames(out) <- c("Mean", "Lower_CrI", "Upper_CrI")
  out
}

tab_s1 <- data.frame(
  Species = "African elephant",
  Control_Mean      = get_mci("z.sum.C", sumtab)[1, "Mean"],
  Control_Lower_CrI = get_mci("z.sum.C", sumtab)[1, "Lower_CrI"],
  Control_Upper_CrI = get_mci("z.sum.C", sumtab)[1, "Upper_CrI"],
  Human_Mean        = get_mci("z.sum.H", sumtab)[1, "Mean"],
  Human_Lower_CrI   = get_mci("z.sum.H", sumtab)[1, "Lower_CrI"],
  Human_Upper_CrI   = get_mci("z.sum.H", sumtab)[1, "Upper_CrI"]
)

tab_s1
# write.csv(tab_s1, "Table_S1_elephant.csv", row.names = FALSE)

psi_pars <- c("b0","b1","b2","b3")
psi_lab  <- c("Intercept (β0)", "Treatment (β1)", "Experimental site (β2)", "Treatment × site (β3)")

tab_s2a <- do.call(rbind, lapply(psi_pars, \(p) get_mci(p, sumtab)))
tab_s2a <- data.frame(
  Term = psi_lab,
  tab_s2a,
  row.names = NULL
)

# flag CrI that does not cross zero
tab_s2a$NotCross0 <- with(tab_s2a, Lower_CrI > 0 | Upper_CrI < 0)

tab_s2a

lam_pars <- c("a0","a1","a2","a3")
lam_lab  <- c("Intercept (α0)", "Treatment (α1)", "Experimental site (α2)", "Treatment × site (α3)")

tab_s2b <- do.call(rbind, lapply(lam_pars, \(p) get_mci(p, sumtab)))
tab_s2b <- data.frame(
  Term = lam_lab,
  tab_s2b,
  row.names = NULL
)

tab_s2b$NotCross0 <- with(tab_s2b, Lower_CrI > 0 | Upper_CrI < 0)

tab_s2b

fmt_ital <- function(x, italic = FALSE, digits = 2){
  s <- format(round(x, digits), nsmall = digits)
  ifelse(italic, paste0("*", s, "*"), s)
}

tab_s2a_fmt <- within(tab_s2a, {
  Mean      <- fmt_ital(Mean,      NotCross0)
  Lower_CrI <- fmt_ital(Lower_CrI, NotCross0)
  Upper_CrI <- fmt_ital(Upper_CrI, NotCross0)
  NotCross0 <- NULL
})

tab_s2b_fmt <- within(tab_s2b, {
  Mean      <- fmt_ital(Mean,      NotCross0)
  Lower_CrI <- fmt_ital(Lower_CrI, NotCross0)
  Upper_CrI <- fmt_ital(Upper_CrI, NotCross0)
  NotCross0 <- NULL
})

tab_s2a_fmt
tab_s2b_fmt


