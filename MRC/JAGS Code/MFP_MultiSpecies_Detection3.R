####################################################################################
#
#                       MFP Playback - Camera Trap Data
#                       Detection Rate Bayesian Models
#                         Simplified Process Model
#
####################################################################################
library(R2jags); library(plyr); library(dplyr); library(ggplot2); library(abind)

load(file = "MFP_MultiSpecies_Detection_Data.rda", verbose = T) 


# MULTI-SPECIES VISITATION - Z Allowed to vary between treatments
sink("MFP_MultiSpecies_SimpProc3.jags")
cat("
    model{
    
    # Prior for gamma distribution from which rho is drawn
    theta <- exp(logtheta)
    logtheta ~ dunif(-5,5)
    
    ############
    # Priors for means and variances of community-level distributions from which parameters on lambda are drawn
    mu.a0 ~ dnorm(0, 0.001)
    mu.a1 ~ dnorm(0, 0.001)
    mu.a2 ~ dnorm(0, 0.001)
    mu.a3 ~ dnorm(0, 0.001)
    
    sig.a0 ~ dunif(0,100); tau.a0 <- 1/sig.a0
    sig.a1 ~ dunif(0,100); tau.a1 <- 1/sig.a1
    sig.a2 ~ dunif(0,100); tau.a2 <- 1/sig.a2
    sig.a3 ~ dunif(0,100); tau.a3 <- 1/sig.a3
    
    ############
    # Priors for means and variances of community-level distributions from which parameters on psi are drawn
    mu.b0 ~ dnorm(0, 0.001)
    mu.b1 ~ dnorm(0, 0.001)
    mu.b2 ~ dnorm(0, 0.001)
    mu.b3 ~ dnorm(0, 0.001)
    
    sig.b0 ~ dunif(0,100); tau.b0 <- 1/sig.b0
    sig.b1 ~ dunif(0,100); tau.b1 <- 1/sig.b1
    sig.b2 ~ dunif(0,100); tau.b2 <- 1/sig.b2
    sig.b3 ~ dunif(0,100); tau.b3 <- 1/sig.b3

    
    ############
    # Create separate parameters (for modeling lambda and psi) for each species 
    for(i in 1:Sp){                    # Loop through species (1 to 6)
    a0[i] ~ dnorm(mu.a0, tau.a0)     # Intercept
    a1[i] ~ dnorm(mu.a1, tau.a1)     # Treatment effect
    a2[i] ~ dnorm(mu.a2, tau.a2)     # Site effect 
    a3[i] ~ dnorm(mu.a3, tau.a3)     # Treatment x site interaction
    
    b0[i] ~ dnorm(mu.b0, tau.b0)     # Intercept
    b1[i] ~ dnorm(mu.b1, tau.b1)     # Treatment effect    
    b2[i] ~ dnorm(mu.b2, tau.b2)     # Site effect
    b3[i] ~ dnorm(mu.b3, tau.b3)     # Treatment x site interaction

    }
    
    ############
    t.track<-c(1,1,1,1,1,2,2,2,2,2)    # Vector to match each of the k survey periods to the appropriate treatment priod
    
    ############
    # Logistic regression submodel for site use
    for(i in 1:Sp){                    # Loop through species (1 to 6)
    for(j in 1:M){                   # Loop through camera sites (1 to 24)
    for(k in 1:K){                 # Loop through week (1 to 10)
    
    z[j,k,i] ~ dbern(psi[j,k,i])
    logit(psi[j,k,i]) <- b0[i] + b1[i]*Trt[j,k] + b2[i]*S[j,k] + b3[i]*Trt[j,k]*S[j,k]
    
    }
    }
    }
    
    ############
    # Negative binomial submodel for detection frequency
    for(i in 1:Sp){                    # Loop through species (1 to 6)
    for(j in 1:M){                   # Loop through camera sites (1 to 24)
    for(k in 1:K){                 # Loop through survey periods (1 to 10)
    
    log(lambda[j,k,i]) <- a0[i] + a1[i]*Trt[j,k] + a2[i]*S[j,k] + a3[i]*Trt[j,k]*S[j,k]
    
    rho[j,k,i] ~ dgamma(theta, theta)
    
    mu[j,k,i] <- rho[j,k,i]*lambda[j,k,i]
    
    y[j,k,i] ~ dpois(z[j,k,i] * mu[j,k,i])
    
    # Create new data for calculating Bayesian p-values
    y_new[j,k,i] ~ dpois(z[j,k,i] * mu[j,k,i])
    
    # Calculate statistics from both observed and new data sets for Bayesian p-values
    eval[j,k,i] <- z[j,k,i] * mu[j,k,i]
    
    # Freeman-Tukey Residual (Royle et al. 2013, Sect. 8.4.2)
    Terr[j,k,i] <- pow(pow(y[j,k,i],.5) - pow(eval[j,k,i],.5),2)
    Terrnew[j,k,i] <- pow(pow(y_new[j,k,i],.5) - pow(eval[j,k,i],.5),2)
    
    # Chi-squared stat (Kery & Schaub 2012, Sect. 12.3)
    ch.err[j,k,i] <- pow((y[j,k,i] - eval[j,k,i]),2)/ (eval[j,k,i] + 0.5)
    ch.errnew[j,k,i] <- pow((y_new[j,k,i] - eval[j,k,i]),2)/ (eval[j,k,i] + 0.5)
    
    }
    }
    }
    
    # Posterior predictive checks - Sum statistics above for calculating Bayesian p-values
    Tobs <- sum(Terr[,,])
    Tnew <- sum(Terrnew[,,])
    Chisq.obs <- sum(ch.err[,,])
    Chisq.new <- sum(ch.errnew[,,])
    
    ############
    # Estimate number of camera sites used as a derived parameter
    
    # First sum within experimental sites for both treatments
    for(i in 1:Sp){
    z.S2.c[i] <- mean(z[1:12,6:10,i])
    z.S2.h[i] <- mean(z[1:12,1:5,i])
    z.N2.c[i] <- mean(z[13:24,1:5,i])
    z.N2.h[i] <- mean(z[13:24,6:10,i])
    }
            
    # Then calculate average across sites for each treatment
    for(i in 1:Sp){
    z.ave.C[i] <- mean(c(z.S2.c[i], z.N2.c[i]))
    z.ave.H[i] <- mean(c(z.S2.h[i], z.N2.h[i]))
    }
    
    # Get average value of Rho
    rho.ave <- mean(rho[,,])
    
    # Estimate mean detection rate for elephants
    eleLam.N2.c <- mean(lambda[13:24,1:5,1])
    eleLam.N2.h <- mean(lambda[13:24,6:10,1])
    eleLam.S2.h <- mean(lambda[1:12,1:5,1])
    eleLam.S2.c <- mean(lambda[1:12,6:10,1])
    
    eleLam.C<-mean(c(eleLam.N2.c, eleLam.S2.c))
    eleLam.H<-mean(c(eleLam.N2.h, eleLam.S2.h))
    
    
    }
    ",fill = TRUE)
sink()

dat<-species.det.S.array[,,c(1,2)]
M<-dim(dat)[1]
K<-dim(dat)[2]
Sp<-dim(dat)[3]

data<-list(y = dat, K = K, M = M, Trt = T.hc, S = S.hc, Sp = Sp)


# Initial values for z
zst<-array(dim = c(M,K,Sp))
for(i in 1:Sp){
  for(j in 1:M){
    for(k in 1:K){
      # Initialize z based on whether any detections occurred in that week
      zst[j,k,i] <- ifelse(dat[j,k,i] > 0, 1, 0)
    }
  }
}


inits<-list(
  
  list(mu.a0 = runif(1,-2,2), mu.a1 = runif(1,-2,2), mu.a2 = runif(1,-2,2), mu.a3 = runif(1,-2,2), 
       mu.b0 = runif(1,-0.5,0.5), mu.b1 = runif(1,-0.5,0.5), mu.b2 = runif(1,-0.5,0.5), mu.b3 = runif(1,-0.5,0.5),
       sig.a0 = runif(1,1,5), sig.a1 = runif(1,1,5), sig.a2 = runif(1,1,5), sig.a3 = runif(1,1,5), 
       sig.b0 = runif(1,1,5), sig.b1 = runif(1,1,5), sig.b2 = runif(1,1,5), sig.b3 = runif(1,1,5), 
       z = zst),
  list(mu.a0 = runif(1,-2,2), mu.a1 = runif(1,-2,2), mu.a2 = runif(1,-2,2), mu.a3 = runif(1,-2,2), 
       mu.b0 = runif(1,-0.5,0.5), mu.b1 = runif(1,-0.5,0.5), mu.b2 = runif(1,-0.5,0.5), mu.b3 = runif(1,-0.5,0.5),
       sig.a0 = runif(1,1,5), sig.a1 = runif(1,1,5), sig.a2 = runif(1,1,5), sig.a3 = runif(1,1,5), 
       sig.b0 = runif(1,1,5), sig.b1 = runif(1,1,5), sig.b2 = runif(1,1,5), sig.b3 = runif(1,1,5), 
       z = zst),
  list(mu.a0 = runif(1,-2,2), mu.a1 = runif(1,-2,2), mu.a2 = runif(1,-2,2), mu.a3 = runif(1,-2,2), 
       mu.b0 = runif(1,-0.5,0.5), mu.b1 = runif(1,-0.5,0.5), mu.b2 = runif(1,-0.5,0.5), mu.b3 = runif(1,-0.5,0.5),
       sig.a0 = runif(1,1,5), sig.a1 = runif(1,1,5), sig.a2 = runif(1,1,5), sig.a3 = runif(1,1,5), 
       sig.b0 = runif(1,1,5), sig.b1 = runif(1,1,5), sig.b2 = runif(1,1,5), sig.b3 = runif(1,1,5), 
       z = zst)
  
)

# Parameters to monitor
parameters <- c('a0','a1','a2','a3','b0','b1','b2','b3','Tobs','Tnew','Chisq.obs','Chisq.new','z.S2.c','z.S2.h','z.N2.c','z.N2.h', 'rho.ave', 'eleLam.C','eleLam.H','eleLam.N2.c','eleLam.N2.h','eleLam.S2.c','eleLam.S2.h','z.ave.C','z.ave.H') 


ni=50000 # iterations
nb=10000 # burn-in
nthin=50 # thining
nc=3 # number of chains

det.out.wZT3<-jags(data = data, inits = inits, parameters.to.save = parameters, model.file = "MFP_MultiSpecies_SimpProc3.jags", n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nthin, progress.bar = "text")

###################################################################################################
###################################################################################################
# Save models

save(list = c("visit.out.wZT"), file = 'SCMP_Cam_VisitRate_FINAL_Model_20180912.rda')

###################################################################################################
###################################################################################################
# Posterior Sampling
samples <- det.out.wZT3$BUGSoutput$sims.list

# Percent change in detection frequency under human treatment vs control
ele_percent_change <- 100 * (samples$eleLam.H - samples$eleLam.C) / samples$eleLam.C

# Summarize
mean(ele_percent_change)
quantile(ele_percent_change, c(0.025, 0.5, 0.975))
