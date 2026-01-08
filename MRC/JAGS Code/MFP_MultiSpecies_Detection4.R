####################################################################################
#
#                       MFP Playback - Camera Trap Data
#                       Detection Rate Bayesian Models
#                         Simplified Process Model
#
####################################################################################
library(R2jags); library(plyr); library(dplyr); library(ggplot2); library(abind)

load(file = "MFP_MultiSpecies_Detection_Data.rda", verbose = T) 

# MULTI-SPECIES VISITATION WITH SPATIAL GAUSSIAN PROCESS
sink("MFP_MultiSpecies_SimpProc4.jags")
cat("
model{
    # Prior for gamma distribution from which rho is drawn
    theta <- exp(logtheta)
    logtheta ~ dunif(-5,5)
    
    ############
    # Spatial Parameters for Gaussian Process (Species-specific)
    # Range parameters (controls spatial decay) - separate for each species
    for(i in 1:Sp){
        rho.psi[i] ~ dunif(0.1, 3)     # Constrain to avoid very small values
        rho.lambda[i] ~ dunif(0.1, 3)  # Constrain to avoid very small values
        
        # Marginal variance parameters - separate for each species
        sigma2.psi[i] ~ dunif(0.01, 3)     # Avoid zero variance
        sigma2.lambda[i] ~ dunif(0.01, 3)  # Avoid zero variance
    }
    
    ############
    # Priors for means and variances of community-level distributions
    mu.a0 ~ dnorm(0, 0.001)
    mu.a1 ~ dnorm(0, 0.001)
    mu.a2 ~ dnorm(0, 0.001)
    mu.a3 ~ dnorm(0, 0.001)
    
    sig.a0 ~ dunif(0,100); tau.a0 <- 1/(sig.a0*sig.a0)
    sig.a1 ~ dunif(0,100); tau.a1 <- 1/(sig.a1*sig.a1)
    sig.a2 ~ dunif(0,100); tau.a2 <- 1/(sig.a2*sig.a2)
    sig.a3 ~ dunif(0,100); tau.a3 <- 1/(sig.a3*sig.a3)
    
    ############
    # Priors for means and variances of community-level distributions
    mu.b0 ~ dnorm(0, 0.001)
    mu.b1 ~ dnorm(0, 0.001)
    mu.b2 ~ dnorm(0, 0.001)
    mu.b3 ~ dnorm(0, 0.001)
    
    sig.b0 ~ dunif(0,100); tau.b0 <- 1/(sig.b0*sig.b0)
    sig.b1 ~ dunif(0,100); tau.b1 <- 1/(sig.b1*sig.b1)
    sig.b2 ~ dunif(0,100); tau.b2 <- 1/(sig.b2*sig.b2)
    sig.b3 ~ dunif(0,100); tau.b3 <- 1/(sig.b3*sig.b3)
    
    ############
    # Create covariance matrices for spatial effects (Species-specific)
    for(i in 1:Sp){
    for(j in 1:M){
        for(k in 1:M){
            Sigma.psi[j,k,i]    <- sigma2.psi[i] *
                                   exp(-pow(D[j,k],2)/(2*pow(rho.psi[i],2))) +
                                   1.0E-6 * equals(j,k)

            Sigma.lambda[j,k,i] <- sigma2.lambda[i] *
                                   exp(-pow(D[j,k],2)/(2*pow(rho.lambda[i],2))) +
                                   1.0E-6 * equals(j,k)
        }
    }

    Omega.psi[1:M,1:M,i]    <- inverse(Sigma.psi[,,i])
    Omega.lambda[1:M,1:M,i] <- inverse(Sigma.lambda[,,i])
    }
    ############
    # Create separate parameters for each species 
    for(i in 1:Sp){
        a0[i] ~ dnorm(mu.a0, tau.a0)
        a1[i] ~ dnorm(mu.a1, tau.a1)
        a2[i] ~ dnorm(mu.a2, tau.a2)
        a3[i] ~ dnorm(mu.a3, tau.a3)
        
        b0[i] ~ dnorm(mu.b0, tau.b0)
        b1[i] ~ dnorm(mu.b1, tau.b1)
        b2[i] ~ dnorm(mu.b2, tau.b2)
        b3[i] ~ dnorm(mu.b3, tau.b3)
        
        # Spatial random effects for each species
        phi.psi[1:M,i] ~ dmnorm(zeros[], Omega.psi[,,i])
        phi.lambda[1:M,i] ~ dmnorm(zeros[], Omega.lambda[,,i])
    }
    
    # Vector of zeros for MVN mean
    for(i in 1:M){
        zeros[i] <- 0
    }
    
    ############
    t.track<-c(1,1,1,1,1,2,2,2,2,2)
    
    ############
    # Logistic regression submodel for site use
    for(i in 1:Sp){
        for(j in 1:M){
            for(t in 1:2){
                z[j,t,i] ~ dbern(psi[j,t,i])
                logit(psi[j,t,i]) <- b0[i] + b1[i]*Trt[j,t] + b2[i]*S.psi[j,t] + b3[i]*Trt[j,t]*S.psi[j,t] + phi.psi[j,i]
            }
        }
    }
    
    ############
    # Negative binomial submodel for detection frequency
    for(i in 1:Sp){
        for(j in 1:M){
            for(k in 1:K){
                log(lambda[j,k,i]) <- a0[i] + a1[i]*Trt[j,k] + a2[i]*S[j,k] + a3[i]*Trt[j,k]*S[j,k] + phi.lambda[j,i]
                rho[j,k,i] ~ dgamma(theta, theta)
                mu[j,k,i] <- rho[j,k,i]*lambda[j,k,i]
                y[j,k,i] ~ dpois(z[j,t.track[k],i] * mu[j,k,i])
                
                # Create new data for calculating Bayesian p-values
                y_new[j,k,i] ~ dpois(z[j,t.track[k],i] * mu[j,k,i])
                
                # Calculate statistics for Bayesian p-values
                eval[j,k,i] <- z[j,t.track[k],i] * mu[j,k,i]
                
                # Freeman-Tukey Residual
                Terr[j,k,i] <- pow(pow(y[j,k,i],.5) - pow(eval[j,k,i],.5),2)
                Terrnew[j,k,i] <- pow(pow(y_new[j,k,i],.5) - pow(eval[j,k,i],.5),2)
                
                # Chi-squared stat
                ch.err[j,k,i] <- pow((y[j,k,i] - eval[j,k,i]),2)/ (eval[j,k,i] + 0.5)
                ch.errnew[j,k,i] <- pow((y_new[j,k,i] - eval[j,k,i]),2)/ (eval[j,k,i] + 0.5)
            }
        }
    }
    
    # Posterior predictive checks
    Tobs <- sum(Terr[,,])
    Tnew <- sum(Terrnew[,,])
    Chisq.obs <- sum(ch.err[,,])
    Chisq.new <- sum(ch.errnew[,,])
    
    ############
    # Derived parameters
    
    # Sum within experimental sites for both treatments
    for(t in 1:2){
        for(i in 1:Sp){
            z.sum.S2[t,i] <- sum(z[1:12,t,i])
            z.sum.N2[t,i] <- sum(z[13:24,t,i])
        }
    }
    
    # Calculate average across sites for each treatment
    for(i in 1:Sp){
        z.sum.C[i] <- mean(c(z.sum.S2[2,i], z.sum.N2[1,i]))
        z.sum.H[i] <- mean(c(z.sum.S2[1,i], z.sum.N2[2,i]))
    }
    
    # Get average value of Rho
    rho.ave <- mean(rho[,,])
    
    # Estimate mean detection rate for African Elephant
    eleLam.N2.c <- mean(lambda[13:24,1:5,1])
    eleLam.N2.h <- mean(lambda[13:24,6:10,1])
    eleLam.S2.h <- mean(lambda[1:12,1:5,1])
    eleLam.S2.c <- mean(lambda[1:12,6:10,1])
    eleLam.C<-mean(c(eleLam.N2.c, eleLam.S2.c))
    eleLam.H<-mean(c(eleLam.N2.h, eleLam.S2.h))
    
    # Estimate mean detection rate for Reticulated Giraffe
    girLam.N2.c <- mean(lambda[13:24,1:5,2])
    girLam.N2.h <- mean(lambda[13:24,6:10,2])
    girLam.S2.h <- mean(lambda[1:12,1:5,2])
    girLam.S2.c <- mean(lambda[1:12,6:10,2])
    girLam.C<-mean(c(girLam.N2.c, girLam.S2.c))
    girLam.H<-mean(c(girLam.N2.h, girLam.S2.h))

    # Estimate mean detection rate for Grevy's Zebra
    gvyLam.N2.c <- mean(lambda[13:24,1:5,3])
    gvyLam.N2.h <- mean(lambda[13:24,6:10,3])
    gvyLam.S2.h <- mean(lambda[1:12,1:5,3])
    gvyLam.S2.c <- mean(lambda[1:12,6:10,3])
    gvyLam.C<-mean(c(gvyLam.N2.c, gvyLam.S2.c))
    gvyLam.H<-mean(c(gvyLam.N2.h, gvyLam.S2.h))

    # Estimate mean detection rate for Black-Backed Jackal
    bbjLam.N2.c <- mean(lambda[13:24,1:5,4])
    bbjLam.N2.h <- mean(lambda[13:24,6:10,4])
    bbjLam.S2.h <- mean(lambda[1:12,1:5,4])
    bbjLam.S2.c <- mean(lambda[1:12,6:10,4])
    bbjLam.C<-mean(c(bbjLam.N2.c, bbjLam.S2.c))
    bbjLam.H<-mean(c(bbjLam.N2.h, bbjLam.S2.h))

    # Estimate mean detection rate for Guenther's Dik-Dik
    gddLam.N2.c <- mean(lambda[13:24,1:5,5])
    gddLam.N2.h <- mean(lambda[13:24,6:10,5])
    gddLam.S2.h <- mean(lambda[1:12,1:5,5])
    gddLam.S2.c <- mean(lambda[1:12,6:10,5])
    gddLam.C<-mean(c(gddLam.N2.c, gddLam.S2.c))
    gddLam.H<-mean(c(gddLam.N2.h, gddLam.S2.h))

    # Estimate mean detection rate for Impala
    impLam.N2.c <- mean(lambda[13:24,1:5,6])
    impLam.N2.h <- mean(lambda[13:24,6:10,6])
    impLam.S2.h <- mean(lambda[1:12,1:5,6])
    impLam.S2.c <- mean(lambda[1:12,6:10,6])
    impLam.C<-mean(c(impLam.N2.c, impLam.S2.c))
    impLam.H<-mean(c(impLam.N2.h, impLam.S2.h))

    # Estimate mean detection rate for Striped Hyena
    shyLam.N2.c <- mean(lambda[13:24,1:5,7])
    shyLam.N2.h <- mean(lambda[13:24,6:10,7])
    shyLam.S2.h <- mean(lambda[1:12,1:5,7])
    shyLam.S2.c <- mean(lambda[1:12,6:10,7])
    shyLam.C<-mean(c(shyLam.N2.c, shyLam.S2.c))
    shyLam.H<-mean(c(shyLam.N2.h, shyLam.S2.h))

    # Estimate mean detection rate for Spotted Hyena
    hnaLam.N2.c <- mean(lambda[13:24,1:5,8])
    hnaLam.N2.h <- mean(lambda[13:24,6:10,8])
    hnaLam.S2.h <- mean(lambda[1:12,1:5,8])
    hnaLam.S2.c <- mean(lambda[1:12,6:10,8])
    hnaLam.C<-mean(c(hnaLam.N2.c, hnaLam.S2.c))
    hnaLam.H<-mean(c(hnaLam.N2.h, hnaLam.S2.h))

    # Estimate mean detection rate for Hares
    lepLam.N2.c <- mean(lambda[13:24,1:5,9])
    lepLam.N2.h <- mean(lambda[13:24,6:10,9])
    lepLam.S2.h <- mean(lambda[1:12,1:5,9])
    lepLam.S2.c <- mean(lambda[1:12,6:10,9])
    lepLam.C<-mean(c(lepLam.N2.c, lepLam.S2.c))
    lepLam.H<-mean(c(lepLam.N2.h, lepLam.S2.h))

    # Estimate mean detection rate for Slender Mongoose
    smoLam.N2.c <- mean(lambda[13:24,1:5,10])
    smoLam.N2.h <- mean(lambda[13:24,6:10,10])
    smoLam.S2.h <- mean(lambda[1:12,1:5,10])
    smoLam.S2.c <- mean(lambda[1:12,6:10,10])
    smoLam.C<-mean(c(smoLam.N2.c, smoLam.S2.c))
    smoLam.H<-mean(c(smoLam.N2.h, smoLam.S2.h))
    
    # Estimate mean detection rate for White-tailed mongoose
    wtmLam.N2.c <- mean(lambda[13:24,1:5,11])
    wtmLam.N2.h <- mean(lambda[13:24,6:10,11])
    wtmLam.S2.h <- mean(lambda[1:12,1:5,11])
    wtmLam.S2.c <- mean(lambda[1:12,6:10,11])
    wtmLam.C<-mean(c(wtmLam.N2.c, wtmLam.S2.c))
    wtmLam.H<-mean(c(wtmLam.N2.h, wtmLam.S2.h))
     
    # Estimate mean detection rate for Genets
    genLam.N2.c <- mean(lambda[13:24,1:5,12])
    genLam.N2.h <- mean(lambda[13:24,6:10,12])
    genLam.S2.h <- mean(lambda[1:12,1:5,12])
    genLam.S2.c <- mean(lambda[1:12,6:10,12])
    genLam.C<-mean(c(genLam.N2.c, genLam.S2.c))
    genLam.H<-mean(c(genLam.N2.h, genLam.S2.h))

    # Derived parameters for spatial correlation assessment
    for(i in 1:Sp){
        eff.range.psi[i] <- rho.psi[i] * sqrt(-2 * log(0.05))
        eff.range.lambda[i] <- rho.lambda[i] * sqrt(-2 * log(0.05))
    }
}
",fill = TRUE)
sink()

#species array order: ele.det, gir.det, gvy.det, bbj.det, gdd.det, imp.det, shy.det, hna.det, lep.det, smo.det, wtm.det, gen.det

# Prepare data (assuming your existing objects are available)
dat <- species.det.S.array[,,c(1,2,3,4,5,6,7,8,9)]
M <- dim(dat)[1]  # Number of camera sites (24)
K <- dim(dat)[2]  # Number of survey periods (10)
Sp <- dim(dat)[3] # Number of species (2)

# Data list for JAGS (including distance matrix)
data <- list(
  y = dat, 
  K = K, 
  M = M, 
  Trt = T.hc, 
  S = S.hc, 
  S.psi = S.psi, 
  Sp = Sp,
  D = D  # Distance matrix
)

# Initial values for z (occupancy states)
zst <- array(dim = c(24, 2, Sp))
for(i in 1:Sp){
  z1 <- dat[,1:5,i]; z1 <- apply(z1,1,sum); z1 <- ifelse(z1>0,1,0)
  z2 <- dat[,6:10,i]; z2 <- apply(z2,1,sum, na.rm = T); z2 <- ifelse(z2>0,1,0)
  zz <- cbind(z1,z2)
  zst[,,i] <- zz
}

# Updated initial values with better constraints
inits <- list(
  list(
    mu.a0 = runif(1,-1,1), mu.a1 = runif(1,-1,1), mu.a2 = runif(1,-1,1), mu.a3 = runif(1,-1,1),
    mu.b0 = runif(1,-0.5,0.5), mu.b1 = runif(1,-0.5,0.5), mu.b2 = runif(1,-0.5,0.5), mu.b3 = runif(1,-0.5,0.5),
    sig.a0 = runif(1,0.5,2), sig.a1 = runif(1,0.5,2), sig.a2 = runif(1,0.5,2), sig.a3 = runif(1,0.5,2),
    sig.b0 = runif(1,0.5,2), sig.b1 = runif(1,0.5,2), sig.b2 = runif(1,0.5,2), sig.b3 = runif(1,0.5,2),
    rho.psi = runif(Sp, 0.5, 1.2), rho.lambda = runif(Sp, 0.5, 1.2),
    sigma2.psi = runif(Sp, 0.4, 1.5), sigma2.lambda = runif(Sp, 0.4, 1.5),
    logtheta = runif(1, -2, 2),
    z = zst
  ),
  list(
    mu.a0 = runif(1,-1,1), mu.a1 = runif(1,-1,1), mu.a2 = runif(1,-1,1), mu.a3 = runif(1,-1,1),
    mu.b0 = runif(1,-0.5,0.5), mu.b1 = runif(1,-0.5,0.5), mu.b2 = runif(1,-0.5,0.5), mu.b3 = runif(1,-0.5,0.5),
    sig.a0 = runif(1,0.5,2), sig.a1 = runif(1,0.5,2), sig.a2 = runif(1,0.5,2), sig.a3 = runif(1,0.5,2),
    sig.b0 = runif(1,0.5,2), sig.b1 = runif(1,0.5,2), sig.b2 = runif(1,0.5,2), sig.b3 = runif(1,0.5,2),
    rho.psi = runif(Sp, 0.5, 1.2), rho.lambda = runif(Sp, 0.5, 1.2),
    sigma2.psi = runif(Sp, 0.4, 1.5), sigma2.lambda = runif(Sp, 0.4, 1.5),
    logtheta = runif(1, -2, 2),
    z = zst
  ),
  list(
    mu.a0 = runif(1,-1,1), mu.a1 = runif(1,-1,1), mu.a2 = runif(1,-1,1), mu.a3 = runif(1,-1,1),
    mu.b0 = runif(1,-0.5,0.5), mu.b1 = runif(1,-0.5,0.5), mu.b2 = runif(1,-0.5,0.5), mu.b3 = runif(1,-0.5,0.5),
    sig.a0 = runif(1,0.5,2), sig.a1 = runif(1,0.5,2), sig.a2 = runif(1,0.5,2), sig.a3 = runif(1,0.5,2),
    sig.b0 = runif(1,0.5,2), sig.b1 = runif(1,0.5,2), sig.b2 = runif(1,0.5,2), sig.b3 = runif(1,0.5,2),
    rho.psi = runif(Sp, 0.5, 1.2), rho.lambda = runif(Sp, 0.5, 1.2),
    sigma2.psi = runif(Sp, 0.4, 1.5), sigma2.lambda = runif(Sp, 0.4, 1.5),
    logtheta = runif(1, -2, 2),
    z = zst
  )
)

#paramters to monitor
parameters <- c("a0", "a1", "a2", "a3","b0","b1","b2","b3","Tobs","Tnew","Chisq.obs","Chisq.new","z.sum.C","z.sum.H","eleLam.C","eleLam.H","girLam.C","girLam.H","gvyLam.C","gvyLam.H","bbjLam.C","bbjLam.H","gddLam.C","gddLam.H","impLam.C","impLam.H","shyLam.C","shyLam.H","hnaLam.C","hnaLam.H","lepLam.C","lepLam.H","smoLam.C","smoLam.H","wtmLam.C","wtmLam.H","genLam.C","genLam.H","rho.psi","rho.lambda","sigma2.psi","sigma2.lambda","eff.range.psi","eff.range.lambda")

ni=300000 # iterations
nb=250000 # burn-in
nthin=50 # thinning
nc=3 # number of chains

det.out.wZT4<-jags(data = data, inits = inits, parameters.to.save = parameters, model.file = "MFP_MultiSpecies_SimpProc4.jags", n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nthin, progress.bar = "text")

###################################################################################################
###################################################################################################
# Save models

save(list = c("visit.out.wZT"), file = 'SCMP_Cam_VisitRate_FINAL_Model_20180912.rda')

###################################################################################################
###################################################################################################
# Posterior Sampling
samples <- det.out.wZT4$BUGSoutput$sims.list

# Percent change in detection frequency under human treatment vs control
ele_percent_change <- 100 * (samples$eleLam.H - samples$eleLam.C) / samples$eleLam.C

# Summarize
mean(ele_percent_change)
quantile(ele_percent_change, c(0.025, 0.5, 0.975))
