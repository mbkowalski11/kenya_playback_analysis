####################################################################################
#
#                       MFP Playback - Camera Trap Data
#                       Detection Rate Bayesian Models
#                         Simplified Process Model
#
####################################################################################
library(R2jags); library(plyr); library(dplyr); library(ggplot2); library(abind); 
library(purrr)

load(file = "MFP_MultiSpecies_Detection_Data.rda", verbose = T) 
species_det_array_24x10x19 <- readRDS("species_det_array_24x10x19.rds")

# MULTI-SPECIES VISITATION WITH SPATIAL GAUSSIAN PROCESS
sink("MFP_MultiSpecies_SimpProc5.jags")
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
    mu.a0 ~ dnorm(0, 0.01)
    mu.a1 ~ dnorm(0, 0.01)
    mu.a2 ~ dnorm(0, 0.01)
    mu.a3 ~ dnorm(0, 0.01)
    
    sig.a0 ~ dunif(0,10); tau.a0 <- 1/(sig.a0*sig.a0)
    sig.a1 ~ dunif(0,10); tau.a1 <- 1/(sig.a1*sig.a1)
    sig.a2 ~ dunif(0,10); tau.a2 <- 1/(sig.a2*sig.a2)
    sig.a3 ~ dunif(0,10); tau.a3 <- 1/(sig.a3*sig.a3)
    
    ############
    # Priors for means and variances of community-level distributions
    mu.b0 ~ dnorm(0, 0.01)
    mu.b1 ~ dnorm(0, 0.01)
    mu.b2 ~ dnorm(0, 0.01)
    mu.b3 ~ dnorm(0, 0.01)
    
    sig.b0 ~ dunif(0,10); tau.b0 <- 1/(sig.b0*sig.b0)
    sig.b1 ~ dunif(0,10); tau.b1 <- 1/(sig.b1*sig.b1)
    sig.b2 ~ dunif(0,10); tau.b2 <- 1/(sig.b2*sig.b2)
    sig.b3 ~ dunif(0,10); tau.b3 <- 1/(sig.b3*sig.b3)
    
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
                loglik.lambda[j,k,i] <- logdensity.pois(y[j,k,i], z[j,t.track[k],i] * mu[j,k,i])

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
        total_loglik[i] <- sum(loglik.lambda[,,i])
    }
    
    # Posterior predictive checks (pooled)
    Tobs <- sum(Terr[,,])
    Tnew <- sum(Terrnew[,,])
    Chisq.obs <- sum(ch.err[,,])
    Chisq.new <- sum(ch.errnew[,,])

    # Posterior predictive checks (per-species)
    for(i in 1:Sp){
    Tobs_sp[i]   <- sum(Terr[,,i])
    Tnew_sp[i]   <- sum(Terrnew[,,i])
    Chisq_obs_sp[i] <- sum(ch.err[,,i])
    Chisq_new_sp[i] <- sum(ch.errnew[,,i])
    }
    
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
    for(i in 1:Sp){

      # Calculate average detections across sites for each treatment
      lamS2h[i] <- mean(lambda[1:12,  1:5,  i])   # human - S2
      lamS2c[i] <- mean(lambda[1:12,  6:10, i])   # control - S2
      lamN2h[i] <- mean(lambda[13:24, 6:10, i])   # human - N2
      lamN2c[i] <- mean(lambda[13:24, 1:5,  i])   # control - N2
    
      lam.mean.H[i] <- mean( c(lamS2h[i], lamN2h[i]) )
      lam.mean.C[i] <- mean( c(lamS2c[i], lamN2c[i]) )
    }
    
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
dat <- species_det_array_24x10x19
M <- dim(dat)[1]  # Number of camera sites (24)
K <- dim(dat)[2]  # Number of survey periods (10)
Sp <- dim(dat)[3] # Number of species (19)

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
parameters <- c(
  "a0","a1","a2","a3",
  "b0","b1","b2","b3",
  "Tobs","Tnew","Chisq.obs","Chisq.new",
  "Tobs_sp","Tnew_sp","Chisq_obs_sp","Chisq_new_sp",
  "z.sum.C","z.sum.H", "lam.mean.C","lam.mean.H",
  "rho.psi","rho.lambda","sigma2.psi","sigma2.lambda",
  "eff.range.psi","eff.range.lambda",
  "loglik.lambda","total_loglik")

ni=500000 # iterations
nb=400000 # burn-in
nthin=50 # thinning
nc=3 # number of chains

det.out.MultiSpec <- R2jags::jags(data = data, inits = inits, parameters.to.save = parameters, model.file = "MFP_MultiSpecies_SimpProc5.jags", n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nthin)

###################################################################################################
###################################################################################################
# Save models

save(list = c("det.out.MultiSpec"), file = 'MFP_DetRate_Model6.all.rda')

###################################################################################################
###################################################################################################

##########################stuff for multispecies and trophic groups

# your list:
# species_matrices  # named list of matrices

# the 19 species in the screenshot, in the order you want
desired <- c(
  "Hippopotamus","African Elephant","Common Warthog","Zorilla",
  "Reticulated Giraffe","Black-backed Jackal","Striped Hyena","Genets",
  "Spotted Hyena","Hares","Guenther's Dik-Dik","Leopard",
  "White-tailed Mongoose","Grevy's Zebra","Lion","Ground Squirrels",
  "Slender Mongoose","Impala","Plains Zebra"
)

# helper to coerce any found matrix to 24x10 numeric, else a 24x10 zero matrix
mk24x10 <- function(x, nm) {
  if (is.null(x)) {
    message("Missing species in list: ", nm, " — filling with zeros.")
    return(matrix(0, nrow = 24, ncol = 10))
  }
  x <- as.matrix(x)
  if (!all(dim(x) == c(24,10))) {
    stop(sprintf("Matrix for '%s' has dim %s, expected 24x10.", nm,
                 paste(dim(x), collapse = "x")))
  }
  storage.mode(x) <- "numeric"
  x[is.na(x)] <- 0
  x
}

# build a list of 24x10 matrices in the desired order
mats <- map(desired, ~ mk24x10(species_matrices[[.x]], .x))

# combine to a 24 x 10 x 19 array
arr <- abind(mats, along = 3)

# set nice dimnames
dimnames(arr) <- list(
  camera = as.character(1:24),
  period = as.character(1:10),
  species = desired
)

# quick checks
dim(arr)          # should be 24 10 19
dimnames(arr)$species

# optional: save for JAGS/NIMBLE/BUGS
saveRDS(arr, file = "species_det_array_24x10x19.rds")
# or .RData: save(arr, file = "species_array_24x10x19.RData")

############
#Trophic Groups

# Build a 3-D array along the 3rd dim from objects like "ele.det", "gir.det", ...
make_group_array <- function(var_names) {
  # pull the matrices from the workspace
  mats <- mget(var_names, inherits = TRUE)
  # quick checks
  stopifnot(all(sapply(mats, is.matrix)))
  stopifnot(length(unique(sapply(mats, nrow))) == 1)  # 24
  stopifnot(length(unique(sapply(mats, ncol))) == 1)  # 10
  
  arr <- abind(mats, along = 3)
  # give species codes as dimnames on 3rd dim
  dimnames(arr) <- list(NULL, NULL, sub("\\.det$", "", var_names))
  arr
}

#labels with sample sizes (n = total detections)
make_labels_with_n <- function(arr) {
  totals <- apply(arr, 3, sum, na.rm = TRUE)
  codes  <- dimnames(arr)[[3]]
  pretty <- c(
    hip = "Hippopotamus", gir = "Reticulated Giraffe", ele = "African Elephant",
    war = "Common Warthog", lep = "Hares", gdd = "Guenther's Dik-Dik",
    pze = "Plains Zebra", gvy = "Grevy's Zebra", imp = "Impala",
    gsq = "Ground Squirrels",
    lio = "Lion", par = "Leopard",
    hna = "Spotted Hyena", shy = "Striped Hyena", gen = "Genets",
    bbj = "Black-backed Jackal", wtm = "White-tailed Mongoose",
    smo = "Slender Mongoose", zor = "Zorilla"
  )
  setNames(
    paste0(pretty[codes], " (n=", totals, ")"),
    codes
  )
}

# Groups (order within each group is preserved)
megaherbivores <- c("hip.det", "gir.det", "ele.det")

herbivores <- c("war.det", "pze.det", "gvy.det", "imp.det")

smallherbivores <- c("lep.det", "gdd.det", "gsq.det")

large_carnivores <- c("lio.det", "par.det")

mesocarnivores <- c("hna.det", "shy.det", "gen.det", "bbj.det",
                    "wtm.det", "smo.det", "zor.det")

mega_array  <- make_group_array(megaherbivores)
herb_array  <- make_group_array(herbivores)
sherb_array <- make_group_array(smallherbivores) 
largec_array <- make_group_array(large_carnivores)
meso_array  <- make_group_array(mesocarnivores)

mega_labels  <- make_labels_with_n(mega_array)
herb_labels  <- make_labels_with_n(herb_array)
sherb_labels  <- make_labels_with_n(sherb_array)
largec_labels <- make_labels_with_n(largec_array)
meso_labels  <- make_labels_with_n(meso_array)

saveRDS(mega_array,  "array_megaherbivores.rds")
saveRDS(herb_array,  "array_herbivores.rds")
saveRDS(sherb_array,  "array_smallherbivores.rds")
saveRDS(largec_array,"array_large_carnivores.rds")
saveRDS(meso_array,  "array_mesocarnivores.rds")

## ---- prerequisites (already in your env) ----
## Required objects:
## - species_det_array_24x10x19  (24 x 10 x 19)
## - full_codes (length 19, in the array’s true order)
## - data (list passed to JAGS), inits (list of chains), parameters, nc, ni, nb, nthin
## - model file "MFP_MultiSpecies_SimpProc5.jags"

# Full species order used in the 24x10x19 array
full_codes <- c(
  "hip","ele","war","zor","gir","bbj","shy","gen","hna","lep",
  "gdd","par","wtm","gvy","lio","gsq","smo","imp","pze"
)

# Group definitions by short codes (order matters)
groups <- list(
  megaherbivores   = c("hip","gir","ele"),
  herbivores       = c("war","pze","gvy","imp"),
  smallherbivores  = c("lep","gdd","gsq"),
  large_carnivores = c("lio","par"),
  mesocarnivores   = c("hna","shy","gen","bbj","wtm","smo","zor")
)
# --- helpers -------------------------------------------

labels_with_n <- function(arr3, name_map = NULL){
  # arr3: 24 x 10 x Sp_g
  n_per_sp <- apply(arr3, 3, function(a) sum(a, na.rm = TRUE))
  codes <- dimnames(arr3)[[3]]
  if (!is.null(name_map)) {
    # use provided names when available, else fall back to codes
    species_names <- ifelse(codes %in% names(name_map), name_map[codes], codes)
  } else {
    species_names <- codes
  }
  paste0(species_names, " (n=", n_per_sp, ")")
}

run_group <- function(group_name, codes, Y3d, full_codes, data_base, inits, parameters,
                      model.file, nc, ni, nb, nthin, name_map = NULL,
                      save_rda = TRUE){
  # slice array
  idx <- match(codes, full_codes)
  if (any(is.na(idx))) stop("Unknown species code(s): ", paste(codes[is.na(idx)], collapse=", "))
  arr <- Y3d[, , idx, drop = FALSE]
  dimnames(arr) <- list(
    site   = dimnames(Y3d)[[1]],
    period = dimnames(Y3d)[[2]],
    species= codes
  )
  
  # trim inits per chain to these species
  trim_inits_to_idx <- function(init, keep_idx){
    if (!is.null(init$rho.psi))       init$rho.psi       <- init$rho.psi[keep_idx]
    if (!is.null(init$rho.lambda))    init$rho.lambda    <- init$rho.lambda[keep_idx]
    if (!is.null(init$sigma2.psi))    init$sigma2.psi    <- init$sigma2.psi[keep_idx]
    if (!is.null(init$sigma2.lambda)) init$sigma2.lambda <- init$sigma2.lambda[keep_idx]
    if (!is.null(init$z))             init$z <- init$z[ , , keep_idx, drop = FALSE]
    init
  }
  keep_idx <- idx
  inits_g  <- lapply(inits, trim_inits_to_idx, keep_idx = keep_idx)
  
  # assemble data list
  data_g <- data_base
  data_g$y  <- arr
  data_g$Sp <- dim(arr)[3]
  
  # run JAGS
  fit <- R2jags::jags(
    data = data_g,
    inits = inits_g,
    parameters.to.save = parameters,
    model.file = model.file,
    n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nthin
  )
  
  # save (optional)
  if (save_rda){
    save_name <- paste0("MFP_DetRate_Model6.", group_name, ".rda")
    obj_name  <- paste0("det.out.", group_name)
    assign(obj_name, fit, envir = .GlobalEnv)
    save(list = obj_name, file = save_name)
  }
  
  # build labels from codes (or mapped names if given)
  lbl <- labels_with_n(arr, name_map)
  list(fit = fit, array = arr, labels = lbl)
}

# ---- multi-group loop ------------------------------------------

trophic.group.out <- lapply(names(groups), function(g) {
  run_group(
    group_name = g,
    codes      = groups[[g]],
    Y3d        = species_det_array_24x10x19,
    full_codes = full_codes,
    data_base  = data,
    inits      = inits,
    parameters = parameters,
    model.file = "MFP_MultiSpecies_SimpProc5.jags",
    nc = nc, ni = ni, nb = nb, nthin = nthin,
    save_rda = TRUE
  )
})

names(trophic.group.out) <- names(groups)
save(trophic.group.out, file = "MFP_DetRate_Model6_TrophicGroupsP-G5.rda")

trophic.group.out$megaherbivores$fit
trophic.group.out$herbivores$fit
trophic.group.out$large_carnivores$fit
trophic.group.out$mesocarnivores$fit