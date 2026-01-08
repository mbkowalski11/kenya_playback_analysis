####################################################################################
#
#                       MFP Playback - Camera Trap Data
#                       Detection Rate Bayesian Models
#                         Model Output Calculations
#
####################################################################################
library(tidyverse);library(purrr); library(extrafont)

##########
#Load Model Objects - both gaussian
load("MFP_DetRate_Model6ele.rda")   
load("MFP_DetRate_Model6gir.rda")
load("MFP_DetRate_Model6gvy.rda")   
load("MFP_DetRate_Model6bbj.rda")
load("MFP_DetRate_Model6gdd.rda")   
load("MFP_DetRate_Model6imp.rda")
load("MFP_DetRate_Model6shy.rda")   
load("MFP_DetRate_Model6hna.rda")
load("MFP_DetRate_Model6lep.rda")   
load("MFP_DetRate_Model6smo.rda")
load("MFP_DetRate_Model6wtm.rda")   
load("MFP_DetRate_Model6gen.rda")
load("MFP_DetRate_Model6pze.rda")   
load("MFP_DetRate_Model6hip.rda")
load("MFP_DetRate_Model6zor.rda")   
load("MFP_DetRate_Model6war.rda")
load("MFP_DetRate_Model6lio.rda")
load("MFP_DetRate_Model6par.rda")   
load("MFP_DetRate_Model6gsq.rda")

model_list <- list(
  ele = det.out.wZT6.ele,
  gir = det.out.wZT6.gir,
  gvy = det.out.wZT6.gvy,
  bbj = det.out.wZT6.bbj,
  gdd = det.out.wZT6.gdd,
  imp = det.out.wZT6.imp,
  shy = det.out.wZT6.shy,
  hna = det.out.wZT6.hna,
  lep = det.out.wZT6.lep,
  smo = det.out.wZT6.smo,
  wtm = det.out.wZT6.wtm,
  gen = det.out.wZT6.gen,
  pze = det.out.wZT6.pze,
  hip = det.out.wZT6.hip,
  zor = det.out.wZT6.zor,
  war = det.out.wZT6.war,
  lio = det.out.wZT6.lio,
  par = det.out.wZT6.par,
  gsq = det.out.wZT6.gsq
  )

#non-spatial
load("MFP_DetRate_Model7ele.rda")   
load("MFP_DetRate_Model7gir.rda")
load("MFP_DetRate_Model7gvy.rda")   
load("MFP_DetRate_Model7bbj.rda")
load("MFP_DetRate_Model7gdd.rda")   
load("MFP_DetRate_Model7imp.rda")
load("MFP_DetRate_Model7shy.rda")   
load("MFP_DetRate_Model7hna.rda")
load("MFP_DetRate_Model7lep.rda")   
load("MFP_DetRate_Model7smo.rda")
load("MFP_DetRate_Model7wtm.rda")   
load("MFP_DetRate_Model7gen.rda")
load("MFP_DetRate_Model7pze.rda")   
load("MFP_DetRate_Model7hip.rda")
load("MFP_DetRate_Model7zor.rda")   
load("MFP_DetRate_Model7war.rda")
load("MFP_DetRate_Model7lio.rda")
load("MFP_DetRate_Model7par.rda")   
load("MFP_DetRate_Model7gsq.rda")

model_list <- list(
  ele = det.out.wZT7.ele,
  gir = det.out.wZT7.gir,
  gvy = det.out.wZT7.gvy,
  bbj = det.out.wZT7.bbj,
  gdd = det.out.wZT7.gdd,
  imp = det.out.wZT7.imp,
  shy = det.out.wZT7.shy,
  hna = det.out.wZT7.hna,
  lep = det.out.wZT7.lep,
  smo = det.out.wZT7.smo,
  wtm = det.out.wZT7.wtm,
  gen = det.out.wZT7.gen,
  pze = det.out.wZT7.pze,
  hip = det.out.wZT7.hip,
  zor = det.out.wZT7.zor,
  war = det.out.wZT7.war,
  lio = det.out.wZT7.lio,
  par = det.out.wZT7.par,
  gsq = det.out.wZT7.gsq
)

#both exponential
# Load Model 9 objects
load("MFP_DetRate_Model9ele.rda")   
load("MFP_DetRate_Model9gir.rda")
load("MFP_DetRate_Model9gvy.rda")   
load("MFP_DetRate_Model9bbj.rda")
load("MFP_DetRate_Model9gdd.rda")   
load("MFP_DetRate_Model9imp.rda")
load("MFP_DetRate_Model9shy.rda")   
load("MFP_DetRate_Model9hna.rda")
load("MFP_DetRate_Model9lep.rda")   
load("MFP_DetRate_Model9smo.rda")
load("MFP_DetRate_Model9wtm.rda")   
load("MFP_DetRate_Model9gen.rda")
load("MFP_DetRate_Model9pze.rda")   
load("MFP_DetRate_Model9hip.rda")
load("MFP_DetRate_Model9zor.rda")   
load("MFP_DetRate_Model9war.rda")
load("MFP_DetRate_Model9lio.rda")
load("MFP_DetRate_Model9par.rda")   
load("MFP_DetRate_Model9gsq.rda")

# Assemble into list
model_list <- list(
  ele = det.out.wZT9.ele,
  gir = det.out.wZT9.gir,
  gvy = det.out.wZT9.gvy,
  bbj = det.out.wZT9.bbj,
  gdd = det.out.wZT9.gdd,
  imp = det.out.wZT9.imp,
  shy = det.out.wZT9.shy,
  hna = det.out.wZT9.hna,
  lep = det.out.wZT9.lep,
  smo = det.out.wZT9.smo,
  wtm = det.out.wZT9.wtm,
  gen = det.out.wZT9.gen,
  pze = det.out.wZT9.pze,
  hip = det.out.wZT9.hip,
  zor = det.out.wZT9.zor,
  war = det.out.wZT9.war,
  lio = det.out.wZT9.lio,
  par = det.out.wZT9.par,
  gsq = det.out.wZT9.gsq
)

load("MFP_DetRate_Model6ele.test.rda") 
load("MFP_DetRate_Model7ele.rda") 
load("MFP_DetRate_Model8ele.rda") 
load("MFP_DetRate_Model9ele.rda") 
load("MFP_DetRate_Model10ele.rda") 

model_list <- list(
  ele6 = det.out.wZT6.ele.test, #gausian both
  ele7 = det.out.wZT7.ele, #non spatial
  ele8 = det.out.wZT8.ele, #gausian lambda
  ele9 = det.out.wZT9.ele, #exponential both
  ele10 = det.out.wZT10.ele #xponential lambda
)

##########
# Convert Each Species Model Outputs into a Single Row 
extract_metrics <- function(mod) {
  sims <- mod$BUGSoutput$sims.list
  
  # Bayesian P-Values
  p_T      <- mean(sims$Tnew       >= sims$Tobs)
  p_Chisq  <- mean(sims$Chisq.new  >= sims$Chisq.obs)
  
  # Effect Sizes
  a1 <- sims$a1   
  a2 <- sims$a2    
  b1 <- sims$b1    
  b2 <- sims$b2     
  
  # Percent Change  
  pc <- function(num, den) 100 * (num - den) / den
  
  #Tibble
  tibble(
    #Bayesian P-Values
    p_T                = p_T,
    p_Chisq            = p_Chisq,
    
    #Model Effect Sizes
    T_eff_abun_mean    = mean(a1),
    T_eff_abun_lo      = quantile(a1, 0.025),
    T_eff_abun_med     = quantile(a1, 0.500),
    T_eff_abun_hi      = quantile(a1, 0.975),
    
    S_eff_abun_mean    = mean(a2),
    S_eff_abun_lo      = quantile(a2, 0.025),
    S_eff_abun_med     = quantile(a2, 0.500),
    S_eff_abun_hi      = quantile(a2, 0.975),
    
    T_eff_occ_mean     = mean(b1),
    T_eff_occ_lo       = quantile(b1, 0.025),
    T_eff_occ_med      = quantile(b1, 0.500),
    T_eff_occ_hi       = quantile(b1, 0.975),
    
    S_eff_occ_mean     = mean(b2),
    S_eff_occ_lo       = quantile(b2, 0.025),
    S_eff_occ_med      = quantile(b2, 0.500),
    S_eff_occ_hi       = quantile(b2, 0.975),
    
    #Percent Change Derived Values
    pc_abun_T_mean     = mean(pc(sims$det.H,     sims$det.C)),
    pc_abun_T_lo       = quantile(pc(sims$det.H, sims$det.C), 0.025),
    pc_abun_T_med      = quantile(pc(sims$det.H, sims$det.C), 0.500),
    pc_abun_T_hi       = quantile(pc(sims$det.H, sims$det.C), 0.975),
    
    pc_occ_T_mean      = mean(pc(sims$z.sum.H,   sims$z.sum.C)),
    pc_occ_T_lo        = quantile(pc(sims$z.sum.H, sims$z.sum.C), 0.025),
    pc_occ_T_med       = quantile(pc(sims$z.sum.H, sims$z.sum.C), 0.500),
    pc_occ_T_hi        = quantile(pc(sims$z.sum.H, sims$z.sum.C), 0.975),
    
    pc_abun_S_mean     = mean(pc(sims$det.D,     sims$det.W)),
    pc_abun_S_lo       = quantile(pc(sims$det.D, sims$det.W), 0.025),
    pc_abun_S_med      = quantile(pc(sims$det.D, sims$det.W), 0.500),
    pc_abun_S_hi       = quantile(pc(sims$det.D, sims$det.W), 0.975),
    
    pc_occ_S_mean      = mean(pc(sims$z.sum.D,   sims$z.sum.W)),
    pc_occ_S_lo        = quantile(pc(sims$z.sum.D, sims$z.sum.W), 0.025),
    pc_occ_S_med       = quantile(pc(sims$z.sum.D, sims$z.sum.W), 0.500),
    pc_occ_S_hi        = quantile(pc(sims$z.sum.D, sims$z.sum.W), 0.975),
    
    #Model Selection
    WAIC               = calculate_waic(mod)$waic,
    WAIC_se            = calculate_waic(mod)$se_waic,
    WAIC_p             = calculate_waic(mod)$p_waic
  )
}

# Bind Rows for Each Species
species_summary <- imap_dfr(model_list, ~extract_metrics(.x) %>% 
                            mutate(species = .y), .id = NULL) %>% 
                            relocate(species)
extract_metrics()
#Write to CSV
print(species_summary, width = Inf)
write.csv(species_summary, "model6.outputs.csv", row.names = FALSE)

library(purrr)
library(dplyr)
library(tibble)

get_rhat_all <- function(model_list) {
  imap_dfr(model_list, function(mod, species) {
    S <- as.data.frame(mod$BUGSoutput$summary)
    S$param <- rownames(mod$BUGSoutput$summary)
    keep <- Reduce(`|`, lapply(
      c("^a(\\d+)?$", "^b(\\d+)?$", "^alpha(\\[\\d+\\])?$", "^beta(\\[\\d+\\])?$"),
      grepl, x = S$param
    ))
    S <- S[keep, , drop = FALSE]
    tibble(
      species = species,
      param   = S$param,
      Rhat    = S$Rhat,
      n_eff   = S[["n.eff"]],
      mean    = S$mean,
      sd      = S$sd,
      lo95    = S[["2.5%"]],
      med     = S[["50%"]],
      hi95    = S[["97.5%"]]
    )
  })
}

rhat_table <- get_rhat_all(model_list)
view(rhat_table)
##########
#Summary Plots
df <- read.csv("model6.outputs.csv")

#Species Labels
species_map <- c(
  ele = "African Elephant (n=242)",
  gir = "Reticulated Giraffe (n=127)",
  gvy = "Grevy's Zebra (n=389)",
  bbj = "Black-backed Jackal (n=219)",
  gdd = "Gunther's Dik-Dik (n=5547)",
  imp = "Impala (n=2153)",
  shy = "Striped Hyena (n=224)",
  hna = "Spotted Hyena (n=765)",
  lep = "Hares (n=2310)",
  smo = "Slender Mongoose (n=150)",
  wtm = "White-tailed Mongoose (n=426)",
  gen = "Genets (n=311)",
  pze = "Plains Zebra (n=68)",
  hip = "Hippopotamus (n=54)",
  zor = "Zorilla (n=63)",
  war = "Common Warthog (n=60)",
  gsq = "Ground Squirrels (n=377)",
  par = "Leopard (n=35)",
  lio = "Lion (n=46)"
)

df$species_full <- species_map[as.character(df$species)]
df$species_full[is.na(df$species_full)] <- as.character(df$species)

#Percent Change Plots
df <- df %>% # Order species by treatment abundance median (top = most negative)
  arrange(pc_abun_T_med) %>%
  mutate(species_full = factor(species_full, levels = rev(unique(species_full))))

plot_data <- df %>%
  select(species_full,
         abun_T_med = pc_abun_T_med, abun_T_lo = pc_abun_T_lo, abun_T_hi = pc_abun_T_hi,
         occ_T_med = pc_occ_T_med, occ_T_lo = pc_occ_T_lo, occ_T_hi = pc_occ_T_hi,
         abun_S_med = pc_abun_S_med, abun_S_lo = pc_abun_S_lo, abun_S_hi = pc_abun_S_hi,
         occ_S_med = pc_occ_S_med, occ_S_lo = pc_occ_S_lo, occ_S_hi = pc_occ_S_hi) %>%
  pivot_longer(
    cols = -species_full,
    names_to = c("metric", "effect", ".value"),
    names_pattern = "([a-z]+)_([TS])_(med|lo|hi)"
  ) %>%
  mutate(
    metric = recode(metric, abun = "Abundance", occ = "Occupancy"),
    effect = recode(effect, T = "Treatment", S = "Season")
  )

plot_data$effect <- factor(plot_data$effect, levels = c("Treatment", "Season")) # Fix facet order

ggplot(plot_data, aes(x = med, y = species_full, color = metric)) +
  geom_point(position = position_dodge(width = 0.6)) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.3, position = position_dodge(width = 0.6)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~effect, scales = "fixed") +
  scale_color_manual(values = c("Abundance" = "#7A0403FF", "Occupancy" = "darkgreen")) +
  labs(
    x = "Percent Change (%)",
    y = NULL,
    color = NULL,
    title = "Posterior Median and 95% CI of Percent Change",
    subtitle = "Effects of Treatment and Season on Occupancy and Abundance"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

#Effect Size Plots
df <- df %>%
  mutate(code = tolower(trimws(species)),
         species_full = coalesce(species_map[code], species)) %>%
  arrange(T_eff_abun_med) %>%
  mutate(species_full = factor(species_full, levels = rev(unique(species_full))))

plot_data <- df %>%
  select(species_full,
         T_eff_abun_med, T_eff_abun_lo, T_eff_abun_hi,
         S_eff_abun_med, S_eff_abun_lo, S_eff_abun_hi) %>%
  pivot_longer(
    -species_full,
    names_to = c("effect", ".value"),
    names_pattern = "([TS])_eff_abun_(med|lo|hi)"
  ) %>%
  mutate(effect = recode(effect, T = "Treatment", S = "Season"))

# Force correct facet order
plot_data$effect <- factor(plot_data$effect, levels = c("Treatment", "Season"))

ggplot(plot_data, aes(x = med, y = species_full)) +
  geom_point(color = "blue", position = position_dodge(width = 0.6)) +
  geom_errorbarh(aes(xmin = lo, xmax = hi),
                 height = 0.3,
                 position = position_dodge(width = 0.6),
                 color = "blue") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~effect, scales = "fixed", labeller = label_value) + 
  labs(x = "Abundance Effect Size (Link Scale)",
       y = NULL,
       title = "Posterior Median ± 95% CI of Abundance Coefficients") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

#Detection Frequency Plot
plot_data <- df %>%
  select(species_full,
         abun_T_med = pc_abun_T_med, abun_T_lo = pc_abun_T_lo, abun_T_hi = pc_abun_T_hi) %>%
  pivot_longer(
    cols = -species_full,
    names_to = c("metric", "effect", ".value"),
    names_pattern = "([a-z]+)_([TS])_(med|lo|hi)"
  )

ggplot(plot_data, aes(x = med, y = species_full, color = metric)) +
  geom_point(position = position_dodge(width = 0.6), size = 24) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, position = position_dodge(width = 0.6), size = 8) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 5) +
  scale_color_manual(values = "#7A0403FF", guide = "none") +
  labs(
    x = "Percent Change (%)",
    y = NULL,
    color = NULL,
  ) +
  theme_minimal() +
  theme(text = element_text(family = "Palatino", size = 180),
        axis.text = element_text(size = 180, color = "black"),
        axis.text.y = element_text(size = 180, color = "black", margin = margin(r = 250)),
        axis.text.x = element_text(size = 180, color = "black"),
        axis.title = element_text(size = 180),
        axis.title.x = element_text(margin = margin(t = 150)),
        plot.subtitle = element_text(size = 180),
        panel.grid = element_line(color = "grey60")
  )

####################################################################################
# Multi-species percent-change summaries + plots from jagsUI object det.out.MultiSpec
####################################################################################
library(tidyverse)

# ---- 1) species names in the SAME order as your 3rd array dim (Sp = 19) ----
species_names <- c(
  "Hippopotamus","African Elephant","Common Warthog","Zorilla",
  "Reticulated Giraffe","Black-backed Jackal","Striped Hyena","Genet",
  "Spotted Hyena","Hares","Guenther's Dik-Dik","Leopard",
  "White-tailed Mongoose","Grevy's Zebra","Lion","Ground Squirrels",
  "Slender Mongoose","Impala","Plains Zebra"
)

# Optional: pretty labels with (n=) counts for plotting (edit counts if needed)
species_labels <- c(
  "Hippopotamus (n=54)","African Elephant (n=242)","Common Warthog (n=60)","Zorilla (n=63)",
  "Reticulated Giraffe (n=127)","Black-backed Jackal (n=219)","Striped Hyena (n=224)","Genets (n=311)",
  "Spotted Hyena (n=765)","Hares (n=2310)","Guenther's Dik-Dik (n=5547)","Leopard (n=35)",
  "White-tailed Mongoose (n=426)","Grevy's Zebra (n=389)","Lion (n=46)","Ground Squirrels (n=377)",
  "Slender Mongoose (n=150)","Impala (n=2153)","Plains Zebra (n=68)"
)

stopifnot(length(species_names) == length(species_labels))

# ---- 2) pull posterior draws ----
sims <- det.out.MultiSpec$BUGSoutput$sims.list

# Expected elements from your multi-species model:
# a1[,i], a2[,i]          : detection (abundance) effects
# b1[,i], b2[,i]          : occupancy effects
# lam.mean.H[,i], lam.mean.C[,i] : treatment-specific mean detection rate per species
# z.sum.H[,i],  z.sum.C[,i]      : treatment-specific occupied-site totals per species
# (If names differ, adjust the lines below accordingly.)

eps <- 1e-10  # guard against divide-by-zero

make_pc <- function(num, den) 100 * (num - den) / (den + eps)

Sp <- length(species_names)
draws_to_summary <- function(x) {
  c(mean = mean(x), lo = quantile(x, .025), med = median(x), hi = quantile(x, .975))
}

# Build a tidy table of summaries per species
summ_list <- vector("list", Sp)
for(i in seq_len(Sp)){
  # percent change: detection (λ means) and occupancy (z sums)
  pc_abun_T <- make_pc(sims$lam.mean.H[,i], sims$lam.mean.C[,i])
  pc_occ_T  <- make_pc(sims$z.sum.H[,i],  sims$z.sum.C[,i])
  
  # coefficient summaries (link scale)
  a1_sum <- draws_to_summary(sims$a1[,i])
  a2_sum <- draws_to_summary(sims$a2[,i])
  b1_sum <- draws_to_summary(sims$b1[,i])
  b2_sum <- draws_to_summary(sims$b2[,i])
  
  summ_list[[i]] <- tibble(
    species              = species_names[i],
    # % change (Treatment: Human vs Control)
    pc_abun_T_mean = mean(pc_abun_T), pc_abun_T_lo = quantile(pc_abun_T,.025),
    pc_abun_T_med  = median(pc_abun_T), pc_abun_T_hi = quantile(pc_abun_T,.975),
    pc_occ_T_mean  = mean(pc_occ_T),  pc_occ_T_lo  = quantile(pc_occ_T,.025),
    pc_occ_T_med   = median(pc_occ_T), pc_occ_T_hi  = quantile(pc_occ_T,.975),
    
    # effect sizes on link scale
    T_eff_abun_mean = a1_sum["mean"], T_eff_abun_lo = a1_sum["lo"],
    T_eff_abun_med  = a1_sum["med"],  T_eff_abun_hi = a1_sum["hi"],
    
    S_eff_abun_mean = a2_sum["mean"], S_eff_abun_lo = a2_sum["lo"],
    S_eff_abun_med  = a2_sum["med"],  S_eff_abun_hi = a2_sum["hi"],
    
    T_eff_occ_mean  = b1_sum["mean"], T_eff_occ_lo  = b1_sum["lo"],
    T_eff_occ_med   = b1_sum["med"],  T_eff_occ_hi  = b1_sum["hi"],
    
    S_eff_occ_mean  = b2_sum["mean"], S_eff_occ_lo  = b2_sum["lo"],
    S_eff_occ_med   = b2_sum["med"],  S_eff_occ_hi  = b2_sum["hi"]
  )
}
species_summary <- bind_rows(summ_list) |>
  mutate(species_full = species_labels[match(species, species_names)])

# Save/inspect
print(species_summary, width = Inf)
write_csv(species_summary, "multispecies_percent_change_summary.csv")

############ The one that works for multispecies!

summ <- function(x) c(mean = mean(x), lo = quantile(x,.025), med = median(x), hi = quantile(x,.975))
pct_change <- function(num, den, eps=1e-10) 100 * (num - den) / (den + eps)

# unify access to posterior draws (matrix with columns like a1[1], lam.mean.H[3], etc.)
# --- safer converters ---------------------------------------------------------
as_posterior_matrix <- function(fit){
  # jagsUI path
  if (!is.null(fit$samples)) {
    mats <- lapply(fit$samples, function(x) {
      m <- as.matrix(x)
      # keep as-is; jagsUI already has proper colnames
      m
    })
    return(do.call(rbind, mats))
  }
  
  # R2jags/bugs path
  sl <- fit$BUGSoutput$sims.list
  if (is.null(sl)) stop("No samples found in object.")
  
  mats <- list()
  for (nm in names(sl)) {
    arr <- sl[[nm]]
    if (is.null(arr)) next
    
    # vector (one value per draw)
    if (is.atomic(arr) && is.null(dim(arr))) {
      it <- length(arr)
      out <- matrix(arr, ncol = 1)
      colnames(out) <- nm
      mats[[nm]] <- out
      next
    }
    
    # matrix (iter x J)
    if (is.matrix(arr)) {
      out <- arr
      if (is.null(colnames(out))) {
        colnames(out) <- paste0(nm, "[", seq_len(ncol(out)), "]")
      } else {
        colnames(out) <- paste0(nm, "[", colnames(out), "]")
      }
      mats[[nm]] <- out
      next
    }
    
    # array (iter x i x j x k ...)
    if (length(dim(arr)) >= 3) {
      iters <- dim(arr)[1]
      idx   <- as.matrix(expand.grid(lapply(dim(arr)[-1], seq_len)))
      cols  <- vector("list", nrow(idx))
      coln  <- character(nrow(idx))
      for (r in seq_len(nrow(idx))) {
        sel <- cbind(seq_len(iters),
                     matrix(idx[r,], ncol = ncol(idx), nrow = iters, byrow = TRUE))
        v <- arr[sel]
        cols[[r]] <- v
        coln[r]   <- paste0(nm, "[", paste(idx[r,], collapse=","), "]")
      }
      out <- do.call(cbind, cols)
      colnames(out) <- coln
      mats[[nm]] <- out
      next
    }
  }
  
  do.call(cbind, mats)
}

get_family <- function(post, prefix){
  idx <- grep(paste0("^", prefix, "\\["), colnames(post))
  if (!length(idx)) stop("No columns found for '", prefix, "'. Did you monitor it?")
  fam <- post[, idx, drop = FALSE]
  ord <- as.numeric(sub(".*\\[(\\d+)\\].*", "\\1", colnames(fam)))
  fam[, order(ord), drop = FALSE]
}

# ---- compute the table from your fitted object ----
make_multispecies_summary <- function(fit, species_names = NULL){
  post <- as_posterior_matrix(fit)
  
  A1 <- get_family(post, "a1"); A2 <- get_family(post, "a2")
  B1 <- get_family(post, "b1"); B2 <- get_family(post, "b2")
  H  <- get_family(post, "lam.mean.H"); C <- get_family(post, "lam.mean.C")
  ZH <- get_family(post, "z.sum.H");    ZC <- get_family(post, "z.sum.C")
  
  Sp <- ncol(A1)
  if (is.null(species_names)) species_names <- paste0("sp", seq_len(Sp))
  
  rows <- vector("list", Sp)
  for(i in seq_len(Sp)){
    pc_abun_T <- pct_change(H[,i],  C[,i])
    pc_occ_T  <- pct_change(ZH[,i], ZC[,i])
    
    a1s <- summ(A1[,i]); a2s <- summ(A2[,i])
    b1s <- summ(B1[,i]); b2s <- summ(B2[,i])
    
    rows[[i]] <- tibble(
      species = species_names[i],
      pc_abun_T_mean = mean(pc_abun_T), pc_abun_T_lo = quantile(pc_abun_T,.025),
      pc_abun_T_med  = median(pc_abun_T), pc_abun_T_hi = quantile(pc_abun_T,.975),
      pc_occ_T_mean  = mean(pc_occ_T),  pc_occ_T_lo  = quantile(pc_occ_T,.025),
      pc_occ_T_med   = median(pc_occ_T), pc_occ_T_hi  = quantile(pc_occ_T,.975),
      
      T_eff_abun_mean = a1s["mean"], T_eff_abun_lo = a1s["lo"],
      T_eff_abun_med  = a1s["med"],  T_eff_abun_hi = a1s["hi"],
      
      S_eff_abun_mean = a2s["mean"], S_eff_abun_lo = a2s["lo"],
      S_eff_abun_med  = a2s["med"],  S_eff_abun_hi = a2s["hi"],
      
      T_eff_occ_mean  = b1s["mean"], T_eff_occ_lo  = b1s["lo"],
      T_eff_occ_med   = b1s["med"],  T_eff_occ_hi  = b1s["hi"],
      
      S_eff_occ_mean  = b2s["mean"], S_eff_occ_lo  = b2s["lo"],
      S_eff_occ_med   = b2s["med"],  S_eff_occ_hi  = b2s["hi"]
    )
  }
  bind_rows(rows)
}

# ---- run it on your fit ----
species_labels <- c(
  "Hippopotamus (n=54)",
  "African Elephant (n=242)",
  "Common Warthog (n=60)",
  "Zorilla (n=63)",
  "Reticulated Giraffe (n=127)",
  "Black-backed Jackal (n=219)",
  "Striped Hyena (n=224)",
  "Genets (n=311)",
  "Spotted Hyena (n=765)",
  "Hares (n=2310)",
  "Guenther's Dik-Dik (n=5547)",
  "Leopard (n=35)",
  "White-tailed Mongoose (n=426)",
  "Grevy's Zebra (n=389)",
  "Lion (n=46)",
  "Ground Squirrels (n=377)",
  "Slender Mongoose (n=150)",
  "Impala (n=2153)",
  "Plains Zebra (n=68)"
)
ms_summary <- make_multispecies_summary(det.out.MultiSpec, species_labels)

print(ms_summary, width = Inf)
# write.csv(ms_summary, "multispecies_percent_change_summary.csv", row.names = FALSE)


library(forcats)

df <- ms_summary %>%
  mutate(
    species_full = if (exists("species_labels")) {
      coalesce(unname(species_labels[species]), species)
    } else species
  ) %>%
  arrange(pc_abun_T_med) %>%
  mutate(
    # order by median percent change, most negative at top
    species_full = fct_reorder(species_full, pc_abun_T_med),
    species_full = fct_rev(species_full)
  )

ggplot(df, aes(x = pc_abun_T_med, y = species_full)) +
  geom_point(color = "#7A0403") +
  geom_errorbarh(aes(xmin = pc_abun_T_lo, xmax = pc_abun_T_hi), height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "Percent change in detection (Human vs Control)",
    y = NULL,
    title = "Posterior Median and 95% CI of Percent Change (Detection)"
  ) +
  theme_minimal(base_size = 12)

####################################################################################
# Trophic group species percent-change summaries + plots
####################################################################################
library(tibble); library(dplyr); library(purrr); library(forcats); library(ggplot2); library(r2jags)

load("MFP_DetRate_Model6_TrophicGroupsP-G5.rda")

# --- helper functions ---
summ <- function(x) c(
  mean = mean(x, na.rm = TRUE),
  lo   = quantile(x, 0.025, na.rm = TRUE, names = FALSE),
  med  = median(x, na.rm = TRUE),
  hi   = quantile(x, 0.975, na.rm = TRUE, names = FALSE)
)

pct_change <- function(num, den, eps = 1e-10) {
  x <- 100 * (num - den) / (den + eps)
  x[!is.finite(x)] <- NA_real_
  x
}

as_mat <- function(x){
  if (is.null(x)) return(NULL)
  if (is.null(dim(x))) return(cbind(x))
  x
}

extract_metrics_multi <- function(mod, species_names = NULL){
  sims <- mod$BUGSoutput$sims.list
  
  A1 <- as_mat(sims$a1);  A2 <- as_mat(sims$a2)
  B1 <- as_mat(sims$b1);  B2 <- as_mat(sims$b2)
  
  H  <- as_mat(if (!is.null(sims$lam.mean.H)) sims$lam.mean.H else sims$det.H)
  C  <- as_mat(if (!is.null(sims$lam.mean.C)) sims$lam.mean.C else sims$det.C)
  
  ZH <- as_mat(if (!is.null(sims$z.sum.H)) sims$z.sum.H else sims$z.sum.t2)
  ZC <- as_mat(if (!is.null(sims$z.sum.C)) sims$z.sum.C else sims$z.sum.t1)
  
  Sp <- ncol(A1)
  if (is.null(species_names)) species_names <- paste0("sp", seq_len(Sp))
  
  rows <- vector("list", Sp)
  for(i in seq_len(Sp)){
    pc_abun_T <- pct_change(H[,i],  C[,i])
    pc_occ_T  <- pct_change(ZH[,i], ZC[,i])
    
    a1s <- summ(A1[,i]); a2s <- summ(A2[,i])
    b1s <- summ(B1[,i]); b2s <- summ(B2[,i])
    
    rows[[i]] <- tibble(
      species = species_names[i],
      
      pc_abun_T_mean = mean(pc_abun_T),
      pc_abun_T_lo   = quantile(pc_abun_T, .025),
      pc_abun_T_med  = median(pc_abun_T),
      pc_abun_T_hi   = quantile(pc_abun_T, .975),
      
      pc_occ_T_mean  = mean(pc_occ_T),
      pc_occ_T_lo    = quantile(pc_occ_T, .025),
      pc_occ_T_med   = median(pc_occ_T),
      pc_occ_T_hi    = quantile(pc_occ_T, .975),
      
      T_eff_abun_mean = a1s["mean"], T_eff_abun_lo = a1s["lo"],
      T_eff_abun_med  = a1s["med"],  T_eff_abun_hi = a1s["hi"],
      
      S_eff_abun_mean = a2s["mean"], S_eff_abun_lo = a2s["lo"],
      S_eff_abun_med  = a2s["med"],  S_eff_abun_hi = a2s["hi"],
      
      T_eff_occ_mean  = b1s["mean"], T_eff_occ_lo  = b1s["lo"],
      T_eff_occ_med   = b1s["med"],  T_eff_occ_hi  = b1s["hi"],
      
      S_eff_occ_mean  = b2s["mean"], S_eff_occ_lo  = b2s["lo"],
      S_eff_occ_med   = b2s["med"],  S_eff_occ_hi  = b2s["hi"]
    )
  }
  bind_rows(rows)
}

#y-axis labels with sample size
codes_for_labels <- c(
  "hip","ele","gir","zor","bbj","shy","war","gen",
  "hna","par","lep","gdd","wtm","gvy","smo","lio","gsq","imp","pze"
)

species_labels <- c(
  "Hippopotamus (n=54)","African Elephant (n=242)","Reticulated Giraffe (n=127)", "Zorilla (n=63)",
  "Black-backed Jackal (n=219)","Striped Hyena (n=224)","Common Warthog (n=60)", "Genets (n=311)",
  "Spotted Hyena (n=765)","Leopard (n=35)","Hares (n=2310)","Gunther's Dik-Dik (n=5547)",
  "White-tailed Mongoose (n=426)","Grevy's Zebra (n=389)","Slender Mongoose (n=150)","Lion (n=46)",
  "Ground Squirrels (n=377)", "Impala (n=2153)","Plains Zebra (n=68)"
)

label_map <- setNames(species_labels, codes_for_labels)

# trophic.group.out[[group]] has: $fit, $array, $labels
all_species_summary <- imap_dfr(
  trophic.group.out,
  ~ {
    # get the species codes actually present in this group's array
    codes_g <- dimnames(.x$array)[[3]]           # e.g., c("hip","gir","ele")
    # map to your pretty labels in the correct order for THIS group
    labels_g <- unname(label_map[codes_g])
    
    extract_metrics_multi(.x$fit, species_names = labels_g) %>%
      mutate(group = .y)
  },
  .id = NULL
)

# order species by median % change (most negative at top) and plot in a single color
df <- all_species_summary %>%
  mutate(
    species_full = species,
    species_full = fct_reorder(species_full, pc_abun_T_med),
    species_full = fct_rev(species_full)
  )

#Detection Frequency Plot
plot_data <- df %>%
  select(species_full,
         abun_T_med = pc_abun_T_med, abun_T_lo = pc_abun_T_lo, abun_T_hi = pc_abun_T_hi) %>%
  pivot_longer(
    cols = -species_full,
    names_to = c("metric", "effect", ".value"),
    names_pattern = "([a-z]+)_([TS])_(med|lo|hi)"
  )

ggplot(plot_data, aes(x = med, y = species_full, color = metric)) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 5, color = "#FABA39FF") +
  geom_point(position = position_dodge(width = 0.6), size = 24) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, position = position_dodge(width = 0.6), size = 8) +
  scale_color_manual(values = "#7A0403FF", guide = "none") +
  labs(
    x = "Percent Change (%)",
    y = NULL,
    color = NULL,
  ) +
  theme_minimal() +
  theme(text = element_text(family = "Palatino", size = 180),
        axis.text = element_text(size = 180, color = "black"),
        axis.text.y = element_text(size = 180, color = "black", margin = margin(r = 250)),
        axis.text.x = element_text(size = 180, color = "black"),
        axis.title = element_text(size = 180),
        axis.title.x = element_text(margin = margin(t = 250)),
        plot.subtitle = element_text(size = 180),
        panel.grid = element_line(color = "grey60")
  )

#Detection Frequency and Occupancy Plot
# 1) Make long df with BOTH detection (abundance) and occupancy percent changes
plot_data2 <- df %>%
  select(
    species_full,
    det_med = pc_abun_T_med, det_lo = pc_abun_T_lo, det_hi = pc_abun_T_hi,
    occ_med = pc_occ_T_med,  occ_lo = pc_occ_T_lo,  occ_hi = pc_occ_T_hi
  ) %>%
  pivot_longer(
    cols = -species_full,
    names_to = c("metric", ".value"),
    names_pattern = "(det|occ)_(med|lo|hi)"
  ) %>%
  mutate(
    metric = case_when(
      metric == "det" ~ "Detection frequency",
      metric == "occ" ~ "Occupancy"
    ),
    metric = factor(metric, levels = c("Occupancy", "Detection frequency"))
  )

# numeric y positions for dodging vertically
plot_data2 <- plot_data2 %>%
  mutate(
    y_num = as.numeric(species_full),
    y_num = ifelse(metric == "Detection frequency", y_num + 0.18, y_num - 0.18)
  )

light_red <- "#D47A7AFF"
dark_red  <- "#7A0403FF"

ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 5, color = "#FABA39FF") +

  ## Error bars (draw separately so colors match)
  geom_errorbarh(
    data = subset(plot_data2, metric == "Detection frequency"),
    aes(xmin = lo, xmax = hi, y = y_num),
    height = 0, linewidth = 8, color = dark_red
  ) +
  geom_errorbarh(
    data = subset(plot_data2, metric == "Occupancy"),
    aes(xmin = lo, xmax = hi, y = y_num),
    height = 0, linewidth = 8, color = light_red
  ) +

  ## Points (detection above; occupancy below)
  geom_point(
    data = subset(plot_data2, metric == "Detection frequency"),
    aes(x = med, y = y_num),
    shape = 16, size = 24, color = dark_red
  ) +
  geom_point(
    data = subset(plot_data2, metric == "Occupancy"),
    aes(x = med, y = y_num),
    shape = 18, size = 30, color = light_red
  ) +

  ## Put factor labels back on the y axis
  scale_y_continuous(
    breaks = seq_along(levels(df$species_full)),
    labels = levels(df$species_full),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(x = "Percent Change (%)", y = NULL) +
  theme_minimal() +
  theme(
    text = element_text(family = "Palatino", size = 180),
    axis.text = element_text(size = 180, color = "black"),
    axis.text.y = element_text(size = 180, color = "black", margin = margin(r = 250)),
    axis.text.x = element_text(size = 180, color = "black"),
    axis.title = element_text(size = 180),
    axis.title.x = element_text(margin = margin(t = 250)),
    panel.grid = element_line(color = "grey60")
  )



dark_red  <- "#7A0403FF"
light_red <- "#D47A7AFF"
yellow    <- "#FABA39FF"

# axis ranges
xr <- range(c(plot_data2$lo, plot_data2$hi, 0), na.rm = TRUE)
yr <- range(plot_data2$y_num, na.rm = TRUE)

# legend placement (outside panel, aligned with y-axis text)
ly  <- yr[2] + 1
x0  <- xr[1] - 1.1 * diff(xr)     # push into left margin to align with y labels
seg <- 0.12 * diff(xr)
pad <- 0.03 * diff(xr)

x1 <- x0
x2 <- x0 + 0.80 * diff(xr)  
x3 <- x0 + 1.70 * diff(xr)  

p <- ggplot() +
  
  ## reference line
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 5, color = yellow) +
  
  ## detection frequency (site use intensity)
  geom_segment(data = subset(plot_data2, metric == "Detection frequency"),
               aes(x = lo, xend = hi, y = y_num, yend = y_num),
               linewidth = 8, color = dark_red) +
  geom_point(data = subset(plot_data2, metric == "Detection frequency"),
             aes(x = med, y = y_num),
             shape = 16, size = 24, color = dark_red) +
  
  ## occupancy
  geom_segment(data = subset(plot_data2, metric == "Occupancy"),
               aes(x = lo, xend = hi, y = y_num, yend = y_num),
               linewidth = 8, color = light_red) +
  geom_point(data = subset(plot_data2, metric == "Occupancy"),
             aes(x = med, y = y_num),
             shape = 18, size = 34, color = light_red) +
  
  ## ---- HORIZONTAL LEGEND ABOVE PANEL ----
annotate("segment", x = x1, xend = x1 + seg, y = ly, yend = ly,
         linewidth = 8, color = dark_red) +
  annotate("point", x = x1 + 0.5*seg, y = ly,
           shape = 16, size = 24, color = dark_red) +
  annotate("text", x = x1 + seg + pad, y = ly,
           label = "Site Use Intensity", hjust = 0,
           family = "Palatino", size = 50) +
  
  annotate("segment", x = x2, xend = x2 + seg, y = ly, yend = ly,
           linewidth = 8, color = light_red) +
  annotate("point", x = x2 + 0.5*seg, y = ly,
           shape = 18, size = 34, color = light_red) +
  annotate("text", x = x2 + seg + pad, y = ly,
           label = "Occurrence Probability", hjust = 0,
           family = "Palatino", size = 50) +
  
annotate("segment", x = x3, xend = x3 + seg, y = ly, yend = ly,
         linewidth = 8, color = yellow, linetype = "dashed") +
annotate("text", x = x3 + seg + pad, y = ly,
         label = "Control", hjust = 0,
         family = "Palatino", size = 50) +
  
  ## axes
  scale_y_continuous(
    breaks = seq_along(levels(df$species_full)),
    labels = levels(df$species_full),
    expand = expansion(mult = c(0.02, 0.02))
  ) +

  coord_cartesian(xlim = xr, ylim = yr, clip = "off") +
  
  labs(x = "Percent Change (%)", y = NULL) +
  
  theme_minimal() +
  theme(
    legend.position = "none",
    text = element_text(family = "Palatino", size = 180),
    axis.text.x = element_text(size = 180, color = "black"),
    axis.text.y = element_text(size = 180, color = "black"),
    axis.title = element_text(size = 180),
    axis.title.x = element_text(margin = margin(t = 250)),
    panel.grid = element_line(color = "grey60"),
    plot.margin = margin(t = 260, r = 60, b = 60, l = 60)
  )

library(grid)
library(gtable)

g <- ggplotGrob(p)

# where is the panel located?
panel_cols <- g$layout[g$layout$name == "panel", c("l","r")]
panel_l <- panel_cols$l[1]

# add a spacer column immediately to the left of the panel
g <- gtable_add_cols(g, unit(6.0, "cm"), pos = panel_l - 1)  # <- change 2.0 cm as needed

grid.newpage()
grid.draw(g)

####################################################################################
# Multi-species percent-change summaries + plots (dependencies from trophic)
####################################################################################
# 1) labels in the exact 3rd-dimension order you printed earlier
species_labels_19 <- c(
  "Hippopotamus (n=54)",
  "African Elephant (n=242)",
  "Common Warthog (n=60)",
  "Zorilla (n=63)",
  "Reticulated Giraffe (n=127)",
  "Black-backed Jackal (n=219)",
  "Striped Hyena (n=224)",
  "Genets (n=311)",
  "Spotted Hyena (n=765)",
  "Hares (n=2310)",
  "Gunther's Dik-Dik (n=5547)",
  "Leopard (n=35)",
  "White-tailed Mongoose (n=426)",
  "Grevy's Zebra (n=389)",
  "Lion (n=46)",
  "Ground Squirrels (n=377)",
  "Slender Mongoose (n=150)",
  "Impala (n=2153)",
  "Plains Zebra (n=68)"
)

# 2) one-shot summary for the 19-species model
ms_19 <- extract_metrics_multi(det.out.MultiSpec, species_labels_19)

# 3) plot ordered by median % change (H vs C), single color
library(dplyr); library(forcats); library(ggplot2)

df19 <- ms_19 %>%
  mutate(species_full = fct_reorder(species, pc_abun_T_med) |> fct_rev())

ggplot(df19, aes(x = pc_abun_T_med, y = species_full)) +
  geom_point(color = "#7A0403", size = 2.5) +
  geom_errorbarh(aes(xmin = pc_abun_T_lo, xmax = pc_abun_T_hi),
                 height = 0, color = "#7A0403") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#FABA39FF") +
  labs(x = "Percent change in detection (Human vs Control)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())

# Bayesian p-values (global / model-level)
ppc_global <- function(fit){
  s <- fit$BUGSoutput$sims.list
  
  data.frame(
    p_FT    = mean(s$Tnew      > s$Tobs),        # Freeman–Tukey
    p_Chisq = mean(s$Chisq.new > s$Chisq.obs)    # Pearson-type chi-square
  )
}

ppc_by_group <- bind_rows(
  lapply(names(trophic.group.out), function(g){
    p <- ppc_global(trophic.group.out[[g]]$fit)
    cbind(group = g, p)
  })
)

ppc_by_group

ppc_summaries <- function(fit){
  s <- fit$BUGSoutput$sims.list
  
  data.frame(
    Tobs_mean = mean(s$Tobs),   Tnew_mean = mean(s$Tnew),
    Chisq_obs_mean = mean(s$Chisq.obs), Chisq_new_mean = mean(s$Chisq.new),
    p_FT    = mean(s$Tnew      > s$Tobs),
    p_Chisq = mean(s$Chisq.new > s$Chisq.obs)
  )
}

ppc_summary_by_group <- bind_rows(
  lapply(names(trophic.group.out), function(g){
    p <- ppc_summaries(trophic.group.out[[g]]$fit)
    cbind(group = g, p)
  })
)

ppc_summary_by_group

#########
#Elephant only outputs and plots
extract_metrics_ele <- function(mod) {
  sims <- mod$BUGSoutput$sims.list
  
  # Bayesian p-values
  p_T     <- mean(sims$Tnew      >= sims$Tobs)
  p_Chisq <- mean(sims$Chisq.new >= sims$Chisq.obs)
  
  # coefficients
  a1 <- sims$a1
  a2 <- sims$a2
  b1 <- sims$b1
  b2 <- sims$b2
  
  # percent change helper (draw-by-draw)
  pc <- function(num, den) 100 * (num - den) / den
  
  tibble(
    species = "elephant",
    
    p_T     = p_T,
    p_Chisq = p_Chisq,
    
    # Effects (posterior summaries)
    T_eff_abun_lo   = as.numeric(quantile(a1, 0.025)),
    T_eff_abun_med  = as.numeric(quantile(a1, 0.500)),
    T_eff_abun_hi   = as.numeric(quantile(a1, 0.975)),
    
    T_eff_occ_lo    = as.numeric(quantile(b1, 0.025)),
    T_eff_occ_med   = as.numeric(quantile(b1, 0.500)),
    T_eff_occ_hi    = as.numeric(quantile(b1, 0.975)),
    
    # Derived percent changes (posterior summaries)
    pc_abun_T_lo   = as.numeric(quantile(pc(sims$det.H, sims$det.C), 0.025)),
    pc_abun_T_med  = as.numeric(quantile(pc(sims$det.H, sims$det.C), 0.500)),
    pc_abun_T_hi   = as.numeric(quantile(pc(sims$det.H, sims$det.C), 0.975)),
    
    pc_occ_T_lo    = as.numeric(quantile(pc(sims$z.sum.H, sims$z.sum.C), 0.025)),
    pc_occ_T_med   = as.numeric(quantile(pc(sims$z.sum.H, sims$z.sum.C), 0.500)),
    pc_occ_T_hi    = as.numeric(quantile(pc(sims$z.sum.H, sims$z.sum.C), 0.975)),
    
  )
}

ele_summary <- extract_metrics_ele(det.out.wZT6.ele)
ele_stack <- bind_rows(
  tibble(
    metric = "Site Use Intensity",
    med = ele_summary$pc_abun_T_med,
    lo  = ele_summary$pc_abun_T_lo,
    hi  = ele_summary$pc_abun_T_hi
  ),
  tibble(
    metric = "Occurrence Probability",
    med = ele_summary$pc_occ_T_med,
    lo  = ele_summary$pc_occ_T_lo,
    hi  = ele_summary$pc_occ_T_hi
  )
) %>%
  mutate(metric = factor(metric, levels = c("Site Use Intensity","Occurrence Probability")))

ggplot(ele_stack, aes(x = med, y = metric)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 8, color = "#FABA39FF") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 5, color = "#7A0403FF") +
  geom_point(size = 16, color = "#7A0403FF") +
  coord_cartesian(xlim = c(-80, 30)) +
  scale_x_continuous(breaks = seq(-75, 25, by = 25)) +
  labs(x = "Percent Change (%)", y = NULL) +
  theme_bw(base_size = 140, base_family = "Palatino") +
  theme(
    axis.ticks.y      = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y = element_text(margin = margin(r = 50)),
    axis.text.x = element_text(margin = margin(t = 50)),
    panel.grid = element_line(color = "grey90")
  )

