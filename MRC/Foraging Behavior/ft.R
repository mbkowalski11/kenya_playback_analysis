library(tidyverse)
library(lubridate)
library(stringr)
library(purrr)
library(circular)
#library(MASS)
library(glmmTMB)
library(gridExtra)
library(ordinal)
library(broom)
library(broom.mixed)
library(lme4)
library(overlap)
library(emmeans)

# Importing + Cleaning ----------------------------------------------------
ft <- read.csv("forage_4.csv")
str(ft)

ft <- ft %>%
  mutate(
    grid_num   = as.integer(str_extract(grid, "\\d+")),       # extract grid number
    camera_num = as.integer(str_extract(camera, "(?<=T)\\d+")) # extract trap number after "T"
  ) %>%
  arrange(week, grid_num, camera_num, start) %>%
  mutate(ID = row_number()) %>%  # renumber IDs sequentially
  select(-grid_num, -camera_num) # drop helper cols if not needed

ft <- ft %>%
  mutate(
    start = ymd_hms(start),
    end   = ymd_hms(end),
    duration = as.numeric(difftime(end, start, units = "mins"))
  ) %>% relocate(duration, .after = end)

ft <- ft %>%
  mutate(
    start = ymd_hms(start),
    end   = ymd_hms(end),
    hour  = hour(start)+ minute(start)/60 + second(start)/3600,  # hour of day from `start` (0–23) in decimals
    bold  = rowSums(across(c(proximity, investigate, contact, forage)))  # sum of 4 cols
  ) %>%
  relocate(hour, .after = end) %>%      # place `hour` right after `end`
  relocate(bold, .after = forage)       # place `bold` right after `forage`

#Canonicalize species and treatment labels
ft <- ft %>%
  mutate(
    genus = genus %>%
      str_squish() %>% str_trim() %>% tolower() %>%
      str_replace_all(" ", "_") %>%
      str_replace_all("[^a-z0-9_]", ""),
    species = species %>%
      str_squish() %>% str_trim() %>% tolower() %>%
      str_replace_all(" ", "_") %>%
      str_replace_all("[^a-z0-9_]", ""),
    common = common %>%
      str_squish() %>% str_trim() %>% tolower() %>%
      str_replace_all(" ", "_") %>%
      str_replace_all("[^a-z0-9_]", ""),
    treatment = treatment %>%
      str_squish() %>% str_trim() %>% tolower(),
    treatment = case_when(
      str_detect(treatment_clean, "^human$")   ~ "human",
      str_detect(treatment_clean, "^control$") ~ "control",
      TRUE ~ NA_character_
    ),
    food_item = food_item %>%
      str_squish() %>% str_trim() %>% tolower() %>%
      str_replace_all(" ", "_") %>%
      str_replace_all("[^a-z0-9_]", "")
  ) %>%
  filter(!is.na(treatment))

#label check
#ft %>% distinct(common, common_clean)
#ft %>% distinct(treatment, treatment_clean)
#ft %>% distinct(genus, genus_clean)
#ft %>% distinct(food_item, food_item_clean)

ft <- ft %>%
  mutate(
    common = case_when(
      common %in% c("african_civet", "genet") ~ "genets",
      common %in% c("striped_ground_squirrel", "unstriped_ground_squirrel") ~ "ground_squirrels",
      TRUE ~ common
    )
  )

# 2) Summary of foraging data (for elephants)
summary_ft <- ft %>%
  filter(common == "african_elephant") %>%
  filter(target == 1) %>%
  group_by(common = common, treatment = treatment) %>%
  summarise(
    n             = n(),
    sum_forage    = sum(forage == 1, na.rm = TRUE),
    sum_light     = sum(light == 1,  na.rm = TRUE),
    avg_hour      = mean(hour,     na.rm = TRUE),
    avg_duration  = mean(duration_mins, na.rm = TRUE),
    avg_light     = mean(light,    na.rm = TRUE),
    avg_bold      = mean(bold,     na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(common, factor(treatment, levels = c("control", "human")))

write.csv(ft, "forage_2.csv", row.names = FALSE)
# Data Exploration --------------------------------------------------------
summary.ft <- ft %>%
  filter(target == 1) %>% 
  group_by(common, treatment) %>%
  summarise(
    n = n(),
    across(c(group_size, hour, duration_mins, light, bold),
           ~ mean(.x, na.rm = TRUE),
           .names = "avg_{.col}")
  ) %>%
  arrange(common, treatment)
print(summary.ft, n =100)
save(summary.ft, file = "summary_ft_nontarget.RData")
summary(ft)

# Statistics --------------------------------------------------------------
ft_target_filtered <- ft %>%
  filter(target == 1) %>%                     
  group_by(common) %>%                         
  filter(n() >= 20) %>%                        
  ungroup()

###Testing for Normality 
# Safe wrapper to handle identical or too few values
safe_shapiro <- possibly(function(x) {
  if (length(unique(x)) < 3) return(NA_real_)
  shapiro.test(x)$p.value
}, otherwise = NA_real_)
# Group size
normality_group <- ft_target_filtered %>%
  group_by(common) %>%
  summarise(
    p_value = safe_shapiro(group_size),
    .groups = "drop"
  ) %>%
  mutate(variable = "group_size")
# Duration
normality_duration <- ft_target_filtered %>%
  group_by(common) %>%
  summarise(
    p_value = safe_shapiro(duration_mins),
    .groups = "drop"
  ) %>%
  mutate(variable = "duration_mins")
# Combine + interpret
normality_tests <- bind_rows(normality_group, normality_duration) %>%
  mutate(
    p_value = round(p_value, 4),
    normal = case_when(
      is.na(p_value) ~ "Skipped (no variance / too few obs)",
      p_value > 0.05 ~ "Yes (≈ normal)",
      TRUE ~ "No (non-normal)"
    )
  )
print(normality_tests, n = 50)
#Log Tranform since Non-Normal
ft_target_filtered <- ft_target_filtered %>%
  mutate(
    log_group_size   = log1p(group_size),
    log_duration_min = log1p(duration_mins)
  )
#Log Normal Normality Test
normality_log <- ft_target_filtered %>%
  group_by(common) %>%
  summarise(
    p_group = if (length(unique(log_group_size)) > 2) shapiro.test(log_group_size)$p.value else NA_real_,
    p_dur   = if (length(unique(log_duration_min)) > 2) shapiro.test(log_duration_min)$p.value else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    normal_group = case_when(
      is.na(p_group) ~ "Skipped",
      p_group > 0.05 ~ "Yes (≈ normal)",
      TRUE ~ "No"
    ),
    normal_dur = case_when(
      is.na(p_dur) ~ "Skipped",
      p_dur > 0.05 ~ "Yes (≈ normal)",
      TRUE ~ "No"
    )
  )
print(normality_log, n = 50)

#------------------------------------------------------------
# 0) Bout Frequency  (Chi Square)
#------------------------------------------------------------
species_list <- unique(ft$common)

ft_target <- ft %>% filter(target == 1)
df <- ft_target

chisq_results <- df %>%
  filter(treatment %in% c("control", "human")) %>%
  group_by(species = common, treatment) %>%
  summarise(bouts = n(), .groups = "drop") %>%
  pivot_wider(names_from = treatment, values_from = bouts, values_fill = 0) %>%
  filter(control > 0, human > 0) %>%
  mutate(
    total_bouts = control + human,
    prop_human  = human / total_bouts,
    test        = map2(control, human, ~ chisq.test(c(.x, .y))),
    chi_squared = map_dbl(test, ~ unname(.x$statistic)),
    p_value     = map_dbl(test, ~ .x$p.value),
    sig = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ ""
    )
  ) %>%
  select(species,
         bouts_control = control,
         bouts_human   = human,
         total_bouts, prop_human, chi_squared, p_value, sig) %>%
  arrange(prop_human)

print(chisq_results, n=50)

table_grob <- gridExtra::tableGrob(chisq_results)
ggsave("chisq_results_full_table.png", table_grob, width = 10, height = 8, dpi = 500)

library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)

label_table <- c(
  african_elephant     = "African Elephant (n=20)",
  spotted_hyena        = "Spotted Hyena (n=69)",
  grevys_zebra         = "Grevy’s Zebra (n=36)",
  impala               = "Impala (n=266)",
  gunthers_dik_dik     = "Gunther's Dik-Dik (n=476)",
  hare                 = "Hares (n=101)",
  ground_squirrels     = "Ground Squirrels (n=34)",
  slender_mongoose     = "Slender Mongoose (n=21)",
  striped_hyena        = "Striped Hyena (n=31)",
  black_backed_jackal  = "Black-backed Jackal (n=26)",
  white_tailed_mongoose= "White-tailed Mongoose (n=76)",
  genets               = "Genets (n=20)"
)

# --- BH flags at the chisq table level (one row per species) ---
chisq_results2 <- chisq_results %>%
  mutate(
    p_adj  = p.adjust(p_value, method = "BH"),
    status = case_when(
      p_value < 0.05 & p_adj >= 0.05 ~ "raw_only",
      TRUE ~ "other"
    )
  )

focal_species <- names(label_table)  

status_tbl <- chisq_results %>%
  filter(species %in% focal_species) %>%
  transmute(species, p_raw = p_value, p_adj = p.adjust(p_raw, method = "BH"),
    status = case_when(
      p_raw < 0.05 & p_adj >= 0.05 ~ "raw_only",
      p_adj < 0.05                 ~ "bh_sig",
      TRUE                         ~ "ns"))

# --- build plotting df (keep status, and keep wide counts for ordering) ---
dfp <- chisq_results2 %>%
  filter(total_bouts >= 20) %>%
  mutate(
    species_lab = label_table[species],
    species_lab = forcats::fct_reorder(species_lab, bouts_human, .desc = TRUE)
  ) %>%
  tidyr::pivot_longer(
    cols = c(bouts_control, bouts_human),
    names_to = "treatment",
    values_to = "bouts"
  ) %>%
  mutate(
    treatment = if_else(treatment == "bouts_control", "Control", "Human"),
    treatment = factor(treatment, levels = c("Human","Control"))
  )

master_y <- master_y <- levels(dfp$species_lab)

dark_red <- "#7A0403FF"; yellow <- "#FABA39FF"

p <- ggplot(dfp, aes(bouts, species_lab, group = species_lab)) +
  geom_line(color = "black", linewidth = 8) +
  geom_point(aes(color = treatment, shape = treatment), size = 30) +
  geom_point(
    data = subset(dfp, status == "raw_only" & treatment == "Human"), aes(bouts, species_lab),
    inherit.aes = FALSE, shape = 21, fill = "white", color = "white", size = 12) +
  geom_point(
    data = subset(dfp, status == "raw_only" & treatment == "Control"),
    aes(bouts, species_lab),inherit.aes = FALSE, shape = 24, fill = "white", color = "white", size = 14) +
  scale_color_manual(breaks = c("Human","Control"), values = c(Human = dark_red, Control = yellow)) +
  scale_shape_manual(breaks = c("Human","Control"), values = c(Human = 16, Control = 17)) +
  guides(
    shape = "none",
    color = guide_legend(override.aes = list(shape = c(16,17), size = 30),
      nrow = 1, byrow = TRUE)) +
  labs(x = "Foraging Bouts", y = NULL) +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text  = element_text(margin = margin(r = 120)),
    text = element_text(family = "Palatino", size = 180),
    axis.text.x = element_text(size = 180, color = "black"),
    axis.text.y = element_text(size = 180, color = "black"),
    axis.title  = element_text(size = 180),
    axis.title.x = element_text(margin = margin(t = 50)),
    panel.grid = element_line(color = "grey70")
  )

# --- optional: add space between y labels and panel (your gtable trick) ---
library(grid); library(gtable)
g <- ggplotGrob(p)
panel_l <- g$layout[g$layout$name == "panel", "l"][1]
g <- gtable_add_cols(g, unit(6.0, "cm"), pos = panel_l - 1)  # adjust gap
grid.newpage(); grid.draw(g)

###for triple plot with maaster y axis
# 1) add BH + status (raw sig but fails BH => hollow)
dfp2 <- dfp %>%
  left_join(status_tbl, by = "species") %>%
  mutate(status_bh12 = status.y) %>%                    
  tidyr::complete(species_lab = master_y, treatment   = levels(dfp$treatment)) %>%
  mutate(species_lab  = factor(species_lab, levels = master_y),
    treatment    = factor(treatment, levels = c("Human","Control")),
    status_bh12  = tidyr::replace_na(status_bh12, "ns"))

# 2) plot: filled points for all, then hollow overlay ONLY for raw_only
p_chi <- ggplot(dfp2, aes(x = bouts, y = species_lab, group = species_lab)) +
  geom_line(color = "black", linewidth = 8, na.rm = TRUE) +
  geom_point(aes(shape = treatment, fill = treatment, color = treatment),
    size = 30, stroke = 4, na.rm = TRUE) +
  geom_point(data = subset(dfp2, status_bh12 == "raw_only"),
    aes(shape = treatment, color = treatment),
    fill = "white", size = 20,stroke = 9,na.rm = TRUE) +
  scale_shape_manual(values = c(Human = 21, Control = 24), breaks = c("Human","Control")) +
  scale_fill_manual(values  = c(Human = dark_red, Control = yellow), breaks = c("Human","Control")) +
  scale_color_manual(values = c(Human = dark_red, Control = yellow), breaks = c("Human","Control")) +
  scale_y_discrete(drop = FALSE, limits = master_y) +
  guides(shape = guide_legend(override.aes = list(size = 30, stroke = 4), 
                              nrow = 1, byrow = TRUE), fill = "none", color = "none") +
  labs(x = "Foraging Bouts", y = NULL) +
  theme_minimal() +
  theme(
    legend.position="none",
    legend.title=element_blank(),
    legend.text=element_text(margin=margin(r=120)),
    text=element_text(family="Palatino", size=180),
    axis.text.x=element_text(size=180, color="black"),
    axis.text.y=element_text(size=180, color="black"),
    axis.title=element_text(size=180),
    axis.title.x=element_text(margin=margin(t=50)),
    panel.grid=element_line(color="grey70"))

library(grid); library(gtable)
g_chi <- ggplotGrob(p_chi)
panel_l <- g_chi$layout[g_chi$layout$name == "panel", "l"][1]
g_chi <- gtable_add_cols(g_chi, unit(6.0, "cm"), pos = panel_l - 1)  # adjust gap
grid.newpage(); grid.draw(g_chi)

#------------------------------------------------------------
# 1) Group size  (Wilcoxon)
#------------------------------------------------------------
hist(ft$group, breaks = 30)

Wilcoxon_group <- ft_target_filtered %>%
  group_by(common) %>%
  summarise({
    w <- wilcox.test(group_size ~ treatment, exact = FALSE)
    tibble(
      statistic = unname(w$statistic),
      parameter = NA_real_,
      p_value   = w$p.value
    )
  }, .groups = "drop")

#------------------------------------------------------------
# 2) Hour  (Watson–Wheeler circular)
#------------------------------------------------------------
WW_hour <- ft %>%
  group_by(common) %>%
  group_modify(~{
    df <- .x %>%
      filter(treatment %in% c("control","human"), !is.na(hour))
    
    # need at least 2 obs per treatment
    if (min(table(df$treatment)) < 2) {
      return(tibble(statistic = NA_real_, parameter = NA_real_, p_value = NA_real_))
    }
    
    # one circular vector + one grouping factor
    theta <- circular((df$hour %% 24) / 24 * 2*pi, units = "radians", modulo = "2pi")
    grp   <- factor(df$treatment, levels = c("control","human"))
    
    ww <- circular::watson.wheeler.test(theta, grp)
    
    tibble(
      statistic = unname(ww$statistic),
      parameter = unname(ww$parameter),  # should be 2 for two groups
      p_value   = ww$p.value
    )
  }) %>%
  ungroup()
print(WW_hour, n=25)

overlap_table <- ft %>%
  filter(treatment %in% c("control","human")) %>%
  group_by(common) %>%
  summarise(
    n_control = sum(treatment=="control"),
    n_human   = sum(treatment=="human"),
    overlap   = tryCatch({
      overlapEst((hour[treatment=="control"]/24)*2*pi,
                 (hour[treatment=="human"]/24)*2*pi,
                 type = ifelse(min(n_control, n_human) > 50, "Dhat4", "Dhat1"))
    }, error = function(e) NA_real_)
  ) %>%
  arrange(overlap)
print(overlap_table, n=50)

species_keep <- c(
  "striped_hyena",
  "grevys_zebra",
  "african_elephant",
  "black_backed_jackal",
  "genets",
  "slender_mongoose",
  "impala",
  "ground_squirrels",
  "spotted_hyena",
  "white_tailed_mongoose",
  "hare",
  "gunthers_dik_dik"
)

# Filter & prep
df12 <- ft_target %>%
  filter(common %in% species_keep,
         treatment %in% c("control","human"),
         !is.na(hour)) %>%
  mutate(
    treatment = factor(treatment, levels = c("control","human")),
    common    = factor(common, levels = species_keep),    # facet order
    hour      = (hour %% 24)                              # keep in [0,24)
  )

# Faceted density + rugs
p <- ggplot(df12, aes(x = hour, color = treatment, fill = treatment)) +
  geom_density(alpha = 0.30, adjust = 1) +
  geom_rug(aes(color = treatment), sides = "b", alpha = 0.35) +
  scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, 6)) +
  labs(
    x = "Hour of day",
    y = "Density",
    title = "Activity hour density by treatment (control vs human)"
  ) +
  facet_wrap(~ common, ncol = 3, scales = "fixed") +
  theme_bw(base_size = 13) +
  ylim(c(0,.3)) +
  theme(legend.position = "top")

print(p)
#------------------------------------------------------------
# 3) Duration  (Wilcoxon)
#------------------------------------------------------------
ft_filtered_bold <- ft %>%
  filter(bold > 0)

summary(ft_filtered_bold$duration_mins)
hist(ft_filtered_bold$duration_mins[ft_filtered_bold$duration_mins < 60], breaks = 30)
str(ft$duration_mins)

Wilcoxon_duration <- ft_target %>%
  filter(common == "african_elephant") %>%
  summarise({
    w <- wilcox.test(duration_mins ~ treatment, exact = FALSE)
    tibble(
      statistic = unname(w$statistic),
      parameter = NA_real_,
      p_value   = w$p.value
    )
  }, .groups = "drop")

library(broom)

###log-linear model for duration min
fit_lognorm_one <- function(sp) {
  
  dat <- ft_target %>%
    filter(common == sp) %>%
    filter(is.finite(duration_mins), duration_mins > 0) %>%
    mutate(treatment = factor(treatment, levels = c("control","human")))
  
  if (nrow(dat) < 5 || n_distinct(dat$treatment) < 2) return(NULL)
  
  mod <- lm(log(duration_mins) ~ treatment, data = dat)
  
  # treatment effect on log scale (human vs control)
  td <- broom::tidy(mod, conf.int = TRUE) %>%
    filter(term == "treatmenthuman")
  if (nrow(td) == 0) return(NULL)
  
  # mu = E[log(Y)] by treatment
  emm_df <- as.data.frame(emmeans(mod, ~ treatment))  # emmean = mu_hat (log scale)
  
  mu_control <- emm_df$emmean[emm_df$treatment == "control"]
  mu_human   <- emm_df$emmean[emm_df$treatment == "human"]
  
  # geometric means (aka model-based medians) on original scale
  geo_control <- exp(mu_control)
  geo_human   <- exp(mu_human)
  
  tibble(
    species        = sp,
    n              = nrow(dat),
    
    geo_control    = geo_control,   # <-- geometric mean / median minutes (model-based)
    geo_human      = geo_human,     # <-- geometric mean / median minutes (model-based)
    
    estimate_log   = td$estimate[1],   # log ratio
    se_log         = td$std.error[1],
    t              = td$statistic[1],
    p_value        = td$p.value[1],
    
    ratio          = geo_human / geo_control,      # = exp(estimate_log)
    ratio_l95      = exp(td$conf.low[1]),
    ratio_u95      = exp(td$conf.high[1]),
    percent_change = (geo_human / geo_control - 1) * 100
  )
}

species_list <- unique(ft_target$common)
view(lognorm_results)
lognorm_results <- species_list %>%
  map(fit_lognorm_one) %>%
  compact() %>%
  bind_rows() %>%
  mutate(sig = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE ~ ""
  )) %>%
  arrange(percent_change) %>%
  filter(n > 19)

lognorm_results
table_grob <- gridExtra::tableGrob(lognorm_results)
ggsave("lognormal_results_target_table.png", table_grob, width = 14, height = 6, dpi = 500)

dur_df <- lognorm_results %>% 
  mutate(
    species_lab = factor(label_table[species], levels = label_table),
    species_lab = forcats::fct_reorder(species_lab, ratio))

dur_df <- dur_df %>%
  mutate(
    p_adj  = p.adjust(p_value, method = "BH"),
    status = dplyr::case_when(
      p_value < 0.05 & p_adj >= 0.05 ~ "raw_only",
      TRUE ~ "other"))

dur_df2 <- dur_df %>%
  mutate(species_lab = factor(unname(label_table[species]), levels = master_y)) %>%
  tidyr::complete(species_lab = master_y)

p_duration <- p_duration <- ggplot(dur_df2, aes(x = ratio, y = species_lab)) +
  geom_vline(xintercept = 1, linetype="dashed", linewidth=6, color=yellow) +
  geom_errorbarh(aes(xmin = ratio_l95, xmax = ratio_u95), height=0, linewidth=8, color=dark_red, na.rm=TRUE) +
  geom_point(size=30, stroke=4, color=dark_red, na.rm=TRUE) +
  scale_y_discrete(drop = FALSE, limits = master_y) +
  scale_x_log10(breaks=c(0.1,1,10), labels=c("0.1","1.0","10.0"))+
  geom_point(
    data = subset(dur_df, status == "raw_only"),
    aes(ratio, species_lab),
    inherit.aes = FALSE,
    shape = 21, fill = "white", color = dark_red,
    size =20, stroke = 9) +
  labs(x = "Human / Control Bout Duration", y = NULL) +
  theme_minimal() +
  theme(
    legend.position = "none",
    text = element_text(family = "Palatino", size = 180),
    axis.text.x = element_text(size = 180, color = "black"),
    axis.text.y = element_text(size = 180, color = "black"),
    axis.title.x = element_text(size = 180, margin = margin(t = 50)),
    panel.grid = element_line(color = "grey60")
  )
p_duration
# > # Histogram and QQ plot
#   > hist(residuals(lm_log_dur), main = "Residuals of log-linear model", xlab = "Residuals")
# > qqnorm(residuals(lm_log_dur))
# > qqline(residuals(lm_log_dur))
# > # Shapiro–Wilk test (for small n)
#   > shapiro.test(residuals(lm_log_dur))
# 
# Shapiro-Wilk normality test
# 
# data:  residuals(lm_log_dur)
# W = 0.98853, p-value = 0.9957


fit_gamma_one <- function(sp) {
  dat <- ft_filtered_bold %>% filter(common == sp)
  if (nrow(dat) < 5 || dplyr::n_distinct(dat$treatment) < 2) return(NULL)
  
  mod <- try(
    glm(duration_mins ~ treatment,
        data = dat,
        family = Gamma(link = "log")),
    silent = TRUE
  )
  if (inherits(mod, "try-error")) return(NULL)
  
  # Treatment effect (log scale)
  td <- broom::tidy(mod) %>% filter(term == "treatmenthuman")
  if (nrow(td) == 0) return(NULL)
  
  est <- td$estimate[1]
  se  <- td$std.error[1]
  ci_log <- est + c(-1.96, 1.96) * se
  
  # --- Mean duration per treatment (response scale) ---
  emm <- emmeans::emmeans(mod, ~ treatment, type = "response")
  emm_df <- as.data.frame(emm)
  
  mean_control <- emm_df$response[emm_df$treatment == "control"]
  mean_human   <- emm_df$response[emm_df$treatment == "human"]
  
  tibble(
    species         = sp,
    n               = nrow(dat),
    
    mean_control    = mean_control,
    mean_human      = mean_human,
    
    estimate        = est,               # log ratio
    se              = se,
    z               = td$statistic[1],
    p_value         = td$p.value[1],
    
    ratio           = exp(est),           # human / control
    ratio_l95       = exp(ci_log[1]),
    ratio_u95       = exp(ci_log[2]),
    percent_change  = (exp(est) - 1) * 100
  )
}

gamma_results <- species_list %>%
  map(fit_gamma_one) %>%
  purrr::compact() %>%
  bind_rows() %>%
  mutate(sig = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE ~ ""
  )) %>%
  arrange(percent_change)

gamma_results <- gamma_results %>% filter(n > 20)
gamma_results

table_grob <- gridExtra::tableGrob(gamma_results)
ggsave("gamma_results_full_table.png", table_grob, width = 14, height = 6, dpi = 500)

#------------------------------------------------------------
# 4) Light  (Logistic regression)
#------------------------------------------------------------
m_light_comm <- glmer(
  light ~ treatment + (1 | common),
  data = ft,
  family = binomial(link = "logit")
)
summary(m_light_comm)

m_light_comm_rs <- glmer(
  light ~ treatment + (treatment | common),
  data = ft,
  family = binomial("logit")
)
summary(m_light_comm_rs)

anova(m_light_comm, m_light_comm_rs) 

ranef_df <- coef(m_light_comm_rs)$common %>%
  tibble::rownames_to_column("species") %>%
  select(species, `(Intercept)`, treatmenthuman) %>%
  mutate(OR_species = exp(treatmenthuman)) %>%
  arrange(desc(OR_species))
print(ranef_df)

# Species-specific intercepts and slopes
coefs <- coef(m_light_comm_rs)$common %>%
  tibble::rownames_to_column("species") %>%
  select(species, `(Intercept)`, treatmenthuman)

# Calculate probabilities
pred_sp_trt <- coefs %>%
  mutate(
    p_control = plogis(`(Intercept)`),
    p_human   = plogis(`(Intercept)` + treatmenthuman),
    diff      = p_human - p_control
  ) %>%
  arrange(desc(diff))

# Now your plotting code will work:
p1 <- ggplot(pred_sp_trt, aes(y = reorder(species, diff))) +
  geom_point(aes(x = p_control), color = "blue", size = 2) +
  geom_point(aes(x = p_human), color = "red", size = 2) +
  geom_segment(aes(x = p_control, xend = p_human, y = species, yend = species),
               arrow = arrow(length = unit(0.1, "cm")), color = "grey40") +
  labs(
    x = "Predicted P(light = 1)",
    y = "Species",
    title = "Species-specific light probability (control vs human)"
  ) +
  scale_x_continuous(limits = c(0, 1)) +
  theme_bw(base_size = 13)

ggsave("LRMElight_full.png", p1, width = 10, height = 6, dpi = 500)

# Per-species profile of what’s missing
profile <- ft_target %>%
  dplyr::group_by(common) %>%
  dplyr::summarise(
    n          = dplyr::n(),
    n_treat    = dplyr::n_distinct(treatment),
    n_light    = dplyr::n_distinct(light),
    has_both_treat  = n_treat >= 2,
    has_both_light  = n_light >= 2,
    .groups = "drop"
  ) %>%
  dplyr::arrange(!has_both_treat, !has_both_light, desc(n))
profile

## 1) Use your profile to define valid species
valid_species <- profile %>%
  filter(has_both_treat, has_both_light) %>%
  pull(common)

## 2) Per-species logistic GLM (human vs control)
fit_logit_one <- function(sp) {
  d <- ft_target %>%
    filter(common == sp) %>%
    droplevels() %>%
    mutate(treatment = factor(treatment, levels = c("control","human")))
  
  m  <- glm(light ~ treatment, data = d, family = binomial("logit"))
  td <- tidy(m) %>% filter(term == "treatmenthuman")
  est <- td$estimate[1]; se <- td$std.error[1]; ci <- est + c(-1.96, 1.96)*se
  
  tibble(
    species = sp,
    n       = nrow(d),
    estimate= est, se = se, z = td$statistic[1], p_value = td$p.value[1],
    OR      = exp(est),
    sig     = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ ""
    )
  )
}

## 3) Fit only valid species and sort
logit_light <- valid_species %>%
  map(fit_logit_one) %>%
  bind_rows() %>%
  arrange(p_value)

## 4) (Optional) add per-treatment prevalence for context
prevs <- ft_target %>%
  filter(common %in% valid_species) %>%
  group_by(common, treatment) %>%
  summarise(n = n(), prop_light = mean(light), .groups = "drop") %>%
  pivot_wider(names_from = treatment, values_from = c(n, prop_light),
              names_sep = "_")

logit_light <- logit_light %>%
  left_join(prevs, by = c("species" = "common"))

print(logit_light)

table_grob <- gridExtra::tableGrob(logit_light)
ggsave("LRlight_results_target_table.png", table_grob, width = 14, height = 6, dpi = 500)

library(dplyr)
library(emmeans)
library(tibble)

fit_diurnal_one <- function(sp) {
  
  d <- ft_target %>%
    filter(common == sp) %>%
    mutate(treatment = factor(treatment, levels = c("control","human")))
  
  if (nrow(d) < 5 || dplyr::n_distinct(d$treatment) < 2) return(NULL)
  
  m <- glm(light ~ treatment, data = d, family = binomial)
  
  em_prob <- regrid(emmeans(m, ~ treatment), transform = "response")
  emdf <- as.data.frame(em_prob)
  
  pC <- emdf$prob[emdf$treatment == "control"]
  pH <- emdf$prob[emdf$treatment == "human"]
  
  ct <- as.data.frame(summary(
    contrast(em_prob, method = list("Human - Control" = c(-1, 1))),
    infer = c(TRUE, TRUE)
  ))
  
  tibble(
    species   = sp,
    n         = nrow(d),
    p_control = pC,
    p_human   = pH,
    delta_p   = ct$estimate[1],
    lo        = ct$asymp.LCL[1],
    hi        = ct$asymp.UCL[1],
    p_value   = ct$p.value[1]
  )
}
fit_diurnal_one("impala")
# --- build diurnal_df (your code) ---
diurnal_df <- purrr::map_dfr(valid_species, fit_diurnal_one) %>%
  mutate(p_adj  = p.adjust(p_value, method = "BH"),
    status = case_when(p_value < 0.05 & p_adj >= 0.05 ~ "raw_only", TRUE ~ "other"),
    species_lab = factor(label_table[species], levels = label_table),
    species_lab = forcats::fct_reorder(species_lab, delta_p))

table(diurnal_df$status)  # should show how many raw_only
table_grob <- gridExtra::tableGrob(diurnal_df)
ggsave("LRlight_results_response_target_table.png", table_grob, width = 14, height = 6, dpi = 500)

diurnal_df2 <- diurnal_df %>%
  mutate(species_lab = factor(unname(label_table[species]), levels = master_y)) %>%
  tidyr::complete(species_lab = master_y)

# --- plot with overlay ---
p_diurnal <- p_diurnal <- ggplot(diurnal_df2, aes(x = delta_p, y = species_lab)) +
  geom_vline(xintercept = 0, linetype="dashed", linewidth=6, color=yellow) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height=0, linewidth=8, color=dark_red, na.rm=TRUE) +
  geom_point(size=30, color=dark_red, na.rm=TRUE) +
  scale_y_discrete(drop = FALSE, limits = master_y) +
  geom_point(
    data = subset(diurnal_df, status == "raw_only"), aes(x = delta_p, y = species_lab),
    inherit.aes = FALSE, shape = 21, fill = "white", color = "white", size = 14) +
  scale_x_continuous(
    limits = range(c(diurnal_df$lo, diurnal_df$hi), na.rm = TRUE),
    expand = expansion(mult = 0.05)) +
  labs(x = "Δ Diurnal Proportion (Human − Control)", y = NULL) +
  theme_minimal() +
  theme(
    legend.position = "none",
    text = element_text(family = "Palatino", size = 180),
    axis.text.x = element_text(size = 180, color = "black"),
    axis.text.y = element_text(size = 180, color = "black"),
    axis.title.x = element_text(size = 180, margin = margin(t = 50)),
    panel.grid = element_line(color = "grey60")
  )
p_diurnal

####Binomial Tests and Fisher's Exact with Paired Plot
binom_summary_one <- function(sp) {
  d <- ft_target %>% filter(common == sp) %>%
    mutate(treatment = factor(treatment, levels = c("control", "human")))
  if (nrow(d) < 5 || n_distinct(d$treatment) < 2) return(NULL)
  # summarize successes / trials
  sumdat <- d %>% group_by(treatment) %>% 
    summarise(successes = sum(light), trials = n(), .groups = "drop")
  # exact binomial CIs
  ci_df <- sumdat %>% rowwise() %>%
    mutate(bt = list(binom.test(successes, trials)),
      prop = successes / trials, lo = bt$conf.int[1], hi = bt$conf.int[2]) %>%
    ungroup() %>% select(-bt)
  # Fisher's exact test
  mat <- matrix(c(sumdat$successes[sumdat$treatment == "control"],
      sumdat$trials[sumdat$treatment == "control"] -
        sumdat$successes[sumdat$treatment == "control"],
      sumdat$successes[sumdat$treatment == "human"],
      sumdat$trials[sumdat$treatment == "human"] -
        sumdat$successes[sumdat$treatment == "human"]),
    nrow = 2, byrow = TRUE)
  
  ft <- fisher.test(mat)
  
  ci_df %>% mutate(species = sp, p_value = ft$p.value, odds_ratio = ft$estimate)}

binom_df <- map_dfr(valid_species, binom_summary_one)
binom_df

species_stats <- binom_df %>%
  distinct(species, p_value, odds_ratio) %>%
  mutate(sig = case_when(p_value < 0.001 ~ "***", p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*", TRUE ~ ""))
species_stats

yellow   <- "#FABA39FF"
dark_red <- "#7A0403FF"

plot_df <- binom_df %>%
  mutate(p_adj = p.adjust(p_value, "BH"),
         status = if_else(p_value < .05 & p_adj >= .05, "raw_only", "other"),
         trt = factor(treatment, c("control","human")),
         species_lab = factor(species_lab, levels = master_y)) %>%
  tidyr::complete(species_lab = factor(master_y, levels = master_y),
                  trt = factor(c("control","human"), levels=c("control","human"))) %>%
  mutate(y0 = as.numeric(species_lab),
         y  = y0 + if_else(trt=="control", +0.22, -0.22))

p_binom <- ggplot(plot_df, aes(y = y)) +
  geom_errorbarh(aes(x = prop, xmin = lo, xmax = hi, color = trt),
                 height=0, linewidth=8, na.rm=TRUE) +
  geom_point(aes(x = prop, shape = trt, fill = trt, color = trt),
             size=30, stroke=4, na.rm=TRUE) +
  geom_point(data = subset(plot_df, status=="raw_only"),
             aes(x = prop, y = y, shape = trt, color = trt),
             inherit.aes = FALSE, size=20, fill="white", stroke=9, na.rm=TRUE) +
  scale_shape_manual(values=c(control=24, human=21)) +
  scale_fill_manual(values=c(control=yellow, human=dark_red)) +
  scale_color_manual(values=c(control=yellow, human=dark_red)) +
  scale_y_continuous(breaks = seq_along(master_y), labels = master_y) +
  scale_x_continuous(limits = c(-.05, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("0", "0.25", "0.5", "0.75", "1")) +
  labs(x="Proportion of Diurnal Bouts", y=NULL) +
  theme_minimal() +
  theme(legend.position="none",
        text=element_text(family="Palatino", size=180),
        axis.text=element_text(size=180, color="black"),
        axis.title.x=element_text(size=180, margin=margin(t=50)),
        panel.grid=element_line(color="grey60"),
        panel.grid.minor = element_blank())

p_binom
####Triple plot (bout frequency, diurnal, and duration)
library(patchwork)
std_theme <- theme(
  plot.margin  = margin(20, 20, 20, 20),
  axis.title.x = element_text(margin = margin(t = 50)),
  axis.text.x  = element_text(margin = margin(t = 10))
)

p_left  <- p_chi 
p_mid   <- p_binom + theme(axis.text.y = element_blank(), axis.title.y = element_blank())
p_right <- p_duration + theme(axis.text.y = element_blank(), axis.title.y = element_blank())

(p_left | p_mid | p_right) +
  patchwork::plot_layout(widths = c(1.1, 1, 1))
#------------------------------------------------------------
# 5) Boldness  (Ordinal logistic regression)
#------------------------------------------------------------
library(MASS)

species_keep <- ft_filtered_bold %>%
  group_by(common) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 19) %>%
  pull(common)

OLR_bold <- ft %>%
  filter(common %in% species_keep) %>%
  mutate(
    treatment = factor(treatment, levels = c("control", "human")),
    bold      = factor(bold, levels = 0:4, ordered = TRUE)
  ) %>%
  group_by(common) %>%
  group_modify(~{
    dat <- .x %>%
      filter(!is.na(bold), !is.na(treatment)) %>%
      droplevels()
    
    if (n_distinct(dat$treatment) < 2 || n_distinct(dat$bold) < 2) return(tibble())
    
    # counts per treatment (still useful)
    n_w <- dat %>%
      count(treatment, name = "n") %>%
      pivot_wider(names_from = treatment, values_from = n, names_prefix = "n_")
    
    fit <- try(polr(bold ~ treatment, data = dat, Hess = TRUE, method = "logistic"), silent = TRUE)
    if (inherits(fit, "try-error")) return(tibble())
    
    # --- model-based expected ordinal score by treatment ---
    # predicted category probabilities for each treatment
    pr_control <- predict(fit, newdata = data.frame(treatment = factor("control", levels = c("control","human"))),
                          type = "probs")
    pr_human   <- predict(fit, newdata = data.frame(treatment = factor("human", levels = c("control","human"))),
                          type = "probs")
    
    # numeric scores corresponding to ordered factor levels
    k <- as.numeric(levels(dat$bold))
    
    expected_control <- sum(k * as.numeric(pr_control))
    expected_human   <- sum(k * as.numeric(pr_human))
    
    # --- inference on treatment effect (same as you did) ---
    ct <- coef(summary(fit))
    ridx <- grep("^treatment", rownames(ct))
    if (length(ridx) == 0) return(tibble())
    
    est <- unname(ct[ridx[1], "Value"])
    se  <- unname(ct[ridx[1], "Std. Error"])
    z   <- est / se
    p   <- 2 * pnorm(abs(z), lower.tail = FALSE)
    
    tibble(
      species          = unique(dat$common),
      n_control        = n_w$n_control,
      n_human          = n_w$n_human,
      
      expected_control = expected_control,
      expected_human   = expected_human,
      
      estimate         = est,
      se               = se,
      z                = z,
      p_value          = p,
      OR               = exp(est)
    )
  }) %>%
  ungroup() %>%
  mutate(sig = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE ~ ""
  )) %>%
  arrange(OR)

OLR_bold <- OLR_bold %>%
  mutate(
    p_adj  = p.adjust(p_value, method = "BH"),
    status = dplyr::case_when(
      p_value < 0.05 & p_adj >= 0.05 ~ "raw_only",
      TRUE ~ "other"
    )
  )

OLR_bold

table_grob <- gridExtra::tableGrob(OLR_bold)
ggsave(filename = "OLR_bold_results_full_table.png", plot = table_grob, width = 16, height = 6, dpi = 500)

library(dplyr); library(tidyr); library(ggplot2); library(forcats)

dark_red <- "#7A0403FF"
yellow   <- "#FABA39FF"

# label_table should map common-name keys -> "Nice Name (n=...)"
label_table <- c(
  african_elephant     = "African Elephant (n=54)",
  spotted_hyena        = "Spotted Hyena (n=146)",
  grevys_zebra         = "Grevy’s Zebra (n=102)",
  impala               = "Impala (n=388)",
  gunthers_dik_dik     = "Gunther's Dik-Dik (n=1208)",
  hare                 = "Hares (n=272)",
  ground_squirrels     = "Ground Squirrels (n=70)",
  slender_mongoose     = "Slender Mongoose (n=35)",
  striped_hyena        = "Striped Hyena (n=60)",
  black_backed_jackal  = "Black-backed Jackal (n=43)",
  white_tailed_mongoose= "White-tailed Mongoose (n=143)",
  genets               = "Genets (n=36)"
)

bold_plot_df <- OLR_bold %>%
  mutate(
    species_lab = unname(label_table[common]),
    species_lab = forcats::fct_reorder(species_lab, expected_human, .desc = FALSE)
  ) %>%
  tidyr::pivot_longer(
    cols = c(expected_control, expected_human),
    names_to = "treatment",
    values_to = "expected"
  ) %>%
  mutate(
    treatment = ifelse(treatment == "expected_human", "Human", "Control"),
    treatment = factor(treatment, levels = c("Human", "Control"))
  )

p_bold <- ggplot(bold_plot_df, aes(x = expected, y = species_lab, group = species_lab)) +
  geom_line(color = "black", linewidth = 8) +
  geom_point(aes(color = treatment, shape = treatment), size = 30) +
  geom_point(
    data = subset(bold_plot_df, status == "raw_only" & treatment == "Human"),
    aes(expected, species_lab),
    inherit.aes = FALSE,
    shape = 21, fill = "white", color = "white",
    size = 20
  ) +
  geom_point(
    data = subset(bold_plot_df, status == "raw_only" & treatment == "Control"),
    aes(expected, species_lab),
    inherit.aes = FALSE,
    shape = 24, fill = "white", color = "white",
    size = 20
  ) +
  scale_color_manual(values = c(Human = dark_red, Control = yellow),
                     breaks = c("Human","Control")) +
  scale_shape_manual(values = c(Human = 16, Control = 17),
                     breaks = c("Human","Control")) +
  guides(shape = "none",
         color = guide_legend(override.aes = list(shape = c(16,17), size = 30),
                              nrow = 1, byrow = TRUE)) +
  scale_x_continuous(limits = c(0, 4.15), breaks = seq(0, 4, by = 1),
                     expand = expansion(mult = 0)) +
  labs(x = "Expected Bout Boldness (0–4)", y = NULL) +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(margin = margin(r = 120)),
    text = element_text(family = "Palatino", size = 180),
    axis.text.x = element_text(size = 180, color = "black"),
    axis.text.y = element_text(size = 180, color = "black"),
    axis.title  = element_text(size = 180),
    axis.title.x = element_text(margin = margin(t = 50)),
    panel.grid = element_line(color = "grey70")
  )
p_bold

#------------------------------------------------------------
# Optional: round and preview each table
#------------------------------------------------------------
res_group_size <- res_group_size %>% mutate(across(c(statistic, p_value), round, 3))
res_hour        <- res_hour        %>% mutate(across(c(statistic, p_value), round, 3))
res_duration    <- res_duration    %>% mutate(across(c(statistic, p_value), round, 3))
res_light       <- res_light       %>% mutate(across(c(statistic, p_value), round, 3))
res_bold        <- res_bold        %>% mutate(across(c(statistic, p_value), round, 3))

# Inspect
print(Wilcoxon_group)
print(WW_hour)
print(Wilcoxon_duration)
print(LR_light)
print(OLR_bold)

# Save individually if you like
# write.csv(res_group_size, "res_group_size.csv", row.names = FALSE)
# write.csv(res_hour,        "res_hour.csv",        row.names = FALSE)
# write.csv(res_duration,    "res_duration.csv",    row.names = FALSE)
# write.csv(res_light,       "res_light.csv",       row.names = FALSE)
# write.csv(res_bold,        "res_bold.csv",        row.names = FALSE)

###########
#elephant plots - duration
ele_log <- lognorm_results %>%
  filter(species == "african_elephant")

ggplot(ele_log, aes(y = "Bout Duration")) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 8, colour = "#FABA39FF") +
  geom_errorbarh(
    aes(xmin = ratio_l95, xmax = ratio_u95),
    height = 0,
    linewidth = 5,
    colour = "#7A0403FF"
  ) +
  geom_point(aes(x = ratio), size = 16, colour = "#7A0403FF") +
  scale_x_log10(
    breaks = c(0.01, 0.1, 1),
    labels = c("0.01", "0.1", "1"),
    limits = c(0.01, 2)
  ) +
  labs(
    x = "Human Treatment / Control Treatment",
    y = NULL
  ) +
  theme_bw(base_size = 140, base_family = "Palatino") +
  theme(
    axis.ticks.y       = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )

# 2) fit the SAME model
ele_mod <- lm(log(duration_mins) ~ trt, data = ele_dat)

# 3) model-based geometric means (medians) + CI on original scale
# emmean is on log scale; type="response" back-transforms via exp()
emm_raw <- as.data.frame(emmeans(ele_mod, ~ trt, type = "response"))
ele_emm <- emm_raw %>%
  mutate(
    trt_lab = case_when(
      trt == "control" ~ "C",
      trt == "human"   ~ "H",
      TRUE ~ NA_character_
    ),
    trt_lab = factor(trt_lab, levels = c("C","H"))
  ) %>%
  transmute(
    trt_lab,
    geo = response,
    l95 = lower.CL,
    u95 = upper.CL
  )

# 4) plot
ggplot(ele_dat, aes(x = trt_lab, y = duration_mins)) +
  geom_point(position = position_jitter(width = 0.08), size = 5, alpha = 0.8) +
  geom_errorbar(
    data = ele_emm,
    aes(x = trt_lab, ymin = l95, ymax = u95),
    width = 0.12, linewidth = 2,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = ele_emm,
    aes(x = trt_lab, y = geo),
    size = 10,
    inherit.aes = FALSE
  ) +
  scale_y_log10() +
  labs(y = "Bout duration (minutes, log scale)", x = NULL) +
  theme_bw(base_size = 140, base_family = "Palatino") +
  theme(legend.position = "none")

#diurnal
forage_diel <- tibble(
  treatment = rep(c("Control", "Human"), each = 2),
  period    = rep(c("Day", "Night"), 2),
  bouts     = c(12, 5, 0, 3)
)

ggplot(forage_diel, aes(x = treatment, y = bouts, fill = period)) +
  geom_col(position = "fill", width = 0.6) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("Day" = "gray70", "Night" = "black")) +
  labs(y = "Proportion of bouts", x = NULL, fill = NULL) +
  theme_bw(base_size = 140, base_family = "Palatino")

#chi square
forage_counts <- tibble(
  treatment = c("Control", "Human"),
  bouts = c(17, 3)
)

ggplot(forage_counts, aes(x = treatment, y = bouts, fill = treatment)) +
  geom_col(width = 0.4) +
  scale_fill_manual(values = cols, guide = "none") +
  scale_x_discrete(
    limits = c("Control", "Human"),
    labels = c(Control = "C", Human = "H")
  ) +
  labs(y = "Foraging Bouts (n)", x = NULL) +
  theme_bw(base_size = 35) +
  theme(
    text = element_text(family = "Palatino"),
    axis.text = element_text(size = 180, color = "black"),
    axis.text.x = element_text(size = 180, color = "black", margin = margin(t = 25)),                         
    axis.title = element_text(size = 180),
    panel.grid = element_line(color = "grey90"),
    panel.border = element_rect(linewidth = 5, color = "black")
  )
