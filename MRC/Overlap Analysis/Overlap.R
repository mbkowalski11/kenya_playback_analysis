library(tidyverse); library(suncalc); library(circular); library(overlap); library(lme4); 
library(patchwork); library(broom.mixed); library(CircStats); library(activity); 
library(purrr); library(broom); library(stringr); library(lubridate); library(hms); library(rlang)

overlap <- read.csv("solaroverlap.csv")
str(overlap)
# Parsing Timestamps ------------------------------------------------------
elephants <- overlap %>% filter(common_name == "African Elephant") %>% group_by(treatment) %>% summarize(n = n())
####Parsing Hell
tz_ke <- "Africa/Nairobi"

## 1) Build a single raw character vector, falling back to start_time if needed
raw_dt <- as.character(overlap$datetime)
raw_st <- as.character(overlap$start_time)
use_st <- is.na(raw_dt) | raw_dt == "" | tolower(raw_dt) == "na"
raw_dt[use_st] <- raw_st[use_st]

## 2) Normalize nasty whitespace and AM/PM artifacts
# replace non-breaking & unicode spaces with normal space, collapse multiples, trim
clean <- raw_dt %>%
  str_replace_all("\\p{Z}+", " ") %>%  # any unicode space
  str_replace_all("[\\u200B\\u200C\\u200D\\uFEFF]", "") %>% # zero-width chars
  str_replace_all("\\s+", " ") %>%
  str_trim()

# remove stray AM/PM (we’ll parse 24h)
clean <- str_replace(clean, "\\s*(AM|PM)\\b", "")

## 3) Split into date part and time part
# Expect either "date time" or just "date"; allow multiple spaces
date_part <- str_match(clean, "^([^ ]+)")[,2]
time_part <- str_match(clean, " ([0-9]{1,2}:[0-9]{2}(:[0-9]{2})?)$")[,2]

# If time missing, default to 00:00:00
time_part[is.na(time_part)] <- "00:00:00"

# Normalize single-digit hour to two digits (e.g., "0:00" -> "00:00")
time_part <- str_replace(time_part, "^(\\d):", "0\\1:")

# If time is HH:MM, add :00 seconds
time_part <- ifelse(str_detect(time_part, "^[0-2]?\\d:\\d{2}$"),
                    paste0(time_part, ":00"), time_part)
time_part <- str_match(clean, " ([0-9]{1,2}:[0-9]{2}(:[0-9]{2})?)$")[,2]
time_part <- as.character(time_part)                # <-- ensure character
time_part[is.na(time_part)] <- "00:00:00"
time_part <- str_replace(time_part, "^(\\d):", "0\\1:")
time_part <- dplyr::if_else(                       # <-- type-stable
  str_detect(time_part, "^\\d{1,2}:\\d{2}$"),
  paste0(time_part, ":00"),
  time_part
)
## 4) Parse date part leniently:
# heuristic: if it has '-' assume YMD; if it has '/' assume MDY (your screenshot shows MDY like 8/3/24)
date_parsed <- ifelse(str_detect(date_part, "-"),
                      as.character(ymd(date_part, quiet = TRUE)),
                      as.character(mdy(date_part, quiet = TRUE)))
date_parsed <- ymd(date_parsed)  # coerce to Date safely

## 5) Parse time part as hms
# hm() handles "HH:MM"; hms() handles "HH:MM:SS"
time_parsed <- suppressWarnings(hms::as_hms(time_part))  # returns NA if truly bad

## 6) Combine to POSIXct in Africa/Nairobi
# Extract components
dt <- make_datetime(
  year = year(date_parsed),
  month = month(date_parsed),
  day = day(date_parsed),
  hour = hour(time_parsed),
  min = minute(time_parsed),
  sec = second(time_parsed),
  tz = tz_ke
)

## 7) Guard: print anything that still failed
if (anyNA(dt)) {
  bad_rows <- which(is.na(dt))
  message("Unparsed rows (showing up to 10):")
  print(head(data.frame(raw = raw_dt[bad_rows], clean = clean[bad_rows]), 10))
  stop("There are still unparsed datetimes. The preview above will show the offending strings.")
}

## 8) Commit and derive your sun columns (same row only)
overlap$datetime     <- dt
overlap$date         <- as_date(dt)

##########
#Solar Times
lat <- 0.2923; lon <- 36.8985

dates <- seq(ymd("2024-07-24"), ymd("2024-10-13"), by = "day")

sun_times <- getSunlightTimes(date = dates, lat = lat,lon = lon,
                              keep = c("sunrise","dawn", "sunset", "dusk"), tz = "Africa/Nairobi")

overlap <- overlap %>% left_join(sun_times, by = "date")

overlap <- overlap %>%
  mutate(
    # 1) Binary daylight (between sunrise and sunset)
    daylight = as.integer(datetime >= dawn & datetime <= dusk),
    # 2) Binary night (between dusk and next dawn — here we treat it as same-day night)
    night = as.integer(datetime < dawn | datetime > dusk),
    # 3) Solar midnight = midpoint between same-row dusk and next dawn
    solar_midnight = dusk + ( (dawn + days(1)) - dusk ) / 2,
    # 4) Hours from solar midnight (signed, decimal)
    # shortest circular distance to solar_midnight in hours (0..12)
    .delta_h = as.numeric(difftime(datetime, solar_midnight, units = "hours")),
    hrs_from_solar_midnight = abs(((.delta_h + 12) %% 24) - 12)) %>%
    select(-.delta_h)

summary(overlap$hrs_from_solar_midnight)

#write.csv(overlap, file = "solaroverlap.csv")

##########
#Calculating Average Sunrise, Sunset, and Solar Midnight at Mpala Research Centre (24/07/25 - 13/10/25)
#Average numeric (hours) sunrise and sunset
# sun_avg <- sun_times %>%
#   mutate(
#     sunrise_hour = hour(sunrise) + minute(sunrise)/60,
#     sunset_hour = hour(sunset) + minute(sunset)/60
#   ) %>%
#   summarise(
#     avg_sunrise = mean(sunrise_hour), #average sunrise = 6.469309 numeric hours
#     avg_sunset = mean(sunset_hour) #average sunset = 18.58679 numeric hours
#   )
# 
# #average solar midnight
# sun_times <- sun_times %>%
#   mutate(sunrise_next = lead(sunrise)) %>%
#   filter(!is.na(sunrise_next))  
# 
# sun_times <- sun_times %>%
#   mutate(
#     solar_midnight = sunset + (sunrise_next - sunset) / 2 )
# 
# avg_solar_midnight <- mean(sun_times$solar_midnight)
# format(avg_solar_midnight, "%H:%M") #solar midnight - 00:32

#Convert detection time to radians
#overlap <- read.csv("recovery.csv")

# overlap <- overlap %>%
#   mutate(start_time = mdy_hm(start_time),
#          hour_decimal = hour(start_time) + minute(start_time) / 60)
# avg_sunrise <- 6.469309 #for study period
# avg_sunset  <- 18.58679 #for study period
# overlap <- overlap %>% 
#   mutate(daylight = if_else(hour_decimal >= avg_sunrise & hour_decimal <= avg_sunset, 1, 0))
# 
# overlap <- overlap %>%
#    mutate(
#      # If no ":" is in start_time, assign midnight (00:00)
#      start_time = if_else(
#        !str_detect(start_time, ":"), 
#        paste0(start_time, " 00:00"), 
#        start_time
#      ),
#      
#      # Parse datetime
#      datetime = parse_date_time(
#        start_time,
#        orders = c("mdy HM", "mdy HMS", "dmy HM", "dmy HMS", "ymd HM", "ymd HMS"),
#        tz = "Africa/Nairobi"
#      ),
#      
#      # Convert to radians for activity analysis
#      time_rad = ((hour(datetime) + minute(datetime)/60 + second(datetime)/3600) / 24) * 2 * pi
#    )

##########
#Two-tailed Fisher Exact Test
# ----- Focal species list -----
focal_species <- c(
  "African Elephant", "Guenther's Dik-Dik", "Hares", "Hippopotamus",
  "Spotted Hyena", "Genets", "Impala", "Leopard", "White-tailed Mongoose",
  "Black-backed Jackal", "Grevy's Zebra", "Reticulated Giraffe", "Zorilla", 
  "Ground Squirrels", "Plains Zebra", "Striped Hyena", "Common Warthog", "Lion",
  "Slender Mongoose"
)

# ----- Pick species column name automatically -----
species_col <- if ("common_name" %in% names(overlap)) "common_name" else "common"

# ----- Clean up fields -----
overlap2 <- overlap %>%
  mutate(
    # standardize treatment to control/human
    treatment = tolower(as.character(treatment)),
    treatment = ifelse(treatment %in% c("control","human"), treatment, treatment),
    treatment = factor(treatment, levels = c("control","human")),
    # coerce daylight to 0/1 numeric
    daylight = case_when(
      is.numeric(daylight) ~ as.integer(daylight > 0),
      is.logical(daylight) ~ as.integer(daylight),
      TRUE ~ as.integer(tolower(as.character(daylight)) %in% c("1","day","daylight","true","yes"))
    )
  ) %>%
  filter(.data[[species_col]] %in% focal_species)

# ----- Per-species Fisher's test (two-sided) -----
fisher_daylight <- overlap2 %>%
  group_split(.data[[species_col]], .keep = TRUE) %>%
  map_df(function(df) {
    sp <- df[[species_col]][1]
    
    # Build a 2x2 contingency table: rows=treatment (control/human), cols=daylight (0/1)
    tab <- table(
      factor(df$treatment, levels = c("control","human")),
      factor(df$daylight,   levels = c(0,1))
    )
    
    # Validate table: must be 2x2, each row/column has >0 total
    valid <- identical(dim(tab), c(2L, 2L)) &&
      all(rowSums(tab) > 0) &&
      all(colSums(tab) > 0)
    
    if (valid) {
      ft <- fisher.test(tab, alternative = "two.sided")
      tibble(
        !!species_col := sp,
        control_night = as.integer(tab["control","0"]),
        control_day   = as.integer(tab["control","1"]),
        human_night   = as.integer(tab["human","0"]),
        human_day     = as.integer(tab["human","1"]),
        odds_ratio    = unname(ft$estimate),
        p_value       = ft$p.value
      )
    } else {
      tibble(
        !!species_col := sp,
        control_night = NA_integer_, control_day = NA_integer_,
        human_night   = NA_integer_, human_day   = NA_integer_,
        odds_ratio    = NA_real_,    p_value     = NA_real_
      )
    }
  }) %>%
  arrange(p_value)

print(fisher_daylight)

##########
#Overlap Analysis
overlap <- read.csv("solaroverlap.csv")

focal_species <- c("African Elephant", "Guenther's Dik-Dik", "Hares", "Hippopotamus",
                   "Spotted Hyena", "Genets", "Impala", "Leopard", "White-tailed Mongoose",
                   "Black-backed Jackal", "Grevy's Zebra", "Reticulated Giraffe", "Zorilla", 
                   "Ground Squirrels", "Plains Zebra", "Striped Hyena", "Common Warthog", "Lion",  
                   "Slender Mongoose")

overlap.data <- overlap %>%
  filter(common_name %in% focal_species)

compute_overlap <- function(species, df, n_boot = 10000) {
  
  cat("Processing:", species, "... ")
  start_time <- Sys.time()
  
  dat <- df %>% 
    filter(common_name == species)
  
  t1 <- dat %>% filter(treatment == "control") %>% pull(time_rad)
  t2 <- dat %>% filter(treatment == "human")   %>% pull(time_rad)
  
  n1 <- length(t1); n2 <- length(t2)
  if (n1 == 0 || n2 == 0) {
    cat("SKIPPED (no data)\n")
    return(tibble(
      species, n_control = n1, n_human = n2,
      estimator = NA_character_, Dhat = NA_real_,
      CI_low = NA_real_, CI_high = NA_real_,
      temporal_shift = NA
    ))
  }
  
  type <- if (min(n1, n2) < 75) "Dhat1" else "Dhat4"
  
  Dhat <- overlapEst(t1, t2, type = type)
  CI   <- quantile(bootstrap(t1, t2, nb = n_boot, type = type),
                   c(0.025, 0.975))
  
  end_time <- Sys.time()
  elapsed <- round(as.numeric(end_time - start_time), 1)
  cat("DONE (", elapsed, "s, n1=", n1, ", n2=", n2, ")\n", sep="")
  
  return(tibble(
    species,
    n_control = n1,
    n_human   = n2,
    estimator = type,
    Dhat      = round(Dhat, 3),
    CI_low    = round(CI[1], 3),
    CI_high   = round(CI[2], 3),
    temporal_shift = Dhat < 0.90
  ))
}

overlap.results <- map_dfr(focal_species, compute_overlap, df = overlap.data)

write.csv(overlap.results, "overlap.results.csv")

###########THIS IS WHERE THE PERMUTATION TEST IS!
perm_test_overlap <- function(species, df, reps = 10000) {
  dat <- df %>% filter(common_name == species)
  t_ctrl  <- dat %>% filter(treatment == "control") %>% pull(time_rad)
  t_human <- dat %>% filter(treatment == "human")   %>% pull(time_rad)
  n1 <- length(t_ctrl); n2 <- length(t_human)
  f_ctrl  <- fitact(t_ctrl); f_human <- fitact(t_human)
  ck <- compareCkern(f_ctrl, f_human, reps = reps)
  tibble(
    species,
    n_control = n1,
    n_human   = n2,
    Dhat_obs  = round(as.numeric(ck["obs"      ]), 3),
    null_mean = round(as.numeric(ck["null"     ]), 3),
    se_null   = round(as.numeric(ck["seNull"   ]), 3),
    p_overlap = round(as.numeric(ck["pNull"    ]), 4)
  )
}

p.results <- map_dfr(focal_species, perm_test_overlap, df = overlap.data)
write.csv(p.results, "compareCkern.csv")

overlap_full <- overlap.results %>% 
  left_join(p.results %>% select(species, p_overlap),
            by = "species")


##########Group Overlap and CompareCKern
df <- overlap_subset
focal_groups <- unique(df$group)

compute_overlap_group <- function(grp, nb = 10000) {
  cat("Processing:", grp, "... ", flush = TRUE)
  dat <- df %>% filter(group == grp)
  
  t1 <- dat %>% filter(treatment == "control") %>% pull(time_rad)
  t2 <- dat %>% filter(treatment == "human")   %>% pull(time_rad)
  
  n1 <- length(t1); n2 <- length(t2)
  if (!n1 || !n2) {
    cat("SKIPPED (no data)\n", flush = TRUE)
    return(tibble(group = grp, n_control=n1, n_human=n2,
                  estimator = NA, Dhat = NA, CI_low = NA, CI_high = NA,
                  temporal_shift = NA))
  }
  
  type <- if (min(n1, n2) < 75) "Dhat1" else "Dhat4"
  Dhat <- overlap::overlapEst(t1, t2, type = type)
  bs   <- overlap::bootstrap(t1, t2, nb = nb, type = type)
  CI   <- quantile(bs, c(0.025, 0.975), na.rm = TRUE)
  
  cat("DONE (n1=", n1, ", n2=", n2, ")\n", sep = "", flush = TRUE)
  
  tibble(
    group = grp, n_control = n1, n_human = n2,
    estimator = type,
    Dhat = round(Dhat, 3),
    CI_low  = round(unname(CI[1]), 3),
    CI_high = round(unname(CI[2]), 3),
    temporal_shift = Dhat < 0.90
  )
}

overlap.results.group <- map_dfr(focal_groups, compute_overlap_group)
write.csv(overlap.results.group, "overlap.results.group.csv", row.names = FALSE)


# -------------------------------------------------------------------
#  Permutation (compareCkern) test
# -------------------------------------------------------------------
perm_test_overlap <- function(species, df, reps = 10000) {
  dat <- df %>% filter(group == species)
  
  t_ctrl  <- dat %>% filter(treatment == "control") %>% pull(time_rad)
  t_human <- dat %>% filter(treatment == "human")   %>% pull(time_rad)
  
  n1 <- length(t_ctrl)
  n2 <- length(t_human)
  
  if (n1 == 0 || n2 == 0) {
    return(tibble(
      species,
      group = unique(dat$group)[1] %||% NA_character_,
      n_control = n1, n_human = n2,
      Dhat_obs  = NA_real_,
      null_mean = NA_real_,
      se_null   = NA_real_,
      p_overlap = NA_real_
    ))
  }
  
  f_ctrl  <- fitact(t_ctrl)
  f_human <- fitact(t_human)
  ck      <- compareCkern(f_ctrl, f_human, reps = reps)
  
  tibble(
    species,
    group     = unique(dat$group)[1] %||% NA_character_,
    n_control = n1,
    n_human   = n2,
    Dhat_obs  = round(as.numeric(ck["obs"     ]), 3),
    null_mean = round(as.numeric(ck["null"    ]), 3),
    se_null   = round(as.numeric(ck["seNull"  ]), 3),
    p_overlap = round(as.numeric(ck["pNull"   ]), 4)
  )
}

p.results <- map_dfr(focal_groups, perm_test_overlap, df = df)
write.csv(p.results, "compareCkern.group.csv", row.names = FALSE)


##########New overlap loop
one_overlap <- function(species, df,
                        nb_boot   = 10000,   # bootstrap reps for CI
                        reps_perm = 10000,   # permutation reps for p-value
                        min_n     = 10) {    # skip if < min_n detections per group
  
  cat("Processing:", species, "... ")
  start_time <- Sys.time()
  
  dat <- df %>% filter(common_name == species)
  
  t_ctrl  <- dat %>% filter(treatment == "control") %>% pull(time_rad)
  t_human <- dat %>% filter(treatment == "human")   %>% pull(time_rad)
  
  n1 <- length(t_ctrl)
  n2 <- length(t_human)
  
  if (n1 < min_n | n2 < min_n) {
    cat("SKIPPED (n1 =", n1, ", n2 =", n2, ")\n")
    return(tibble(
      species, n_control = n1, n_human = n2,
      estimator = NA_character_, Dhat = NA_real_,
      CI_low = NA_real_, CI_high = NA_real_,
      null_mean = NA_real_, se_null = NA_real_,
      p_overlap = NA_real_
    ))
  }
  
  estimator <- "Dhat4"
  
  # point estimate
  Dhat <- overlapEst(t_ctrl, t_human, type = estimator)
  
  # bootstrap CI
  CI <- quantile(
    bootstrap(t_ctrl, t_human, nb = nb_boot, type = estimator),
    c(0.025, 0.975))
  
  # permutation test (left-tailed: expecting overlap to decrease)
  pooled <- c(t_ctrl, t_human)
  perm_vals <- replicate(reps_perm, {
    idx <- sample.int(n1 + n2)
    overlapEst(pooled[idx[1:n1]],
               pooled[idx[(n1 + 1):(n1 + n2)]],
               type = estimator)
  })
  
  null_mean <- mean(perm_vals)
  se_null   <- sd(perm_vals)
  p_overlap <- mean(perm_vals != Dhat)            # two-tailed
  
  cat("DONE (", round(as.numeric(Sys.time() - start_time), 1),
      "s, n1=", n1, ", n2=", n2, ")\n", sep = "")
  
  tibble(
    species,
    n_control = n1,
    n_human   = n2,
    estimator = estimator,
    Dhat      = round(Dhat, 3),
    CI_low    = round(CI[1], 3),
    CI_high   = round(CI[2], 3),
    null_mean = round(null_mean, 3),
    se_null   = round(se_null, 3),
    p_overlap = round(p_overlap, 4)
  )
}

twotailed.results <- map_dfr(focal_species, one_overlap, df = overlap.data)
write.csv(left.results, "all.KDEperm.csv")

##########
#KDE Bootstrapping 
# -------- settings --------
min_n   <- 10          # minimum detections in *each* treatment
reps    <- 10000       # bootstrap reps for fitact()
overlap_reps <- 10000  # bootstrap reps for overlap analysis

# helper to convert a fitact object → tidy tibble
fit_to_df <- function(fit, treatment, species) {
  tibble(time_rad = fit@pdf[, 1],
         density  = fit@pdf[, 2]) |>
    mutate(time_hr  = (time_rad %% (2 * pi)) * 24 / (2 * pi),
           treatment = treatment,
           species   = species)
}

#bootstrap loop
all_kde <- map_dfr(focal_species, function(sp) {
  dat <- overlap |>
    filter(common_name == sp)
  
  # must have control & human, each ≥ min_n records
  if (!all(c("control", "human") %in% dat$treatment)) return(NULL)
  if (min(table(dat$treatment)) < min_n)              return(NULL)
  
  ctrl  <- dat |> filter(treatment == "control") |> pull(time_rad)
  human <- dat |> filter(treatment == "human")   |> pull(time_rad)
  
  fc <- fitact(ctrl,  sample = "data", reps = reps)
  fh <- fitact(human, sample = "data", reps = reps)
  
  bind_rows(
    fit_to_df(fc, "Control", sp),
    fit_to_df(fh, "Human",   sp)
  )
})
write.csv(all_kde, "all.kde.csv")

all_kde <- read.csv(("all.kde.csv"))

cols <- c("Control" = "#FABA39FF", "Human"   = "#7A0403FF")

spec_order <- c(
  # Group 1
  "African Elephant","Lion", "Leopard", "Black-backed Jackal", "Zorilla", 
  # Group 2
  "Hippopotamus", "Reticulated Giraffe", "Plains Zebra", "Grevy's Zebra", "Common Warthog",  
  # Group 3
  "Spotted Hyena", "Striped Hyena", "Genets", "White-tailed Mongoose", "Hares",
  # Group 4
  "Ground Squirrels", "Impala", "Guenther's Dik-Dik",  "Slender Mongoose"
)

all_kde$species <- factor(all_kde$species, levels = spec_order)

rug_ctrl  <- overlap |>
  filter(common_name %in% spec_order, treatment == "control") |>
  transmute(
    species   = factor(common_name, levels = spec_order),
    time_hr   = (time_rad %% (2*pi))*24/(2*pi),
    y_start   = 0,
    y_end     = -0.04         
  )

rug_human <- overlap |>
  filter(common_name %in% spec_order, treatment == "human") |>
  transmute(
    species   = factor(common_name, levels = spec_order),
    time_hr   = (time_rad %% (2*pi))*24/(2*pi),
    y_start   = -0.041,
    y_end     = -0.08           
  )

#Composite Plot
ggplot(all_kde,
       aes(x = time_hr, y = density, group = treatment)) +
  geom_ribbon(aes(ymin = 0, ymax = density, fill = treatment),
              alpha = 0.3, colour = NA) +
  geom_line(aes(colour = treatment), linewidth = 8) +
  geom_segment(
    data = rug_ctrl,
    aes(x = time_hr,  xend = time_hr,
        y = y_start, yend = y_end),
    inherit.aes = FALSE,
    linewidth   = 2,
    colour      = cols["Control"],
    alpha       = 0.9
  ) +
  geom_segment(
    data = rug_human,
    aes(x = time_hr,  xend = time_hr,
        y = y_start, yend = y_end),
    inherit.aes = FALSE,
    linewidth   = 2,
    colour      = cols["Human"],  
    alpha       = 0.9
  ) +
  scale_fill_manual(values = cols) +
  scale_colour_manual(values = cols) +
  scale_x_continuous(breaks = c(0, 6, 12, 18, 24), limits = c(0, 24)) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    limits = c(-0.1, NA)) +
  labs(x = "Time (Hours)", y = "Kernel Density") +
  facet_wrap(~ species, ncol = 4, dir = "v") +  
  theme_bw(base_size = 180, base_family = "Palatino") +
  theme(
    legend.position   = "none",
    strip.background  = element_rect(fill = "white", color = NA),
    strip.text        = element_text(face = "bold", size = 100),
    axis.title        = element_text(size = 180, colour = "black"),
    axis.text         = element_text(size = 140, colour = "black"),
    panel.grid.major  = element_line(color = "gray90", linewidth = 4),
    panel.grid.minor  = element_blank(),
    axis.text.x = element_text(size = 140),
    axis.ticks.x = element_line(),
  )


overlap <- read.csv("solaroverlap.csv")
all_kde_group <- read.csv(("all.kde.group.csv"))

cols <- c("Control" = "#FABA39FF", "Human"   = "#7A0403FF")

spec_to_group <- tribble(
  ~common_name,              ~group,
  "African Elephant",        "African Elephant",
  "Lion",                    "Large Felids",
  "Leopard",                 "Large Felids",
  "Plains Zebra",            "Zebras",
  "Grevy's Zebra",           "Zebras",
  "Black-backed Jackal",     "Mesocarnivores",
  "Zorilla",                 "Mesocarnivores",
  "Hippopotamus",            "Hippopotamus"
)

spec_order_grp <- c("Large Felids", "Zebras","Mesocarnivores", "African Elephant", "Hippopotamus")

rug_ctrl <- overlap %>%
  left_join(spec_to_group, by = "common_name") %>%
  filter(!is.na(group), treatment == "control") %>%
  transmute(
    group   = factor(group, levels = spec_order_grp),
    time_hr = (time_rad %% (2*pi)) * 24/(2*pi),
    y_start = 0,
    y_end   = -0.04
  )

rug_human <- overlap %>%
  left_join(spec_to_group, by = "common_name") %>%
  filter(!is.na(group), treatment == "human") %>%
  transmute(
    group   = factor(group, levels = spec_order_grp),
    time_hr = (time_rad %% (2*pi)) * 24/(2*pi),
    y_start = -0.041,
    y_end   = -0.08
  )

all_kde_group$group <- factor(all_kde_group$group, levels = spec_order_grp)

ggplot(all_kde_group,
       aes(x = time_hr, y = density, group = treatment)) +
  geom_ribbon(aes(ymin = 0, ymax = density, fill = treatment),
              alpha = 0.3, colour = NA) +
  geom_line(aes(colour = treatment), linewidth = 8) +
  geom_segment(
    data = rug_ctrl,
    aes(x = time_hr,  xend = time_hr,
        y = y_start, yend = y_end),
    inherit.aes = FALSE,
    linewidth   = 2,
    colour      = cols["Control"],
    alpha       = 0.9
  ) +
  geom_segment(
    data = rug_human,
    aes(x = time_hr,  xend = time_hr,
        y = y_start, yend = y_end),
    inherit.aes = FALSE,
    linewidth   = 2,
    colour      = cols["Human"],  
    alpha       = 0.9
  ) +
  scale_fill_manual(values = cols) +
  scale_colour_manual(values = cols) +
  scale_x_continuous(breaks = c(0, 6, 12, 18, 24), limits = c(0, 24)) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    limits = c(-0.1, NA)) +
  labs(x = "Time (Hours)", y = "Kernel Density") +
  facet_wrap(~ group, ncol = 2, dir = "v") +  
  theme_bw(base_size = 180, base_family = "Palatino") +
  theme(
    legend.position   = "none",
    strip.background  = element_rect(fill = "white", color = NA),
    strip.text        = element_text(face = "bold", size = 100),
    axis.title        = element_text(size = 180, colour = "black"),
    axis.text         = element_text(size = 160, colour = "black"),
    panel.grid.major  = element_line(color = "gray90", linewidth = 4),
    panel.grid.minor  = element_blank(),
    axis.text.x = element_text(size = 140),
    axis.ticks.x = element_line(),
  )

########
#Elephant only plot
overlap <- read.csv("solaroverlap.csv")
all_kde_group <- read.csv(("all.kde.group.csv"))

cols <- c("Control" = "#FABA39FF", "Human"   = "#7A0403FF")

spec_to_group <- tribble(
  ~common_name,              ~group,
  "African Elephant",        "African Elephant",
  "Lion",                    "Large Felids",
  "Leopard",                 "Large Felids",
  "Plains Zebra",            "Zebras",
  "Grevy's Zebra",           "Zebras",
  "Black-backed Jackal",     "Mesocarnivores",
  "Zorilla",                 "Mesocarnivores",
  "Hippopotamus",            "Hippopotamus"
)

spec_order_grp <- c("Large Felids", "Zebras","Mesocarnivores", "African Elephant", "Hippopotamus")

rug_ctrl <- overlap %>%
  left_join(spec_to_group, by = "common_name") %>%
  filter(!is.na(group), treatment == "control") %>%
  transmute(
    group   = factor(group, levels = spec_order_grp),
    time_hr = (time_rad %% (2*pi)) * 24/(2*pi),
    y_start = 0,
    y_end   = -0.04
  )

rug_human <- overlap %>%
  left_join(spec_to_group, by = "common_name") %>%
  filter(!is.na(group), treatment == "human") %>%
  transmute(
    group   = factor(group, levels = spec_order_grp),
    time_hr = (time_rad %% (2*pi)) * 24/(2*pi),
    y_start = -0.041,
    y_end   = -0.08
  )

all_kde_group$group <- factor(all_kde_group$group, levels = spec_order_grp)

ggplot(all_kde_group,
       aes(x = time_hr, y = density, group = treatment)) +
  geom_ribbon(aes(ymin = 0, ymax = density, fill = treatment),
              alpha = 0.3, colour = NA) +
  geom_line(aes(colour = treatment), linewidth = 8) +
  geom_segment(
    data = rug_ctrl,
    aes(x = time_hr,  xend = time_hr,
        y = y_start, yend = y_end),
    inherit.aes = FALSE,
    linewidth   = 2,
    colour      = cols["Control"],
    alpha       = 0.9
  ) +
  geom_segment(
    data = rug_human,
    aes(x = time_hr,  xend = time_hr,
        y = y_start, yend = y_end),
    inherit.aes = FALSE,
    linewidth   = 2,
    colour      = cols["Human"],  
    alpha       = 0.9
  ) +
  scale_fill_manual(values = cols) +
  scale_colour_manual(values = cols) +
  scale_x_continuous(breaks = c(0, 6, 12, 18, 24), limits = c(0, 24)) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    limits = c(-0.1, NA)) +
  labs(x = "Time (Hours)", y = "Kernel Density") +
  facet_wrap(~ group, ncol = 2, dir = "v") +  
  theme_bw(base_size = 180, base_family = "Palatino") +
  theme(
    legend.position   = "none",
    strip.background  = element_rect(fill = "white", color = NA),
    strip.text        = element_text(face = "bold", size = 100),
    axis.title        = element_text(size = 180, colour = "black"),
    axis.text         = element_text(size = 160, colour = "black"),
    panel.grid.major  = element_line(color = "gray90", linewidth = 4),
    panel.grid.minor  = element_blank(),
    axis.text.x = element_text(size = 140),
    axis.ticks.x = element_line(),
  )
##########
#elephant only plot
overlap <- read.csv("solaroverlap.csv")
all_kde_group <- read.csv("all.kde.group.csv")

cols <- c("Control" = "#FABA39FF", "Human" = "#7A0403FF")

# Elephant-only rugs
rug_ctrl_ele <- overlap %>%
  filter(common_name == "African Elephant", treatment == "control") %>%
  transmute(
    time_hr = (time_rad %% (2*pi)) * 24/(2*pi),
    y_start = 0,
    y_end   = -0.02   # was -0.04
  )

rug_human_ele <- overlap %>%
  filter(common_name == "African Elephant", treatment == "human") %>%
  transmute(
    time_hr = (time_rad %% (2*pi)) * 24/(2*pi),
    y_start = -0.021,
    y_end   = -0.04   # was -0.08
  )

# Elephant-only KDE
kde_ele <- all_kde_group %>%
  filter(group == "African Elephant")

ggplot(kde_ele, aes(x = time_hr, y = density, group = treatment)) +
  geom_ribbon(aes(ymin = 0, ymax = density, fill = treatment),
              alpha = 0.3, colour = NA) +
  geom_line(aes(colour = treatment), linewidth = 8) +
  geom_segment(
    data = rug_ctrl_ele,
    aes(x = time_hr, xend = time_hr, y = y_start, yend = y_end),
    inherit.aes = FALSE,
    linewidth = 2,
    colour = cols["Control"],
    alpha = 0.9
  ) +
  geom_segment(
    data = rug_human_ele,
    aes(x = time_hr, xend = time_hr, y = y_start, yend = y_end),
    inherit.aes = FALSE,
    linewidth = 2,
    colour = cols["Human"],
    alpha = 0.9
  ) +
  #geom_hline(yintercept = 0, linewidth = 8) +
  scale_fill_manual(values = cols) +
  scale_colour_manual(values = cols) +
  scale_x_continuous(breaks = c(0, 6, 12, 18, 24), limits = c(0, 24)) +
  scale_y_continuous(
    limits = c(-0.05, NA),                 # room for rugs
    breaks = scales::pretty_breaks(n = 4), # nice ticks
    labels = function(x) ifelse(x < 0, "", x),  # hide negatives
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(x = "Time (Hours)", y = "Kernel Density") +
  theme_bw(base_size = 180, base_family = "Palatino") +
  theme(
    legend.position  = "none",
    axis.title       = element_text(size = 180, colour = "black"),
    axis.text        = element_text(size = 160, colour = "black"),
    panel.grid.major = element_line(color = "gray90", linewidth = 4),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(size = 140),
    axis.ticks.x     = element_line()
  )

##############
#Group KDE Bootstrapping 
# -------- settings (reuse yours) --------
reps          <- 10000

groups <- list(
  "African Elephant" = c("African Elephant"),
  "Hippopotamus"     = c("Hippopotamus"),
  "Large Felids"     = c("Lion","Leopard"),
  "Zebras"           = c("Plains Zebra","Grevy's Zebra"),
  "Mesocarnivores"   = c("Black-backed Jackal","Zorilla")
)
group_levels <- names(groups)

reps <- 10000

fit_to_df <- function(fit, treatment, grp) {
  tibble(time_rad = fit@pdf[,1], density = fit@pdf[,2]) |>
    mutate(time_hr = (time_rad %% (2*pi))*24/(2*pi),
           treatment = treatment,
           group = grp)
}

get_group_times <- function(group_species, data = overlap) {
  dat <- dplyr::filter(data, common_name %in% group_species)
  list(
    control = dplyr::filter(dat, treatment == "control") |> dplyr::pull(time_rad),
    human   = dplyr::filter(dat, treatment == "human")   |> dplyr::pull(time_rad)
  )
}

# -------------------- KDE build (from the SAME mapping) --------------------
all_kde_groups <- purrr::map_dfr(group_levels, function(grp){
  sp_vec <- groups[[grp]]
  times  <- get_group_times(sp_vec, overlap)
  ctrl_times  <- times$control
  human_times <- times$human
  if (!length(ctrl_times) || !length(human_times)) return(NULL)
  
  fc <- activity::fitact(ctrl_times,  sample = "data", reps = reps)
  fh <- activity::fitact(human_times, sample = "data", reps = reps)
  
  dplyr::bind_rows(
    fit_to_df(fc, "Control", grp),
    fit_to_df(fh, "Human",   grp)
  )
})

# ---- write/read if you want parity with your species workflow ----
write.csv(all_kde_groups, "all.kde.group.csv", row.names = FALSE)
# all_kde_groups <- read.csv("all.kde.groups.csv")

# -------- rugs for grouped plot (pooled across species) --------
# 1) Rename to avoid conflict with dplyr::group_map()
grp_map <- list(
  `African Elephant`   = c("African Elephant"),
  `Hippopotamus`       = c("Hippopotamus"),
  `Large Felids`       = c("Lion","Leopard"),
  `Zebras`    = c("Common Warthog","Plains Zebra","Grevy's Zebra"),
  `Mesocarnivores`     = c("Black-backed Jackal","Zorilla")
)
stopifnot(is.list(grp_map))

# 2) Tidy species→group lookup (one row per species)
sp2grp <- tibble::enframe(grp_map, name = "group", value = "species") |>
  tidyr::unnest_longer(species) |>      # use unnest(species) if tidyr<1.3
  mutate(species = as.character(species))

# 3) Attach groups to overlap once
overlap_g <- overlap |>
  mutate(common_name = as.character(common_name)) |>
  inner_join(sp2grp, by = c("common_name" = "species"))

group_levels <- names(grp_map)

# 4) Rugs (now simple & robust)
rug_ctrl_g <- overlap_g |>
  filter(treatment == "control") |>
  transmute(
    group   = factor(group, levels = group_levels),
    time_hr = (time_rad %% (2*pi))*24/(2*pi),
    y_start = 0,
    y_end   = -0.04
  )

rug_human_g <- overlap_g |>
  filter(treatment == "human") |>
  transmute(
    group   = factor(group, levels = group_levels),
    time_hr = (time_rad %% (2*pi))*24/(2*pi),
    y_start = -0.041,
    y_end   = -0.08
  )

# -------- plot (facet by group) --------
cols <- c("Control" = "#FABA39FF", "Human" = "#7A0403FF")
group_levels <- c("Large Felids",  "Zebras", "Mesocarnivores", "African Elephant","Hippopotamus")

ggplot(all_kde_groups, aes(x = time_hr, y = density, group = treatment)) +
  geom_ribbon(aes(ymin = 0, ymax = density, fill = treatment), alpha = 0.3, colour = NA) +
  geom_line(aes(colour = treatment), linewidth = 8) +
  geom_segment(
    data = rug_ctrl_g,
    aes(x = time_hr, xend = time_hr, y = y_start, yend = y_end),
    inherit.aes = FALSE, linewidth = 2, colour = cols["Control"], alpha = 0.9
  ) +
  geom_segment(
    data = rug_human_g,
    aes(x = time_hr, xend = time_hr, y = y_start, yend = y_end),
    inherit.aes = FALSE, linewidth = 2, colour = cols["Human"], alpha = 0.9
  ) +
  scale_fill_manual(values = cols) +
  scale_colour_manual(values = cols) +
  scale_x_continuous(breaks = c(0, 6, 12, 18, 24), limits = c(0, 24)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(-0.1, NA)) +
  labs(x = "Time (Hours)", y = "Kernel Density") +
  facet_wrap(~ factor(group, levels = group_levels), ncol = 2, dir = "v") +
  theme_bw(base_size = 180, base_family = "Palatino") +
  theme(
    legend.position   = "none",
    strip.background  = element_rect(fill = "white", color = NA),
    strip.text        = element_text(face = "bold", size = 100),
    axis.title        = element_text(size = 180, colour = "black"),
    axis.text         = element_text(size = 140, colour = "black"),
    panel.grid.major  = element_line(color = "gray90", linewidth = 4),
    panel.grid.minor  = element_blank(),
    axis.text.x       = element_text(size = 140),
    axis.ticks.x      = element_line()
  )
##########
#Circular Statistics
# Vector of species to test
focal_species <- c("African Elephant", "Guenther's Dik-Dik", "Hares",
                   "Hippopotamus", "Spotted Hyena", "Genets", "Impala",
                   "Leopard", "White-tailed Mongoose", "Black-backed Jackal",
                   "Grevy's Zebra", "Reticulated Giraffe", "Zorilla",
                   "Ground Squirrels", "Plains Zebra", "Striped Hyena",
                   "Common Warthog", "Lion", "Slender Mongoose")

# Function to calculate Watson U² p-value
watson_pvalue <- function(U2, n1, n2) {
  # For large samples, Watson U² follows an asymptotic distribution
  # The critical values and p-value calculation are based on:
  # Watson, G.S. (1962). "Goodness-of-fit tests on a circle"
  
  # Effective sample size
  n_eff <- (n1 * n2) / (n1 + n2)
  
  # Standardized test statistic
  # For large samples, use the asymptotic distribution
  if (n_eff >= 17) {
    # Large sample approximation
    # P-value based on asymptotic distribution
    if (U2 < 0.051) {
      pval <- 1.0
    } else if (U2 < 0.084) {
      pval <- 0.5
    } else if (U2 < 0.136) {
      pval <- 0.2
    } else if (U2 < 0.187) {
      pval <- 0.1
    } else if (U2 < 0.267) {
      pval <- 0.05
    } else if (U2 < 0.347) {
      pval <- 0.02
    } else if (U2 < 0.461) {
      pval <- 0.01
    } else if (U2 < 0.743) {
      pval <- 0.001
    } else {
      pval <- 0.0001  # Very small p-value
    }
  } else {
    # Small sample - use more conservative approach
    # These are approximate critical values
    if (U2 < 0.1) {
      pval <- 0.5
    } else if (U2 < 0.15) {
      pval <- 0.2
    } else if (U2 < 0.2) {
      pval <- 0.1
    } else if (U2 < 0.3) {
      pval <- 0.05
    } else if (U2 < 0.4) {
      pval <- 0.01
    } else {
      pval <- 0.001
    }
  }
  
  return(pval)
}

# Updated function with manual p-value calculation
watson_test_loop <- function(species, df) {
  cat("Running Watson test for:", species, "\n")
  
  dat <- df %>% filter(.data[["common_name"]] == species)
  t_ctrl  <- dat %>% filter(treatment == "control") %>% pull(time_rad)
  t_human <- dat %>% filter(treatment == "human")   %>% pull(time_rad)
  
  n1 <- length(t_ctrl)
  n2 <- length(t_human)
  cat("  n_control =", n1, "; n_human =", n2, "\n")
  
  t_ctrl_circ  <- circular(t_ctrl,  units = "radians", modulo = "2pi")
  t_human_circ <- circular(t_human, units = "radians", modulo = "2pi")
  
  if (!is.null(out)) {
    stat <- out$statistic
    # Calculate p-value manually
    pval <- watson_pvalue(stat, n1, n2)
    cat("  ✅ Watson test succeeded: U2 =", round(stat, 4), ", p ≈", pval, "\n")
  } else {
    cat("  ❌ Watson test failed for", species, "\n")
    stat <- NA_real_
    pval <- NA_real_
  }
  
  data.frame(
    species = species,
    n_control = n1,
    n_human = n2,
    U2_stat = round(stat, 4),
    p_value = pval,
    stringsAsFactors = FALSE
  )
}

circular_tests <- function(sp, df, reps = 10000) {
  print(paste("Processing species:", sp))
  
  # Filter data for the species
  dat <- df %>% 
    filter(common_name == sp) %>%
    filter(!is.na(time_rad))
  
  # Create circular objects
  ctrl_data <- dat$time_rad[dat$treatment == "control"]
  hum_data <- dat$time_rad[dat$treatment == "human"]
  
  # Skip if one group is empty
  if (length(ctrl_data) == 0 | length(hum_data) == 0) {
    return(tibble(
      species = sp, 
      watson_W = NA, 
      watson_p = NA, 
      ww_F = NA, 
      ww_p = NA, 
      overlap_est = NA,
      overlap_p = NA
    ))
  }
  
  ctrl <- circular(ctrl_data)
  hum <- circular(hum_data)
  
  # 1. Watson's U² test (fix: extract correct values)
  watson_result <- tryCatch({
    result <- watson.two.test(ctrl, hum)
    # Watson's test doesn't return p.value directly, calculate significance
    # Based on critical values: >0.187 (p<0.05), >0.267 (p<0.01)
    p_val <- if (result$statistic > 0.267) {
      "< 0.01"
    } else if (result$statistic > 0.187) {
      "< 0.05"
    } else {
      "> 0.05"
    }
    
    list(statistic = result$statistic, p.value = p_val)
  }, error = function(e) {
    print(paste("Watson two-sample test error:", e$message))
    # Try alternative: Kuiper's test
    tryCatch({
      result <- kuiper.test(ctrl, hum)
      list(statistic = result$statistic, p.value = result$p.value)
    }, error = function(e2) {
      print(paste("Kuiper test error:", e2$message))
      list(statistic = NA, p.value = NA)
    })
  })
  
  # 2. Watson-Williams test (fix: combine into single vector with grouping)
  ww_result <- tryCatch({
    # Combine data and create grouping factor
    combined_data <- circular(c(ctrl_data, hum_data))
    groups <- factor(c(rep("control", length(ctrl_data)), rep("human", length(hum_data))))
    
    result <- watson.williams.test(combined_data, groups)
    list(statistic = result$statistic, p.value = result$p.value)
  }, error = function(e) {
    print(paste("Watson-Williams test error:", e$message))
    list(statistic = NA, p.value = NA)
  })
  
  # 3. Overlap coefficient and bootstrap test (simplified)
  overlap_result <- tryCatch({
    # Convert to numeric (0 to 2π)
    ctrl_num <- as.numeric(ctrl)
    hum_num <- as.numeric(hum)
    
    # Calculate overlap coefficient - use simpler approach
    overlap_coef <- overlapEst(ctrl_num, hum_num)
    
    # Try bootstrap test with simpler parameters
    bootstrap_result <- tryCatch({
      overlapTrue(ctrl_num, hum_num)
    }, error = function(e) {
      print(paste("Bootstrap test failed:", e$message))
      NA
    })
    
    list(overlap = overlap_coef, p.value = bootstrap_result)
  }, error = function(e) {
    print(paste("Overlap test error:", e$message))
    list(overlap = NA, p.value = NA)
  })
  
  # Return results (ensure single row)
  result_tibble <- tibble(
    species = sp,
    watson_W = as.numeric(watson_result$statistic)[1],  # Take first value only
    watson_p = watson_result$p.value,                   # Keep as character for significance
    ww_F = as.numeric(ww_result$statistic)[1],          # Take first value only
    ww_p = as.numeric(ww_result$p.value)[1],            # Take first value only
    overlap_est = as.numeric(overlap_result$overlap)[1], # Take first value only
    overlap_p = as.numeric(overlap_result$p.value)[1]    # Take first value only
  )
  
  print("Results:")
  print(result_tibble)
  return(result_tibble)
}

# Test with a single species first
print("=== TESTING FIXED FUNCTION ===")
test_result <- circular_tests(focal_species[1], overlap.data)
print("Final test result:")
print(test_result)

# If that works, run for all species
circ_results <- map_dfr(focal_species, circular_tests, df = overlap.data)


##########
# Detection Density Plots
## -------------------------------------------------------------------
## helper: convert a fitted actmod object to a tidy data frame
## -------------------------------------------------------------------
fit_to_df <- function(fit, treatment) {
  tibble(time_rad = fit@pdf[, 1],          # column 1 = evaluation grid
         density  = fit@pdf[, 2]) %>%      # column 2 = KDE values
    mutate(time_hr = (time_rad %% (2 * pi)) * 24 / (2 * pi),
           treatment = treatment)
}

## -------------------------------------------------------------------
## main plotting function
## -------------------------------------------------------------------
plot_species_density <- function(species_name, df,
                                 reps = 500,          # bootstrap reps for fitact
                                 colours = c(control = "#FABA39FF",
                                             human   = "#7A0403FF")) {

  sp_dat <- df %>% filter(common_name == species_name)

  ## require records in *both* treatments to plot
  if (!all(c("control", "human") %in% sp_dat$treatment)) return(NULL)

  ctrl   <- sp_dat %>% filter(treatment == "control") %>% pull(time_rad)
  human  <- sp_dat %>% filter(treatment == "human")   %>% pull(time_rad)

  if (length(ctrl)  < 10 || length(human) < 10) return(NULL)   # skip thin data

  ## fit the two kernels (same settings for each)
  fc <- fitact(ctrl,  sample = "data", reps = reps)
  fh <- fitact(human, sample = "data", reps = reps)

  ## convert to tidy data frames
  plot_df <- bind_rows(fit_to_df(fc, "control"),
                       fit_to_df(fh, "human"))

  ggplot(plot_df, aes(x = time_hr, y = density,
                      colour = treatment)) +
    geom_line(size = 1.2) +
    scale_colour_manual(values = colours,
                        breaks = c("control", "human"),
                        labels = c("Control", "Human"),
                        name = "Treatment") +
    scale_x_continuous(breaks = c(0, 6, 12, 18, 24),
                       limits = c(0, 24)) +
    labs(title = species_name,
         x = "Time (Hours)",
         y = "Kernel density") +
    theme_minimal(base_size = 12) +
    theme(plot.title      = element_text(face = "bold",
                                         hjust = 0.5, size = 13),
          legend.position = "none")
}

# Generate plots
plot_list <- vector("list", length(focal_species))
names(plot_list) <- focal_species

for (species in focal_species) {
  cat("Creating plot for:", species, "\n")
  plot_list[[species]] <- plot_species_density(species, overlap.data)
}

# Display individual plot (example)
plot_list[[3]]

# Create combined plot grid
valid_plots <- plot_list[!sapply(plot_list, is.null)]
combined_plot_grid <- wrap_plots(valid_plots, ncol = 3)
print(combined_plot_grid)

##########
#Linear Model for Diurnality vs Nocturnality
solar.midnight <- 00:32

overlap <- overlap %>%
  mutate(
    datetime = parse_date_time(start_time, 
                               orders = c("mdy HM", "mdy HMS"),
                               tz = "Africa/Nairobi"),
    
    hour_decimal = hour(datetime) + minute(datetime)/60 + second(datetime)/3600,
    
    # Diurnal activity = absolute distance from solar midnight
    diurnal_metric = abs(hour_decimal - solar.midnight),
    
    # Correct wrap-around for detections near midnight (e.g. 23:00 and 02:00)
    diurnal_metric = if_else(diurnal_metric > 12, 24 - diurnal_metric, diurnal_metric)
  )

# Linear Model for Diurnality vs Nocturnality
solar.midnight <- 00:32  # Note: This should be 0.533 (32 minutes past midnight as decimal)

# Calculate diurnal metric
overlap <- overlap %>%
  mutate(
    datetime = parse_date_time(start_time, 
                               orders = c("mdy HM", "mdy HMS"),
                               tz = "Africa/Nairobi"),
    
    hour_decimal = hour(datetime) + minute(datetime)/60 + second(datetime)/3600,
    
    # Diurnal activity = absolute distance from solar midnight
    diurnal_metric = abs(hour_decimal - 0.533),  # 00:32 = 0.533 hours
    
    # Correct wrap-around for detections near midnight (e.g. 23:00 and 02:00)
    diurnal_metric = if_else(diurnal_metric > 12, 24 - diurnal_metric, diurnal_metric)
  )

# Function to run LMM for a single species
run_species_lmm <- function(species_name, df) {
  cat("Running LMM for:", species_name, "\n")
  
  # Filter data for the specific species
  species_data <- df %>%
    filter(common_name == species_name)
  
  # Check if we have enough data
  if (nrow(species_data) < 10) {
    cat("  Warning: Insufficient data for", species_name, "(n =", nrow(species_data), ")\n")
    return(NULL)
  }
  
  # Fit the linear mixed model
  tryCatch({
    model <- lmer(diurnal_metric ~ treatment + grid_id + treatment * grid_id + (1 | camera_id), 
                  data = species_data)
    
    # Extract model summary
    model_summary <- tidy(model, effects = "fixed")
    model_summary$species <- species_name
    model_summary$n_obs <- nrow(species_data)
    model_summary$n_cameras <- length(unique(species_data$camera_id))
    
    cat("  Successfully fitted model for", species_name, "\n")
    return(list(model = model, summary = model_summary))
    
  }, error = function(e) {
    cat("  Error fitting model for", species_name, ":", e$message, "\n")
    return(NULL)
  })
}

# Run LMMs for all focal species
lmm_results <- vector("list", length(focal_species))
names(lmm_results) <- focal_species

for (species in focal_species) {
  lmm_results[[species]] <- run_species_lmm(species, overlap)
}

# Extract successful model summaries
successful_models <- lmm_results[!sapply(lmm_results, is.null)]
model_summaries <- map_dfr(successful_models, ~ .x$summary)

# Display results
cat("\n=== LMM Results Summary ===\n")
print(model_summaries , n=100)

# Create a more readable summary table
results_table <- model_summaries %>%
  select(species, term, estimate, std.error, statistic, p.value, n_obs, n_cameras) %>%
  mutate(
    p_value_formatted = case_when(
      p.value < 0.001 ~ "< 0.001",
      p.value < 0.01 ~ "< 0.01",
      p.value < 0.05 ~ "< 0.05",
      TRUE ~ as.character(round(p.value, 3))
    ),
    significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  arrange(species, term)

print(results_table, n=100)

# Save results to CSV (optional)
# write.csv(results_table, "species_lmm_results.csv", row.names = FALSE)

##########
#Activity Package exploration
# Subset by species and treatment
ele_control <- subset(overlap, common_name == "African Elephant" & treatment == "control")
ele_human   <- subset(overlap, common_name == "African Elephant" & treatment == "human")

# Fit activity models with bootstrapping
f_control <- fitact(ele_control$time_rad, sample = "data", reps = 1000)
f_human   <- fitact(ele_human$time_rad, sample = "data", reps = 1000)

# Plot
plot(f_control, tline = list(col = "blue"))
plot(f_human, add = TRUE, tline = list(col = "red"))
legend("topright", legend = c("Control", "Human"), col = c("blue", "red"), lty = 1)

compareAct(list(f_control, f_human))
compareCkern(f_control, f_human, reps = 999)

cmean(ele_control$time_rad)
cmean(ele_human$time_rad)

# Plot KDEs for control (blue) and human playback (red)
plot(f_control, yunit = "density", data = "none", tline = list(col = "blue"))
plot(f_human, add = TRUE, data = "none", tline = list(col = "red"))
legend("topright", legend = c("Control", "Human"), col = c("blue", "red"), lty = 1)

# helper to convert radians to 24-h clock (0–24)
rad2hour <- function(rad) (rad %% (2*pi)) / (2*pi) * 24

results <- lapply(focal_species, function(sp) {
  sub_control <- subset(overlap, common_name == sp & treatment == "control")
  sub_human   <- subset(overlap, common_name == sp & treatment == "human")
  
  if (nrow(sub_control) > 10 & nrow(sub_human) > 10) {  # arbitrary threshold
    
    ## 1. Fit activity models
    fc <- fitact(sub_control$time_rad, sample = "data", reps = 500)
    fh <- fitact(sub_human$time_rad,   sample = "data", reps = 500)
    
    ## 2. Compare activity & overlap
    act_comp  <- compareAct(list(fc, fh))
    kern_comp <- compareCkern(fc, fh, reps = 999)
    ovl_val   <- ovl4(fc, fh)
    
    ## 3. Circular means
    mean_ctrl_rad <- cmean(sub_control$time_rad)
    mean_hum_rad  <- cmean(sub_human$time_rad)
    
    data.frame(
      species         = sp,
      # activity comparison
      diff_act        = act_comp[1],              # control – human
      p_act           = act_comp[4],
      # overlap comparison
      Dhat4           = kern_comp["obs"],
      p_overlap       = kern_comp["pNull"],
      # raw overlap index (same as Dhat4 but kept for clarity)
      overlap         = ovl_val,
      # circular means
      mean_ctrl_rad   = mean_ctrl_rad,
      mean_hum_rad    = mean_hum_rad,
      mean_ctrl_hour  = rad2hour(mean_ctrl_rad),
      mean_hum_hour   = rad2hour(mean_hum_rad),
      mean_diff_hour  = rad2hour(mean_ctrl_rad - mean_hum_rad)  # signed shift
    )
  }
})

# combine into one data frame
activity_results <- bind_rows(results)
write.csv(activity_results, "activity.results.csv")

## -------------------------------------------------------------------
## main plotting function
## -------------------------------------------------------------------
plot_species_density <- function(species_name, df,
                                 reps = 500,          # bootstrap reps for fitact
                                 colours = c(control = "#FABA39FF",
                                             human   = "#7A0403FF")) {
  
  sp_dat <- df %>% filter(common_name == species_name)
  
  ## require records in *both* treatments to plot
  if (!all(c("control", "human") %in% sp_dat$treatment)) return(NULL)
  
  ctrl   <- sp_dat %>% filter(treatment == "control") %>% pull(time_rad)
  human  <- sp_dat %>% filter(treatment == "human")   %>% pull(time_rad)
  
  if (length(ctrl)  < 10 || length(human) < 10) return(NULL)   # skip thin data
  
  ## fit the two kernels (same settings for each)
  fc <- fitact(ctrl,  sample = "data", reps = reps)
  fh <- fitact(human, sample = "data", reps = reps)
  
  ## convert to tidy data frames
  plot_df <- bind_rows(fit_to_df(fc, "control"),
                       fit_to_df(fh, "human"))
  
  ggplot(plot_df, aes(x = time_hr, y = density,
                      colour = treatment)) +
    geom_line(size = 1.2) +
    scale_colour_manual(values = colours,
                        breaks = c("control", "human"),
                        labels = c("Control", "Human"),
                        name = "Treatment") +
    scale_x_continuous(breaks = c(0, 6, 12, 18, 24),
                       limits = c(0, 24)) +
    labs(title = species_name,
         x = "Time (Hours)",
         y = "Kernel density") +
    theme_minimal(base_size = 12) +
    theme(plot.title      = element_text(face = "bold",
                                         hjust = 0.5, size = 13),
          legend.position = "none")
}





##########
#making a frequency table of counts by week
# 1. Identify all unique factors
all_combos <- expand_grid(
  common_name = unique(overlap$common_name),
  camera_id   = unique(overlap$camera_id),
  week_id     = unique(overlap$week_id))

#Bout counts
det_counts <- overlap %>%
  group_by(common_name, camera_id, week_id, treatment, grid_id) %>%
  summarise(
    count = n(),   # or n() if no foraging column
    .groups = "drop"
  )
# 2. Merge with actual detections
det_full <- all_combos %>%
  left_join(det_counts, by = c("common_name","camera_id","week_id")) %>%
  mutate(count = tidyr::replace_na(count, 0))

# camera-week lookup (source of truth)
week_lkp <- overlap %>%
  distinct(camera_id, week_id, treatment, grid_id)

# fill NAs in bout_full$treatment / grid_id from the lookup
det_full <- det_full %>%
  left_join(week_lkp, by = c("camera_id","week_id"),
            suffix = c("", ".lkp")) %>%
  mutate(
    treatment = coalesce(treatment, treatment.lkp),
    grid_id   = coalesce(grid_id,   grid_id.lkp)
  ) %>%
  select(-ends_with(".lkp"))

# quick check
colSums(is.na(det_full[c("treatment","grid_id")]))

str(bout_full)

###### now for daylight
# 0) Build det_tgd in a controlled way (types set explicitly)
det_tgd <- overlap %>%
  mutate(
    common_name = as.character(common_name),
    treatment   = as.character(treatment),
    grid_id     = as.character(grid_id),
    # ensure daylight is 0/1 integer
    daylight    = case_when(
      is.logical(daylight) ~ as.integer(daylight),         # TRUE/FALSE -> 1/0
      is.numeric(daylight) ~ as.integer(daylight != 0),    # any nonzero -> 1
      TRUE                 ~ as.integer(trimws(as.character(daylight)) %in% c("1","day","Day","TRUE"))
    )
  ) %>%
  count(common_name, treatment, grid_id, daylight, name = "count")

# 1) Camera-week → (treatment, grid_id) lookup, with types aligned
week_lkp <- overlap %>%
  distinct(camera_id, week_id, treatment, grid_id) %>%
  mutate(
    treatment = as.character(treatment),
    grid_id   = as.character(grid_id)
  )

# 2) Base table of *observed* (treatment, grid) pairs, then cross with species & daylight=0/1
base_tgd <- week_lkp %>%
  distinct(treatment, grid_id) %>%
  mutate(
    treatment = as.character(treatment),
    grid_id   = as.character(grid_id)
  ) %>%
  crossing(
    common_name = unique(overlap$common_name) |> as.character(),
    daylight    = c(0L, 1L)   # <-- integer 0/1 to match det_tgd
  )

# 3) Now the join will match cleanly
det_tgd_full <- base_tgd %>%
  left_join(det_tgd,
            by = c("common_name","treatment","grid_id","daylight")) %>%
  mutate(count = replace_na(count, 0L))

# 4) Quick diagnostics (should be zero rows)
unmatched <- det_tgd_full %>%
  filter(is.na(count)) %>%
  select(common_name, treatment, grid_id, daylight) %>%
  distinct()

nrow(unmatched)  # expect 0

# fit negative binomial GLMM with interaction
m_nightday <- glmer.nb(
  count ~ treatment * daylight + (1 | common_name) + (1 | grid_id),
  data = det_tgd_full
)

summary(m_nightday)



##########
#von mises
library(dplyr)
library(purrr)
library(brms)
library(loo)
library(ggplot2)
library(tidyr)


overlap <- read.csv("solaroverlap.csv")


species_left <- c(
  "African Elephant","Hippopotamus","Lion","Reticulated Giraffe","Leopard",
  "Plains Zebra","Black-backed Jackal","Grevy's Zebra","Zorilla","Common Warthog"
)

# Singles vs mixtures
family_map <- list(
  "African Elephant"     = von_mises(),
  "Hippopotamus"         = von_mises(),
  "Lion"                 = mixture(von_mises(), von_mises(), order = "mu"),
  "Reticulated Giraffe"  = von_mises(),
  "Leopard"              = von_mises(),
  "Plains Zebra"         = mixture(von_mises(), von_mises(), order = "mu"),
  "Black-backed Jackal"  = mixture(von_mises(), von_mises(), order = "mu"),
  "Grevy's Zebra"        = mixture(von_mises(), von_mises(), order = "mu"),
  "Zorilla"              = von_mises(),
  "Common Warthog"       = mixture(von_mises(), von_mises(), order = "mu")
)

wrap_pi <- function(x) ((x + pi) %% (2*pi)) - pi
to_rad  <- function(hr) (hr / 24) * 2*pi
to_hour <- function(rad) ((rad %% (2*pi)) / (2*pi)) * 24
rad2hr  <- function(r) ((r %% (2*pi)) / (2*pi)) * 24

# Empirical circular KDE in hours (boundary corrected)
kde_circ_hours <- function(x_hours, grid = seq(0, 24, length.out = 241), adjust = 1) {
  x <- x_hours %% 24
  x_aug <- c(x - 24, x, x + 24)
  d <- density(x_aug, from = -24, to = 48, adjust = adjust, n = 4096)
  y <- approx(d$x, d$y, xout = grid)$y
  y <- y / sum(y) / mean(diff(grid))
  tibble(hour = grid, density = y)
}

make_dat <- function(sp) {
  overlap %>%
    filter(common_name == sp,
           treatment %in% c("control","human"),
           !is.na(time_rad),
           !is.na(hrs_from_solar_midnight)) %>%
    mutate(
      treatment = factor(treatment, levels = c("control","human")),
      # rotate into solar time
      angle_rot = wrap_pi( to_rad(hrs_from_solar_midnight %% 24) ),
      # keep offset to map predictions back to clock time later if needed
      offset_rad = to_rad((24 - (hrs_from_solar_midnight %% 24)) %% 24)
      # note: angle_rot + offset_rad ≈ time_rad (mod 2π)
    )
}

fit_mix_full <- function(dat) {
  brm(
    bf(
      angle_rot ~ 1,
      mu1    ~ treatment,  mu2    ~ treatment,  mu3    ~ treatment,
      kappa1 ~ treatment,  kappa2 ~ treatment,  kappa3 ~ treatment,
      theta1 ~ treatment,  theta2 ~ treatment
    ),
    data   = dat,
    family = mixture(von_mises(), von_mises(), von_mises(), order = "mu"),
    prior  = c(
      # Intercepts for each distributional parameter
      prior(student_t(3,0,2.5), class="Intercept", dpar="mu1"),
      prior(student_t(3,0,2.5), class="Intercept", dpar="mu2"),
      prior(student_t(3,0,2.5), class="Intercept", dpar="mu3"),
      prior(student_t(3,0,2.5), class="Intercept", dpar="kappa1"),
      prior(student_t(3,0,2.5), class="Intercept", dpar="kappa2"),
      prior(student_t(3,0,2.5), class="Intercept", dpar="kappa3"),
      prior(student_t(3,0,2.5), class="Intercept", dpar="theta1"),
      prior(student_t(3,0,2.5), class="Intercept", dpar="theta2"),
      
      # Treatment slopes
      prior(normal(0,1), class="b", dpar="mu1"),
      prior(normal(0,1), class="b", dpar="mu2"),
      prior(normal(0,1), class="b", dpar="mu3"),
      prior(normal(0,1), class="b", dpar="kappa1"),   # on log(kappa)
      prior(normal(0,1), class="b", dpar="kappa2"),
      prior(normal(0,1), class="b", dpar="kappa3"),
      prior(normal(0,1), class="b", dpar="theta1"),   # logit weights
      prior(normal(0,1), class="b", dpar="theta2")
    ),
    chains = 4, iter = 6000, warmup = 2000, seed = 123,
    control = list(adapt_delta = 0.997, max_treedepth = 13)
  )
}


# Posterior predictive density for ONE treatment, mapped back to CLOCK time
ppd_summaries_clock <- function(fit, dat, trt,
                                grid = seq(0, 24, length.out = 241),
                                draws = 300, nrows = 300, adjust = 1) {
  sub <- dat %>% dplyr::filter(.data$treatment == !!trt)
  if (nrow(sub) == 0) return(NULL)
  sub <- dplyr::slice_sample(sub, n = min(nrows, nrow(sub)))
  
  # simulate in SOLAR space (angle_rot is the response)
  yrep <- posterior_predict(fit, newdata = sub, draws = draws)  # draws x nrows
  
  # map each row’s draws back to CLOCK radians via its offset
  y_clock <- wrap_pi(yrep + matrix(sub$offset_rad, nrow(yrep), ncol(yrep), byrow = TRUE))
  
  # per-draw KDE on clock hours, then summarise to mean & 80% band
  rows <- lapply(seq_len(nrow(y_clock)), function(i) y_clock[i, ])
  dens_list <- lapply(rows, function(rads) kde_circ_hours(rad2hr(rads), grid = grid, adjust = adjust)$density)
  dens_mat <- do.call(cbind, dens_list)
  
  tibble(
    hour = grid,
    mean = rowMeans(dens_mat),
    lo80 = apply(dens_mat, 1, quantile, 0.10),
    hi80 = apply(dens_mat, 1, quantile, 0.90),
    treatment = trt
  )
}

# choose your fitted models list:
species_left <- c(
  "African Elephant","Hippopotamus","Lion","Reticulated Giraffe","Leopard",
  "Plains Zebra","Black-backed Jackal","Grevy's Zebra","Zorilla","Common Warthog"
)

fits_full <- list()

for (sp in species_left) {
  dat <- make_dat(sp)
  have <- dat %>% count(treatment, name="n")
  if (nrow(dat) == 0 || any(have$n < 10)) { message("Skipping ", sp, " (low n)"); next }
  
  message("Fitting full mixture for ", sp, "  (n=", nrow(dat), ")...")
  fit <- fit_mix_full(dat)
  fits_full[[sp]] <- fit
  
  # empirical KDE (clock time) by treatment
  grid <- seq(0, 24, length.out = 241)
  obs <- dat %>% transmute(treatment, hour = rad2hr(wrap_pi(time_rad)))
  emp <- obs %>%
    group_by(treatment) %>%
    summarize(kde = list(kde_circ_hours(hour, grid = grid, adjust = 1)),
              .groups = "drop") %>%
    unnest(kde) %>% rename(empirical = density)
  
  # posterior predictive by treatment (mapped back to CLOCK)
  pp_ctrl <- ppd_summaries_clock(fit, dat, trt = "control", grid = grid)
  pp_hum  <- ppd_summaries_clock(fit, dat, trt = "human",   grid = grid)
  ppd <- bind_rows(pp_ctrl, pp_hum)
  
  dfp <- left_join(emp, ppd, by = c("treatment","hour"))
  
  p <- ggplot(dfp, aes(hour)) +
    geom_ribbon(aes(ymin = lo80, ymax = hi80, fill = treatment), alpha = 0.25, na.rm = TRUE) +
    geom_line(aes(y = mean, color = treatment), linewidth = 0.9, na.rm = TRUE) +
    geom_line(aes(y = empirical, linetype = treatment), linewidth = 0.7, na.rm = TRUE) +
    scale_x_continuous(limits = c(0,24), breaks = seq(0,24,6)) +
    labs(title = sp,
         subtitle = "Posterior predictive (mean ±80%) vs empirical KDE, by treatment",
         x = "Hour of day (clock time)", y = "Density") +
    theme_bw(12) + theme(legend.position = "top")
  
  print(p)
  
  # quick convergence glance in console
  s <- summary(fit)
  cat("\n", sp, ": max Rhat = ", max(s$fixed$Rhat, na.rm = TRUE), 
      " | min bulk_ESS = ", min(s$fixed$ESS_Bulk, na.rm = TRUE), "\n\n", sep = "")
}





walk(names(model_set), function(sp) {
  fit <- model_set[[sp]]
  if (is.null(fit)) return(invisible())
  
  dat <- make_dat(sp)
  if (nrow(dat) == 0) return(invisible())
  
  # --- Empirical KDE on CLOCK time (to compare apples-to-apples) ---
  grid <- seq(0, 24, length.out = 241)
  obs <- dat %>%
    transmute(treatment, hour = rad2hr(wrap_pi(time_rad)))  # clock hours from raw time_rad
  emp <- obs %>%
    group_by(treatment) %>%
    summarize(kde = list(kde_circ_hours(hour, grid = grid, adjust = 1)),
              .groups = "drop") %>%
    unnest(kde) %>% rename(empirical = density)
  
  # --- Posterior predictive by treatment (mapped back to CLOCK time) ---
  ppd_ctrl <- ppd_summaries_clock(fit, dat, trt = "control", grid = grid)
  ppd_hum  <- ppd_summaries_clock(fit, dat, trt = "human",   grid = grid)
  ppd <- bind_rows(ppd_ctrl, ppd_hum)
  
  dfp <- left_join(emp, ppd, by = c("treatment","hour"))
  
  p <- ggplot(dfp, aes(hour)) +
    geom_ribbon(aes(ymin = lo80, ymax = hi80, fill = treatment), alpha = 0.25) +
    geom_line(aes(y = mean, color = treatment), linewidth = 0.9) +
    geom_line(aes(y = empirical, linetype = treatment), linewidth = 0.7) +
    scale_x_continuous(limits = c(0,24), breaks = seq(0,24,6)) +
    labs(title = sp,
         x = "Hour of day (clock time)", y = "Density",
         subtitle = "Posterior predictive (mean ±80%) vs empirical KDE, by treatment") +
    theme_bw(12) + theme(legend.position = "top")
  
  print(p)
})




wrap_pi <- function(x) ((x + pi) %% (2*pi)) - pi
is_mixture <- function(fam) identical(fam$family, "mixture")  # robust check

fit_one <- function(sp) {
  dat <- overlap %>%
    filter(common_name == sp,
           treatment %in% c("control","human"),
           !is.na(time_rad)) %>%
    mutate(angle = wrap_pi(time_rad))
  
  have <- dat %>% count(treatment, name = "n")
  if (nrow(dat) == 0 || n_distinct(dat$treatment) < 2 || any(have$n < 2)) {
    message(sprintf("Skipping '%s': insufficient data (n=%d; control=%d, human=%d)",
                    sp, nrow(dat),
                    have$n[have$treatment=="control"] %||% 0,
                    have$n[have$treatment=="human"] %||% 0))
    return(NULL)
  }
  
  fam <- family_map[[sp]]
  ctrl <- list(adapt_delta = if (is_mixture(fam)) 0.99 else 0.95,
               max_treedepth = 12)
  
  message(sprintf("Fitting %s (%s) n=%d",
                  sp, if (is_mixture(fam)) "mixture" else "single", nrow(dat)))
  
  # Keep it simple: no kappa priors/formulas; rely on defaults
  brm(
    angle ~ treatment,
    data    = dat,
    family  = fam,
    prior   = c(
      prior(normal(0,1), class = "b"),
      prior(student_t(3,0,2.5), class = "Intercept")
    ),
    chains  = 4, iter = 4000, warmup = 1000, seed = 123,
    control = ctrl
  )
}

fits <- map(species_left, ~ tryCatch(fit_one(.x), error = function(e) { 
  message("Error: ", .x, " -> ", e$message); NULL 
})) %>% set_names(species_left)

# Quick summary of treatment effect per species (if present)
diag_table <- map_dfr(names(fits), function(sp) {
  fit <- fits[[sp]]
  if (is.null(fit)) return(tibble(species = sp, Rhat = NA, Bulk_ESS = NA, Tail_ESS = NA))
  
  summ <- as.data.frame(summary(fit)$fixed)
  tibble(
    species = sp,
    Rhat    = mean(summ$Rhat, na.rm = TRUE),
    Bulk_ESS = mean(summ$Bulk_ESS, na.rm = TRUE),
    Tail_ESS = mean(summ$Tail_ESS, na.rm = TRUE)
  )
})
print(diag_table)

walk(names(fits), function(sp) {
  f <- fits[[sp]]
  if (is.null(f)) return(invisible())
  print(pp_check(f) + ggtitle(sp))
})

summ <- purrr::imap_dfr(fits, function(f, sp){
  if (is.null(f)) return(tibble(species = sp, est = NA_real_, l95 = NA_real_, u95 = NA_real_))
  ps <- posterior_summary(f, variable = "b_treatmenthuman")
  tibble(species = sp, est = ps[,"Estimate"], l95 = ps[,"Q2.5"], u95 = ps[,"Q97.5"])
})
summ %>%
  mutate(est_hr  = est  * 24 / (2*pi),
         l95_hr  = l95  * 24 / (2*pi),
         u95_hr  = u95  * 24 / (2*pi))
print(summ, n = Inf)

wrap_pi  <- function(x) ((x + pi) %% (2*pi)) - pi
rad2hr   <- function(r) ((r %% (2*pi)) / (2*pi)) * 24
kde_circ_hours <- function(x_hours, grid = seq(0, 24, length.out = 241), adjust = 1) {
  x <- x_hours %% 24
  x_aug <- c(x - 24, x, x + 24)
  d <- density(x_aug, from = -24, to = 48, adjust = adjust, n = 4096)
  y <- approx(d$x, d$y, xout = grid)$y
  y <- y / sum(y) / mean(diff(grid))
  tibble(hour = grid, density = y)
}
ppd_summaries <- function(fit, treatment, draws = 300, n_per_treat = 300,
                          grid = seq(0, 24, length.out = 241), adjust = 1) {
  nd <- tibble(treatment = factor(rep(treatment, n_per_treat),
                                  levels = c("control","human")))
  yrep <- posterior_predict(fit, newdata = nd, draws = draws)  # draws x n_per_treat (matrix)
  
  # turn rows into a list so do.call(cbind, ...) works
  rows <- asplit(yrep, 1)  # base R: list of numeric vectors (one per draw)
  dens_list <- lapply(rows, function(rads) {
    xh <- rad2hr(wrap_pi(rads))
    kde_circ_hours(xh, grid = grid, adjust = adjust)$density
  })
  
  # bind columns: |grid| x draws
  dens_mat <- do.call(cbind, dens_list)
  
  tibble(
    hour = grid,
    mean = rowMeans(dens_mat),
    lo80 = apply(dens_mat, 1, quantile, 0.10),
    hi80 = apply(dens_mat, 1, quantile, 0.90),
    treatment = treatment
  )
}

# Print posterior predictive vs. empirical KDE inline
walk(names(fits), function(sp) {
  fit <- fits[[sp]]
  if (is.null(fit)) return(invisible())
  
  obs <- overlap %>%
    filter(common_name == sp,
           treatment %in% c("control","human"),
           !is.na(time_rad)) %>%
    transmute(treatment = factor(treatment, levels = c("control","human")),
              hour = rad2hr(wrap_pi(time_rad)))
  
  grid <- seq(0, 24, length.out = 241)
  emp <- obs %>%
    group_by(treatment) %>%
    summarize(kde = list(kde_circ_hours(hour, grid = grid, adjust = 1)), .groups = "drop") %>%
    unnest(kde) %>%
    rename(empirical = density)
  
  ppd_ctrl <- ppd_summaries(fit, "control", grid = grid)
  ppd_hum  <- ppd_summaries(fit, "human", grid = grid)
  ppd <- bind_rows(ppd_ctrl, ppd_hum)
  
  dfp <- left_join(emp, ppd, by = c("treatment","hour"))
  
  p <- ggplot(dfp, aes(hour)) +
    geom_ribbon(aes(ymin = lo80, ymax = hi80, fill = treatment), alpha = 0.25) +
    geom_line(aes(y = mean, color = treatment), linewidth = 0.9) +
    geom_line(aes(y = empirical, linetype = treatment), linewidth = 0.6) +
    scale_x_continuous(limits = c(0,24), breaks = seq(0,24,6)) +
    labs(title = sp, x = "Hour of day", y = "Density",
         subtitle = "Posterior predictive (mean ±80%) vs empirical KDE") +
    theme_bw(12) + theme(legend.position = "top")
  
  print(p)
})



#####################
#Bayesian GAMs
library(tidyverse)
library(purrr)
library(brms)
library(rlang)
library(posterior)
library(tibble)

overlap <- read.csv("solaroverlap.csv")

species_left <- c("African Elephant","Hippopotamus","Lion","Reticulated Giraffe",
                  "Leopard", "Plains Zebra","Black-backed Jackal", "Grevy's Zebra",      
                  "Zorilla", "Common Warthog")

# keep expo in the data
make_bins_1h <- function(dat){
  dat %>%
    dplyr::filter(treatment %in% c("control","human"), !is.na(time_rad)) %>%
    dplyr::transmute(
      treatment = factor(treatment, levels = c("control","human")),
      hour = ((time_rad %% (2*pi)) / (2*pi)) * 24,
      bin  = floor(pmin(23, hour))
    ) %>%
    dplyr::count(treatment, bin, name = "y") %>%
    tidyr::complete(treatment, bin = 0:23, fill = list(y = 0)) %>%
    dplyr::mutate(
      hour_ctr = bin + 0.5,
      expo = 1.0                           # <-- ensure numeric
    )
}

# IMPORTANT: offset inside the formula, not as an argument
fit_brm_spline <- function(bdat, K = 10, family = poisson()){
  stopifnot(is.factor(bdat$treatment))
  
  fml <- as.formula(
    sprintf(
      "y ~ 0 + treatment + s(hour_ctr, by = treatment, bs = 'cc', k = %d) + offset(log(expo))",
      as.integer(K)
    )
  )
  
  brms::brm(
    formula = fml,
    family  = family,
    data    = bdat,
    knots   = list(hour_ctr = c(0, 24)),   # needed for cyclic spline
    prior   = c(
      brms::prior(exponential(1), class = "sds"),
      brms::prior(normal(0, 2), class = "b")
    ),
    chains  = 3, iter = 10000, warmup = 1000, seed = 123,
    control = list(adapt_delta = 0.99)
  )
}

predict_spline_density <- function(fit, grid = seq(0, 24, length.out = 241)){
  new <- expand.grid(
    treatment = factor(c("control","human"), levels = c("control","human")),
    hour_ctr  = grid,
    expo = 1
  )
  mu <- posterior_epred(fit, newdata = new)
  dens <- apply(mu, 1, function(m){
    tibble(treatment = new$treatment, hour = new$hour_ctr, val = m) %>%
      group_by(treatment) %>%
      mutate(val = val / sum(val) / mean(diff(unique(hour)))) %>%
      ungroup()
  })
  bind_rows(dens, .id = "draw") %>%
    group_by(treatment, hour) %>%
    summarise(mean = mean(val),
              lo80 = quantile(val, 0.10),
              hi80 = quantile(val, 0.90),
              .groups = "drop")
}

plot_spline_ppd <- function(dat, fit, sp){
  grid <- seq(0, 24, length.out = 241)
  ppd  <- predict_spline_density(fit, grid = grid)
  
  emp <- dat %>%
    filter(treatment %in% c("control","human"), !is.na(time_rad)) %>%
    transmute(treatment = factor(treatment, levels=c("control","human")),
              hour = ((time_rad %% (2*pi)) / (2*pi))*24) %>%
    group_by(treatment) %>%
    summarise(n_raw = n(),
              kde = list(kde_hours(hour, grid = grid, adjust = 1)),
              .groups = "drop") %>%
    tidyr::unnest(kde) %>%
    rename(empirical = density)
  
  ggplot(left_join(emp, ppd, by = c("treatment","hour")), aes(hour)) +
    geom_ribbon(aes(ymin = lo80, ymax = hi80, fill = treatment),
                alpha = 0.25, na.rm = TRUE) +
    geom_line(aes(y = mean, color = treatment), linewidth = 0.9, na.rm = TRUE) +
    geom_line(aes(y = empirical, linetype = treatment), linewidth = 0.7, na.rm = TRUE) +
    scale_x_continuous(limits = c(0,24), breaks = seq(0,24,6)) +
    labs(x = "Hour of day", y = "Density",
         title = sp, subtitle = "Cyclic spline posterior vs empirical KDE") +
    theme_bw(12) + theme(legend.position = "top")
}


# ----------------------------
# Single-species runner
# ----------------------------
run_species_spline <- function(
    sp,
    overlap,
    return_fit = TRUE,
    save_plot  = FALSE,
    out_dir    = "plots",
    quiet      = FALSE
){
  # minimal, ASCII-only messages
  msg <- function(...) if (!quiet) cat(sprintf(...), "\n")
  
  # 1) subset + bin
  dat  <- overlap[overlap$common_name == sp, , drop = FALSE]
  if (nrow(dat) == 0L) {
    return(list(species = sp, status = "skipped", n_raw = 0L,
                note = "no rows for species", fit = NULL))
  }
  
  bdat <- make_bins_1h(dat)
  
  # require both treatments present after binning
  have <- unique(as.character(bdat$treatment))
  need <- c("control","human")
  if (!all(need %in% have)) {
    return(list(species = sp, status = "skipped",
                n_raw = nrow(dat),
                note = "missing a treatment after binning",
                fit = NULL))
  }
  
  # 2) fit and plot with guardrails
  fit  <- NULL
  plt  <- NULL
  err  <- NULL
  
  tryCatch({
    fit <- fit_brm_spline(bdat)
    ppd <- ppd_from_fit(fit, level = 0.90)
    sig <- ppd_significance(ppd)
    plt <- plot_ppd_no_kde(ppd, sig, species_name = sp)
    print(plt)
    
    if (save_plot) {
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      fn <- file.path(out_dir, paste0(gsub("[^A-Za-z0-9_-]+","_", sp), ".png"))
      ggplot2::ggsave(filename = fn, plot = plt, width = 7, height = 4.5, dpi = 300)
      msg("saved plot: %s", fn)
    }
  }, error = function(e){
    err <<- conditionMessage(e)
  })
  
  if (!is.null(err)) {
    return(list(species = sp, status = "error", n_raw = nrow(dat),
                note = err, fit = if (return_fit) fit else NULL))
  }
  
  list(species = sp, status = "ok", n_raw = nrow(dat),
       note = "", fit = if (return_fit) fit else NULL)
}

# ----------------------------
# Multi-species batch runner
# ----------------------------
run_species_batch <- function(
    species,
    overlap,
    return_fits = TRUE,
    save_plots  = FALSE,
    out_dir     = "plots",
    quiet       = FALSE
){
  results <- vector("list", length(species))
  names(results) <- species
  
  for (i in seq_along(species)) {
    sp <- species[i]
    if (!quiet) cat(sprintf("[%d/%d] %s ...", i, length(species), sp))
    
    # filter + fit
    dat  <- overlap %>% dplyr::filter(common_name == sp)
    bdat <- make_bins_1h(dat)
    if (nrow(bdat) == 0 || !all(c("control","human") %in% bdat$treatment)) {
      results[[i]] <- list(species = sp, status = "skipped", n_raw = nrow(dat), note = "insufficient data")
      if (!quiet) cat("skipped\n")
      next
    }
    
    # fit spline model
    fit <- tryCatch(
      fit_brm_spline(bdat),
      error = function(e) { 
        list(error = conditionMessage(e))
      }
    )
    
    if (inherits(fit, "list") && "error" %in% names(fit)) {
      results[[i]] <- list(species = sp, status = "error", n_raw = nrow(dat), note = fit$error)
      if (!quiet) cat("error\n")
      next
    }
    
    # posterior predictive density + significance
    ppd <- ppd_from_fit(fit, level = 0.90)
    sig <- ppd_significance(ppd)
    
    # plot without KDE
    p <- plot_ppd_no_kde(ppd, sig, species_name = sp)
    
    if (save_plots) {
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
      ggplot2::ggsave(
        file.path(out_dir, paste0(gsub(" ", "_", sp), "_ppd.png")),
        p, width = 7, height = 4, dpi = 300
      )
    }
    
    results[[i]] <- list(
      species = sp,
      status  = "ok",
      n_raw   = nrow(dat),
      note    = "",
      fit     = if (return_fits) fit else NULL
    )
    
    if (!quiet) cat("ok\n")
  }
  
  summary <- data.frame(
    species = vapply(results, `[[`, "", "species"),
    status  = vapply(results, `[[`, "", "status"),
    n_raw   = vapply(results, function(x) as.integer(x$n_raw), 1L),
    note    = vapply(results, `[[`, "", "note"),
    stringsAsFactors = FALSE
  )
  
  list(summary = summary, results = results)
}

 out <- run_species_batch(species_left, overlap, return_fits = TRUE,
                          save_plots = TRUE, out_dir = "ppd_plots", quiet = FALSE)

 ####Model diagnostics
 diag_per_species <- function(out){
   tibble(
     species = names(out$results),
     status  = vapply(out$results, \(x) x$status, ""),
     max_rhat = vapply(
       out$results,
       \(x) if (!is.null(x$fit) && x$status == "ok")
         max(posterior::rhat(x$fit), na.rm = TRUE) else NA_real_,
       numeric(1)
     ),
     min_ess = vapply(
       out$results,
       \(x) if (!is.null(x$fit) && x$status == "ok")
         min(posterior::ess_bulk(x$fit), na.rm = TRUE) else NA_real_,
       numeric(1)
     )
   )
 }
 
 diags <- diag_per_species(out)
 diags

# -----------------------------
# 1) Posterior predictive density
# -----------------------------
ppd_from_fit <- function(fit, level = 0.90,
                         grid = seq(0, 24, length.out = 241)) {
  # newdata for both treatments across the hour grid
  new <- expand.grid(
    treatment = factor(c("control","human"), levels = c("control","human")),
    hour_ctr  = grid,
    expo      = 1
  )
  
  # posterior draws of expected counts
  mu <- brms::posterior_epred(fit, newdata = new)  # draws x nrow(new)
  
  # convert each draw to a density (normalize within treatment)
  dx <- mean(diff(unique(new$hour_ctr)))
  dens_draws <- apply(mu, 1, function(m) {
    tibble::tibble(treatment = new$treatment, hour = new$hour_ctr, val = m) |>
      dplyr::group_by(treatment) |>
      dplyr::mutate(val = val / sum(val) / dx) |>
      dplyr::ungroup()
  })
  
  # summarise across draws
  ppd <- dplyr::bind_rows(dens_draws, .id = "draw") |>
    dplyr::group_by(treatment, hour) |>
    dplyr::summarise(
      mean = mean(val),
      lo   = quantile(val, (1 - level) / 2),
      hi   = quantile(val, 1 - (1 - level) / 2),
      .groups = "drop"
    )
  ppd
}

# -----------------------------
# 2) Credible-interval overlap test (per-hour)
# -----------------------------
ppd_significance <- function(ppd) {
  wide <- tidyr::pivot_wider(
    ppd,
    id_cols = hour,
    names_from  = treatment,
    values_from = c(mean, lo, hi)
  )
  
  wide |>
    dplyr::mutate(
      diff_mean = mean_human - mean_control,
      # "human > control" if the *lower* CI of human exceeds the *upper* CI of control
      sig = dplyr::case_when(
        lo_human > hi_control ~ "higher under human",
        hi_human < lo_control ~ "lower under human",
        TRUE ~ "no diff"
      )
    )
}

# -----------------------------
# 3) Plot without KDE + highlight significant periods
# -----------------------------
plot_ppd_no_kde <- function(ppd, sig_tbl, species_name = "",
                            fill_alpha = 0.25) {
  # collapse contiguous significant runs into rectangles
  runs <- sig_tbl |>
    dplyr::mutate(flg = sig != "no diff",
                  grp = cumsum(dplyr::lag(flg, default = FALSE) != flg)) |>
    dplyr::filter(flg) |>
    dplyr::group_by(sig, grp) |>
    dplyr::summarise(xmin = min(hour), xmax = max(hour), .groups = "drop")
  
  ggplot2::ggplot(ppd, ggplot2::aes(hour)) +
    # light background stripes for significant regions
    ggplot2::geom_rect(
      data = runs,
      ggplot2::aes(xmin = xmin, xmax = xmax,
                   ymin = -Inf, ymax = Inf, fill = sig),
      inherit.aes = FALSE, alpha = 0.07
    ) +
    # ribbons + means for the two treatments
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lo, ymax = hi, fill = treatment),
      alpha = fill_alpha
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = mean, color = treatment),
      linewidth = 0.9
    ) +
    ggplot2::scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, 6)) +
    ggplot2::labs(
      x = "Hour of day", y = "Density",
      title = species_name,
      subtitle = "Cyclic spline posterior (GAM); shaded bands = credible differences"
    ) +
    ggplot2::theme_bw(12) +
    ggplot2::theme(legend.position = "top") +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(alpha = 0.25)))
}

#######
 #support for species becoming more nocturnal or diurnal based on KDE estimates
 
 # ---- trapezoid integration ----
 trapz <- function(x, y) {
   o <- order(x)
   x <- x[o]; y <- y[o]
   sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2, na.rm = TRUE)
 }
 
 # ---- integrate KDE between [a, b] hours ----
 integrate_kde_window <- function(time_hr, density, a, b) {
   o <- order(time_hr)
   x <- time_hr[o]
   y <- density[o]
   if (length(x) < 2) return(NA_real_)
   
   ya <- approx(x, y, xout = a, rule = 2)$y
   yb <- approx(x, y, xout = b, rule = 2)$y
   
   keep <- (x > a & x < b)
   xx <- c(a, x[keep], b)
   yy <- c(ya, y[keep], yb)
   
   trapz(xx, yy)
 }
 
 # ---- 1) define species -> group mapping (edit if your group definitions differ) ----
 group_map <- list(
   `African Elephant` = c("African Elephant"),
   `Hippopotamus`     = c("Hippopotamus"),
   `Large Felids`     = c("Lion", "Leopard"),
   `Zebras`           = c("Plains Zebra", "Grevy's Zebra"),
   `Mesocarnivores`   = c("Black-backed Jackal", "Zorilla")
 )
 
 species_to_group <- enframe(group_map, name = "group", value = "common_name") %>%
   unnest(common_name)
 
 # ---- 2) compute dawn/dusk hours per GROUP from overlap2 ----
 dawn_dusk_by_group <- overlap2 %>%
   inner_join(species_to_group, by = "common_name") %>%
   mutate(
     dawn_dt = ymd_hms(dawn, tz = "UTC"),
     dusk_dt = ymd_hms(dusk, tz = "UTC"),
     dawn_hr = hour(dawn_dt) + minute(dawn_dt)/60 + second(dawn_dt)/3600,
     dusk_hr = hour(dusk_dt) + minute(dusk_dt)/60 + second(dusk_dt)/3600
   ) %>%
   group_by(group) %>%
   summarise(
     dawn_hr = median(dawn_hr, na.rm = TRUE),
     dusk_hr = median(dusk_hr, na.rm = TRUE),
     .groups = "drop"
   )
 
 # ---- 3) integrate KDE daylight density for each GROUP x treatment ----
 kde_daylight_group <- all_kde_group %>%
   left_join(dawn_dusk_by_group, by = "group") %>%
   group_by(group, treatment) %>%
   summarise(
     dawn_hr = mean(dawn_hr),
     dusk_hr = mean(dusk_hr),
     
     total_area    = trapz(time_hr, density),
     daylight_area = integrate_kde_window(time_hr, density, dawn_hr, dusk_hr),
     
     prop_daylight = daylight_area / total_area,
     prop_night    = 1 - prop_daylight,
     .groups = "drop"
   )
 
 kde_daylight_group
write.csv(kde_daylight_group, file = "kde_daylight_group.csv")
