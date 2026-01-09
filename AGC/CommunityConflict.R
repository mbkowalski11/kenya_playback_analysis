######Community Playback Models
library(glmmTMB)
library(DHARMa)
library(broom.mixed)
library(tidyverse)
library(emmeans)
library(scales)
library(patchwork)

#######
#Data Wrangling
farm.effort <- tribble(
  ~Farm, ~night_human, ~night_control,
  "Farm1", 28, 34,
  "Farm2", 28, 24,
  "Farm3", 28, 24,
  "Farm4", 28, 24,
  # Farm5 excluded
  "Farm6", 28, 24,
  "Farm7", 28, 24,
  "Farm8", 28, 25,
  "Farm9", 28, 25,
  "Farm10", 28, 25,
  "Farm11", 28, 25,
  "Farm12", 28, 25,
  "Farm13", 28, 25,
  "Farm14", 28, 25,
  "Farm15", 28, 25,
  "Farm16", 28, 33,
  "Farm17", 28, 33,
  "Farm18", 28, 33,
  "Farm19", 28, 32,
)

write.csv(farm.effort, "farmeffort.csv")

boma.effort <- tribble(
  ~BomaID, ~night_human, ~night_control,
  "Boma1", 30, 29,
  "Boma2", 30, 29,
  "Boma3", 30, 29,
  "Boma4", 30, 29,
  "Boma5", 29, 30,
  "Boma6", 29, 30,
  "Boma7", 34, 26,
  "Boma8", 34, 26,
  "Boma9", 34, 26,
  "Boma10", 34, 26,
  "Boma11", 34, 26,
  "Boma12", 34, 26,
  "Boma13", 40, 22,
  # Boma14 excluded
  "Boma15", 40, 22,
  "Boma16", 40, 22,
  "Boma17", 40, 24,
  "Boma18", 40, 24,
  "Boma19", 40, 24,
  "Boma20", 40, 24
)

write.csv(boma.effort, "bomaeffort.csv")

####Elephant Data
raw_text <- "
Farm	Ele-present	Proximity(m)	Cropdamaged	Cropspeciedamaged	Area(m2)	Treatment
F1	Yes	200	No	NA	NA	Humans
F2	Yes	400	No	NA	NA	Humans
F3	Yes	200	No	NA	NA	Humans
F6	Yes	5	Yes	Cabbage	100	Humans
F7	Yes	100	No	NA	NA	Humans
F8	Yes	5	Yes	Beans	200	Humans
F9	Yes	20	No	NA	NA	Humans
F11	Yes	100	No	NA	NA	Humans
F12	Yes	20	Yes	Maize	150	Humans
F18	Yes	75	Yes	Maize	1500	Humans
F19	Yes	50	Yes	Maize	50	Humans
F1	Yes	100	No	NA	NA	Humans
F2	Yes	300	No	NA	NA	Humans
F3	Yes	100	No	NA	NA	Humans
F4	Yes	200	No	NA	NA	Humans
F5	Yes	200	No	NA	NA	Humans
F6	Yes	200	No	NA	NA	Humans
F7	Yes	100	No	NA	NA	Humans
F8	Yes	5	No	NA	NA	Humans
F9	yes	20	No	NA	NA	Humans
F10	Yes	5	Yes	Snow peas	100	Humans
F15	Yes	15	No	NA	NA	Humans
F3	Yes	100	No	NA	NA	Humans
F4	Yes	20	Yes	Cabbage	50	Humans
F1	Yes	150	No	NA	NA	Humans
F2	Yes	30	No	NA	NA	Humans
F3	Yes	100	No	NA	NA	Humans
F4	Yes	20	Yes	cabbage	15	Humans
F7	Yes	80	Yes	Snow peas, sweet potatoes	150	Humans
F8	Yes	20	No	NA	NA	Humans
F9	Yes	80	No	NA	NA	Humans
F10	Yes	150	Yes	NA	NA	Humans
F16	Yes	50	No	NA	NA	Humans
F1	Yes	20	No	NA	NA	Humans
F2	Yes	20	No	NA	NA	Crickets
F3	Yes	10	Yes	Irish potatoes	8	Crickets
F3	Yes	5	Yes	Irish potatoes	20	Crickets
F4	Yes	4	Yes	Maize and Kales	30	Crickets
F4	Yes	5	Yes	Maize and Kales	12	Crickets
F6	Yes	10	No	NA	NA	Crickets
F6	Yes	20	No	NA	NA	Crickets
F7	Yes	15	Yes	Kales	5	Crickets
F8	Yes	10	Yes	Beans	30	Crickets
F9	Yes	5	Yes	NA	NA	Crickets
F10	Yes	5	No	NA	NA	Crickets
F11	Yes	5	Yes	NA	NA	Crickets
F12	Yes	5	No	NA	NA	Crickets
F16	Yes	10	Yes	cabbage	50	Crickets
F17	Yes	5	Yes	Maize	20	Crickets
F18	Yes	5	Yes	Kales	15	Crickets
F1	Yes	50	No	NA	NA	Crickets
F2	Yes	50	No	NA	NA	Crickets
F3	Yes	10	Yes	Irish potatoes	25	Crickets
F3	Yes	5	Yes	Irish potatoes	8	Crickets
F3	Yes	10	Yes	Irish potatoes	5	Crickets
F4	Yes	5	Yes	Beans	6	Crickets
F4	Yes	15	Yes	Beans	6	Crickets
F4	Yes	10	Yes	Beans	6	Crickets
F6	Yes	3	Yes	Maize	80	Crickets
F6	Yes	5	Yes	Maize	20	Crickets
F7	Yes	20	Yes	Sweet potatoes	4	Crickets
F8	Yes	10	Yes	Snow peas	20	Crickets
F9	Yes	10	Yes	Snow peas	20	Crickets
F10	Yes	10	Yes	Beans	25	Crickets
F10	Yes	5	Yes	Irish potatoes	30	Crickets
F10	Yes	2	Yes	Kales	300	Crickets
F11	Yes	5	Yes	Sweet potatoes	20	Crickets
F11	Yes	50	No	NA	NA	Crickets
F11	Yes	70	No	NA	NA	Crickets
F12	Yes	70	No	NA	NA	Crickets
F12	Yes	50	No	NA	NA	Crickets
F12	Yes	5	Yes	Sweet potatoes	20	Crickets
F13	Yes	2	Yes	Kales	50	Crickets
F13	Yes	7	Yes	Kales	50	Crickets
F14	Yes	200	No	NA	NA	Crickets
F15	Yes	150	No	NA	NA	Crickets
F16	Yes	1	No	NA	NA	Crickets
F18	Yes	100	No	NA	NA	Crickets
F1	Yes	20	Yes	Maize	20	Crickets
F1	Yes	50	No	NA	NA	Crickets
F2	Yes	50	No	NA	NA	Crickets
F3	Yes	10	Yes	Irish potatoes	25	Crickets
F3	Yes	5	Yes	Irish potatoes	8	Crickets
F3	Yes	10	Yes	Irish potatoes	5	Crickets
F4	Yes	5	Yes	Beans	20	Crickets
F4	Yes	5	Yes	Beans	6	Crickets
F4	Yes	5	Yes	Beans	6	Crickets
F4	Yes	5	Yes	Beans	25	Crickets
F6	Yes	20	Yes	Maize	2	Crickets
F6	Yes	5	Yes	Maize	80	Crickets
F6	Yes	20	No	NA	NA	Crickets
F6	Yes	20	No	NA	NA	Crickets
F6	Yes	20	No	NA	NA	Crickets
F7	Yes	15	Yes	Kales	5	Crickets
F8	Yes	1	Yes	Butternut	300	Crickets
F8	Yes	1	Yes	Beans	400	Crickets
F8	Yes	1	Yes	Kales	300	Crickets
F8	Yes	1	Yes	Butternut	300	Crickets
F8	Yes	1	Yes	Butternut	100	Crickets
F9	Yes	5	Yes	Beans	60	Crickets
F9	Yes	1	Yes	Beans	20	Crickets
F9	Yes	10	Yes	Beans	20	Crickets
F9	Yes	20	Yes	Beans	60	Crickets
F9	Yes	20	Yes	Beans	60	Crickets
F10	Yes	1	Yes	Irish potatoes	300	Crickets
F10	Yes	1	Yes	Kales	150	Crickets
F10	Yes	1	Yes	Beans	50	Crickets
F10	Yes	1	Yes	Kales	6	Crickets
F10	Yes	50	No	NA	NA	Crickets
F11	Yes	1	Yes	Irish potatoes	25	Crickets
F11	Yes	5	No	NA	NA	Crickets
F11	Yes	1	Yes	Irish potatoes	25	Crickets
F11	Yes	5	No	NA	NA	Crickets
F12	Yes	1	Yes	Maize	20	Crickets
F12	Yes	1	Yes	Maize	20	Crickets
F12	Yes	1	Yes	Maize	20	Crickets
F12	Yes	1	Yes	Maize	20	Crickets
F13	Yes	1	Yes	Kales	225	Crickets
F13	Yes	1	Yes	Kales	225	Crickets
F14	Yes	5	No	NA	NA	Crickets
F15	Yes	50	No	NA	NA	Crickets
F19	Yes	20	Yes	Sweet potatoes	100	Crickets
F19	Yes	10	Yes	Sweet potatoes	100	Crickets
F19	Yes	20	No	NA	NA	Crickets
F2	Yes	20	Yes	Maize	4	Crickets
F3	Yes	10	Yes	Irish potatoes	30	Crickets
F3	Yes	5	Yes	Irish potatoes	8	Crickets
F3	Yes	10	Yes	Irish potatoes	5	Crickets
F4	Yes	5	Yes	Beans	300	Crickets
F4	Yes	5	Yes	Beans	20	Crickets
F4	Yes	5	Yes	Beans	25	Crickets
F4	Yes	5	Yes	Beans	25	Crickets
F6	Yes	20	Yes	Maize	5	Crickets
F6	Yes	5	Yes	Maize	45	Crickets
F6	Yes	20	Yes	Maize	20	Crickets
F6	Yes	20	Yes	Maize	20	Crickets
F6	Yes	20	No	NA	NA	Crickets
F7	Yes	15	Yes	Kales	4	Crickets
F8	Yes	1	Yes	Butternut	100	Crickets
F8	Yes	1	Yes	Beans	100	Crickets
F8	Yes	1	Yes	Kales	25	Crickets
F8	Yes	1	Yes	Butternut	50	Crickets
F8	Yes	1	Yes	Butternut	25	Crickets
F9	Yes	5	Yes	Beans	100	Crickets
F9	Yes	1	Yes	Beans	100	Crickets
F9	Yes	10	Yes	NA	NA	Crickets
F9	Yes	20	Yes	None	60	Crickets
F9	Yes	20	Yes	None	60	Crickets
F10	Yes	1	Yes	Irish potatoes	300	Crickets
F10	Yes	1	Yes	Kales	150	Crickets
F10	Yes	1	Yes	Beans	50	Crickets
F10	Yes	1	No	NA	NA	Crickets
F10	Yes	50	No	NA	NA	Crickets
F11	Yes	1	Yes	Irish potatoes	100	Crickets
F11	Yes	5	Yes	Irish potatoes	50	Crickets
F11	Yes	1	Yes	Irish potatoes	25	Crickets
F11	Yes	5	No	NA	NA	Crickets
F12	Yes	1	Yes	Maize	30	Crickets
F12	Yes	1	Yes	Maize	40	Crickets
F12	Yes	1	Yes	Maize	50	Crickets
F12	Yes	1	Yes	Maize	20	Crickets
F13	Yes	1	Yes	Kales	100	Crickets
F13	Yes	1	Yes	Kales	150	Crickets
F14	Yes	5	No	NA	NA	Crickets
F15	Yes	50	No	NA	NA	Crickets
F19	Yes	20	Yes	Sweet potatoes	25	Crickets
F19	Yes	10	Yes	Sweet potatoes	25	Crickets
"

# Read as a data frame
elecon <- read.table(
  text       = raw_text,
  header     = TRUE,
  sep        = "\t",
  na.strings = c("NA", ""),
  check.names = TRUE,        
  stringsAsFactors = FALSE
)

# Convert F1 → "Farm 1", etc., and clean types
elecon <- elecon %>%
  mutate(
    Farm = stringr::str_replace(Farm, "^F(\\d+)$", "Farm \\1"),
    Ele.present = factor(`Ele.present`, levels = c("No", "Yes")),
    Treatment   = factor(Treatment),
    Proximity_m = as.numeric(`Proximity.m.`),
    Area_m2     = as.numeric(`Area.m2.`)
  ) %>%
  dplyr::select(Farm, Ele.present, Proximity_m, Cropdamaged,
         Cropspeciedamaged, Area_m2, Treatment)

elecon <- elecon %>%
  mutate(Farm = gsub("^Farm\\s+(\\d+)$", "Farm\\1", Farm))

write.csv(elecon, "elecon.csv")

####Carnivore Data
raw_text <- "
Boma	CarnPresent	CarnSpecies	Proximity	LivestockInjured	LivestockSpeciesInjured	LivestockNumberInjured	LivestockKilled	LivestockSpeciesKilled	LivestockNumberKilled	Treatment
Leruso lopiro	Yes	Spotted hyena	200	No	NA	NA	No	NA	NA	Humans
Leseketeti	Yes	Leopard	100	No	NA	NA	No	NA	NA	Humans
Lenaiptari	Yes	Leopard	1	No	NA	NA	Yes	Sheep	1	Humans
Lengolos	Yes	Spotted hyena	400	No	NA	NA	No	NA	NA	Humans
Lekisima	Yes	Spotted hyena	100	No	NA	NA	No	NA	NA	Humans
Lenaino2	Yes	Spotted hyena	100	No	NA	NA	No	NA	NA	Humans
Nakoyia Lesamaita	Yes	Spotted hyena	NA	No	NA	NA	No	NA	NA	Humans
Kokoroi Telaita	Yes	Spotted hyena	100	No	NA	NA	No	NA	NA	Humans
Long'oria	Yes	Spotted hyena	200	No	NA	NA	No	NA	NA	Humans
Ang'orinyang	Yes	Spotted hyena	300	No	NA	NA	No	NA	NA	Humans
Susan	Yes	Spotted hyena	200	No	NA	NA	No	NA	NA	Humans
Leruso lopiro	Yes	Spotted hyena	200	No	NA	NA	No	NA	NA	Humans
Lesangurikuri	Yes	Spotted hyena	50	No	NA	NA	No	NA	NA	Humans
Lekulal	Yes	Spotted hyena	100	No	NA	NA	No	NA	NA	Humans
Lengolos	Yes	Spotted hyena	400	No	NA	NA	No	NA	NA	Humans
Lekisima	Yes	Spotted hyena	100	No	NA	NA	No	NA	NA	Humans
Lenaino2	Yes	Spotted hyena	100	No	NA	NA	No	NA	NA	Humans
Nakoyia Lesamaita	Yes	Spotted hyena	NA	No	NA	NA	No	NA	NA	Humans
Lesangurikuri2	Yes	Spotted hyena	5	No	NA	NA	No	NA	NA	Humans
Kokoroi Telaita	Yes	Spotted hyena	100	No	NA	NA	No	NA	NA	Humans
Long'oria	Yes	Spotted hyena	150	No	NA	NA	No	NA	NA	Humans
Ang'orinyang	Yes	Spotted hyena	200	No	NA	NA	No	NA	NA	Humans
Susan	Yes	Spotted hyena	200	No	NA	NA	No	NA	NA	Humans
Leruso lopiro	Yes	Spotted hyena	100	No	NA	NA	No	NA	NA	Humans
Lewarani	Yes	Spotted hyena	80	No	NA	NA	No	NA	NA	Humans
Lenaino2	Yes	Spotted hyena	50	No	NA	NA	No	NA	NA	Humans
Nakoyia Lesamaita	Yes	Spotted hyena	50	No	NA	NA	No	NA	NA	Humans
Long'oria	Yes	Spotted hyena	NA	No	NA	NA	No	NA	NA	Humans
Susan	Yes	Leopard	200	No	NA	NA	No	NA	NA	Humans
Leseketeti	Yes	Leopard	10	No	NA	NA	No	NA	NA	Humans
Leseketeti	Yes	Leopard	10	No	NA	NA	No	NA	NA	Humans
Lenaiptari	Yes	Leopard	10	No	NA	NA	No	NA	NA	Humans
Lenaiptari	Yes	Leopard	10	No	NA	NA	No	NA	NA	Humans
Lengolos	Yes	Spotted hyenae	100	No	NA	NA	No	NA	NA	Humans
Lewarani	Yes	Spotted hyenae	100	No	NA	NA	No	NA	NA	Humans
Lekisima	Yes	Spotted hyenae	50	No	NA	NA	No	NA	NA	Humans
Lenaino2	Yes	Spotted hyenae	100	No	NA	NA	No	NA	NA	Humans
Long'oria	Yes	Spotted hyenae	300	No	NA	NA	No	NA	NA	Humans
Ang'orinyang	Yes	Spotted hyenae	400	No	NA	NA	No	NA	NA	Humans
Susan	Yes	Spotted hyenae	200	No	NA	NA	No	NA	NA	Humans
Leseketeti	Yes	Leopard	1	No	NA	NA	Yes	Sheep	1	Crickets
Lenaiptari	Yes	Leopard	1	No	NA	NA	Yes	Sheep	1	Crickets
Lesangurikuri	Yes	Spotted hyenae	100	No	NA	NA	No	NA	NA	Humans
Lesangurikuri	Yes	Spotted hyenae	100	No	NA	NA	No	NA	NA	Humans
Leseketeti	Yes	Leopard	50	No	NA	NA	No	NA	NA	Crickets
Leseketeti	Yes	Leopard	20	No	NA	NA	No	NA	NA	Crickets
Lenaiptari	Yes	Leopard	50	No	NA	NA	No	NA	NA	Crickets
Long'oria	Yes	Spotted hyenae	300	No	NA	NA	No	NA	NA	Humans
Susan	Yes	Leopard	300	No	NA	NA	No	NA	NA	Humans
Leruso lopiro	Yes	Spotted hyenae	10	No	NA	NA	No	NA	NA	Crickets
Lesangurikuri	Yes	Spotted hyenae	100	No	NA	NA	No	NA	NA	Crickets
Lekulal	Yes	Spotted hyenae	1	No	NA	NA	No	NA	NA	Crickets
Lenaiptari	Yes	Leopard	10	No	NA	NA	No	NA	NA	Crickets
Lenaiptari	Yes	Leopard	15	No	NA	NA	No	NA	NA	Crickets
Lewarani	Yes	Spotted hyenae	10	No	NA	NA	No	NA	NA	Crickets
Lenaino2	Yes	Spotted hyenae	20	No	NA	NA	No	NA	NA	Crickets
Lenaino2	Yes	Spotted hyenae	20	No	NA	NA	No	NA	NA	Crickets
Nakoyia Lesamaita	Yes	Spotted hyenae	20	No	NA	NA	No	NA	NA	Crickets
Nakoyia Lesamaita	Yes	Spotted hyenae	15	No	NA	NA	No	NA	NA	Crickets
Lesangurikuri2	Yes	Spotted hyenae	40	No	NA	NA	No	NA	NA	Crickets
Long'oria	Yes	Leopard	200	No	NA	NA	No	NA	NA	Humans
Susan	Yes	Leopard	400	No	NA	NA	No	NA	NA	Humans
Leruso lopiro	Yes	Spotted hyenae	15	No	NA	NA	No	NA	NA	Crickets
Lesangurikuri	Yes	Spotted hyenae	10	No	NA	NA	No	NA	NA	Crickets
Lekulal	Yes	Spotted hyenae	5	No	NA	NA	No	NA	NA	Crickets
Leseketeti	Yes	Leopard	1	No	NA	NA	Yes	NA	NA	Crickets
Leseketeti	Yes	Leopard	10	No	NA	NA	No	NA	NA	Crickets
Lenaiptari	Yes	Leopard	20	No	NA	NA	No	NA	NA	Crickets
Lenaiptari	Yes	Leopard	15	No	NA	NA	No	NA	NA	Crickets
Lewarani	Yes	Spotted hyenae	50	No	NA	NA	No	NA	NA	Crickets
Lekisima	Yes	Spotted hyenae	5	No	NA	NA	No	NA	NA	Crickets
Lenaino2	Yes	Spotted hyenae	70	No	NA	NA	No	NA	NA	Crickets
Nakoyia Lesamaita	Yes	Spotted hyenae	50	No	NA	NA	No	NA	NA	Crickets
Nakoyia Lesamaita	Yes	Spotted hyenae	20	No	NA	NA	No	NA	NA	Crickets
Lesangurikuri2	Yes	Spotted hyenae	30	No	NA	NA	No	NA	NA	Crickets
Leruso lopiro	Yes	Spotted hyenae	20	No	NA	NA	No	NA	NA	Crickets
Lesangurikuri	Yes	Spotted hyenae	10	No	NA	NA	No	NA	NA	Crickets
Lekulal	Yes	Spotted hyenae	20	No	NA	NA	No	NA	NA	Crickets
Lenaino	Yes	Spotted hyenae	30	No	NA	NA	No	NA	NA	Crickets
Leseketeti	Yes	Leopard	5	No	NA	NA	No	NA	NA	Crickets
Leseketeti	Yes	Leopard	20	No	NA	NA	No	NA	NA	Crickets
Lenaiptari	Yes	Leopard	25	No	NA	NA	No	NA	NA	Crickets
Lenaiptari	Yes	Leopard	50	No	NA	NA	No	NA	NA	Crickets
Lengolos	Yes	Spotted hyenae	100	No	NA	NA	No	NA	NA	Crickets
Lewarani	Yes	Spotted hyenae	50	No	NA	NA	No	NA	NA	Crickets
Lekisima	Yes	Spotted hyenae	5	No	NA	NA	No	NA	NA	Crickets
Lenaino2	Yes	Spotted hyenae	80	No	NA	NA	No	NA	NA	Crickets
Nakoyia Lesamaita	Yes	Spotted hyenae	15	No	NA	NA	No	NA	NA	Crickets
Long'oria	Yes	Spotted hyenae	100	No	NA	NA	No	NA	NA	Crickets
Ang'orinyang	Yes	Spotted hyenae	1	Yes	Sheep	1	No	NA	NA	Crickets
Susan	Yes	Spotted hyenae	200	No	NA	NA	No	NA	NA	Crickets
Leruso lopiro	Yes	Spotted hyenae	10	No	NA	NA	No	NA	NA	Crickets
Leruso lopiro	Yes	Spotted hyenae	20	No	NA	NA	No	NA	NA	Crickets
Leruso lopiro	Yes	Spotted hyenae	5	No	NA	NA	No	NA	NA	Crickets
Lesangurikuri	Yes	Spotted hyenae	5	No	NA	NA	No	NA	NA	Crickets
Lesangurikuri	Yes	Spotted hyenae	5	No	NA	NA	No	NA	NA	Crickets
Lekulal	Yes	Spotted hyenae	10	No	NA	NA	No	NA	NA	Crickets
Lenaino	Yes	Spotted hyenae	5	No	NA	NA	No	NA	NA	Crickets
Leseketeti	Yes	Leopard	1	Yes	Dog	1	No	NA	NA	Crickets
Leseketeti	Yes	Leopard	5	No	NA	NA	No	NA	NA	Crickets
Lenaiptari	Yes	Leopard	1	Yes	Sheep	1	Yes	Sheep	2	Crickets
Lenaiptari	Yes	Leopard	10	No	NA	NA	No	NA	NA	Crickets
Lengolos	Yes	Spotted hyenae	5	No	NA	NA	No	NA	NA	Crickets
Lengolos	Yes	Spotted hyenae	1	No	NA	NA	No	NA	NA	Crickets
Lengolos	Yes	Spotted hyenae	1	No	NA	NA	No	NA	NA	Crickets
Lengolos	Yes	Spotted hyenae	2	No	NA	NA	No	NA	NA	Crickets
Lekisima	Yes	Spotted hyenae	2	No	NA	NA	No	NA	NA	Crickets
Lenaino2	Yes	Spotted hyenae	2	No	NA	NA	No	NA	NA	Crickets
Nakoyia Lesamaita	Yes	Spotted hyenae	15	No	NA	NA	No	NA	NA	Crickets
Lesangurikuri2	Yes	Spotted hyenae	5	No	NA	NA	No	NA	NA	Crickets
Lesangurikuri2	Yes	Spotted hyenae	10	No	NA	NA	No	NA	NA	Crickets
Lesangurikuri2	Yes	Spotted hyenae	1	No	NA	NA	No	NA	NA	Crickets
Long'oria	Yes	Spotted hyenae	20	No	NA	NA	No	NA	NA	Crickets
Ang'orinyang	Yes	Leopard	10	No	NA	NA	No	NA	NA	Crickets
Susan	Yes	Leopard	1	Yes	Goat	1	No	NA	NA	Crickets
"

carncon <- read.table(
  text       = raw_text,
  header     = TRUE,
  sep        = "\t",
  quote      = "\"",       
  na.strings = c("NA", ""),
  stringsAsFactors = FALSE
)

boma_lookup <- tribble(
  ~Boma, ~BomaID,
  "Leruso lopiro",      "Boma1",
  "Lesangurikuri",    "Boma2",
  "Lekulal",            "Boma3",
  "Lenaino",          "Boma4",
  "Leseketeti",         "Boma5",
  "Lenaiptari",         "Boma6",
  "Lengolos",           "Boma7",
  "Lewarani",           "Boma8",
  "Lekisima",           "Boma9",
  "Lenaino2",          "Boma10",
  "Nakoyia Lesamaita",  "Boma11",
  "Lesangurikuri2",    "Boma12",
  "Kokoroi Telaita",    "Boma13",
  "Kasirtei Lotulia",   "Boma14",
  "Limasiya Akale",     "Boma15",
  "Kachichi Akale",     "Boma16",
  "Chepkomol",          "Boma17",
  "Long'oria",          "Boma18",
  "Ang'orinyang",       "Boma19",
  "Susan",              "Boma20"
)

carncon <- carncon %>%
  left_join(boma_lookup, by = "Boma") %>%
  relocate(BomaID, .before = CarnPresent)

carncon <- carncon %>%
  mutate(CarnAttack = rowSums(across(c(LivestockNumberInjured, LivestockNumberKilled)), na.rm = TRUE)) %>%
  relocate(CarnAttack, .before = Treatment)

write.csv(carncon, "carncon.csv")

# Diagnostic Functions
# Overdispersion check function
check_overdispersion <- function(model) {
  rdf <- df.residual(model)
  rp <- resid(model, type = "pearson")
  ratio <- sum(rp^2) / rdf
  pval <- pchisq(sum(rp^2), df = rdf, lower.tail = FALSE)
  data.frame(
    chisq = sum(rp^2),
    ratio = ratio,
    rdf = rdf,
    p_value = pval
  )
}


##########
#Elephant Models
elencon <- read.csv("elecon.csv")
farm.effort <- read.csv("farmeffort.csv")
#Detection Rate - Farm Level
elecon_detect <- elecon %>%
  group_by(Farm, Treatment) %>%
  summarise(n_conflicts = n()) %>%
  left_join(farm.effort, by = "Farm") %>%
  mutate(nights = ifelse(Treatment == "Humans", night_human, night_control))

ele_rate <- glmmTMB(
  n_conflicts ~ Treatment + offset(log(nights)) + (1 | Farm),
  family = poisson,
  data = elecon_detect
)
check_overdispersion(ele_rate) #below 1 so can use Poisson
#chisq     ratio  rdf   p_value
#22.118 0.7372667  30 0.8497172
summary(ele_rate)
#Conditional model:
#                Estimate Std. Error z value Pr(>|z|)    
#(Intercept)      -1.4108     0.1712  -8.240  < 2e-16 ***
#TreatmentHumans  -1.3696     0.1966  -6.965 3.28e-12 ***
exp(-1.3696) #.254 or a 74.6% reduction in conflict due to human treatment

coef <- fixef(ele_rate)$cond
rate_control <- exp(coef["(Intercept)"])
rate_human <- exp(coef["(Intercept)"] + coef["TreatmentHumans"])
convert <- function(x) {
  tibble(per_night = x, per_month = x * 30, per_year  = x * 365)
  }
eledetrates_table <- bind_rows(
  Control = convert(rate_control),Humans  = convert(rate_human),.id = "Treatment")
eledetrates_table
#Treatment per_night per_month per_year
#Control      0.244       7.32     89.0
#Humans       0.0620      1.86     22.6

#Proximity
eleprox <- glmmTMB(
  Proximity_m ~ Treatment + (1|Farm),
  family = Gamma(link = "log"),
  data = elecon)
summary(eleprox)
#Conditional model:
#                 Estimate Std. Error z value Pr(>|z|)    
#(Intercept)       2.6462     0.2177  12.154   <2e-16 ***
#TreatmentHumans   1.8822     0.2115   8.898   <2e-16 ***
exp(2.646)       # Control mean proximity - 14.1m
exp(1.882)       # Multiplicative effect of Humans - 92.7m (6.57x greater distance than controls)


#Elephant Conversion Rate
ele_visit <- elecon %>%
  filter(Ele.present == "Yes") %>%
  mutate(conflict = as.integer(Cropdamaged == "Yes"))

ele_con_conv <- glmmTMB(
  conflict ~ Treatment + (1 | Farm),
  family = binomial(link = "logit"),
  data   = ele_visit
)

summary(ele_conf_conv)
#Conditional model:
#Estimate Std. Error z value Pr(>|z|)    
#(Intercept)       1.1201     0.4192   2.672  0.00754 ** 
#TreatmentHumans  -2.3837     0.5408  -4.407 1.05e-05 ***
emmeans(ele_conf_conv, ~ Treatment, type = "response")
#Treatment  prob     SE  df asymp.LCL asymp.UCL
#Crickets  0.754 0.0778 Inf    0.5741     0.875
#Humans    0.220 0.0990 Inf    0.0837     0.466
#So 75.4% probability of conflict vs 22% probability of conflict, given presence


#Elephant Damage Rate
elecon_damage <- elecon %>%
  group_by(Farm, Treatment) %>%
  summarise(sum_damage = sum(Area_m2, na.rm = TRUE)) %>%
  left_join(farm.effort, by = "Farm") %>%
  mutate(nights = ifelse(Treatment == "Humans", night_human, night_control))

ele_damrate <- glmmTMB(
  sum_damage ~ Treatment + offset(log(nights)) + (1 | Farm),
  family = nbinom2(link = "log"),
  data = elecon_damage2
)
check_overdispersion(ele_damrate) #0 so need to use NB, which corrects the overdispersion
#chisq     ratio    rdf   p_value ######Poison
#5950.822 198.3607  30       0
#chisq     ratio     rdf   p_value ######NB
#26.75128 0.9224578  29   0.5851224
summary(ele_damrate)
#Conditional model:
#                Estimate Std. Error z value Pr(>|z|)    
#(Intercept)       2.6226     0.4992   5.254 1.49e-07 ***
#TreatmentHumans  -0.9156     0.7405  -1.236    0.216   
exp(-0.9156) #.400 or a 60% reduction in damage due to human treatment
###Raw numbers: 6239m2 for crickets, 2315m2 for humans

#Conditional model (outlier excluded):
#                Estimate Std. Error z value Pr(>|z|)    
#(Intercept)       2.6225     0.4641   5.651  1.6e-08 ***
#TreatmentHumans  -1.7907     0.6890  -2.599  0.00935 **  
exp(-1.7907) #.167 or a 83% reduction in damage due to human treatment

coef <- fixef(ele_damrate)$cond
rate_control <- exp(coef["(Intercept)"])
rate_human <- exp(coef["(Intercept)"] + coef["TreatmentHumans"])
convert <- function(x) {
  tibble(per_night = x, per_month = x * 30, per_year  = x * 365)
}
eledamrates_table <- bind_rows(
  Control = convert(rate_control),Humans  = convert(rate_human),.id = "Treatment")
eledamrates_table
#Treatment per_night per_month per_year
#Control       13.8       413.    5026.
#Humans         5.51      165.    2012.  (outlier removed: 2.30      68.9     839)
##These farms are about 1,000-2,500m2 



##########
#Carnivore Models
carncon <- read.csv("carncon.csv")
boma.effort <- read.csv("bomaeffort.csv")
#Detection Rate - Boma Level
carncon_detect <- carncon %>%
  group_by(BomaID, Treatment) %>%
  summarise(n_detections = n()) %>%
  left_join(boma.effort, by = "BomaID") %>%
  mutate(nights = ifelse(Treatment == "Humans", night_human, night_control))

carn_detrate <- glmmTMB(
  n_detections ~ Treatment + offset(log(nights)) + (1 | BomaID),
  family = poisson,
  data = carncon_detect
)
check_overdispersion(carn_detrate) #below 1 so can use Poisson
#chisq     ratio  rdf   p_value
#21.98068 0.8140993  27 0.7383583
summary(carn_detrate)
#Conditional model:
#                Estimate Std. Error z value Pr(>|z|)    
#(Intercept)      -1.7673     0.1204 -14.680  < 2e-16 ***
#TreatmentHumans  -0.6424     0.1903  -3.375 0.000739 ***
exp(-.6424) #.526 or a 48.4% reduction in detection due to human treatment

coef <- fixef(carn_detrate)$cond
rate_control <- exp(coef["(Intercept)"])
rate_human <- exp(coef["(Intercept)"] + coef["TreatmentHumans"])
convert <- function(x) {
  tibble(per_night = x, per_month = x * 30, per_year  = x * 365)
}
cardetrates_table <- bind_rows(
  Control = convert(rate_control),Humans  = convert(rate_human),.id = "Treatment")
cardetrates_table
#Treatment per_night per_month per_year
#Control      0.171       5.12     62.3
#Humans       0.0898      2.70     32.8


#Proximity
carnprox <- glmmTMB(
  Proximity ~ Treatment + (1|BomaID),
  family = Gamma(link = "log"),
  data = carncon)
summary(carnprox)
#Conditional model:
#                 Estimate Std. Error z value Pr(>|z|)    
#(Intercept)       3.1369     0.1630   19.25  < 2e-16 ***
#TreatmentHumans   1.7205     0.2288    7.52 5.49e-14 ***
exp(3.1369)       # Control proximity (meters) - 23.03m
exp(1.7205)       # Multiplicative effect of Humans - 128.7m (5.59x greater distance than controls)

#Carnivore Conversion Rate
carn_visit <- carncon %>%
  filter(CarnPresent == "Yes") %>%
  mutate(conflict = as.integer(
    LivestockInjured == "Yes" | LivestockKilled == "Yes"
  ))
carn_conf_conv <- glmmTMB(
  conflict ~ Treatment + (1|BomaID),
  family = binomial,
  data = carn_visit
)
summary(carn_conf_conv)
#Conditional model:
#                 Estimate Std. Error z value Pr(>|z|)   
#(Intercept)       -2.935      1.023  -2.868  0.00413 **
#TreatmentHumans   -1.635      1.137  -1.437  0.15068   
emmeans(carn_conf_conv, ~ Treatment, type = "response")
#Treatment   prob     SE  df asymp.LCL asymp.UCL
#Crickets  0.0505 0.0490 Inf  0.007103     0.283
#Humans    0.0103 0.0146 Inf  0.000617     0.148
#So 5.05% change of attacks vs 1.03% change of attacks when carnivores present


#Conflict Damage Rate
carncon_damage <- carncon %>%
  group_by(BomaID, Treatment) %>%
  summarise(n_conflicts = sum(CarnAttack, na.rm = TRUE)) %>%
  left_join(boma.effort, by = "BomaID") %>%
  mutate(  nights = ifelse(Treatment == "Humans", night_human, night_control))

carn_damrate <- glmmTMB(
  n_conflicts ~ Treatment + offset(log(nights)) + (1 | BomaID),
  family = poisson,
  data = carncon_damage
)
check_overdispersion(carn_damrate) #below 1 so can use Poisson
#chisq     ratio  rdf   p_value
#4.26984 0.1581422  27 0.9999998
summary(carn_damrate)
#Conditional model:
#                Estimate Std. Error z value Pr(>|z|)    
#(Intercept)       -5.325      1.196  -4.453 8.47e-06 ***
#TreatmentHumans   -2.229      1.067  -2.089   0.0367 *  
exp(-2.229) #.107 or a 89.3% reduction in damage due to human treatment

coef <- fixef(carn_damrate)$cond
rate_control <- exp(coef["(Intercept)"])
rate_human <- exp(coef["(Intercept)"] + coef["TreatmentHumans"])
convert <- function(x) {
  tibble(per_night = x, per_month = x * 30, per_year  = x * 365)
}
cardamrates_table <- bind_rows(
  Control = convert(rate_control),Humans  = convert(rate_human),.id = "Treatment")
cardamrates_table
#Treatment per_night per_month per_year
#Control      0.00487    0.146     1.78
#Humans       0.000524   0.0157    0.191

#########
# PLOTS
#########

## Shared colours, labels, and theme ----
treat_cols <- c("Crickets" = "#FABA39FF",
                "Humans"   = "#7A0403FF")

treat_labels <- c("Crickets" = "C",
                  "Humans"   = "H")

base_theme <- theme_bw() +
  theme(text = element_text(family = "Palatino", size = 180),
        axis.text = element_text(size = 180, color = "black"),
        axis.text.y = element_text(size = 180, color = "black"),
        axis.text.x = element_text(size = 180, color = "black", margin = margin(t = 20)),
        axis.title = element_text(size = 180),
        panel.grid = element_line(color = "grey90"),
        panel.border = element_rect(linewidth = 5, color = "black"),
        legend.position = "none"
  )


## Convenience jitter layer ----
jitter_layer <- function(raw_df, yvar) {
  geom_jitter(
    data  = raw_df,
    aes(x = Treatment, y = {{ yvar }}, fill = Treatment),
    width  = 0.15,
    height = 0,  
    shape  = 21,
    size   = 20,
    alpha  = 0.35,
    colour = "black",
    stroke = 1,
    show.legend = FALSE
  )
}

############################
## ELEPHANT RAW DATA
############################

# Detection rate: conflicts per night per farm
ele_detect_raw <- elecon_detect %>%
  mutate(rate_night = n_conflicts / nights) %>%
  filter(is.finite(rate_night))

# Proximity: individual measurements (CORRECT AS IS)
ele_prox_raw <- elecon %>%
  filter(!is.na(Proximity_m))

# Conversion rate: proportion of visits resulting in damage per farm
ele_conv_raw <- ele_visit %>%
  group_by(Farm, Treatment) %>%
  summarise(
    n_visits = n(),
    n_damage = sum(conflict),
    prop_damage = n_damage / n_visits,
    .groups = "drop"
  )

# Damage rate: m² per night per farm
ele_dam_raw <- elecon_damage2 %>%
  mutate(rate_night = sum_damage / nights) %>%
  filter(is.finite(rate_night))

############################
## CARNIVORE RAW DATA 
############################

# Detection rate: detections per night per boma
carn_detect_raw <- carncon_detect %>%
  mutate(rate_night = n_detections / nights) %>%
  filter(is.finite(rate_night))

# Proximity: individual measurements (CORRECT AS IS)
carn_prox_raw <- carncon %>%
  filter(!is.na(Proximity))

# Conversion rate: proportion of visits resulting in attacks per boma
carn_conv_raw <- carn_visit %>%
  group_by(BomaID, Treatment) %>%
  summarise(
    n_visits = n(),
    n_attacks = sum(conflict),
    prop_attack = n_attacks / n_visits,
    .groups = "drop"
  )

# Damage rate: attacks per night per boma
carn_dam_raw <- carncon_damage %>%
  mutate(rate_night = n_conflicts / nights) %>%
  filter(is.finite(rate_night))

############################
## ELEPHANT PLOTS
############################

## 1. Elephant detection rate (per night) ----
emm_ele_rate <- emmeans(
  ele_rate, ~ Treatment,
  type   = "response",
  offset = 0
)

ele_rate_df <- as.data.frame(emm_ele_rate)
# Standardize column names
colnames(ele_rate_df)[colnames(ele_rate_df) %in% c("rate", "response")] <- "response"
ele_rate_df <- ele_rate_df %>%
  rename(LCL = asymp.LCL,
         UCL = asymp.UCL)

p_ele_detect <- ggplot(ele_rate_df,
                       aes(x = Treatment, y = response)) +
  jitter_layer(ele_detect_raw, rate_night) +
  geom_point(aes(fill = Treatment), size = 24, shape = 21, colour = "black", stroke = 1, show.legend = FALSE) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL, colour = Treatment),
                width = 0.05, linewidth = 8, show.legend = FALSE) +
  scale_fill_manual(values = treat_cols) +
  scale_colour_manual(values = treat_cols) +
  scale_x_discrete(labels = treat_labels) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6), limits = c(0, .7)) +
  labs(x = NULL,
       y = "Detections/Night") +
  base_theme

## 2. Elephant proximity (m) ----
emm_ele_prox <- emmeans(eleprox, ~ Treatment, type = "link")
ele_prox_df  <- as.data.frame(emm_ele_prox) %>%
  mutate(
    response = exp(emmean),
    LCL      = exp(asymp.LCL),
    UCL      = exp(asymp.UCL)
  )

p_ele_prox <- ggplot(ele_prox_df,
                     aes(x = Treatment, y = response)) +
  jitter_layer(ele_prox_raw, Proximity_m) +
  geom_point(aes(fill = Treatment), size = 24, shape = 21, colour = "black", stroke = 1, show.legend = FALSE) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL, colour = Treatment),
                width = 0.05, linewidth = 8, show.legend = FALSE) +
  scale_fill_manual(values = treat_cols) +
  scale_colour_manual(values = treat_cols) +
  scale_x_discrete(labels = treat_labels) +
  labs(x = NULL,
       y = "Mean Proximity (m)") +
  base_theme

## 3. Elephant conversion (P(crop damage | present)) ----
emm_ele_conv <- emmeans(ele_conf_conv, ~ Treatment, type = "response")
ele_conv_df  <- as.data.frame(emm_ele_conv) %>%
  rename(
    response = prob,
    LCL      = asymp.LCL,
    UCL      = asymp.UCL
  )

p_ele_conv <- ggplot(ele_conv_df,
                     aes(x = Treatment, y = response)) +
  jitter_layer(ele_conv_raw, prop_damage) +
  geom_point(aes(fill = Treatment), size = 24, shape = 21, colour = "black", stroke = 1, show.legend = FALSE) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL, colour = Treatment),
                width = 0.05, linewidth = 8, show.legend = FALSE) +
  scale_fill_manual(values = treat_cols) +
  scale_colour_manual(values = treat_cols) +
  scale_x_discrete(labels = treat_labels) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(x = NULL,
       y = "P(Damage | Presence)") +
  base_theme

## 4. Elephant crop damage rate (m² per night) ----
emm_ele_dam <- emmeans(
  ele_damrate, ~ Treatment,
  type   = "response",
  offset = 0
)

ele_dam_df <- as.data.frame(emm_ele_dam)
# Standardize column names
colnames(ele_dam_df)[colnames(ele_dam_df) %in% c("rate", "response")] <- "response"
ele_dam_df <- ele_dam_df %>%
  rename(LCL = asymp.LCL,
         UCL = asymp.UCL)

p_ele_dam <- ggplot(ele_dam_df,
                    aes(x = Treatment, y = response)) +
  jitter_layer(ele_dam_raw, rate_night) +
  geom_point(aes(fill = Treatment), size = 24, shape = 21, colour = "black", stroke = 1, show.legend = FALSE) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL, colour = Treatment),
                width = 0.05, linewidth = 8, show.legend = FALSE) +
  scale_fill_manual(values = treat_cols) +
  scale_colour_manual(values = treat_cols) +
  scale_x_discrete(labels = treat_labels) +
  labs(
    x = NULL,
    y = expression("Damage/Night (m"^2*")")
  ) +
  base_theme

############################
## CARNIVORE PLOTS
############################

## 5. Carnivore detection rate (per night) ----
emm_carn_rate <- emmeans(
  carn_detrate, ~ Treatment,
  type   = "response",
  offset = 0
)

carn_rate_df <- as.data.frame(emm_carn_rate)
# Standardize column names
colnames(carn_rate_df)[colnames(carn_rate_df) %in% c("rate", "response")] <- "response"
carn_rate_df <- carn_rate_df %>%
  rename(LCL = asymp.LCL,
         UCL = asymp.UCL)

p_carn_detect <- ggplot(carn_rate_df,
                        aes(x = Treatment, y = response)) +
  jitter_layer(carn_detect_raw, rate_night) +
  geom_point(aes(fill = Treatment), size = 24, shape = 21, colour = "black", stroke = 1, show.legend = FALSE) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL, colour = Treatment),
                width = 0.05, linewidth = 8, show.legend = FALSE) +
  scale_fill_manual(values = treat_cols) +
  scale_colour_manual(values = treat_cols) +
  scale_x_discrete(labels = treat_labels) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3), limits = c(0, NA)) +
  labs(x = NULL,
       y = "Detections/Night") +
  base_theme

## 6. Carnivore proximity (m) ----
emm_carn_prox <- emmeans(carnprox, ~ Treatment, type = "link")
carn_prox_df  <- as.data.frame(emm_carn_prox) %>%
  mutate(
    response = exp(emmean),
    LCL      = exp(asymp.LCL),
    UCL      = exp(asymp.UCL)
  )

p_carn_prox <- ggplot(carn_prox_df,
                      aes(x = Treatment, y = response)) +
  jitter_layer(carn_prox_raw, Proximity) +
  geom_point(aes(fill = Treatment), size = 24, shape = 21, colour = "black", stroke = 1, show.legend = FALSE) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL, colour = Treatment),
                width = 0.05, linewidth = 8, show.legend = FALSE) +
  scale_fill_manual(values = treat_cols) +
  scale_colour_manual(values = treat_cols) +
  scale_x_discrete(labels = treat_labels) +
  labs(x = NULL,
       y = "Mean Proximity (m)") +
  base_theme

## 7. Carnivore conversion (P(attack | present)) ----
emm_carn_conv <- emmeans(carn_conf_conv, ~ Treatment, type = "response")
carn_conv_df  <- as.data.frame(emm_carn_conv) %>%
  rename(
    response = prob,
    LCL      = asymp.LCL,
    UCL      = asymp.UCL
  )

p_carn_conv <- ggplot(carn_conv_df,
                      aes(x = Treatment, y = response)) +
  jitter_layer(carn_conv_raw, prop_attack) +
  geom_point(aes(fill = Treatment), size = 24, shape = 21, colour = "black", stroke = 1, show.legend = FALSE) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL, colour = Treatment),
                width = 0.05, linewidth = 8, show.legend = FALSE) +
  scale_fill_manual(values = treat_cols) +
  scale_colour_manual(values = treat_cols) +
  scale_x_discrete(labels = treat_labels) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits  = c(0, 1)) +
  labs(x = NULL,
       y = "P(Attack | Presence)") +
  base_theme

## 8. Carnivore damage rate (attacks per night) ----
emm_carn_dam <- emmeans(
  carn_damrate, ~ Treatment,
  type   = "response",
  offset = 0
)

carn_dam_df <- as.data.frame(emm_carn_dam)
# Standardize column names
colnames(carn_dam_df)[colnames(carn_dam_df) %in% c("rate", "response")] <- "response"
carn_dam_df <- carn_dam_df %>%
  rename(LCL = asymp.LCL,
         UCL = asymp.UCL)

p_carn_dam <- ggplot(carn_dam_df,
                     aes(x = Treatment, y = response)) +
  jitter_layer(carn_dam_raw, rate_night) +
  geom_point(aes(fill = Treatment), size = 24, shape = 21, colour = "black", stroke = 1, show.legend = FALSE) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL, colour = Treatment),
                width = 0.05, linewidth = 8, show.legend = FALSE) +
  scale_fill_manual(values = treat_cols) +
  scale_colour_manual(values = treat_cols) +
  scale_x_discrete(labels = treat_labels) +
  labs(x = NULL,
       y = "Attacks/Night") +
  base_theme 

############################
## PATCHWORK GRIDS
############################

# Elephants: 2×2 grid
elephant_panel <- (p_ele_detect | p_ele_prox) /
  (p_ele_conv   | p_ele_dam)

# Carnivores: 2×2 grid
carnivore_panel <- (p_carn_detect | p_carn_prox) /
  (p_carn_conv   | p_carn_dam)

# Poster: 3x2 grid
poster_panel <- (p_carn_detect | p_ele_detect) / 
                (p_carn_prox  | p_ele_prox) /
                (p_carn_dam  | p_ele_dam) 

poster_panel <- (p_carn_detect) / 
                (p_carn_prox) /
                (p_carn_dam) 

poster_panel <- (p_ele_detect) / 
                (p_ele_prox) /
                (p_ele_dam) 
# Display panels
elephant_panel
carnivore_panel
poster_panel

#############
#Habituation
#############
mrc_detections <- read.csv("solaroverlap.csv")
mrc_ele_det <- mrc_detections %>%
  filter(common_name == "African Elephant") %>%
  group_by(treatment, grid_id,week_id) %>% summarise(detections = n(), .groups = "drop")
full_weeks <- expand.grid(treatment = c("control", "human"),grid_id   = c("N2", "S2"), week_id   = 1:10)
mrc_ele_det <- full_weeks %>%
  left_join(mrc_ele_det, by = c("treatment", "grid_id", "week_id")) %>%
  mutate( detections = tidyr::replace_na(detections, 0))
mrc_ele_hab <- glmmTMB(
  detections ~ treatment + week_id + grid_id, family = nbinom2, data = mrc_ele_det)
summary(mrc_ele_hab)
res <- simulateResiduals(mrc_ele_hab)
plot(res)
testDispersion(res)
# Conditional model:
#                 Estimate  Std. Error  z value Pr(>|z|)   
# (Intercept)     2.15589    0.83498   2.582  0.00982 **
# treatmenthuman -0.73737    0.66530  -1.108  0.26772   
# week_id         0.01507    0.12188   0.124  0.90162   
# grid_idS2      -0.28841    0.68088  -0.424  0.67187   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

mrc_ele_human <- mrc_ele_det %>%
  filter(treatment == "human")

mrc_ele_hab_human <- glmmTMB(
  detections ~ week_id + (1 | grid_id),
  family = nbinom2,
  data = mrc_ele_human
)

summary(mrc_ele_hab_human)
# Conditional model:
# Groups  Name        Variance  Std.Dev. 
# grid_id (Intercept) 1.596e-09 3.995e-05
# Number of obs: 20, groups:  grid_id, 2
# 
# Dispersion parameter for nbinom2 family (): 0.241 
# 
# Conditional model:
# Estimate Std. Error z value Pr(>|z|)
# (Intercept)  1.09790    1.24577   0.881    0.378
# week_id      0.04905    0.20924   0.234    0.815

agc_hab <- read.csv("elecon2.csv")

agc_hab <- agc_hab %>%
  mutate(Date = dmy(Date),Farm = factor(Farm),Treatment = factor(Treatment)) %>%
  arrange(Farm, Date)

agc_hab <- agc_hab %>% group_by(Farm) %>%
  mutate(exposure_night = row_number()) %>% ungroup()

agc_hab <- agc_hab %>% mutate(ele_present_bin = as.integer(Ele.present == "Yes"))

agc_hab <- agc_hab %>% filter(Treatment == "Humans")

m_hab_pres <- glmmTMB(ele_present_bin ~ exposure_night + (1 | Farm),
                      family = binomial(link = "logit"), data = agc_hab)

summary(m_hab_pres)
# Conditional model:
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)     -1.455216   0.552893  -2.632  0.00849 ** 
# exposure_night   0.008584   0.011725   0.732  0.46410    
# TreatmentHumans -1.697449   0.385963  -4.398 1.09e-05 ***

# Conditional model:
# Groups Name        Variance  Std.Dev. 
# Farm   (Intercept) 1.066e-08 0.0001033
# Number of obs: 553, groups:  Farm, 19
# 
# Conditional model:
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    -2.32169    0.33834  -6.862 6.79e-12 ***
# exposure_night -0.03090    0.02178  -1.419    0.156  

agc_damage_area <- agc_hab %>%
  mutate(Area_m2 = as.numeric(Area.m2.)) %>% filter(Area_m2 > 0)

m_hab_damage_area <- glmmTMB(
  Area_m2 ~ exposure_night + Treatment + (1 | Farm),
  family = Gamma(link = "log"),  data = agc_damage_area)

summary(m_hab_damage_area)

m_hab_damage_area_human <- glmmTMB(
  Area_m2 ~ exposure_night + (1 | Farm),
  family = Gamma(link = "log"),
  data = agc_damage_area %>% filter(Treatment == "Humans")
)

summary(m_hab_damage_area_human)
# Conditional model:
# Groups Name        Variance Std.Dev.
# Farm   (Intercept) 0.1346   0.3668  
# Number of obs: 9, groups:  Farm, 8
# 
# Dispersion estimate for Gamma family (sigma^2): 0.266 
# 
# Conditional model:
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)     4.91710    0.33335   14.75   <2e-16 ***
# exposure_night -0.02639    0.02563   -1.03    0.303 


