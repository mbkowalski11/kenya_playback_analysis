#millet
millet <- read.csv("~/Desktop/UCSC/Research/Data/Kenya Playback Pilot/Food Measurements/Millet.csv")
head(millet)
library(tidyverse)
millet <- millet %>% select(Camera, Grid, Treatment, Week, Trial, Weight.Start,Weight.End, Consumed, Chewed, Poo)
millet2 <- millet %>% select(Camera, Grid, Treatment, Consumed)
ggplot(millet2, aes(x=Consumed, y=Treatment)) + 
  geom_violin()
bartlett.test(x = millet2$Consumed, g = millet2$Treatment) #p-value is above .05 so not enough evidence to reject the null that variances are equal (normal)
library(car)
leveneTest(Consumed ~ Treatment, millet2) #pvalue above .05 so we do not have enough evidence to reject the null of equal variances (non-normal)
qqPlot(millet2$Consumed)
qqnorm(millet2$Consumed)
qqline(millet2$Consumed)

millethuman <- millet2 %>% filter(Treatment=="human")
milletcontrol <- millet2 %>% filter(Treatment=="control")
wilcox.test(x = millethuman$Consumed, y = milletcontrol$Consumed, paired = TRUE, alternative = "less")
t.test(x = millethuman$Consumed, 
       y = milletcontrol$Consumed,
       alternative = "less",
       mu = 0, 
       paired = TRUE,   
       var.equal = TRUE,
       conf.level = 0.95)

sample <- c(colors()[375],colors()[142])

mplot <- ggplot(millet2, aes(x=Treatment, y=Consumed, fill = Treatment)) + 
  geom_jitter(width = .2, height = .2, aes(fill=Treatment), size = 3, shape = 21, alpha = 0.4, show.legend = FALSE) +
  geom_boxplot(width = 0.5, alpha=.6, outliers = FALSE, show.legend = FALSE) +
  scale_fill_manual(values = c("#FABA39FF", "#7A0403FF")) + 
  scale_x_discrete(labels = c("human" = "Human", "control" = "Control")) +
  labs(y = "Millet Consumed (g/tray)") +
  theme_bw(base_size=35)

diffs <- millethuman$Consumed - milletcontrol$Consumed
# Or density
plot(density(na.omit(diffs)), main = "Density of Differences", xlab = "Difference (Human - Control)")
abline(v = 0, lty = 2)

percent_change <- (millethuman$Consumed - milletcontrol$Consumed) / milletcontrol$Consumed * 100
percent_change <- percent_change[!is.na(percent_change)]
mean(percent_change)
log_ratio <- log1p(millethuman$Consumed) - log1p(milletcontrol$Consumed)
mean(log_ratio, na.rm = TRUE)
exp(mean(log_ratio, na.rm = TRUE))  # Average multiplicative change


#eggs
eggs <- read.csv("~/Desktop/UCSC/Research/Data/Kenya Playback Pilot/Food Measurements/Eggs.csv")
head(eggs)
ggplot(eggs, aes(x=Consumed, y=Treatment)) + 
  geom_violin()
leveneTest(Consumed ~ Treatment, eggs) #pvalue above .05 so we do not have enough evidence to reject the null of equal variances (non-normal)
qqPlot(eggs$Consumed)
eggshuman <- eggs %>% filter(Treatment=="human")
eggscontrol <- eggs %>% filter(Treatment=="control")
wilcox.test(x = eggshuman$Consumed, y = eggscontrol$Consumed, paired = TRUE, alternative = "less")
t.test(x = eggshuman$Consumed, 
       y = eggscontrol$Consumed,
       alternative = "less",
       mu = 0, 
       paired = TRUE,   
       var.equal = TRUE,
       conf.level = 0.95)

percent_change <- (eggshuman$Consumed - eggscontrol$Consumed) / eggscontrol$Consumed * 100
percent_change <- percent_change[!is.na(percent_change)]
log_ratio <- log1p(eggshuman$Consumed) - log1p(eggscontrol$Consumed)
mean(log_ratio, na.rm = TRUE)
exp(mean(log_ratio, na.rm = TRUE))

eplot <- ggplot(eggs, aes(x=Treatment, y=Consumed, fill = Treatment)) + 
  geom_jitter(width = .2, height = .2, aes(fill=Treatment), size = 3, shape = 21, alpha = 0.4, show.legend = FALSE) +
  geom_boxplot(width = 0.5, alpha=.6, outliers = FALSE, show.legend = FALSE) +
  scale_fill_manual(values = c("#FABA39FF", "#7A0403FF")) + 
  scale_x_discrete(labels = c("human" = "Human", "control" = "Control")) +
  scale_y_continuous(breaks = c(0,1,2), n.breaks = 3) +
  labs(y = "Eggs Consumed (n)") +
  theme_bw(base_size=35)

ggplot(eggs, aes(x=Consumed, y=Treatment, fill = Treatment)) + 
  geom_jitter(width = .2, height = .2, aes(fill=Treatment), size = 3, shape = 21, alpha = 0.4, show.legend = FALSE) +
  geom_violin(alpha=.6, show.legend = FALSE) +
  scale_fill_manual(values = c("#FABA39FF", "#7A0403FF")) + 
  scale_y_discrete(labels = c("human" = "Human", "control" = "Control")) +
  scale_x_continuous(breaks = c(0,1,2), n.breaks = 3) +
  labs(x = "Eggs Consumed (n)") +
  theme_bw(base_size=35)

ggplot(eggs, aes(x=as.factor(Consumed), fill = Treatment)) + 
  geom_bar(alpha=.6, width = .5, stat = "count", position = "dodge", color = "black") +
  scale_fill_manual(values = c("#FABA39FF", "#7A0403FF"),
                    labels = c("Control", "Human")) +
  scale_x_discrete() +
  labs(x = "Eggs Consumed", y = "Frequency (n)") +
  theme_bw(base_size=35) +
  theme(
    legend.position = c(0.8, 0.8), # Relative position inside the plot (x, y)
    legend.background = element_rect(fill = "white", color = "black", size = 0.3), # Optional styling
    legend.title = element_text(face = "bold") # Bold legend title
  )

#meso eggs
meggs <- read.csv("~/Desktop/UCSC/Research/Data/Kenya Playback Pilot/Meggs.csv")
head(meggs)
ggplot(meggs, aes(x=Consumed, y=Treatment)) + 
  geom_violin()
leveneTest(Consumed ~ Treatment, meggs) #pvalue above .05 so we do not have enough evidence to reject the null of equal variances (non-normal)
qqPlot(meggs$Consumed)
meggshuman <- meggs %>% filter(Treatment=="human")
meggscontrol <- meggs %>% filter(Treatment=="control")
wilcox.test(x = meggshuman$Consumed, y = meggscontrol$Consumed, paired = TRUE, alternative = "less")
t.test(x = meggshuman$Consumed, 
       y = meggscontrol$Consumed,
       alternative = "less",
       mu = 0, 
       paired = TRUE,   
       var.equal = TRUE,
       conf.level = 0.95)

#eles
eles <- read.csv("~/Desktop/UCSC/Research/Data/Kenya Playback Pilot/eleact.csv")
head(eles)
ggplot(eles, aes(x=Detections, y=Treatment)) + 
  geom_violin() + geom_jitter()
leveneTest(Detections ~ Treatment, eles) #pvalue above .05 so we do not have enough evidence to reject the null of equal variances (non-normal)
qqPlot(eles$Detections)
eleshuman <- eles %>% filter(Treatment=="human")
elescontrol <- eles %>% filter(Treatment=="control")
wilcox.test(x = eleshuman$Detections, y = elescontrol$Detections, paired = TRUE, alternative = "less")
t.test(x = eleshuman$Detections, 
       y = elescontrol$Detections,
       alternative = "less",
       mu = 0, 
       paired = TRUE,   
       var.equal = TRUE,
       conf.level = 0.95)
eleshuman
elescontrol
#Remove Week 5 for Lion Party
eles9 <- eles %>% filter(!Week==5)
eles9human <- eleshuman %>% filter(!Week==5)
eles9control <- elescontrol %>% filter(!Week==5)

ggplot(eles9, aes(x=Detections, y=Treatment)) + 
  geom_violin() + geom_jitter()
leveneTest(Detections ~ Treatment, eles9) #pvalue above .05 so we do not have enough evidence to reject the null of equal variances (non-normal)
qqPlot(eles$Detections)
wilcox.test(x = eles9human$Detections, y = eles9control$Detections, paired = TRUE, alternative = "less")
t.test(x = eles9human$Detections, 
       y = eles9control$Detections,
       alternative = "less",
       mu = 0, 
       paired = TRUE,   
       var.equal = TRUE,
       conf.level = 0.95)

sum(eles9human$Detections)
sum(eles9control$Detections)

ggplot(eles9, aes(x=Detections, y=Treatment, fill = Treatment)) + 
  geom_jitter(width = .2, height = .2, aes(fill=Treatment), size = 3, shape = 21, alpha = 0.4, show.legend = FALSE) +
  geom_boxplot(width = 0.5, alpha=.6, outliers = FALSE, show.legend = FALSE) +
  scale_fill_manual(values = c("gold", "indianred3")) + 
  scale_y_discrete(labels = c("human" = "Human", "control" = "Control")) +
  labs(x = "Detections (>30 min apart)") +
  theme_bw(base_size=35)

#camera numbers for space use
leveneTest(Camera ~ Treatment, eles9) #pvalue above .05 so we do not have enough evidence to reject the null of equal variances (non-normal)
qqPlot(eles$Camera)
wilcox.test(x = eles9human$Camera, y = eles9control$Camera, paired = TRUE, alternative = "less")
t.test(x = eles9human$Camera, 
       y = eles9control$Camera,
       alternative = "less",
       mu = 0, 
       paired = TRUE,   
       var.equal = FALSE,
       conf.level = 0.95)

sum(eles9human$Camera)
sum(eles9control$Camera)

ggplot(eles9, aes(x=Camera, y=Treatment, fill = Treatment)) + 
  geom_jitter(width = .2, height = .2, aes(fill=Treatment), size = 3, shape = 21, alpha = 0.4, show.legend = FALSE) +
  geom_boxplot(width = 0.5, alpha=.6, outliers = FALSE, show.legend = FALSE) +
  scale_fill_manual(values = c("gold", "indianred3")) + 
  scale_y_discrete(labels = c("human" = "Human", "control" = "Control")) +
  labs(x = "Unique Cameras Per Week") +
  theme_bw(base_size=35)

#tree damage
trees <- read.csv("~/Desktop/UCSC/Research/Data/Kenya Playback Pilot/eledamage2.csv")
head(trees)
str(trees)
as.factor(trees$damage)
as.factor(trees$species)
as.factor(trees$treatment)
as.factor(trees$grid)
ggplot(trees, aes(x=percent, y=treatment)) + 
  geom_violin()
leveneTest(percent ~ treatment, trees) #pvalue below .05 so normal
qqPlot(trees$percent)
treeshuman <- trees %>% filter(treatment=="human")
treescontrol <- trees %>% filter(treatment=="control")
wilcox.test(x = treeshuman$percent, y = treescontrol$percent, paired = FALSE, alternative = "two.sided")
t.test(x = treeshuman$percent, 
       y = treescontrol$percent,
       alternative = "less",
       mu = 0, 
       paired = FALSE,   
       var.equal = FALSE,
       conf.level = 0.95)

bigtrees <- trees %>% filter(percent>=3) #trees damaged more than 2%
ggplot(bigtrees, aes(x=percent, y=treatment)) + 
  geom_violin()
leveneTest(percent ~ treatment, bigtrees) #pvalue above .05 so not-normal
qqPlot(bigtrees$percent)
bigtreeshuman <- bigtrees %>% filter(treatment=="human")
bigtreescontrol <- bigtrees %>% filter(treatment=="control")
wilcox.test(x = bigtreeshuman$percent, y = bigtreescontrol$percent, paired = FALSE, alternative = "less")
t.test(x = bigtreeshuman$percent, 
       y = bigtreescontrol$percent,
       alternative = "less",
       mu = 0, 
       paired = FALSE,   
       var.equal = FALSE,
       conf.level = 0.95)

bigtrees.pf <- bigtrees %>%
  group_by(treatment) %>%
  summarise(mean = mean(percent))
bigtrees.pf

ggplot(bigtrees, aes(x=treatment, y=percent, fill = treatment)) + 
  geom_jitter(width = .2, height = .2, aes(fill=treatment), size = 3, shape = 21, alpha = 0.4, show.legend = FALSE) +
  geom_boxplot(width = 0.5, alpha=.6, outliers = FALSE, show.legend = FALSE) +
  scale_fill_manual(values = c("#FABA39FF", "#7A0403FF")) + 
  scale_x_discrete(labels = c("human" = "Human", "control" = "Control")) +
  labs(x = "Treatment", y = "Tree Damage (%)") +
  theme_bw(base_size=35)

# Calculate frequencies
trees.f <- trees %>%
  group_by(treatment, damage) %>%
  summarise(count = n(), .groups = "drop")

bigtrees.f <- bigtrees %>%
  group_by(treatment, damage) %>%
  summarise(count = n(), .groups = "drop")

trees.gf <- trees %>%
  group_by(treatment, grid, transect) %>%
  summarise(count = n(), .groups = "drop")
trees.gf.human <- trees.gf %>% filter(treatment=="human")
trees.gf.control <- trees.gf %>% filter(treatment=="control")
leveneTest(count ~ treatment, trees.gf) #normal data
qqPlot(trees.gf$count)

bigtrees.gf <- bigtrees %>%
  group_by(treatment, grid, transect) %>%
  summarise(count = n())
new_rows <- data.frame(
  treatment = c("human", "human"),
  grid = c("N2", "S2"),
  transect = c("b", "a"),
  count = c(0, 0))
bigtrees.gf <- bigtrees.gf %>% 
  bind_rows(new_rows)
bigtrees.gf.human <- trees.gf %>% filter(treatment=="human")
bigtrees.gf.control <- trees.gf %>% filter(treatment=="control")
leveneTest(count ~ treatment, bigtrees.gf) #normal data
qqPlot(bigtrees.gf$count)

wilcox.test(x = trees.gf.human$count, y = trees.gf.control$count, paired = TRUE, alternative = "less")
wilcox.test(x = bigtrees.gf.human$count, y = bigtrees.gf.control$count, paired = TRUE, alternative = "less")

t.test(x = trees.gf.human$count, 
       y = trees.gf.control$count,
       alternative = "less",
       mu = 0, 
       paired = TRUE,   
       var.equal = FALSE,
       conf.level = 0.95)

t.test(x = bigtrees.gf.human$count, 
       y = bigtrees.gf.control$count,
       alternative = "less",
       mu = 0, 
       paired = TRUE,   
       var.equal = FALSE,
       conf.level = 0.95)

#vertical plot
trees.gf$treatment <- factor(trees.gf$treatment, levels = c("control", "human"))
tplot <- ggplot(trees.gf, aes(y=count, x=treatment, fill = treatment)) + 
  geom_jitter(width = .2, height = .2, aes(fill=treatment), size = 3, shape = 21, alpha = 0.4, show.legend = FALSE) +
  geom_boxplot(width = 0.5, alpha=.6, outliers = FALSE, show.legend = FALSE) +
  scale_fill_manual(values = c("#FABA39FF", "#7A0403FF")) + 
  scale_x_discrete(labels = c("human" = "Human", "control" = "Control")) +
  scale_y_continuous(n.breaks = 6) +
  labs(y = "Damaged Trees (n)", x = "Treatment") +
  theme_bw(base_size=35)

#A different way to do the same repeated measures paired t-test for browsing frequency:
eledamage <- read.csv("eledamage2.csv")
transect_sum <- eledamage %>%
  group_by(transect, treatment) %>%
  summarise(
    n_browsed = sum(percent > 0, na.rm = TRUE),
    .groups = "drop"
  )

transect_wide_sum <- transect_sum %>%
  pivot_wider(names_from = treatment, values_from = n_browsed) %>%
  drop_na(control, human)

# Paired t-test: same transects, before (control) vs after (human)
t_test_sum <- t.test(transect_wide_sum$human, transect_wide_sum$control,
                     paired = TRUE, alternative = "less")

transect_wide_sum
t_test_sum



# Calculate relative frequencies
trees.rf <- trees %>%
  group_by(treatment, damage) %>%
  summarise(count = n()) %>%  # Count the number of occurrences
  mutate(relative.frequency = count / sum(count))  # Calculate relative frequency
trees.rf

bigtrees.rf <- bigtrees %>%
  group_by(treatment, damage) %>%
  summarise(count = n()) %>%  # Count the number of occurrences
  mutate(relative.frequency = count / sum(count))  # Calculate relative frequency
bigtrees.rf

# Function to run t-test for each damage type
damage.ttest <- function(damagetype) {
  # Subset data for the given damage type
  damagedata <- trees %>% filter(damage == damagetype)
  # Run t-test between 'human' and 'control' treatments for percent damage
  t_test_result <- t.test(percent ~ treatment, data = damagedata, var.equal = TRUE) # assuming equal variances
  print(paste("T-test for damage type:", damagetype))
  print(t_test_result)
}
damage.ttest(damagetype = "lb")

bigdamage.ttest <- function(damagetype) {
  # Subset data for the given damage type
  damagedata <- bigtrees %>% filter(damage == damagetype)
  # Run t-test between 'human' and 'control' treatments for percent damage
  t_test_result <- t.test(percent ~ treatment, data = damagedata, var.equal = TRUE) # assuming equal variances
  print(paste("T-test for damage type:", damagetype))
  print(t_test_result)
}
bigdamage.ttest(damagetype = "lb")





##############Plotting summaries
# Millet plot with mean and 95% CI
mplot <- ggplot(millet2, aes(x=Treatment, y=Consumed, fill = Treatment)) + 
  geom_jitter(width = .2, height = .2, aes(fill=Treatment), size = 24, shape = 21, alpha = 0.3, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "point", size = 24, shape = 21, color = "black", show.legend = FALSE) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.05, size = 4, aes(color = Treatment), show.legend = FALSE) +
  scale_fill_manual(values = c("#FABA39FF", "#7A0403FF")) + 
  scale_color_manual(values = c("#FABA39FF", "#7A0403FF")) +
  scale_y_continuous(breaks = c(0, 4, 8, 12)) +
  coord_cartesian(ylim = c(0, 12)) +
  scale_x_discrete(labels = c("human" = "H", "control" = "C")) +
  labs(y = "Millet Consumed (g/tray)", x = NULL) +
  theme_bw(base_size=35) +
  theme(text = element_text(family = "Palatino", size = 180),
        axis.text = element_text(size = 180, color = "black"),
        axis.text.y = element_text(size = 180, color = "black", margin = margin(r = 5)),
        axis.text.x = element_text(size = 180, color = "black", margin = margin(t = 25)),
        axis.title = element_text(size = 180),
        panel.grid = element_line(color = "grey90"),
        panel.border = element_rect(linewidth = 5, color = "black")
  )

# Eggs plot with mean and 95% CI
eplot <- ggplot(eggs, aes(x=Treatment, y=Consumed, fill = Treatment)) + 
  geom_jitter(width = .2, height = .1, aes(fill=Treatment), size = 24, shape = 21, alpha = 0.3, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "point", size = 24, shape = 21, color = "black", show.legend = FALSE) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.05, size = 4, aes(color = Treatment), show.legend = FALSE) +
  scale_fill_manual(values = c("#FABA39FF", "#7A0403FF")) + 
  scale_color_manual(values = c("#FABA39FF", "#7A0403FF")) +
  scale_x_discrete(labels = c("human" = "H", "control" = "C")) +
  scale_y_continuous(breaks = c(0,1,2)) +
  coord_cartesian(ylim = c(-.05, 2.05)) +
  labs(y = "Eggs Consumed (n)", x = NULL) +
  theme_bw(base_size=35) +
  theme(text = element_text(family = "Palatino", size = 180),
        axis.text = element_text(size = 180, color = "black"),
        axis.text.y = element_text(size = 180, color = "black"),
        axis.text.x = element_text(size = 180, color = "black", margin = margin(t = 25)),
        axis.title = element_text(size = 180),
        panel.grid = element_line(color = "grey90"),
        panel.border = element_rect(linewidth = 5, color = "black")
  )

# Trees plot with mean and 95% CI
tplot <- ggplot(trees.gf, aes(y=count, x=treatment, fill = treatment)) + 
  geom_jitter(width = .2, height = .2, aes(fill=treatment), size = 24, shape = 21, alpha = 0.3, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "point", size = 24, shape = 21, color = "black", show.legend = FALSE) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.05, size = 4, aes(color = treatment), show.legend = FALSE) +
  scale_fill_manual(values = c("#FABA39FF", "#7A0403FF")) + 
  scale_color_manual(values = c("#FABA39FF", "#7A0403FF")) +
  scale_x_discrete(labels = c("human" = "H", "control" = "C")) +
  scale_y_continuous(breaks = c(0,15,30,45,60)) +
  coord_cartesian(ylim = c(0, 60)) +
  labs(y = "Damaged Trees (n)", x = NULL) +
  theme_bw(base_size=35) + 
  theme(text = element_text(family = "Palatino", size = 180),
        axis.text = element_text(size = 180, color = "black"),
        axis.text.y = element_text(size = 180, color = "black"),
        axis.text.x = element_text(size = 180, color = "black", margin = margin(t = 25)),
        axis.title = element_text(size = 180),
        panel.grid = element_line(color = "grey90"),
        panel.border = element_rect(linewidth = 5, color = "black")
  )

# Trees plot with mean and 95% CI
dplot <- ggplot(bigtrees, aes(y=percent, x=treatment, fill = treatment)) + 
  geom_jitter(width = .2, height = .2, aes(fill=treatment), size = 24, shape = 21, alpha = 0.3, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "point", size = 24, shape = 21, color = "black", show.legend = FALSE) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.05, size = 4, aes(color = treatment), show.legend = FALSE) +
  scale_fill_manual(values = c("#FABA39FF", "#7A0403FF")) + 
  scale_color_manual(values = c("#FABA39FF", "#7A0403FF")) +
  scale_x_discrete(labels = c("human" = "H", "control" = "C")) +
  scale_y_continuous(breaks = c(0,25,50,75,100)) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(y = "Tree Damage (%)", x = NULL) +
  theme_bw(base_size=35) + 
  theme(text = element_text(family = "Palatino", size = 180),
        axis.text = element_text(size = 180, color = "black"),
        axis.text.y = element_text(size = 180, color = "black"),
        axis.text.x = element_text(size = 180, color = "black", margin = margin(t = 25)),
        axis.title = element_text(size = 180),
        panel.grid = element_line(color = "grey90"),
        panel.border = element_rect(linewidth = 5, color = "black")
  )

library(gridExtra)
splot <- grid.arrange(mplot, eplot, tplot, ncol = 3)
splot <- grid.arrange(tplot, dplot, ncol = 2)
splot <- grid.arrange(mplot, eplot, tplot, dplot, ncol = 4)


################
#Whitney-U Test for humans
anthro_camera <- combined %>%
  group_by(grid_id, camera_id) %>%
  summarise(
    total_occurrences = n(),
    .groups = "drop"
  )
anthro_camera <- anthro_camera %>%
  mutate(
    camera_nights = 10 * 7,   # 70 camera nights
    rate_per_night = total_occurrences / camera_nights
  )
table(anthro_camera$grid_id)
wilcox.test(
  rate_per_night ~ grid_id,
  data  = anthro_camera,
  exact = FALSE
)

