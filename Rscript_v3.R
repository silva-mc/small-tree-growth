###
### Code for "Small understory trees increase growth following sustained drought in the Amazon"
###

###
### Load packages
###

library(tidyverse)
library(lme4)
library(MuMIn)
library(lmerTest)
library(ggpubr)
library(factoextra)

###
### Load datasets
###

data = read.csv("Cax_Silva_M_data_v3.csv") # Small tree DBH timeseries
data_larg = read.csv("Cax_Silva_M_data_larg_v3.csv") # Large tree DBH timeseries
data_rgr = read.csv("Cax_Silva_M_data_rgr_v3.csv") # Individual-level stem increment
data_rgr_avg = read.csv("Cax_Silva_M_data_rgr_avg_v3.csv") # Species-averaged stem increment
data_subp = read.csv("Cax_Silva_M_data_subp_v3.csv") # Subplot-averaged stem increment
data_subp_long = read.csv("Cax_Silva_M_data_subp_long_v3.csv") # Subplot species richness and invidual density
data_trai_pca = read.csv("Cax_Silva_M_data_trai_pca_v3.csv") # Species-level stem increment, traits and PCA axes
data_trai_pca_delt = read.csv("Cax_Silva_M_data_trai_pca_delt_v3.csv") # Species-level delta stem increment, traits and PCA axes
pca_scores = read.csv("Cax_Silva_M_pca_scores_v3.csv") # PCA scores
pca_loadings = read.csv("Cax_Silva_M_pca_loadings_v3.csv") # PCA loadings
pca = readRDS("Cax_Silva_M_pca_v3.rds") # PCA

### 
### Process data
###

# Remove individuals outside the inclusion criteria
data = data %>%
  filter((dbh >= 1 & dbh <= 10) | is.na(dbh))
data_larg = data_larg %>%
  filter(dbh >= 10 | is.na(dbh))

# List shared species
shar_spec = intersect(unique(data_rgr[data_rgr$plot == "control", "search.str"]), # Species in the control plot
                      unique(data_rgr[data_rgr$plot == "tfe", "search.str"])) # Species in the tfe plot

# Dataset with RGR for shared species
shar_rgr = data_rgr_avg %>%
  filter(search.str %in% shar_spec) %>%
  spread(key = plot, value = rgr)

###
### Analyses
###

# Number of trees per plot
nrow(data_rgr[data_rgr$plot == "control",]) # Control
nrow(data_rgr[data_rgr$plot == "tfe",]) # TFE

# Number of dead individuals
nrow(data %>% filter(date == "2024-04-19", plot == "control", obs == "dead")) # Control
nrow(data %>% filter(date == "2024-04-19", plot == "tfe", obs == "dead")) # TFE

# Species richness and individual density vs time * plot
summary(aov(stem_numb ~ year * plot, data_subp_long))
summary(aov(rich ~ year * plot, data_subp_long))

# Stem increment vs plot
fit_rgr = lmer(rgr ~ plot + (1|genus/search.str), data_rgr)
summary(fit_rgr)
r.squaredGLMM(fit_rgr)

# Stem increment stats
mean(data_rgr[data_rgr$plot == "tfe", "rgr"], na.rm = T) # TFE
mean(data_rgr[data_rgr$plot == "control", "rgr"], na.rm = T) # Control
sd(data_rgr[data_rgr$plot == "tfe", "rgr"], na.rm = T) # TFE
sd(data_rgr[data_rgr$plot == "control", "rgr"], na.rm = T) # Control
quantile(data_rgr[data_rgr$plot == "tfe", "rgr"], probs = .95, na.rm = T) # TFE
quantile(data_rgr[data_rgr$plot == "control", "rgr"], probs = .95, na.rm = T) # Control

# Stem increment vs initial DBH
summary(lmer(rgr ~ dbh_init * plot + (1|genus/search.str), data_rgr))

# Stem increment vs plot (generalist species)
t.test(x = shar_rgr[shar_rgr$search.str %in% shar_spec, ][["control"]],
       y = shar_rgr[shar_rgr$search.str %in% shar_spec, ][["tfe"]],
       paired = T)

# Mean increment (generalist species)
summary(shar_rgr)

# Stem increment generalist vs specialist (Control)
t.test(rgr ~ occur,
       data_rgr_avg %>%
         mutate(occur = if_else(search.str %in% shar_spec, "Generalist", "Specialist")) %>%
         filter(plot == "control"), var.equal = F)

# Stem increment generalist vs specialist (TFE)
t.test(rgr ~ occur,
       data_rgr_avg %>%
         mutate(occur = if_else(search.str %in% shar_spec, "Generalist", "Specialist")) %>%
         filter(plot == "tfe"), var.equal = F)

# Check PCA summary
summary(pca)

# Check PCA scores
pca_scores

# Check PCA loading
pca_loadings

# Trait vs PCs
print(
  `dimnames<-`(
    outer(
      c("vcmax", "jmax", "rleaf", "gsmin", "nmass", "pmass", "lma", "lth", 
        "wd", "pdwp", "mdwp", "p50", "p88", "ksmax", "plc", "lasa"),
      c("PC1", "PC2", "PC3"),
      Vectorize(function(trait, pc) {
        test <- cor.test(data_trai_pca[[trait]], data_trai_pca[[pc]], use = "pairwise.complete.obs")
        paste0("r = ", formatC(test$estimate, digits = 2, format = "f"),
               ", p = ", format.pval(test$p.value, digits = 3, eps = 0.0001))
      })
    ),
    list(
      c("vcmax", "jmax", "rleaf", "gsmin", "nmass", "pmass", "lma", "lth", 
        "wd", "pdwp", "mdwp", "p50", "p88", "ksmax", "plc", "lasa"),
      c("PC1", "PC2", "PC3")
    )
  )
)

# Stem increment vs PCs
fit_glob = lm(rgr ~ (PC1 + PC2 + PC3) * plot, data_trai_pca, na.action = na.fail) # Global model
dred = dredge(fit_glob) # Model selection
fit_best = get.models(dred, subset = 1)[[1]] # Best model
summary(fit_best)

# Delta stem increment vs delta PCs
data_long = data_trai_pca %>% # Pivot long
  rename(PC_rgr = rgr) %>%
  pivot_longer(cols = starts_with("PC"),
               names_to = "PC",
               values_to = "value") %>%
  mutate(PC_plot = paste0(PC, "_", plot))
data_wide = data_long %>% pivot_wider(id_cols = c(family, genus, search.str), # Pivot wide
                                      names_from = PC_plot,
                                      values_from = value) %>% 
  drop_na() 
data_trai_pca_delt = data_wide %>% # Calculate delta
  rename("rgr_control" = "PC_rgr_control",
         "rgr_tfe" = "PC_rgr_tfe") %>%
  mutate(across(
    .cols = ends_with("_tfe"),
    .fns = ~ .x - get(str_replace(cur_column(), "_tfe", "_control")),
    .names = "{str_replace(.col, '_tfe', '_delt')}"))
fit_glob_delt = lm(rgr_delt ~ (PC1_delt + PC2_delt + PC3_delt), # Global model
                   data_trai_pca_delt, na.action = na.fail)
dred_delt = dredge(fit_glob_delt) # Model selection
fit_best_delt = get.models(dred_delt, subset = 1)[[1]] # Best model
summary(fit_best_delt)
rm(data_wide, data_long) # Remove intermediary objects

# Stem increment stats (subplot)
mean(data_subp[data_subp$plot == "control", "rgr_mean"]) # Control
mean(data_subp[data_subp$plot == "tfe", "rgr_mean"]) # TFE
sd(data_subp[data_subp$plot == "control", "rgr_mean"]) # Control
sd(data_subp[data_subp$plot == "tfe", "rgr_mean"]) # TFE

# Stem increment vs small tree density * plot
summary(lm(rgr_mean ~ sapli_numbe * plot, data_subp))

# Stem increment vs large tree density * plot
summary(lm(rgr_mean ~ tree_numbe * plot, data_subp))

# Stem increment vs (small + large tree density) * plot
summary(lm(rgr_mean ~ (sapli_numbe + tree_numbe) * plot, data_subp))

# Stem increment vs small tree basal area * plot
summary(lm(rgr_mean ~ sapli_ba_m2_m2 * plot, data_subp %>% mutate(sapli_ba_m2_m2 = sapli_ba_m2_m2 * 10000)))

# Stem increment vs large tree basal area * plot
summary(lm(rgr_mean ~ tree_ba_m2_m2 * plot, data_subp %>% mutate(tree_ba_m2_m2 = tree_ba_m2_m2 * 10000)))

# Large tree density vs small tree density (Control)
cor.test(~ tree_densi + sapli_densi, data_subp %>% filter(plot == "control"))

# Large tree density vs small tree density (TFE)
cor.test(~ tree_densi + sapli_densi, data_subp %>% filter(plot == "tfe"))

# Large tree basal area vs small tree basal area (Control)
cor.test(~ tree_ba_m2_m2 + sapli_ba_m2_m2, data_subp %>% filter(plot == "control"))

# Large tree basal area vs small tree basal area (Control)
cor.test(~ tree_ba_m2_m2 + sapli_ba_m2_m2, data_subp %>% filter(plot == "tfe"))

###
### Main text figures
###

# Figure 1

# Initial year (2017)
plot1 = ggplot() +
  geom_histogram(data = data[data$date == "2017-05-22" & data$plot == "control", ], 
                 aes(x = dbh, fill = "Control", y = ..count..), colour = NA, alpha = .8, breaks = seq(1, 10, .5)) +
  geom_histogram(data = data[data$date == "2017-05-19" & data$plot == "tfe", ], 
                 aes(x = dbh, fill = "TFE", y = -..count..), colour = NA, alpha = .8, breaks = seq(1, 10, .5)) +
  scale_fill_viridis_d(option = "turbo", begin = .2, end = .8, labels = c("Control", "TFE")) +
  scale_y_continuous(labels = abs, limits = c(-35, 115)) +
  labs(x = expression("DBH" ~ ("cm")),
       y = expression("Count"),
       title = "2017") +
  coord_flip() +
  theme_classic() +
  theme(legend.title = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12, colour = "black"))

# Final year (2024)
plot2 = ggplot() +
  geom_histogram(data = data[data$date == "2024-04-19" & data$plot == "control", ], 
                 aes(x = dbh, fill = "Control", y = ..count..), colour = NA, alpha = .8, breaks = seq(1, 10, .5)) +
  geom_histogram(data = data[data$date == "2024-04-19" & data$plot == "tfe", ], 
                 aes(x = dbh, fill = "TFE", y = -..count..), colour = NA, alpha = .8, breaks = seq(1, 10, .5)) +
  scale_fill_viridis_d(option = "turbo", begin = .2, end = .8, labels = c("Control", "TFE")) +
  scale_y_continuous(labels = abs, limits = c(-35, 115)) +
  labs(x = expression("DBH" ~ ("cm")),
       y = expression("Count"),
       title = "2024") +
  coord_flip() +
  theme_classic() +
  theme(legend.title = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12, colour = "black"))

# Dead final year (2024)
plot3 = ggplot() +
  geom_histogram(data = data[data$date == "2024-04-19" & data$plot == "control" & data$obs == "dead", ], 
                 aes(x = dbh_init, fill = "Control", y = ..count..), colour = NA, alpha = .8, breaks = seq(1, 10, .5)) +
  geom_histogram(data = data[data$date == "2024-04-19" & data$plot == "tfe" & data$obs == "dead", ], 
                 aes(x = dbh_init, fill = "TFE", y = -..count..), colour = NA, alpha = .8, breaks = seq(1, 10, .5)) +
  scale_fill_viridis_d(option = "turbo", begin = .2, end = .8, labels = c("Control", "TFE")) +
  scale_y_continuous(labels = abs) +
  labs(x = expression("DBH" ~ ("cm")),
       y = expression("Count"),
       title = "Dead") +
  coord_flip() +
  theme_classic() +
  theme(legend.title = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12, colour = "black"))

jpeg("Fig 1.jpg", width = 30, height = 10, units = "cm", res = 1000)

ggarrange(plot1, plot2, plot3,
          ncol = 3, nrow = 1,
          labels = "auto", common.legend = T, legend = "right")

dev.off()

# Clean environment
rm(plot1, plot2, plot3)

# Figure 2

# Tree density
plot1 = ggplot(data_subp_long, aes(x = factor(year), y = stem_numb, fill = plot)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  scale_fill_viridis_d(option = "turbo", begin = .2, end = .8, labels = c("Control", "TFE")) +
  labs(x = element_blank(),
       y = expression("Tree density " ~ ("ind." ~ subplot^-1))) +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text = element_text(size = 10, colour = "black")) +
  annotate("text", x = Inf, y = 52, hjust = 1, vjust = 0, size = 3.5,
           label = "Year~'ns'", parse = T) +
  annotate("text", x = Inf, y = 49, hjust = 1, vjust = 0, size = 3.5,
           label = "bold('Plot'[TFE])~bold('***')", parse = T) +
  annotate("text", x = Inf, y = 46, hjust = 1, vjust = 0, size = 3.5,
           label = "'Year:Plot'[TFE]~'ns'", parse = T)

# Tree richness
plot2 = ggplot(data_subp_long, aes(x = factor(year), y = rich, fill = plot)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  scale_fill_viridis_d(option = "turbo", begin = .2, end = .8, labels = c("Control", "TFE")) +
  labs(x = element_blank(),
       y = expression("Species richness " ~ ("spp." ~ subplot^-1))) +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text = element_text(size = 10, colour = "black")) +
  annotate("text", x = Inf, y = 38, hjust = 1, vjust = 0, size = 3.5,
           label = "Year~'ns'", parse = T) +
  annotate("text", x = Inf, y = 35.5, hjust = 1, vjust = 0, size = 3.5,
           label = "bold('Plot'[TFE])~bold('***')", parse = T) +
  annotate("text", x = Inf, y = 33, hjust = 1, vjust = 0, size = 3.5,
           label = "'Year:Plot'[TFE]~'ns'", parse = T)

# Export it
jpeg("Fig 2.jpg", width = 22.5, height = 10, units = "cm", res = 1000)

ggarrange(plot1, plot2,
          ncol = 2, nrow = 1,
          labels = "auto",
          common.legend = T, legend = "right")

dev.off()

# Clean environment
rm(plot1, plot2)

# Figure 3

# General
plot1 = ggplot(data_rgr,
               aes(x = plot, y = rgr, fill = plot)) +
  geom_boxplot(outlier.alpha = .1, colour = "black") +
  scale_fill_viridis_d(option = "turbo", begin = .2, end = .8, labels = c("control", "tfe")) +
  scale_colour_viridis_d(option = "turbo", begin = .2, end = .8, labels = c("control", "tfe")) +
  scale_x_discrete(labels = c("Control", "TFE")) +
  ylim(0, 1) +
  labs(x = element_blank(),
       y = expression("Stem increment" ~ (cm ~ year^-1)),
       fill = element_blank()) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text = element_text(size = 12, colour = "black")) +
  annotate("text", x = -Inf, y = Inf, hjust = 0, vjust = 1, size = 5,
           label = " ***", fontface = "bold")

# RGR generalist species
plot2 = ggplot(data_rgr_avg %>%
                 filter(search.str %in% shar_spec), 
               aes(x = plot, y = rgr, fill = plot)) +
  geom_boxplot(outlier.alpha = 0, colour = "black") +
  geom_line(aes(group = search.str), alpha = .05) +
  geom_point(shape = 21, alpha = .05) +
  scale_fill_viridis_d(option = "turbo", begin = .2, end = .8, labels = c("Control", "TFE")) +
  scale_x_discrete(labels = c("Control", "TFE")) +
  ylim(0, 1) +
  labs(x = element_blank(),
       y = expression("Stem increment" ~ (cm ~ year^-1)),
       fill = element_blank()) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text = element_text(size = 12, colour = "black")) +
  annotate("text", x = -Inf, y = Inf, hjust = 0, vjust = 1, size = 5,
           label = " ***", fontface = "bold")

# Specialist vs. generalist RGR
plot3 = ggplot(data_rgr_avg %>%
                 mutate(occur = if_else(search.str %in% shar_spec, "Shared", "Exclusive")),
               aes(x = occur, y = rgr, fill = plot)) +
  geom_boxplot(outlier.alpha = .1) +
  scale_fill_viridis_d(option = "turbo", begin = .2, end = .8, labels = c("Control", "TFE")) +
  scale_x_discrete(limits = c("Shared", "Exclusive")) +
  ylim(0, 1) +
  labs(x = element_blank(),
       y = expression("Stem increment" ~ (cm ~ year^-1))) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12, colour = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size = 16)) +
  facet_grid(. ~ plot, labeller = as_labeller(c(
    "control" = "Control",
    "tfe" = "TFE"
  ))) +
  annotate("text", x = -Inf, y = Inf, hjust = 0, vjust = 1, size = 4,
           label = " ns")

# Export it
jpeg("Fig 3.jpg", width = 30, height = 9, units = "cm", res = 1000)

ggarrange(plot1, plot2, plot3,
          ncol = 3, nrow = 1,
          labels = "auto")

dev.off()

# Clean environment
rm(plot1, plot2, plot3)

# Figure 4

# Create expression labels
labe = c(vcmax  = "V[cmax]",
         jmax   = "J[max]",
         rleaf  = "R[dark]",
         gsmin  = "g[min]",
         nmass  = "N[mass]",
         pmass  = "P[mass]",
         lma    = "LMA",
         lth    = "L[th]",
         wd     = "rho",
         pdwp   = "Psi[pd]",
         mdwp   = "Psi[md]",
         p50    = "P50",
         p88    = "P88",
         ksmax  = "K[smax]",
         plc    = "PLC",
         lasa   = "A[l]:A[s]")

# PC1 vs PC2 
plot1 = fviz_pca_var(pca,
                     col.var = "contrib", # color by contribution to PCs
                     geom = c("arrow", "text"),
                     col.circle = NA,
                     axes.linetype = "blank",
                     labelsize = 5,
                     repel = T) + # avoid overlapping text labels
  scale_colour_viridis_c(option = "rocket", direction = -1, limits = c(0, 16.6)) +
  labs(x = "Photosynthetic potential (PC1 22.5%)",
       y = "Nutrient content (PC2 14.6%)",
       colour = "Contribution (%)",
       title = NULL) +
  theme_classic() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12, colour = "black"))

# Replace labels in the data layer used for text
text_p1 = which(sapply(plot1$layers, function(l) "label" %in% names(l$mapping)))
plot1$layers[[text_p1]]$mapping$label <- 
  aes(label = labe[as.character(plot1$data$name)])$label
plot1$layers[[text_p1]]$geom_params$parse <- TRUE

# PC1 vs PC3
plot2 = fviz_pca_var(pca,
                     axes = c(1, 3),
                     col.var = "contrib", # color by contribution to PCs
                     geom = c("arrow", "text"),
                     col.circle = NA,
                     axes.linetype = "blank",
                     labelsize = 5,
                     repel = T) + # avoid overlapping text labels
  scale_colour_viridis_c(option = "rocket", direction = -1, limits = c(0, 16.6)) +
  labs(x = "Photosynthetic potential (PC1 22.5%)",
       y = "Embolism resistance (PC3 13.1%)",
       colour = "Contribution (%)",
       title = NULL) +
  theme_classic() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12, colour = "black"))

# Replace labels in the data layer used for text
text_p2 <- which(sapply(plot2$layers, function(l) "label" %in% names(l$mapping)))
plot2$layers[[text_p2]]$mapping$label <- 
  aes(label = labe[as.character(plot1$data$name)])$label
plot2$layers[[text_p2]]$geom_params$parse <- T

# PC2 vs PC3
plot3 = fviz_pca_var(pca,
                     axes = c(2, 3),
                     col.var = "contrib", # color by contribution to PCs
                     geom = c("arrow", "text"),
                     col.circle = NA,
                     axes.linetype = "blank",
                     labelsize = 5,
                     repel = T) + # avoid overlapping text labels
  scale_colour_viridis_c(option = "rocket", direction = -1, limits = c(0, 16.6)) +
  labs(x = "Nutrient content (PC2 14.6%)",
       y = "Embolism resistance (PC3 13.1%)",
       colour = "Contribution (%)",
       title = NULL) +
  theme_classic() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12, colour = "black"))

# Replace labels in the data layer used for text
text_p3 <- which(sapply(plot3$layers, function(l) "label" %in% names(l$mapping)))
plot3$layers[[text_p3]]$mapping$label <- 
  aes(label = labe[as.character(plot1$data$name)])$label
plot3$layers[[text_p3]]$geom_params$parse <- TRUE

# Export it
jpeg("Fig 4.jpg", width = 37.5, height = 12.5, units = "cm", res = 1000)

ggarrange(plot1, plot2, plot3, 
          ncol = 3, nrow = 1,
          labels = "auto",
          common.legend = T, legend = "right")

dev.off()

rm(plot1, plot2, plot3,
   text_p1, text_p2, text_p3, labe)

# Figure 5

# Absolute growth
plot1 = ggplot(data = data_trai_pca, aes(x = PC2, y = rgr, colour = plot, fill = plot)) + 
  geom_smooth(data = subset(data_trai_pca, plot == "tfe"), 
              method = "lm", alpha = 0.5, show.legend = F) +
  geom_point(alpha = .5, size = 3) +
  labs(x = "Nutrient content (unitless)",
       y = expression("Stem increment" ~ (cm ~ year^-1)),
       colour = element_blank(),
       fill = element_blank()) +
  scale_colour_viridis_d(option = "turbo", begin = .2, end = .8, labels = c("Control", "TFE")) +
  scale_fill_viridis_d(option = "turbo", begin = .2, end = .8, labels = c("Control", "TFE")) +
  annotate("text", x = -Inf, y = .42, hjust = 0, vjust = 1, size = 3.5,
           label = "' Nutrient content'~'ns'", parse = T) +
  annotate("text", x = -Inf, y = .39, hjust = 0, vjust = 1, size = 3.5,
           label = "bold(' Plot'[TFE])~bold('**')", parse = T) +
  annotate("text", x = -Inf, y = .36, hjust = 0, vjust = 1, size = 3.5,
           label = "bold(' Nutrient content:Plot'[TFE])~bold('*')", parse = T) +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12, colour = "black"))

# Delta growth
plot2 = ggplot(data = data_trai_pca_delt, aes(x = PC2_delt, y = rgr_delt)) + 
  geom_smooth(method = "lm", alpha = .5, colour = "black") +
  geom_point(size = 3) +
  labs(x = expression(Delta * "Nutrient content (unitless)"),
       y = expression(Delta * "Stem increment" ~ (cm ~ year^-1))) +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12, colour = "black"))  +
  annotate("text", x = -Inf, y = Inf, hjust = 0, vjust = 1, size = 5,
           label = " **", fontface = "bold")

# Export it
jpeg("Fig 5.jpg", width = 22.5, height = 10, units = "cm", res = 1000)

ggarrange(plot1, plot2,
          ncol = 2, nrow = 1,
          labels = "auto",
          common.legend = F, widths = c(5, 4))

dev.off()

rm(plot1, plot2)

# Figure 6

# Sapling density
plot1 = ggplot(data = data_subp, aes(x = sapli_numbe, y = rgr_mean, colour = plot, fill = plot)) + 
  geom_abline(slope = 0, intercept = 0, linetype = "dashed") +
  geom_smooth(data = subset(data_subp, plot == "tfe"), 
              method = "lm", alpha = 0.5, show.legend = F) +
  geom_point(alpha = .5, size = 3) +
  labs(x = "Small tree density" ~ ("ind." ~ subplot^-1),
       y = expression("Stem increment" ~ (cm ~ year^-1)),
       colour = element_blank(),
       fill = element_blank()) +
  scale_colour_viridis_d(option = "turbo", begin = .2, end = .8, labels = c("Control", "TFE")) +
  scale_fill_viridis_d(option = "turbo", begin = .2, end = .8, labels = c("Control", "TFE")) +
  annotate("text", x = Inf, y = 0.43, hjust = 1, vjust = 0, size = 3.5,
           label = "Density~'ns'", parse = T) +
  annotate("text", x = Inf, y = 0.405, hjust = 1, vjust = 0, size = 3.5,
           label = "bold('Plot'[TFE])~bold('***')", parse = T) +
  annotate("text", x = Inf, y = 0.38, hjust = 1, vjust = 0, size = 3.5,
           label = "bold('Density:Plot'[TFE])~bold('**')", parse = T) +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12, colour = "black"))

# Tree density
plot2 = ggplot(data = data_subp, aes(x = tree_numbe, y = rgr_mean, colour = plot, fill = plot)) + 
  geom_abline(slope = 0, intercept = 0, linetype = "dashed") +
  # geom_smooth(method = "lm", alpha = .5) +
  geom_point(alpha = .5, size = 3) +
  labs(x = "Large tree density" ~ ("ind." ~ subplot^-1),
       y = expression("Stem increment" ~ (cm ~ year^-1)),
       colour = element_blank(),
       fill = element_blank()) +
  scale_colour_viridis_d(option = "turbo", begin = .2, end = .8, labels = c("Control", "TFE")) +
  scale_fill_viridis_d(option = "turbo", begin = .2, end = .8, labels = c("Control", "TFE")) +
  annotate("text", x = Inf, y = 0.43, hjust = 1, vjust = 0, size = 3.5,
           label = "Density~'ns'", parse = T) +
  annotate("text", x = Inf, y = 0.405, hjust = 1, vjust = 0, size = 3.5,
           label = "bold('Plot'[TFE])~bold('**')", parse = T) +
  annotate("text", x = Inf, y = 0.38, hjust = 1, vjust = 0, size = 3.5,
           label = "'Density:Plot'[TFE]~'ns'", parse = T) +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12, colour = "black"))

# Export it
jpeg("Fig 6.jpg", width = 22.5, height = 10, units = "cm", res = 1000)

ggarrange(plot1, plot2, 
          ncol = 2, nrow = 1,
          labels = "auto",
          common.legend = T, legend = "right")

dev.off()

rm(plot1, plot2)

####
#### SI Figures
####

# Figure S4

# Large tree size distribution
plot1 = ggplot() +
  geom_histogram(data = data_larg[data_larg$plot == "control", ], 
                 aes(x = dbh, fill = "Control", y = ..count..), colour = NA, alpha = .8, breaks = seq(10, 180, 2)) +
  geom_histogram(data = data_larg[data_larg$plot == "tfe", ], 
                 aes(x = dbh, fill = "TFE", y = -..count..), colour = NA, alpha = .8, breaks = seq(10, 180, 2)) +
  scale_fill_viridis_d(option = "turbo", begin = .2, end = .8, labels = c("Control", "TFE")) +
  scale_y_continuous(labels = abs) +
  labs(x = expression("DBH" ~ ("cm")),
       y = expression("Count")) +
  coord_flip() +
  theme_classic() +
  theme(legend.title = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12, colour = "black"))

# Combined size distribution
data_weig = bind_rows(data %>% filter(date == "2024-04-19") %>% select(plot, dbh) %>%
                        mutate(area_ha = 20 * 10 * 10 / 10000, # 20 subplots x 100 m2 each / 10000 m2/ha
                               weight = 1 / area_ha, size = "Small"),
                      data_larg %>% select(plot, dbh) %>%
                        mutate(area_ha = 1, # 1 ha plot
                               weight = 1 / area_ha, size = "Large"))
plot2 = ggplot(data_weig) +
  geom_histogram(data = data_weig %>% filter(plot == "control"),
                 aes(x = dbh, fill = "Control", weight = weight),
                 colour = NA, alpha = 0.8, breaks = seq(0, 180, 2)) +
  geom_histogram(data = data_weig %>% filter(plot == "tfe"),
                 aes(x = dbh, fill = "TFE", weight = -weight),
                 colour = NA, alpha = 0.8, breaks = seq(0, 180, 2)) +
  scale_fill_viridis_d(option = "turbo", begin = 0.2, end = 0.8, labels = c("Control", "TFE")) +
  scale_y_continuous(labels = abs) +
  labs(x = expression("DBH" ~ ("cm")),
       y = expression("Count")) +
  coord_flip() +
  theme_classic() +
  theme(legend.title = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12, colour = "black"))

# Export it
jpeg("Fig S4.jpg", width = 20, height = 10, units = "cm", res = 1000)

ggarrange(plot1, plot2, 
          ncol = 2, nrow = 1,
          labels = "auto",
          common.legend = T, legend = "right")

dev.off()

rm(plot1, plot2, data_weig)

# Figure S5

jpeg("Fig S5.jpg", width = 12.5, height = 10, units = "cm", res = 1000)

ggplot(data_rgr, aes(x = dbh_init, y = rgr, colour = plot, fill = plot)) +
  geom_point(alpha = .5, size = 3) +
  geom_smooth(data = data_rgr, aes(x = dbh_init, y = rgr), colour = "black", fill = "black", method = "lm") +
  labs(x = expression("DBH"["2017"] ~ ("cm")),
       y = expression("Stem increment" ~ (cm ~ year^-1))) +
  scale_colour_viridis_d(option = "turbo", begin = .2, end = .8, labels = c("Control", "TFE")) +
  scale_fill_viridis_d(option = "turbo", begin = .2, end = .8, labels = c("Control", "TFE")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text = element_text(size = 10, colour = "black")) +
  annotate("text", x = -Inf, y = 1.1, hjust = 0, vjust = 1, size = 3.5,
           label = "~bold(DBH[2017]) ~ bold('***')", parse = T) +
  annotate("text", x = -Inf, y = 1, hjust = 0, vjust = 1, size = 3.5,
           label = "~bold(Plot[TFE]) ~ bold('***')", parse = T) +
  annotate("text", x = -Inf, y = 0.9, hjust = 0, vjust = 1, size = 3.5,
           label = "~DBH[2017] * ':Plot'[TFE] ~ 'ns'", parse = T)

dev.off()

# Figure S6

# Sapling basal area
plot1 = ggplot(data = data_subp, aes(x = sapli_ba_m2_m2 * 10000, y = rgr_mean, colour = plot, fill = plot)) + 
  # geom_smooth(method = "lm", alpha = .5) +
  geom_point(alpha = .5, size = 3) +
  labs(x = expression("Small tree basal area" ~ (m^2 ~ ha^-1)),
       y = expression("Stem increment" ~ (cm ~ year^-1)),
       colour = element_blank(),
       fill = element_blank()) +
  scale_colour_viridis_d(option = "turbo", begin = .2, end = .8, labels = c("Control", "TFE")) +
  scale_fill_viridis_d(option = "turbo", begin = .2, end = .8, labels = c("Control", "TFE")) +
  annotate("text", x = Inf, y = .43, hjust = 1, vjust = 0, size = 3.5,
           label = "Basal~area~'ns'", parse = T) +
  annotate("text", x = Inf, y = .405, hjust = 1, vjust = 0, size = 3.5,
           label = "bold('Plot'[TFE])~bold('***')", parse = T) +
  annotate("text", x = Inf, y = .38, hjust = 1, vjust = 0, size = 3.5,
           label = "Basal~area:Plot[TFE]~'ns'", parse = T) +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12, colour = "black"))

# Tree basal area
plot2 = ggplot(data = data_subp, aes(x = tree_ba_m2_m2 * 10000, y = rgr_mean, colour = plot, fill = plot)) + 
  # geom_smooth(method = "lm", alpha = .5) +
  geom_point(alpha = .5, size = 3) +
  labs(x = expression("Large tree basal area" ~ (m^2 ~ ha^-1)),
       y = expression("Stem increment" ~ (cm ~ year^-1)),
       colour = element_blank(),
       fill = element_blank()) +
  scale_colour_viridis_d(option = "turbo", begin = .2, end = .8, labels = c("Control", "TFE")) +
  scale_fill_viridis_d(option = "turbo", begin = .2, end = .8, labels = c("Control", "TFE")) +
  annotate("text", x = Inf, y = .43, hjust = 1, vjust = 0, size = 3.5,
           label = "Basal~area~'ns'", parse = T) +
  annotate("text", x = Inf, y = .405, hjust = 1, vjust = 0, size = 3.5,
           label = "bold('Plot'[TFE])~bold('***')", parse = T) +
  annotate("text", x = Inf, y = .38, hjust = 1, vjust = 0, size = 3.5,
           label = "Basal~area:Plot[TFE]~'ns'", parse = T) +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12, colour = "black"))

# Export it
jpeg("Fig S6.jpg", width = 22.5, height = 10, units = "cm", res = 1000)

ggarrange(plot1, plot2, 
          ncol = 2, nrow = 1,
          labels = "auto",
          common.legend = T, legend = "right")

dev.off()

rm(plot1, plot2)

