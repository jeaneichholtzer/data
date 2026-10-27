# TITLE: Leaf and stem extracts of Mexican sunflower differentially affect growth, defense and herbivory in domesticated and wild lima bean
# AUTHORS: Jean Eichholtzer, Roberto Rafael Ruiz-Santiago, Alejandra del Socorro González-Moreno, Esaú Ruiz-Sánchez, Germán Giácoman Vallejos, Horacio Salomón Ballina-Gómez*
# CORREPONDING: horacio.bg@conkal.tecnm.mx

# 0 Set Working Directory

#setwd("yourpath/foldername")

### I Load packages ############################################################################################################################################################################################################

load_packages <- function(pkgs) 
  {missing_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
  if (length(missing_pkgs) > 0) 
    {install.packages(missing_pkgs, dependencies = TRUE)}
invisible(lapply(pkgs, library, character.only = TRUE))}

load_packages(c(
  "bestNormalize",
  "car",
  "DHARMa",
  "emmeans",
  "ggtext",
  "glmmTMB",
  "multcomp",
  "patchwork",
  "performance",
  "plotrix",
  "readxl",
  "scales",
  "tidyverse"))

### II Load data ############################################################################################################################################################################################################

defenses <- read_excel(path = "Measurements.xlsx", sheet = "Defense_traits") %>%
  mutate(Sampling_time = as.factor(Sampling_time),
         Accession = as.factor(Accession),
         Accession_id = as.factor(Accession_id),
         Dom_status = as.factor(Domestication_status),
         Row = as.factor(Row),
         Leaf_id = as.factor(Leaf_id),
         Treatment = factor(Treatment, levels = c("Control", "Leaf extract", "Stem extract")))

growth <- read_excel(path = "Measurements.xlsx", sheet = "Growth_traits") %>%
  mutate(Sampling_time = as.factor(Sampling_time),
         Date = as.Date(Date, format = "%m/%d/%Y"),
         Accession = as.factor(Accession),
         Accession_id = as.factor(str_sub(Accession_id)),
         Dom_status = as.factor(Domestication_status),
         Row = as.factor(Row),
         Treatment = factor(Treatment, levels = c("Control", "Leaf extract", "Stem extract")))

herbivory <- read_excel(path = "Measurements.xlsx", sheet = "Herbivory") %>%
  mutate(Sampling_time = as.factor(Sampling_time),
         Date = as.Date(Date, format = "%m/%d/%Y"),
         Accession = as.factor(Accession),
         Accession_id = as.factor(Accession_id),
         Dom_status = as.factor(Domestication_status),
         Row = as.factor(Row),
         Treatment = factor(Treatment, levels = c("Control", "Leaf extract", "Stem extract")),
         Count = Damaged_obs + Intact_obs)

trichomes <- read_excel(path = "Measurements.xlsx", sheet = "Trichomes") %>%
  mutate(Sampling_time = as.factor(Sampling_time),
         Accession = as.factor(Accession),
         Accession_id = as.factor(Accession_id),
         Dom_status = as.factor(Domestication_status),
         Row = as.factor(Row),
         Leaf_id = as.factor(Leaf_id),
         Treatment = factor(Treatment, levels = c("Control", "Leaf extract", "Stem extract")))

area_biomass <- read_excel(path = "Measurements.xlsx", sheet = "Leaf_area_biomass") %>%
  mutate(Sampling_time = as.factor(Sampling_time),
         Date = as.Date(Date, format = "%m/%d/%Y"),
         Accession = as.factor(Accession),
         Accession_id = as.factor(Accession_id),
         Dom_status = as.factor(Domestication_status),
         Row = as.factor(Row),
         Treatment = factor(Treatment, levels = c("Control", "Leaf extract", "Stem extract")),
         Leaf_area_sample = as.numeric(Leaf_area_sample),
         Leaf_biomass_sample = as.numeric(Leaf_biomass_sample))

nectaries <- read_excel(path = "Measurements.xlsx", sheet = "Nectaries") %>%
  subset(Accession_id %in% c("D1","W1","D5")) %>%
  mutate(Accession = as.factor(Accession),
         Accession_id = as.factor(Accession_id),
         Dom_status = as.factor(Domestication_status),
         Row = as.factor(Row),
         Leaf_id = as.factor(Leaf_id),
         Leaflet_id = as.factor(Leaflet_id),
         Treatment = factor(Treatment, levels = c("Control", "Leaf extract", "Stem extract")))

yield <- read_excel(path = "Measurements.xlsx", sheet = "Yield") %>%
  subset(Accession_id %in% c("D1","D5")) %>%
  mutate(Accession = as.factor(Accession),
         Accession_id = as.factor(Accession_id),
         Dom_status = as.factor(Domestication_status),
         Row = as.factor(Row),
         Treatment = factor(Treatment, levels = c("Control", "Leaf extract", "Stem extract")))

### III Themes ggplot2 ############################################################################################################################################################################################################

theme_1 <- theme(axis.line = element_line(),
                 axis.text.x = element_text(size= 14),
                 axis.text.y = element_text(size= 14),
                 axis.ticks.length = unit(0.2, "cm"),
                 axis.title.x = element_markdown(size= 16, margin = margin(t = 10)),
                 axis.title.y = element_text(size= 16, margin = margin(r = 10)),
                 legend.position = "none",
                 panel.grid = element_blank(),
                 panel.background = element_blank(),
                 panel.grid.major.y = element_line(colour = "grey85", linewidth = 0.3),
                 panel.grid.minor.y = element_line(colour = "grey93", linewidth = 0.2),
                 plot.title = element_text(size= 16, face="bold"),
                 strip.text.x = element_text(hjust = 0.05, size = 16),
                 strip.background = element_blank())

theme_legend_left <- theme(legend.position = c(0.02, 0.98),
                           legend.direction = "horizontal",
                           legend.justification = "left",
                           legend.box.just = "left",
                           legend.text = element_text(size = 12),
                           legend.title = element_blank(),
                           legend.key.height = unit(0.3, "cm"),
                           legend.key.width = unit(0.6, "cm"),
                           legend.background = element_rect(fill = alpha("white", 0.8), colour = NA))

theme_legend <- theme(legend.position = c(0.98, 0.98),
                      legend.direction = "horizontal",
                      legend.justification = "right",
                      legend.box.just = "right",
                      legend.text = element_text(size = 12),
                      legend.title = element_blank(),
                      legend.key.height = unit(0.3, "cm"),
                      legend.key.width = unit(0.6, "cm"),
                      legend.background = element_rect(fill = alpha("white", 0.8), colour = NA))

theme_legend_vertical <- theme(legend.position = c(0.98, 0.85),
                               legend.direction = "vertical",
                               legend.justification = "right",
                               legend.box.just = "right",
                               legend.text = element_text(size = 12),
                               legend.title = element_blank(),
                               legend.key.height = unit(0.3, "cm"),
                               legend.key.width = unit(0.6, "cm"),
                               legend.key.spacing.y = unit(0.1, "cm"),
                               legend.background = element_rect(fill = alpha("white", 0.8), colour = NA))

theme_legend_vertical_left <- theme(legend.position = c(0.02, 0.90),
                                    legend.direction = "vertical",
                                    legend.justification = "left",
                                    legend.box.just = "left",
                                    legend.text = element_text(size = 12),
                                    legend.title = element_blank(),
                                    legend.key.height = unit(0.3, "cm"),
                                    legend.key.width = unit(0.6, "cm"),
                                    legend.key.spacing.y = unit(0.1, "cm"),
                                    legend.background = element_rect(fill = alpha("white", 0.8), colour = NA))

trat_colors <- c("Control" = "white",
                 "Leaf extract" = "#4eb934",
                 "Stem extract" = "#e97034")

trat_labels_herbivory <- c("Control"  = "Control",
                           "Leaf extract" = "Leaf",
                           "Stem extract" = "Stem")

Sampling_time_labels = c("1" = "35",
                         "2" = "55",
                         "3" = "75")

### IV Calculations ############################################################################################################################################################################################################

###### 1 Leaf_area_sample#####

area_sample_raw <- trichomes %>%
  group_by(Sampling_time, Accession_id, Row, Treatment) %>%
  summarise(Leaf_area_sample = sum(Trifoliate_leaf_area, na.rm = TRUE))

###### 2 SLA #####

# Initial LA and LM means at Sampling time = 0

means_sampling_0 <- area_biomass %>%
  filter(Sampling_time == "0") %>%
  group_by(Accession_id) %>%
  summarise(mean_prev_area = mean(Leaf_area_sample, na.rm = TRUE),
            mean_prev_biomass = mean(Leaf_biomass_sample, na.rm = TRUE),
            .groups = "drop")

# LA and LM means at Sampling time = 1 and 2

means_previous_sampling <- area_biomass %>%
  filter(Sampling_time %in% c("1", "2")) %>%
  group_by(Accession_id, Treatment, Sampling_time) %>%
  summarise(mean_prev_area = mean(Leaf_area_sample, na.rm = TRUE),
            mean_prev_biomass = mean(Leaf_biomass_sample, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(Sampling_time = as.character(as.numeric(as.character(Sampling_time)) + 1))

# SLA and LMA calculations for Sampling_time = 1

sla_s1 <- area_biomass %>%
  filter(Sampling_time == "1") %>%
  left_join(means_sampling_0,
            by = "Accession_id") %>%
  mutate(SLA = ((mean_prev_area / mean_prev_biomass) + (Leaf_area_sample / Leaf_biomass_sample)) / 2)

lma_s1 <- area_biomass %>%
  filter(Sampling_time == "1") %>%
  left_join(means_sampling_0,
            by = "Accession_id") %>%
  mutate(LMA_2 = ((mean_prev_biomass / mean_prev_area) + (Leaf_biomass_sample / Leaf_area_sample)) / 2)

# SLA and LMA calculations for Sampling_time = 2 and 3

sla_s23 <- area_biomass %>%
  filter(Sampling_time %in% c("2", "3")) %>%
  left_join(means_previous_sampling,
            by = c("Accession_id", "Treatment", "Sampling_time")) %>%
  mutate(SLA = ((mean_prev_area / mean_prev_biomass) +  (Leaf_area_sample / Leaf_biomass_sample)) / 2)

lma_s23 <- area_biomass %>%
  filter(Sampling_time %in% c("2", "3")) %>%
  left_join(means_previous_sampling,
            by = c("Accession_id", "Treatment", "Sampling_time")) %>%
  mutate(LMA_2 = ((mean_prev_biomass / mean_prev_area) +  (Leaf_biomass_sample / Leaf_area_sample)) / 2)

# Merge data

sla <- bind_rows(area_biomass %>% filter(Sampling_time == "0"),
                 full_join(sla_s1, lma_s1),
                 full_join(sla_s23, lma_s23)) %>%
  arrange(Accession_id, Row, Treatment, Sampling_time)

###### 3 RGR #####

# Means at Sampling time = 0

means_sampling_0 <- growth %>%
  filter(Sampling_time == "0") %>%
  group_by(Accession_id) %>%
  summarise(mean_log_prev_diameter = mean(log(Diameter), na.rm = TRUE),
            mean_log_prev_no_leaves = mean(log(No_leaves), na.rm = TRUE),
            mean_prev_no_leaves = mean(No_leaves, na.rm = TRUE),
            mean_log_prev_length   = mean(log(Length), na.rm = TRUE),
            date_prev = first(Date),
            .groups = "drop")

# Means at Sampling time = 1 and 2

means_previous_sampling <- growth %>%
  filter(Sampling_time %in% c("1", "2")) %>%
  group_by(Accession_id, Treatment, Sampling_time) %>%
  summarise(mean_log_prev_diameter = mean(log(Diameter), na.rm = TRUE),
            mean_log_prev_no_leaves = mean(log(No_leaves), na.rm = TRUE),
            mean_prev_no_leaves = mean(No_leaves, na.rm = TRUE),
            mean_log_prev_length   = mean(log(Length), na.rm = TRUE),
            date_prev = first(Date),
            .groups = "drop") %>%
  mutate(Sampling_time = as.character(as.numeric(as.character(Sampling_time)) + 1))

# RGR and LPR calculations for Sampling_time = 1

rgr_s1 <- growth %>%
  filter(Sampling_time == "1") %>%
  left_join(means_sampling_0,
            by = "Accession_id") %>%
  mutate(days_diff = as.numeric(Date - date_prev),
         RGR_Diameter = (log(Diameter) - mean_log_prev_diameter) / days_diff,
         RGR_No_leaves = (log(No_leaves) - mean_log_prev_no_leaves) / days_diff,
         LPR = (No_leaves - mean_prev_no_leaves) / days_diff,
         RGR_Length = (log(Length) - mean_log_prev_length) / days_diff)

# RGR and LPR calculations for Sampling_time 2 and 3

rgr_s23 <- growth %>%
  filter(Sampling_time %in% c("2", "3")) %>%
  left_join(means_previous_sampling,
            by = c("Accession_id", "Treatment", "Sampling_time")) %>%
  mutate(days_diff = as.numeric(Date - date_prev),
         RGR_Diameter = (log(Diameter) - mean_log_prev_diameter) / days_diff,
         RGR_No_leaves = (log(No_leaves) - mean_log_prev_no_leaves) / days_diff,
         LPR = (No_leaves - mean_prev_no_leaves) / days_diff,
         RGR_Length = (log(Length) - mean_log_prev_length) / days_diff)

# Merge data and apply filter

growth_rgr <- bind_rows(growth %>% filter(Sampling_time == "0"), rgr_s1, rgr_s23) %>%
  filter(Sampling_time %in% c("1","2","3")) %>%
  mutate(Sampling_time = as.factor(Sampling_time)) %>%
  arrange(Accession_id, Row, Treatment, Sampling_time)

### V Statistical analyses and plots ############################################################################################################################################################################################################
###### 1 Diameter #####

model_diameter <- glmmTMB(Diameter ~ Sampling_time * Dom_status * Treatment + (1|Accession_id) + (1|Row),
                          family = gaussian(),
                          data = growth_rgr,
                          dispformula = ~ Sampling_time * Dom_status)

part1_model_diameter <- glmmTMB(Diameter ~ Sampling_time * Dom_status * Treatment + (1|Row),
                                family = gaussian(),
                                data = growth_rgr,
                                dispformula = ~ Sampling_time * Dom_status)

part2_model_diameter <- glmmTMB(Diameter ~ Sampling_time * Dom_status * Treatment + (1|Accession_id),
                                family = gaussian(),
                                data = growth_rgr,
                                dispformula = ~ Sampling_time * Dom_status)

plot(simulateResiduals(model_diameter))
check_model(model_diameter)

summary(model_diameter)
Anova(model_diameter, type="III")
anova(part1_model_diameter, model_diameter)
anova(part2_model_diameter, model_diameter)

ls1 <- emmeans(model_diameter, ~ Treatment | Sampling_time | Dom_status, type="response")
cld1 <- as.data.frame(cld(ls1, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld1$.group <- str_trim(cld1$.group)
cld1$Treatment <- factor(cld1$Treatment,levels = c("Control", "Leaf extract", "Stem extract"))
cld1 <- cld1 %>% arrange(Dom_status, Sampling_time, Treatment)

df_raw1 <- growth_rgr %>%
  group_by(Sampling_time, Treatment, Dom_status) %>%
  summarise(mean = mean(Diameter, na.rm = TRUE),
            se = std.error(Diameter, na.rm = TRUE),
            .groups = "drop") %>% 
  left_join(cld1 %>% dplyr::select(Sampling_time, Treatment, Dom_status, .group), by = c("Sampling_time", "Treatment", "Dom_status"))

plot_diameter <- ggplot(df_raw1, aes(x = Sampling_time, y = mean, fill = Treatment)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se,
                    ymax = mean+se),
                position = position_dodge(width = 0.8),
                width = 0.2,
                linewidth = 0.5) +
  facet_wrap(~ Dom_status) +
  geom_text(aes(label = .group, y = mean+se + 0.05 * max(mean)), position = position_dodge(width = 0.8), size = 5) +
  labs(x = "Time (days after sowing)",
       y = "Stem diameter (mm)",
       fill = "Treatment") +
  scale_x_discrete(labels = Sampling_time_labels) +
  scale_y_continuous(limits=c(0,15)) +
  scale_fill_manual(values = trat_colors) +
  theme_1

plot_diameter

###### 2 Hardness #####

set.seed(123)
BNobject <- bestNormalize(defenses$Hardness)
defenses$Hardness_tr <- predict(BNobject)

model_hardness <- glmmTMB(Hardness_tr ~ Sampling_time * Dom_status * Treatment + (1|Accession_id) + (1|Row),
                  family = t_family(),
                  data = defenses,
                  dispformula = ~ Sampling_time)

part1_model_hardness <- glmmTMB(Hardness_tr ~ Sampling_time * Dom_status * Treatment + (1|Row),
                        family = t_family(),
                        data = defenses,
                        dispformula = ~ Sampling_time)

part2_model_hardness <- glmmTMB(Hardness_tr ~ Sampling_time * Dom_status * Treatment + (1|Accession_id),
                        family = t_family(),
                        data = defenses,
                        dispformula = ~ Sampling_time)

plot(simulateResiduals(model_hardness))
check_model(model_hardness)

summary(model_hardness)
Anova(model_hardness, type="III")
anova(part1_model_hardness, model_hardness)
anova(part2_model_hardness, model_hardness)

ls1 <- emmeans(model_hardness, ~ Treatment | Sampling_time | Dom_status, type="response")
cld1 <- as.data.frame(cld(ls1, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld1$.group <- str_trim(cld1$.group)
cld1$Treatment <- factor(cld1$Treatment,levels = c("Control", "Leaf extract", "Stem extract"))
cld1 <- cld1 %>% arrange(Dom_status, Sampling_time, Treatment)
cld1

df_raw1 <- defenses %>%
  group_by(Sampling_time, Treatment, Dom_status) %>%
  summarise(mean = mean(Hardness, na.rm = TRUE),
            se = std.error(Hardness, na.rm = TRUE),
            .groups = "drop") %>% 
  left_join(cld1 %>% dplyr::select(Sampling_time, Treatment, Dom_status, .group), by = c("Sampling_time", "Treatment", "Dom_status"))

plot_hardness <- ggplot(df_raw1, aes(x = Sampling_time, y = mean, fill = Treatment)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se,
                    ymax = mean+se),
                position = position_dodge(width = 0.8),
                width = 0.2,
                linewidth = 0.5) +
  facet_wrap(~ Dom_status) +
  geom_text(aes(label = .group, y = mean+se + 0.05 * max(mean)), position = position_dodge(width = 0.8), size = 5) +
  labs(x = "Time (days after sowing)",
       y = expression(Hardness~(kg~cm^{-2})),
       fill = "Treatment") +
  scale_x_discrete(labels = Sampling_time_labels) +
  scale_fill_manual(values = trat_colors) +
  theme_1

plot_hardness

###### 3 Herbivory #####

model_herbivory <- glmmTMB(cbind(Damaged_obs, Intact_obs) ~ Sampling_time * Dom_status * Treatment + (1|Accession_id) + (1|Row),
                           family=binomial(link="logit"),
                           data = herbivory)

part1_model_herbivory <- glmmTMB(cbind(Damaged_obs, Intact_obs) ~ Sampling_time * Dom_status * Treatment + (1|Row),
                                 family=binomial(link="logit"),
                                 data = herbivory)

part2_model_herbivory <- glmmTMB(cbind(Damaged_obs, Intact_obs) ~ Sampling_time * Dom_status * Treatment + (1|Accession_id),
                                 family=binomial(link="logit"),
                                 data = herbivory)

plot(simulateResiduals(model_herbivory))
check_model(model_herbivory)

summary(model_herbivory)
Anova(model_herbivory, type="II")
anova(part1_model_herbivory, model_herbivory)
anova(part2_model_herbivory, model_herbivory)

ls1 <- emmeans(model_herbivory, ~ Sampling_time, type="response")
cld1 <- as.data.frame(cld(ls1, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld1$.group <- str_trim(cld1$.group)
cld1 <- cld1 %>% arrange(Sampling_time)

df_raw1 <- herbivory %>%
  group_by(Sampling_time) %>%
  summarise(mean = mean(Herbivory_incidence, na.rm = TRUE),
            se = std.error(Herbivory_incidence, na.rm = TRUE),
            .groups = "drop") %>%
  left_join(cld1 %>% dplyr::select(Sampling_time, .group), by = c("Sampling_time"))

plot_herbivory_1 <- ggplot(df_raw1, aes(x = Sampling_time, y = mean)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6, colour = "black", fill="grey95", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se,
                    ymax = mean+se),
                position = position_dodge(width = 0.8),
                width = 0.2,
                linewidth = 0.5) +
  geom_text(aes(label = .group, y = mean+se + 0.05 * max(mean)), position = position_dodge(width = 0.8), size = 5) +
  labs(x = "Time (days after sowing)",
       y = "Herbivory (%)") +
  scale_x_discrete(labels = Sampling_time_labels) +
  scale_y_continuous(labels = scales::percent, limits=c(0,0.80)) +
  theme_1

plot_herbivory_1

ls2 <- emmeans(model_herbivory, ~ Dom_status, type="response")
cld2 <- as.data.frame(cld(ls2, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld2$.group <- str_trim(cld2$.group)
cld2 <- cld2 %>% arrange(Dom_status)

df_raw2 <- herbivory %>%
  group_by(Dom_status) %>%
  summarise(mean = mean(Herbivory_incidence, na.rm = TRUE),
            se = std.error(Herbivory_incidence, na.rm = TRUE),
            .groups = "drop") %>% 
  left_join(cld2 %>% dplyr::select(Dom_status, .group), by = c("Dom_status"))

plot_herbivory_2 <- ggplot(df_raw2, aes(x = Dom_status, y = mean)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6, colour = "black", fill="grey95", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se,
                    ymax = mean+se),
                position = position_dodge(width = 0.8),
                width = 0.2,
                linewidth = 0.5) +
  geom_text(aes(label = .group, y = mean+se + 0.06 * max(mean)), position = position_dodge(width = 0.8), size = 5) +
  labs(x = "Status",
       y = "Herbivory (%)") +
  scale_y_continuous(labels = scales::percent, limits=c(0,0.80)) +
  theme_1

plot_herbivory_2

ls3 <- emmeans(model_herbivory, ~ Treatment, type="response")
cld3 <- as.data.frame(cld(ls3, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld3$.group <- str_trim(cld3$.group)
cld3$Treatment <- factor(cld3$Treatment,levels = c("Control", "Leaf extract", "Stem extract"))
cld3 <- cld3 %>% arrange(Treatment)

df_raw3 <- herbivory %>%
  group_by(Treatment) %>%
  summarise(mean = mean(Herbivory_incidence, na.rm = TRUE),
            se = std.error(Herbivory_incidence, na.rm = TRUE),
            .groups = "drop") %>%
  left_join(cld3 %>% dplyr::select(Treatment, .group), by = c("Treatment"))

plot_herbivory_3 <- ggplot(df_raw3, aes(x = Treatment, y = mean, fill = Treatment)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se,
                    ymax = mean+se),
                position = position_dodge(width = 0.8),
                width = 0.2,
                linewidth = 0.5) +
  geom_text(aes(label = .group, y = mean+se + 0.06 * max(mean)), position = position_dodge(width = 0.8), size = 5) +
  labs(x = "Treatment",
       y = "Herbivory (%)",
       fill = "Treatment") +
  scale_x_discrete(labels = trat_labels_herbivory) +
  scale_fill_manual(values = trat_colors) +
  scale_y_continuous(labels = scales::percent, limits=c(0,0.80)) +
  theme_1

plot_herbivory_3

ls4 <- emmeans(model_herbivory, ~ Treatment | Sampling_time | Dom_status, type="response")
cld4 <- as.data.frame(cld(ls4, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld4$.group <- str_trim(cld4$.group)
cld4$Treatment <- factor(cld4$Treatment,levels = c("Control", "Leaf extract", "Stem extract"))
cld4 <- cld4 %>% arrange(Dom_status, Sampling_time, Treatment)
cld4
pairs(ls4)

df_raw4 <- herbivory %>%
  group_by(Sampling_time, Treatment, Dom_status) %>%
  summarise(mean = mean(Herbivory_incidence, na.rm = TRUE),
            se = std.error(Herbivory_incidence, na.rm = TRUE),
            .groups = "drop") %>% 
  left_join(cld4 %>% dplyr::select(Sampling_time, Treatment, Dom_status, .group), by = c("Sampling_time", "Treatment", "Dom_status"))

plot_herbivory_4 <- ggplot(df_raw4, aes(x = Sampling_time, y = mean, fill = Treatment)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se,
                    ymax = mean+se),
                position = position_dodge(width = 0.8),
                width = 0.2,
                linewidth = 0.5) +
  facet_wrap(~ Dom_status) +
  geom_text(aes(label = .group, y = mean+se + 0.05 * max(mean)), position = position_dodge(width = 0.8), size = 5) +
  labs(x = "Treatment",
       y = "Herbivory (%)",
       fill = "Treatment") +
  scale_x_discrete(labels = Sampling_time_labels) +
  scale_fill_manual(values = trat_colors, labels = trat_labels_herbivory) +
  theme_1

plot_herbivory_4

###### 4 Length #####

set.seed(123)
BNobject <- bestNormalize(growth_rgr$Length)
growth_rgr$Length_tr <- predict(BNobject)

model_length <- glmmTMB(Length_tr ~ Sampling_time * Dom_status * Treatment + (1|Accession_id) + (1|Row),
                        family = gaussian(),
                        dispformula = ~ Sampling_time,
                        data = growth_rgr)

part1_model_length <- glmmTMB(Length_tr ~ Sampling_time * Dom_status * Treatment + (1|Row),
                              family = gaussian(),
                              dispformula = ~ Sampling_time,
                              data = growth_rgr)

part2_model_length <- glmmTMB(Length_tr ~ Sampling_time * Dom_status * Treatment + (1|Accession_id),
                              family = gaussian(),
                              dispformula = ~ Sampling_time,
                              data = growth_rgr)

plot(simulateResiduals(model_length))
check_model(model_length)

summary(model_length)
Anova(model_length, type="III")
anova(part1_model_length, model_length)
anova(part2_model_length, model_length)

ls1 <- emmeans(model_length, ~ Treatment | Sampling_time, type="response")
cld1 <- as.data.frame(cld(ls1, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld1$.group <- str_trim(cld1$.group)
cld1$Treatment <- factor(cld1$Treatment,levels = c("Control", "Leaf extract", "Stem extract"))
cld1 <- cld1 %>% arrange(Sampling_time, Treatment)

df_raw1 <- growth_rgr %>%
  group_by(Sampling_time, Treatment) %>%
  summarise(mean = mean(Length, na.rm = TRUE),
            se = std.error(Length, na.rm = TRUE),
            .groups = "drop") %>% 
  left_join(cld1 %>% dplyr::select(Sampling_time, Treatment, .group), by = c("Sampling_time", "Treatment"))

plot_length <- ggplot(df_raw1, aes(x = Sampling_time, y = mean, fill = Treatment)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se,
                    ymax = mean+se),
                position = position_dodge(width = 0.8),
                width = 0.2,
                linewidth = 0.5) +
  geom_text(aes(label = .group, y = mean+se + 0.05 * max(mean)), position = position_dodge(width = 0.8), size = 5) +
  labs(x = "Time (days after sowing)",
       y = "Length (cm)",
       fill = "Treatment") +
  scale_x_discrete(labels = Sampling_time_labels) +
  scale_fill_manual(values = trat_colors) +
  theme_1

plot_length

###### 5 Nectaries #####

model_nectaries <- glmmTMB(Degrees_brix ~ Dom_status * Treatment + (1|Row/Leaf_id),
                           family = gaussian(),
                           data = nectaries)

part1_model_nectaries <- glmmTMB(Degrees_brix ~ Dom_status * Treatment,
                                 family = gaussian(),
                                 data = nectaries)

plot(simulateResiduals(model_nectaries))
check_model(model_nectaries)

summary(model_nectaries)
Anova(model_nectaries, type = "III")
anova(part1_model_nectaries, model_nectaries)

ls1 <- emmeans(model_nectaries, ~ Treatment | Dom_status, type="response")
cld1 <- as.data.frame(cld(ls1, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld1$.group <- str_trim(cld1$.group)
cld1$Treatment <- factor(cld1$Treatment,levels = c("Control", "Leaf extract", "Stem extract"))
cld1 <- cld1 %>% arrange(Dom_status, Treatment)

df_raw1 <- nectaries %>%
  group_by(Treatment, Dom_status) %>%
  summarise(mean = mean(Degrees_brix, na.rm = TRUE),
            se = std.error(Degrees_brix, na.rm = TRUE),
            .groups = "drop") %>% 
  left_join(cld1 %>% dplyr::select(Treatment, Dom_status, .group), by = c("Treatment", "Dom_status"))

plot_nectaries <- ggplot(df_raw1, aes(x=Dom_status, y=mean, fill=Treatment)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.5, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position = position_dodge(width = 0.6),
                width = 0.2) +
  geom_text(aes(label = .group, y = mean+se + 0.15 * max(mean)), position = position_dodge(width = 0.6), size = 5) +
  labs(x = "Status",
       y = "EFN (°Bx)",
       fill = "Treatment") +
  scale_fill_manual(values = trat_colors) +
  theme_1

plot_nectaries

###### 6 Number of leaves #####

set.seed(123)
BNobject <- bestNormalize(growth_rgr$No_leaves)
growth_rgr$No_leaves_tr <- predict(BNobject)

model_no_leaves <- glmmTMB(No_leaves_tr ~ Sampling_time * Dom_status * Treatment + (1|Accession_id) + (1|Row),
                           family = gaussian(),
                           data = growth_rgr,
                           dispformula = ~ Sampling_time * Dom_status * Treatment)

part1_model_no_leaves <- glmmTMB(No_leaves_tr ~ Sampling_time * Dom_status * Treatment + (1|Row),
                                 family = gaussian(),
                                 data = growth_rgr,
                                 dispformula = ~ Sampling_time * Dom_status * Treatment)

part2_model_no_leaves <- glmmTMB(No_leaves_tr ~ Sampling_time * Dom_status * Treatment + (1|Accession_id),
                                 family = gaussian(),
                                 data = growth_rgr,
                                 dispformula = ~ Sampling_time * Dom_status * Treatment)

plot(simulateResiduals(model_no_leaves))
check_model(model_no_leaves)

summary(model_no_leaves)
Anova(model_no_leaves, type="III")
anova(part1_model_no_leaves, model_no_leaves)
anova(part2_model_no_leaves, model_no_leaves)

plot(simulateResiduals(model_no_leaves))

ls1 <- emmeans(model_no_leaves, ~ Treatment | Dom_status, type="response")
cld1 <- as.data.frame(cld(ls1, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld1$.group <- str_trim(cld1$.group)
cld1$Treatment <- factor(cld1$Treatment,levels = c("Control", "Leaf extract", "Stem extract"))
cld1 <- cld1 %>% arrange(Dom_status, Treatment)

df_raw1 <- growth_rgr %>%
  group_by(Treatment, Dom_status) %>%
  summarise(mean = mean(No_leaves, na.rm = TRUE),
            se = std.error(No_leaves, na.rm = TRUE),
            .groups = "drop") %>% 
  left_join(cld1 %>% dplyr::select(Treatment, Dom_status, .group), by = c("Treatment", "Dom_status"))

plot_no_leaves <- ggplot(df_raw1, aes(x = Dom_status, y = mean, fill = Treatment)) +
  geom_col(position = position_dodge(width = 0.5), width = 0.3, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se,
                    ymax = mean+se),
                position = position_dodge(width = 0.5),
                width = 0.12,
                linewidth = 0.5) +
  geom_text(aes(label = .group, y = mean+se + 0.05 * max(mean)), position = position_dodge(width = 0.5), size = 5) +
  labs(x = "Status",
       y = "Number of leaves",
       fill = "Treatment") +
  scale_x_discrete(labels = Sampling_time_labels) +
  scale_fill_manual(values = trat_colors) +
  theme_1

plot_no_leaves

###### 7 RGR Diameter #####

set.seed(123)
BNobject <- bestNormalize(growth_rgr$RGR_Diameter)
growth_rgr$RGR_Diameter_tr <- predict(BNobject)

model_rgr_diameter <- glmmTMB(RGR_Diameter_tr ~ Sampling_time * Dom_status * Treatment + (1|Accession_id) + (1|Row),
                              family = t_family(),
                              dispformula = ~ Sampling_time,
                              data = growth_rgr)

part1_model_rgr_diameter <- glmmTMB(RGR_Diameter_tr ~ Sampling_time * Dom_status * Treatment + (1|Row),
                                    family = t_family(),
                                    dispformula = ~ Sampling_time,
                                    data = growth_rgr)

part2_model_rgr_diameter <- glmmTMB(RGR_Diameter_tr ~ Sampling_time * Dom_status * Treatment + (1|Accession_id),
                                    family = t_family(),
                                    dispformula = ~ Sampling_time,
                                    data = growth_rgr)

plot(simulateResiduals(model_rgr_diameter))
check_model(model_rgr_diameter)

summary(model_rgr_diameter)
Anova(model_rgr_diameter, type="III")
anova(part1_model_rgr_diameter, model_rgr_diameter)
anova(part2_model_rgr_diameter, model_rgr_diameter)

ls1 <- emmeans(model_rgr_diameter, ~ Treatment | Sampling_time | Dom_status, type="response")
cld1 <- as.data.frame(cld(ls1, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld1$.group <- str_trim(cld1$.group)
cld1$Treatment <- factor(cld1$Treatment,levels = c("Control", "Leaf extract", "Stem extract"))
cld1 <- cld1 %>% arrange(Dom_status, Sampling_time, Treatment)

df_raw1 <- growth_rgr %>%
  group_by(Sampling_time, Treatment, Dom_status) %>%
  summarise(mean = mean(RGR_Diameter, na.rm = TRUE),
            se = std.error(RGR_Diameter, na.rm = TRUE),
            .groups = "drop") %>% 
  left_join(cld1 %>% dplyr::select(Sampling_time, Treatment, Dom_status, .group), by = c("Sampling_time", "Treatment", "Dom_status"))

plot_rgr_diameter <- ggplot(df_raw1, aes(x = Sampling_time, y = mean, fill = Treatment)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se,
                    ymax = mean+se),
                position = position_dodge(width = 0.8),
                width = 0.2,
                linewidth = 0.5) +
  facet_wrap(~ Dom_status) +
  geom_text(aes(label = .group, y = mean+se + 0.1 * max(mean)), position = position_dodge(width = 0.8), size = 5) +
  labs(x = "Time (days after sowing)",
       y = expression(RGR[italic(D)]~(mm~mm^{-1}~day^{-1})),
       fill = "Treatment") +
  scale_x_discrete(labels = Sampling_time_labels) +
  scale_fill_manual(values = trat_colors) +
  theme_1

plot_rgr_diameter

###### 8 RGR Length #####

set.seed(123)
BNobject <- bestNormalize(growth_rgr$RGR_Length)
growth_rgr$RGR_Length_tr <- predict(BNobject)

model_rgr_length <- glmmTMB(RGR_Length_tr ~ Sampling_time * Dom_status * Treatment + (1|Accession_id) + (1|Row),
                            family = t_family(),
                            data = growth_rgr,
                            dispformula = ~ Sampling_time * Treatment * Dom_status)

part1_model_rgr_length <- glmmTMB(RGR_Length_tr ~ Sampling_time * Dom_status * Treatment + (1|Row),
                                  family = t_family(),
                                  data = growth_rgr,
                                  dispformula = ~ Sampling_time * Treatment * Dom_status)

part2_model_rgr_length <- glmmTMB(RGR_Length_tr ~ Sampling_time * Dom_status * Treatment + (1|Accession_id),
                                  family = t_family(),
                                  data = growth_rgr,
                                  dispformula = ~ Sampling_time * Treatment * Dom_status)

plot(simulateResiduals(model_rgr_length))
check_model(model_rgr_length)

summary(model_rgr_length)
Anova(model_rgr_length, type="III")
anova(part1_model_rgr_length, model_rgr_length)
anova(part2_model_rgr_length, model_rgr_length)

ls1 <- emmeans(model_rgr_length, ~ Treatment | Sampling_time | Dom_status, type="response")
cld1 <- as.data.frame(cld(ls1, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld1$.group <- str_trim(cld1$.group)
cld1$Treatment <- factor(cld1$Treatment,levels = c("Control", "Leaf extract", "Stem extract"))
cld1 <- cld1 %>% arrange(Dom_status, Sampling_time, Treatment)

df_raw1 <- growth_rgr %>%
  group_by(Sampling_time, Treatment, Dom_status) %>%
  summarise(mean = mean(RGR_Length, na.rm = TRUE),
            se = std.error(RGR_Length, na.rm = TRUE),
            .groups = "drop") %>% 
  left_join(cld1 %>% dplyr::select(Sampling_time, Treatment, Dom_status, .group), by = c("Sampling_time", "Treatment", "Dom_status"))

plot_rgr_height <- ggplot(df_raw1, aes(x = Sampling_time, y = mean, fill = Treatment)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se,
                    ymax = mean+se),
                position = position_dodge(width = 0.8),
                width = 0.2,
                linewidth = 0.5) +
  facet_wrap(~ Dom_status) +
  geom_text(aes(label = .group, y = mean+se + 0.1 * max(mean)), position = position_dodge(width = 0.8), size = 5) +
  labs(x = "Time (days after sowing)",
       y = expression(RGR[italic(L)]~(cm~cm^{-1}~day^{-1})),
       fill = "Treatment") +
  scale_x_discrete(labels = Sampling_time_labels) +
  scale_fill_manual(values = trat_colors) +
  theme_1

plot_rgr_height

###### 9 RGR Number of leaves #####

set.seed(123)
BNobject <- bestNormalize(growth_rgr$RGR_No_leaves)
growth_rgr$RGR_No_leaves_tr <- predict(BNobject)

model_rgr_no_leaves <- glmmTMB(RGR_No_leaves_tr ~ Sampling_time * Dom_status * Treatment + (1|Accession_id) + (1|Row),
                               family = gaussian(),
                               dispformula = ~ Sampling_time * Dom_status * Treatment,
                               data = growth_rgr)

part1_model_rgr_no_leaves <- glmmTMB(RGR_No_leaves_tr ~ Sampling_time * Dom_status * Treatment + (1|Row),
                                     family = gaussian(),
                                     dispformula = ~ Sampling_time * Dom_status * Treatment,
                                     data = growth_rgr)

part2_model_rgr_no_leaves <- glmmTMB(RGR_No_leaves_tr ~ Sampling_time * Dom_status * Treatment + (1|Accession_id),
                                     family = gaussian(),
                                     dispformula = ~ Sampling_time * Dom_status * Treatment,
                                     data = growth_rgr)

plot(simulateResiduals(model_rgr_no_leaves))
check_model(model_rgr_no_leaves)

summary(model_rgr_no_leaves)
Anova(model_rgr_no_leaves, type="III")
anova(part1_model_rgr_no_leaves, model_rgr_no_leaves)
anova(part2_model_rgr_no_leaves, model_rgr_no_leaves)

ls1 <- emmeans(model_rgr_no_leaves, ~ Treatment | Sampling_time | Dom_status, type="response")
cld1 <- as.data.frame(cld(ls1, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld1$.group <- str_trim(cld1$.group)
cld1$Treatment <- factor(cld1$Treatment,levels = c("Control", "Leaf extract", "Stem extract"))
cld1 <- cld1 %>% arrange(Dom_status, Sampling_time, Treatment)

df_raw1 <- growth_rgr %>%
  group_by(Sampling_time, Treatment, Dom_status) %>%
  summarise(mean = mean(RGR_No_leaves, na.rm = TRUE),
            se = std.error(RGR_No_leaves, na.rm = TRUE),
            .groups = "drop") %>% 
  left_join(cld1 %>% dplyr::select(Sampling_time, Treatment, Dom_status, .group), by = c("Sampling_time", "Treatment", "Dom_status"))

plot_rgr_no_leaves <- ggplot(df_raw1, aes(x = Sampling_time, y = mean, fill = Treatment)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se,
                    ymax = mean+se),
                position = position_dodge(width = 0.8),
                width = 0.2,
                linewidth = 0.5) +
  facet_wrap(~ Dom_status) +
  geom_text(aes(label = .group, y = mean+se + 0.1 * max(mean)), position = position_dodge(width = 0.8), size = 5) +
  labs(x = "Time (days after sowing)",
       y = expression(RGR[italic(NL)]~(no.~no.^{-1}~day^{-1})),
       fill = "Treatment") +
  scale_x_discrete(labels = Sampling_time_labels) +
  scale_fill_manual(values = trat_colors) +
  theme_1

plot_rgr_no_leaves

###### 10 SLA #####

set.seed(123)

model_sla <- glmmTMB(SLA ~ Sampling_time * Dom_status * Treatment + (1|Accession_id) + (1|Row),
                     family = Gamma(link="log"),
                     data = sla)

part1_model_sla <- glmmTMB(SLA ~ Sampling_time * Dom_status * Treatment + (1|Row),
                           family = Gamma(link="log"),
                           data = sla)

part2_model_sla <- glmmTMB(SLA ~ Sampling_time * Dom_status * Treatment + (1|Accession_id),
                           family = Gamma(link="log"),
                           data = sla)

plot(simulateResiduals(model_sla))
check_model(model_sla)

summary(model_sla)
Anova(model_sla, type="II")
anova(part1_model_sla, model_sla)
anova(part2_model_sla, model_sla)

###### 11 Thickness #####

model_thickness <- glmmTMB(Thickness ~ Sampling_time * Dom_status * Treatment + (1|Accession_id) + (1|Row),
                           family =  gaussian(),
                           data = defenses,
                           dispformula = ~ Treatment + Dom_status)

part1_model_thickness <- glmmTMB(Thickness ~ Sampling_time * Dom_status * Treatment + (1|Row),
                                 family =  gaussian(),
                                 data = defenses,
                                 dispformula = ~ Treatment + Dom_status)

part2_model_thickness <- glmmTMB(Thickness ~ Sampling_time * Dom_status * Treatment + (1|Accession_id),
                                 family =  gaussian(),
                                 data = defenses,
                                 dispformula = ~ Treatment + Dom_status)

plot(simulateResiduals(model_thickness))
check_model(model_thickness)

summary(model_thickness)
Anova(model_thickness, type="III")
anova(part1_model_thickness, model_thickness)
anova(part2_model_thickness, model_thickness)

ls1 <- emmeans(model_thickness, ~ Sampling_time | Dom_status, type="response")
cld1 <- as.data.frame(cld(ls1, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld1$.group <- str_trim(cld1$.group)
cld1 <- cld1 %>% arrange(Dom_status, Sampling_time)
cld1

df_raw1 <- defenses %>%
  group_by(Sampling_time, Dom_status) %>%
  summarise(mean = mean(Thickness, na.rm = TRUE),
            se = std.error(Thickness, na.rm = TRUE),
            .groups = "drop") %>%
  left_join(cld1 %>% dplyr::select(Sampling_time, Dom_status, .group), by = c("Sampling_time", "Dom_status"))

plot_thickness <- ggplot(df_raw1, aes(x = Sampling_time, y = mean)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se,
                    ymax = mean+se),
                position = position_dodge(width = 0.8),
                width = 0.2,
                linewidth = 0.5) +
  facet_wrap(~ Dom_status) +
  geom_text(aes(label = .group, y = mean+se + 0.05 * max(mean)), position = position_dodge(width = 0.8), size = 5) +
  labs(x = "Time (days after sowing)",
       y = "Leaf thickness (mm)") +
  scale_x_discrete(labels = Sampling_time_labels) +
  theme_1

plot_thickness

###### 12 Trichome Density #####

model_trichomes <- glmmTMB(sqrt(Trichomes_cm2) ~ Sampling_time * Dom_status * Treatment + (1|Accession_id) + (1|Row),
                           family = gaussian(),
                           dispformula = ~ Dom_status,
                           data = trichomes)

part1_model_trichomes <- glmmTMB(sqrt(Trichomes_cm2) ~ Sampling_time * Dom_status * Treatment + (1|Row),
                                 family = gaussian(),
                                 dispformula = ~ Dom_status,
                                 data = trichomes)

part2_model_trichomes <- glmmTMB(sqrt(Trichomes_cm2) ~ Sampling_time * Dom_status * Treatment + (1|Accession_id),
                                 family = gaussian(),
                                 dispformula = ~ Dom_status,
                                 data = trichomes)

plot(simulateResiduals(model_trichomes))
check_model(model_trichomes)

summary(model_trichomes)
Anova(model_trichomes, type="III")
anova(part1_model_trichomes, model_trichomes)
anova(part2_model_trichomes, model_trichomes)

ls1 <- emmeans(model_trichomes, ~ Treatment | Dom_status, type="response")
cld1 <- as.data.frame(cld(ls1, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld1$.group <- str_trim(cld1$.group)
cld1$Treatment <- factor(cld1$Treatment,levels = c("Control", "Leaf extract", "Stem extract"))
cld1 <- cld1 %>% arrange(Dom_status, Treatment)

df_raw1 <- trichomes %>%
  group_by(Treatment, Dom_status) %>%
  summarise(mean = mean(Trichomes_cm2, na.rm = TRUE),
            se = std.error(Trichomes_cm2, na.rm = TRUE),
            .groups = "drop") %>% 
  left_join(cld1 %>% dplyr::select(Treatment, Dom_status, .group), by = c("Treatment", "Dom_status"))

plot_trichomes <- ggplot(df_raw1, aes(x = Dom_status, y = mean, fill = Treatment)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se,
                    ymax = mean+se),
                position = position_dodge(width = 0.8),
                width = 0.2,
                linewidth = 0.5) +
  geom_text(aes(label = .group, y = mean+se + 0.15 * max(mean)), position = position_dodge(width = 0.8), size = 5) +
  labs(x = "Status",
       y = expression(TD~(no.~cm^{-2})),
       fill = "Treatment") +
  scale_fill_manual(values = trat_colors) +
  theme_1

plot_trichomes

###### 13 Trifoliate LA #####

model_trifoliate_la <- glmmTMB(log(Trifoliate_leaf_area) ~ Sampling_time * Dom_status * Treatment + (1|Accession_id) + (1|Row),
                               family = gaussian(),
                               dispformula = ~ Dom_status,
                               data = trichomes)

part1_model_trifoliate_la <- glmmTMB(log(Trifoliate_leaf_area) ~ Sampling_time * Dom_status * Treatment + (1|Row),
                               family = gaussian(),
                               dispformula = ~ Dom_status,
                               data = trichomes)

part2_model_trifoliate_la <- glmmTMB(log(Trifoliate_leaf_area) ~ Sampling_time * Dom_status * Treatment + (1|Accession_id),
                               family = gaussian(),
                               dispformula = ~ Dom_status,
                               data = trichomes)

plot(simulateResiduals(model_trifoliate_la))
check_model(model_trifoliate_la)

summary(model_trifoliate_la)
Anova(model_trifoliate_la, type="II")
anova(part1_model_trifoliate_la, model_trifoliate_la)
anova(part2_model_trifoliate_la, model_trifoliate_la)

ls1 <- emmeans(model_trifoliate_la, ~ Treatment | Sampling_time, type="response")
cld1 <- as.data.frame(cld(ls1, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld1$.group <- str_trim(cld1$.group)
cld1$Treatment <- factor(cld1$Treatment,levels = c("Control", "Leaf extract", "Stem extract"))
cld1 <- cld1 %>% arrange(Sampling_time, Treatment)

df_raw1 <- trichomes %>%
  group_by(Sampling_time, Treatment) %>%
  summarise(mean = mean(Trifoliate_leaf_area, na.rm = TRUE),
            se = std.error(Trifoliate_leaf_area, na.rm = TRUE),
            .groups = "drop") %>% 
  left_join(cld1 %>% dplyr::select(Sampling_time, Treatment, .group), by = c("Sampling_time", "Treatment"))

plot_leaf_area <- ggplot(df_raw1, aes(x = Sampling_time, y = mean, fill = Treatment)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se,
                    ymax = mean+se),
                position = position_dodge(width = 0.8),
                width = 0.2,
                linewidth = 0.5) +
  geom_text(aes(label = .group, y = mean+se + 0.05 * max(mean)), position = position_dodge(width = 0.8), size = 5) +
  labs(x = "Time (days after sowing)",
       y = expression(LA~(cm^{2})),
       fill = "Treatment") +
  scale_x_discrete(labels = Sampling_time_labels) +
  scale_fill_manual(values = trat_colors) +
  theme_1

plot_leaf_area

###### 14 Yield #####

model_yield <- glmmTMB(Yield_per_plant ~ Treatment * Accession_id + (1|Row),
                       family= Gamma(link="log"),
                       data = yield)

part1_model_yield <- glmmTMB(Yield_per_plant ~ Treatment * Accession_id,
                        family= Gamma(link="log"),
                        data = yield)

plot(simulateResiduals(model_yield))
check_model(model_yield)

summary(model_yield)
Anova(model_yield, type ="III")
anova(part1_model_yield, model_yield)

ls1 <- emmeans(model_yield, ~ Accession_id, type="response")
cld1 <- as.data.frame(cld(ls1, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld1$.group <- str_trim(cld1$.group)
cld1 <- cld1 %>% arrange(Accession_id)

df_raw1 <- yield %>%
  group_by(Accession_id) %>%
  summarise(mean = mean(Yield_per_plant, na.rm = TRUE),
            se = std.error(Yield_per_plant, na.rm = TRUE),
            .groups = "drop") %>% 
  left_join(cld1 %>% dplyr::select(Accession_id, .group), by = c("Accession_id"))

plot_yield_1 <- ggplot(df_raw1, aes(x=Accession_id, y=mean)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.5, colour = "black", fill="gray95", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position = position_dodge(width = 0.6),
                width = 0.2) +
  geom_text(aes(label = .group, y = mean+se + 0.05 * max(mean)), position = position_dodge(width = 0.6), size = 5) +
  labs(x = "Accession",
       y = expression("Seeds plant"^{-1}~"(g)")) +
  scale_y_continuous(limits = c(0,250)) +
  theme_1

plot_yield_1

ls2 <- emmeans(model_yield, ~ Treatment, type="response")
cld2 <- as.data.frame(cld(ls2, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld2$.group <- str_trim(cld2$.group)
cld2$Treatment <- factor(cld2$Treatment,levels = c("Control", "Leaf extract", "Stem extract"))
cld2 <- cld2 %>% arrange(Treatment)

df_raw2 <- yield %>%
  group_by(Treatment) %>%
  summarise(mean = mean(Yield_per_plant, na.rm = TRUE),
            se = std.error(Yield_per_plant, na.rm = TRUE),
            .groups = "drop") %>% 
  left_join(cld2 %>% dplyr::select(Treatment, .group), by = c("Treatment"))

plot_yield_2 <- ggplot(df_raw2, aes(x=Treatment, y=mean, fill=Treatment)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.5, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position = position_dodge(width = 0.6),
                width = 0.2) +
  geom_text(aes(label = .group, y = mean+se + 0.05 * max(mean)), position = position_dodge(width = 0.6), size = 5) +
  labs(x = "Accession",
       y = expression("Seeds plant"^{-1}~"(g)"),
       fill = "Treatment") +
  scale_y_continuous(limits = c(0,250)) +
  scale_fill_manual(values = trat_colors) +
  theme_1

plot_yield_2

ls3 <- emmeans(model_yield, ~ Treatment | Accession_id, type="response")
cld3 <- as.data.frame(cld(ls3, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld3$.group <- str_trim(cld3$.group)
cld3 <- cld3 %>% arrange(Accession_id, Treatment)

df_raw3 <- yield %>%
  group_by(Treatment, Accession_id) %>%
  summarise(mean = mean(Yield_per_plant, na.rm = TRUE),
            se = std.error(Yield_per_plant, na.rm = TRUE),
            .groups = "drop") %>% 
  left_join(cld3 %>% dplyr::select(Treatment, Accession_id, .group), by = c("Treatment", "Accession_id"))

plot_yield_3 <- ggplot(df_raw3, aes(x=Accession_id, y=mean, fill=Treatment)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.5, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position = position_dodge(width = 0.6),
                width = 0.2) +
  geom_text(aes(label = .group, y = mean+se + 0.05 * max(mean)), position = position_dodge(width = 0.6), size = 5) +
  labs(x = "Accession",
       y = expression("Seeds plant"^{-1}~"(g)"),
       fill = "Treatment") +
  scale_y_continuous(limits = c(0,250)) +
  scale_fill_manual(values = trat_colors) +
  theme_1

plot_yield_3

### VI. Figure construction and export #####
######  Figure 1 #####

p1 <- plot_rgr_diameter + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

p2 <- plot_rgr_height +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank()) +
  theme_legend_vertical

p3 <- plot_rgr_no_leaves +
  theme(strip.text.x = element_blank())

figure_1 <- p1 / p2 / p3 +
  plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)"))) &
  theme(plot.tag = element_text(size = 18))

figure_1

ggsave("Plots/Figure 1.svg",
       figure_1,
       width = 174,
       height = 234,
       units = "mm")

######  Figure 2 #####

p1 <- plot_leaf_area

p2 <- plot_nectaries +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

p3 <- plot_trichomes +
  theme(strip.text.x = element_blank())

p4 <- plot_hardness +
  theme_legend_left

figure_2 <- ((p1 + p2 / p3) + plot_layout(widths = c(1.4, 1))) / p4 +
  plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)", "(d)"))) &
  theme(plot.tag = element_text(size = 18))

figure_2

ggsave("Plots/Figure 2.svg",
       figure_2,
       width = 174,
       height = 234,
       units = "mm")

######  Figure 3 #####

figure_3 <- plot_yield_3 +
  theme_legend_vertical_left

figure_3

ggsave("Plots/Figure 3.svg",
       figure_3,
       width = 94,
       height = 124,
       units = "mm")

######  Figure 4 #####

p1 <- plot_herbivory_1

p2 <- plot_herbivory_2 +
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())

p3 <- plot_herbivory_3 +
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())

figure_4 <- p1 + p2 + p3 + plot_layout(widths = c(1, 0.7, 1)) +
  plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)"))) &
  theme(plot.tag = element_text(size = 18))

figure_4

ggsave("Plots/Figure 4.svg",
       figure_4,
       width = 210,
       height = 126,
       units = "mm")
