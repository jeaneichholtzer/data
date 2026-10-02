setwd("C:/Users/94/OneDrive - IT de Conkal/Tithonia_Phaseolus_project/data/Datos_brutos")
#setwd("D:/Jean/OneDrive - IT de Conkal/Tithonia_Phaseolus_project/data/Datos_brutos")

# Packages ######################################################################################################

#install.packages("writexl")

#paquetes <- c("factoextra","mice")

#install.packages(paquetes, dependencies = TRUE)

library(tidyverse)
library(ggsignif)
library(scales)
library(readxl)
library(writexl)
# Stats
library(car)
library(emmeans)
library(glmmTMB)
library(lme4)
library(lmerTest)
library(multcomp)
library(MuMIn)
library(rcompanion)
library(DHARMa)
library(SuppDists)
library(plotrix)
library(performance)
library(patchwork)
library(ggtext)
library(factoextra)
library(corrplot)
library(mice)
library(piecewiseSEM)
library(DiagrammeR)
library(bestNormalize)
set.seed(123)

# Themes ###################################################

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
                 #strip.text = element_text(margin = margin(2, 0, 0, 0), size = 16, face = "bold"),
                 strip.background = element_blank())

theme_legend_left <- theme(legend.position = c(0.02, 0.98),
                           legend.direction = "horizontal",
                           legend.justification = "left",
                           legend.box.just = "left",
                           legend.text = element_text(size = 12),
                           legend.title = element_blank(), #element_text(size = 13, face = "bold"),
                           legend.key.height = unit(0.3, "cm"),
                           legend.key.width = unit(0.6, "cm"),
                           legend.background = element_rect(fill = alpha("white", 0.8),
                                                            colour = NA))

theme_legend <- theme(legend.position = c(0.98, 0.98),
                      legend.direction = "horizontal",
                      legend.justification = "right",
                      legend.box.just = "right",
                      legend.text = element_text(size = 12),
                      legend.title = element_blank(), #element_text(size = 13, face = "bold"),
                      legend.key.height = unit(0.3, "cm"),
                      legend.key.width = unit(0.6, "cm"),
                      legend.background = element_rect(fill = alpha("white", 0.8),
                                                       colour = NA))

theme_legend_vertical <- theme(legend.position = c(0.98, 0.85),
                               legend.direction = "vertical",
                               legend.justification = "right",
                               legend.box.just = "right",
                               legend.text = element_text(size = 12),
                               legend.title = element_blank(), #element_text(size = 13, face = "bold"),
                               legend.key.height = unit(0.3, "cm"),
                               legend.key.width = unit(0.6, "cm"),
                               legend.key.spacing.y = unit(0.1, "cm"),
                               legend.background = element_rect(fill = alpha("white", 0.8),
                                                                colour = NA))

theme_legend_vertical_left <- theme(legend.position = c(0.02, 0.90),
                                    legend.direction = "vertical",
                                    legend.justification = "left",
                                    legend.box.just = "left",
                                    legend.text = element_text(size = 12),
                                    legend.title = element_blank(), #element_text(size = 13, face = "bold"),
                                    legend.key.height = unit(0.3, "cm"),
                                    legend.key.width = unit(0.6, "cm"),
                                    legend.key.spacing.y = unit(0.1, "cm"),
                                    legend.background = element_rect(fill = alpha("white", 0.8),
                                                                     colour = NA))

trat_labels <- c(
  "H2O"  = "Control", #expression(H[2]*O),
  "H100" = "Leaf extract",
  "T100" = "Stem extract"
)

#trat_colors <- c(
#  "H2O" = "#FFFFFF",
#  "H100" = "grey",
#  "T100" = "black"
#)

trat_colors <- c(
  "H2O" = "white",
  "H100" = "#4eb934",
  "T100" = "#e97034"
)

cosecha_labels = c("1" = "35",
                   "2" = "55",
                   "3" = "75")

# Load data ######################################################################################################################################################################
### BIOMASA SLA #####
###### LOAD DATA #####

biomasa <- read_excel("Biomasa_SLA.xlsx")
biomasa <- biomasa %>%
  mutate(
    Cosecha = as.factor(Cosecha),
    Fecha = as.Date(Fecha, format = "%m/%d/%Y"),
    Genotype = as.factor(Genotype),
    Origin = as.factor(Origin),
    Gen_id = as.factor(str_sub(Genotype, 1, 1)),
    Tratamiento = factor(Tratamiento, levels = c("H2O", "H25", "T25", "H50", "T50", "H75", "T75", "H100", "T100")),
    Linea = as.factor(Linea),
    Linea_2 = as.factor(Linea_2),
    Area_foliar_muestra = as.numeric(Area_foliar_muestra)
  )

# Calculations for SLA_2 and LMA_2 

# ---- Means from Cosecha 0

means_cosecha_0 <- biomasa %>%
  filter(Cosecha == "0") %>%
  group_by(Genotype) %>%
  summarise(mean_prev_area = mean(Area_foliar_muestra, na.rm = TRUE),
            mean_prev_biomasa = mean(Biomasa_hojas_muestra, na.rm = TRUE),
            .groups = "drop")

# ---- Means from previous cosechas

means_previous_cosecha <- biomasa %>%
  filter(Cosecha %in% c("1", "2")) %>%
  group_by(Genotype, Tratamiento, Cosecha) %>%
  summarise(mean_prev_area = mean(Area_foliar_muestra, na.rm = TRUE),
            mean_prev_biomasa = mean(Biomasa_hojas_muestra, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(Cosecha = as.character(as.numeric(as.character(Cosecha)) + 1))

# ---- SLA2 for Cosecha 1

sla_c1 <- biomasa %>%
  filter(Cosecha == "1") %>%
  left_join(means_cosecha_0,
            by = "Genotype") %>%
  mutate(SLA_2 = ((mean_prev_area / mean_prev_biomasa) + (Area_foliar_muestra / Biomasa_hojas_muestra)) / 2)

# ---- SLA2 for Cosecha 2 and 3

sla_c23 <- biomasa %>%
  filter(Cosecha %in% c("2", "3")) %>%
  left_join(means_previous_cosecha,
            by = c("Genotype", "Tratamiento", "Cosecha")) %>%
  mutate(SLA_2 = ((mean_prev_area / mean_prev_biomasa) +  (Area_foliar_muestra / Biomasa_hojas_muestra)) / 2)

# ---- LMA2 for Cosecha 1

lma_c1 <- biomasa %>%
  filter(Cosecha == "1") %>%
  left_join(means_cosecha_0,
            by = "Genotype") %>%
  mutate(LMA_2 = ((mean_prev_biomasa / mean_prev_area) + (Biomasa_hojas_muestra / Area_foliar_muestra)) / 2)

# ---- SLA2 for Cosecha 2 and 3

lma_c23 <- biomasa %>%
  filter(Cosecha %in% c("2", "3")) %>%
  left_join(means_previous_cosecha,
            by = c("Genotype", "Tratamiento", "Cosecha")) %>%
  mutate(LMA_2 = ((mean_prev_biomasa / mean_prev_area) +  (Biomasa_hojas_muestra / Area_foliar_muestra)) / 2)

# ---- Merge all together

biomasa_sla <- bind_rows(
  biomasa %>% filter(Cosecha == "0"),
  full_join(sla_c1,  lma_c1),
  full_join(sla_c23, lma_c23)
) %>%
  arrange(Genotype, Linea, Tratamiento, Cosecha)

biomasa_filtrado <- biomasa_sla %>%
  filter(Cosecha %in% c("1","2","3"),
         !Tratamiento %in% c("H25", "T25", "H50", "T50", "H75", "T75"),
         !Gen_id %in% c("I")) %>%
  mutate(Cosecha = factor(Cosecha),
         Origin = factor(Origin),
         Tratamiento = factor(Tratamiento, levels = c("H2O", "H100", "T100")),
         Genotype = factor(Genotype),
         Linea = factor(Linea),
         Gen_id = factor(Gen_id))

biomasa_domesticados <- biomasa_filtrado %>%
  filter(Origin=="Domesticated") %>%
  mutate(Linea = factor(Linea),
         Gen_id = factor(Gen_id))

biomasa_silvestres <- biomasa_filtrado %>%
  filter(Origin=="Wild") %>%
  mutate(Linea = factor(Linea),
         Gen_id = factor(Gen_id))

###### SLA #####

set.seed(123)
BNobject <- bestNormalize(biomasa_filtrado$SLA_2)
biomasa_filtrado$SLA_2_tr <- predict(BNobject)

model1 <- glmmTMB(SLA_2 ~ Cosecha * Origin * Tratamiento + (1|Gen_id) + (1|Linea),
                  family = Gamma(link="log"),
                  data = biomasa_filtrado)

part1_model1 <- glmmTMB(SLA_2 ~ Cosecha * Origin * Tratamiento + (1|Linea),
                        family = Gamma(link="log"),
                        data = biomasa_filtrado)

part2_model1 <- glmmTMB(SLA_2 ~ Cosecha * Origin * Tratamiento + (1|Gen_id),
                        family = Gamma(link="log"),
                        data = biomasa_filtrado)

summary(model1)
Anova(model1, type="II")
anova(part1_model1, model1)
anova(part2_model1, model1)

### DEFENSE #####
###### LOAD DATA #####

defensas <- read_excel("Defensas.xlsx")
defensas <- defensas %>%
  mutate(
    Cosecha = as.factor(Cosecha),
    Origin = as.factor(Origin),
    Genotype = as.factor(Genotype),
    Gen_id = as.factor(str_sub(Genotype, 1, 1)),
    Tratamiento = factor(Tratamiento, levels = c("H2O", "H25", "T25", "H50", "T50", "H75", "T75", "H100", "T100")),
    Linea = as.factor(Linea),
    Linea_2 = as.factor(Linea_2),
    Muestreador = as.factor(Muestreador),
    Hoja_id = as.factor(Hoja_id)
  )

defensas_filtrado <- defensas %>%
  filter(Cosecha %in% c("1","2","3"),
         !Tratamiento %in% c("H25", "T25", "H50", "T50", "H75", "T75"),
         !Gen_id %in% c("I")) %>%
  mutate(Cosecha = factor(Cosecha),
         Origin = factor(Origin),
         Tratamiento = factor(Tratamiento, levels = c("H2O", "H100", "T100")),
         Genotype = factor(Genotype),
         Gen_id = factor(Gen_id))

defensas_domesticados <- defensas_filtrado %>%
  filter(Origin=="Domesticated")

defensas_silvestres <- defensas_filtrado %>%
  filter(Origin=="Wild")

#write.csv(defensas_filtrados,"defensas_filtrados.csv")

###### THICKNESS #####

model1 <- glmmTMB(Grosor ~ Cosecha * Origin * Tratamiento + (1|Gen_id) + (1|Linea),
                  family =  gaussian(),
                  data = defensas_filtrado,
                  dispformula = ~ Tratamiento + Origin)

part1_model1 <- glmmTMB(Grosor ~ Cosecha * Origin * Tratamiento + (1|Linea),
                        family =  gaussian(),
                        data = defensas_filtrado,
                        dispformula = ~ Tratamiento + Origin)

part2_model1 <- glmmTMB(Grosor ~ Cosecha * Origin * Tratamiento + (1|Gen_id),
                        family =  gaussian(),
                        data = defensas_filtrado,
                        dispformula = ~ Tratamiento + Origin)

summary(model1)
Anova(model1, type="III")
anova(part1_model1, model1)
anova(part2_model1, model1)

ls1 <- emmeans(model1, ~ Tratamiento | Cosecha | Origin, type="response")
cld1 <- as.data.frame(cld(ls1, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld1$.group <- str_trim(cld1$.group) # Remove space
cld1$Tratamiento <- factor(cld1$Tratamiento,levels = c("H2O", "H100", "T100"))
cld1 <- cld1 %>% arrange(Origin, Cosecha, Tratamiento)

df_raw1 <- defensas_filtrado %>%
  group_by(Cosecha, Tratamiento, Origin) %>%
  summarise(
    mean = mean(Grosor, na.rm = TRUE),
    se   = std.error(Grosor, na.rm = TRUE),
    .groups = "drop") %>% 
  left_join(cld1 %>% dplyr::select(Cosecha, Tratamiento, Origin, .group), by = c("Cosecha", "Tratamiento", "Origin"))

plot1 <- ggplot(df_raw1, aes(x = Cosecha, y = mean, fill = Tratamiento)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se,
                    ymax = mean+se),
                position = position_dodge(width = 0.8),
                width = 0.2,
                size = 0.5) +
  facet_wrap(~ Origin) +
  geom_text(aes(label = .group, y = mean+se + 0.05 * max(mean)), position = position_dodge(width = 0.8), size = 5) +
  labs(x = "Time (days after sowing)",
       y = "Leaf thickness (mm)",
       fill = "Treatment") +
  scale_x_discrete(labels = cosecha_labels) +
  scale_fill_manual(values = trat_colors, labels = trat_labels) +
  theme_1 +
  theme_legend_left +
  theme(legend.position=c(0.02,0.95))

plot_thickness <- plot1

ls1 <- emmeans(model1, ~ Cosecha | Origin, type="response")
cld1 <- as.data.frame(cld(ls1, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld1$.group <- str_trim(cld1$.group) # Remove space
cld1 <- cld1 %>% arrange(Origin, Cosecha)
cld1

df_raw1 <- defensas_filtrado %>%
  group_by(Cosecha, Origin) %>%
  summarise(
    mean = mean(Grosor, na.rm = TRUE),
    se   = std.error(Grosor, na.rm = TRUE),
    .groups = "drop") %>% 
  left_join(cld1 %>% dplyr::select(Cosecha, Origin, .group), by = c("Cosecha", "Origin"))

plot1 <- ggplot(df_raw1, aes(x = Cosecha, y = mean)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-0.001,
                    ymax = mean+se),
                position = position_dodge(width = 0.8),
                width = 0.2,
                size = 0.5) +
  facet_wrap(~ Origin) +
  geom_text(aes(label = .group, y = mean+se + 0.05 * max(mean)), position = position_dodge(width = 0.8), size = 5) +
  labs(x = "Time (days after sowing)",
       y = "Leaf thickness (mm)",
       fill = "Treatment") +
  scale_x_discrete(labels = cosecha_labels) +
  #scale_fill_manual(values = trat_colors, labels = trat_labels) +
  theme_1 +
  theme_legend_left +
  theme(legend.position=c(0.02,0.95))

plot_thickness2 <- plot1

###### TOUGHNESS #####

set.seed(123)
BNobject <- bestNormalize(defensas_filtrado$Dureza)
defensas_filtrado$Dureza_tr <- predict(BNobject)

model1 <- glmmTMB(Dureza_tr ~ Cosecha * Origin * Tratamiento + (1|Gen_id) + (1|Linea),
                  family =  t_family(),
                  data = defensas_filtrado,
                  dispformula = ~ Cosecha)

part1_model1 <- glmmTMB(Dureza_tr ~ Cosecha * Origin * Tratamiento + (1|Linea),
                        family =  t_family(),
                        data = defensas_filtrado,
                        dispformula = ~ Cosecha)

part2_model1 <- glmmTMB(Dureza_tr ~ Cosecha * Origin * Tratamiento + (1|Gen_id),
                        family =  t_family(),
                        data = defensas_filtrado,
                        dispformula = ~ Cosecha)

summary(model1)
Anova(model1, type="III")
anova(part1_model1, model1)
anova(part2_model1, model1)

ls1 <- emmeans(model1, ~ Tratamiento | Cosecha | Origin, type="response")
cld1 <- as.data.frame(cld(ls1, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld1$.group <- str_trim(cld1$.group) # Remove space
cld1$Tratamiento <- factor(cld1$Tratamiento,levels = c("H2O", "H100", "T100"))
cld1 <- cld1 %>% arrange(Origin, Cosecha, Tratamiento)
cld1

df_raw1 <- defensas_filtrado %>%
  group_by(Cosecha, Tratamiento, Origin) %>%
  summarise(
    mean = mean(Dureza, na.rm = TRUE),
    se   = std.error(Dureza, na.rm = TRUE),
    .groups = "drop") %>% 
  left_join(cld1 %>% dplyr::select(Cosecha, Tratamiento, Origin, .group), by = c("Cosecha", "Tratamiento", "Origin"))

plot1 <- ggplot(df_raw1, aes(x = Cosecha, y = mean, fill = Tratamiento)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se,
                    ymax = mean+se),
                position = position_dodge(width = 0.8),
                width = 0.2,
                size = 0.5) +
  facet_wrap(~ Origin) +
  geom_text(aes(label = .group, y = mean+se + 0.05 * max(mean)), position = position_dodge(width = 0.8), size = 5) +
  labs(x = "Time (days after sowing)",
       y = expression(Leaf~hardness~(kg~cm^{-2})),
       fill = "Treatment") +
  scale_x_discrete(labels = cosecha_labels) +
  scale_fill_manual(values = trat_colors, labels = trat_labels) +
  theme_1

plot_toughness <- plot1

### GROWTH #####
###### LOAD DATA #####

crecimiento <- read_excel("Crecimiento.xlsx")
crecimiento <- crecimiento %>%
  mutate(
    Cosecha = as.factor(Cosecha),
    Fecha = as.Date(Fecha, format = "%m/%d/%Y"),
    Origin = as.factor(Origin),
    Genotype = as.factor(Genotype),
    Gen_id = as.factor(str_sub(Genotype, 1, 1)),
    Tratamiento = factor(Tratamiento, levels = c("H2O", "H25", "T25", "H50", "T50", "H75", "T75", "H100", "T100")),
    Linea = as.factor(Linea),
    Linea_2 = as.factor(Linea_2)
  )

# Cálculo de los RGRs

# ---- Means from Cosecha 0

means_cosecha_0 <- crecimiento %>%
  filter(Cosecha == "0") %>%
  group_by(Genotype) %>%
  summarise(mean_log_prev_diametro = mean(log(Diametro), na.rm = TRUE),
            mean_log_prev_no_hojas = mean(log(No_hojas), na.rm = TRUE),
            mean_prev_no_hojas = mean(No_hojas, na.rm = TRUE),
            mean_log_prev_altura   = mean(log(Altura), na.rm = TRUE),
            fecha_prev = first(Fecha),
            .groups = "drop")

# ---- Means from previous cosechas

means_previous_cosecha <- crecimiento %>%
  filter(Cosecha %in% c("1", "2")) %>%
  group_by(Genotype, Tratamiento, Cosecha) %>%
  summarise(mean_log_prev_diametro = mean(log(Diametro), na.rm = TRUE),
            mean_log_prev_no_hojas = mean(log(No_hojas), na.rm = TRUE),
            mean_prev_no_hojas = mean(No_hojas, na.rm = TRUE),
            mean_log_prev_altura   = mean(log(Altura), na.rm = TRUE),
            fecha_prev = first(Fecha),
            .groups = "drop") %>%
  mutate(Cosecha = as.character(as.numeric(as.character(Cosecha)) + 1))

# ---- RGR for Cosecha 1

rgr_c1 <- crecimiento %>%
  filter(Cosecha == "1") %>%
  left_join(means_cosecha_0,
            by = "Genotype") %>%
  mutate(days_diff = as.numeric(Fecha - fecha_prev),
         RGR_Diametro = (log(Diametro) - mean_log_prev_diametro) / days_diff,
         RGR_No_hojas = (log(No_hojas) - mean_log_prev_no_hojas) / days_diff,
         LPR = (No_hojas - mean_prev_no_hojas) / days_diff,
         RGR_Altura = (log(Altura) - mean_log_prev_altura) / days_diff)

# ---- RGR for Cosecha 2 and 3

rgr_c23 <- crecimiento %>%
  filter(Cosecha %in% c("2", "3")) %>%
  left_join(means_previous_cosecha,
            by = c("Genotype", "Tratamiento", "Cosecha")) %>%
  mutate(days_diff = as.numeric(Fecha - fecha_prev),
         RGR_Diametro = (log(Diametro) - mean_log_prev_diametro) / days_diff,
         RGR_No_hojas = (log(No_hojas) - mean_log_prev_no_hojas) / days_diff,
         LPR = (No_hojas - mean_prev_no_hojas) / days_diff,
         RGR_Altura = (log(Altura) - mean_log_prev_altura) / days_diff)

# ---- Merge all together

crecimiento_rgr <- bind_rows(crecimiento %>% filter(Cosecha == "0"), rgr_c1, rgr_c23) %>%
  arrange(Genotype, Linea, Tratamiento, Cosecha)

crecimiento_filtrado <- crecimiento_rgr %>%
  filter(Cosecha %in% c("1","2","3"),
         !Tratamiento %in% c("H25", "T25", "H50", "T50", "H75", "T75"),
         !Gen_id %in% c("I")) %>%
  mutate(Cosecha = factor(Cosecha),
         Origin = factor(Origin),
         Tratamiento = factor(Tratamiento, levels = c("H2O", "H100", "T100")),
         Genotype = factor(Genotype),
         Gen_id = factor(Gen_id))

crecimiento_domesticados <- crecimiento_filtrado %>%
  filter(Origin=="Domesticated")

crecimiento_silvestres <- crecimiento_filtrado %>%
  filter(Origin=="Wild")

###### Diameter #####

model1 <- glmmTMB(Diametro ~ Cosecha * Origin * Tratamiento + (1|Gen_id) + (1|Linea),
                  family = gaussian(),
                  data = crecimiento_filtrado,
                  dispformula = ~ Cosecha * Origin)

part1_model1 <- glmmTMB(Diametro ~ Cosecha * Origin * Tratamiento + (1|Linea),
                        family = gaussian(),
                        data = crecimiento_filtrado,
                        dispformula = ~ Cosecha * Origin)

part2_model1 <- glmmTMB(Diametro ~ Cosecha * Origin * Tratamiento + (1|Gen_id),
                        family = gaussian(),
                        data = crecimiento_filtrado,
                        dispformula = ~ Cosecha * Origin)

summary(model1)
Anova(model1, type="III")
anova(part1_model1, model1)
anova(part2_model1, model1)

plot(simulateResiduals(model1))

ls1 <- emmeans(model1, ~ Tratamiento | Cosecha | Origin, type="response")
cld1 <- as.data.frame(cld(ls1, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld1$.group <- str_trim(cld1$.group) # Remove space
cld1$Tratamiento <- factor(cld1$Tratamiento,levels = c("H2O", "H100", "T100"))
cld1 <- cld1 %>% arrange(Origin, Cosecha, Tratamiento)

df_raw1 <- crecimiento_filtrado %>%
  group_by(Cosecha, Tratamiento, Origin) %>%
  summarise(
    mean = mean(Diametro, na.rm = TRUE),
    se   = std.error(Diametro, na.rm = TRUE),
    .groups = "drop") %>% 
  left_join(cld1 %>% dplyr::select(Cosecha, Tratamiento, Origin, .group), by = c("Cosecha", "Tratamiento", "Origin"))

plot1 <- ggplot(df_raw1, aes(x = Cosecha, y = mean, fill = Tratamiento)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se,
                    ymax = mean+se),
                position = position_dodge(width = 0.8),
                width = 0.2,
                size = 0.5) +
  facet_wrap(~ Origin) +
  geom_text(aes(label = .group, y = mean+se + 0.05 * max(mean)), position = position_dodge(width = 0.8), size = 5) +
  labs(x = "Time (days after sowing)",
       y = "Stem diameter (mm)",
       fill = "Treatment") +
  scale_x_discrete(labels = cosecha_labels) +
  scale_y_continuous(limits=c(0,15)) +
  scale_fill_manual(values = trat_colors, labels = trat_labels) +
  theme_1

plot_diameter <- plot1

###### Height #####

set.seed(123)
BNobject <- bestNormalize(crecimiento_filtrado$Altura)
crecimiento_filtrado$Altura_tr <- predict(BNobject)

model1 <- glmmTMB(Altura_tr ~ Cosecha * Origin * Tratamiento + (1|Gen_id) + (1|Linea),
                  family = gaussian(),
                  dispformula = ~ Cosecha,
                  data = crecimiento_filtrado)

part1_model1 <- glmmTMB(Altura_tr ~ Cosecha * Origin * Tratamiento + (1|Linea),
                        family = gaussian(),
                        dispformula = ~ Cosecha,
                        data = crecimiento_filtrado)

part2_model1 <- glmmTMB(Altura_tr ~ Cosecha * Origin * Tratamiento + (1|Gen_id),
                        family = gaussian(),
                        dispformula = ~ Cosecha,
                        data = crecimiento_filtrado)

summary(model1)
Anova(model1, type="III")
anova(part1_model1, model1)
anova(part2_model1, model1)

ls1 <- emmeans(model1, ~ Tratamiento | Cosecha, type="response")
cld1 <- as.data.frame(cld(ls1, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld1$.group <- str_trim(cld1$.group) # Remove space
cld1$Tratamiento <- factor(cld1$Tratamiento,levels = c("H2O", "H100", "T100"))
cld1 <- cld1 %>% arrange(Cosecha, Tratamiento)

df_raw1 <- crecimiento_filtrado %>%
  group_by(Cosecha, Tratamiento) %>%
  summarise(
    mean = mean(Altura, na.rm = TRUE),
    se   = std.error(Altura, na.rm = TRUE),
    .groups = "drop") %>% 
  left_join(cld1 %>% dplyr::select(Cosecha, Tratamiento, .group), by = c("Cosecha", "Tratamiento"))

plot1 <- ggplot(df_raw1, aes(x = Cosecha, y = mean, fill = Tratamiento)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se,
                    ymax = mean+se),
                position = position_dodge(width = 0.8),
                width = 0.2,
                size = 0.5) +
  #facet_wrap(~ Origin) +
  geom_text(aes(label = .group, y = mean+se + 0.05 * max(mean)), position = position_dodge(width = 0.8), size = 5) +
  labs(x = "Time (days after sowing)",
       y = "Length (cm)",
       fill = "Treatment") +
  scale_x_discrete(labels = cosecha_labels) +
  scale_fill_manual(values = trat_colors, labels = trat_labels) +
  theme_1

plot_height <- plot1
plot_height

###### Leaf number #####

set.seed(123)
BNobject <- bestNormalize(crecimiento_filtrado$No_hojas)
crecimiento_filtrado$No_hojas_tr <- predict(BNobject)

model1 <- glmmTMB(No_hojas_tr ~ Cosecha * Origin * Tratamiento + (1|Gen_id) + (1|Linea),
                  family = gaussian(),
                  data = crecimiento_filtrado,
                  dispformula = ~ Cosecha * Origin * Tratamiento)

part1_model1 <- glmmTMB(No_hojas_tr ~ Cosecha * Origin * Tratamiento + (1|Linea),
                        family = gaussian(),
                        data = crecimiento_filtrado,
                        dispformula = ~ Cosecha * Origin * Tratamiento)

part2_model1 <- glmmTMB(No_hojas_tr ~ Cosecha * Origin * Tratamiento + (1|Gen_id),
                        family = gaussian(),
                        data = crecimiento_filtrado,
                        dispformula = ~ Cosecha * Origin * Tratamiento)

summary(model1)
Anova(model1, type="III")
anova(part1_model1, model1)
anova(part2_model1, model1)

plot(simulateResiduals(model1))

ls1 <- emmeans(model1, ~ Tratamiento | Origin, type="response")
cld1 <- as.data.frame(cld(ls1, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld1$.group <- str_trim(cld1$.group) # Remove space
cld1$Tratamiento <- factor(cld1$Tratamiento,levels = c("H2O", "H100", "T100"))
cld1 <- cld1 %>% arrange(Origin, Tratamiento)

df_raw1 <- crecimiento_filtrado %>%
  group_by(Tratamiento, Origin) %>%
  summarise(
    mean = mean(No_hojas, na.rm = TRUE),
    se   = std.error(No_hojas, na.rm = TRUE),
    .groups = "drop") %>% 
  left_join(cld1 %>% dplyr::select(Tratamiento, Origin, .group), by = c("Tratamiento", "Origin"))

plot1 <- ggplot(df_raw1, aes(x = Origin, y = mean, fill = Tratamiento)) +
  geom_col(position = position_dodge(width = 0.5), width = 0.3, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se,
                    ymax = mean+se),
                position = position_dodge(width = 0.5),
                width = 0.12,
                size = 0.5) +
  #facet_wrap(~ Origin) +
  geom_text(aes(label = .group, y = mean+se + 0.05 * max(mean)), position = position_dodge(width = 0.5), size = 5) +
  labs(x = "Status",
       y = "Number of leaves",
       fill = "Treatment") +
  scale_x_discrete(labels = cosecha_labels) +
  #scale_y_continuous(limits=c(0,15)) +
  scale_fill_manual(values = trat_colors, labels = trat_labels) +
  theme_1

plot_leaf_number <- plot1
plot_leaf_number

###### RGR Diameter #####

set.seed(123)
BNobject <- bestNormalize(crecimiento_filtrado$RGR_Diametro)
crecimiento_filtrado$RGR_Diametro_tr <- predict(BNobject)

model1 <- glmmTMB(RGR_Diametro_tr ~ Cosecha * Origin * Tratamiento + (1|Gen_id) + (1|Linea),
                  family = t_family(),
                  dispformula = ~ Cosecha,
                  data = crecimiento_filtrado)

part1_model1 <- glmmTMB(RGR_Diametro_tr ~ Cosecha * Origin * Tratamiento + (1|Linea),
                        family = t_family(),
                        dispformula = ~ Cosecha,
                        data = crecimiento_filtrado)

part2_model1 <- glmmTMB(RGR_Diametro_tr ~ Cosecha * Origin * Tratamiento + (1|Gen_id),
                        family = t_family(),
                        dispformula = ~ Cosecha,
                        data = crecimiento_filtrado)

#plot(simulateResiduals(model1))

summary(model1)
Anova(model1, type="III")
anova(part1_model1, model1)
anova(part2_model1, model1)

ls1 <- emmeans(model1, ~ Tratamiento | Cosecha | Origin, type="response")
cld1 <- as.data.frame(cld(ls1, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld1$.group <- str_trim(cld1$.group) # Remove space
cld1$Tratamiento <- factor(cld1$Tratamiento,levels = c("H2O", "H100", "T100"))
cld1 <- cld1 %>% arrange(Origin, Cosecha, Tratamiento)

df_raw1 <- crecimiento_filtrado %>%
  group_by(Cosecha, Tratamiento, Origin) %>%
  summarise(
    mean = mean(RGR_Diametro, na.rm = TRUE),
    se   = std.error(RGR_Diametro, na.rm = TRUE),
    .groups = "drop") %>% 
  left_join(cld1 %>% dplyr::select(Cosecha, Tratamiento, Origin, .group), by = c("Cosecha", "Tratamiento", "Origin"))

plot1 <- ggplot(df_raw1, aes(x = Cosecha, y = mean, fill = Tratamiento)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se,
                    ymax = mean+se),
                position = position_dodge(width = 0.8),
                width = 0.2,
                size = 0.5) +
  facet_wrap(~ Origin) +
  geom_text(aes(label = .group, y = mean+se + 0.1 * max(mean)), position = position_dodge(width = 0.8), size = 5) +
  labs(x = "Time (days after sowing)",
       y = expression(RGR[italic(D)]~(mm~mm^{-1}~day^{-1})),
       fill = "Treatment") +
  scale_x_discrete(labels = cosecha_labels) +
  scale_fill_manual(values = trat_colors, labels = trat_labels) +
  theme_1

plot_rgr_diameter <- plot1

plot_test <- ggplot(cld1, aes(x = Cosecha, y = emmean, fill = Tratamiento)) +
  geom_errorbar(aes(ymin = emmean-SE,
                    ymax = emmean+SE),
                position = position_dodge(width = 0.8),
                width = 0.2,
                size = 0.5) +
  geom_point(position = position_dodge(width = 0.8), colour = "black") +
  facet_wrap(~ Origin) +
  geom_text(aes(label = .group, y = emmean+SE + 0.1), position = position_dodge(width = 0.8), size = 5) +
  labs(x = "Time (days after sowing)",
       y = expression(RGR[italic(D)]~(mm~mm^{-1}~day^{-1})),
       fill = "Treatment") +
  scale_x_discrete(labels = cosecha_labels) +
  scale_fill_manual(values = trat_colors, labels = trat_labels) +
  theme_1

###### RGR Height #####

set.seed(123)
BNobject <- bestNormalize(crecimiento_filtrado$RGR_Altura)
crecimiento_filtrado$RGR_Altura_tr <- predict(BNobject)

model1 <- glmmTMB(RGR_Altura_tr ~ Cosecha * Origin * Tratamiento + (1|Gen_id) + (1|Linea),
                  family = t_family(),
                  data = crecimiento_filtrado,
                  dispformula = ~ Cosecha * Tratamiento * Origin)

part1_model1 <- glmmTMB(RGR_Altura_tr ~ Cosecha * Origin * Tratamiento + (1|Linea),
                        family = t_family(),
                        data = crecimiento_filtrado,
                        dispformula = ~ Cosecha * Tratamiento * Origin)

part2_model1 <- glmmTMB(RGR_Altura_tr ~ Cosecha * Origin * Tratamiento + (1|Gen_id),
                        family = t_family(),
                        data = crecimiento_filtrado,
                        dispformula = ~ Cosecha * Tratamiento * Origin)

#plot(simulateResiduals(model1))

summary(model1)
Anova(model1, type="III")
anova(part1_model1, model1)
anova(part2_model1, model1)

ls1 <- emmeans(model1, ~ Tratamiento | Cosecha | Origin, type="response")
cld1 <- as.data.frame(cld(ls1, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld1$.group <- str_trim(cld1$.group) # Remove space
cld1$Tratamiento <- factor(cld1$Tratamiento,levels = c("H2O", "H100", "T100"))
cld1 <- cld1 %>% arrange(Origin, Cosecha, Tratamiento)

df_raw1 <- crecimiento_filtrado %>%
  group_by(Cosecha, Tratamiento, Origin) %>%
  summarise(
    mean = mean(RGR_Altura, na.rm = TRUE),
    se   = std.error(RGR_Altura, na.rm = TRUE),
    .groups = "drop") %>% 
  left_join(cld1 %>% dplyr::select(Cosecha, Tratamiento, Origin, .group), by = c("Cosecha", "Tratamiento", "Origin"))

plot1 <- ggplot(df_raw1, aes(x = Cosecha, y = mean, fill = Tratamiento)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se,
                    ymax = mean+se),
                position = position_dodge(width = 0.8),
                width = 0.2,
                size = 0.5) +
  facet_wrap(~ Origin) +
  geom_text(aes(label = .group, y = mean+se + 0.1 * max(mean)), position = position_dodge(width = 0.8), size = 5) +
  labs(x = "Time (days after sowing)",
       y = expression(RGR[italic(L)]~(cm~cm^{-1}~day^{-1})),
       fill = "Treatment") +
  scale_x_discrete(labels = cosecha_labels) +
  scale_fill_manual(values = trat_colors, labels = trat_labels) +
  theme_1

plot_rgr_height <- plot1

###### RGR Leaf number #####

set.seed(123)
BNobject <- bestNormalize(crecimiento_filtrado$RGR_No_hojas)
crecimiento_filtrado$RGR_No_hojas_tr <- predict(BNobject)

model1 <- glmmTMB(RGR_No_hojas_tr ~ Cosecha * Origin * Tratamiento + (1|Gen_id) + (1|Linea),
                  family = gaussian(),
                  dispformula = ~ Cosecha * Origin * Tratamiento,
                  data = crecimiento_filtrado)

part1_model1 <- glmmTMB(RGR_No_hojas_tr ~ Cosecha * Origin * Tratamiento + (1|Linea),
                        family = gaussian(),
                        dispformula = ~ Cosecha * Origin * Tratamiento,
                        data = crecimiento_filtrado)

part2_model1 <- glmmTMB(RGR_No_hojas_tr ~ Cosecha * Origin * Tratamiento + (1|Gen_id),
                        family = gaussian(),
                        dispformula = ~ Cosecha * Origin * Tratamiento,
                        data = crecimiento_filtrado)

plot(simulateResiduals(model1))

summary(model1)
Anova(model1, type="III")
anova(part1_model1, model1)
anova(part2_model1, model1)

ls1 <- emmeans(model1, ~ Tratamiento | Cosecha | Origin, type="response")
cld1 <- as.data.frame(cld(ls1, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld1$.group <- str_trim(cld1$.group) # Remove space
cld1$Tratamiento <- factor(cld1$Tratamiento,levels = c("H2O", "H100", "T100"))
cld1 <- cld1 %>% arrange(Origin, Cosecha, Tratamiento)
cld1
pairs(ls1)

df_raw1 <- crecimiento_filtrado %>%
  group_by(Cosecha, Tratamiento, Origin) %>%
  summarise(
    mean = mean(RGR_No_hojas, na.rm = TRUE),
    se   = std.error(RGR_No_hojas, na.rm = TRUE),
    .groups = "drop") %>% 
  left_join(cld1 %>% dplyr::select(Cosecha, Tratamiento, Origin, .group), by = c("Cosecha", "Tratamiento", "Origin"))

plot1 <- ggplot(df_raw1, aes(x = Cosecha, y = mean, fill = Tratamiento)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se,
                    ymax = mean+se),
                position = position_dodge(width = 0.8),
                width = 0.2,
                size = 0.5) +
  facet_wrap(~ Origin) +
  geom_text(aes(label = .group, y = mean+se + 0.1 * max(mean)), position = position_dodge(width = 0.8), size = 5) +
  labs(x = "Time (days after sowing)",
       y = expression(RGR[italic(NL)]~(no.~no.^{-1}~day^{-1})),
       fill = "Treatment") +
  scale_x_discrete(labels = cosecha_labels) +
  scale_fill_manual(values = trat_colors, labels = trat_labels) +
  theme_1

plot_rgr_leaf_number <- plot1

### LEAF DAMAGE #####

daño <- read_excel("Daño_foliar.xlsx")

daño <- daño %>%
  mutate(Count = Damaged_obs + Intact_obs,
         #Damaged = round(Incidencia * 5),       # successes
         #Intact   = 5 - Damaged,
         Cosecha = as.factor(Cosecha),
         Genotype = as.factor(Genotype),
         Gen_id = as.factor(Gen_id),
         Origin = as.factor(Origin),
         Tratamiento = factor(Tratamiento, levels = c("H2O", "H25", "T25", "H50", "T50", "H75", "T75", "H100", "T100")),
         Linea = as.factor(Linea),
         Linea_2 = as.factor(Linea_2))

daño_filtrado <- daño %>%
  filter(!Tratamiento %in% c("H25", "T25", "H50", "T50", "H75", "T75"),
         !Gen_id %in% c("I")) %>%
  #!Gen_id %in% c("C","E","F","I")) %>%
  #Gen_id %in% c("A","B","D","G","H","I")) %>%
  mutate(Tratamiento = factor(Tratamiento, levels = c("H2O", "H100", "T100")),
         Genotype = factor(Genotype),
         Gen_id = factor(Gen_id),
         Origin = factor(Origin))

trat_labels_herbivory <- c(
  "H2O"  = "Control", #expression(H[2]*O),
  "H100" = "Leaf",
  "T100" = "Stem"
)

model1 <- glmmTMB(cbind(Damaged_obs, Intact_obs) ~ Cosecha * Origin * Tratamiento + (1|Gen_id) + (1|Linea),
                  family=binomial(link="logit"),
                  data = daño_filtrado)

part1_model1 <- glmmTMB(cbind(Damaged_obs, Intact_obs) ~ Cosecha * Origin * Tratamiento + (1|Linea),
                        family=binomial(link="logit"),
                        data = daño_filtrado)

part2_model1 <- glmmTMB(cbind(Damaged_obs, Intact_obs) ~ Cosecha * Origin * Tratamiento + (1|Gen_id),
                        family=binomial(link="logit"),
                        data = daño_filtrado)

#plot(simulateResiduals(model1))

summary(model1)
Anova(model1, type="II")
anova(part1_model1, model1)
anova(part2_model1, model1)

ls1 <- emmeans(model1, ~ Cosecha, type="response")
cld1 <- as.data.frame(cld(ls1, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld1$.group <- str_trim(cld1$.group) # Remove space
cld1 <- cld1 %>% arrange(Cosecha)
cld1

df_raw1 <- daño_filtrado %>%
  group_by(Cosecha) %>%
  summarise(
    mean = mean(Incidencia, na.rm = TRUE),
    se   = std.error(Incidencia, na.rm = TRUE),
    .groups = "drop") %>% 
  left_join(cld1 %>% dplyr::select(Cosecha, .group), by = c("Cosecha"))

plot1 <- ggplot(df_raw1, aes(x = Cosecha, y = mean)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6, colour = "black", fill="grey95", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se,
                    ymax = mean+se),
                position = position_dodge(width = 0.8),
                width = 0.2,
                size = 0.5) +
  geom_text(aes(label = .group, y = mean+se + 0.05 * max(mean)), position = position_dodge(width = 0.8), size = 5) +
  labs(x = "Time (days after sowing)",
       y = "Herbivory (%)") +
  scale_x_discrete(labels = cosecha_labels) +
  scale_y_continuous(labels = scales::percent, limits=c(0,0.80)) +
  theme_1

plot_herbivory_1 <- plot1
plot_herbivory_1

ls2 <- emmeans(model1, ~ Origin, type="response")
cld2 <- as.data.frame(cld(ls2, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld2$.group <- str_trim(cld2$.group) # Remove space
cld2 <- cld2 %>% arrange(Origin)
cld2

df_raw2 <- daño_filtrado %>%
  group_by(Origin) %>%
  summarise(
    mean = mean(Incidencia, na.rm = TRUE),
    se   = std.error(Incidencia, na.rm = TRUE),
    .groups = "drop") %>% 
  left_join(cld2 %>% dplyr::select(Origin, .group), by = c("Origin"))

plot2 <- ggplot(df_raw2, aes(x = Origin, y = mean)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6, colour = "black", fill="grey95", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se,
                    ymax = mean+se),
                position = position_dodge(width = 0.8),
                width = 0.2,
                size = 0.5) +
  #facet_wrap(~ Origin) +
  geom_text(aes(label = .group, y = mean+se + 0.06 * max(mean)), position = position_dodge(width = 0.8), size = 5) +
  labs(x = "Status",
       y = "Herbivory (%)") +
  scale_y_continuous(labels = scales::percent, limits=c(0,0.80)) +
  theme_1

plot_herbivory_2 <- plot2
plot_herbivory_2

ls3 <- emmeans(model1, ~ Tratamiento, type="response")
cld3 <- as.data.frame(cld(ls3, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld3$.group <- str_trim(cld3$.group) # Remove space
cld3$Tratamiento <- factor(cld3$Tratamiento,levels = c("H2O", "H100", "T100"))
cld3 <- cld3 %>% arrange(Tratamiento)
cld3

df_raw3 <- daño_filtrado %>%
  group_by(Tratamiento) %>%
  summarise(
    mean = mean(Incidencia, na.rm = TRUE),
    se   = std.error(Incidencia, na.rm = TRUE),
    .groups = "drop") %>% 
  left_join(cld3 %>% dplyr::select(Tratamiento, .group), by = c("Tratamiento"))

plot3 <- ggplot(df_raw3, aes(x = Tratamiento, y = mean, fill = Tratamiento)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se,
                    ymax = mean+se),
                position = position_dodge(width = 0.8),
                width = 0.2,
                size = 0.5) +
  geom_text(aes(label = .group, y = mean+se + 0.06 * max(mean)), position = position_dodge(width = 0.8), size = 5) +
  labs(x = "Treatment",
       y = "Herbivory (%)",
       fill = "Treatment") +
  scale_x_discrete(labels = trat_labels_herbivory) +
  scale_fill_manual(values = trat_colors, labels = trat_labels) +
  scale_y_continuous(labels = scales::percent, limits=c(0,0.80)) +
  theme_1

plot_herbivory_3 <- plot3
plot_herbivory_3

ls4 <- emmeans(model1, ~ Tratamiento | Cosecha | Origin, type="response")
cld4 <- as.data.frame(cld(ls4, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld4$.group <- str_trim(cld4$.group) # Remove space
cld4$Tratamiento <- factor(cld4$Tratamiento,levels = c("H2O", "H100", "T100"))
cld4 <- cld4 %>% arrange(Origin, Cosecha, Tratamiento)
cld4
pairs(ls4)

df_raw4 <- daño_filtrado %>%
  group_by(Cosecha, Tratamiento, Origin) %>%
  summarise(
    mean = mean(Incidencia, na.rm = TRUE),
    se   = std.error(Incidencia, na.rm = TRUE),
    .groups = "drop") %>% 
  left_join(cld4 %>% dplyr::select(Cosecha, Tratamiento, Origin, .group), by = c("Cosecha", "Tratamiento", "Origin"))

plot4 <- ggplot(df_raw4, aes(x = Cosecha, y = mean, fill = Tratamiento)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se,
                    ymax = mean+se),
                position = position_dodge(width = 0.8),
                width = 0.2,
                size = 0.5) +
  facet_wrap(~ Origin) +
  geom_text(aes(label = .group, y = mean+se + 0.05 * max(mean)), position = position_dodge(width = 0.8), size = 5) +
  labs(x = "Treatment",
       y = "Herbivory (%)",
       fill = "Treatment") +
  scale_x_discrete(labels = cosecha_labels) +
  scale_fill_manual(values = trat_colors, labels = trat_labels_herbivory) +
  theme_1

plot4

### TRICHOMES #####
###### LOAD DATA #####

tricomas <- read_excel("Tricomas.xlsx")
tricomas <- tricomas %>%
  mutate(
    Cosecha = as.factor(Cosecha),
    Origin = as.factor(Origin),
    Genotype = as.factor(Genotype),
    Gen_id = as.factor(Gen_id),
    Linea = as.factor(Linea),
    Linea_2 = as.factor(Linea_2),
    Tratamiento = factor(Tratamiento, levels = c("H2O", "H25", "T25", "H50", "T50", "H75", "T75", "H100", "T100")),
    Hoja_id = as.factor(Hoja_id),
    Photo = as.factor(Photo)
  )

head(tricomas)

tricomas_hoja_raw <- tricomas %>%
  group_by(Cosecha, Origin, Genotype, Gen_id, Linea, Linea_2, Tratamiento, Hoja_id) %>%
  summarise(Suma_Tricomas = sum(Tricomas, na.rm = TRUE),
            N_Photos = n_distinct(Photo),
            Tricomas_cm2 = Suma_Tricomas / (N_Photos * 20) * 100,
            .groups = "drop")

# Se mide el Area foliar individual en "tricomas_hoja_raw"
write_xlsx(tricomas_hoja_raw, "Tricomas_hoja_raw.xlsx")

tricomas_hoja <- read_excel("Tricomas_hoja.xlsx")
tricomas_hoja <- tricomas_hoja %>%
  mutate(
    Cosecha = as.factor(Cosecha),
    Origin = as.factor(Origin),
    Genotype = as.factor(Genotype),
    Gen_id = as.factor(Gen_id),
    Linea = as.factor(Linea),
    Linea_2 = as.factor(Linea_2),
    Tratamiento = factor(Tratamiento, levels = c("H2O", "H25", "T25", "H50", "T50", "H75", "T75", "H100", "T100")),
    Hoja_id = as.factor(Hoja_id)
  )

tricomas_hoja_filtrado <- tricomas_hoja %>%
  filter(Cosecha %in% c("1","2","3"),
         !Gen_id %in% c("I"), 
         !Tratamiento %in% c("H25", "T25", "H50", "T50", "H75", "T75")) %>%
  mutate(Cosecha = factor(Cosecha),
         Origin = factor(Origin),
         Genotype = factor(Genotype),
         Gen_id = factor(Gen_id),
         Linea = factor(Linea),
         Tratamiento = factor(Tratamiento, levels = c("H2O", "H100", "T100")))

# Se calcula el "Area_sum" (area foliar por hoja) para pasar la variable a Biomasa_SLA

area_planta_raw <- tricomas_hoja %>%
  group_by(Cosecha, Genotype, Linea, Tratamiento) %>%
  summarise(Area_sum = sum(Area_foliar_individual, na.rm = TRUE))

#write_xlsx(area_planta_raw, "Area_planta_raw.xlsx")

tricomas_domesticados <- tricomas_hoja_filtrado %>%
  filter(Origin=="Domesticated") %>%
  mutate(Origin = factor(Origin),
         Cultivar = factor(Genotype))

tricomas_silvestres <- tricomas_hoja_filtrado %>%
  filter(Origin=="Wild") %>%
  mutate(Origin = factor(Origin),
         Cultivar = factor(Genotype))

###### Trifoliate LA #####

model1 <- glmmTMB(log(Area_foliar_individual) ~ Cosecha * Origin * Tratamiento + (1|Gen_id) + (1|Linea),
                  family = gaussian(),
                  dispformula = ~ Origin,
                  data = tricomas_hoja_filtrado)

part1_model1 <- glmmTMB(log(Area_foliar_individual) ~ Cosecha * Origin * Tratamiento + (1|Linea),
                        family = gaussian(),
                        dispformula = ~ Origin,
                        data = tricomas_hoja_filtrado)

part2_model1 <- glmmTMB(log(Area_foliar_individual) ~ Cosecha * Origin * Tratamiento + (1|Gen_id),
                        family = gaussian(),
                        dispformula = ~ Origin,
                        data = tricomas_hoja_filtrado)

plot(simulateResiduals(model1))

summary(model1)
Anova(model1, type="II")
anova(part1_model1, model1)
anova(part2_model1, model1)

ls1 <- emmeans(model1, ~ Tratamiento | Cosecha, type="response")
cld1 <- as.data.frame(cld(ls1, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld1$.group <- str_trim(cld1$.group) # Remove space
cld1$Tratamiento <- factor(cld1$Tratamiento,levels = c("H2O", "H100", "T100"))
cld1 <- cld1 %>% arrange(Cosecha, Tratamiento)

df_raw1 <- tricomas_hoja_filtrado %>%
  group_by(Cosecha, Tratamiento) %>%
  summarise(
    mean = mean(Area_foliar_individual, na.rm = TRUE),
    se   = std.error(Area_foliar_individual, na.rm = TRUE),
    .groups = "drop") %>% 
  left_join(cld1 %>% dplyr::select(Cosecha, Tratamiento, .group), by = c("Cosecha", "Tratamiento"))

plot1 <- ggplot(df_raw1, aes(x = Cosecha, y = mean, fill = Tratamiento)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se,
                    ymax = mean+se),
                position = position_dodge(width = 0.8),
                width = 0.2,
                size = 0.5) +
  #facet_wrap(~ Cosecha, labeller = labeller(Cosecha = cosecha_labels)) +
  geom_text(aes(label = .group, y = mean+se + 0.05 * max(mean)), position = position_dodge(width = 0.8), size = 5) +
  labs(x = "Time (days after sowing)",
       y = expression(Trifoliate~leaf~area~(cm^{2})),
       fill = "Treatment") +
  scale_x_discrete(labels = cosecha_labels) +
  scale_fill_manual(values = trat_colors, labels = trat_labels) +
  theme_1

plot_leaf_area <- plot1

###### Trichome density #####

model1 <- glmmTMB(sqrt(Tricomas_cm2) ~ Cosecha * Origin * Tratamiento + (1|Gen_id) + (1|Linea),
                  family = gaussian(),
                  dispformula = ~ Origin,
                  data = tricomas_hoja_filtrado)

part1_model1 <- glmmTMB(sqrt(Tricomas_cm2) ~ Cosecha * Origin * Tratamiento + (1|Linea),
                        family = gaussian(),
                        dispformula = ~ Origin,
                        data = tricomas_hoja_filtrado)

part2_model1 <- glmmTMB(sqrt(Tricomas_cm2) ~ Cosecha * Origin * Tratamiento + (1|Gen_id),
                        family = gaussian(),
                        dispformula = ~ Origin,
                        data = tricomas_hoja_filtrado)

summary(model1)
Anova(model1, type="III")
anova(part1_model1, model1)
anova(part2_model1, model1)

ls1 <- emmeans(model1, ~ Tratamiento | Origin, type="response")
cld1 <- as.data.frame(cld(ls1, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld1$.group <- str_trim(cld1$.group) # Remove space
cld1$Tratamiento <- factor(cld1$Tratamiento,levels = c("H2O", "H100", "T100"))
cld1 <- cld1 %>% arrange(Origin, Tratamiento)

df_raw1 <- tricomas_hoja_filtrado %>%
  group_by(Tratamiento, Origin) %>%
  summarise(
    mean = mean(Tricomas_cm2, na.rm = TRUE),
    se   = std.error(Tricomas_cm2, na.rm = TRUE),
    .groups = "drop") %>% 
  left_join(cld1 %>% dplyr::select(Tratamiento, Origin, .group), by = c("Tratamiento", "Origin"))

plot1 <- ggplot(df_raw1, aes(x = Origin, y = mean, fill = Tratamiento)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se,
                    ymax = mean+se),
                position = position_dodge(width = 0.8),
                width = 0.2,
                size = 0.5) +
  #facet_wrap(~ Origin) +
  geom_text(aes(label = .group, y = mean+se + 0.15 * max(mean)), position = position_dodge(width = 0.8), size = 5) +
  labs(x = "Status",
       y = expression(TD~(no.~cm^{-2})),
       fill = "Treatment") +
  scale_x_discrete(labels = trat_labels) +
  scale_fill_manual(values = trat_colors, labels = trat_labels) +
  theme_1

plot_trichomes <- plot1

### NECTARIOS #####

nectarios <- read_excel("Nectarios.xlsx")

nectarios_filtrados <- nectarios %>%
  subset(Gen_id %in% c("A","C","H")) %>%
  mutate(Gen_id = as.factor(Gen_id),
         Origin = as.factor(Origin),
         Linea = as.factor(Linea),
         Linea_2 = as.factor(Linea_2),
         Hoja = as.factor(Hoja),
         Hoja_2 = as.factor(Hoja_2),
         Foliolo = as.factor(Foliolo),
         Tratamiento = factor(Tratamiento, levels = c("H2O", "H100", "T100")))

model1 <- glmmTMB(Grados_brix ~ Origin * Tratamiento + (1|Linea/Hoja_2),
                  family = gaussian(),
                  data = nectarios_filtrados)

part1_model1 <- glmmTMB(Grados_brix ~ Origin * Tratamiento,
                        family = gaussian(),
                        data = nectarios_filtrados)

summary(model1)
Anova(model1, type = "III")
anova(part1_model1, model1)

ls1 <- emmeans(model1, ~ Tratamiento | Origin, type="response")
cld1 <- as.data.frame(cld(ls1, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld1$.group <- str_trim(cld1$.group) # Remove space
cld1$Tratamiento <- factor(cld1$Tratamiento,levels = c("H2O", "H100", "T100"))
cld1 <- cld1 %>% arrange(Origin, Tratamiento)

df_raw1 <- nectarios_filtrados %>%
  group_by(Tratamiento, Origin) %>%
  summarise(
    mean = mean(Grados_brix, na.rm = TRUE),
    se   = std.error(Grados_brix, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  left_join(cld1 %>% dplyr::select(Tratamiento, Origin, .group), by = c("Tratamiento", "Origin"))

plot_1 <- ggplot(df_raw1, aes(x=Origin, y=mean, fill=Tratamiento)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.5, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position = position_dodge(width = 0.6),
                width = 0.2) +
  geom_text(aes(label = .group, y = mean+se + 0.15 * max(mean)), position = position_dodge(width = 0.6), size = 5) +
  labs(x = "Status",
       y = "EFN (°Bx)",
       fill = "Treatment") +
  #scale_y_continuous(limits = c(0,250)) +
  scale_fill_manual(values = trat_colors, labels = trat_labels) +
  theme_1

plot_nectarios <- plot_1

### YIELD #####

rendimiento <- read_excel("Rendimiento.xlsx")
rendimiento <- rendimiento %>%
  mutate(Cultivar = as.factor(Cultivar),
         Origin = as.factor(Origin),
         Line = as.factor(Line),
         Treatment = factor(Treatment, levels = c("H2O", "H100", "T100")),
         Yield_total = as.numeric(Yield_total),
         Yield_plant = as.numeric(Yield_plant))

rendimiento_filtrado <- rendimiento %>%
  filter(Cultivar %in% c("D1","D5")) %>%
  mutate(Cultivar = as.factor(Cultivar),
         Line = as.factor(Line))

model1 <- glmmTMB(Yield_plant ~ Treatment * Cultivar + (1|Line),
                  family= Gamma(link="log"),
                  data = rendimiento_filtrado)

part1_model1 <- glmmTMB(Yield_plant ~ Treatment * Cultivar,
                  family= Gamma(link="log"),
                  data = rendimiento_filtrado)

summary(model1)
Anova(model1, type ="III")
anova(part1_model1, model1)

ls1 <- emmeans(model1, ~ Cultivar, type="response")
cld1 <- as.data.frame(cld(ls1, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld1$.group <- str_trim(cld1$.group) # Remove space
cld1 <- cld1 %>% arrange(Cultivar)

df_raw1 <- rendimiento_filtrado %>%
  group_by(Cultivar) %>%
  summarise(
    mean = mean(Yield_plant, na.rm = TRUE),
    se   = std.error(Yield_plant, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  left_join(cld1 %>% dplyr::select(Cultivar, .group), by = c("Cultivar"))

plot_1 <- ggplot(df_raw1, aes(x=Cultivar, y=mean)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.5, colour = "black", fill="gray95", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position = position_dodge(width = 0.6),
                width = 0.2) +
  geom_text(aes(label = .group, y = mean+se + 0.05 * max(mean)), position = position_dodge(width = 0.6), size = 5) +
  labs(x = "Accession",
       y = expression("Seeds plant"^{-1}~"(g)")) +
  scale_y_continuous(limits = c(0,250)) +
  theme_1

plot_yield_1 <- plot_1

ls2 <- emmeans(model1, ~ Treatment, type="response")
cld2 <- as.data.frame(cld(ls2, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld2$.group <- str_trim(cld2$.group) # Remove space
cld2$Tratamiento <- factor(cld2$Treatment,levels = c("H2O", "H100", "T100"))
cld2 <- cld2 %>% arrange(Treatment)

df_raw2 <- rendimiento_filtrado %>%
  group_by(Treatment) %>%
  summarise(
    mean = mean(Yield_plant, na.rm = TRUE),
    se   = std.error(Yield_plant, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  left_join(cld2 %>% dplyr::select(Treatment, .group), by = c("Treatment"))

plot_2 <- ggplot(df_raw2, aes(x=Treatment, y=mean, fill=Treatment)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.5, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position = position_dodge(width = 0.6),
                width = 0.2) +
  geom_text(aes(label = .group, y = mean+se + 0.05 * max(mean)), position = position_dodge(width = 0.6), size = 5) +
  labs(x = "Accession",
       y = expression("Seeds plant"^{-1}~"(g)"),
       fill = "Treatment") +
  scale_y_continuous(limits = c(0,250)) +
  scale_x_discrete(labels = trat_labels) +
  scale_fill_manual(values = trat_colors, labels = trat_labels) +
  theme_1

plot_yield_2 <- plot_2

ls3 <- emmeans(model1, ~ Treatment | Cultivar, type="response")
cld3 <- as.data.frame(cld(ls3, alpha=0.05, Letters=letters, type="response", decreasing=T))
cld3$.group <- str_trim(cld3$.group) # Remove space
cld3 <- cld3 %>% arrange(Cultivar, Treatment)

df_raw3 <- rendimiento_filtrado %>%
  group_by(Treatment, Cultivar) %>%
  summarise(
    mean = mean(Yield_plant, na.rm = TRUE),
    se   = std.error(Yield_plant, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  left_join(cld3 %>% dplyr::select(Treatment, Cultivar, .group), by = c("Treatment", "Cultivar"))

plot_3 <- ggplot(df_raw3, aes(x=Cultivar, y=mean, fill=Treatment)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.5, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position = position_dodge(width = 0.6),
                width = 0.2) +
  geom_text(aes(label = .group, y = mean+se + 0.05 * max(mean)), position = position_dodge(width = 0.6), size = 5) +
  labs(x = "Accession",
       y = expression("Seeds plant"^{-1}~"(g)"),
       fill = "Treatment") +
  scale_y_continuous(limits = c(0,250)) +
  scale_fill_manual(values = trat_colors, labels = trat_labels) +
  theme_1

plot_yield_3 <- plot_3

### Figures EXPORT #####

p1 <- plot_rgr_diameter + theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.x = element_blank())

p2 <- plot_rgr_height + theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.x = element_blank(),
  strip.text.x = element_blank()) +
  theme_legend_vertical

p3 <- plot_rgr_leaf_number + theme(
  strip.text.x = element_blank())

figura_1 <- p1 / p2 / p3 +
  plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)"))) &
  theme(plot.tag = element_text(size = 18))
figura_1

ggsave(
  "plots/Figure 1.svg",
  figura_1,
  width = 174,
  height = 234,
  units = "mm"
)

p1 <- plot_leaf_area
p2 <- plot_nectarios + theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.x = element_blank())
p3 <- plot_trichomes + theme(
  strip.text.x = element_blank())
p4 <- plot_toughness +
  theme_legend_left

figura_2 <- ((p1 + p2/p3) + plot_layout(widths = c(1.4, 1))) / p4 +
  plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)", "(d)"))) &
  theme(plot.tag = element_text(size = 18))
figura_2

ggsave(
  "plots/Figure 2.svg",
  figura_2,
  width = 174,
  height = 234,
  units = "mm"
)

figura_3 <- plot_yield_3 +
  theme_legend_vertical_left
figura_3

ggsave(
  "plots/Figure 3.svg",
  figura_3,
  width = 94,
  height = 124,
  units = "mm"
)

p1 <- plot_herbivory_1

p2 <- plot_herbivory_2 + theme(
  axis.line.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.title.y = element_blank())

p3 <- plot_herbivory_3 + theme(
  axis.line.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.title.y = element_blank())

figura_4 <- p1 + p2 + p3 + plot_layout(widths = c(1, 0.7, 1)) +
  plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)"))) &
  theme(plot.tag = element_text(size = 18))
figura_4

ggsave(
  "plots/Figure 4.svg",
  figura_4,
  width = 210,
  height = 126,
  units = "mm"
)





p1 <- plot_height + theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.x = element_blank()) +
  theme_legend_vertical_left
p2 <- plot_leaf_area
p3 <- plot_diameter
p4 <- plot_leaf_number

figure_2 <- ((p1 / p2) | (p3 / p4)) + plot_layout(widths = c(1, 2)) +
  plot_annotation(tag_levels = list(c("(a)", "(c)", "(b)", "(d)"))) &
  theme(plot.tag = element_text(size = 18))
figure_2

ggsave(
  "plots/Figure 2.png",
  figure_2,
  width = 1700,
  height = 1800,
  units = "px",
  dpi = 200
)

p1 <- plot_toughness +
  theme_legend_vertical_left
p2 <- plot_nectarios + theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.x = element_blank())
p3 <- plot_trichomes + theme(
  strip.text.x = element_blank())

figura_4 <- p1 + (p2/p3) + plot_layout(widths = c(2.4, 1)) +
  plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)"))) &
  theme(plot.tag = element_text(size = 18))
figura_4

ggsave(
  "plots/Figure 4.png",
  figura_4,
  width = 2100,
  height = 1400,
  units = "px",
  dpi = 200
)

#plot_thickness
#plot_thickness2

### DATA ORGANIZATION for SEM #####

summary(biomasa_filtrado)
str(biomasa_filtrado)
#write_xlsx(biomasa_filtrado, "biomasa_filtrado.xlsx")

summary(crecimiento_filtrado)
str(crecimiento_filtrado)
#write_xlsx(crecimiento_filtrado, "crecimiento_filtrado.xlsx")

summary(daño_filtrado)
str(daño_filtrado)
#write_xlsx(daño_filtrado, "daño_filtrado.xlsx")

summary(defensas_filtrado)
str(defensas_filtrado)
#write_xlsx(defensas_filtrado, "defensas_filtrado.xlsx")

summary(tricomas_hoja_filtrado)
str(tricomas_hoja_filtrado)
#write_xlsx(tricomas_hoja_filtrado, "tricomas_hoja_filtrado.xlsx")

crecimiento_master <- crecimiento_filtrado %>%
  group_by(Cosecha, Genotype, Origin, Linea, Tratamiento) %>%
  summarise(Diametro = mean(Diametro, na.rm = TRUE),
            No_hojas = mean(No_hojas, na.rm = TRUE),
            Altura = mean(Altura, na.rm = TRUE),
            RGR_Diametro = mean(RGR_Diametro, na.rm = TRUE),
            RGR_No_hojas = mean(RGR_No_hojas, na.rm = TRUE),
            RGR_Altura = mean(RGR_Altura, na.rm = TRUE),
            LPR = mean(LPR, na.rm = TRUE),
            .groups = "drop")

defensas_master <- defensas_filtrado %>%
  group_by(Cosecha, Genotype, Origin, Linea, Tratamiento) %>%
  summarise(Dureza = mean(Dureza, na.rm = TRUE),
            Grosor = mean(Grosor, na.rm = TRUE),
            .groups = "drop")

tricomas_master <- tricomas_hoja_filtrado %>%
  group_by(Cosecha, Genotype, Origin, Linea, Tratamiento) %>%
  summarise(Tricomas = mean(Tricomas_cm2, na.rm = TRUE),
            Area_ind = mean(Area_foliar_individual, na.rm = TRUE),
            .groups = "drop")

master_df <- biomasa_filtrado %>%
  full_join(daño_filtrado, by = c("Cosecha","Genotype","Origin","Linea","Tratamiento")) %>%
  full_join(crecimiento_master, by = c("Cosecha","Genotype","Origin","Linea","Tratamiento")) %>%
  full_join(defensas_master, by = c("Cosecha","Genotype","Origin","Linea","Tratamiento")) %>%
  full_join(tricomas_master, by = c("Cosecha","Genotype","Origin","Linea","Tratamiento"))
write_xlsx(master_df, "master_df.xlsx")

vars <- master_df %>%
  dplyr::select(Altura, Diametro, No_hojas,
                SLA_2, LMA_2, Area_ind,
                RGR_Altura, RGR_Diametro, RGR_No_hojas, LPR,
                Tricomas, Dureza, Grosor, 
                Damaged_obs, Intact_obs, Incidencia)

M = cor(vars, use = "pairwise.complete.obs")
print(M)
corrplot(M, method = 'circle', type = 'lower', addCoef.col ='black', number.cex = 0.6, diag=FALSE)


### Model selection tools #####

model1 <- glmmTMB(Yield_plant ~ Treatment * Cultivar + (1|Line),
                  family= Gamma(link="log"),
                  data = rendimiento_filtrado)

model2 <- glmmTMB(Yield_plant ~ Treatment * Cultivar + (1|Line),
                  family= Gamma(link="log"),
                  data = rendimiento_filtrado,
                  dispformula = ~ Cultivar)

model3 <- glmmTMB(Yield_plant ~ Treatment * Cultivar + (1|Line),
                  family= Gamma(link="log"),
                  data = rendimiento_filtrado,
                  dispformula = ~ Treatment)

model4 <- glmmTMB(Yield_plant ~ Treatment * Cultivar + (1|Line),
                  family= Gamma(link="log"),
                  data = rendimiento_filtrado,
                  dispformula = ~ Cultivar + Treatment)

model5 <- glmmTMB(Yield_plant ~ Treatment * Cultivar + (1|Line),
                  family= Gamma(link="log"),
                  data = rendimiento_filtrado,
                  dispformula = ~ Cultivar * Treatment)

model1 <- glmmTMB(No_hojas_tr ~ Cosecha * Origin * Tratamiento + (1|Gen_id) + (1|Linea),
                  family = gaussian(),
                  data = crecimiento_filtrado)

model2 <- glmmTMB(No_hojas_tr ~ Cosecha * Origin * Tratamiento + (1|Gen_id) + (1|Linea),
                  family = gaussian(),
                  data = crecimiento_filtrado,
                  dispformula = ~ Cosecha)

model3 <- glmmTMB(No_hojas_tr ~ Cosecha * Origin * Tratamiento + (1|Gen_id) + (1|Linea),
                  family = gaussian(),
                  data = crecimiento_filtrado,
                  dispformula = ~ Origin)

model4 <- glmmTMB(No_hojas_tr ~ Cosecha * Origin * Tratamiento + (1|Gen_id) + (1|Linea),
                  family = gaussian(),
                  data = crecimiento_filtrado,
                  dispformula = ~ Tratamiento)

model5 <- glmmTMB(No_hojas_tr ~ Cosecha * Origin * Tratamiento + (1|Gen_id) + (1|Linea),
                  family = gaussian(),
                  data = crecimiento_filtrado,
                  dispformula = ~ Tratamiento + Origin)

model6 <- glmmTMB(No_hojas_tr ~ Cosecha * Origin * Tratamiento + (1|Gen_id) + (1|Linea),
                  family = gaussian(),
                  data = crecimiento_filtrado,
                  dispformula = ~ Tratamiento + Cosecha)

model7 <- glmmTMB(No_hojas_tr ~ Cosecha * Origin * Tratamiento + (1|Gen_id) + (1|Linea),
                  family = gaussian(),
                  data = crecimiento_filtrado,
                  dispformula = ~ Cosecha + Origin)

model8 <- glmmTMB(No_hojas_tr ~ Cosecha * Origin * Tratamiento + (1|Gen_id) + (1|Linea),
                  family = gaussian(),
                  data = crecimiento_filtrado,
                  dispformula = ~ Cosecha + Tratamiento + Origin)

model9 <- glmmTMB(No_hojas_tr ~ Cosecha * Origin * Tratamiento + (1|Gen_id) + (1|Linea),
                  family = gaussian(),
                  data = crecimiento_filtrado,
                  dispformula = ~ Tratamiento * Origin)

model10 <- glmmTMB(No_hojas_tr ~ Cosecha * Origin * Tratamiento + (1|Gen_id) + (1|Linea),
                   family = gaussian(),
                   data = crecimiento_filtrado,
                   dispformula = ~ Tratamiento * Cosecha)

model11 <- glmmTMB(No_hojas_tr ~ Cosecha * Origin * Tratamiento + (1|Gen_id) + (1|Linea),
                   family = gaussian(),
                   data = crecimiento_filtrado,
                   dispformula = ~ Cosecha * Origin)

model12 <- glmmTMB(No_hojas_tr ~ Cosecha * Origin * Tratamiento + (1|Gen_id) + (1|Linea),
                   family = gaussian(),
                   data = crecimiento_filtrado,
                   dispformula = ~ Cosecha * Tratamiento * Origin)

performance::compare_performance(
  model1, model2, model3, model4, model5, model6, model7, model8, model9, model10, model11, model12
)

plot(simulateResiduals(model5))

Anova(model1, type = "III")
Anova(model2, type = "III")
Anova(model3, type = "III")
Anova(model4, type = "III")
Anova(model5, type = "III")
Anova(model6, type = "III")
Anova(model7, type = "III")#
Anova(model8, type = "III")#
Anova(model9, type = "III")
Anova(model10, type = "III")
Anova(model11, type = "III")#
Anova(model12, type = "III")

Anova(model1, type = "II")
Anova(model2, type = "II")
Anova(model3, type = "II")
Anova(model4, type = "II")
Anova(model5, type = "II")
Anova(model6, type = "II")
Anova(model7, type = "II")
Anova(model8, type = "II")
Anova(model9, type = "II")
Anova(model10, type = "II")
Anova(model11, type = "II")
Anova(model12, type = "II")