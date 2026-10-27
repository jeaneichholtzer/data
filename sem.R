# TITLE: Leaf and stem extracts of Mexican sunflower differentially affect growth, defense and herbivory in domesticated and wild lima bean
# AUTHORS: Jean Eichholtzer, Roberto Rafael Ruiz-Santiago, Alejandra del Socorro González-Moreno, Esaú Ruiz-Sánchez, Germán Giácoman Vallejos, Horacio Salomón Ballina-Gómez*
# CORREPONDING: horacio.bg@conkal.tecnm.mx

# Set Working Directory

#setwd("yourpath/data")

### I Load packages ############################################################################################################################################################################################################

load_packages <- function(pkgs)
{missing_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(missing_pkgs) > 0)
{install.packages(missing_pkgs, dependencies = TRUE)}
invisible(lapply(pkgs, library, character.only = TRUE))}

load_packages(c("corrplot",
                "DiagrammeR",
                "mice",
                "piecewiseSEM",
                "readxl", 
                "tidyverse"))

set.seed(123)

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
         Accession_id = as.factor(Accession_id),
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

### II Calculations ############################################################################################################################################################################################################

###### 1 SLA #####

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
  mutate(LMA = ((mean_prev_biomass / mean_prev_area) + (Leaf_biomass_sample / Leaf_area_sample)) / 2)

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
  mutate(LMA = ((mean_prev_biomass / mean_prev_area) +  (Leaf_biomass_sample / Leaf_area_sample)) / 2)

# Merge data

sla <- bind_rows(area_biomass %>% filter(Sampling_time == "0"),
                 full_join(sla_s1, lma_s1),
                 full_join(sla_s23, lma_s23)) %>%
  filter(Sampling_time %in% c("1","2","3")) %>%
  arrange(Accession_id, Row, Treatment, Sampling_time)

###### 2 RGR #####

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
            mean_log_prev_length   = mean(log(Length), na.rm = TRUE),
            date_prev = first(Date),
            .groups = "drop") %>%
  mutate(Sampling_time = as.character(as.numeric(as.character(Sampling_time)) + 1))

# RGR calculations for Sampling_time = 1

rgr_s1 <- growth %>%
  filter(Sampling_time == "1") %>%
  left_join(means_sampling_0,
            by = "Accession_id") %>%
  mutate(days_diff = as.numeric(Date - date_prev),
         RGR_Diameter = (log(Diameter) - mean_log_prev_diameter) / days_diff,
         RGR_No_leaves = (log(No_leaves) - mean_log_prev_no_leaves) / days_diff,
         RGR_Length = (log(Length) - mean_log_prev_length) / days_diff)

# RGR calculations for Sampling_time 2 and 3

rgr_s23 <- growth %>%
  filter(Sampling_time %in% c("2", "3")) %>%
  left_join(means_previous_sampling,
            by = c("Accession_id", "Treatment", "Sampling_time")) %>%
  mutate(days_diff = as.numeric(Date - date_prev),
         RGR_Diameter = (log(Diameter) - mean_log_prev_diameter) / days_diff,
         RGR_No_leaves = (log(No_leaves) - mean_log_prev_no_leaves) / days_diff,
         RGR_Length = (log(Length) - mean_log_prev_length) / days_diff)

# Merge data and apply filter

growth_rgr <- bind_rows(growth %>% filter(Sampling_time == "0"), rgr_s1, rgr_s23) %>%
  filter(Sampling_time %in% c("1","2","3")) %>%
  mutate(Sampling_time = as.factor(Sampling_time)) %>%
  arrange(Accession_id, Row, Treatment, Sampling_time)

### III Data organization for SEM ############################################################################################################################################################################################################

growth_sem <- growth_rgr %>%
  group_by(Sampling_time, Accession_id, Dom_status, Row, Treatment) %>%
  summarise(Diameter = mean(Diameter, na.rm = TRUE),
            No_leaves = mean(No_leaves, na.rm = TRUE),
            Length = mean(Length, na.rm = TRUE),
            RGR_Diameter = mean(RGR_Diameter, na.rm = TRUE),
            RGR_No_leaves = mean(RGR_No_leaves, na.rm = TRUE),
            RGR_Length = mean(RGR_Length, na.rm = TRUE),
            .groups = "drop")

defenses_sem <- defenses %>%
  group_by(Sampling_time, Accession_id, Dom_status, Row, Treatment) %>%
  summarise(Hardness = mean(Hardness, na.rm = TRUE),
            Thickness = mean(Thickness, na.rm = TRUE),
            .groups = "drop")

trichomes_sem <- trichomes %>%
  group_by(Sampling_time, Accession_id, Dom_status, Row, Treatment) %>%
  summarise(Trichomes_cm2 = mean(Trichomes_cm2, na.rm = TRUE),
            Trifoliate_leaf_area = mean(Trifoliate_leaf_area, na.rm = TRUE),
            .groups = "drop")

master_df <- sla %>%
  full_join(herbivory, by = c("Sampling_time","Accession_id","Dom_status","Row","Treatment")) %>%
  full_join(growth_sem, by = c("Sampling_time","Accession_id","Dom_status","Row","Treatment")) %>%
  full_join(defenses_sem, by = c("Sampling_time","Accession_id","Dom_status","Row","Treatment")) %>%
  full_join(trichomes_sem, by = c("Sampling_time","Accession_id","Dom_status","Row","Treatment"))

# Select variables for SEM environment

sem_df <- master_df %>%
  dplyr::select(Sampling_time, Accession_id, Dom_status, Row, Treatment,
                # Size
                Diameter, Length, No_leaves,
                # Leaf economics
                SLA, Trifoliate_leaf_area,
                # Physical resistance
                Hardness, Thickness, Trichomes_cm2,
                # Growth
                RGR_Diameter, RGR_Length, RGR_No_leaves,
                # Herbivory
                Herbivory_incidence, Damaged_obs, Intact_obs)

# Check missing values
md.pattern(sem_df)

# Check for unusual correlations among variables
M = cor(sem_df[,c(6:19)], use = "pairwise.complete.obs")
corrplot(M, method = 'circle', type = 'lower', addCoef.col ='black', number.cex = 0.6, diag=FALSE)

# Imputation MICE of variables used in SEM
vars_imp <- c("Diameter", "Length", "No_leaves", "SLA", "Trifoliate_leaf_area",
              "Hardness", "Thickness", "Trichomes_cm2", "RGR_Diameter", "RGR_Length", "RGR_No_leaves")
pred <- make.predictorMatrix(sem_df)
pred[,] <- 0
pred[vars_imp, c("Sampling_time", "Accession_id", "Treatment",
                 "Diameter", "Length", "No_leaves", "SLA", "Trifoliate_leaf_area", 
                 "Hardness", "Thickness", "Trichomes_cm2", "RGR_Diameter", "RGR_Length", "RGR_No_leaves",
                 "Herbivory_incidence", "Damaged_obs", "Intact_obs")] <- 1
diag(pred) <- 0
imp <- mice(sem_df, predictorMatrix = pred, method = "pmm", m = 20, maxit = 50, seed = 123)

avg_matrix <- Reduce("+",lapply(1:20, function(i) complete(imp, i)[, vars_imp])) / 20
sem_df_final <- sem_df
sem_df_final[, vars_imp] <- avg_matrix

# Scale all continuous predictors (mean = 0, sd = 1) to get path coefficients directly comparable
sem_df_final <- sem_df_final %>% 
  mutate(DAS = c(-20,0,20)[as.numeric(Sampling_time)], # Centered only
         N_leaves = Damaged_obs + Intact_obs,
         across(c(DAS,
                  Hardness, Thickness, Trichomes_cm2,
                  RGR_Diameter, RGR_Length, RGR_No_leaves,
                  Herbivory_incidence),
                ~ scale(.x)[,1],
                .names = "{.col}_s"))

# Check for unusual correlations among variables
M = cor(sem_df_final[,c(6:19,21)], use = "pairwise.complete.obs")
corrplot(M, method = 'circle', type = 'lower', addCoef.col ='black', number.cex = 0.6, diag=FALSE)

### IV SEM construction ############################################################################################################################################################################################################

# SEM function
fit_sem_treatment <- function(df, trat_label) {
  # Physical resistance
  mod_Hardness <- glmmTMB(Hardness ~ Herbivory_incidence + DAS + (1|Accession_id), family = gaussian(), data = df)
  mod_Thickness <- glmmTMB(Thickness ~ Herbivory_incidence + DAS + (1|Accession_id), family = gaussian(), data = df)
  mod_Trichomes_cm2 <- glmmTMB(Trichomes_cm2 ~ Herbivory_incidence + DAS + (1|Accession_id), family = gaussian(), data = df)
  # Tolerance/compensatory growth
  mod_RGR_No_leaves <- glmmTMB(RGR_No_leaves ~ Herbivory_incidence + DAS + (1|Accession_id), family = gaussian(), data = df)
  mod_RGR_Length <- glmmTMB(RGR_Length ~ Herbivory_incidence + DAS + (1|Accession_id), family = gaussian(), data = df)
  mod_RGR_Diameter <- glmmTMB(RGR_Diameter ~ Herbivory_incidence + DAS + (1|Accession_id), family = gaussian(), data = df)
  # SEM assamblage
  sem <- psem(mod_Hardness, mod_Thickness, mod_Trichomes_cm2, 
              mod_RGR_No_leaves, mod_RGR_Length, mod_RGR_Diameter,
              # Residual correlations within constructs
              RGR_Diameter %~~% Thickness,
              RGR_Diameter %~~% Hardness,
              RGR_Diameter %~~% Trichomes_cm2,
              RGR_Length %~~% Thickness,
              RGR_Length %~~% Hardness,
              RGR_Length %~~% Trichomes_cm2,
              RGR_No_leaves %~~% Thickness,
              RGR_No_leaves %~~% Hardness,
              RGR_No_leaves %~~% Trichomes_cm2,
              # New relationships suggested by d-separation
              Hardness %~~% Thickness,
              Trichomes_cm2 %~~% Thickness,
              RGR_Diameter %~~% RGR_No_leaves,
              RGR_Diameter %~~% RGR_Length,
              RGR_Length %~~% RGR_No_leaves
  )
  return(sem)}

treatment_levels <- c("Control", "Leaf extract", "Stem extract")
status_levels <- c("Domesticated", "Wild")
group_names <- as.vector(outer(treatment_levels, status_levels, paste, sep = "_"))

# Split into named list of subsets
sem_subsets <- setNames(lapply(group_names, function(g){
  parts <- strsplit(g, "_")[[1]]
  sem_df_final %>%
    filter(Treatment == parts[1],
           Dom_status == parts[2])}),
  group_names)

# Fit SEM on each subset
sem_results <- lapply(group_names, function(t) {
  tryCatch(
    fit_sem_treatment(sem_subsets[[t]], t),
    error = function(e) {
      cat("✗ ERROR in", t, ":", conditionMessage(e), "\n")
      return(NULL)
    }
  )
})

names(sem_results) <- group_names

extract_sem_paths <- function(sem_model){
  tab <- coefs(sem_model, standardize = "scale") %>% as.data.frame()
  tab$Response <- case_when(grepl("^cbind\\(", tab$Response) ~ "Herbivory",
                            tab$Response == "Herbivory_incidence" ~ "Herbivory",
                            TRUE ~ tab$Response)
  tab %>% dplyr::select(Response, Predictor, Beta = Std.Estimate, Estimate, P.Value)}

extract_sem_r2 <- function(sem_model) {
  r2 <- rsquared(sem_model) %>% as.data.frame()
  r2 %>% dplyr::select(Response, Marginal) %>%
    dplyr::mutate(R2 = sprintf("%.2f", Marginal))}

node_label <- function(title, r2) {paste0("<",title,
                                          "<BR/><FONT POINT-SIZE='14'><I>R² =</I>",
                                          r2,"</FONT>",">")}

sig_code <- function(p){
  if(is.na(p)) return("")
  if(p < 0.001) return("<FONT POINT-SIZE='10'><SUP>***</SUP></FONT>")
  if(p < 0.01)  return("<FONT POINT-SIZE='10'><SUP>**</SUP></FONT>")
  if(p < 0.05)  return("<FONT POINT-SIZE='10'><SUP>*</SUP></FONT>")
  #if(p < 0.10)  return("<FONT POINT-SIZE='10'><SUP>†</SUP></FONT>")
  return("")
}

make_label <- function(beta, p){
  stars <- sig_code(p)
  if(!is.na(p) && p < 0.05){paste0("<<BR/><B><I>β = </I>", sprintf("%.2f", beta), stars, "</B>>")} 
  else {paste0("<<BR/><I>β = </I>", sprintf("%.2f", beta), stars,">")}}

edge_style <- function(tab, pred, resp){
  row <- tab[tab$Predictor == pred & tab$Response == resp,]
  if(nrow(row) == 0) return("")
  beta <- row$Beta[1]
  pval <- row$P.Value[1]
  lab <- make_label(beta, pval)
  sig <- pval < 0.05
  col <- if(sig){ifelse(beta > 0, "forestgreen", "firebrick")} else {"gray60"}
  lty <- ifelse(sig, "solid", "dashed")
  w <- round(1 + 4 * abs(beta), 2)
  paste0(
    "[label=", lab,
    " color='", col,
    "' fontcolor='black'",
    " fontname='Arial'",
    " style='", lty,
    "' penwidth=", w,
    "]"
  )}

corr_style <- function(){
  paste0(
    "[dir=both",
    " color='gray60'",
    " style='dashed'",
    " constraint=false",
    "]"
  )
}

plot_sem <- function(sem_model){
  
  tab <- extract_sem_paths(sem_model)
  r2  <- extract_sem_r2(sem_model)
  
  # Helper to retrieve marginal R2
  get_r2 <- function(response) {
    x <- r2$R2[r2$Response == response]
    if(length(x) == 0 || is.na(x)) return("NA")
    x[1]
  }
  
  grViz(paste0(
    "
digraph {

graph[
layout = dot
rankdir = LR
splines = true
labelloc='t'
fontsize=24
fontname='Arial'
]

node[
shape = box
style = 'filled'
fontsize=18
fontname='Arial'
]

Herbivory [label='Herbivory', fillcolor='#FFFFFF']

Thickness [
  label=", node_label("Thickness", get_r2("Thickness")), "
  fillcolor='#F2F2F2'
]

Hardness [
  label=", node_label("Toughness", get_r2("Hardness")), "
  fillcolor='#F2F2F2'
]

Trichomes [
  label=", node_label("TD", get_r2("Trichomes_cm2")), "
  fillcolor='#F2F2F2'
]

RGR_D [
  label=<RGR<FONT POINT-SIZE='8'><I>D</I></FONT>
  <BR/><FONT POINT-SIZE='14'><I>R² =</I>",
get_r2("RGR_Diameter"),
"</FONT>>,
  fillcolor='#F2F2F2'
]

RGR_H [
  label=<RGR<FONT POINT-SIZE='8'><I>H</I></FONT>
  <BR/><FONT POINT-SIZE='14'><I>R² =</I>",
get_r2("RGR_Length"),
"</FONT>>,
  fillcolor='#F2F2F2'
]

RGR_N [
  label=<RGR<FONT POINT-SIZE='8'><I>LN</I></FONT>
  <BR/><FONT POINT-SIZE='14'><I>R² =</I>",
get_r2("RGR_No_leaves"),
"</FONT>>,
  fillcolor='#F2F2F2'
]

Herbivory -> Thickness
", edge_style(tab, 'Herbivory_incidence', 'Thickness'), "

Herbivory -> Hardness
", edge_style(tab, 'Herbivory_incidence', 'Hardness'), "

Herbivory -> Trichomes
", edge_style(tab, 'Herbivory_incidence', 'Trichomes_cm2'), "

Herbivory -> RGR_D
", edge_style(tab, 'Herbivory_incidence', 'RGR_Diameter'), "

Herbivory -> RGR_H
", edge_style(tab, 'Herbivory_incidence', 'RGR_Length'), "

Herbivory -> RGR_N
", edge_style(tab, 'Herbivory_incidence', 'RGR_No_leaves'), "

}
"))}

### V Diagram generation ############################################################################################################################################################################################################

domesticated_control <- plot_sem(sem_results[["Control_Domesticated"]])
domesticated_h10 <- plot_sem(sem_results[["Leaf extract_Domesticated"]])
domesticated_t10 <- plot_sem(sem_results[["Stem extract_Domesticated"]])
wild_control <- plot_sem(sem_results[["Control_Wild"]])
wild_h10 <- plot_sem(sem_results[["Leaf extract_Wild"]])
wild_t10 <- plot_sem(sem_results[["Stem extract_Wild"]])

# Diagram visualization

domesticated_control
domesticated_h10
domesticated_t10
wild_control
wild_h10
wild_t10

# Fisher's C statistics

fisherC(sem_results[["Control_Domesticated"]])
fisherC(sem_results[["Leaf extract_Domesticated"]])
fisherC(sem_results[["Stem extract_Domesticated"]])
fisherC(sem_results[["Control_Wild"]])
fisherC(sem_results[["Leaf extract_Wild"]])
fisherC(sem_results[["Stem extract_Wild"]])

# Tests of directed separation

dSep(sem_results[["Control_Domesticated"]])
dSep(sem_results[["Leaf extract_Domesticated"]])
dSep(sem_results[["Stem extract_Domesticated"]])
dSep(sem_results[["Control_Wild"]])
dSep(sem_results[["Leaf extract_Wild"]])
dSep(sem_results[["Stem extract_Wild"]])

# Summary of SEM results

summary(sem_results[["Control_Domesticated"]])
summary(sem_results[["Leaf extract_Domesticated"]])
summary(sem_results[["Stem extract_Domesticated"]])
summary(sem_results[["Control_Wild"]])
summary(sem_results[["Leaf extract_Wild"]])
summary(sem_results[["Stem extract_Wild"]])

# Residual correlations between physical resistance traits and growth responses

res_correlations <- rbind(coefs(sem_results[["Control_Domesticated"]])[c(13:21),],
                          coefs(sem_results[["Leaf extract_Domesticated"]])[c(13:21),],
                          coefs(sem_results[["Stem extract_Domesticated"]])[c(13:21),],
                          coefs(sem_results[["Control_Wild"]])[c(13:21),],
                          coefs(sem_results[["Leaf extract_Wild"]])[c(13:21),],
                          coefs(sem_results[["Stem extract_Wild"]])[c(13:21),])
