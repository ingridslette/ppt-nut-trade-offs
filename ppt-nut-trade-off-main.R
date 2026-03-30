library(tidyverse)
library(lme4)
library(lmerTest)
library(performance)
#library(boot)
#library(MuMIn)
#library(MASS)
#library(broom)
#library(broom.mixed)
#library(emmeans)
#library(ggpubr)
#library(cowplot)
#library(ggeffects)

### Loading, viewing, and filtering precipitation and mass data 

mass <- read.csv('/Users/ingridslette/Library/CloudStorage/GoogleDrive-slett152@umn.edu/Shared drives/NutNet_DRAGNet_Shared_External/NutNet Shared/NutNet Core Data/comb-by-plot-clim-soil-diversity-2025-12-09.csv',
                 na.strings = c("NULL","NA"))

unique(mass$site_code)
unique(mass$trt)
unique(mass$year_trt)

mass1 <- filter(mass, year_trt > 0)
unique(mass1$year_trt)

mass1 <- filter(mass1, trt %in% c("Control", "NPK")) 
unique(mass1$trt)

mass1 <- mass1 %>%
  mutate(live_mass = case_when(
    !is.na(vascular_live_mass) | !is.na(nonvascular_live_mass) ~ 
      rowSums(across(c(vascular_live_mass, nonvascular_live_mass, standing_dead_mass)), na.rm = TRUE),
    is.na(vascular_live_mass) & is.na(nonvascular_live_mass) ~ 
      unsorted_live_mass
  )
  )

site_year_counts <- mass1 %>%
  group_by(site_code, trt) %>%
  filter(!is.na(live_mass)) %>% 
  summarise(year_count = n_distinct(year), .groups = 'drop')

sites_with_6_years <- site_year_counts %>%
  filter(year_count >= 6) %>%
  group_by(site_code) %>% 
  filter(n_distinct(trt) == 2) 

mass6 <- mass1 %>%
  filter(site_code %in% sites_with_6_years$site_code)

sites_with_5_years <- site_year_counts %>%
  filter(year_count >= 5) %>%
  group_by(site_code) %>% 
  filter(n_distinct(trt) == 2) 

mass5 <- mass1 %>%
  filter(site_code %in% sites_with_5_years$site_code)

unique(mass5$site_code) 
unique(mass6$site_code)

## popped over to script "calculate-gs-ppt-pet.R" here, to get growing season ppt for the sites included in mass2
## exported that as csv and now loading it here

ppt_data <- read.csv("/Users/ingridslette/Desktop/NutNet/ppt_pet_annual_gs_only_2025-10-09.csv")

unique(ppt_data$site_code)

ppt_data <- filter(ppt_data, year >= 1983)
ppt_data <- filter(ppt_data, year < 2025) 
unique(ppt_data$year)

ppt_data <- ppt_data %>%
  group_by(site_code) %>%
  mutate(
    avg_ppt = mean(ppt, na.rm = TRUE),
    avg_pet = mean(pet, na.rm = TRUE),
    sd_ppt = sd(ppt, na.rm = TRUE)
  ) %>%
  ungroup()

unique(mass6$site_code)
unique(ppt_data$site_code)

mass_ppt <- inner_join(mass6, ppt_data, by = c("site_code", "year"))

unique(mass_ppt$site_code)

mass_ppt <- mass_ppt %>%
  mutate(log_mass = log10(live_mass),
         log_ppt = log10(ppt))

mass_ppt <- mass_ppt %>%
  group_by(site_code) %>%
  mutate(min_ppt = min(ppt, na.rm = TRUE),
         max_ppt = max(ppt, na.rm = TRUE)) %>%
  ungroup()


# Filter to keep only sites with an observed ppt range that spans at least +- 1 sd of long-term avg
mass_ppt <- mass_ppt %>%
  group_by(site_code) %>%
  filter(min_ppt <= (avg_ppt - sd_ppt), max_ppt >= (avg_ppt + sd_ppt)) %>%
  ungroup()

unique(mass_ppt$site_code)


### Calculate slope of ppt vs mass at each site

site_codes <- unique(mass_ppt$site_code)

site_slopes <- data.frame(site_code = character(),
                      control_slope = numeric(), 
                      npk_slope = numeric(),
                      slope_diff = numeric(),
                      stringsAsFactors = FALSE)

for (site in site_codes) {
  site_data_control <- subset(mass_ppt, site_code == site & trt == "Control")
  site_data_npk <- subset(mass_ppt, site_code == site & trt == "NPK")
  control_model <- lm(log_mass ~ log_ppt, data = site_data_control)
  npk_model <- lm(log_mass ~ log_ppt, data = site_data_npk)
  control_slope <- coef(control_model)["log_ppt"]
  npk_slope <- coef(npk_model)["log_ppt"]
  site_slopes <- rbind(site_slopes, data.frame(
    site_code = site,
    control_slope = control_slope,
    npk_slope = npk_slope,
    slope_diff = npk_slope - control_slope
  ))
}





### Covariate analyses

## Incorporating log response ratio of mass to trt
lrr_df <- mass_ppt %>%
  group_by(site_code, year) %>%
  summarize(
    lrr_mass = log(mean(live_mass[trt == "NPK"], na.rm = TRUE) /
                     mean(live_mass[trt == "Control"], na.rm = TRUE))
  )

## Incorporating C3/C4 and annual/perennial information from cover data
cover <- read.csv("/Users/ingridslette/Desktop/NutNet/full-cover_2025-01-31.csv",
                  na.strings = c("NULL","NA"))

cover <- cover %>%
  filter(site_code %in% site_codes)

dat1 <- subset(cover, is.na(ps_path) == TRUE)%>%
  separate(Taxon, into = c("Genus", "Species"), remove = FALSE, sep = " ")

dat1$ps_path <- ifelse(dat1$Genus == "BINERTIA" | dat1$Genus == "TIDESTROMIA" | dat1$Genus == "PECTIS" | dat1$Genus == "EUPLOCA" | dat1$Genus == "BULBOSTYLIS" | dat1$Genus == "CYPERUS" | dat1$Genus == "FIMBRISTYLIS" | dat1$Genus == "CHAMAESYCE" | dat1$Genus == "ALLIONIA" | dat1$Genus == "CALLIGONUM" | dat1$Genus == "PORTULACA" | dat1$Genus == "EUPHORBIA", "C4",
                       ifelse(dat1$Genus == "BELAPHARIS" | dat1$Genus == "AERVA" | dat1$Genus == "ALTERNANTHERA" | dat1$Genus == "ATRIPLEX" | dat1$Genus == "SUAEDA" | dat1$Genus == "TECTICORNIA" | dat1$Genus == "FLAVERIA" | dat1$Genus == "POLYCARPOREA" | dat1$Genus == "ELEOCHARS" | dat1$Genus == "RHYNCHOSPORA" | dat1$Genus == "EUPHORBIA" | dat1$Genus == "MOLLUGO" | dat1$Genus == "BOERHAVIA" | dat1$Genus == "BASSIA" | dat1$Family == "Poaceae", NA,
                              "C3"))
unique(dat1$ps_path)

dat1 <- dat1 %>%
  rename(ps_path2 = ps_path)

names(dat1)
names(cover)

cover <- cover %>% 
  left_join(dat1, by = c("year", "site_name", "site_code", "block", "plot", "subplot", "year_trt", "trt", 
                         "Family", "Taxon", "live", "local_provenance", "local_lifeform", "local_lifespan", 
                         "functional_group", "max_cover"))

cover <- cover %>%
  mutate(ps_path2 = if_else(is.na(ps_path2) & !is.na(ps_path), ps_path, ps_path2))

cover_by_site_plot_year <- cover %>%
  group_by(site_code, plot, year) %>%
  summarise(
    total_cover = sum(max_cover, na.rm = TRUE),
    c4_cover = if (any(ps_path2 == "C4", na.rm = TRUE)) {
      sum(max_cover[ps_path2 == "C4"], na.rm = TRUE)} else {0},
    c4_proportion = c4_cover / total_cover,
    annual_cover = if (any(local_lifespan == "ANNUAL", na.rm = TRUE)) {
      sum(max_cover[local_lifespan == "ANNUAL"], na.rm = TRUE)} else {0},
    annual_proportion = annual_cover / total_cover,
    .groups = "drop"
  )


# Calculating light interception from PAR

mass_ppt <- mass_ppt %>%
  mutate(light_intercepted = 1 - (Ground_PAR / Ambient_PAR))

# Getting rid of unneeded columns

mass_ppt_edited <- mass_ppt %>%
  dplyr::select(site_code, block, plot, continent, country, region, habitat, trt, 
                year, live_mass, log_mass, ppt, log_ppt, prev_ppt, year_trt,
                proportion_par, avg_ppt, sd_ppt, p05_ppt, p95_ppt, p10_ppt, 
                p90_ppt,rich, MAT_v2, AI, PET, MAP_v2, ppt_pet, light_intercepted)

unique(mass_ppt_edited$site_code)

# Joining LRR and cover data to main dataframe

mass_ppt_edited <- mass_ppt_edited %>% 
  left_join(lrr_df, by = c("site_code", "year"))

mass_ppt_edited <- mass_ppt_edited %>% 
  left_join(cover_by_site_plot_year, by = c("site_code", "plot", "year"))

mass_ppt_edited <- na.omit(mass_ppt_edited)
unique(mass_ppt_edited$site_code)


## Covariate model of mass

full_model <- lmer(log_mass ~ trt * (log_ppt + light_intercepted + AI + rich + prev_ppt
                                     + c4_proportion + annual_proportion + year_trt)
                   + (1 | year) + (1 | site_code/block), 
                   data = mass_ppt_edited, REML = FALSE, na.action = "na.fail")
summary(full_model)


full_model_table <- dredge(full_model, m.lim = c(NA, 6))
full_model_avg <- model.avg(get.models(full_model_table, subset = delta < 50))
summary(full_model_avg); sw(full_model_avg)

r2_full_model <- performance::r2(full_model)
print(r2_full_model)

nrow(full_model_table)


mass_lrr_mass_plot <- ggplot(data = mass_ppt_edited, 
                             aes(x = lrr_mass, y = live_mass, color = trt, fill = trt, shape = trt)) +
  geom_point(alpha = 0.2) + 
  geom_smooth(method = "lm", se = F, alpha = 0.3) +
  labs(x = "Biomass Response \nRatio", y = "Biomass (g/m²)", 
       color = "Treatment", shape = "Treatment", fill = "Treatment") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                     #, labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_fill_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                    #, labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_shape_manual(values = c("Control" = 21, "NPK" = 24)
                     #, labels = c("Control" = "Control", "NPK" = "Fertilized")
  )

mass_par_plot <- ggplot(data = mass_ppt_edited, 
                        aes(x = light_intercepted, y = live_mass, color = trt, fill = trt, shape = trt)) +
  geom_point(alpha = 0.2) + 
  geom_smooth(method = "lm", se = F, alpha = 0.3) +
  labs(x = "Light Interception", y = "Biomass (g/m²)",
       color = "Treatment", shape = "Treatment", fill = "Treatment") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                     #, labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_fill_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                    #, labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_shape_manual(values = c("Control" = 21, "NPK" = 24)
                     #, labels = c("Control" = "Control", "NPK" = "Fertilized")
  )

mass_ai_plot <- ggplot(data = mass_ppt_edited, 
                       aes(x = AI, y = live_mass, color = trt, fill = trt, shape = trt)) +
  geom_point(alpha = 0.3) +
  labs(x = "Aridity Index", y = "Biomass (g/m²)",
       color = "Treatment", shape = "Treatment", fill = "Treatment") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                     #, labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_fill_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                    #, labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_shape_manual(values = c("Control" = 16, "NPK" = 17)
                     #, labels = c("Control" = "Control", "NPK" = "Fertilized")
  )

mass_rich_plot <- ggplot(data = mass_ppt_edited, 
                         aes(x = rich, y = live_mass, color = trt, shape = trt, fill = trt)) +
  geom_point(alpha = 0.3) +
  labs(x = "Species Richness", y = "Biomass (g/m²)",
       color = "Treatment", shape = "Treatment", fill = "Treatment") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                     #, labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_fill_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                    #, labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_shape_manual(values = c("Control" = 21, "NPK" = 24)
                     #, labels = c("Control" = "Control", "NPK" = "Fertilized")
  )

mass_prev_ppt_plot <- ggplot(data = mass_ppt_edited, 
                             aes(x = prev_ppt, y = live_mass, color = trt, shape = trt, fill = trt)) +
  geom_point(alpha = 0.3) +
  labs(x = "Previous Growing Season \nPrecipitation (mm)", y = "Biomass (g/m²)",
       color = "Treatment", shape = "Treatment", fill = "Treatment") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                     #, labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_fill_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                    #, labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_shape_manual(values = c("Control" = 21, "NPK" = 24)
                     #, labels = c("Control" = "Control", "NPK" = "Fertilized")
  )

mass_c4_plot <- ggplot(data = mass_ppt_edited, 
                       aes(x = c4_proportion, y = live_mass, color = trt, shape = trt, fill = trt)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = F, alpha = 0.3) +
  labs(x = expression("Proportion C"[4]), y = "Biomass (g/m²)",
       color = "Treatment", shape = "Treatment", fill = "Treatment") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                     #, labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_fill_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                    #, labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_shape_manual(values = c("Control" = 21, "NPK" = 24)
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  )

mass_annual_plot <- ggplot(data = mass_ppt_edited, 
                           aes(x = annual_proportion, y = live_mass, color = trt, fill = trt, shape = trt)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = F, alpha = 0.3) +
  labs(x = "Proportion Annual \nLifespan", y = "Biomass (g/m²)",
       color = "Treatment", shape = "Treatment", fill = "Treatment") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                     #, labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_fill_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                    #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_shape_manual(values = c("Control" = 21, "NPK" = 24)
                     #, labels = c("Control" = "Control", "NPK" = "Fertilized")
  )

mass_year_trt_plot <- ggplot(data = mass_ppt_edited, 
                             aes(x = year_trt, y = live_mass, color = trt, shape = trt, fill = trt)) +
  geom_point(alpha = 0.3) +
  labs(x = "Treatment Year", y = "Biomass (g/m²)",
       color = "Treatment", shape = "Treatment", fill = "Treatment") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_fill_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                    #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_shape_manual(values = c("Control" = 21, "NPK" = 24)
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  )


mass_covar_fig_all <- ggarrange(mass_ai_plot,
                                
                                mass_par_plot + rremove("ylab") +
                                  theme(axis.text.y = element_blank()),
                                
                                mass_prev_ppt_plot + rremove("ylab") +
                                  theme(axis.text.y = element_blank()),
                                
                                mass_year_trt_plot + rremove("ylab") +
                                  theme(axis.text.y = element_blank()),
                                
                                mass_rich_plot,
                                
                                mass_annual_plot + rremove("ylab") +
                                  theme(axis.text.y = element_blank()),
                                
                                mass_c4_plot + rremove("ylab") +
                                  theme(axis.text.y = element_blank()),
                                
                                align = 'hv',
                                ncol = 4, nrow = 2,
                                common.legend = TRUE, 
                                legend = "bottom")
mass_covar_fig_all


## Covariate models of R2 and slope

results_long <- data.frame(site_code = character(), 
                           trt = character(), 
                           r2 = numeric(), 
                           slope = numeric(),
                           mean = numeric(),
                           stringsAsFactors = FALSE)

for (site in site_codes) {
  site_data_control <- subset(mass_ppt, site_code == site & trt == "Control")
  site_data_npk <- subset(mass_ppt, site_code == site & trt == "NPK")
  control_model <- lm(log_mass ~ log_ppt, data = site_data_control)
  npk_model <- lm(log_mass ~ log_ppt, data = site_data_npk)
  control_r2 <- summary(control_model)$r.squared
  npk_r2 <- summary(npk_model)$r.squared
  control_slope <- coef(control_model)["log_ppt"]
  npk_slope <- coef(npk_model)["log_ppt"]
  control_mean <- mean(site_data_control$log_mass, na.rm = TRUE)
  npk_mean <- mean(site_data_npk$log_mass, na.rm = TRUE)
  results_long <- rbind(results_long, data.frame(
    site_code = site,
    trt = "Control",
    r2 = control_r2,
    slope = control_slope,
    mean = control_mean
  ))
  results_long <- rbind(results_long, data.frame(
    site_code = site,
    trt = "NPK",
    r2 = npk_r2,
    slope = npk_slope,
    mean = npk_mean
  ))
}

averages <- mass_ppt_edited %>%
  group_by(site_code, trt) %>%
  summarise(
    avg_light = mean(light_intercepted, na.rm = TRUE),
    avg_avg_ppt = mean(avg_ppt, na.rm = TRUE),
    avg_mat = mean(MAT_v2, na.rm = TRUE),
    avg_map = mean(MAP_v2, na.rm = TRUE),
    avg_richness = mean(rich, na.rm = TRUE),
    avg_lrr_mass = mean(lrr_mass, na.rm = TRUE),
    avg_c4_proportion = mean(c4_proportion, na.rm = TRUE),
    avg_annual_proportion = mean(annual_proportion, na.rm = TRUE),
    avg_ai = mean(AI, na.rm = TRUE)
  )

results_with_averages <- results_long %>%
  left_join(averages, by = c("site_code", "trt"))

full_r2_model <- lm(r2 ~ trt * (avg_light + avg_ai + avg_richness + avg_lrr_mass
                                + avg_c4_proportion + avg_annual_proportion), 
                    data = results_with_averages, na.action = "na.fail")

summary(full_r2_model)
full_r2_model_table <- dredge(full_r2_model, m.lim = c(NA, 6))
full_r2_model_avg <- model.avg(get.models(full_r2_model_table, subset = delta < 10))
summary(full_r2_model_avg); sw(full_r2_model_avg)

nrow(full_r2_model_table)

ai_r2_model <- lm(r2 ~ trt * avg_ai, data = results_with_averages)
summary(ai_r2_model)

annual_r2_model <- lm(r2 ~ trt * avg_annual_proportion, data = results_with_averages)
summary(annual_r2_model)


full_slope_model <- lm(slope ~ trt * (avg_light  + avg_ai + avg_richness + avg_lrr_mass
                                      + avg_c4_proportion + avg_annual_proportion), 
                       data = results_with_averages, na.action = "na.fail")

summary(full_slope_model)
full_slope_model_table <- dredge(full_slope_model, m.lim = c(NA, 6))
full_slope_model_avg <- model.avg(get.models(full_slope_model_table, subset = delta < 10))
summary(full_slope_model_avg); sw(full_slope_model_avg)

nrow(full_slope_model_table)

ai_slope_model <- lm(slope ~ trt * avg_ai, data = results_with_averages)
summary(ai_slope_model)

model_quad <- lm(slope ~ poly(avg_ai, 2) * trt, data = results_with_averages)

model_log <- lm(slope ~ log(avg_ai) * trt, data = results_with_averages)

AIC(model_quad, model_log, ai_slope_model)

summary(model_quad)$adj.r.squared
summary(model_log)$adj.r.squared
summary(model_log)$adj.r.squared

annual_slope_model <- lm(slope ~ trt * avg_annual_proportion, data = results_with_averages)
summary(annual_slope_model)


avg_slope_r2_by_trt <- results_long %>%
  group_by(trt) %>%
  summarise(
    avg_slope = mean(slope, na.rm = TRUE),
    se_slope = sd(slope, na.rm = TRUE) / sqrt(sum(!is.na(slope))),
    avg_r2 = mean(r2, na.rm = TRUE),
    se_r2 = sd(r2, na.rm = TRUE) / sqrt(sum(!is.na(r2)))
  )

avg_slope_r2_by_trt

control_slope <- avg_slope_r2_by_trt$avg_slope[avg_slope_r2_by_trt$trt == "Control"]
npk_slope <- avg_slope_r2_by_trt$avg_slope[avg_slope_r2_by_trt$trt == "NPK"]
control_r2 <- avg_slope_r2_by_trt$avg_r2[avg_slope_r2_by_trt$trt == "Control"]
npk_r2 <- avg_slope_r2_by_trt$avg_r2[avg_slope_r2_by_trt$trt == "NPK"]

percent_increase_slope <- ((npk_slope - control_slope) / control_slope) * 100
percent_increase_slope

percent_increase_r2 <- ((npk_r2 - control_r2) / control_r2) * 100
percent_increase_r2


## Main graph(s)

fig_2 <- ggplot(predictions_allsites, 
                aes(x = 10^log_ppt, y = predicted_mass, colour = trt)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = mass_lower, ymax = mass_upper, fill = trt), alpha = 0.3, , colour = NA) +
  scale_color_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                     #, labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_fill_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                    #, labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  labs(x = "Growing Season Precipitation (mm)", y = "Biomass (g/m²)",
       color = "Treatment", fill = "Treatment",) +
  scale_y_continuous(limits = c(0, 1000)) +
  theme_bw(base_size = 14) +
  theme(legend.title = element_text(size = 12))

fig_2  


predictions <- predictions %>%
  left_join(results_with_averages, by = c("site_code", "trt"))

fig2_control <- ggplot(subset(predictions, trt == "Control"), aes(x = 10^log_ppt, y = predicted_mass)) +
  geom_line(aes(group = site_code, colour = r2)) +
  geom_ribbon(data = subset(predictions_allsites, trt == "Control"),
              aes(x = 10^log_ppt, ymin = mass_lower, ymax = mass_upper), 
              fill = "#0092E0", alpha = 0.3) +
  geom_line(data = subset(predictions_allsites, trt == "Control"), 
            aes(x = 10^log_ppt, y = predicted_mass),
            color = "#0092E0") +
  labs(x = "Growing Season Precipitation (mm)", y = "Biomass (g/m²)", 
       title = "Control", colour = "Correlation \nStrength") +
  scale_y_continuous(limits = c(0, 2300)) +
  scale_colour_gradient2(low = "#E4D3EE", mid = "#B185DB", high = "#423073",
                         midpoint = 0.3) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 10),
    axis.title.y = element_text(size = 14, margin = margin(r = 20))
  )

fig2_control

fig2_npk <- ggplot(subset(predictions, trt == "NPK"), aes(x = 10^log_ppt, y = predicted_mass)) +
  geom_line(aes(group = site_code, colour = r2)) +
  geom_ribbon(data = subset(predictions_allsites, trt == "NPK"),
              aes(x = 10^log_ppt, ymin = mass_lower, ymax = mass_upper), 
              fill = "#ff924c", alpha = 0.3) +
  geom_line(data = subset(predictions_allsites, trt == "NPK"), 
            aes(x = 10^log_ppt, y = predicted_mass),
            color = "#ff924c") +
  labs(x = "Growing Season Precipitation (mm)", y = "Biomass (g/m²)", 
       title = "NPK", colour = "Correlation \nStrength") +
  scale_y_continuous(limits = c(0, 2300)) +
  scale_colour_gradient2(low = "#E4D3EE", mid = "#B185DB", high = "#423073",
                         midpoint = 0.3) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 14)
  )

fig2_npk

fig2_both <- ggarrange(
  fig2_control + rremove("xlab"),
  fig2_npk + 
    rremove("ylab") + rremove("xlab") +
    theme(
      axis.text.y = element_blank(), 
      axis.ticks.y = element_blank()
    ),
  ncol = 2,
  common.legend = TRUE,
  legend = "right",
  align = 'hv'
)

fig2_both

fig2_both <- annotate_figure(
  fig2_both,
  bottom = text_grob("Growing Season Precipitation (mm)", size = 14)
)

fig2_both


## Calculating and testing trt effect on RUE 

mass_ppt_edited <- mass_ppt_edited %>% 
  mutate(rue = live_mass/ppt)

rue_trt_model <- lmer(rue ~ trt + (1 | site_code) + (1 | year), data = mass_ppt_edited)
summary(rue_trt_model)

AIC(rue_trt_model)

r2_rue_trt_model <- performance::r2(rue_trt_model)
print(r2_rue_trt_model)

ggplot(mass_ppt_edited, aes(x = trt, y = rue)) +
  geom_boxplot() +
  labs(x = "Treatment", y = "Rain Use Efficiency") +
  theme_bw(14)

avg_rue_by_trt <- mass_ppt_edited %>%
  group_by(trt) %>%
  summarise(
    avg_rue = mean(rue, na.rm = TRUE),
    se_rue = sd(rue, na.rm = TRUE) / sqrt(sum(!is.na(rue)))
  )

avg_rue_by_trt

control_rue <- avg_rue_by_trt$avg_rue[avg_rue_by_trt$trt == "Control"]
npk_rue <- avg_rue_by_trt$avg_rue[avg_rue_by_trt$trt == "NPK"]

percent_increase_rue <- ((npk_rue - control_rue) / control_rue) * 100
percent_increase_rue



## Analyzing only data from non-extreme precipitation years

mass_ppt_nominal <- subset(mass_ppt_edited, ppt > p05_ppt & ppt < p95_ppt)

main_model_non_extreme <- lmer(
  log_mass ~ log_ppt * trt + (1 | site_code / block) + (1 | year),
  data = mass_ppt_nominal,
  na.action = na.exclude
)

summary(main_model_non_extreme)

AIC(main_model_non_extreme)

r2_main_non_extreme <- performance::r2(main_model_non_extreme)
print(r2_main_non_extreme)

non_extreme_plot <- ggplot(data = mass_ppt_nominal, 
                           aes(x = ppt, y = live_mass, color = trt, fill = trt, shape = trt)) +
  geom_point(alpha = 0.7) + 
  geom_smooth(method = "lm", alpha = 0.3) +
  labs(x = "Growing Season Precipitation (mm)", y = "Biomass (g/m²)", 
       color = "Treatment", shape = "Treatment", fill = "Treatment",
       title = "Non-Extreme Precipitation Years") +
  theme_bw(base_size = 16) +
  scale_color_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_fill_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                    #, labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_shape_manual(values = c("Control" = 21, "NPK" = 24)
                     #, labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14))

non_extreme_plot


### Comparing control vs. NPK R2 - non-extreme years

site_codes_nominal <- unique(mass_ppt_nominal$site_code)

results_nominal <- data.frame(site_code = character(), 
                              control_r2 = numeric(), 
                              npk_r2 = numeric(),
                              r2_diff = numeric(),
                              control_slope = numeric(), 
                              npk_slope = numeric(),
                              slope_diff = numeric(),
                              stringsAsFactors = FALSE)

for (site in site_codes_nominal) {
  site_data_control <- subset(mass_ppt_nominal, site_code == site & trt == "Control")
  site_data_npk <- subset(mass_ppt_nominal, site_code == site & trt == "NPK")
  control_model <- lm(log_mass ~ log_ppt, data = site_data_control)
  npk_model <- lm(log_mass ~ log_ppt, data = site_data_npk)
  control_r2 <- summary(control_model)$r.squared
  npk_r2 <- summary(npk_model)$r.squared
  control_slope <- coef(control_model)["log_ppt"]
  npk_slope <- coef(npk_model)["log_ppt"]
  results_nominal <- rbind(results_nominal, data.frame(
    site_code = site,
    control_r2 = control_r2,
    npk_r2 = npk_r2,
    r2_diff = npk_r2 - control_r2,
    control_slope = control_slope,
    npk_slope = npk_slope,
    slope_diff = npk_slope - control_slope
  ))
}

paired_t_test_r2_nominal <- t.test(results_nominal$control_r2, results_nominal$npk_r2, paired = TRUE)
paired_t_test_r2_nominal


## Analyzing only data from extreme dry years

mass_ppt_dry <- subset(mass_ppt_edited, ppt<p05_ppt)

main_model_dry <- lmer(log_mass ~ log_ppt * trt + (1 | site_code / block) + (1 | year), 
                       data = mass_ppt_dry)

summary(main_model_dry)

AIC(main_model_dry)

r2_main_dry <- performance::r2(main_model_dry)
r2_main_dry

dry_plot <- ggplot(data = mass_ppt_dry, 
                   aes(x = ppt, y = live_mass)) +
  geom_point(aes(color = trt, fill = trt, shape = trt), alpha = 0.7) + 
  geom_smooth(method = "lm", color = "#6F6F6F", alpha = 0.3) +
  labs(x = "Growing Season Precipitation (mm)", y = "Biomass (g/m²)", 
       color = "Treatment", shape = "Treatment", fill = "Treatment",
       title = "Extreme Dry Years") +
  theme_bw(base_size = 16) +
  scale_color_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_fill_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                    #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_shape_manual(values = c("Control" = 21, "NPK" = 24)
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14))

dry_plot

### Comparing control vs. NPK R2 - extreme dry years

site_codes_dry <- unique(mass_ppt_dry$site_code)

results_dry <- data.frame(site_code = character(), 
                          control_r2 = numeric(), 
                          npk_r2 = numeric(),
                          r2_diff = numeric(),
                          control_slope = numeric(), 
                          npk_slope = numeric(),
                          slope_diff = numeric(),
                          stringsAsFactors = FALSE)

for (site in site_codes_dry) {
  site_data_control <- subset(mass_ppt_dry, site_code == site & trt == "Control")
  site_data_npk <- subset(mass_ppt_dry, site_code == site & trt == "NPK")
  control_model <- lm(log_mass ~ log_ppt, data = site_data_control)
  npk_model <- lm(log_mass ~ log_ppt, data = site_data_npk)
  control_r2 <- summary(control_model)$r.squared
  npk_r2 <- summary(npk_model)$r.squared
  control_slope <- coef(control_model)["log_ppt"]
  npk_slope <- coef(npk_model)["log_ppt"]
  results_dry <- rbind(results_dry, data.frame(
    site_code = site,
    control_r2 = control_r2,
    npk_r2 = npk_r2,
    r2_diff = npk_r2 - control_r2,
    control_slope = control_slope,
    npk_slope = npk_slope,
    slope_diff = npk_slope - control_slope
  ))
}

paired_t_test_r2_dry <- t.test(results_dry$control_r2, results_dry$npk_r2, paired = TRUE)
paired_t_test_r2_dry


mean_se_control_dry <- mean_se(results_dry$control_r2)
mean_se_npk_dry <- mean_se(results_dry$npk_r2)

mean_se_control_dry
mean_se_npk_dry


## Analyzing only data from extreme wet years

main_model_p95 <- lmer(log_mass ~ log_ppt * trt + (1 | site_code / block) + (1 | year), 
                       data = subset(mass_ppt_edited, ppt>p95_ppt))

summary(main_model_p95)

p95_plot <- ggplot(data = subset(mass_ppt_edited, ppt>p95_ppt), 
                   aes(x = ppt, y = live_mass, color = trt, fill = trt, shape = trt)) +
  geom_point(alpha = 0.7) + 
  geom_smooth(method = "lm", alpha = 0.3) +
  labs(x = "Growing Season \nPrecipitation (mm)", y = "Biomass (g/m²)", 
       color = "Treatment", shape = "Treatment", fill = "Treatment",
       title = "Extreme Wet Years") +
  theme_bw(base_size = 16) +
  scale_color_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_fill_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                    #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_shape_manual(values = c("Control" = 21, "NPK" = 24)
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14))

p95_plot


### Calculating and graphing effect sizes

results_graphing <- data.frame(site_code = character(), 
                               trt = character(), 
                               r2 = numeric(), 
                               slope = numeric(),
                               mean = numeric(),
                               stringsAsFactors = FALSE)

for (site in site_codes) {
  site_data_control <- subset(mass_ppt, site_code == site & trt == "Control")
  site_data_npk <- subset(mass_ppt, site_code == site & trt == "NPK")
  control_model <- lm(live_mass ~ ppt, data = site_data_control)
  npk_model <- lm(live_mass ~ ppt, data = site_data_npk)
  control_r2 <- summary(control_model)$r.squared
  npk_r2 <- summary(npk_model)$r.squared
  control_slope <- coef(control_model)["ppt"]
  npk_slope <- coef(npk_model)["ppt"]
  control_mean <- mean(site_data_control$live_mass, na.rm = TRUE)
  npk_mean <- mean(site_data_npk$live_mass, na.rm = TRUE)
  results_graphing <- rbind(results_graphing, data.frame(
    site_code = site,
    trt = "Control",
    r2 = control_r2,
    slope = control_slope,
    mean = control_mean
  ))
  results_graphing <- rbind(results_graphing, data.frame(
    site_code = site,
    trt = "NPK",
    r2 = npk_r2,
    slope = npk_slope,
    mean = npk_mean
  ))
}

mean_model <- lmer(mean ~ trt + (1| site_code), data = results_graphing)
summary(mean_model)

slope_model <- lmer(slope ~ trt + (1| site_code), data = results_graphing)
summary(slope_model)

r2_model <- lmer(r2 ~ trt + (1| site_code), data = results_graphing)
summary(r2_model)

mean_estimate <- 188.55
mean_se <- 20.59
mean_resid_sd <- 86.15
n_mean <- 70

slope_estimate <- 0.3174
slope_se <- 0.1210
slope_resid_sd <- 0.5063
n_slope <- 70

r2_estimate <- -0.001285
r2_se <- 0.018987
r2_resid_sd <- 0.07943
n_r2 <- 70

rue_estimate <- 0.5365
rue_se <- 0.04865
rue_resid_sd <- 1.0923
n_rue <- 2046


calc_cohen_d <- function(estimate, resid_sd, n) {
  d <- estimate / resid_sd
  n1 <- n / 2
  n2 <- n / 2
  SE_d <- sqrt((n1 + n2) / (n1 * n2) + (d^2) / (2 * (n1 + n2))) 
  lower_CI <- d - 1.96 * SE_d
  upper_CI <- d + 1.96 * SE_d
  return(c(d, lower_CI, upper_CI))
}

mean_results <- calc_cohen_d(mean_estimate, mean_resid_sd, n_mean)
slope_results <- calc_cohen_d(slope_estimate, slope_resid_sd, n_slope)
r2_results <- calc_cohen_d(r2_estimate, r2_resid_sd, n_r2)
rue_results <- calc_cohen_d(rue_estimate, rue_resid_sd, n_rue)

cohen_d_df <- data.frame(
  Variable = c("Biomass", "Rain Use \nEfficiency", "Sensitivity", "Correlation \n Strength"),
  Cohen_d = c(mean_results[1], rue_results[1], slope_results[1], r2_results[1]),
  Lower_CI = c(mean_results[2], rue_results[2], slope_results[2], r2_results[2]),
  Upper_CI = c(mean_results[3], rue_results[3], slope_results[3], r2_results[3])
)

cohen_d_df$Variable <- factor(cohen_d_df$Variable, levels = c("Correlation \n Strength", "Sensitivity", "Rain Use \nEfficiency", "Biomass"))

es_fig <- ggplot(cohen_d_df, aes(x = Cohen_d, y = Variable, color = Variable)) +
  geom_point(size = 4) +
  geom_errorbar(aes(xmin = Lower_CI, xmax = Upper_CI), width = 0.2) +
  labs(x = "Effect Size",
       y = "") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("#6F6F6F", "#ff924c", "#ff924c", "#ff924c")) +
  theme_bw(base_size = 14) +
  theme(axis.text.y = element_text(size = 12, color = "black"), legend.position = "none")

es_fig

es_fig_y <- es_fig + 
  labs(y = "variable") + 
  theme(axis.title.y = element_text(color = "white"))
es_fig_y

fig_2_inside <- fig_2 + theme(legend.position = "inside",
                              legend.position.inside = c(0.25, 0.8))

es_fig2 <- ggarrange(fig_2_inside, 
                     es_fig_y,
                     labels = c("a", "b"),
                     label.x = 0.3,
                     label.y = 0.95,
                     font.label = list(size = 12),
                     align = "hv",
                     nrow = 1, ncol = 2
)
es_fig2 


## Graphing covariates against R2 and slope

results_with_averages_graphing <- results_graphing %>%
  left_join(averages, by = c("site_code", "trt"))

r2_lrr_mass_plot <- ggplot(data = results_with_averages, 
                           aes(x = avg_lrr_mass, y = r2, color = trt, shape = trt, fill = trt)) +
  geom_point(alpha = 0.7) +
  labs(x = "Biomass Response \nRatio", y = "Correlation Strength", 
       color = "Treatment", shape = "Treatment", fill = "Treatment") +
  theme_bw(base_size = 14) +
  scale_color_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_fill_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                    #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_shape_manual(values = c("Control" = 21, "NPK" = 24)
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  )

r2_par_plot <- ggplot(data = results_with_averages, 
                      aes(x = avg_light, y = r2, color = trt, shape = trt, fill = trt)) +
  geom_point(alpha = 0.7) +
  labs(x = "Light Interception", y = "Correlation Strength", 
       color = "Treatment", shape = "Treatment", fill = "Treatment") +
  theme_bw(base_size = 14) +
  #  scale_y_continuous(limits = c(0, 0.6)) +
  scale_color_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_fill_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                    #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_shape_manual(values = c("Control" = 21, "NPK" = 24)
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  )


r2_ai_plot <- ggplot(data = results_with_averages, 
                     aes(x = avg_ai, y = r2, color = trt, shape = trt, fill = trt)) +
  geom_point(alpha = 0.7) +
  labs(x = "Aridity Index", y = "Correlation Strength",
       color = "Treatment", shape = "Treatment", fill = "Treatment") +
  theme_bw(base_size = 14) +
  #  scale_y_continuous(limits = c(0, 0.6)) +
  scale_color_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_fill_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                    #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_shape_manual(values = c("Control" = 21, "NPK" = 24)
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  )

r2_ai_plot2 <- ggplot(data = results_with_averages, 
                      aes(x = avg_ai, y = r2)) +
  geom_point(aes(color = trt, shape = trt, fill = trt), alpha = 0.7) + 
  geom_smooth(method = lm, formula = y ~ log(x), se = F, color = "#6F6F6F", linewidth = 0.75) +
  labs(x = "Aridity Index \n\nMore Arid   →   Less Arid", y = "Correlation Strength",
       color = "Treatment", shape = "Treatment", fill = "Treatment") +
  theme_bw(base_size = 14) +
  #  scale_y_continuous(limits = c(0, 0.6)) +
  scale_color_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                     #, labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_fill_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                    #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_shape_manual(values = c("Control" = 21, "NPK" = 24)
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  )

r2_ai_plot2

r2_ai_plot3 <- ggplot(data = results_with_averages, 
                      aes(x = avg_ai, y = r2)) +
  geom_point(aes(color = trt, shape = trt, fill = trt), alpha = 0.7) + 
  geom_smooth(method = lm, formula = y ~ log(x), se = F, color = "#6F6F6F", linewidth = 0.75) +
  labs(x = "Aridity Index", y = "Correlation Strength",
       color = "Treatment", shape = "Treatment", fill = "Treatment") +
  theme_bw(base_size = 14) +
  #  scale_y_continuous(limits = c(0, 0.6)) +
  scale_color_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_fill_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                    #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_shape_manual(values = c("Control" = 21, "NPK" = 24)
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  )

r2_ai_plot3

r2_ai_plot2_trt <- ggplot(data = results_with_averages, 
                          aes(x = avg_ai, y = r2, color = trt, shape = trt, fill = trt)) +
  geom_point(alpha = 0.7) + 
  geom_smooth(method = lm, formula = y ~ log(x), se = F, linewidth = 0.75) +
  labs(x = "Aridity Index \n\nMore Arid   →   Less Arid", y = "Correlation Strength",
       color = "Treatment", shape = "Treatment", fill = "Treatment") +
  theme_bw(base_size = 14) +
  #  scale_y_continuous(limits = c(0, 0.6)) +
  scale_color_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_fill_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                    #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_shape_manual(values = c("Control" = 21, "NPK" = 24)
                     #, labels = c("Control" = "Control", "NPK" = "Fertilized")
  )

r2_ai_plot2_trt

r2_rich_plot <- ggplot(data = results_with_averages, 
                       aes(x = avg_richness, y = r2, color = trt, shape = trt, fill = trt)) +
  geom_point(alpha = 0.7) + 
  labs(x = "Species Richness", y = "Correlation Strength", 
       color = "Treatment", shape = "Treatment", fill = "Treatment") +
  theme_bw(base_size = 14) +
  # scale_y_continuous(limits = c(0, 0.6)) +
  scale_color_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_fill_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                    #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_shape_manual(values = c("Control" = 21, "NPK" = 24)
                     #, labels = c("Control" = "Control", "NPK" = "Fertilized")
  )


r2_c4_plot <- ggplot(data = results_with_averages, 
                     aes(x = avg_c4_proportion, y = r2, color = trt, shape = trt, fill = trt)) +
  geom_point(alpha = 0.7) + 
  labs(x = expression("Proportion C"[4]), y = "Correlation Strength", 
       color = "Treatment", shape = "Treatment", fill = "Treatment") +
  theme_bw(base_size = 14) +
  #  scale_y_continuous(limits = c(0, 0.6)) +
  scale_color_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_fill_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                    #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_shape_manual(values = c("Control" = 21, "NPK" = 24)
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  )


r2_annual_plot <- ggplot(data = results_with_averages, 
                         aes(x = avg_annual_proportion, y = r2, color = trt, shape = trt, fill = trt)) +
  geom_point(alpha = 0.7) + 
  labs(x ="Proportion Annual \nLifespan", y = "Correlation Strength", 
       color = "Treatment", shape = "Treatment", fill = "Treatment") +
  theme_bw(base_size = 14) +
  #  scale_y_continuous(limits = c(0, 0.6)) +
  scale_color_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_fill_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                    #, labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_shape_manual(values = c("Control" = 21, "NPK" = 24)
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  )


r2_covar_figure <- ggarrange(r2_ai_plot2, 
                             
                             r2_lrr_mass_plot + rremove("ylab") +
                               theme(axis.text.y = element_blank()),
                             
                             r2_par_plot + rremove("ylab") +
                               theme(axis.text.y = element_blank()),
                             
                             r2_rich_plot, 
                             
                             r2_annual_plot  + rremove("ylab") +
                               theme(axis.text.y = element_blank()),
                             
                             r2_c4_plot + rremove("ylab") +
                               theme(axis.text.y = element_blank()),
                             
                             ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom", align = 'hv')
r2_covar_figure


slope_lrr_mass_plot <- ggplot(data = results_with_averages, 
                              aes(x = avg_lrr_mass, y = slope, color = trt, shape = trt, fill = trt)) +
  geom_point(alpha = 0.7) + 
  labs(x = "Biomass Response \nRatio", y = "Sensitivity",
       color = "Treatment", shape = "Treatment", fill = "Treatment") +
  theme_bw(base_size = 14) +
  #  scale_y_continuous(limits = c(-1.1, 3.8)) +
  scale_color_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_fill_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                    #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_shape_manual(values = c("Control" = 21, "NPK" = 24)
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  )

slope_par_plot <- ggplot(data = results_with_averages, 
                         aes(x = avg_light, y = slope, color = trt, shape = trt, fill = trt)) +
  geom_point(alpha = 0.7) + 
  labs(x = "Light Interception", y = "Sensitivity",
       color = "Treatment", shape = "Treatment", fill = "Treatment") +
  theme_bw(base_size = 14) +
  #  scale_y_continuous(limits = c(-1.1, 3.8)) +
  scale_color_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_fill_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                    #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_shape_manual(values = c("Control" = 21, "NPK" = 24)
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  )

slope_ai_plot <- ggplot(data = results_with_averages, 
                        aes(x = avg_ai, y = slope, color = trt, shape = trt, fill = trt)) +
  geom_point(alpha = 0.7) + 
  geom_smooth(method = lm, se = FALSE, alpha = 0.2) +
  labs(x = "Aridity Index", y = "Sensitivity", 
       color = "Treatment", shape = "Treatment", fill = "Treatment") +
  theme_bw(base_size = 14) +
  #  scale_y_continuous(limits = c(-1.1, 3.8)) +
  scale_color_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_fill_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                    #, labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_shape_manual(values = c("Control" = 21, "NPK" = 24)
                     #, labels = c("Control" = "Control", "NPK" = "Fertilized")
  )

slope_ai_plot2 <- ggplot(data = results_with_averages, 
                         aes(x = avg_ai, y = slope, color = trt, shape = trt, fill = trt)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = lm, formula = y ~ log(x), se = FALSE, linewidth = 0.75) +
  labs(x = "Aridity Index", y = "Sensitivity",
       color = "Treatment", shape = "Treatment", fill = "Treatment") +
  theme_bw(base_size = 14) +
  #  scale_y_continuous(limits = c(-1.1, 3.8)) +
  scale_color_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                     #, labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_fill_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                    #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_shape_manual(values = c("Control" = 21, "NPK" = 24)
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  )
slope_ai_plot2


slope_rich_plot <- ggplot(data = results_with_averages, 
                          aes(x = avg_richness, y = slope, color = trt, shape = trt, fill = trt)) +
  geom_point(alpha = 0.7) + 
  labs(x = "Species Richness", y = "Sensitivity",
       color = "Treatment", shape = "Treatment", fill = "Treatment") +
  theme_bw(base_size = 14) +
  #  scale_y_continuous(limits = c(-1.1, 3.8)) +
  scale_color_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_fill_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                    #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_shape_manual(values = c("Control" = 21, "NPK" = 24)
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  )

slope_c4_plot <- ggplot(data = results_with_averages, 
                        aes(x = avg_c4_proportion, y = slope, color = trt, shape = trt, fill = trt)) +
  geom_point(alpha = 0.7) +
  labs(x = expression("Proportion C"[4]), y = "Sensitivity",
       color = "Treatment", shape = "Treatment", fill = "Treatment") +
  theme_bw(base_size = 14) +
  #  scale_y_continuous(limits = c(-1.1, 3.8)) +
  scale_color_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_fill_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                    #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_shape_manual(values = c("Control" = 21, "NPK" = 24)
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  )

slope_annual_plot <- ggplot(data = results_with_averages, 
                            aes(x = avg_annual_proportion, y = slope, color = trt, shape = trt, fill = trt)) +
  geom_point(alpha = 0.7) + 
  geom_smooth(method = lm, se = FALSE, linewidth = 0.75) +
  labs(x = "Proportion Annual \nLifespan", y = "Sensitivity",
       color = "Treatment", shape = "Treatment", fill = "Treatment") +
  theme_bw(base_size = 14) +
  #  scale_y_continuous(limits = c(-1.1, 3.8)) +
  scale_color_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_fill_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                    #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_shape_manual(values = c("Control" = 21, "NPK" = 24)
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  )


slope_covar_figure <- ggarrange(slope_ai_plot2, 
                                
                                slope_lrr_mass_plot + rremove("ylab") +
                                  theme(axis.text.y = element_blank()),
                                
                                slope_par_plot + rremove("ylab") +
                                  theme(axis.text.y = element_blank()),
                                
                                slope_rich_plot, 
                                
                                slope_annual_plot  + rremove("ylab") +
                                  theme(axis.text.y = element_blank()),
                                
                                slope_c4_plot + rremove("ylab") +
                                  theme(axis.text.y = element_blank()),
                                
                                ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom", align = 'hv')

slope_covar_figure


slope_r2_sig_fig <- ggarrange(slope_ai_plot2 + rremove("xlab"),
                              slope_annual_plot  + rremove("ylab") +
                                theme(axis.text.y = element_blank()),
                              r2_ai_plot2,
                              ncol = 2, nrow = 2, common.legend = TRUE, 
                              legend = "right", align = 'hv',
                              labels = c("a*", "b", "c*"),
                              label.x = 0.27,          # left aligned
                              label.y = 0.95,          # top aligned
                              font.label = list(size = 12)
)

slope_r2_sig_fig

slope_r2_ai_fig <- ggarrange(slope_ai_plot2 + rremove("xlab"),
                             r2_ai_plot2,
                             ncol = 1, nrow = 2, common.legend = TRUE, 
                             legend = "bottom", align = 'hv',
                             labels = c("a", "b"),
                             label.x = 0.26,          # left aligned
                             label.y = 0.95,          # top aligned
                             font.label = list(size = 12))
slope_r2_ai_fig


## Mapping sites

library(ggmap)
library(viridis)

register_stadiamaps("98dda47f-a35b-4ead-aa46-dda02c95d912", write = FALSE)
stadiamaps_key()
has_stadiamaps_key()

bbox <- c(left = -170, bottom = -60, right = 160, top = 75)
myMap <- get_stadiamap(bbox, zoom = 4, maptype = "stamen_terrain_background")

# Main map with overlaid legend
map <- ggmap(myMap) + 
  geom_point(data = mass_ppt, aes(x = longitude, y = latitude, fill = AI),
             shape = 21, size = 1.2, alpha = 0.2) +
  scale_fill_viridis_c(option = "plasma", direction = 1) +
  labs(x = "", y = "", fill = "Aridity Index") +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.65, 0.07),
    legend.direction = "horizontal",
    legend.title = element_text(color = "white", size = 14, face = "bold"),
    legend.text = element_text(color = "white", size = 10),
    legend.background = element_rect(fill = "transparent", colour = NA)
  )

map

aridity_fig <- ggplot(data = mass_ppt, aes(x = MAP_v2, y = MAT_v2, fill = AI)) +
  geom_point(shape = 21, color = "black", size = 1.2, stroke = 0.3) + 
  labs(x = "Mean Annual \nPrecipitation (mm)",
       y = "Mean Annual \nTemperature (°C)",
       fill = "Aridity \nIndex") +
  scale_fill_viridis_c(option = "plasma", direction = 1) +
  theme_bw() +
  theme(legend.position = "none")

map_aridity_inset <- ggdraw() + 
  draw_plot(map) +
  draw_plot(aridity_fig, x = 0.04, y = 0.07, width = 0.24, height = 0.35)

map_aridity_inset


# SUPPLEMENTAL ANALYSES

### Another approach to comparing control vs. NPK R2 - fit separate models for control and NPK data, calculate and compare z scores

model_control <- lmer(log_mass ~ log_ppt + (1 | site_code / year), 
                      data = subset(mass_ppt, trt == "Control"))
model_npk <- lmer(log_mass ~ log_ppt + (1 | site_code / year), 
                  data = subset(mass_ppt, trt == "NPK"))

summary(model_control)
summary(model_npk)

AIC(model_control, model_npk)

r2_control <- performance::r2(model_control)
r2_npk <- performance::r2(model_npk)

conditional_r2_control <- r2_control$R2_conditional
conditional_r2_npk <- r2_npk$R2_conditional

marginal_r2_control <- r2_control$R2_marginal
marginal_r2_npk <- r2_npk$R2_marginal

# Compare R2 values using Fisher's Z transformation
z_control <- 0.5 * log((1 + sqrt(marginal_r2_control)) / (1 - sqrt(marginal_r2_control)))
z_npk <- 0.5 * log((1 + sqrt(marginal_r2_npk)) / (1 - sqrt(marginal_r2_npk)))

n_control <- length(unique(subset(mass_ppt, trt == "Control")$site_code))
n_npk <- length(unique(subset(mass_ppt, trt == "NPK")$site_code))
se_diff <- sqrt((1 / (n_control - 3)) + (1 / (n_npk - 3)))

# Calculate the Z-score for the difference
z_diff <- (z_control - z_npk) / se_diff

p_value <- 2 * (1 - pnorm(abs(z_diff)))

cat("Z-score for the difference:", z_diff, "\n")
cat("P-value for the difference in marginal R-squared values:", p_value, "\n")


### Yet another approach to comparing control vs. NPK R2 - Bootstrapping to test for difference in R2 between Control and NKP models

r2_diff <- function(data, indices) {
  data_resampled <- data[indices, ]
  model_control <- lmer(log_mass ~ log_ppt + (1 | site_code), 
                        data = data_resampled[data_resampled$trt == "Control", ])
  model_npk <- lmer(log_mass ~ log_ppt + (1 | site_code), 
                    data = data_resampled[data_resampled$trt == "NPK", ])
  r2_control <- r.squaredGLMM(model_control)[2]
  r2_npk <- r.squaredGLMM(model_npk)[2]          
  return(r2_control - r2_npk)
}

set.seed(123)

boot_r2 <- boot(data = mass_ppt, statistic = r2_diff, R = 1000)

print(boot_r2)

boot_ci <- boot.ci(boot_r2, type = "perc")
print(boot_ci)


## site info table, etc.

site_info_table <- mass_ppt %>%
  dplyr::select(site_code, site_name, latitude, longitude, MAP_v2, MAT_v2, AI) %>%
  distinct()

write.csv(site_info_table, file = "/Users/ingridslette/Desktop/NutNet/site-info-table.csv")


site_block_counts <- mass_ppt %>%
  group_by(site_code) %>%
  summarise(block_count = n_distinct(block), .groups = 'drop')

site_plot_counts <- mass_ppt %>%
  group_by(site_code, block, trt) %>%
  summarise(plot_count = n_distinct(plot), .groups = 'drop')

arid_sites <- mass_ppt %>%
  group_by(site_code) %>%
  filter(AI <= 0.65)

unique(arid_sites$site_code)

water_limited_sites <- mass_ppt %>%
  group_by(site_code) %>%
  filter(avg_ppt<avg_pet)

unique(water_limited_sites$site_code)



