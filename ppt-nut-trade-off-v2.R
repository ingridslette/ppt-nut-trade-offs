
library(tidyverse)
library(lme4)
library(lmerTest)
library(performance)

## LOAD AND FILTER BIOMASS DATA
biomass <- read.csv("/Users/ingridslette/Desktop/NutNet/full-biomass-2026-06-04.csv",
                 na.strings = c("NULL","NA"))

str(biomass)
summary(biomass)

biomass <- biomass |> 
  filter(
    live == 1,
    year_trt > 0, 
    trt %in% c("Control", "NPK"))

unique(biomass$trt)
unique(biomass$year_trt)
unique(biomass$year)
unique(biomass$live)

site_year_counts_biomass <- biomass %>%
  group_by(site_code, trt) %>% 
  summarise(year_count = n_distinct(year), .groups = 'drop')

sites_with_6_years_biomass <- site_year_counts_biomass %>%
  filter(year_count >= 6) %>%
  group_by(site_code) %>% 
  filter(n_distinct(trt) == 2) 

biomass6 <- biomass %>%
  filter(site_code %in% sites_with_6_years_biomass$site_code)

unique(biomass6$site_code)


## LOAD AND FILTER PRECIPITATION DATA
precip <- read.csv("/Users/ingridslette/Desktop/NutNet/ppt_pet_annual_gs_only_2025-10-09.csv")

str(precip)
summary(precip)

unique(precip$site_code)

precip <- filter(precip, year >= 1983)
precip <- filter(precip, year < 2025) 
unique(precip$year)

precip <- precip %>%
  group_by(site_code) %>%
  mutate(
    avg_ppt = mean(ppt, na.rm = TRUE),
    sd_ppt = sd(ppt, na.rm = TRUE),
    precip_z = ((ppt - avg_ppt) / sd_ppt)
  ) %>%
  ungroup()

unique(biomass6$site_code)
unique(precip$site_code)

# Join biomass and precip data
biomass_precip <- inner_join(biomass6, precip, by = c("site_code", "year"))
unique(biomass_precip$site_code)

# Apply 6-year filer again to see if that dropped any sites (due to lack of 2025 precip data)
biomass_precip <- biomass_precip %>%
  group_by(site_code, year) %>% filter(n_distinct(trt) == 2) %>%
  group_by(site_code) %>% filter(n_distinct(year) >= 6) %>%
  ungroup()
unique(biomass_precip$site_code)

# Filter to keep only sites with an observed ppt range that spans at least +- 1 sd of long-term avg
biomass_precip <- biomass_precip %>%
  group_by(site_code) %>%
  mutate(min_ppt = min(ppt, na.rm = TRUE),
         max_ppt = max(ppt, na.rm = TRUE)) %>%
  ungroup()

biomass_precip <- biomass_precip %>%
  group_by(site_code) %>%
  filter(min_ppt <= (avg_ppt - sd_ppt), max_ppt >= (avg_ppt + sd_ppt)) %>%
  ungroup()

unique(biomass_precip$site_code)


##### TOTAL BIOMASS ANALYSIS #####

## Create dataframe for total biomass analysis
total_biomass_precip <- biomass_precip %>%
  group_by(site_code, site_name, year, year_trt, trt, block, plot, ppt, avg_ppt, sd_ppt, precip_z) %>%
  summarise(total_biomass = sum(mass), .groups = "drop") %>%
  mutate(log_mass = log(total_biomass))

# site-year-average (sya) dataframe (as opposed to site-plot-year values)

# total_biomass_precip_sya <- total_biomass_precip %>%
#   group_by(site_code, site_name, year, year_trt, trt, ppt, avg_ppt, sd_ppt, precip_z) %>%
#   summarise(log_mass = mean(log_mass),
#             total_biomass = mean(total_biomass),
#             .groups = "drop")

total_biomass_precip_sya <- total_biomass_precip %>%
  group_by(site_code, site_name, year, year_trt, trt, ppt, avg_ppt, sd_ppt, precip_z) %>%
  summarise(total_biomass = mean(total_biomass), .groups = "drop") %>%
  mutate(log_mean_mass = log(total_biomass))


## SLOPE FUNCTION
get_slope_stats <- function(data, xvar, yvar) {
  data <- data[!is.na(data[[xvar]]) & !is.na(data[[yvar]]), ]
  mod <- tryCatch(lm(reformulate(xvar, yvar), data = data),
                  error = function(e) NULL)
  if (is.null(mod) || anyNA(coef(mod))) {
    return(data.frame(slope = NA_real_, se = NA_real_,
                      ci_low = NA_real_, ci_high = NA_real_,
                      p_value = NA_real_, n = nrow(data)))
  }
  s  <- summary(mod)$coefficients
  ci <- suppressWarnings(confint(mod))
  data.frame(
    slope   = s[xvar, "Estimate"],
    se      = s[xvar, "Std. Error"],
    ci_low  = ci[xvar, 1],
    ci_high = ci[xvar, 2],
    p_value = s[xvar, "Pr(>|t|)"],
    n       = nrow(data)
  )
}

## SLOPE CLASSIFICATION FUNCTION - classify slopes as a positive, negative, or none, based on CI 
classify_slope <- function(ci_low, ci_high) {
  case_when(
    is.na(ci_low) | is.na(ci_high) ~ NA_character_,
    ci_low  > 0 ~ "positive",
    ci_high < 0 ~ "negative",
    TRUE        ~ "none"
  )
}

## LOG RESPONSE RATIO VARIANCE FUNCTION (Hedges, Gurevitch & Curtis 1999)
lrr_variance <- function(mean_trt, sd_trt, n_trt, mean_ctrl, sd_ctrl, n_ctrl) {
  (sd_trt^2)  / (n_trt  * mean_trt^2) +
    (sd_ctrl^2) / (n_ctrl * mean_ctrl^2)
}


## TOTAL BIOMASS SLOPES PER SITE
slopes_total_biomass <- total_biomass_precip_sya %>%
  filter(trt == "Control") %>%
  group_by(site_code) %>%
  group_modify(~ get_slope_stats(.x, "precip_z", "log_mean_mass")) %>%
  ungroup() %>%
  rename(control_slope = slope, control_se = se,
         control_ci_low = ci_low, control_ci_high = ci_high,
         control_p = p_value, control_n = n)


## TOTAL BIOMASS LRR PER SITE
site_stats_total_biomass <- total_biomass_precip_sya %>%
  group_by(site_code, trt) %>%
  summarize(
    trt_mean = mean(total_biomass, na.rm = TRUE),
    trt_sd   = sd(total_biomass, na.rm = TRUE),
    trt_n    = sum(!is.na(total_biomass)),
    .groups = "drop"
  )

lrr_total_biomass <- site_stats_total_biomass %>%
  pivot_wider(names_from = trt, values_from = c(trt_mean, trt_sd, trt_n)) %>%
  mutate(
    lrr_mass        = log(trt_mean_NPK / trt_mean_Control),
    lrr_mass_var     = lrr_variance(trt_mean_NPK, trt_sd_NPK, trt_n_NPK,
                                    trt_mean_Control, trt_sd_Control, trt_n_Control),
    lrr_mass_se      = sqrt(lrr_mass_var),
    lrr_mass_ci_low  = lrr_mass - 1.96 * lrr_mass_se,
    lrr_mass_ci_high = lrr_mass + 1.96 * lrr_mass_se
  ) %>%
  select(site_code, lrr_mass, lrr_mass_se, lrr_mass_ci_low, lrr_mass_ci_high)


## TOTAL BIOMASS SLOPES AND LRRS
## join, classify and tally, graph
slopes_lrr_total_biomass <- left_join(slopes_total_biomass, lrr_total_biomass, by = "site_code") %>%
  mutate(
    precip_response = classify_slope(control_ci_low, control_ci_high),
    npk_response    = classify_slope(lrr_mass_ci_low, lrr_mass_ci_high),
    quadrant = paste(npk_response, precip_response, sep = " / ")
  )

response_counts_total_biomass <- slopes_lrr_total_biomass %>%
  count(npk_response, precip_response, name = "n_sites") %>%
  complete(npk_response    = c("positive", "none", "negative"),
           precip_response = c("positive", "none", "negative"),
           fill = list(n_sites = 0))

print(response_counts_total_biomass)

quadrant_counts_total_biomass <- slopes_lrr_total_biomass %>%
  count(quadrant, name = "n_sites") 

print(quadrant_counts_total_biomass)

plot_total_biomass <- ggplot(slopes_lrr_total_biomass,
                        aes(x = lrr_mass, y = control_slope, color = quadrant)) +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Response to NPK (log response ratio)",
       y = "Response to precipitation\n(slope of mass vs. precipitation)",
       title = "Total biomass responses to precipitation and NPK",
       color = "Response type") +
  theme_bw(base_size = 14)

plot_total_biomass

bar_total_biomass <- ggplot(response_counts_total_biomass,
                            aes(x = npk_response, y = n_sites, fill = precip_response)) +
  geom_col() +
  labs(x = "NPK response", y = "Number of sites",
       fill = "Precipitation response") +
  theme_bw(base_size = 14)

bar_total_biomass


bar_quadrant_total_biomass <- ggplot(quadrant_counts_total_biomass,
                            aes(x = quadrant, y = n_sites)) +
  geom_col() +
  labs(x = "", y = "Number of sites") +
  theme_bw(base_size = 14)

bar_quadrant_total_biomass


## TOTAL BIOMASS STATISTICAL MODELS
model_total_biomass <- lm(control_slope ~ lrr_mass, data = slopes_lrr_total_biomass)
summary(model_total_biomass)


mixed_model_total_biomass <- lmer(log_mass ~ trt * precip_z + (1 | site_code/plot) + (1 | site_code:year), 
                         data = total_biomass_precip)
summary(mixed_model_total_biomass)


##### BIOMASS BY FUNTIONAL GROUP ANALYSIS #####

## Create dataframe for biomass analysis by functional group (FG) 
fg_biomass <- biomass_precip %>%
  filter(category %in% c("GRAMINOID", "FORB", "WOODY", "LEGUME"))

# 0-fill categories found in some but not all plots and years at a site
years_by_site_fg_biomass <- fg_biomass %>%
  distinct(site_code, year, year_trt)

plots_fg_biomass <- fg_biomass %>%
  distinct(site_code, site_name, block, plot, trt)

fg_by_site <- fg_biomass %>%
  distinct(site_code, category)

fg_biomass_full_design <- plots_fg_biomass %>%
  left_join(fg_by_site, by = "site_code", relationship = "many-to-many") %>%
  left_join(years_by_site_fg_biomass, by = "site_code", relationship = "many-to-many")

fg_biomass_focal_cols <- fg_biomass %>%
  distinct(site_code, plot, category, year, mass)

fg_biomass_complete <- fg_biomass_full_design %>%
  left_join(fg_biomass_focal_cols,
            by = c("site_code", "plot", "category", "year")) %>%
  mutate(mass = replace_na(mass, 0),
         log_mass = log(mass + 0.01))

# Join precip data to 0-filled fg biomass data
fg_biomass_complete_precip <- fg_biomass_complete %>%
  inner_join(precip, by = c("site_code", "year"))

# Site-year-average (sya) dataframe (as opposed to site-plot-year values)
fg_biomass_complete_precip_sya <- fg_biomass_complete_precip %>%
  group_by(site_code, site_name, year, year_trt, trt, category, ppt, avg_ppt, sd_ppt, precip_z) %>%
  summarise(mass = mean(mass), .groups = "drop") %>%
  mutate(log_mean_mass = log(mass + 0.01))

# Keep only site-FGs where the FG was present in at least 4 years
fg_keep <- fg_biomass_complete_precip_sya %>%
  filter(mass > 0) %>%
  group_by(site_code, category) %>%
  summarise(n_years_present = n_distinct(year), .groups = "drop") %>%
  filter(n_years_present >= 4)

fg_biomass_complete_precip_sya <- fg_biomass_complete_precip_sya %>%
  semi_join(fg_keep, by = c("site_code", "category"))


## FG BIOMASS SLOPES PER SITE 
slopes_fg_biomass <- fg_biomass_complete_precip_sya %>%
  filter(trt == "Control") %>%
  group_by(site_code, category) %>%
  group_modify(~ get_slope_stats(.x, "precip_z", "log_mean_mass")) %>%
  ungroup() %>%
  rename(control_slope = slope, control_se = se,
         control_ci_low = ci_low, control_ci_high = ci_high,
         control_p = p_value, control_n = n)


## FG BIOMASS LRR PER SITE
site_stats_fg_biomass <- fg_biomass_complete_precip_sya %>%
  group_by(site_code, category, trt) %>%
  summarize(
    trt_mean = mean(mass, na.rm = TRUE) + 0.01,
    trt_sd   = sd(mass, na.rm = TRUE),
    trt_n    = sum(!is.na(mass)),
    .groups = "drop"
  ) 

lrr_fg_biomass <- site_stats_fg_biomass %>%
  pivot_wider(names_from = trt, values_from = c(trt_mean, trt_sd, trt_n)) %>%
  mutate(
    lrr_mass        = log(trt_mean_NPK / trt_mean_Control),
    lrr_mass_var    = lrr_variance(trt_mean_NPK, trt_sd_NPK, trt_n_NPK,
                                    trt_mean_Control, trt_sd_Control, trt_n_Control),
    lrr_mass_se     = sqrt(lrr_mass_var),
    lrr_mass_ci_low  = lrr_mass - 1.96 * lrr_mass_se,
    lrr_mass_ci_high = lrr_mass + 1.96 * lrr_mass_se
  ) %>%
  select(site_code, category, lrr_mass, lrr_mass_se, lrr_mass_ci_low, lrr_mass_ci_high)


## FG BIOMASS SLOPES AND LRRS
## join, classify and tally, graph
slopes_lrr_fg_biomass <- left_join(slopes_fg_biomass, lrr_fg_biomass,
                                       by = c("site_code", "category")) %>%
  mutate(
    precip_response = classify_slope(control_ci_low, control_ci_high),
    npk_response    = classify_slope(lrr_mass_ci_low, lrr_mass_ci_high),
    quadrant = paste(npk_response, precip_response, sep = " / ")
  )

response_counts_fg_biomass <- slopes_lrr_fg_biomass %>%
  filter(!is.na(precip_response), !is.na(npk_response)) %>%
  count(category, npk_response, precip_response, name = "n_sites") %>%
  complete(category, 
           npk_response    = c("positive", "none", "negative"),
           precip_response = c("positive", "none", "negative"),
           fill = list(n_sites = 0))

print(response_counts_fg_biomass)

quadrant_counts_fg_biomass <- slopes_lrr_fg_biomass %>%
  filter(!is.na(precip_response), !is.na(npk_response)) %>%
  count(category, quadrant, name = "n_sites")

print(quadrant_counts_fg_biomass)


plot_fg_biomass <- ggplot(slopes_lrr_fg_biomass,
                         aes(x = lrr_mass, y = control_slope, color = quadrant)) +
  geom_point(size = 2, alpha = 0.6) +
  facet_wrap(~ category) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Response to NPK (log response ratio)",
       y = "Response to precipitation\n(slope of cover vs. precipitation)",
       title = "FG Biomass responses to precipitation and NPK",
       color = "Response category") +
  theme_bw(base_size = 14)

plot_fg_biomass


bar_fg_biomass <- ggplot(response_counts_fg_biomass,
                             aes(x = npk_response, y = n_sites, fill = precip_response)) +
  geom_col() +
  facet_wrap(~ category) +
  labs(x = "NPK response", y = "Number of site-fg combinations",
       fill = "Precipitation response") +
  theme_bw(base_size = 14)

bar_fg_biomass


## FG BIOMASS STATISTICAL MODELS
model_fg_biomass <- lmer(control_slope ~ lrr_mass + category + (1 | site_code), data = slopes_lrr_fg_biomass)
summary(model_fg_biomass)

model_fg_biomass_x <- lmer(control_slope ~ lrr_mass * category + (1 | site_code), data = slopes_lrr_fg_biomass)
summary(model_fg_biomass_x)

mixed_model_fg_biomass <- lmer(log_mass ~ trt * precip_z  + category + 
                                 (1 | site_code/plot) + (1 | site_code:year), 
                                  data = fg_biomass_complete_precip)
summary(mixed_model_fg_biomass)















## LOAD AND FILTER COVER DATA
cover <- read.csv("/Users/ingridslette/Desktop/NutNet/full-cover-2026-06-04.csv",
                  na.strings = c("NULL","NA"))

str(cover)
summary(cover)

cover <- cover %>%
  filter(
    live == 1,
    year_trt > 0,
    trt %in% c("Control", "NPK")
  )

site_year_counts_cover <- cover %>%
  group_by(site_code, trt) %>% 
  summarise(year_count = n_distinct(year), .groups = 'drop')

sites_with_6_years_cover <- site_year_counts_cover %>%
  filter(year_count >= 6) %>%
  group_by(site_code) %>% 
  filter(n_distinct(trt) == 2) 

cover6 <- cover %>%
  filter(site_code %in% sites_with_6_years_cover$site_code)

unique(cover6$site_code) ## more sites with >5 years of cover data than with >5 years of biomass data
unique(biomass6$site_code) ## every site with >5 years of cover data has >5 years of biomass data

cover6 <- cover %>%
  filter(site_code %in% biomass_precip$site_code)

unique(cover6$site_code)
unique(cover6$year)


# Add rows for species found in some but not all plots and years at a site
years_by_site <- cover6 %>%
  distinct(site_code, year, year_trt)

plots <- cover6 %>%
  distinct(site_code, site_name, block, plot, subplot, trt)

taxa_by_site <- cover6 %>%
  distinct(site_code, Taxon)

full_design <- plots %>%
  left_join(taxa_by_site, by = "site_code", relationship = "many-to-many") %>%
  left_join(years_by_site, by = "site_code", relationship = "many-to-many")

cover_focal_cols <- cover6 %>%
  distinct(site_code, plot, Taxon, year, max_cover)

cover_complete <- full_design %>%
  left_join(cover_focal_cols,
            by = c("site_code", "plot", "Taxon", "year")) %>%
  mutate(max_cover = replace_na(max_cover, 0))

cover_trait_cols <- cover6 %>%
  distinct(site_code, Taxon, Family, functional_group, local_lifeform, 
           local_lifespan, local_provenance, ps_path)

cover_complete <- cover_complete %>%
  left_join(cover_trait_cols, by = c("site_code", "Taxon"), relationship = "many-to-many")

# Join precip data to 0-filled cover data
cover_precip <- inner_join(cover_complete, precip, by = c("site_code", "year"))


### -------------------------------------------------------------------------
### 5. COVER: precipitation-response slopes (control plots only, per site x Taxon)
### -------------------------------------------------------------------------

cover_slopes_control <- cover_precip %>%
  filter(trt == "Control") %>%
  group_by(site_code, Taxon) %>%
  group_modify(~ get_slope_stats(.x, "ppt", "max_cover")) %>%
  ungroup() %>%
  rename(Control_slope = slope, Control_se = se,
         Control_ci_low = ci_low, Control_ci_high = ci_high,
         Control_p = p_value, Control_n = n)

### 6. COVER: NPK log response ratio + its CI, per site x Taxon
### (keeps your existing +0.01 pseudocount so zero-cover taxa are still defined)

site_taxon_cover_stats <- cover_precip %>%
  group_by(site_code, Taxon, trt) %>%
  summarize(
    trt_mean = mean(max_cover, na.rm = TRUE),
    trt_sd   = sd(max_cover, na.rm = TRUE),
    trt_n    = sum(!is.na(max_cover)),
    .groups = "drop"
  ) %>%
  mutate(trt_mean = trt_mean + 0.01)

site_taxon_lrr_cover <- site_taxon_cover_stats %>%
  pivot_wider(names_from = trt, values_from = c(trt_mean, trt_sd, trt_n)) %>%
  mutate(
    lrr_cover        = log(trt_mean_NPK / trt_mean_Control),
    lrr_cover_var    = lrr_variance(trt_mean_NPK, trt_sd_NPK, trt_n_NPK,
                                    trt_mean_Control, trt_sd_Control, trt_n_Control),
    lrr_cover_se     = sqrt(lrr_cover_var),
    lrr_cover_ci_low  = lrr_cover - 1.96 * lrr_cover_se,
    lrr_cover_ci_high = lrr_cover + 1.96 * lrr_cover_se
  ) %>%
  select(site_code, Taxon, lrr_cover, lrr_cover_se, lrr_cover_ci_low, lrr_cover_ci_high)

### 7. COVER: join, classify, tally

spp_site_slopes_lrr_cover <- left_join(cover_slopes_control, site_taxon_lrr_cover,
                                       by = c("site_code", "Taxon")) %>%
  mutate(
    precip_response = classify_response(Control_ci_low, Control_ci_high),
    npk_response    = classify_response(lrr_cover_ci_low, lrr_cover_ci_high),
    quadrant = paste(npk_response, precip_response, sep = " / ")
  )

cover_quadrant_counts <- spp_site_slopes_lrr_cover %>%
  filter(!is.na(precip_response), !is.na(npk_response)) %>%
  count(npk_response, precip_response, name = "n_site_species") %>%
  complete(npk_response    = c("positive", "none", "negative"),
           precip_response = c("positive", "none", "negative"),
           fill = list(n_site_species = 0))

print(cover_quadrant_counts)

### 8. COVER: figure

site_cover_fig <- ggplot(spp_site_slopes_lrr_cover,
                         aes(x = lrr_cover, y = Control_slope, color = quadrant)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Response to NPK (log response ratio)",
       y = "Response to precipitation\n(slope of cover vs. precipitation)",
       title = "Site-species cover responses to precipitation and NPK",
       color = "Response category") +
  theme_bw(base_size = 14)

site_cover_fig

cover_quadrant_bar <- ggplot(cover_quadrant_counts,
                             aes(x = precip_response, y = n_site_species, fill = npk_response)) +
  geom_col(position = "dodge") +
  labs(x = "Precipitation response", y = "Number of site-species combinations",
       fill = "NPK response",
       title = "Number of site-species combinations by response category (cover)") +
  theme_bw(base_size = 14)

cover_quadrant_bar

