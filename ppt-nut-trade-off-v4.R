
library(tidyverse)
library(lme4)
library(lmerTest)
library(purrr)
library(broom)
library(ggpubr)

##### DEFINITIONS AND FUNCTIONS #####

pc      <- 0.01
fg_list <- c("GRAMINOID", "FORB", "WOODY", "LEGUME")

# Random-effect structure used for every model
rand_form <- "(1 | year) + (1 | block)"

# Fit one model per site (or site x FG or site x species) and return coefficient of interest
#   fixed = fixed-effect predictor ("trt" or "precip_z")
#   term  = coefficient name to extract ("trtNPK" or "precip_z")
# `singular` is TRUE when the fit is on the boundary (e.g., a variance
# component estimated at zero); the estimate is still returned.
fit_site_lmm <- function(data, fixed, term, response = "log_mass") {
  data <- data[!is.na(data[[response]]) & !is.na(data[[fixed]]), ]
  data$year  <- factor(data$year)
  data$block <- factor(data$block)
  
  out_na <- tibble(
    estimate = NA_real_, se = NA_real_, df = NA_real_,
    ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_,
    n_obs = nrow(data), n_years = n_distinct(data$year),
    singular = NA
  )
  
  f <- as.formula(paste(response, "~", fixed, "+", rand_form))
  m <- tryCatch(suppressMessages(lmer(f, data = data)),
                error = function(e) NULL)
  if (is.null(m)) return(out_na)
  
  cf <- summary(m)$coefficients
  if (!term %in% rownames(cf)) return(out_na)
  
  est  <- cf[term, "Estimate"]
  se   <- cf[term, "Std. Error"]
  dfs  <- cf[term, "df"]
  crit <- qt(0.975, dfs)
  
  tibble(
    estimate = est, se = se, df = dfs,
    ci_low = est - crit * se, ci_high = est + crit * se,
    p_value = cf[term, "Pr(>|t|)"],
    n_obs = nrow(data), n_years = n_distinct(data$year),
    singular = isSingular(m)
  )
}

# Classify an estimate as positive, negative, or none, based on its CI
classify_slope <- function(ci_low, ci_high) {
  case_when(
    is.na(ci_low) | is.na(ci_high) ~ NA_character_,
    ci_low  > 0 ~ "positive",
    ci_high < 0 ~ "negative",
    TRUE        ~ "none"
  )
}


#### LOAD AND FILTER BIOMASS DATA ####
biomass <- read.csv("/Users/ingridslette/Desktop/NutNet/full-biomass-2026-06-04.csv",
                    na.strings = c("NULL","NA"))

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


#### LOAD AND FILTER PRECIPITATION DATA ####
precip <- read.csv("/Users/ingridslette/Desktop/NutNet/ppt_pet_annual_gs_only_2025-10-09.csv")

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


#### LOAD AND FILTER COVER DATA ####
cover <- read.csv("/Users/ingridslette/Desktop/NutNet/full-cover-2026-06-04.csv",
                  na.strings = c("NULL","NA"))

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
  filter(site_code %in% biomass6$site_code)

unique(cover6$site_code)


#### SET UP TOTAL BIOMASS DATASET FOR ANALYSIS ####

## Join biomass and precip data
biomass_precip <- inner_join(biomass6, precip, by = c("site_code", "year"))
unique(biomass_precip$site_code)

## Apply 6-year filer again to see if that dropped any sites (due to lack of 2025 precip data)
biomass_precip <- biomass_precip %>%
  group_by(site_code, year) %>% filter(n_distinct(trt) == 2) %>%
  group_by(site_code) %>% filter(n_distinct(year) >= 6) %>%
  ungroup()

unique(biomass_precip$site_code)

## Filter to keep only sites with an observed ppt range that spans at least +- 1 sd of long-term avg
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

biomass_precip <- biomass_precip %>%
  mutate(trt = factor(trt, levels = c("Control", "NPK")))

## Calculate plot-level total biomass
total_biomass_precip <- biomass_precip %>%
  group_by(site_code, site_name, year, year_trt, trt, block, plot,
           ppt, avg_ppt, sd_ppt, precip_z) %>%
  summarise(total_biomass = sum(mass), .groups = "drop")

## Check how many plots per block
site_plot_counts <- biomass_precip %>%
  group_by(site_code, block, trt) %>%
  summarise(plot_count = n_distinct(plot), .groups = 'drop')

print(site_plot_counts, n=248)

## Calculate average total biomass per block
total_biomass_precip_block <- total_biomass_precip %>%
  group_by(site_code, site_name, year, year_trt, trt, block,
           ppt, avg_ppt, sd_ppt, precip_z) %>%
  summarise(total_biomass = mean(total_biomass), .groups = "drop") %>%
  mutate(log_mass = log(total_biomass))

## Check how many blocks per site
site_block_counts <- total_biomass_precip_block %>%
  group_by(site_code) %>%
  summarise(block_count = n_distinct(block), .groups = 'drop')

print(site_block_counts, n=37)


#### SET UP FG BIOMASS DATASET FOR ANALYSIS ####

## Filter to FGs of interest
fg_biomass <- biomass_precip %>%
  filter(category %in% fg_list)

## Check no duplicate FGs per plot
mass_cat_block_counts <- fg_biomass %>% 
  count(site_code, plot, year, category) %>% filter(n > 1)
nrow(mass_cat_block_counts)

## 0-fill the dataset 
## Add rows for FGs that were observed at a site but absent from some plot-years
plot_years <- biomass_precip %>%
  distinct(site_code, site_name, block, plot, trt, year, year_trt)

fg_by_site <- fg_biomass %>%
  distinct(site_code, category)

fg_biomass_complete <- plot_years %>%
  left_join(fg_by_site, by = "site_code", relationship = "many-to-many") %>%
  left_join(select(fg_biomass, site_code, plot, year, category, mass),
            by = c("site_code", "plot", "year", "category")) %>%
  mutate(mass = replace_na(mass, 0))

## Attach precip variables (already computed for biomass_precip)
precip_vars <- biomass_precip %>%
  distinct(site_code, year, ppt, avg_ppt, sd_ppt, precip_z)

fg_biomass_complete_precip <- fg_biomass_complete %>%
  inner_join(precip_vars, by = c("site_code", "year"))

## Keep only site-FGs present (mass > 0 in any plot) in at least 4 years
fg_keep <- fg_biomass_complete_precip %>%
  group_by(site_code, category, year) %>%
  summarise(present = any(mass > 0), .groups = "drop") %>%
  filter(present) %>%
  group_by(site_code, category) %>%
  summarise(n_years_present = n(), .groups = "drop") %>%
  filter(n_years_present >= 4)

fg_biomass_complete_precip <- fg_biomass_complete_precip %>%
  semi_join(fg_keep, by = c("site_code", "category")) 

## Check how many plots per block
site_plot_counts_fg_mass <- fg_biomass_complete_precip %>%
  group_by(site_code, block, trt) %>%
  summarise(plot_count = n_distinct(plot), .groups = 'drop')

print(site_plot_counts_fg_mass, n=234)

## Calculate average biomass per block per category
fg_biomass_precip_block <- fg_biomass_complete_precip %>%
  group_by(site_code, site_name, year, year_trt, trt, block,
           category, ppt, avg_ppt, sd_ppt, precip_z) %>%
  summarise(mass = mean(mass), .groups = "drop") %>%
  mutate(log_mass = log(mass + pc))

## Check how many blocks per site
site_block_counts_fg_mass <- fg_biomass_precip_block %>%
  group_by(site_code) %>%
  summarise(block_count = n_distinct(block), .groups = 'drop')

print(site_block_counts_fg_mass, n=35)


#### SET UP SPECIES COVER DATASET FOR ANALYSIS ####

## Standardize functional group naming
## Biomass doesn't separate grass vs graminoid, cover does
cover6 <- cover6 %>%
  mutate(functional_group = recode(functional_group, "GRASS" = "GRAMINOID"))

cover6 <- cover6 %>%
  filter(functional_group %in% fg_list)

## Check if same site-year combos as for biomass
biomass_site_years <- biomass_precip %>% distinct(site_code, year)

cover6 <- cover6 %>%
  semi_join(biomass_site_years, by = c("site_code", "year"))

# re-check 6-year/both-trt requirement on this reduced set
site_year_counts_cover <- cover6 %>%
  group_by(site_code, trt) %>%
  summarise(year_count = n_distinct(year), .groups = "drop")

sites_ok_cover <- site_year_counts_cover %>%
  filter(year_count >= 6) %>%
  group_by(site_code) %>%
  filter(n_distinct(trt) == 2)

unique(sites_ok_cover$site_code)

cover6 <- cover6 %>% filter(site_code %in% sites_ok_cover$site_code)

setdiff(unique(biomass_precip$site_code), unique(cover6$site_code))

site_year_counts_cover %>%
  anti_join(sites_ok_cover, by = "site_code") %>%
  filter(site_code %in% biomass_precip$site_code)


## 0-fill the dataset 
## Add rows for FGs that were observed at a site but absent from some plot-years
unique(cover6$functional_group)

plot_years <- cover6 %>%
  distinct(site_code, site_name, block, plot, trt, year, year_trt)

taxa_by_site <- cover6 %>%
  distinct(site_code, Taxon)

full_design <- plot_years %>%
  left_join(taxa_by_site, by = "site_code", relationship = "many-to-many")

cover_focal_cols <- cover6 %>%
  distinct(site_code, plot, Taxon, year, max_cover)

cover_complete <- full_design %>%
  left_join(cover_focal_cols,
            by = c("site_code", "plot", "Taxon", "year")) %>%
  mutate(max_cover = replace_na(max_cover, 0))

cover_trait_cols <- cover6 %>%
  distinct(site_code, Taxon, Family, functional_group, local_lifeform, 
           local_lifespan, local_provenance, ps_path)

## Check for inconsistent trait rows per site-taxon
dup_traits <- cover_trait_cols %>% count(site_code, Taxon) %>% filter(n > 1)
nrow(dup_traits)
print(dup_traits)

# Collapse to one row per site-taxon: first non-NA value in each trait column
first_non_na <- function(x) {
  x_clean <- x[!is.na(x)]
  if (length(x_clean) == 0) NA_character_ else x_clean[1]
}

cover_trait_cols <- cover_trait_cols %>%
  group_by(site_code, Taxon) %>%
  summarise(across(c(Family, functional_group, local_lifeform,
                     local_lifespan, local_provenance, ps_path),
                   first_non_na),
            .groups = "drop")

dup_traits <- cover_trait_cols %>% count(site_code, Taxon) %>% filter(n > 1)
nrow(dup_traits)

cover_complete <- cover_complete %>%
  left_join(cover_trait_cols, by = c("site_code", "Taxon"))

## Join precip data to 0-filled cover data
cover_precip <- inner_join(cover_complete, precip, by = c("site_code", "year"))

## Keep only site x species combos present (cover > 0 in any plot) in >= 4 years
species_keep <- cover_precip %>%
  group_by(site_code, Taxon, year) %>%
  summarise(present = any(max_cover > 0), .groups = "drop") %>%
  filter(present) %>%
  group_by(site_code, Taxon) %>%
  summarise(n_years_present = n(), .groups = "drop") %>%
  filter(n_years_present >= 4)

cover_precip_keep <- cover_precip %>%
  semi_join(species_keep, by = c("site_code", "Taxon")) 

## Check how many plots per block
site_plot_counts_cover <- cover_precip_keep %>%
  group_by(site_code, block, trt) %>%
  summarise(plot_count = n_distinct(plot), .groups = 'drop')

print(site_plot_counts_cover, n=242)

## Calculate average cover per block per species
cover_precip_block <- cover_precip_keep %>%
  group_by(site_code, site_name, year, year_trt, trt, block, Taxon,
           Family, functional_group, local_lifeform, local_lifespan,
           local_provenance, ppt, avg_ppt, sd_ppt, precip_z) %>%
  summarise(cover = mean(max_cover), .groups = "drop") %>%
  mutate(log_cover = log(cover + pc),
         trt = factor(trt, levels = c("Control", "NPK")))

## Check how many blocks per site
site_block_counts_cover <- cover_precip_block %>%
  group_by(site_code) %>%
  summarise(block_count = n_distinct(block), .groups = 'drop')

print(site_block_counts_cover, n=36)


##### TOTAL BIOMASS ANALYSIS #####

## Calculate precip effect (precip-biomass slope) per site for control plots and for NPK plots
slopes_total_biomass_control <- total_biomass_precip_block %>%
  filter(trt == "Control") %>%
  group_by(site_code) %>%
  group_modify(~ fit_site_lmm(.x, fixed = "precip_z", term = "precip_z")) %>%
  ungroup() %>%
  rename_with(~ paste0("ppt_", .x), -site_code) %>%
  rename(ppt_slope = ppt_estimate, ppt_p = ppt_p_value)

slopes_total_biomass_npk <- total_biomass_precip_block %>%
  filter(trt == "NPK") %>%
  group_by(site_code) %>%
  group_modify(~ fit_site_lmm(.x, fixed = "precip_z", term = "precip_z")) %>%
  ungroup() %>%
  rename_with(~ paste0("ppt_", .x), -site_code) %>%
  rename(ppt_slope = ppt_estimate, ppt_p = ppt_p_value)

## Calculate NPK effect (log response ratio) per site: trt coefficient
lrr_total_biomass <- total_biomass_precip_block %>%
  group_by(site_code) %>%
  group_modify(~ fit_site_lmm(.x, fixed = "trt", term = "trtNPK")) %>%
  ungroup() %>%
  rename_with(~ paste0("npk_", .x), -site_code) %>%
  rename(npk_lrr = npk_estimate, npk_p = npk_p_value)

## Check for singularity 
table(slopes_total_biomass_control$ppt_singular, useNA = "ifany")
table(slopes_total_biomass_npk$ppt_singular, useNA = "ifany")
table(lrr_total_biomass$npk_singular, useNA = "ifany")

## Join, classify and tally
slopes_lrr_total_biomass_ctrl <- left_join(slopes_total_biomass_control, lrr_total_biomass,
                                      by = "site_code") %>%
  mutate(
    precip_response = classify_slope(ppt_ci_low, ppt_ci_high),
    npk_response    = classify_slope(npk_ci_low, npk_ci_high),
    quadrant = paste(npk_response, precip_response, sep = " / ")
  )

slopes_lrr_total_biomass_npk <- left_join(slopes_total_biomass_npk, lrr_total_biomass,
                                           by = "site_code") %>%
  mutate(
    precip_response = classify_slope(ppt_ci_low, ppt_ci_high),
    npk_response    = classify_slope(npk_ci_low, npk_ci_high),
    quadrant = paste(npk_response, precip_response, sep = " / ")
  )

response_counts_total_biomass_ctrl <- slopes_lrr_total_biomass_ctrl %>%
  count(npk_response, precip_response, name = "n_sites") %>%
  complete(npk_response    = c("positive", "none", "negative"),
           precip_response = c("positive", "none", "negative"),
           fill = list(n_sites = 0))%>%
  mutate(quadrant = paste(npk_response, precip_response, sep = " / "))

print(response_counts_total_biomass_ctrl)

quadrant_counts_total_biomass_ctrl <- slopes_lrr_total_biomass_ctrl %>%
  count(quadrant, name = "n_sites")

print(quadrant_counts_total_biomass_ctrl)

response_counts_total_biomass_npk <- slopes_lrr_total_biomass_npk %>%
  count(npk_response, precip_response, name = "n_sites") %>%
  complete(npk_response    = c("positive", "none", "negative"),
           precip_response = c("positive", "none", "negative"),
           fill = list(n_sites = 0))%>%
  mutate(quadrant = paste(npk_response, precip_response, sep = " / "))

print(response_counts_total_biomass_npk)

quadrant_counts_total_biomass_npk <- slopes_lrr_total_biomass_npk %>%
  count(quadrant, name = "n_sites")

print(quadrant_counts_total_biomass_npk)


## Analyze
model_total_biomass_ctrl <- lm(ppt_slope ~ npk_lrr, data = slopes_lrr_total_biomass_ctrl)
summary(model_total_biomass_ctrl)

model_total_biomass_npk <- lm(ppt_slope ~ npk_lrr, data = slopes_lrr_total_biomass_npk)
summary(model_total_biomass_npk)

mixed_model_total_biomass <- lmer(log_mass ~ trt * precip_z +
                                    (1 | site_code/block) + (1 | site_code:year),
                                  data = total_biomass_precip_block)
summary(mixed_model_total_biomass)


##### BIOMASS BY FUNCTIONAL GROUP ANALYSIS #####

## Calculate precip effect (precip-biomass slope) per site for control and for NPK
slopes_fg_biomass_ctrl <- fg_biomass_precip_block %>%
  filter(trt == "Control") %>%
  group_by(site_code, category) %>%
  group_modify(~ fit_site_lmm(.x, fixed = "precip_z", term = "precip_z")) %>%
  ungroup() %>%
  rename_with(~ paste0("ppt_", .x), -c(site_code, category)) %>%
  rename(ppt_slope = ppt_estimate, ppt_p = ppt_p_value)

slopes_fg_biomass_npk <- fg_biomass_precip_block %>%
  filter(trt == "NPK") %>%
  group_by(site_code, category) %>%
  group_modify(~ fit_site_lmm(.x, fixed = "precip_z", term = "precip_z")) %>%
  ungroup() %>%
  rename_with(~ paste0("ppt_", .x), -c(site_code, category)) %>%
  rename(ppt_slope = ppt_estimate, ppt_p = ppt_p_value)

## Calculate NPK effect (log response ratio) per site: trt coefficient
lrr_fg_biomass <- fg_biomass_precip_block %>%
  group_by(site_code, category) %>%
  group_modify(~ fit_site_lmm(.x, fixed = "trt", term = "trtNPK")) %>%
  ungroup() %>%
  rename_with(~ paste0("npk_", .x), -c(site_code, category)) %>%
  rename(npk_lrr = npk_estimate, npk_p = npk_p_value, )

# Check for singularity
table(slopes_fg_biomass_ctrl$ppt_singular, useNA = "ifany")
table(slopes_fg_biomass_npk$ppt_singular, useNA = "ifany")
table(lrr_fg_biomass$npk_singular, useNA = "ifany")

## Join, classify and tally
slopes_lrr_fg_biomass_ctrl <- left_join(slopes_fg_biomass_ctrl, lrr_fg_biomass,
                                   by = c("site_code", "category")) %>%
  mutate(
    precip_response = classify_slope(ppt_ci_low, ppt_ci_high),
    npk_response    = classify_slope(npk_ci_low, npk_ci_high),
    quadrant = paste(npk_response, precip_response, sep = " / ")
  )

slopes_lrr_fg_biomass_npk <- left_join(slopes_fg_biomass_npk, lrr_fg_biomass,
                                        by = c("site_code", "category")) %>%
  mutate(
    precip_response = classify_slope(ppt_ci_low, ppt_ci_high),
    npk_response    = classify_slope(npk_ci_low, npk_ci_high),
    quadrant = paste(npk_response, precip_response, sep = " / ")
  )

response_counts_fg_biomass_ctrl <- slopes_lrr_fg_biomass_ctrl %>%
  filter(!is.na(precip_response), !is.na(npk_response)) %>%
  count(category, npk_response, precip_response, name = "n_sites") %>%
  complete(category,
           npk_response    = c("positive", "none", "negative"),
           precip_response = c("positive", "none", "negative"),
           fill = list(n_sites = 0)) %>%
  mutate(quadrant = paste(npk_response, precip_response, sep = " / "))

print(response_counts_fg_biomass_ctrl, n=36)

quadrant_counts_fg_biomass_ctrl <- slopes_lrr_fg_biomass_ctrl %>%
  filter(!is.na(precip_response), !is.na(npk_response)) %>%
  count(category, quadrant, name = "n_sites")

print(quadrant_counts_fg_biomass_ctrl)

response_counts_fg_biomass_npk <- slopes_lrr_fg_biomass_npk %>%
  filter(!is.na(precip_response), !is.na(npk_response)) %>%
  count(category, npk_response, precip_response, name = "n_sites") %>%
  complete(category,
           npk_response    = c("positive", "none", "negative"),
           precip_response = c("positive", "none", "negative"),
           fill = list(n_sites = 0)) %>%
  mutate(quadrant = paste(npk_response, precip_response, sep = " / "))

print(response_counts_fg_biomass_npk, n=36)

quadrant_counts_fg_biomass_npk <- slopes_lrr_fg_biomass_npk %>%
  filter(!is.na(precip_response), !is.na(npk_response)) %>%
  count(category, quadrant, name = "n_sites")

print(quadrant_counts_fg_biomass_npk, n=22)

## Analyze
model_fg_biomass_ctrl <- lmer(ppt_slope ~ npk_lrr + category + (1 | site_code),
                         data = slopes_lrr_fg_biomass_ctrl)
summary(model_fg_biomass_ctrl)

model_fg_biomass_npk <- lmer(ppt_slope ~ npk_lrr + category + (1 | site_code),
                              data = slopes_lrr_fg_biomass_npk)
summary(model_fg_biomass_npk)

fg_model_list_ctrl <- slopes_lrr_fg_biomass_ctrl %>%
  split(.$category) %>%
  map(~ lm(ppt_slope ~ npk_lrr, data = .x))

map(fg_model_list_ctrl, summary)

fg_model_list_npk <- slopes_lrr_fg_biomass_npk %>%
  split(.$category) %>%
  map(~ lm(ppt_slope ~ npk_lrr, data = .x))

map(fg_model_list_npk, summary)


##### SPECIES COVER ANALYSIS #####

## How many site x species models will be fit?
n_distinct(cover_precip_block %>% distinct(site_code, Taxon))

## Precip slope per site x species for Control plots and for NPK plots
slopes_species_cover_ctrl <- cover_precip_block %>%
  filter(trt == "Control") %>%
  group_by(site_code, Taxon) %>%
  group_modify(~ fit_site_lmm(.x, fixed = "precip_z", term = "precip_z",
                              response = "log_cover")) %>%
  ungroup() %>%
  rename_with(~ paste0("ppt_", .x), -c(site_code, Taxon)) %>%
  rename(ppt_slope = ppt_estimate, ppt_p = ppt_p_value)

warnings()

slopes_species_cover_npk <- cover_precip_block %>%
  filter(trt == "NPK") %>%
  group_by(site_code, Taxon) %>%
  group_modify(~ fit_site_lmm(.x, fixed = "precip_z", term = "precip_z",
                              response = "log_cover")) %>%
  ungroup() %>%
  rename_with(~ paste0("ppt_", .x), -c(site_code, Taxon)) %>%
  rename(ppt_slope = ppt_estimate, ppt_p = ppt_p_value)

warnings()

## NPK effect (log response ratio) per site x species
lrr_species_cover <- cover_precip_block %>%
  group_by(site_code, Taxon) %>%
  group_modify(~ fit_site_lmm(.x, fixed = "trt", term = "trtNPK",
                              response = "log_cover")) %>%
  ungroup() %>%
  transmute(site_code, Taxon,
            npk_lrr = estimate, npk_se = se,
            npk_ci_low = ci_low, npk_ci_high = ci_high,
            npk_df = df, npk_singular = singular)

## Check for singularity
table(slopes_species_cover_ctrl$ppt_singular, useNA = "ifany")
table(slopes_species_cover_npk$ppt_singular, useNA = "ifany")
table(lrr_species_cover$npk_singular, useNA = "ifany")

## Join, classify, tally
slopes_lrr_species_cover_ctrl <- left_join(slopes_species_cover_ctrl, lrr_species_cover,
                                      by = c("site_code", "Taxon")) %>%
  left_join(distinct(cover_precip_block, site_code, Taxon, functional_group),
            by = c("site_code", "Taxon")) %>%
  mutate(
    precip_response = classify_slope(ppt_ci_low, ppt_ci_high),
    npk_response    = classify_slope(npk_ci_low, npk_ci_high),
    quadrant = paste(npk_response, precip_response, sep = " / ")
  )

slopes_lrr_species_cover_npk <- left_join(slopes_species_cover_npk, lrr_species_cover,
                                           by = c("site_code", "Taxon")) %>%
  left_join(distinct(cover_precip_block, site_code, Taxon, functional_group),
            by = c("site_code", "Taxon")) %>%
  mutate(
    precip_response = classify_slope(ppt_ci_low, ppt_ci_high),
    npk_response    = classify_slope(npk_ci_low, npk_ci_high),
    quadrant = paste(npk_response, precip_response, sep = " / ")
  )

response_counts_species_cover_ctrl <- slopes_lrr_species_cover_ctrl %>%
  filter(!is.na(precip_response), !is.na(npk_response)) %>%
  count(npk_response, precip_response, name = "n_site_species") %>%
  complete(npk_response    = c("positive", "none", "negative"),
           precip_response = c("positive", "none", "negative"),
           fill = list(n_site_species = 0))  %>%
  mutate(quadrant = paste(npk_response, precip_response, sep = " / "))

print(response_counts_species_cover_ctrl)

quadrant_counts_species_cover_ctrl <- slopes_lrr_species_cover_ctrl %>%
  filter(!is.na(precip_response), !is.na(npk_response)) %>%
  count(quadrant, name = "n_site_species")

print(quadrant_counts_species_cover_ctrl)

response_counts_fg_cover_ctrl <- slopes_lrr_species_cover_ctrl %>%
  filter(!is.na(precip_response), !is.na(npk_response)) %>%
  count(functional_group, npk_response, precip_response, name = "n_site_species") %>%
  complete(functional_group,
           npk_response    = c("positive", "none", "negative"),
           precip_response = c("positive", "none", "negative"),
           fill = list(n_site_species = 0)) %>%
  mutate(quadrant = paste(npk_response, precip_response, sep = " / "))

print(response_counts_fg_cover_ctrl, n=36)

quadrant_counts_fg_cover_ctrl <- slopes_lrr_species_cover_ctrl %>%
  filter(!is.na(precip_response), !is.na(npk_response)) %>%
  count(functional_group, quadrant, name = "n_site_species")

print(quadrant_counts_fg_cover_ctrl, n=25)

response_counts_species_cover_npk <- slopes_lrr_species_cover_npk %>%
  filter(!is.na(precip_response), !is.na(npk_response)) %>%
  count(npk_response, precip_response, name = "n_site_species") %>%
  complete(npk_response    = c("positive", "none", "negative"),
           precip_response = c("positive", "none", "negative"),
           fill = list(n_site_species = 0))  %>%
  mutate(quadrant = paste(npk_response, precip_response, sep = " / "))

print(response_counts_species_cover_npk)

quadrant_counts_species_cover_npk <- slopes_lrr_species_cover_npk %>%
  filter(!is.na(precip_response), !is.na(npk_response)) %>%
  count(quadrant, name = "n_site_species")

print(quadrant_counts_species_cover_npk)

response_counts_fg_cover_npk <- slopes_lrr_species_cover_npk %>%
  filter(!is.na(precip_response), !is.na(npk_response)) %>%
  count(functional_group, npk_response, precip_response, name = "n_site_species") %>%
  complete(functional_group,
           npk_response    = c("positive", "none", "negative"),
           precip_response = c("positive", "none", "negative"),
           fill = list(n_site_species = 0)) %>%
  mutate(quadrant = paste(npk_response, precip_response, sep = " / "))

print(response_counts_fg_cover_npk, n=36)

quadrant_counts_fg_cover_npk <- slopes_lrr_species_cover_npk %>%
  filter(!is.na(precip_response), !is.na(npk_response)) %>%
  count(functional_group, quadrant, name = "n_site_species")

print(quadrant_counts_fg_cover_npk, n=28)


## Analyze
model_species_cover_ctrl <- lmer(ppt_slope ~ npk_lrr + (1 | site_code) + (1 | Taxon),
                            data = slopes_lrr_species_cover_ctrl)
summary(model_species_cover_ctrl)

model_species_cover_npk <- lmer(ppt_slope ~ npk_lrr + (1 | site_code) + (1 | Taxon),
                                 data = slopes_lrr_species_cover_npk)
summary(model_species_cover_npk)

fg_cover_model_list_ctrl <- slopes_lrr_species_cover_ctrl %>%
  split(.$functional_group) %>%
  map(~ lm(ppt_slope ~ npk_lrr, data = .x))

map(fg_cover_model_list_ctrl, summary)

fg_cover_model_list_npk <- slopes_lrr_species_cover_npk %>%
  split(.$functional_group) %>%
  map(~ lm(ppt_slope ~ npk_lrr, data = .x))

map(fg_cover_model_list_npk, summary)

site_cover_model_list_ctrl <- slopes_lrr_species_cover_ctrl %>%
  split(.$site_code) %>%
  map(~ lm(ppt_slope ~ npk_lrr, data = .x))

site_cover_results_ctrl <- site_cover_model_list_ctrl %>%
  map_dfr(tidy, .id = "site_code")

site_cover_slopes_ctrl <- site_cover_results_ctrl %>%
  filter(term == "npk_lrr")

View(site_cover_slopes_ctrl)
print(site_cover_slopes_ctrl, n=36)

site_cover_model_list_npk <- slopes_lrr_species_cover_npk %>%
  split(.$site_code) %>%
  map(~ lm(ppt_slope ~ npk_lrr, data = .x))

site_cover_results_npk <- site_cover_model_list_npk %>%
  map_dfr(tidy, .id = "site_code")

site_cover_slopes_npk <- site_cover_results_npk %>%
  filter(term == "npk_lrr")

View(site_cover_slopes_npk)
print(site_cover_slopes_ctrl, n=36)


### Common Species Analysis ###

sp_site_counts_ctrl <- slopes_lrr_species_cover_ctrl %>%
  group_by(Taxon) %>% 
  summarise(site_count = n_distinct(site_code), .groups = 'drop')

sp_site_6_ctrl <- sp_site_counts_ctrl %>%
  filter(site_count >= 6)

slopes_lrr_species_cover6_ctrl <- slopes_lrr_species_cover_ctrl %>%
  filter(Taxon %in% sp_site_6_ctrl$Taxon)

cover6_model_list_ctrl <- slopes_lrr_species_cover6_ctrl %>%
  split(.$Taxon) %>%
  map(~ lm(ppt_slope ~ npk_lrr, data = .x))

map(cover6_model_list_ctrl, summary)


sp_site_counts_npk <- slopes_lrr_species_cover_npk %>%
  group_by(Taxon) %>% 
  summarise(site_count = n_distinct(site_code), .groups = 'drop')

sp_site_6_npk <- sp_site_counts_npk %>%
  filter(site_count >= 6)

slopes_lrr_species_cover6_npk <- slopes_lrr_species_cover_npk %>%
  filter(Taxon %in% sp_site_6_npk$Taxon)

cover6_model_list_npk <- slopes_lrr_species_cover6_npk %>%
  split(.$Taxon) %>%
  map(~ lm(ppt_slope ~ npk_lrr, data = .x))

map(cover6_model_list_npk, summary)


### Species Traits Analysis ###

# Fill in missing ps_path info

ps_path_fill <- subset(cover_trait_cols, is.na(ps_path) == TRUE) %>%
  separate(Taxon, into = c("Genus", "Species"), remove = FALSE, sep = " ")

ps_path_fill$ps_path <- ifelse(ps_path_fill$Genus == "BINERTIA" | ps_path_fill$Genus == "TIDESTROMIA" | ps_path_fill$Genus == "PECTIS" | 
                         ps_path_fill$Genus == "EUPLOCA" | ps_path_fill$Genus == "BULBOSTYLIS" | ps_path_fill$Genus == "CYPERUS" | 
                         ps_path_fill$Genus == "FIMBRISTYLIS" | ps_path_fill$Genus == "CHAMAESYCE" | ps_path_fill$Genus == "ALLIONIA" | 
                         ps_path_fill$Genus == "CALLIGONUM" | ps_path_fill$Genus == "PORTULACA" | ps_path_fill$Genus == "EUPHORBIA", 
                       "C4",
                       ifelse(ps_path_fill$Genus == "BELAPHARIS" | ps_path_fill$Genus == "AERVA" | ps_path_fill$Genus == "ALTERNANTHERA" | 
                                ps_path_fill$Genus == "ATRIPLEX" | ps_path_fill$Genus == "SUAEDA" | ps_path_fill$Genus == "TECTICORNIA" | 
                                ps_path_fill$Genus == "FLAVERIA" | ps_path_fill$Genus == "POLYCARPOREA" | ps_path_fill$Genus == "ELEOCHARS" | 
                                ps_path_fill$Genus == "RHYNCHOSPORA" | ps_path_fill$Genus == "EUPHORBIA" | ps_path_fill$Genus == "MOLLUGO" | 
                                ps_path_fill$Genus == "BOERHAVIA" | ps_path_fill$Genus == "BASSIA" | ps_path_fill$Family == "Poaceae", NA,
                              "C3"))
unique(ps_path_fill$ps_path)

ps_path_fill <- ps_path_fill %>%
  rename(ps_path2 = ps_path)

names(ps_path_fill)
colnames(cover_trait_cols)

cover_trait_cols <- cover_trait_cols %>% 
  left_join(ps_path_fill, by = c("site_code", "Taxon", "Family", "functional_group",
                                 "local_provenance", "local_lifeform", "local_lifespan"))

cover_trait_cols <- cover_trait_cols %>%
  mutate(ps_path2 = if_else(is.na(ps_path2) & !is.na(ps_path), ps_path, ps_path2))


# Join species traits to slopes_lrr dataframes

slopes_lrr_species_cover_ctrl <- slopes_lrr_species_cover_ctrl %>%
  inner_join(cover_trait_cols, by=c("site_code", "Taxon", "functional_group"))

slopes_lrr_species_cover_npk <- slopes_lrr_species_cover_npk %>%
  inner_join(cover_trait_cols, by=c("site_code", "Taxon", "functional_group"))


annual_species_model_ctrl <- lm(ppt_slope ~ npk_lrr, 
                           data = subset(slopes_lrr_species_cover_ctrl, 
                                         local_lifespan == "ANNUAL"))
summary(annual_species_model_ctrl)

perennial_species_model_ctrl <- lm(ppt_slope ~ npk_lrr, 
                           data = subset(slopes_lrr_species_cover_ctrl, 
                                         local_lifespan == "PERENNIAL"))
summary(perennial_species_model_ctrl)

annual_species_model_npk <- lm(ppt_slope ~ npk_lrr, 
                                data = subset(slopes_lrr_species_cover_npk, 
                                              local_lifespan == "ANNUAL"))
summary(annual_species_model_npk)

perennial_species_model_npk <- lm(ppt_slope ~ npk_lrr, 
                                   data = subset(slopes_lrr_species_cover_npk, 
                                                 local_lifespan == "PERENNIAL"))
summary(perennial_species_model_npk)


int_species_model_ctrl <- lm(ppt_slope ~ npk_lrr, 
                           data = subset(slopes_lrr_species_cover_ctrl, 
                                         local_provenance == "INT"))
summary(int_species_model_ctrl)

nat_species_model_ctrl <- lm(ppt_slope ~ npk_lrr, 
                             data = subset(slopes_lrr_species_cover_ctrl, 
                                           local_provenance == "NAT"))
summary(nat_species_model_ctrl)


int_species_model_npk <- lm(ppt_slope ~ npk_lrr, 
                             data = subset(slopes_lrr_species_cover_npk, 
                                           local_provenance == "INT"))
summary(int_species_model_npk)

nat_species_model_npk <- lm(ppt_slope ~ npk_lrr, 
                             data = subset(slopes_lrr_species_cover_npk, 
                                           local_provenance == "NAT"))
summary(nat_species_model_npk)


c3_species_model_ctrl <- lm(ppt_slope ~ npk_lrr, 
                             data = subset(slopes_lrr_species_cover_ctrl, 
                                           ps_path2 == "C3"))
summary(c3_species_model_ctrl)

c4_species_model_ctrl <- lm(ppt_slope ~ npk_lrr, 
                             data = subset(slopes_lrr_species_cover_ctrl, 
                                           ps_path2 == "C4"))
summary(c4_species_model_ctrl)


c3_species_model_npk <- lm(ppt_slope ~ npk_lrr, 
                            data = subset(slopes_lrr_species_cover_npk, 
                                          ps_path2 == "C3"))
summary(c3_species_model_npk)

c4_species_model_npk <- lm(ppt_slope ~ npk_lrr, 
                            data = subset(slopes_lrr_species_cover_npk, 
                                          ps_path2 == "C4"))
summary(c4_species_model_npk)


## Family-level Analysis ##

fam_site_counts_ctrl <- slopes_lrr_species_cover_ctrl %>%
  group_by(Family) %>% 
  summarise(site_count = n_distinct(site_code), .groups = 'drop')

fam_site_6 <- fam_site_counts_ctrl %>%
  filter(site_count >= 6)

slopes_lrr_species_cover_fam6_ctrl <- slopes_lrr_species_cover_ctrl %>%
  filter(Family %in% fam_site_6$Family)

fam_site_counts_npk <- slopes_lrr_species_cover_npk %>%
  group_by(Family) %>% 
  summarise(site_count = n_distinct(site_code), .groups = 'drop')

fam_site_6_npk <- fam_site_counts_npk %>%
  filter(site_count >= 6)

slopes_lrr_species_cover_fam6_npk <- slopes_lrr_species_cover_npk %>%
  filter(Family %in% fam_site_6_npk$Family)


fam_cover_model_list_ctrl <- slopes_lrr_species_cover_fam6_ctrl %>%
  split(.$Family) %>%
  map(~ lm(ppt_slope ~ npk_lrr, data = .x))

map(fam_cover_model_list_ctrl, summary)

fam_cover_model_list_npk <- slopes_lrr_species_cover_fam6_npk %>%
  split(.$Family) %>%
  map(~ lm(ppt_slope ~ npk_lrr, data = .x))

map(fam_cover_model_list_npk, summary)


##### Graphing #####

quad_levels <- c(
  "positive / positive", "positive / none", "positive / negative",
  "none / positive",     "none / none",     "none / negative",
  "negative / positive", "negative / none", "negative / negative"
)

quad_colors <- setNames(
  c("#D55E00","#0072B2","#F0E442","#E69F00","#999999","#CC79A7","#6A3D9A","#009E73","#56B4E9"),
  quad_levels
)

scale_quad_color <- scale_color_manual(values = quad_colors, breaks = quad_levels)
scale_quad_fill  <- scale_fill_manual(values = quad_colors,  breaks = quad_levels)

## Total Biomass Graphs

graph_total_biomass_ctrl <- ggplot(slopes_lrr_total_biomass_ctrl,
                              aes(x = npk_lrr, y = ppt_slope, color = quadrant)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(limits = c(-0.3, 1.1)) +
  scale_quad_color +
  labs(x = "Response to NPK",
       y = "Response to precipitation \n Control Plots",
       color = "NPK Response / \nPrecip Resoponse") +
  theme_bw(base_size = 14)

graph_total_biomass_ctrl

graph_total_biomass_npk <- ggplot(slopes_lrr_total_biomass_npk,
                                   aes(x = npk_lrr, y = ppt_slope, color = quadrant)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(limits = c(-0.3, 1.1)) +
  scale_quad_color +
  labs(x = "Response to NPK",
       y = "Response to precipitation \n NPK Plots",
       color = "NPK Response / \nPrecip Resoponse") +
  theme_bw(base_size = 14)

graph_total_biomass_npk


graph_total_biomass_both <- ggarrange(
  graph_total_biomass_ctrl + rremove("xlab"),
  graph_total_biomass_npk,
  ncol = 1,
  common.legend = TRUE,
  legend = "right",
  align = 'hv'
)

graph_total_biomass_both


bargraph_total_biomass_ctrl <- ggplot(response_counts_total_biomass_ctrl,
                                 aes(x = npk_response, y = n_sites, fill = quadrant)) +
  geom_col() +
  scale_quad_fill +
  labs(x = "NPK response", y = "Number of sites",
       fill = "NPK Response / \nPrecip Resoponse") +
  theme_bw(base_size = 14)

bargraph_total_biomass_ctrl

bargraph_quadrant_total_biomass_ctrl <- ggplot(quadrant_counts_total_biomass_ctrl,
                                          aes(x = quadrant, y = n_sites, fill = quadrant)) +
  geom_col() +
  scale_quad_fill +
  guides(fill = "none") +
  labs(x = "NPK Response / \nPrecip Resoponse", y = "Number of sites") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

bargraph_quadrant_total_biomass_ctrl


## Biomass FG Graphs

graph_fg_biomass_ctrl <- ggplot(slopes_lrr_fg_biomass_ctrl,
                           aes(x = npk_lrr, y = ppt_slope, color = quadrant)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_y_continuous(limits = c(-2, 3.1)) +
  scale_quad_color +
  facet_wrap(~ category) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Response to NPK",
       y = "Response to precipitation \nControl Plots",
       color = "NPK Response / \nPrecip Resoponse") +
  theme_bw(base_size = 14)

graph_fg_biomass_ctrl

graph_fg_biomass_npk <- ggplot(slopes_lrr_fg_biomass_npk,
                                aes(x = npk_lrr, y = ppt_slope, color = quadrant)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_y_continuous(limits = c(-2, 3.1)) +
  scale_quad_color +
  facet_wrap(~ category) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Response to NPK",
       y = "Response to precipitation \nNPK Plots",
       color = "NPK Response / \nPrecip Resoponse") +
  theme_bw(base_size = 14)

graph_fg_biomass_npk

graph_fg_biomass_both <- ggarrange(
  graph_fg_biomass_ctrl + rremove("xlab"),
  graph_fg_biomass_npk,
  ncol = 1,
  common.legend = TRUE,
  legend = "right",
  align = 'hv'
)

graph_fg_biomass_both


bargraph_fg_biomass_ctrl <- ggplot(response_counts_fg_biomass_ctrl,
                              aes(x = npk_response, y = n_sites, fill = quadrant)) +
  geom_col() +
  scale_quad_fill +
  facet_wrap(~ category) +
  labs(x = "NPK response", y = "Number of site-FG combinations",
       fill = "NPK Response / \nPrecip Resoponse") +
  theme_bw(base_size = 14)

bargraph_fg_biomass_ctrl

bargraph_fg_biomass_npk <- ggplot(response_counts_fg_biomass_npk,
                                   aes(x = npk_response, y = n_sites, fill = quadrant)) +
  geom_col() +
  scale_quad_fill +
  facet_wrap(~ category) +
  labs(x = "NPK response", y = "Number of site-FG combinations",
       fill = "NPK Response / \nPrecip Resoponse") +
  theme_bw(base_size = 14)

bargraph_fg_biomass_npk

bargraph_fg_biomass_both <- ggarrange(
  bargraph_fg_biomass_ctrl + rremove("xlab"),
  bargraph_fg_biomass_npk,
  ncol = 1,
  common.legend = TRUE,
  legend = "right",
  align = 'hv'
)

bargraph_fg_biomass_both


bargraph_quadrant_fg_biomass_ctrl <- ggplot(quadrant_counts_fg_biomass_ctrl,
                                       aes(x = quadrant, y = n_sites, fill = quadrant)) +
  geom_col() +
  scale_quad_fill +
  guides(fill = "none") +
  facet_wrap(~ category) +
  labs(x = "NPK Response / Precip Resoponse", y = "Number of sites") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

bargraph_quadrant_fg_biomass_ctrl

bargraph_quadrant_fg_biomass_npk <- ggplot(quadrant_counts_fg_biomass_npk,
                                            aes(x = quadrant, y = n_sites, fill = quadrant)) +
  geom_col() +
  scale_quad_fill +
  guides(fill = "none") +
  facet_wrap(~ category) +
  labs(x = "NPK Response / Precip Resoponse", y = "Number of sites") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

bargraph_quadrant_fg_biomass_npk

bargraph_quadrant_fg_biomass_both <- ggarrange(
  bargraph_quadrant_fg_biomass_ctrl + rremove("xlab"),
  bargraph_quadrant_fg_biomass_npk,
  ncol = 1,
  common.legend = TRUE,
  align = 'hv'
)

bargraph_quadrant_fg_biomass_both


## Species Cover Graphs
plot_species_cover_ctrl <- ggplot(slopes_lrr_species_cover_ctrl,
                             aes(x = npk_lrr, y = ppt_slope, color = quadrant)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_y_continuous(limits = c(-2.8, 4)) +
  scale_quad_color +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Response to NPK",
       y = "Response to precipitation \nControl Plots",
       color = "NPK Response / \nPrecip Resoponse") +
  theme_bw(base_size = 14)

plot_species_cover_ctrl

plot_species_cover_npk <- ggplot(slopes_lrr_species_cover_npk,
                                  aes(x = npk_lrr, y = ppt_slope, color = quadrant)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_y_continuous(limits = c(-2.8, 4)) +
  scale_quad_color +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Response to NPK",
       y = "Response to precipitation \nNPK plots",
       color = "NPK Response / \nPrecip Resoponse") +
  theme_bw(base_size = 14)

plot_species_cover_npk

plot_species_cover_both <- ggarrange(
  plot_species_cover_ctrl + rremove("xlab"),
  plot_species_cover_npk,
  ncol = 1,
  common.legend = TRUE,
  legend = "right",
  align = 'hv'
)

plot_species_cover_both


barplot_species_cover_ctrl <- ggplot(response_counts_species_cover_ctrl,
                              aes(x = npk_response, y = n_site_species, fill = quadrant)) +
  geom_col() +
  scale_quad_fill +
  labs(x = "NPK response", y = "Number of site-species combinations \nControl Plots",
       fill = "NPK Response / \nPrecip Resoponse") +
  theme_bw(base_size = 14)

barplot_species_cover_ctrl

barplot_species_cover_npk <- ggplot(response_counts_species_cover_npk,
                                      aes(x = npk_response, y = n_site_species, fill = quadrant)) +
  geom_col() +
  scale_quad_fill +
  labs(x = "NPK response", y = "Number of site-species combinations \nNPK Plots",
       fill = "NPK Response / \nPrecip Resoponse") +
  theme_bw(base_size = 14)

barplot_species_cover_npk

barplot_species_cover_both <- ggarrange(
  barplot_species_cover_ctrl + rremove("xlab"),
  barplot_species_cover_npk,
  ncol = 1,
  common.legend = TRUE,
  legend = "right",
  align = 'hv'
)

barplot_species_cover_both


barplot_quadrant_species_cover_ctrl <- ggplot(quadrant_counts_species_cover_ctrl,
                                              aes(x = quadrant, y = n_site_species, fill = quadrant)) +
  geom_col() +
  scale_quad_fill +
  guides(fill = "none") +
  labs(x = "NPK Response / Precip Resoponse", 
       y = "Number of site-species combinations \nControl Plots") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

barplot_quadrant_species_cover_ctrl

barplot_quadrant_species_cover_npk <- ggplot(quadrant_counts_species_cover_npk,
                                       aes(x = quadrant, y = n_site_species, fill = quadrant)) +
  geom_col() +
  scale_quad_fill +
  guides(fill = "none") +
  labs(x = "NPK Response / Precip Resoponse", 
       y = "Number of site-species combinations \nNPK Plots") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

barplot_quadrant_species_cover_npk

barplot_quadrant_species_cover_both <- ggarrange(
  barplot_quadrant_species_cover_ctrl + rremove("xlab"),
  barplot_quadrant_species_cover_npk,
  ncol = 1,
  common.legend = TRUE,
  legend = "right",
  align = 'hv'
)

barplot_quadrant_species_cover_both


## FG Cover Graphs

plot_fg_cover_ctrl <- ggplot(slopes_lrr_species_cover_ctrl,
                                aes(x = npk_lrr, y = ppt_slope, color = quadrant)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_y_continuous(limits = c(-2.8, 4)) +
  scale_quad_color +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Response to NPK",
       y = "Response to precipitation \nControl Plots",
       color = "NPK Response / \nPrecip Resoponse") +
  facet_wrap(~ functional_group) +
  theme_bw(base_size = 14)

plot_fg_cover_ctrl

plot_fg_cover_npk <- ggplot(slopes_lrr_species_cover_npk,
                                     aes(x = npk_lrr, y = ppt_slope, color = quadrant)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_y_continuous(limits = c(-2.8, 4)) +
  scale_quad_color +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Response to NPK",
       y = "Response to precipitation \nNPK Plots",
       color = "NPK Response / \nPrecip Resoponse") +
  facet_wrap(~ functional_group) +
  theme_bw(base_size = 14)

plot_fg_cover_npk

plot_fg_cover_both <- ggarrange(
  plot_species_cover_ctrl + rremove("xlab"),
  plot_species_cover_npk,
  ncol = 1,
  common.legend = TRUE,
  legend = "right",
  align = 'hv'
)

plot_fg_cover_both


barplot_fg_cover_ctrl <- ggplot(response_counts_fg_cover_ctrl,
                                 aes(x = npk_response, y = n_site_species, fill = quadrant)) +
  geom_col() +
  scale_quad_fill +
  facet_wrap(~ functional_group) +
  labs(x = "NPK response", 
       y = "Number of site-species combinations \nControl Plots",
       fill = "NPK Response / \nPrecip Resoponse") +
  theme_bw(base_size = 14)

barplot_fg_cover_ctrl

barplot_fg_cover_npk <- ggplot(response_counts_fg_cover_npk,
                                 aes(x = npk_response, y = n_site_species, fill = quadrant)) +
  geom_col() +
  scale_quad_fill +
  facet_wrap(~ functional_group) +
  labs(x = "NPK response", 
       y = "Number of site-species combinations \nNPK Plots",
       fill = "NPK Response / \nPrecip Resoponse") +
  theme_bw(base_size = 14)

barplot_fg_cover_npk

barplot_fg_cover_both <- ggarrange(
  barplot_fg_cover_ctrl + rremove("xlab"),
  barplot_fg_cover_npk,
  ncol = 1,
  common.legend = TRUE,
  legend = "right",
  align = 'hv'
)

barplot_species_cover_both


barplot_quadrant_fg_cover_ctrl <- ggplot(response_counts_fg_cover_ctrl,
                                          aes(x = quadrant, y = n_site_species, fill = quadrant)) +
  geom_col() +
  scale_quad_fill +
  guides(fill = "none") +
  facet_wrap(~ functional_group) +
  labs(x = "NPK Response / Precip Resoponse", 
       y = "Number of site-species combinations \nControl Plots") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

barplot_quadrant_fg_cover_ctrl

barplot_quadrant_fg_cover_npk <- ggplot(response_counts_fg_cover_npk,
                                         aes(x = quadrant, y = n_site_species, fill = quadrant)) +
  geom_col() +
  scale_quad_fill +
  guides(fill = "none") +
  facet_wrap(~ functional_group) +
  labs(x = "NPK Response / Precip Resoponse", 
       y = "Number of site-species combinations \nControl Plots") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

barplot_quadrant_fg_cover_npk

barplot_quadrant_fg_cover_both <- ggarrange(
  barplot_quadrant_fg_cover_ctrl + rremove("xlab"),
  barplot_quadrant_fg_cover_npk,
  ncol = 1,
  common.legend = TRUE,
  legend = "right",
  align = 'hv'
)

barplot_quadrant_fg_cover_both


## Site Graphs

plot_species_cover_sites_ctrl <- ggplot(slopes_lrr_species_cover_ctrl,
                             aes(x = npk_lrr, y = ppt_slope, color = quadrant)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_quad_color +
  facet_wrap(~ site_code) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Response to NPK",
       y = "Response to precipitation \nControl Plots",
       color = "NPK Response / \nPrecip Resoponse") +
  theme_bw(base_size = 14)

plot_species_cover_sites_ctrl

plot_species_cover_sites_npk <- ggplot(slopes_lrr_species_cover_npk,
                                        aes(x = npk_lrr, y = ppt_slope, color = quadrant)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_quad_color +
  facet_wrap(~ site_code) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Response to NPK",
       y = "Response to precipitation \nNPK Plots",
       color = "NPK Response / \nPrecip Resoponse") +
  theme_bw(base_size = 14)

plot_species_cover_sites_npk


## Common Species Graphs

plot_species_cover_taxon_ctrl <- ggplot(slopes_lrr_species_cover6_ctrl,
                                   aes(x = npk_lrr, y = ppt_slope, color = quadrant)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_quad_color +
  facet_wrap(~ Taxon) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Response to NPK",
       y = "Response to precipitation \nControl Plots",
       color = "NPK Response / \nPrecip Resoponse") +
  theme_bw(base_size = 14)

plot_species_cover_taxon_ctrl

plot_species_cover_taxon_npk <- ggplot(slopes_lrr_species_cover6_npk,
                                        aes(x = npk_lrr, y = ppt_slope, color = quadrant)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_quad_color +
  facet_wrap(~ Taxon) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Response to NPK",
       y = "Response to precipitation \nNPK Plots",
       color = "NPK Response / \nPrecip Resoponse") +
  theme_bw(base_size = 14)

plot_species_cover_taxon_npk


## Other species trait graphs

plot_species_cover_lifespan_ctrl <- ggplot(
  data = subset(slopes_lrr_species_cover_ctrl, local_lifespan %in% c("ANNUAL", "PERENNIAL")),
                                aes(x = npk_lrr, y = ppt_slope, color = quadrant)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_quad_color +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Response to NPK",
       y = "Response to precipitation \nControl Plots",
       color = "NPK Response / \nPrecip Resoponse") +
  facet_wrap(~ local_lifespan) +
  theme_bw(base_size = 14)

plot_species_cover_lifespan_ctrl

plot_species_cover_lifespan_npk <- ggplot(
  data = subset(slopes_lrr_species_cover_npk, local_lifespan %in% c("ANNUAL", "PERENNIAL")),
  aes(x = npk_lrr, y = ppt_slope, color = quadrant)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_quad_color +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Response to NPK",
       y = "Response to precipitation \nNPK Plots",
       color = "NPK Response / \nPrecip Resoponse") +
  facet_wrap(~ local_lifespan) +
  theme_bw(base_size = 14)

plot_species_cover_lifespan_npk

plot_species_cover_lifespan_both <- ggarrange(
  plot_species_cover_lifespan_ctrl + rremove("xlab"),
  plot_species_cover_lifespan_npk,
  ncol = 1,
  common.legend = TRUE,
  legend = "right",
  align = 'hv'
)

plot_species_cover_lifespan_both


plot_species_cover_provenance_ctrl <- ggplot(
  data = subset(slopes_lrr_species_cover_ctrl, local_provenance %in% c("INT", "NAT")),
  aes(x = npk_lrr, y = ppt_slope, color = quadrant)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_quad_color +
  scale_y_continuous(limits = c(-2.8, 4)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Response to NPK",
       y = "Response to precipitation \nControl Plots",
       color = "NPK Response / \nPrecip Resoponse") +
  facet_wrap(~ local_provenance) +
  theme_bw(base_size = 14)

plot_species_cover_provenance_ctrl

plot_species_cover_provenance_npk <- ggplot(
  data = subset(slopes_lrr_species_cover_npk, local_provenance %in% c("INT", "NAT")),
  aes(x = npk_lrr, y = ppt_slope, color = quadrant)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_quad_color +
  scale_y_continuous(limits = c(-2.8, 4)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Response to NPK",
       y = "Response to precipitation \nNPK Plots",
       color = "NPK Response / \nPrecip Resoponse") +
  facet_wrap(~ local_provenance) +
  theme_bw(base_size = 14)

plot_species_cover_provenance_npk

plot_species_cover_provenance_both <- ggarrange(
  plot_species_cover_provenance_ctrl + rremove("xlab"),
  plot_species_cover_provenance_npk,
  ncol = 1,
  common.legend = TRUE,
  legend = "right",
  align = 'hv'
)

plot_species_cover_provenance_both


plot_species_cover_ps_ctrl <- ggplot(
  data = subset(slopes_lrr_species_cover_ctrl, ps_path %in% c("C3", "C4")),
  aes(x = npk_lrr, y = ppt_slope, color = quadrant)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_quad_color +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Response to NPK",
       y = "Response to precipitation \nControl Plots",
       color = "NPK Response / \nPrecip Resoponse") +
  facet_wrap(~ ps_path) +
  theme_bw(base_size = 14)

plot_species_cover_ps_ctrl


plot_species_cover_ps_npk <- ggplot(
  data = subset(slopes_lrr_species_cover_npk, ps_path %in% c("C3", "C4")),
  aes(x = npk_lrr, y = ppt_slope, color = quadrant)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_quad_color +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Response to NPK",
       y = "Response to precipitation \nNPK Plots",
       color = "NPK Response / \nPrecip Resoponse") +
  facet_wrap(~ ps_path) +
  theme_bw(base_size = 14)

plot_species_cover_ps_npk


plot_species_cover_ps_both <- ggarrange(
  plot_species_cover_ps_ctrl + rremove("xlab"),
  plot_species_cover_ps_npk,
  ncol = 1,
  common.legend = TRUE,
  legend = "right",
  align = 'hv'
)

plot_species_cover_ps_both


## Family Graph

plot_species_cover_family_ctrl <- ggplot(
  data = slopes_lrr_species_cover_fam6,
  aes(x = npk_lrr, y = ppt_slope, color = quadrant)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_quad_color +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Response to NPK",
       y = "Response to precipitation \nControl Plots",
       color = "NPK Response / \nPrecip Resoponse") +
  facet_wrap(~ Family) +
  theme_bw(base_size = 14)

plot_species_cover_family_ctrl

plot_species_cover_family_npk <- ggplot(
  data = slopes_lrr_species_cover_fam6_npk,
  aes(x = npk_lrr, y = ppt_slope, color = quadrant)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_quad_color +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Response to NPK",
       y = "Response to precipitation \nNPK Plots",
       color = "NPK Response / \nPrecip Resoponse") +
  facet_wrap(~ Family) +
  theme_bw(base_size = 14)

plot_species_cover_family_npk


