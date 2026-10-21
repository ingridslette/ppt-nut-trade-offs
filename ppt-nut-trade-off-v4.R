
library(tidyverse)
library(lme4)
library(lmerTest)



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

## Calculate precip effect (precip-biomass slope) per site (control plots only)
slopes_total_biomass <- total_biomass_precip_block %>%
  filter(trt == "Control") %>%
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
table(slopes_total_biomass$ppt_singular, useNA = "ifany")
table(lrr_total_biomass$npk_singular, useNA = "ifany")

## Join, classify and tally
slopes_lrr_total_biomass <- left_join(slopes_total_biomass, lrr_total_biomass,
                                      by = "site_code") %>%
  mutate(
    precip_response = classify_slope(ppt_ci_low, ppt_ci_high),
    npk_response    = classify_slope(npk_ci_low, npk_ci_high),
    quadrant = paste(npk_response, precip_response, sep = " / ")
  )

response_counts_total_biomass <- slopes_lrr_total_biomass %>%
  count(npk_response, precip_response, name = "n_sites") %>%
  complete(npk_response    = c("positive", "none", "negative"),
           precip_response = c("positive", "none", "negative"),
           fill = list(n_sites = 0))%>%
  mutate(quadrant = paste(npk_response, precip_response, sep = " / "))

print(response_counts_total_biomass)

quadrant_counts_total_biomass <- slopes_lrr_total_biomass %>%
  count(quadrant, name = "n_sites")

print(quadrant_counts_total_biomass)

## Analyze
model_total_biomass <- lm(ppt_slope ~ npk_lrr, data = slopes_lrr_total_biomass)
summary(model_total_biomass)

mixed_model_total_biomass <- lmer(log_mass ~ trt * precip_z +
                                    (1 | site_code/block) + (1 | site_code:year),
                                  data = total_biomass_precip_block)
summary(mixed_model_total_biomass)


##### BIOMASS BY FUNCTIONAL GROUP ANALYSIS #####

## Calculate precip effect (precip-biomass slope) per site (control plots only)
slopes_fg_biomass <- fg_biomass_precip_block %>%
  filter(trt == "Control") %>%
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
table(slopes_fg_biomass$ppt_singular, useNA = "ifany")
table(lrr_fg_biomass$npk_singular, useNA = "ifany")

## Join, classify and tally
slopes_lrr_fg_biomass <- left_join(slopes_fg_biomass, lrr_fg_biomass,
                                   by = c("site_code", "category")) %>%
  mutate(
    precip_response = classify_slope(ppt_ci_low, ppt_ci_high),
    npk_response    = classify_slope(npk_ci_low, npk_ci_high),
    quadrant = paste(npk_response, precip_response, sep = " / ")
  )

response_counts_fg_biomass <- slopes_lrr_fg_biomass %>%
  filter(!is.na(precip_response), !is.na(npk_response)) %>%
  count(category, npk_response, precip_response, name = "n_sites") %>%
  complete(category,
           npk_response    = c("positive", "none", "negative"),
           precip_response = c("positive", "none", "negative"),
           fill = list(n_sites = 0)) %>%
  mutate(quadrant = paste(npk_response, precip_response, sep = " / "))

print(response_counts_fg_biomass, n=36)

quadrant_counts_fg_biomass <- slopes_lrr_fg_biomass %>%
  filter(!is.na(precip_response), !is.na(npk_response)) %>%
  count(category, quadrant, name = "n_sites")

print(quadrant_counts_fg_biomass)

## Analyze
model_fg_biomass <- lmer(ppt_slope ~ npk_lrr + category + (1 | site_code),
                         data = slopes_lrr_fg_biomass)
summary(model_fg_biomass)

model_fg_biomass_x <- lmer(ppt_slope ~ npk_lrr * category + (1 | site_code),
                           data = slopes_lrr_fg_biomass)
summary(model_fg_biomass_x)

fg_model_list <- slopes_lrr_fg_biomass %>%
  split(.$category) %>%
  map(~ lm(ppt_slope ~ npk_lrr, data = .x))

map(fg_model_list, summary)


##### SPECIES COVER ANALYSIS #####

## How many site x species models will be fit?
n_distinct(cover_precip_block %>% distinct(site_code, Taxon))

## Precip slope per site x species (control plots only)
slopes_species_cover <- cover_precip_block %>%
  filter(trt == "Control") %>%
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
table(slopes_species_cover$ppt_singular, useNA = "ifany")
table(lrr_species_cover$npk_singular, useNA = "ifany")

## Join, classify, tally
slopes_lrr_species_cover <- left_join(slopes_species_cover, lrr_species_cover,
                                      by = c("site_code", "Taxon")) %>%
  left_join(distinct(cover_precip_block, site_code, Taxon, functional_group),
            by = c("site_code", "Taxon")) %>%
  mutate(
    precip_response = classify_slope(ppt_ci_low, ppt_ci_high),
    npk_response    = classify_slope(npk_ci_low, npk_ci_high),
    quadrant = paste(npk_response, precip_response, sep = " / ")
  )

response_counts_species_cover <- slopes_lrr_species_cover %>%
  filter(!is.na(precip_response), !is.na(npk_response)) %>%
  count(npk_response, precip_response, name = "n_site_species") %>%
  complete(npk_response    = c("positive", "none", "negative"),
           precip_response = c("positive", "none", "negative"),
           fill = list(n_site_species = 0))  %>%
  mutate(quadrant = paste(npk_response, precip_response, sep = " / "))

print(response_counts_species_cover)

quadrant_counts_species_cover <- slopes_lrr_species_cover %>%
  filter(!is.na(precip_response), !is.na(npk_response)) %>%
  count(quadrant, name = "n_site_species")

print(quadrant_counts_species_cover)

response_counts_fg_cover <- slopes_lrr_species_cover %>%
  filter(!is.na(precip_response), !is.na(npk_response)) %>%
  count(functional_group, npk_response, precip_response, name = "n_site_species") %>%
  complete(functional_group,
           npk_response    = c("positive", "none", "negative"),
           precip_response = c("positive", "none", "negative"),
           fill = list(n_site_species = 0)) %>%
  mutate(quadrant = paste(npk_response, precip_response, sep = " / "))

print(response_counts_fg_cover, n=36)

quadrant_counts_fg_cover <- slopes_lrr_species_cover %>%
  filter(!is.na(precip_response), !is.na(npk_response)) %>%
  count(functional_group, quadrant, name = "n_site_species")

print(quadrant_counts_fg_cover)


## Analyze
model_species_cover <- lmer(ppt_slope ~ npk_lrr + (1 | site_code) + (1 | Taxon),
                            data = slopes_lrr_species_cover)
summary(model_species_cover)

fg_cover_model_list <- slopes_lrr_species_cover %>%
  split(.$functional_group) %>%
  map(~ lm(ppt_slope ~ npk_lrr, data = .x))

map(fg_cover_model_list, summary)


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

graph_total_biomass <- ggplot(slopes_lrr_total_biomass,
                              aes(x = npk_lrr, y = ppt_slope, color = quadrant)) +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_quad_color +
  labs(x = "Response to NPK",
       y = "Response to precipitation",
       color = "NPK Response / \nPrecip Resoponse") +
  theme_bw(base_size = 14)

graph_total_biomass

bargraph_total_biomass <- ggplot(response_counts_total_biomass,
                                 aes(x = npk_response, y = n_sites, fill = quadrant)) +
  geom_col() +
  scale_quad_fill +
  labs(x = "NPK response", y = "Number of sites",
       fill = "NPK Response / \nPrecip Resoponse") +
  theme_bw(base_size = 14)

bargraph_total_biomass

bargraph_quadrant_total_biomass <- ggplot(quadrant_counts_total_biomass,
                                          aes(x = quadrant, y = n_sites, fill = quadrant)) +
  geom_col() +
  scale_quad_fill +
  guides(fill = "none") +
  labs(x = "NPK Response / \nPrecip Resoponse", y = "Number of sites") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

bargraph_quadrant_total_biomass


## Biomass FG Graphs

graph_fg_biomass <- ggplot(slopes_lrr_fg_biomass,
                           aes(x = npk_lrr, y = ppt_slope, color = quadrant)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_quad_color +
  facet_wrap(~ category) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Response to NPK",
       y = "Response to precipitation",
       color = "NPK Response / \nPrecip Resoponse") +
  theme_bw(base_size = 14)

graph_fg_biomass

bargraph_fg_biomass <- ggplot(response_counts_fg_biomass,
                              aes(x = npk_response, y = n_sites, fill = quadrant)) +
  geom_col() +
  scale_quad_fill +
  facet_wrap(~ category) +
  labs(x = "NPK response", y = "Number of site-FG combinations",
       fill = "NPK Response / \nPrecip Resoponse") +
  theme_bw(base_size = 14)

bargraph_fg_biomass

bargraph_quadrant_fg_biomass <- ggplot(quadrant_counts_fg_biomass,
                                       aes(x = quadrant, y = n_sites, fill = quadrant)) +
  geom_col() +
  scale_quad_fill +
  guides(fill = "none") +
  facet_wrap(~ category) +
  labs(x = "NPK Response / Precip Resoponse", y = "Number of sites") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

bargraph_quadrant_fg_biomass


## Species Cover Graphs

plot_species_cover <- ggplot(slopes_lrr_species_cover,
                             aes(x = npk_lrr, y = ppt_slope, color = quadrant)) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_quad_color +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Response to NPK",
       y = "Response to precipitation",
       color = "NPK Response / \nPrecip Resoponse") +
  theme_bw(base_size = 14)

plot_species_cover

bargraph_species_cover <- ggplot(response_counts_species_cover,
                              aes(x = npk_response, y = n_site_species, fill = quadrant)) +
  geom_col() +
  scale_quad_fill +
  labs(x = "NPK response", y = "Number of site-species combinations",
       fill = "NPK Response / \nPrecip Resoponse") +
  theme_bw(base_size = 14)

bargraph_species_cover

bargraph_quadrant_species_cover <- ggplot(quadrant_counts_species_cover,
                                       aes(x = quadrant, y = n_site_species, fill = quadrant)) +
  geom_col() +
  scale_quad_fill +
  guides(fill = "none") +
  labs(x = "NPK Response / Precip Resoponse", y = "Number of site-species combinations") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

bargraph_quadrant_species_cover


## FG Cover Graphs

plot_species_cover_fg <- ggplot(slopes_lrr_species_cover,
                                aes(x = npk_lrr, y = ppt_slope, color = quadrant)) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_quad_color +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Response to NPK",
       y = "Response to precipitation",
       color = "NPK Response / \nPrecip Resoponse") +
  facet_wrap(~ functional_group) +
  theme_bw(base_size = 14)

plot_species_cover_fg

bargraph_fg_cover <- ggplot(response_counts_fg_cover,
                                 aes(x = npk_response, y = n_site_species, fill = quadrant)) +
  geom_col() +
  scale_quad_fill +
  facet_wrap(~ functional_group) +
  labs(x = "NPK response", y = "Number of site-species combinations",
       fill = "NPK Response / \nPrecip Resoponse") +
  theme_bw(base_size = 14)

bargraph_fg_cover

bargraph_quadrant_fg_cover <- ggplot(response_counts_fg_cover,
                                          aes(x = quadrant, y = n_site_species, fill = quadrant)) +
  geom_col() +
  scale_quad_fill +
  guides(fill = "none") +
  facet_wrap(~ functional_group) +
  labs(x = "NPK Response / Precip Resoponse", y = "Number of site-species combinations") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

bargraph_quadrant_fg_cover


## Site Graphs

