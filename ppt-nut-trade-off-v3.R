

library(tidyverse)
library(lme4)
library(lmerTest)

## LOAD AND FILTER BIOMASS DATA
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


## LOAD AND FILTER PRECIPITATION DATA
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


##############################################################################
# NutNet: biomass sensitivity to precipitation vs. NPK fertilization
# Site-level mixed-model version
#
# KEEP your existing code above this point, unchanged: package loading,
# biomass loading/filtering, precip loading/filtering, the join, and the
# 6-year and +/-1 SD filters. This script starts where your
# "TOTAL BIOMASS ANALYSIS" section began and assumes `biomass_precip` exists
# (columns used: site_code, site_name, year, year_trt, trt, block, plot,
#  category, mass, ppt, avg_ppt, sd_ppt, precip_z).
#
# What changed relative to your original:
#   * NPK effect (LRR) per site (and site x FG) = the `trt` coefficient from a
#     mixed model on plot-level log biomass, with random effects for year
#     (shared weather), block, and plot (repeated measures). Replaces the
#     Hedges et al. variance formula.
#   * Precip slope per site (and site x FG) = the `precip_z` coefficient from
#     a mixed model on CONTROL plots, with the same random effects. Replaces
#     lm() on site-year means, so both classifications use the same framework.
#   * CIs are t-based (Satterthwaite df from lmerTest) for both estimates.
#   * FG 0-filling is built from observed plot-year combinations, and the
#     "present in >= 4 years" filter is applied to the plot-level data used
#     for ALL FG models (including the global mixed models).
##############################################################################

##### SETTINGS #####
pc      <- 0.01   # pseudocount added before log-transforming FG biomass
# (choose relative to your units; try 0.1 and 1 as sensitivity checks)
fg_list <- c("GRAMINOID", "FORB", "WOODY", "LEGUME")


##### FUNCTIONS #####

# Random-effect structures, tried in order until one gives a non-singular fit.
# (1 | year) is always kept: it absorbs shared weather, and for the precip
# slope it is what makes the slope's SE reflect year-level replication.
rand_ladder <- c(
  full      = "(1 | year) + (1 | block) + (1 | plot)",
  no_block  = "(1 | year) + (1 | plot)",
  no_plot   = "(1 | year) + (1 | block)",
  year_only = "(1 | year)"
)

# Fit one site (or site x FG) model and return the coefficient of interest.
#   fixed = fixed-effect predictor ("trt" or "precip_z")
#   term  = coefficient name to extract ("trtNPK" or "precip_z")
# If every structure in the ladder is singular, the full model is used and
# `singular` is flagged TRUE. `model` records which structure was used.
fit_site_lmm <- function(data, fixed, term, response = "log_mass") {
  data <- data[!is.na(data[[response]]) & !is.na(data[[fixed]]), ]
  data$year  <- factor(data$year)
  data$block <- factor(data$block)
  data$plot  <- factor(data$plot)
  
  out_na <- tibble(
    estimate = NA_real_, se = NA_real_, df = NA_real_,
    ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_,
    n_obs = nrow(data), n_plots = n_distinct(data$plot),
    n_years = n_distinct(data$year),
    singular = NA, model = NA_character_
  )
  
  chosen <- NULL; chosen_nm <- NA_character_
  first_fit <- NULL; first_nm <- NA_character_
  
  for (nm in names(rand_ladder)) {
    f <- as.formula(paste(response, "~", fixed, "+", rand_ladder[[nm]]))
    m <- tryCatch(suppressMessages(lmer(f, data = data)),
                  error = function(e) NULL)
    if (is.null(m)) next
    if (is.null(first_fit)) { first_fit <- m; first_nm <- nm }
    if (!isSingular(m)) { chosen <- m; chosen_nm <- nm; break }
  }
  if (is.null(chosen)) { chosen <- first_fit; chosen_nm <- first_nm }
  if (is.null(chosen)) return(out_na)
  
  cf <- summary(chosen)$coefficients
  if (!term %in% rownames(cf)) return(out_na)
  
  est  <- cf[term, "Estimate"]
  se   <- cf[term, "Std. Error"]
  dfs  <- cf[term, "df"]
  crit <- qt(0.975, dfs)
  
  tibble(
    estimate = est, se = se, df = dfs,
    ci_low = est - crit * se, ci_high = est + crit * se,
    p_value = cf[term, "Pr(>|t|)"],
    n_obs = nrow(data), n_plots = n_distinct(data$plot),
    n_years = n_distinct(data$year),
    singular = isSingular(chosen), model = chosen_nm
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


##### TOTAL BIOMASS ANALYSIS #####

## Plot-level total biomass
total_biomass_precip <- biomass_precip %>%
  group_by(site_code, site_name, year, year_trt, trt, block, plot,
           ppt, avg_ppt, sd_ppt, precip_z) %>%
  summarise(total_biomass = sum(mass), .groups = "drop") %>%
  filter(total_biomass > 0) %>%   # log(0) = -Inf; check that this drops ~nothing
  mutate(log_mass = log(total_biomass),
         trt = factor(trt, levels = c("Control", "NPK")))

## Precip slope per site (control plots only)
slopes_total_biomass <- total_biomass_precip %>%
  filter(trt == "Control") %>%
  group_by(site_code) %>%
  group_modify(~ fit_site_lmm(.x, fixed = "precip_z", term = "precip_z")) %>%
  ungroup() %>%
  rename_with(~ paste0("control_", .x), -site_code) %>%
  rename(control_slope = control_estimate, control_p = control_p_value)

## NPK effect (log response ratio) per site: trt coefficient
lrr_total_biomass <- total_biomass_precip %>%
  group_by(site_code) %>%
  group_modify(~ fit_site_lmm(.x, fixed = "trt", term = "trtNPK")) %>%
  ungroup() %>%
  transmute(site_code,
            lrr_mass = estimate, lrr_mass_se = se,
            lrr_mass_ci_low = ci_low, lrr_mass_ci_high = ci_high,
            lrr_df = df, lrr_singular = singular, lrr_model = model)

# Diagnostics: which random-effect structures were used, any failures?
table(slopes_total_biomass$control_model, useNA = "ifany")
table(lrr_total_biomass$lrr_model, useNA = "ifany")
table(slopes_total_biomass$control_singular, useNA = "ifany")
table(lrr_total_biomass$lrr_singular, useNA = "ifany")


## TOTAL BIOMASS SLOPES AND LRRS
## join, classify and tally, graph
slopes_lrr_total_biomass <- left_join(slopes_total_biomass, lrr_total_biomass,
                                      by = "site_code") %>%
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
  labs(x = "Response to NPK (mean log response ratio)",
       y = "Response to precipitation\n(slope of log biomass vs. precip z-score)",
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


## TOTAL BIOMASS STATISTICAL MODELS (unchanged)
model_total_biomass <- lm(control_slope ~ lrr_mass, data = slopes_lrr_total_biomass)
summary(model_total_biomass)

mixed_model_total_biomass <- lmer(log_mass ~ trt * precip_z +
                                    (1 | site_code/plot) + (1 | site_code:year),
                                  data = total_biomass_precip)
summary(mixed_model_total_biomass)


##### BIOMASS BY FUNCTIONAL GROUP ANALYSIS #####

## Observed FG biomass
fg_biomass <- biomass_precip %>%
  filter(category %in% fg_list)

# Sanity check: one row per site-plot-year-FG (stops if not)
stopifnot(!anyDuplicated(fg_biomass[c("site_code", "plot", "year", "category")]))

## 0-fill FGs that were observed at a site but absent from particular plot-years.
## Design is built from the plot-years that were actually sampled (all
## categories), so we only create zeros for plot-years that really exist.
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
  semi_join(fg_keep, by = c("site_code", "category")) %>%
  mutate(log_mass = log(mass + pc),
         trt = factor(trt, levels = c("Control", "NPK")))

## Precip slope per site x FG (control plots only)
slopes_fg_biomass <- fg_biomass_complete_precip %>%
  filter(trt == "Control") %>%
  group_by(site_code, category) %>%
  group_modify(~ fit_site_lmm(.x, fixed = "precip_z", term = "precip_z")) %>%
  ungroup() %>%
  rename_with(~ paste0("control_", .x), -c(site_code, category)) %>%
  rename(control_slope = control_estimate, control_p = control_p_value)

## NPK effect (log response ratio) per site x FG
lrr_fg_biomass <- fg_biomass_complete_precip %>%
  group_by(site_code, category) %>%
  group_modify(~ fit_site_lmm(.x, fixed = "trt", term = "trtNPK")) %>%
  ungroup() %>%
  transmute(site_code, category,
            lrr_mass = estimate, lrr_mass_se = se,
            lrr_mass_ci_low = ci_low, lrr_mass_ci_high = ci_high,
            lrr_df = df, lrr_singular = singular, lrr_model = model)

# Diagnostics
table(slopes_fg_biomass$control_model, useNA = "ifany")
table(lrr_fg_biomass$lrr_model, useNA = "ifany")
table(slopes_fg_biomass$control_singular, useNA = "ifany")
table(lrr_fg_biomass$lrr_singular, useNA = "ifany")


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
  labs(x = "Response to NPK (mean log response ratio)",
       y = "Response to precipitation\n(slope of log biomass vs. precip z-score)",
       title = "FG biomass responses to precipitation and NPK",
       color = "Response category") +
  theme_bw(base_size = 14)

plot_fg_biomass

bar_fg_biomass <- ggplot(response_counts_fg_biomass,
                         aes(x = npk_response, y = n_sites, fill = precip_response)) +
  geom_col() +
  facet_wrap(~ category) +
  labs(x = "NPK response", y = "Number of site-FG combinations",
       fill = "Precipitation response") +
  theme_bw(base_size = 14)

bar_fg_biomass


## FG BIOMASS STATISTICAL MODELS
# Cross-site models (unchanged formulas)
model_fg_biomass <- lmer(control_slope ~ lrr_mass + category + (1 | site_code),
                         data = slopes_lrr_fg_biomass)
summary(model_fg_biomass)

model_fg_biomass_x <- lmer(control_slope ~ lrr_mass * category + (1 | site_code),
                           data = slopes_lrr_fg_biomass)
summary(model_fg_biomass_x)

# Plot-level global model (formula unchanged; now fit to the fg_keep-filtered data).
# Consider adding trt/precip_z x category interactions and a
# (1 | site_code:plot:category) term, as discussed.
mixed_model_fg_biomass <- lmer(log_mass ~ trt * precip_z + category +
                                 (1 | site_code/plot) + (1 | site_code:year),
                               data = fg_biomass_complete_precip)
summary(mixed_model_fg_biomass)


