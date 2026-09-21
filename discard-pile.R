#library(boot)
#library(MuMIn)
#library(MASS)
#library(broom)
#library(broom.mixed)
#library(emmeans)
#library(ggpubr)
#library(cowplot)
#library(ggeffects)


## DIAGNOSTICS ON BIOMASS DATASET


# Which sites/years are missing ppt?
biomass_precip %>% filter(is.na(ppt)) %>% count(site_code, year)

# Which sites report which biomass categories?
biomass_precip %>%
  group_by(site_code) %>%
  summarise(cats = paste(sort(unique(category)), collapse = ", ")) %>%
  print(n = Inf)

# Count of biomass categories
table(biomass_precip$category, useNA = "ifany")

# How many sites use each category?
bio %>% distinct(site_code, category) %>%
  count(category, name = "n_sites") %>% arrange(desc(n_sites))

# Where do uncommon labels occur?
bio %>% filter(category %in% c("LIVE", "VASCULAR", "FORB + PHLOX DIFFUSA", "CACTUS")) %>%
  distinct(site_code, year, category) %>% arrange(site_code, year)



# Add column to indicate if row is total biomass
bio <- biomass_precip %>%
  filter(!is.na(ppt)) %>%
  mutate(is_total_row = category %in% c("LIVE"))

# Which site-years report totals only, categories only, or both?
sy_type <- bio %>%
  group_by(site_code, year, plot) %>%
  summarise(has_total = any(is_total_row),
            has_cat   = any(!is_total_row), .groups = "drop")
count(sy_type, has_total, has_cat)


# Double-counting check: plot-years with VASCULAR/LIVE AND other vascular categories
bio %>%
  group_by(site_code, year, plot) %>%
  summarise(pooled = any(category %in% c("VASCULAR", "LIVE")),
            split  = any(category %in% c("GRAMINOID", "FORB", "LEGUME", "WOODY")),
            .groups = "drop") %>%
  count(pooled, split)

# Are plots missing categories that other plots at the same site-year have?
# (Missing rows may mean zero biomass, or "not measured".)
bio %>%
  group_by(site_code, year) %>%
  mutate(n_plots = n_distinct(plot)) %>%
  group_by(site_code, year, category) %>%
  summarise(frac_plots = n_distinct(plot) / first(n_plots), .groups = "drop") %>%
  filter(frac_plots < 1) %>% count(category)


vasc_labels <- c("GRAMINOID", "FORB", "FORB + PHLOX DIFFUSA", "LEGUME", "WOODY",
                 "CACTUS", "PTERIDOPHYTE", "VASCULAR", "LIVE")




core <- c("GRAMINOID", "FORB", "LEGUME", "WOODY")

# 1. Check that Phlox labelling doesn't double-count
phlox_sites <- bio %>% filter(category == "FORB + PHLOX DIFFUSA") %>% distinct(site_code)
bio %>% semi_join(phlox_sites, by = "site_code") %>%
  group_by(site_code, year, plot) %>%
  summarise(both = all(c("FORB", "FORB + PHLOX DIFFUSA") %in% category), .groups = "drop") %>%
  count(both)   # want all FALSE

# 2. Recode and sum any duplicates
bio2 <- bio %>%
  mutate(category = recode(category, "FORB + PHLOX DIFFUSA" = "FORB")) %>%
  filter(category %in% c(core, "VASCULAR", "LIVE")) %>%
  group_by(site_code, year, year_trt, trt, block, plot, ppt, avg_ppt, sd_ppt, category) %>%
  summarise(mass = sum(mass), .groups = "drop")

# 3. Site-years with split data only (drop pooled site-years)
pooled_sy <- bio2 %>% filter(category %in% c("VASCULAR", "LIVE")) %>% distinct(site_code, year)
split_sy  <- bio2 %>% filter(category %in% core) %>% distinct(site_code, year) %>%
  anti_join(pooled_sy, by = c("site_code", "year"))

# 4. Which categories each site ever recorded
site_cat <- bio2 %>% filter(category %in% core) %>% distinct(site_code, category)

# 5. Plot skeleton x recorded categories, then fill zeros
plots <- bio2 %>%
  semi_join(split_sy, by = c("site_code", "year")) %>%
  distinct(site_code, year, year_trt, trt, block, plot, ppt, avg_ppt, sd_ppt)

bio_cat <- plots %>%
  inner_join(site_cat, by = "site_code", relationship = "many-to-many") %>%
  left_join(bio2 %>% filter(category %in% core) %>%
              select(site_code, year, plot, category, mass),
            by = c("site_code", "year", "plot", "category")) %>%
  mutate(mass = tidyr::replace_na(mass, 0),
         precip_z = (ppt - avg_ppt) / sd_ppt,
         trt = factor(trt, levels = c("Control", "NPK")))


bio_cat %>%
  group_by(site_code, category) %>%
  summarise(n_years = n_distinct(year[mass > 0]),
            n_years_total = n_distinct(year),
            prop_zero = mean(mass == 0), .groups = "drop") %>%
  filter(prop_zero > 0.8) %>% arrange(desc(prop_zero))



bio %>%
  semi_join(phlox_sites, by = "site_code") %>%
  filter(category %in% c("FORB", "FORB + PHLOX DIFFUSA")) %>%
  group_by(site_code, year, plot) %>% filter(n() > 1) %>%
  arrange(site_code, year, plot, category) %>%
  select(site_code, year, plot, trt, category, mass)


bio %>%
  filter(site_code == "bnch.us", category %in% c("FORB", "FORB + PHLOX DIFFUSA")) %>%
  group_by(year, category) %>%
  summarise(n_plots = n_distinct(plot), mean_mass = mean(mass), .groups = "drop") %>%
  arrange(year, category)

# Site-year-category combos that are absent from ALL plots,
# at sites that recorded the category in other years
site_year_cat <- bio2 %>%
  filter(category %in% core) %>%
  distinct(site_code, year, category)

gaps <- split_sy %>%
  inner_join(site_cat, by = "site_code", relationship = "many-to-many") %>%
  anti_join(site_year_cat, by = c("site_code", "year", "category"))

gaps %>% count(category)     # how common is it?
gaps %>% arrange(site_code, year, category) %>% print(n = 50)



n_sy <- split_sy %>% count(site_code, name = "n_sy_total")

prev <- bio_cat_f %>%
  group_by(site_code, category) %>%
  summarise(n_nz_years = n_distinct(year[mass > 0]), .groups = "drop") %>%
  left_join(n_sy, by = "site_code") %>%
  mutate(prev_nz = n_nz_years / n_sy_total)

gaps_flag <- gaps %>%
  semi_join(keep, by = c("site_code", "category")) %>%   # only combos that survived the filter
  left_join(prev, by = c("site_code", "category")) %>%
  mutate(suspicious = category %in% c("FORB", "GRAMINOID") | prev_nz >= 0.5)

gaps_flag %>% count(category, suspicious)

# Look at them by eye: are they isolated years or runs of years?
gaps_flag %>% filter(suspicious) %>% arrange(site_code, category, year) %>% print(n = Inf)

bio_cat_final <- bio_cat_f %>%
  anti_join(filter(gaps_flag, suspicious), by = c("site_code", "year", "category"))







### MAIN MODEL

main_model <- lmer(log_mass ~ log_ppt * trt + (1 | site_code / block) + (1 | year), 
                   data = mass_ppt, na.action = na.exclude)

summary(main_model)

r2_main_model <- performance::r2(main_model)
r2_main_model

avg_live_mass_by_trt <- mass_ppt %>%
  group_by(trt) %>%
  summarise(
    avg_live_mass = mean(live_mass, na.rm = TRUE),
    se_live_mass = sd(live_mass, na.rm = TRUE) / sqrt(sum(!is.na(live_mass)))
  )

avg_live_mass_by_trt

control_mean <- avg_live_mass_by_trt$avg_live_mass[avg_live_mass_by_trt$trt == "Control"]
npk_mean <- avg_live_mass_by_trt$avg_live_mass[avg_live_mass_by_trt$trt == "NPK"]

percent_increase <- ((npk_mean - control_mean) / control_mean) * 100

percent_increase

# Model assumptions check 
plot(main_model)
resid <- residuals(main_model)
hist(resid, breaks = 30, main = "Histogram of Residuals")
qqnorm(resid)
qqline(resid)
plot(fitted(main_model), resid, main = "Residuals vs Fitted")

summary(mass_ppt$live_mass)
hist(mass_ppt$live_mass)


## Trying a PPT-PET model

main_model_pet <- lmer(log_mass ~ ppt_pet * trt + (1 | site_code / block / plot) + (1 | year), 
                       data = mass_ppt, na.action = na.exclude)

summary(main_model_pet)

plot(main_model_pet)
resid_pet <- residuals(main_model_pet)
hist(resid_pet, breaks = 30, main = "Histogram of Residuals")
qqnorm(resid_pet)
qqline(resid_pet)
plot(fitted(main_model_pet), resid, main = "Residuals vs Fitted")

AIC(main_model, main_model_pet)
BIC(main_model, main_model_pet)

r2_ppt_pet_model <- performance::r2(main_model_pet)
print(r2_ppt_pet_model)


### Back-transforming data for graphing - allows for non-linear curves on linear scale

# Back transform from log-log scale and graph
fit_model_and_predict_allsites <- function(data) {
  model <- lm(log_mass ~ log_ppt, data = data)
  new_data <- data.frame(log_ppt = seq(min(data$log_ppt, na.rm = TRUE),
                                       max(data$log_ppt, na.rm = TRUE),
                                       length.out = 100))
  preds <- predict(model, newdata = new_data, se.fit = TRUE)
  new_data$predicted_log_mass <- preds$fit
  new_data$se_log_mass <- preds$se.fit
  new_data$predicted_mass <- 10^new_data$predicted_log_mass
  new_data$mass_lower <- 10^(new_data$predicted_log_mass - 1.96 * new_data$se_log_mass)
  new_data$mass_upper <- 10^(new_data$predicted_log_mass + 1.96 * new_data$se_log_mass)
  new_data$trt <- unique(data$trt)
  return(new_data)
}

predictions_allsites <- mass_ppt %>%
  group_by(trt) %>%
  group_modify(~ fit_model_and_predict_allsites(.x)) %>%
  ungroup()

ggplot(data = mass_ppt, aes(x = ppt, y = live_mass, color = trt, shape = trt)) +
  geom_point(alpha = 0.7) + 
  geom_line(data = predictions_allsites, aes(x = 10^log_ppt, y = predicted_mass), linewidth = 1) +
  xlab("Growing Season Precipitation (mm)") + ylab("Biomass (g/m²)") +
  labs(color = "Treatment", shape = "Treatment") +
  scale_color_manual(values = c("#0092E0", "#ff924c")) +
  theme_bw(base_size = 14)


fit_model_and_predict <- function(data) {
  model <- lm(log_mass ~ log_ppt, data = data)
  p_value <- summary(model)$coefficients["log_ppt", "Pr(>|t|)"]
  
  new_data <- data.frame(
    log_ppt = seq(min(data$log_ppt, na.rm = TRUE),
                  max(data$log_ppt, na.rm = TRUE),
                  length.out = 100)
  )
  new_data$predicted_log_mass <- predict(model, newdata = new_data)
  new_data$predicted_mass <- 10^new_data$predicted_log_mass
  new_data$p_value <- p_value
  return(new_data)
}

predictions <- mass_ppt %>%
  group_by(site_code, trt) %>%
  group_modify(~ fit_model_and_predict(.x)) %>%
  ungroup()

predictions_sig <- predictions %>%
  filter(!is.na(p_value) & p_value < 0.05)

ggplot(mass_ppt, aes(x = ppt, y = live_mass, color = trt, shape = trt, fill = trt)) +
  geom_point(alpha = 0.7) +
  geom_line(data = predictions_sig,
            aes(x = 10^log_ppt, y = predicted_mass), 
            linewidth = 1) +
  labs(x = "Growing Season Precipitation (mm)", y = "Biomass (g/m²)",
       color = "Treatment", shape = "Treatment", fill = "Treatment") +
  facet_wrap(~ site_code, scales = "free") +
  theme_bw(base_size = 12) +
  scale_color_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                     #, labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_fill_manual(values = c("Control" = "#0092E0", "NPK" = "#ff924c")
                    #, labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  scale_shape_manual(values = c("Control" = 21, "NPK" = 24)
                     #,labels = c("Control" = "Control", "NPK" = "Fertilized")
  ) +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18))


ggplot(mass_ppt, aes(x = ppt, y = live_mass, color = site_code)) +
  geom_line(data = predictions, aes(x = 10^log_ppt, y = predicted_mass), linewidth = 1) +
  geom_line(data = predictions_allsites, aes(x = 10^log_ppt, y = predicted_mass), 
            linewidth = 1, color = "black") +
  labs(x = "Growing Season Precipitation (mm)", y = "Biomass (g/m²)") +
  facet_wrap(~ trt) +
  theme_bw(base_size = 14)




paired_t_test_r2 <- t.test(results$control_r2, results$npk_r2, paired = TRUE)
paired_t_test_r2


# Function to calculate mean ± SE of R2
mean_se <- function(x) {
  m <- mean(x, na.rm = TRUE)
  se <- sd(x, na.rm = TRUE) / sqrt(length(x))
  return(c(mean = m, se = se))
}

mean_se_control <- mean_se(results$control_r2)
mean_se_npk     <- mean_se(results$npk_r2)

mean_se_control
mean_se_npk


# comparing R2s for ppt-pet model

results_ppt_pet <- data.frame(site_code = character(), 
                              control_r2 = numeric(), 
                              npk_r2 = numeric(),
                              r2_diff = numeric(),
                              control_slope = numeric(), 
                              npk_slope = numeric(),
                              slope_diff = numeric(),
                              stringsAsFactors = FALSE)

for (site in site_codes) {
  site_data_control <- subset(mass_ppt, site_code == site & trt == "Control")
  site_data_npk <- subset(mass_ppt, site_code == site & trt == "NPK")
  control_model <- lm(log_mass ~ ppt_pet, data = site_data_control)
  npk_model <- lm(log_mass ~ ppt_pet, data = site_data_npk)
  control_r2 <- summary(control_model)$r.squared
  npk_r2 <- summary(npk_model)$r.squared
  control_slope <- coef(control_model)["ppt_pet"]
  npk_slope <- coef(npk_model)["ppt_pet"]
  results_ppt_pet <- rbind(results_ppt_pet, data.frame(
    site_code = site,
    control_r2 = control_r2,
    npk_r2 = npk_r2,
    r2_diff = npk_r2 - control_r2,
    control_slope = control_slope,
    npk_slope = npk_slope,
    slope_diff = npk_slope - control_slope
  ))
}

paired_t_test_r2_ppt_pet <- t.test(results_ppt_pet$control_r2, results_ppt_pet$npk_r2, paired = TRUE)
paired_t_test_r2_ppt_pet



# Joining LRR and cover data to main dataframe

mass_ppt_edited <- mass_ppt_edited %>% 
  left_join(lrr_df, by = c("site_code", "year"))

mass_ppt_edited <- mass_ppt_edited %>% 
  left_join(cover_by_site_plot_year, by = c("site_code", "plot", "year"))

mass_ppt_edited <- na.omit(mass_ppt_edited)
unique(mass_ppt_edited$site_code)



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


