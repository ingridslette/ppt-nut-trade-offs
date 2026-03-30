
### Main model

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


