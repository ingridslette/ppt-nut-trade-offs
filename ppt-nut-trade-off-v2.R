
library(tidyverse)
library(lme4)
library(lmerTest)
library(performance)

biomass <- read.csv("/Users/ingridslette/Desktop/NutNet/full-biomass-2025-12-09.csv",
                 na.strings = c("NULL","NA"))

str(biomass)
summary(biomass)

biomass <- biomass |> 
  filter(year_trt > 0, trt %in% c("Control", "NPK"))

unique(biomass$trt)
unique(biomass$year_trt)
unique(biomass$year)

site_year_counts_biomass <- biomass %>%
  group_by(site_code, trt) %>% 
  summarise(year_count = n_distinct(year), .groups = 'drop')

sites_with_6_biomass_years <- site_year_counts_biomass %>%
  filter(year_count >= 6) %>%
  group_by(site_code) %>% 
  filter(n_distinct(trt) == 2) 

biomass6 <- biomass %>%
  filter(site_code %in% sites_with_6_biomass_years$site_code)

unique(biomass6$site_code)



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
    avg_pet = mean(pet, na.rm = TRUE),
    sd_ppt = sd(ppt, na.rm = TRUE)
  ) %>%
  ungroup()

unique(biomass6$site_code)
unique(precip$site_code)

biomass_precip <- inner_join(biomass6, precip, by = c("site_code", "year"))

unique(biomass_precip$site_code)

biomass_precip <- biomass_precip %>%
  group_by(site_code) %>%
  mutate(min_ppt = min(ppt, na.rm = TRUE),
         max_ppt = max(ppt, na.rm = TRUE)) %>%
  ungroup()


# Filter to keep only sites with an observed ppt range that spans at least +- 1 sd of long-term avg
biomass_precip <- biomass_precip %>%
  group_by(site_code) %>%
  filter(min_ppt <= (avg_ppt - sd_ppt), max_ppt >= (avg_ppt + sd_ppt)) %>%
  ungroup()

unique(biomass_precip$site_code)



cover <- read.csv("/Users/ingridslette/Desktop/full-cover-2025-12-09.csv",
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

sites_with_6_cover_years <- site_year_counts_cover %>%
  filter(year_count >= 6) %>%
  group_by(site_code) %>% 
  filter(n_distinct(trt) == 2) 

cover6 <- cover %>%
  filter(site_code %in% sites_with_6_cover_years$site_code)

unique(cover6$site_code) ## more sites with >5 years of cover data than with >5 years of biomass data
unique(biomass6$site_code) ## every site with >5 years of cover data has >5 years of biomass data

cover6 <- cover %>%
  filter(site_code %in% sites_with_6_biomass_years$site_code)

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
  left_join(taxa_by_site, by = "site_code") %>%
  left_join(years_by_site, by = "site_code")

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
  left_join(cover_trait_cols, by = c("site_code", "Taxon"))

cover_precip <- inner_join(cover_complete, precip, by = c("site_code", "year"))


str(biomass_precip)
summary(biomass_precip)
str(cover_precip)
summary(cover_precip)
View(cover_precip)

