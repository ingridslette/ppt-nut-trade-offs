library(tidyverse)

local <- read.csv("/Users/ingridslette/Desktop/Weather_monthly_20220615.csv")

unique(local$site_code)

#local <- local %>% 
#  filter(site_code %in% c("bayr.de","ukul.za"))

local <- local %>% 
  dplyr::select(site_code, month, year, pre)


ghcn <- read.csv("/Users/ingridslette/Desktop/NutNet/nutnet-ghcn-ppt-data.csv")

unique(ghcn$STATION)

ghcn <- ghcn %>% 
  mutate(site_code = case_when(
    STATION == "USS0020K04S" ~ "sage.us",
    STATION == "ASN00010044" ~ "mtca.au",
    STATION == "USC00350652" ~ "bnch.us",
    STATION == "USC00130203" ~ "cbgb.us",
    STATION == "USC00252020" ~ "doane.us",
    STATION == "USC00212881" ~ "cdcr.us",
    STATION == "USC00353692" ~ "hart.us",
    STATION == "GME00125074" ~ "badlau.de",
    STATION == "ASN00083084" ~ "bogong.au",
    STATION == "USC00157049" ~ "hall.us",
    STATION == "USW00003131" ~ "elliot.us",
    STATION == "USC00080236" ~ "arch.us",
    STATION == "USW00023275" ~ "hopl.us",
    STATION == "GM000004204" ~ "jena.de",
    STATION == "USC00153194" ~ "spin.us",
    STATION == "USC00145628" ~ "saline.us",
    STATION == "ASN00082165" ~ "nilla.au",
    STATION == "USW00094074" ~ "sgs.us",
    STATION == "USC00203504" ~ "kbs.us",
    STATION == "ASN00080002" ~ "kiny.au",
    STATION == "ASN00014847" ~ "kidman.au",
    STATION == "USC00418646" ~ "temple.us",
    STATION == "USC00247894" ~ "msla.us",
    STATION == "ASN00010626" ~ "ping.au",
    STATION == "USC00111743" ~ "trel.us",
    STATION == "ASN00067105" ~ "yarra.au",
    STATION == "NLE00152497" ~ "veluwe.nl",
    STATION == "CAW00064757" ~ "koffler.ca"
  ))

unique(ghcn$site_code)

#ghcn <- ghcn %>% 
#  filter(site_code %in% c("sage.us","mtca.au","bnch.us","elliot.us"))

ghcn_monthly <- ghcn %>% 
  group_by(site_code, month, year) %>% 
  summarise(PRCP_monthly = sum(PRCP))


cru <- read.csv('/Users/ingridslette/Desktop/CRU/CRU-monthly-pre-pet-1901-2024.csv')

unique(cru$site_code)

cru <- cru %>%
  rename(month_abb = month)

cru$month <- match(cru$month_abb, month.abb)

cru_pre <- cru %>% 
  dplyr::select(site_code, month, year, pre_mm.month)

cru_pet <- cru %>% 
  dplyr::select(site_code, month, year, pet_mm.month)

#cru_ppt <- cru %>% 
#  filter(site_code == "bldr.us")

mswep <- read.csv("/Users/ingridslette/Library/CloudStorage/GoogleDrive-slett152@umn.edu/Shared drives/NutNet_DRAGNet_Shared_External/NutNet Shared/NutNet Non-Core Data/weather/MSWEP/precip-daily-mswep-2025-02-19.csv")

unique(mswep$site_code)

#mswep <- mswep %>% 
#  filter(site_code %in% c("doane.us","marc.ar","sedg.us","bnbt.us","spin.us","yarra.au","cdcr.us",
                          # "cbgb.us","sier.us","ping.au","koffler.ca","kbs.us","smith.us",
                          # "look.us","trel.us","frue.ch","saana.fi","cowi.ca","nilla.au","arch.us",
                          # "konz.us","hall.us","chilcas.ar","cedr.us","mcla.us","valm.ch",
                          # "badlau.de","jena.de","pape.de","comp.pt","kiny.au","bogong.au",
                          # "temple.us","cdpt.us","ethass.au","hopl.us","kidman.au",
                          # "lancaster.uk","lagoas.br","ethamc.au","burrawan.au",
                          # "lubb.us","msla.us","kilp.fi","sgs.us","veluwe.nl", "saline.us",
                          # "sevi.us","potrok.ar","hero.uk","msla_2.us","msla_3.us"))

mswep_monthly <- mswep %>% 
  group_by(site_code, month, year) %>% 
  summarise(precip_monthly = sum(precip))

unique(local$site_code)
unique(cru_pre$site_code)
unique(ghcn_monthly$site_code)
unique(mswep_monthly$site_code)

mswep_monthly <- mswep_monthly %>%
  rename(ppt = precip_monthly) %>%
  mutate(source = "mswep")

ghcn_monthly <- ghcn_monthly %>%
  rename(ppt = PRCP_monthly) %>%
  mutate(source = "ghcn")

local <- local %>%
  rename(ppt = pre) %>%
  mutate(source = "local")

cru_pre <- cru_pre %>%
  rename(ppt = pre_mm.month) %>%
  mutate(source = "cru")

ppt_data <- bind_rows(mswep_monthly, ghcn_monthly, local, cru_pre)

unique(ppt_data$site_code)
unique(ppt_data$source)


# Join in PET data from CRU for all sites
ppt_pet_data <- left_join(ppt_data, cru_pet, by = c("site_code", "month", "year"))
unique(ppt_pet_data$site_code)

# for sites where the growing season spans multiple calendar years:
ppt_pet_data$gs_year <- ppt_pet_data$year 

ppt_pet_data$gs_year <- as.numeric(ppt_pet_data$gs_year)
ppt_pet_data$year <- as.numeric(ppt_pet_data$year)
ppt_pet_data$month <- as.numeric(ppt_pet_data$month)

unique(ppt_pet_data$site_code)

ppt_pet_data <- ppt_pet_data %>%
  mutate(gs_year = case_when(
    site_code == "bogong.au" & month %in% c(10, 11, 12) ~ year + 1,
    site_code == "burrawan.au" & month %in% c(10, 11, 12) ~ year + 1,
    site_code == "chilcas.ar" & month %in% c(08, 09, 10, 11, 12) ~ year + 1,
    site_code == "comp.pt" & month %in% c(10, 11, 12) ~ year + 1,
    site_code == "elliot.us" & month %in% c(11, 12) ~ year + 1,
    site_code == "ethamc.au" & month %in% c(5, 6, 7, 8, 9, 10, 11, 12) ~ year + 1,
    site_code == "ethass.au" & month %in% c(5, 6, 7, 8, 9, 10, 11, 12) ~ year + 1,
    site_code == "gilb.za" & month %in% c(9, 10, 11, 12) ~ year + 1,
    site_code == "hart.us" & month %in% c(10, 11, 12) ~ year + 1,
    site_code == "hopl.us" & month %in% c(11, 12) ~ year + 1,
    site_code == "kidman.au" & month %in% c(11, 12) ~ year + 1,
    site_code == "lagoas.br" & month %in% c(7, 8, 9, 10, 11, 12) ~ year + 1,
    site_code == "mcla.us" & month %in% c(11, 12) ~ year + 1,
    site_code == "nilla.au" & month %in% c(2, 3 ,4, 5, 6, 7, 8, 9, 10, 11, 12) ~ year + 1,
    site_code == "sedg.us" & month %in% c(11, 12) ~ year + 1,
    site_code == "sier.us" & month %in% c(11, 12) ~ year + 1,
    site_code == "smith.us" & month %in% c(10, 11, 12) ~ year + 1,
    site_code == "spin.us" & month %in% c(6, 7, 8, 9, 10, 11, 12) ~ year + 1,
    site_code == "ukul.za" & month %in% c(9, 10, 11, 12) ~ year + 1,
    site_code == "yarra.au" & month %in% c(9, 10, 11, 12) ~ year + 1,
    TRUE ~ year
  ))

# keep only growing season months at each site (there is definitely a better way to do this...)
ppt_pet_gs_only <- filter(ppt_pet_data, site_code =="arch.us" & month %in% c(5, 6, 7, 8, 9, 10) |
                            site_code =="badlau.de" & month %in% c(4, 5, 6, 7, 8, 9, 10) |
                            site_code =="bayr.de" & month %in% c(3, 4, 5, 6, 7, 8, 9) |
                            site_code =="bldr.us" & month %in% c(3, 4, 5, 6, 7, 8) |
                            site_code =="bnbt.us" & month %in% c(5, 6, 7, 8, 9) |
                            site_code =="bnch.us" & month %in% c(4, 5, 6, 7, 8) |
                            site_code =="bogong.au" & month %in% c(10, 11, 12, 1) |
                            site_code =="burrawan.au" & month %in% c(10, 11, 12, 1, 2, 3, 4, 5) |
                            site_code =="cbgb.us" & month %in% c(5, 6, 7, 8, 9, 10) |
                            site_code =="cdcr.us" & month %in% c(4, 5, 6, 7, 8) |
                            site_code =="cedr.us" & month %in% c(4, 5, 6, 7, 8) |
                            site_code =="cdpt.us" & month %in% c(4, 5, 6, 7) |
                            site_code =="chilcas.ar" & month %in% c(8, 9, 10, 11, 12, 1, 2, 3) |
                            site_code =="comp.pt" & month %in% c(10, 11, 12, 1, 2, 3, 4, 5) |
                            site_code =="cowi.ca" & month %in% c(4, 5, 6, 7) |
                            site_code =="doane.us" & month %in% c(5, 6, 7, 8, 9, 10, 11) |
                            site_code =="elliot.us" & month %in% c(11, 12, 1, 2, 3, 4) |
                            site_code =="ethamc.au" & month %in% c(5, 6, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4) |
                            site_code =="ethass.au" & month %in% c(5, 6, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4) |
                            site_code =="frue.ch" & month %in% c(4, 5, 6, 7, 8, 9) |
                            site_code =="hall.us" & month %in% c(4, 5, 6, 7, 8, 9) |
                            site_code =="hero.uk" & month %in% c(4, 5, 6, 7, 8, 9, 10) |
                            site_code =="hopl.us" & month %in% c(11, 12, 1, 2, 3, 4) |
                            site_code =="jena.de" & month %in% c(3, 4, 5, 6, 7, 8, 9, 10) |
                            site_code =="kbs.us" & month %in% c(4, 5, 6, 7, 8, 9) |
                            site_code =="kidman.au" & month %in% c(11, 12, 1, 2, 3, 4) |
                            site_code =="kilp.fi" & month %in% c(6, 7, 8) |
                            site_code =="kiny.au" & month %in% c(5, 6, 7, 8, 9, 10) |
                            site_code =="koffler.ca" & month %in% c(4, 5, 6, 7, 8) |
                            site_code =="konz.us" & month %in% c(5, 6, 7, 8, 9) |
                            site_code =="lagoas.br" & month %in% c(7, 8, 9, 10, 11, 12, 1, 2) |
                            site_code =="lancaster.uk" & month %in% c(3, 4, 5, 6, 7, 8) |
                            site_code =="look.us" & month %in% c(3, 4, 5, 6, 7, 8) | 
                            site_code =="lubb.us" & month %in% c(3, 4, 5, 6, 7, 8, 9, 10) |
                            site_code =="marc.ar" & month %in% c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12) |
                            site_code =="mcla.us" & month %in% c(11, 12, 1, 2, 3, 4) |
                            site_code =="msla_3.us" & month %in% c(4, 5, 6, 7) |
                            site_code =="msla.us" & month %in% c(4, 5, 6, 7) |
                            site_code =="msla_2.us" & month %in% c(4, 5, 6, 7) |
                            site_code =="mtca.au" & month %in% c(8, 9, 10) |
                            site_code =="nilla.au" & month %in% c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12) |
                            site_code =="pape.de" & month %in% c(6, 7, 8, 9) |
                            site_code =="ping.au" & month %in% c(4, 5, 6, 7, 8, 9, 10) |
                            site_code =="potrok.ar" & month %in% c(1, 2, 3, 4) |
                            site_code =="saana.fi" & month %in% c(6, 7, 8) |
                            site_code =="sage.us" & month %in% c(4, 5, 6, 7) |
                            site_code =="saline.us" & month %in% c(5, 6, 7, 8, 9) |
                            site_code =="sedg.us" & month %in% c(11, 12, 1, 2, 3, 4, 5, 6, 7) |
                            site_code =="sevi.us" & month %in% c(4, 5, 6, 7, 8, 9, 10, 11) |
                            site_code =="sgs.us" & month %in% c(4, 5, 6, 7, 8) |
                            site_code =="sier.us" & month %in% c(11, 12, 1, 2, 3, 4) |
                            site_code =="smith.us" & month %in% c(10, 11, 12, 1, 2, 3, 4, 5, 6) |
                            site_code =="spin.us" & month %in% c(1, 2, 3, 4, 5) |
                            site_code =="temple.us" & month %in% c(3, 4, 5, 6, 7, 8, 9, 10) |
                            site_code =="trel.us" & month %in% c(1, 2, 3, 4, 5, 6, 7, 8, 9) |
                            site_code =="ukul.za" & month %in% c(9, 10, 11, 12, 1, 2, 3, 4) |
                            site_code =="valm.ch" & month %in% c(6, 7, 8) |
                            site_code =="veluwe.nl" & month %in% c(3, 4, 5, 6, 7, 8) |
                            site_code =="yarra.au" & month %in% c(9, 10, 11, 12, 1, 2, 3)
                          
)

unique(ppt_pet_gs_only$site_code)

# sum to get total growing season precip per year per site

ppt_pet_annual_gs_only <- ppt_pet_gs_only %>% 
  group_by(site_code, gs_year) %>% 
  summarise(ags_ppt = sum(ppt), ags_pet = sum(pet_mm.month))

ppt_pet_annual_gs_only <- ppt_pet_annual_gs_only %>%
  rename(year = gs_year, ppt = ags_ppt, pet = ags_pet)

unique(ppt_pet_annual_gs_only$site_code)

ppt_pet_annual_gs_only <- ppt_pet_annual_gs_only %>% 
  mutate(ppt_pet = ppt-pet)


write.csv(ppt_pet_annual_gs_only, file = "/Users/ingridslette/Desktop/NutNet/ppt_pet_annual_gs_only_2025-10-09.csv", row.names=FALSE)




