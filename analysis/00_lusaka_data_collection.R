
# from https://www.bmj.com/content/372/bmj.n334
BMJdata <- read.csv("analysis/data/raw/lusaka_bmj_data.csv")

# different methods to understand sampling effort
BMJdata <- BMJdata %>%
  mutate(total_burial_registry_comp = CovSamDeaths*(3676/364)/0.8, # 3676 deaths during this period
         total_burial_registry_comp_strict = CovSamDeaths_Strict*(3676/364)/0.8, # 3676 deaths during this period
         Sampling = c(5,5,5,5,3,3,3*1/10+2*9/10,2),  # Sampling was 1 in 5 for June, July, 1 in 3 for Aug, 1 in 2 for Sep.
         Cap = c(1,1,1,1,1,1,1/10,0)) # Daily tests were capped at 5 or 6 in June, July and Aug. No cap in Sep.
BMJdata <- BMJdata %>%
  mutate(total_sample_effort = (CovSamDeaths*Sampling + CovSamDeaths/TotSamDeaths * Cap*(3676-sum(TotSamDeaths*Sampling))/(sum(Cap)))/0.8,
         total_sample_effort_strict = (CovSamDeaths_Strict*Sampling + CovSamDeaths_Strict/TotSamDeaths * Cap*(3676-sum(TotSamDeaths*Sampling))/(sum(Cap)))/0.8)

# Explanation of High Estimation:
# MinDeaths <- sum(BMJdata$TotSamDeaths*BMJdata$Sampling) # Min number of deaths, estimated from sampling:
# UnaccountedDeaths <- round(3676 - MinDeaths) # Number of deaths unaccounted for
# Unacc_Deaths_Attrib <- BMJdata$Cap*UnaccountedDeaths/(sum(BMJdata$Cap)) # Attribute unaccounted deaths into time periods according to days that had caps
# cvd_death_prev <- BMJdata$CovSamDeaths/BMJdata$TotSamDeaths # Covid death prevalence:
# Min_cvd_deaths <- BMJdata$CovSamDeaths*BMJdata$Sampling # Min cvd deaths from sampling
# sum(Min_cvd_deaths + cvd_death_prev * Unacc_Deaths_Attrib) # Add min cvd deaths with estimated

# and now merge with the lusaka official data
df_off_zambia <- read.csv("analysis/data/raw/zambia_official_covid.csv")
df_off_lusaka <- df_off_zambia %>%
  mutate(date = as.Date(Date), deaths = Total_Deaths) %>%
  filter(Province=="Lusaka" & date < "2020-11-01") %>%
  select(date, deaths) %>% na.omit() %>%
  group_by(date) %>%
  summarise(deaths = sum(deaths)) %>%
  complete(date = seq.Date(min(date), max(date), 1)) %>%
  mutate(deaths = replace_na(deaths, 0)) %>%
  filter(date <= "2020-09-27")

saveRDS(df_off_lusaka, "analysis/data/derived/deaths_time_series/lusaka_official_deaths.rds")

df_off_lusaka <- df_off_lusaka %>%
  filter(date >= "2020-06-08") %>%
  mutate(deaths = cumsum(deaths))

BMJdata_to_fit <- df_off_lusaka[seq(1, 106, 7),] %>%
  mutate(total_burial_registry_comp = rep((BMJdata$total_burial_registry_comp)/2, each = 2),
         total_burial_registry_comp_strict = rep((BMJdata$total_burial_registry_comp_strict)/2, each = 2),
         total_sample_effort = rep((BMJdata$total_sample_effort)/2, each = 2),
         total_sample_effort_strict = rep((BMJdata$total_sample_effort_strict)/2, each = 2)) %>%
  mutate(across(total_burial_registry_comp:total_sample_effort_strict, round)) %>%
  mutate(official_covid_deaths = c(0, diff(deaths))) %>%
  select(-deaths)

BMJdata_to_fit <- BMJdata_to_fit %>%
  mutate(total_burial_registry_comp_mid = round(total_burial_registry_comp + (total_burial_registry_comp_strict - total_burial_registry_comp)/2)) %>%
  mutate(total_sample_effort_mid = round(total_sample_effort + (total_sample_effort_strict - total_sample_effort)/2))

saveRDS(BMJdata_to_fit, "analysis/data/derived/deaths_time_series/lusaka_postmortem_deaths.rds")
