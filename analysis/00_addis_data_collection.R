library(readxl)
library(tidyverse)
library(aweek)
library(haven)

##### COVID-19 deaths #####

# deaths data taken from EPHI weekly reports (https://ephi.gov.et/download/ephi-pheoc-covid-19/)

df_ephi <- data.frame(read_excel("analysis/data/raw/addis_deaths_ephi.xlsx",sheet="EPHI") %>%
                             select(!notes) %>%
                             mutate(date_begin=as.Date(date_begin),
                                    date_end=as.Date(date_end)))

# nationally reported deaths for Ethiopia

df_who_eth <- read.csv("https://covid19.who.int/WHO-COVID-19-global-data.csv") %>%
  rename(date = Date_reported) %>% filter(Country=="Ethiopia") %>% mutate(date=as.Date(date)) %>%
  filter(date<=as.Date("2021-01-07")) %>%
  mutate(epiweek = date2week(date, floor_day = TRUE)) %>% group_by(epiweek) %>%
  mutate(epiweek_number=as.numeric(substring(epiweek,7,8)),
         weekly_new_deaths=sum(New_deaths)) %>% data.frame()

# split into three periods for which we have different data

set.seed(2904)

# period 1: 6 April - 21 June
# 63 deaths in Addis Ababa out of 72 total in Ethiopia as of 21 June

period_1 <- data.frame(df_who_eth %>% filter(date>=as.Date("2020-04-01")&date<=as.Date("2020-06-21")) %>%
                            select(date,epiweek_number,New_deaths) %>%
                            mutate(New_deaths_aa=round(63/72 * New_deaths)))

period_1 %>% select(New_deaths_aa) %>% sum() # this gives 66 due to rounding but there should only be 63

#randomly sample 3 dates and remove 1 death from these time points from 7 April (as 6 April confirmed as in Addis)

period_1_sample <- period_1 %>% filter(New_deaths_aa!=0&date>=as.Date("2020-04-07")) %>% select(date) %>%
  do(sample_n(.,3,replace=FALSE))

# remove one death from each of the sampled dates
period_1$New_deaths_aa[which(period_1$date %in% period_1_sample$date)] <- period_1$New_deaths_aa[which(period_1$date %in% period_1_sample$date)]-1

period_1 %>% select(New_deaths_aa) %>% sum() #63 deaths

# period 2: 22 June - 9 August
# 295 deaths in Addis Ababa by this date
# implies 232 deaths in Addis during this period compared to 318 at national level

period_2 <- data.frame(df_who_eth %>% filter(date>as.Date("2020-06-21")&date<=as.Date("2020-08-09")) %>%
                            select(date,epiweek_number,New_deaths) %>%
                            mutate(New_deaths_aa=round(232/318 * New_deaths)))

period_2 %>% select(New_deaths_aa) %>% sum() # this gives 230 due to rounding but there should be 232

#randomly sample 2 dates and add 1 death from these time points

period_2_sample <- period_2 %>% filter(New_deaths_aa!=0) %>% select(date) %>%
  do(sample_n(.,2,replace=FALSE))


# add one death to each of the sampled dates
period_2$New_deaths_aa[which(period_2$date %in% period_2_sample$date)] <- period_2$New_deaths_aa[which(period_2$date %in% period_2_sample$date)]+1

period_2 %>% select(New_deaths_aa) %>% sum() #232 deaths


# period 3: 10 August - 31 December
# EPHI publishes weekly number of deaths in Addis so we allocate these according to the national trend

df_who_eth_daily <- data.frame(df_who_eth %>% filter(date>as.Date("2020-08-09")) %>% group_by(epiweek_number) %>%
                                 mutate(daily_death_prop_per_epiweek=New_deaths/weekly_new_deaths))


period_3 <- data.frame(merge(df_who_eth_daily,df_ephi %>% select(epiweek_number,addis_deaths),
                             by="epiweek_number",all.x=TRUE) %>% arrange(date) %>%
                         mutate(New_deaths_aa=round(daily_death_prop_per_epiweek*addis_deaths)) %>%
                         group_by(epiweek_number) %>%
                         mutate(aa_deaths_check=sum(New_deaths_aa)) %>%
                         filter(!is.na(addis_deaths)))

period_3 %>% select(New_deaths_aa) %>% sum()
period_3 %>% select(epiweek_number,aa_deaths_check) %>% unique() %>% select(aa_deaths_check) %>% sum()


# bring together to get daily time series for Addis

df_covid <- rbind(period_1 %>% select(date,New_deaths_aa),
                  period_2 %>% select(date,New_deaths_aa),
                  period_3 %>% select(date,New_deaths_aa)) %>% rename(deaths = New_deaths_aa)

saveRDS(df_covid,"analysis/data/derived/deaths_time_series/addis_covid_deaths.rds")

##### excess deaths #####

# Loading in graveyard mortality data
df_graves <- haven::read_dta("analysis/data/raw/addis_grave_data.dta")

# Wrangling dataframe
deaths_by_date <- df_graves %>%
  filter(!is.na(Cemetery_code), !is.na(dateburialgc)) %>%
  mutate(yr = as.numeric(yr)) %>%
  group_by(Cemetery_code) %>%
  complete(dateburialgc = seq.Date(min(dateburialgc), max(dateburialgc), by = "days"),
           fill = list(burial_death = 0)) %>%
  fill(Cemetery_code, .direction = "up") %>%
  mutate(month_year = format(as.Date(dateburialgc, format="%Y-%m-%d"), "%Y-%m")) %>%
  mutate(yr = ifelse(is.na(yr),
                     as.numeric(format(as.Date(dateburialgc, format="%Y-%m-%d"),"%Y")), yr)) %>%
  mutate(month = ifelse(is.na(month),
                        as.numeric(format(as.Date(dateburialgc, format="%Y-%m-%d"),"%m")), month)) %>%
  group_by(dateburialgc, yr) %>%
  summarise(sum = sum(burial_death),
            month = median(month))

# Estimate baseline death negative binomial distributions for each month for 2015-2019 or 2019 only
baseline_all_years <- deaths_by_date %>%
  filter(dateburialgc >= "2015-01-01" & dateburialgc < "2020-01-01") %>%
  group_by(month) %>%
  summarise(nb = list(fitdistrplus::fitdist(sum, "nbinom")))

baseline_2019_only <- deaths_by_date %>%
  filter(dateburialgc >= "2019-01-01" & dateburialgc < "2020-01-01") %>%
  group_by(month) %>%
  summarise(nb = list(fitdistrplus::fitdist(sum, "nbinom")))


# Generate 100 draws from monthly negative binomials for each period definition and calculation of excess mortality
source("R/excess_utils.R")

draws_all_years <- lapply(1:100, create_realisation, deaths_df = deaths_by_date
                          %>% filter(dateburialgc > "2019-12-31"),
                          baseline = baseline_all_years)
draws_all_years <- do.call(rbind, draws_all_years)
draws_all_years_overall <- draws_all_years %>% group_by(dateburialgc) %>%
  summarise(deaths = mean(deaths)) %>%
  mutate(weekly = zoo::rollmean(deaths, 7, na.pad = TRUE))

draws_2019_only  <- lapply(1:100, create_realisation, deaths_df = deaths_by_date %>%
                             filter(dateburialgc > "2019-12-31"),
                           baseline = baseline_2019_only )
draws_2019_only <- do.call(rbind, draws_2019_only)
draws_2019_only_overall <- draws_2019_only %>% group_by(dateburialgc) %>%
  summarise(deaths = mean(deaths)) %>%
  mutate(weekly = zoo::rollmean(deaths, 7, na.pad = TRUE))

df_excess <- rbind(draws_all_years_overall %>% mutate("baseline"="all_years"),
                   draws_2019_only_overall %>% mutate("baseline"="only_2019")) %>% data.frame()

saveRDS(df_excess,"analysis/data/derived/deaths_time_series/addis_excess_deaths.rds")

