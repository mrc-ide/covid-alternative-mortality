##trends in mortality from graveyard surveillance

df <- haven::read_dta("analysis/data/raw/addis_grave_data.dta")
end_date <- as.Date("2020-12-31") ### keep this as is for now as is the annual trends

df_year <- df %>% group_by(yr) %>% mutate(deaths_yr=sum(burial_death,na.rm=TRUE)) %>%
  dplyr::select(yr,deaths_yr) %>% unique() %>% filter(yr!=2021,!is.na(yr))

ggplot(df_year,aes(x=yr,y=deaths_yr))+
  geom_col(fill="black")+
  theme_bw()+
  labs(x="year",y="number of deaths")+
  scale_y_continuous(n.breaks=6)
#ggsave("scripts/unfinished_scripts/Addis/mortality_trend.png",height=3,width=6)



##test for significant differences

# this was months 5 and 6 due to the blip
base_rate_all_years <- df %>% filter(month==5|month==6) %>% group_by(yr) %>%
  mutate(deaths_month_yr=sum(burial_death,na.rm=TRUE)) %>%
  dplyr::select(yr,deaths_month_yr) %>% unique() %>% data.frame()

base_rate_all_years <- df %>% filter(!is.na(yr)) %>% group_by(yr) %>%
  mutate(deaths_month_yr=sum(burial_death,na.rm=TRUE)) %>%
  dplyr::select(yr,deaths_month_yr) %>% unique() %>% data.frame()

# poisson.test(base_rate_all_years %>% filter(yr %in% c(2015,2016,2017,2018)) %>% select(deaths_month_yr) %>%
#                 summarise(mean(deaths_month_yr)) %>% round() %>% as.numeric(),
#               base_rate_all_years %>% filter(yr==2019) %>% select(deaths_month_yr) %>%
#                 summarise(mean(deaths_month_yr)) %>% round() %>% as.numeric(),
#               alternative="two.sided")


### do chi squared test instead
chisq.test(x=c(base_rate_all_years %>% filter(yr %in% c(2015,2016,2017,2018)) %>% select(deaths_month_yr) %>%
                 summarise(mean(deaths_month_yr)) %>% round() %>% as.numeric(),
               base_rate_all_years %>% filter(yr==2019) %>% select(deaths_month_yr) %>%
                 summarise(mean(deaths_month_yr)) %>% round() %>% as.numeric()))

## manual check
exp <- base_rate_all_years %>% filter(yr %in% c(2015,2016,2017,2018)) %>% select(deaths_month_yr) %>%
  summarise(mean(deaths_month_yr)) %>% round() %>% as.numeric()
obs <- base_rate_all_years %>% filter(yr==2019) %>% select(deaths_month_yr) %>%
  summarise(mean(deaths_month_yr)) %>% round() %>% as.numeric()

test_stat <- (obs-exp)^2/exp
1-pchisq(test_stat,df=1)

#### correlation between covid and excess death time series


addis_covid <- readRDS("analysis/data/derived/deaths_time_series/addis_covid_deaths.RDS") %>%
  filter(date >= as.Date("2020-04-05")&date<as.Date("2020-11-01"))

addis_excess <- readRDS("analysis/data/derived/deaths_time_series/addis_excess_deaths.RDS") %>%
  select(dateburialgc,weekly,baseline) %>% rename(date=dateburialgc,deaths=weekly) %>%
  mutate(deaths=round(deaths),date=as.Date(date,format="%Y-%m-%d")) %>%
  filter(date >= as.Date("2020-04-01")&date<as.Date("2020-11-01"))

addis_deaths <- merge(addis_excess %>% rename(excess_deaths=deaths),
                      addis_covid %>% rename(covid_deaths=deaths),
                      by="date") %>%
  mutate(baseline_label = ifelse(baseline=="only_2019","2019 baseline","2015 - 2019 baseline"))

addis_deaths %>% group_by(baseline_label) %>% summarise(cor(excess_deaths,covid_deaths,method="spearman"))
addis_deaths %>% filter(date>=as.Date("2020-06-01")) %>%
  group_by(baseline_label) %>% summarise(cor(excess_deaths,covid_deaths,method="spearman"))



#
# #all years
# poisson.test(base_rate_all_years %>% filter(yr!=2020) %>% mutate(avg_baseline=mean(deaths_month_yr)) %>% dplyr::select(avg_baseline) %>% unique() %>% as.numeric() %>% round(),
#              base_rate_all_years %>% filter(yr==2020) %>% mutate(avg_baseline=mean(deaths_month_yr)) %>% dplyr::select(avg_baseline) %>% unique() %>% as.numeric() %>% round(),
#              alternative="two.sided")
#
# #2015 - 2018
# poisson.test(base_rate_all_years %>% filter(yr!=2020,yr!=2019) %>% mutate(avg_baseline=mean(deaths_month_yr)) %>% dplyr::select(avg_baseline) %>% unique() %>% as.numeric() %>% round(),
#              base_rate_all_years %>% filter(yr==2020) %>% mutate(avg_baseline=mean(deaths_month_yr)) %>% dplyr::select(avg_baseline) %>% unique() %>% as.numeric() %>% round(),
#              alternative="two.sided")
#
#
# #2019
# poisson.test(base_rate_all_years %>% filter(yr==2019) %>% mutate(avg_baseline=mean(deaths_month_yr)) %>% dplyr::select(avg_baseline) %>% unique() %>% as.numeric() %>% round(),
#              base_rate_all_years %>% filter(yr==2020) %>% mutate(avg_baseline=mean(deaths_month_yr)) %>% dplyr::select(avg_baseline) %>% unique() %>% as.numeric() %>% round(),
#              alternative="two.sided")
#



## addis total excess mortality under the different baselines

addis_excess <- readRDS("analysis/data/derived/deaths_time_series/addis_excess_deaths.RDS")

#%>%
# filter(baseline==baseline_select) %>%  select(dateburialgc,weekly) %>% rename(date=dateburialgc,deaths=weekly) %>%
#   mutate(deaths=round(deaths),date=as.Date(date,format="%Y-%m-%d")) %>%
#   filter(date >= as.Date(date_filter)&date <= as.Date("2021-01-01")) %>%
#   mutate(deaths=ifelse(deaths<0,0,deaths))


addis_excess %>% filter(baseline=="all_years",dateburialgc>=as.Date("2020-04-05"),dateburialgc<=as.Date("2020-10-31"),weekly>=0)  %>% dplyr::select(weekly) %>% sum() %>% round() #1549

addis_excess %>% filter(baseline=="only_2019",dateburialgc>=as.Date("2020-04-05"),dateburialgc<=as.Date("2020-10-31"),weekly>=0)  %>% dplyr::select(weekly) %>% sum() %>% round() #2535

addis_excess %>% filter(baseline=="all_years",dateburialgc>=as.Date("2020-06-05"),dateburialgc<=as.Date("2020-10-31"),weekly>=0)  %>% dplyr::select(weekly) %>% sum() %>% round() # 1217

addis_excess %>% filter(baseline=="only_2019",dateburialgc>=as.Date("2020-06-05"),dateburialgc<=as.Date("2020-10-31"),weekly>=0)  %>% dplyr::select(weekly) %>% sum() %>% round() # 2137


addis_deaths <- readRDS("analysis/data/derived/deaths_time_series/addis_covid_deaths.RDS") %>%
  filter(date >= as.Date("05/04/2020"),date<=as.Date("2020-10-31"))
addis_deaths$deaths %>% sum() ## 1064



