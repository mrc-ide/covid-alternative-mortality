library(rtweet)

##### COVID-19 deaths #####

# get the tweet database
tmls <- rtweet::get_timelines(c("YSNECCOVID19"), n = 3200)

# get the media links
links <- tmls$ext_media_url

# now we scrape the text from these images
library(tesseract)
eng <- tesseract("eng")

text_extr <- lapply(links, function(x) {

  if(length(x) == 1) {
    if(is.na(x)) {
      text <- NULL
    } else {
      text <- tesseract::ocr(x, engine = eng)
    }
  } else {
    text <- NULL
  }

  return(text)

})

# which tweet texts include Aden
adens <- lapply(text_extr, function(x) {
  if(length(x) == 1) {
    grepl("Aden",x)
  } else {
    FALSE
  }
})

# scrape the cumulative deaths from these
deaths <- lapply(text_extr[which(unlist(adens))], function(text) {

  if(length(text) == 1) {

as.numeric(strsplit(strsplit(text, "\n")[[1]][grep("Aden",strsplit(text, "\n")[[1]])+2], " ")[[1]][4])

  } else {
    return(NULL)
  }
})

# build into a data.frame
dates <- as.Date(tmls$created_at[which(unlist(adens))])
df_covid <- data.frame("date" = dates,
                       "deaths" = unlist(deaths)) %>% filter(date<=as.Date("2020-12-31"))

# and now tidy for where it has gone wrong
df_covid$deaths[df_covid$deaths>max(df_covid$deaths[1])] <- NA
df_covid$deaths[df_covid$deaths==0] <- NA

df2_covid <- fill(df_covid, deaths, .direction = "up")
df2_covid$deaths <- rev(c(tail(df2_covid$deaths,1), diff(rev(df2_covid$deaths))))
saveRDS(df2_covid, "analysis/data/derived/deaths_time_series/aden_covid_deaths.rds")



##### excess deaths #####

library(zoo)

df_excess <- read.csv("analysis/data/raw/aden_excess_deaths.csv") %>%
  mutate(date=as.Date(date,format="%Y-%m-%d"),
         deaths = c(0,diff(deaths)),
         deaths = rollmean(deaths,k=7,na.pad=TRUE,align="center")) %>%
  filter(date>="2020-03-01",!is.na(deaths))

saveRDS(df_excess,"analysis/data/derived/deaths_time_series/aden_excess_deaths.rds")


### seroprevalence scenarios to consider

## all combinations including IFR

ifr <- c(0.2,0.3,0.4,0.5)

igg_scale_vals <- c(109.5,132,155,176.5,199,219.67,239,266.67,287.5,310)
igg_scenario <- data.frame("igg_scenario" = c("100 days","120 days","140 days","160 days","180 days","200 days",
                                              "220 days","240 days","260 days","280 days"),
                           "igg_serorev"= igg_scale_vals)

combos_igg <- merge(expand.grid("antibody_type"="igg",
                                "ifr"=ifr,
                                "igg_serorev"=igg_scale_vals,
                                "igm_serorev"=54,
                                "igm_scenario"="50 days"
),
igg_scenario,
by="igg_serorev") %>% mutate()


igm_scale_vals <- c(4,9,15,20,26,31.5,37,42.5,48,54,59,65,70,76)
igm_scenario <- data.frame("igm_scenario" = c("5 days","10 days","15 days","20 days","25 days","30 days",
                                              "35 days","40 days","45 days","50 days","55 days","60 days","65 days","70 days"),
                           "igm_serorev"= igm_scale_vals)

combos_igm <- merge(expand.grid("antibody_type"="igm",
                                "ifr"=ifr,
                                "igg_serorev"=199,
                                "igm_serorev"=igm_scale_vals,
                                "igg_scenario"="180 days"
),
igm_scenario,
by="igm_serorev")

combos_iggm <- merge(merge(expand.grid("antibody_type"="iggm",
                                       "ifr"=ifr,
                                       "igg_serorev"=igg_scale_vals,
                                       "igm_serorev"=igm_scale_vals),
                           igg_scenario,
                           by="igg_serorev"),
                     igm_scenario,
                     by="igm_serorev")


combos_all <- rbind(combos_igg %>% dplyr::select(antibody_type,ifr,igg_serorev,igg_scenario,igm_serorev,igm_scenario),
                    combos_igm %>% dplyr::select(antibody_type,ifr,igg_serorev,igg_scenario,igm_serorev,igm_scenario),
                    combos_iggm %>% dplyr::select(antibody_type,ifr,igg_serorev,igg_scenario,igm_serorev,igm_scenario))

saveRDS(combos_all,"analysis/data/derived/seroprevalence/Aden/sero_scenarios_with_IFR.RDS")


### now without IFR for use in model fitting function

combos_no_ifr <- combos_all %>% dplyr::select(-ifr) %>% unique()
saveRDS(combos_no_ifr,"analysis/data/derived/seroprevalence/Aden/sero_scenarios_without_IFR.RDS")



