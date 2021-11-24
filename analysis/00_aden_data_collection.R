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
saveRDS(df2_covid, "analysis/data/derived/aden_covid_deaths.rds")



##### excess deaths #####

library(zoo)

df_excess <- read.csv("analysis/data/raw/aden_excess_deaths.csv") %>%
  mutate(date=as.Date(date,format="%Y-%m-%d"),
         deaths = c(0,diff(deaths)),
         deaths = rollmean(deaths,k=7,na.pad=TRUE,align="center")) %>%
  filter(date>="2020-03-01",!is.na(deaths))

saveRDS(df_excess,"analysis/data/derived/aden_excess_deaths.rds")





