library(zoo)

## -----------------------------------------------------------------------------
## Step 1: Incoming Data
## -----------------------------------------------------------------------------
data <- read.csv("analysis/data/raw/khartoum_deaths.csv")
data$date <- as.Date(as.character(data$date), "%d/%m/%Y")
names(data)[2:3] <- c("cumu", "deaths")

data <- data[order(data$date), ]
data$cumu <- cumsum(data$deaths)

# New data added since last report
# http://sho.gov.sd/corona/index.php

# nov
date <- seq.Date(as.Date("2020-11-27"), as.Date("2020-11-30"), 1)
deaths <- c(6,7,6,10)
to_add <- data.frame(date, cumu = 0, deaths)

# dec - 10 assumed over period of nor reporting to get to total
deaths <- c(6, 9,10,4,6,4,1,10,16,4,3,7,14,12,12,13,9,11,7,5,2, rep(1,10))
date <- seq.Date(as.Date("2020-12-01"), as.Date("2020-12-31"), 1)
to_add <- rbind(to_add, data.frame(date, cumu = 0, deaths))

# jan
deaths <- c(3,7,3)
date <- seq.Date(as.Date("2021-01-01"), as.Date("2021-01-03"), 1)
to_add <- rbind(to_add, data.frame(date, cumu = 0, deaths))

deaths <- rep(5,6) + c(0,0,1,0,1,0) # 32 across these days
date <- seq.Date(as.Date("2021-01-04"), as.Date("2021-01-09"), 1)
to_add <- rbind(to_add, data.frame(date, cumu = 0, deaths))

date <- seq.Date(as.Date("2021-01-10"), as.Date("2021-01-31"), 1)
deaths <- c(7,4,3,6,8,3,7,6,0,0,7,0,4,0,0,6,0,2,0,2,0,0)
to_add <- rbind(to_add, data.frame(date, cumu = 0, deaths))

# feb 2021
date <- seq.Date(as.Date("2021-02-01"), as.Date("2021-02-28"), 1)
deaths <- c(0,0,0,0,0,0,0,0,2,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,2,1,0)
to_add <- rbind(to_add, data.frame(date, cumu = 0, deaths))

# mar 2021
date <- seq.Date(as.Date("2021-03-01"), as.Date("2021-03-31"), 1)
deaths <- c(0,0,0,4,4,0,3,2,4,9,5,3,0,0,0,4,0,7,6,0,9,5,2,5,6,0,9,0,7,3,3)
to_add <- rbind(to_add, data.frame(date, cumu = 0, deaths))

# add to the original data
new <- rbind(data, to_add) %>%
  complete(date = seq.Date(min(date), max(date), 1)) %>%
  mutate(deaths = replace_na(deaths, 0)) %>%
  mutate(cumu = cumsum(deaths))


# By Jul 3 2020 273 # prior to the adjustment made on the 8th November
# By Sep 8 2020 314 # prior to the adjustment made on the 8th November
new %>% filter(date %in% as.Date(c("2020-07-03","2020-09-08")))

# We then need to add the 99 that the government announced on 8th November (cf health cluster email)
ps <- new$deaths[new$date < "2020-11-08"]
rs <- zoo::rollapply(ps, width = 7, FUN = mean, align = "center", partial = TRUE)
p1 <- rs/sum(rs) * 99
r1 <- round(p1)

# 16 left to allocate
99 - sum(r1)

# highest 16 remainders used (and no ties in rank stretching over 16)
incr <- which(order(abs(p1 - round(p1)), decreasing = TRUE) <= 16)
r1[incr] <- r1[incr] + 1

# add these to our data
new$deaths[new$date < "2020-11-08"] <- new$deaths[new$date < "2020-11-08"] + r1

# 46 from private labs betwen nov 8 and nov 17
# https://m.facebook.com/story.php?story_fbid=2764654530474174&id=182728912754405#

orig <- new$deaths[new$date %in% seq.Date(as.Date("2020-11-09"), as.Date("2020-11-17"), 1)]
new$deaths[new$date %in% seq.Date(as.Date("2020-11-09"), as.Date("2020-11-17"), 1)] <- orig + round(exp((1:9)*0.282))
new <- new %>% mutate(cumu = cumsum(deaths))

# By Nov 8 2020 414 # after to the adjustment made on the 8th November
# By Nov 17 2020 464 # after to the adjustment made on the 8th November and the private lab deaths
# By Nov 19 2020 477 # after to the adjustment made on the 8th November and the private lab deaths
# https://www.dabangasudan.org/en/all-news/article/khartoum-state-covid-19-cases-on-the-rise
# By Jan 8 2021 778 deaths
c(273, 315,414,464,477,778) - (new %>% filter(date %in% as.Date(c("2020-07-03","2020-11-07","2020-11-08","2020-11-17","2020-11-19","2021-01-12"))) %>% pull(cumu))

# quick check to see if seems sensible
ggplot(new, aes(date, deaths)) + geom_point() +
  scale_x_date(date_breaks = "1 month") +
  geom_line(aes(y = zoo::rollmean(deaths, 7, na.pad = TRUE)))

# and save this dataset for fitting to in next script
saveRDS(new, "analysis/data/derived/deaths_time_series/khartoum_covid_deaths.rds")
