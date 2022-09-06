### compare excess deaths to reported deaths

library(tidyverse)
library(readxl)


who_rep <- read.csv("https://covid19.who.int/WHO-COVID-19-global-data.csv")
colnames(who_rep)[1] <- "date_reported"

who_rep_countries <- who_rep %>% filter(Country=="Ethiopia"|Country=="Yemen"|Country=="Sudan") %>%
  mutate(month = as.numeric(format.Date(date_reported,"%m")),
         year = as.numeric(format.Date(date_reported,"%Y"))) %>%
  group_by(Country,month,year) %>%
  summarise(reported_deaths = sum(New_deaths,na.rm=TRUE))


who_excess <- read_excel("analysis/data/raw/who_excess.xlsx",sheet="Country by year and month",skip=12) %>%
  filter(country=="Ethiopia"|country=="Yemen"|country=="Sudan")


who_df <- merge(who_rep_countries %>% rename(country=Country),
                who_excess %>% dplyr::select(country, year, month, excess.mean,
                                             cumul.excess.mean,cumul.excess.low,cumul.excess.high),all=TRUE) %>%
  rename(excess_deaths = excess.mean) %>%
  mutate(date = as.Date(paste0("01-",month,"-",year),"%d-%m-%Y"))

ggplot(who_df,aes(x=date))+
  geom_line(aes(y=reported_deaths,group=country,col="Reported COVID-19 \ndeaths"),lwd=1)+
  geom_line(aes(y=excess_deaths,group=country,col="Estimated excess \ndeaths"),lwd=1)+
  facet_wrap(~country,ncol=1,scales="free_y")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position="bottom",
        axis.text.x = element_text(angle=45,hjust=1),
        panel.grid.minor = element_blank())+
  scale_x_date(limits=c(as.Date("2020-04-01"),as.Date("2021-01-01")),
               date_labels="%b-%Y",
               date_breaks = "1 month")+
  scale_colour_manual(values=c("skyblue2","royalblue"))+
  labs(x="Date",y="Monthly Deaths",col="")

cowplot::plot_grid(ggplot(who_df %>% filter(country=="Ethiopia"),aes(x=date))+
                     geom_line(aes(y=reported_deaths,group=country,col="Reported COVID-19 \ndeaths"),lwd=1)+
                     geom_line(aes(y=excess_deaths,group=country,col="Estimated excess \ndeaths"),lwd=1)+
                     theme_bw()+
                     theme(legend.position="none",
                           axis.text.x = element_blank(),
                           axis.ticks.x = element_blank(),
                           panel.grid.minor = element_blank())+
                     scale_x_date(limits=c(as.Date("2020-04-01"),as.Date("2021-01-01")),
                                  date_labels="%b-%Y",
                                  date_breaks = "1 month")+
                     scale_colour_manual(values=c("skyblue2","royalblue"))+
                     labs(x="",y="Monthly Deaths",col="",tag="A"),
                   ggplot(who_df %>% filter(country=="Yemen"),aes(x=date))+
                     geom_line(aes(y=reported_deaths,group=country,col="Reported COVID-19 \ndeaths"),lwd=1)+
                     geom_line(aes(y=excess_deaths,group=country,col="Estimated excess \ndeaths"),lwd=1)+
                     theme_bw()+
                     theme(legend.position="none",
                           axis.text.x = element_blank(),
                           axis.ticks.x = element_blank(),
                           panel.grid.minor = element_blank())+
                     scale_x_date(limits=c(as.Date("2020-04-01"),as.Date("2021-01-01")),
                                  date_labels="%b-%Y",
                                  date_breaks = "1 month")+
                     scale_colour_manual(values=c("skyblue2","royalblue"))+
                     labs(x="",y="Monthly Deaths",col="",tag="B"),
                   ggplot(who_df %>% filter(country=="Sudan"),aes(x=date))+
                     geom_line(aes(y=reported_deaths,group=country,col="Reported COVID-19 \ndeaths"),lwd=1)+
                     geom_line(aes(y=excess_deaths,group=country,col="Estimated excess \ndeaths"),lwd=1)+
                     theme_bw()+
                     theme(legend.position="bottom",
                           panel.grid.minor = element_blank(),
                           axis.text.x = element_text(angle=45,hjust=1))+
                     scale_x_date(limits=c(as.Date("2020-04-01"),as.Date("2021-01-01")),
                                  date_labels="%b-%Y",
                                  date_breaks = "1 month")+
                     scale_colour_manual(values=c("skyblue2","royalblue"))+
                     labs(x="",y="Monthly Deaths",col="",tag="C"),
                   nrow=3,rel_heights = c(1,1,1.5),align="v")
ggsave("analysis/figures/who_covid_excess.png",height=5,width=4)

## rough reporting fractions

# who_df %>% filter(country=="Ethiopia",date<as.Date("2020-10-31"),date>=as.Date("2020-04-01")) %>%
#   mutate(excess = ifelse(excess_deaths<0,0,excess_deaths),
#          rf = sum(reported_deaths)/sum(excess)) # 8.7%

who_df %>% filter(country=="Ethiopia",date<as.Date("2020-10-31"),date>=as.Date("2020-04-01")) %>%
  mutate(excess = ifelse(excess_deaths<0,0,excess_deaths)) %>%
  summarise(covid = sum(reported_deaths),excess = sum(excess),
            rf = covid/excess)


# who_excess %>% filter(country=="Ethiopia",year==2020,month>=4,excess.mean>0) %>% data.frame()  %>% summarise(sum(excess.mean)) ## this gets the old total

# who_df %>% filter(country=="Yemen",date<as.Date("2020-10-01"),date>=as.Date("2020-03-01")) %>%
#   mutate(excess = ifelse(excess_deaths<0,0,excess_deaths),
#          rf = sum(reported_deaths)/sum(excess)) # 5.6%


who_df %>% filter(country=="Yemen",date<as.Date("2020-10-01"),date>=as.Date("2020-03-01")) %>%
  mutate(excess = ifelse(excess_deaths<0,0,excess_deaths)) %>%
  summarise(covid = sum(reported_deaths),excess = sum(excess),
            rf = covid/excess)


# who_df %>% filter(country=="Sudan",date<as.Date("2020-10-01"),date>=as.Date("2020-03-01")) %>%
#   mutate(excess = ifelse(excess_deaths<0,0,excess_deaths),
#          rf = sum(reported_deaths)/sum(excess)) # 9.1%


who_df %>%filter(country=="Sudan",date<as.Date("2020-10-01"),date>=as.Date("2020-03-01")) %>%
  mutate(excess = ifelse(excess_deaths<0,0,excess_deaths)) %>%
  summarise(covid = sum(reported_deaths),excess = sum(excess),
            rf = covid/excess)
