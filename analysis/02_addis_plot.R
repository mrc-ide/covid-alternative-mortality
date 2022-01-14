# Plot results of Addis Ababa analysis

library(tidyverse)
library(cowplot)

# (A) Estimated excess vs reported COVID-19 deaths over time

addis_covid <- readRDS("analysis/data/derived/deaths_time_series/addis_covid_deaths.RDS") %>%
  filter(date >= as.Date("05/04/2020"))

addis_excess <- readRDS("analysis/data/derived/deaths_time_series/addis_excess_deaths.RDS") %>%
  select(dateburialgc,weekly,baseline) %>% rename(date=dateburialgc,deaths=weekly) %>%
  mutate(deaths=round(deaths),date=as.Date(date,format="%Y-%m-%d")) %>%
  filter(date >= as.Date("2020-04-01")&date <= as.Date("2021-01-01"))

addis_deaths <- merge(addis_excess %>% rename(excess_deaths=deaths),
                      addis_covid %>% rename(covid_deaths=deaths),
                      by="date") %>%
  mutate(baseline_label = ifelse(baseline=="only_2019","2019 baseline","2015 - 2019 baseline"))

addis_plot_a <- ggplot(addis_deaths)+
  geom_rect(aes(xmin=as.Date("2020-05-09"),xmax=as.Date("2020-06-09"),ymin=-Inf,ymax=Inf),fill="#E8E8E8")+
  geom_line(aes(x=date,y=excess_deaths,group=baseline_label,col="Estimated excess deaths"))+
  facet_wrap(~baseline_label,nrow=1)+
  geom_line(aes(x=date,y=covid_deaths,col="Reported COVID-19 deaths"))+
  theme_bw()+
  theme(legend.position="bottom",strip.background = element_rect(fill="white"),
        axis.text.x = element_text(angle=45,hjust=1),
        panel.grid.minor = element_blank())+
  scale_colour_manual(values=c("cornflowerblue","#0818A8"))+
  scale_x_date(date_labels = "%b-%y",breaks="1.5 months")+
  geom_hline(yintercept=0,linetype="dashed")+
  labs(x="Date",y="Daily deaths",col="",tag="A")


# (B) Estimated seroprevalence over time from model fit to excess mortality with 2015 - 2019 baseline

excess_all_years_April_sero <- readRDS("analysis/data/derived/seroprevalence/Addis/addis_excess_sero_all_years_April.RDS")

abdella <- data.frame("date"=rep(as.Date("2020-07-31"),2),
                      "med"=c(1.9,3.5),
                      "antibody"=factor(c("IgG","Combined IgG/IgM"),
                                        levels=c("IgG","Combined IgG/IgM")))

addis_plot_b <- ggplot(excess_all_years_April_sero)+
  theme_bw()+
  geom_ribbon(aes(x=date,ymin=min*100,ymax=max*100,group=antibody,fill=antibody),alpha=0.2)+
  geom_line(aes(x=date,y=med*100,group=antibody,col=antibody),lwd=1)+
  geom_segment(aes(x=as.Date("2020-07-31"),xend=as.Date("2020-07-31"),y=0.004*100,yend=0.037*100))+
  geom_segment(aes(x=as.Date("2020-07-22"),xend=as.Date("2020-08-10"),y=0.037*100,yend=0.037*100))+
  geom_segment(aes(x=as.Date("2020-07-22"),xend=as.Date("2020-08-10"),y=0.004*100,yend=0.004*100))+
  geom_segment(aes(x=as.Date("2020-07-31"),xend=as.Date("2020-07-31"),y=0.017*100,yend=0.054*100))+
  geom_segment(aes(x=as.Date("2020-07-22"),xend=as.Date("2020-08-10"),y=0.054*100,yend=0.054*100))+
  geom_segment(aes(x=as.Date("2020-07-22"),xend=as.Date("2020-08-10"),y=0.017*100,yend=0.017*100))+
  geom_point(data=abdella,aes(x=date,y=med,shape=antibody),size=5)+
  labs(fill="Model \npredicted",x="Date",y="Seroprevalence (%)",shape="Observed \n(Abdella et al.)",
       col="Model \npredicted",tag="B")+
  scale_shape_manual(values=c(18,15))+
  guides(shape=guide_legend(nrow=2,byrow=TRUE,title.position = "top"),
         col=guide_legend(nrow=2,byrow=TRUE,title.position = "top"),
         fill=guide_legend(nrow=2,byrow=TRUE,title.position = "top"))+
  scale_x_date(date_labels = "%b-%y",breaks="1.5 months")+
  scale_y_continuous(breaks=c(0,2,4,6,8,10),limits=c(0,10))+
  scale_colour_viridis_d(begin=0.6,end=0.8)+
  scale_fill_viridis_d(begin=0.6,end=0.8)+
  theme(axis.text.x = element_text(angle=45,hjust=1),
        legend.position="bottom",
        panel.grid.minor = element_blank())

# (c) Model-estimated seroprevalence under difference scenarios

model_estimates <- readRDS("analysis/data/derived/seroprevalence/Addis/addis_seroprevalence_model_estimates.RDS")

model_estimates <- rbind(model_estimates,
                         setNames(data.frame("IgG",1.9,0.4,3.7,"Observed (Abdella et al.)"),
                                  names(model_estimates)),
                         setNames(data.frame("Combined IgG/IgM",3.5,1.7,5.4,"Observed (Abdella et al.)"),
                                  names(model_estimates))) #%>%
 # mutate(median=median*100,lower=lower*100,upper=upper*100)

model_estimates$model <- factor(model_estimates$model,
                                levels=c("Observed (Abdella et al.)","Model predicted (COVID-19)",
                                         "Model predicted (2015 - 19 baseline)", "Model predicted (2019 baseline)",
                                         "Model predicted (2015 - 19 baseline without 1st peak)",
                                         "Model predicted (2019 baseline without 1st peak)"))
model_estimates$antibody <- factor(model_estimates$antibody,levels=c("IgG","Combined IgG/IgM"))


addis_plot_c <- ggplot(model_estimates %>% filter(model!="Model predicted (2015 - 19 baseline)"))+
  geom_point(aes(x=model,y=median,shape=antibody,group=antibody,col=model),position=position_dodge(width=0.5),size=5)+
  geom_errorbar(aes(x=model,ymin=lower,ymax=upper,group=antibody,col=model),position=position_dodge(width=0.5))+
  theme_bw()+
  scale_shape_manual(values=c(18,15))+
  scale_colour_viridis_d(option="magma",labels = function(x) str_wrap(x, width = 30),end=0.9)+
  labs(x="Model",y="Seroprevalence (%)",tag="C",col="",shape="")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 16))+
  scale_y_continuous(breaks=c(0,2,4,6,8,10),limits=c(0,10))+
  guides(shape=guide_legend(order=2),col="none")+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle=90,vjust=0.5),
        panel.grid.minor = element_blank())

# combine
plot_grid(addis_plot_a,plot_grid(addis_plot_b,addis_plot_c,ncol=2,align="hv"),nrow=2,rel_heights = c(1,1.5))
ggsave("analysis/figures/addis_figure_main.png",height=8,width=7)







