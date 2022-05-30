# Plot results of Aden analysis

library(zoo)
library(tidyverse)
library(cowplot)

# (A) Estimated excess vs reported COVID-19 deaths over time

aden_covid <- readRDS("analysis/data/derived/deaths_time_series/aden_covid_deaths.rds")

aden_excess <- readRDS("analysis/data/derived/deaths_time_series/aden_excess_deaths.rds")

aden_deaths <- rbind(aden_covid %>% mutate(data="Reported COVID-19 deaths"),
                     aden_excess %>% mutate(data="Estimated excess deaths"))


aden_plot_a <- ggplot(aden_deaths %>% filter(date>as.Date("2020-03-01"),
                                             date<max(aden_excess$date)),aes(x=date,y=deaths,col=data))+
  geom_line()+
  theme_bw()+
  labs(col="",x="Date",y="Daily deaths",tag="A")+
  theme(legend.position="bottom",strip.background = element_rect(fill="white"),
        axis.text.x = element_text(angle=45,hjust=1),
        panel.grid.minor = element_blank())+
  scale_colour_manual(values=c("cornflowerblue","#0818A8"))+
  scale_x_date(date_labels = "%b-%y",breaks="1.5 months")+
  geom_hline(yintercept=0,linetype="dashed")


# (B) Estimated seroprevalence over time under different seroreversion half-lifes

excess_sero <- readRDS("analysis/data/derived/seroprevalence/Aden/aden_excess_complete_sero.RDS")

sero_scen_join <- excess_sero %>% filter(date >= as.Date("2020-06-15"),date<=as.Date("2020-12-31")) %>%
  mutate(min=min*100,med=med*100,max=max*100)

sero_base <- excess_sero %>% filter(scenario=="160 days",date <= as.Date("2020-06-15")) %>%
  mutate(min=min*100,med=med*100,max=max*100)

aden_plot_b <- ggplot()+
  theme_bw()+
  geom_ribbon(data=sero_base,aes(x=date,ymin=min,ymax=max),fill="lightgrey",alpha=0.8)+
  geom_line(data=sero_base,aes(x=date,y=med))+
  geom_vline(xintercept = as.Date("2020-06-15"),linetype="dotted")+
  geom_ribbon(data=sero_scen_join %>% filter(scenario %in% c("140 days","160 days","180 days","200 days")),
              aes(x=date,ymin=min,ymax=max,group=scenario,fill=scenario),alpha=0.2)+
  geom_line(data=sero_scen_join %>% filter(scenario %in% c("140 days","160 days","180 days","200 days")),
            aes(x=date,med,group=scenario,col=scenario))+
  theme(strip.background = element_rect(fill="white"),
        axis.text.x = element_text(angle=45,hjust=1),
        panel.grid.minor = element_blank(),
        legend.position="bottom")+
  scale_x_date(date_labels = "%b-%y",breaks="1.5 months",limits=c(as.Date("2020-03-01"),as.Date("2020-12-31")))+
  labs(x="Date",y="Seroprevalence (%)",col="IgG seroreversion \nhalf-life",
       fill="IgG seroreversion \nhalf-life",shape="Observed",tag="B")+
  geom_point(aes(x=as.Date("2020-12-06"),y=0.274*100,shape="IgG + IgM"),size=3)+
  geom_segment(aes(x=as.Date("2020-11-28"),xend=as.Date("2020-12-13"),y=0.274*100,yend=0.274*100))+
  geom_segment(aes(x=as.Date("2020-12-06"),xend=as.Date("2020-12-06"),y=0.2*100,yend=0.363*100))+
  scale_colour_viridis_d(option="plasma",end=0.75,begin=0.2)+
  scale_fill_viridis_d(option="plasma",end=0.75,begin=0.2)+
  guides(colour = guide_legend(order=1,nrow=2,title.position = "top"),shape=guide_legend(order=2,title.position="top"),
         fill = guide_legend(order=1,nrow=2))+
  scale_shape_manual(values=15)


# (C) Heatmap of log likelihood of models under different IFRs and seroreversion half-lifes

ll_IFR_seroreversion <- rbind(
  readRDS("analysis/data/derived/seroprevalence/Aden/aden_ll_IFR02.rds"),
  readRDS("analysis/data/derived/seroprevalence/Aden/aden_ll_IFR03.rds"),
  readRDS("analysis/data/derived/seroprevalence/Aden/aden_ll_IFR04.rds"),
  readRDS("analysis/data/derived/seroprevalence/Aden/aden_ll_IFR05.rds")
)

aden_plot_c <- ggplot(ll_IFR_seroreversion)+
  geom_tile(aes(y=IFR,x=half_life,fill=ll))+
  theme_bw()+
  scale_fill_viridis_c(begin=0,end=0.90,breaks=c(-1200,-600,-5))+
  theme(legend.position = "bottom")+
  guides(fill=guide_colorbar(ticks=FALSE))+
  labs(x="IgG seroreversion half-life",y="Infection Fatality Ratio (IFR)",fill="Model log-likelihood",tag="C")+
  scale_x_continuous(expand=c(0,0),breaks=c(100,120,140,160,180,200,220,240,260,280))+
  scale_y_continuous(expand=c(0,0),breaks=c(0.2,0.3,0.4,0.5))

# combine
plot_grid(aden_plot_a,plot_grid(aden_plot_b,aden_plot_c,ncol=2,align="hv"),align="b",nrow=2,rel_heights = c(1,1.5))
ggsave("analysis/figures/aden_figure_main.png",height=6,width=7)


