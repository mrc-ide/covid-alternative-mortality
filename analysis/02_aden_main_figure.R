# Plot results of Aden analysis

library(zoo)
library(tidyverse)
library(cowplot)

# (A) Estimated excess vs reported COVID-19 deaths over time

aden_covid <- readRDS("analysis/data/derived/deaths_time_series/aden_covid_deaths.rds") %>%
  arrange(date) %>% complete(date = seq.Date(as.Date("2020-03-01"),max(date),by="day")) %>% data.frame() %>%
  mutate(deaths = ifelse(is.na(deaths),0,deaths)) ## just extending so plots 0 from beginning of excess time series

aden_excess <- readRDS("analysis/data/derived/deaths_time_series/aden_excess_deaths.rds")

aden_deaths <- rbind(aden_covid %>% mutate(data="Reported COVID-19 deaths"),
                     aden_excess %>% mutate(data="Estimated excess deaths"))


aden_plot_a <- ggplot(aden_deaths %>% filter(date>as.Date("2020-03-01"),
                                             date<max(aden_excess$date)),aes(x=date,y=deaths,col=data))+
  geom_line()+
  theme_bw()+
  labs(col="",x="Date",y="Daily \ndeaths",tag="A")+
  theme(legend.position="top",strip.background = element_rect(fill="white"),
        axis.text.x = element_text(angle=45,hjust=1),
        panel.grid.minor = element_blank())+
  scale_colour_manual(values=c("cornflowerblue","#0818A8"))+
  scale_x_date(date_labels = "%b-%y",breaks="1.5 months")#+
#geom_hline(yintercept=0,linetype="dashed")


# (B) Estimated seroprevalence over time under different seroreversion half-lifes for IgG

excess_sero_iggm <- readRDS("analysis/data/derived/seroprevalence/Aden/aden_sero_excess_time_series.rds") %>%
  filter(ifr==0.3,antibody=="iggm",igm_scenario=="50 days",
         igg_scenario %in% c("140 days","160 days","180 days","200 days"))

sero_scen_join_iggm <- excess_sero_iggm %>% filter(date >= as.Date("2020-06-15"),date<=as.Date("2020-12-31")) %>%
  mutate(min=min*100,med=med*100,max=max*100)

sero_base_iggm <- excess_sero_iggm %>% filter(igg_scenario=="160 days",date <= as.Date("2020-06-15")) %>%
  mutate(min=min*100,med=med*100,max=max*100)

aden_plot_b <- ggplot()+
  theme_bw()+
  geom_rect(aes(xmin=as.Date("2020-11-28"),xmax=as.Date("2020-12-13"),ymin=-Inf,ymax=Inf),fill="#E8E8E8")+
  geom_ribbon(data=sero_base_iggm,aes(x=date,ymin=min,ymax=max),fill="grey",alpha=0.8)+
  geom_line(data=sero_base_iggm,aes(x=date,y=med))+
  geom_vline(xintercept = as.Date("2020-06-15"),linetype="dashed")+
  geom_ribbon(data=sero_scen_join_iggm,
              aes(x=date,ymin=min,ymax=max,group=igg_scenario,fill=igg_scenario),alpha=0.2)+
  geom_line(data=sero_scen_join_iggm,
            aes(x=date,med,group=igg_scenario,col=igg_scenario))+
  theme(strip.background = element_rect(fill="white"),
        axis.text.x = element_text(angle=45,hjust=1),
        panel.grid.minor = element_blank(),
        legend.position="right",
        legend.text = element_text(size=10),
        #legend.title = element_text(size=10)
        )+
  scale_x_date(date_labels = "%b-%y",breaks="1.5 months",limits=c(as.Date("2020-03-01"),as.Date("2020-12-31")))+
  labs(x="Date",y="Combined IgG/IgM \nseroprevalence (%)",col="IgG seroreversion \nhalf-life",
       fill="IgG seroreversion \nhalf-life",shape=" ",tag="B")+
  geom_point(aes(x=as.Date("2020-12-06"),y=0.274*100,shape="Reported \ncombined IgG/IgM"),size=2.5)+
  #geom_segment(aes(x=as.Date("2020-11-28"),xend=as.Date("2020-12-13"),y=0.2*100,yend=0.2*100))+
  #geom_segment(aes(x=as.Date("2020-11-28"),xend=as.Date("2020-12-13"),y=0.363*100,yend=0.363*100))+
  geom_linerange(aes(x=as.Date("2020-12-06"),ymin=0.256*100,ymax=0.293*100))+
  scale_colour_viridis_d(option="plasma",end=0.75,begin=0.2)+
  scale_fill_viridis_d(option="plasma",end=0.75,begin=0.2)+
  guides(colour = guide_legend(order=1,nrow=2,
                               title.position = "top"),
         shape=guide_legend(order=2,title.position="top"),
         fill = guide_legend(order=1,nrow=2
                             ))+
  scale_shape_manual(values=15)


# (C) Estimated seroprevalence over time under different seroreversion half-lifes for IgG

excess_sero_igg <- readRDS("analysis/data/derived/seroprevalence/Aden/aden_sero_excess_time_series.rds") %>%
  filter(ifr==0.3,antibody=="igg",igm_scenario=="50 days",
         igg_scenario %in% c("140 days","160 days","180 days","200 days"))

sero_scen_join_igg <- excess_sero_igg %>% filter(date >= as.Date("2020-06-15"),date<=as.Date("2020-12-31")) %>%
  mutate(min=min*100,med=med*100,max=max*100)

sero_base_igg <- excess_sero_igg %>% filter(igg_scenario=="160 days",date <= as.Date("2020-06-15")) %>%
  mutate(min=min*100,med=med*100,max=max*100)

aden_plot_c <- ggplot()+
  theme_bw()+
  geom_rect(aes(xmin=as.Date("2020-11-28"),xmax=as.Date("2020-12-13"),ymin=-Inf,ymax=Inf),fill="#E8E8E8")+
  geom_ribbon(data=sero_base_igg,aes(x=date,ymin=min,ymax=max),fill="grey",alpha=0.8)+
  geom_line(data=sero_base_igg,aes(x=date,y=med))+
  geom_vline(xintercept = as.Date("2020-06-15"),linetype="dashed")+
  geom_ribbon(data=sero_scen_join_igg,
              aes(x=date,ymin=min,ymax=max,group=igg_scenario,fill=igg_scenario),alpha=0.2)+
  geom_line(data=sero_scen_join_igg,
            aes(x=date,med,group=igg_scenario,col=igg_scenario))+
  theme(strip.background = element_rect(fill="white"),
        axis.text.x = element_text(angle=45,hjust=1),
        panel.grid.minor = element_blank(),
        legend.position="bottom",
        legend.text = element_text(size=10))+
  scale_x_date(date_labels = "%b-%y",breaks="1.5 months",limits=c(as.Date("2020-03-01"),as.Date("2020-12-31")))+
  labs(x="Date",y="IgG seroprevalence (%)",col="IgG seroreversion \nhalf-life",
       fill="IgG seroreversion \nhalf-life",shape=" ",tag="C")+
  geom_point(aes(x=as.Date("2020-12-06"),y=0.25*100,shape="Reported IgG"),size=2.5)+
  #geom_segment(aes(x=as.Date("2020-11-28"),xend=as.Date("2020-12-13"),y=0.2*100,yend=0.2*100))+
  #geom_segment(aes(x=as.Date("2020-11-28"),xend=as.Date("2020-12-13"),y=0.363*100,yend=0.363*100))+
  geom_linerange(aes(x=as.Date("2020-12-06"),ymin=0.232*100,ymax=0.269*100))+
  scale_colour_viridis_d(option="plasma",end=0.75,begin=0.2)+
  scale_fill_viridis_d(option="plasma",end=0.75,begin=0.2)+
  guides(#colour = guide_legend(order=1,nrow=2,title.position = "top"),
    shape=guide_legend(order=2,title.position="top"),
    colour = "none",
    fill = "none"
         #fill = guide_legend(order=1,nrow=2)
    )+
  scale_shape_manual(values=18)


# (D) Heatmap of log likelihood of models under different IFRs and seroreversion half-lifes

## igm constant at 45 days and then can sum together
ll_IFR_seroreversion <- readRDS("analysis/data/derived/seroprevalence/Aden/aden_sero_excess_ll.rds") %>%
  filter(antibody!="igm",
         igm_scenario=="50 days") %>% group_by(igg_scenario,ifr) %>%
  summarise(ll = sum(ll)) %>%
  mutate(igg_halflife = as.numeric(substring(igg_scenario,1,3)))



aden_plot_d <- ggplot(ll_IFR_seroreversion )+
  geom_tile(aes(y=ifr,x=igg_halflife,fill=-ll))+
  theme_bw()+
  scale_fill_viridis_c(begin=0.9,end=0,trans="log",breaks=c(5,25,200,1200))+
  theme(legend.position = "right")+
  guides(fill=guide_colorbar(ticks=FALSE))+
  labs(x="IgG seroreversion half-life (days)",y="Infection Fatality \nRatio (IFR) (%)",fill="Negative \nmodel log-likelihood",tag="D")+
  scale_x_continuous(expand=c(0,0),breaks=c(100,120,140,160,180,200,220,240,260,280))+
  scale_y_continuous(expand=c(0,0),breaks=c(0.2,0.3,0.4,0.5))


plot_grid(aden_plot_a,
          plot_grid(aden_plot_b+theme(legend.position = "none"),aden_plot_c+theme(legend.position = "none"),
                    plot_grid(get_legend(aden_plot_b),
                              get_legend(aden_plot_c),
                              ncol=1,rel_heights = c(4,1)),
                    nrow=1,rel_widths = c(1,1,0.6)),
          plot_grid(aden_plot_d+theme(legend.position = "none"),get_legend(aden_plot_d),rel_widths = c(2,0.6)),
          nrow=3,align="h",axis="l",
          rel_heights=c(1.2,1.3,1))
ggsave("analysis/figures/aden_figure_main_updated.png",height=7,width=8)
### then used paint to fix the legend & align top one

