# plots of model fits and estimated seroprevalence for supplementary material for Aden

library(tidyverse)
library(cowplot)

source("R/mcmc_utils.R")

# covid fit

aden_covid_fit_raw <- readRDS("analysis/data/derived/model_fits/Aden/aden_covid_fit.RDS")
aden_covid_fit <- generate_draws(aden_covid_fit_raw,generate_parameters(aden_covid_fit_raw,draws=100),parallel=TRUE)
aden_covid_sero <- readRDS("analysis/data/derived/seroprevalence/Aden/aden_sero_covid_time_series.RDS")

aden_covid_deaths_plot <- plot(aden_covid_fit,var_select="deaths",particle_fit=TRUE)+
  scale_x_date(date_breaks="1.5 months",date_labels = "%b-%y",limits=c(as.Date("2020-04-01"),as.Date("2020-12-31")))+
  theme(axis.text.x = element_text(angle=45,hjust=1),legend.position="none",
        panel.grid.minor = element_blank())+labs(tag="A")

aden_covid_sero_plot_igg <- ggplot(aden_covid_sero %>% filter(antibody=="igg",
                                                              igm_scenario=="50 days",
                                                              igg_scenario=="100 days"|igg_scenario=="280 days"),
                               aes(x=date,y=med*100,ymin=min*100,ymax=max*100))+
  geom_rect(aes(xmin=as.Date("2020-11-28"),xmax=as.Date("2020-12-13"),ymin=-Inf,ymax=Inf),fill="#E8E8E8")+
  geom_line(aes(col=igg_scenario))+
  geom_ribbon(alpha=0.2,aes(fill=igg_scenario))+
  theme_bw()+
  geom_point(aes(x=as.Date("2020-12-06"),y=0.25*100,shape="Reported IgG"),size=3)+
  geom_segment(aes(x=as.Date("2020-12-06"),xend=as.Date("2020-12-06"),y=0.232*100,yend=0.269*100))+
  scale_x_date(date_breaks="1.5 months",date_labels = "%b-%y",limits=c(as.Date("2020-04-01"),as.Date("2020-12-31")))+
  theme(axis.text.x = element_text(angle=45,hjust=1),legend.position="bottom",
        panel.grid.minor = element_blank())+
  guides(colour = "none",shape=guide_legend(order=2),
         fill = "none")+
  scale_colour_viridis_d(option="plasma",end=0.75,begin=0.3)+
  scale_fill_viridis_d(option="plasma",end=0.75,begin=0.3)+
  ylim(c(0,30))+
  labs(x="Date",y="IgG seroprevalence (%)",col="IgG seroreversion half-life",
       fill="IgG seroreversion half-life",shape=" ",tag="D")+
  scale_shape_manual(values=18)


aden_covid_sero_plot_iggm <- ggplot(aden_covid_sero %>% filter(antibody=="iggm",
                                                               igm_scenario=="50 days",
                                                               igg_scenario=="100 days"|igg_scenario=="280 days"),
                                   aes(x=date,y=med*100,ymin=min*100,ymax=max*100))+
  geom_rect(aes(xmin=as.Date("2020-11-28"),xmax=as.Date("2020-12-13"),ymin=-Inf,ymax=Inf),fill="#E8E8E8")+
  geom_line(aes(col=igg_scenario))+
  geom_ribbon(alpha=0.2,aes(fill=igg_scenario))+
  theme_bw()+
  geom_point(aes(x=as.Date("2020-12-06"),y=0.274*100,shape="Reported combined IgG/IgM"),size=3)+
  geom_segment(aes(x=as.Date("2020-12-06"),xend=as.Date("2020-12-06"),y=0.256*100,yend=0.293*100))+
  scale_x_date(date_breaks="1.5 months",date_labels = "%b-%y",limits=c(as.Date("2020-04-01"),as.Date("2020-12-31")))+
  theme(axis.text.x = element_text(angle=45,hjust=1),legend.position="bottom",
        panel.grid.minor = element_blank())+
  guides(colour = "none",shape=guide_legend(order=2),
         fill = "none")+
  scale_colour_viridis_d(option="plasma",end=0.75,begin=0.3)+
  scale_fill_viridis_d(option="plasma",end=0.75,begin=0.3)+
  labs(x="Date",y="Combined IgG/IgM \nseroprevalence (%)",col="IgG seroreversion half-life",
       fill="IgG seroreversion half-life",shape=" ",tag="C")+
  scale_shape_manual(values=15)



aden_covid_sero_plot_igm <- ggplot(aden_covid_sero %>% filter(antibody=="igm",
                                                              igm_scenario=="30 days"|igm_scenario=="70 days",
                                                              igg_scenario=="180 days") %>%
                                     mutate(igm_scenario=factor(igm_scenario,
                                                                levels=c("30 days","70 days"))),
                                    aes(x=date,y=med*100,ymin=min*100,ymax=max*100))+
  geom_rect(aes(xmin=as.Date("2020-11-28"),xmax=as.Date("2020-12-13"),ymin=-Inf,ymax=Inf),fill="#E8E8E8")+
  geom_line(aes(col=igm_scenario))+
  geom_ribbon(alpha=0.2,aes(fill=igm_scenario))+
  theme_bw()+
  geom_point(aes(x=as.Date("2020-12-06"),y=0.002*100,shape="Reported IgM"),size=3)+
  geom_segment(aes(x=as.Date("2020-12-06"),xend=as.Date("2020-12-06"),y=0.001*100,yend=0.004*100))+
  scale_x_date(date_breaks="1.5 months",date_labels = "%b-%y",limits=c(as.Date("2020-04-01"),as.Date("2020-12-31")))+
  theme(axis.text.x = element_text(angle=45,hjust=1),legend.position="bottom",
        panel.grid.minor = element_blank(),
        legend.box="vertical",
        legend.title.align=0.5)+
  guides(colour=guide_legend(title.position="top"),
    shape=guide_legend(order=2))+
  scale_colour_viridis_d(end=0.9,begin=0.3)+
  scale_fill_viridis_d(end=0.9,begin=0.3)+
  scale_y_continuous(limits=c(0,2))+
  labs(x="Date",y="IgM seroprevalence (%)",col="IgM seroreversion half-life",
       fill="IgM seroreversion half-life",shape=" ",tag="B")+
  scale_shape_manual(values=17)




## get the legend from somewhere and add that in

aden_covid_sero_legend <- ggplot(aden_covid_sero %>% filter(antibody=="iggm",
                                                            igm_scenario=="50 days",
                                                            igg_scenario=="100 days"|igg_scenario=="280 days"),
                                 aes(x=date,y=med*100,ymin=min*100,ymax=max*100))+
  geom_rect(aes(xmin=as.Date("2020-11-28"),xmax=as.Date("2020-12-13"),ymin=-Inf,ymax=Inf),fill="#E8E8E8")+
  geom_line(aes(col=igg_scenario))+
  geom_ribbon(alpha=0.2,aes(fill=igg_scenario))+
  theme_bw()+
  #geom_point(aes(x=as.Date("2020-12-06"),y=0.274*100,shape="Reported combined IgG/IgM"),size=3)+
  #geom_segment(aes(x=as.Date("2020-12-06"),xend=as.Date("2020-12-06"),y=0.256*100,yend=0.293*100))+
  scale_x_date(date_breaks="1.5 months",date_labels = "%b-%y",limits=c(as.Date("2020-04-01"),as.Date("2020-12-31")))+
  theme(axis.text.x = element_text(angle=45,hjust=1),legend.position="bottom",
        panel.grid.minor = element_blank(),
        legend.title.align=0.5)+
  guides(colour = guide_legend(title.position="top"))+
  scale_colour_viridis_d(option="plasma",end=0.75,begin=0.3)+
  scale_fill_viridis_d(option="plasma",end=0.75,begin=0.3)+
  labs(x="Date",y="Combined IgG/IgM \nseroprevalence (%)",col="IgG seroreversion half-life",
       fill="IgG seroreversion half-life",shape=" ",tag="C")



plot_grid(plot_grid(aden_covid_deaths_plot,aden_covid_sero_plot_igm,
          aden_covid_sero_plot_iggm,aden_covid_sero_plot_igg,
          nrow=2,#rel_heights = c(1,1.1),
          align="hv",axis="tb"),
          get_legend(aden_covid_sero_legend),nrow=2,rel_heights=c(1,0.1))
ggsave("analysis/figures/supplementary/aden_covid.png",height=8,width=7)


#####

aden_excess_sero_igm <- readRDS("analysis/data/derived/seroprevalence/Aden/aden_sero_excess_time_series.RDS") %>% filter(antibody=="igm",
                                                                                                                         ifr==0.3) %>%
  mutate(igm_scenario=factor(igm_scenario,
                             levels=c("5 days","10 days","15 days","20 days","25 days","30 days","35 days","40 days","45 days",
                                      "50 days","55 days","60 days","65 days","70 days")))

plot_excess_sero_a <- ggplot(aden_excess_sero_igm %>% filter(igm_scenario %in% c("30 days","35 days","40 days","45 days",
                                                                                 "50 days","55 days","60 days","65 days","70 days")),
                             aes(x=date,y=med*100,ymin=min*100,ymax=max*100))+
  geom_rect(aes(xmin=as.Date("2020-11-28"),xmax=as.Date("2020-12-13"),ymin=-Inf,ymax=Inf),fill="#E8E8E8")+
  geom_line(aes(col=igm_scenario))+
  geom_ribbon(alpha=0.2,aes(fill=igm_scenario))+
  theme_bw()+
  geom_point(aes(x=as.Date("2020-12-06"),y=0.002*100,shape="Reported IgM"),size=3)+
  geom_segment(aes(x=as.Date("2020-12-06"),xend=as.Date("2020-12-06"),y=0.001*100,yend=0.004*100))+
  scale_x_date(date_breaks="1.5 months",date_labels = "%b-%y",limits=c(as.Date("2020-04-01"),as.Date("2020-12-31")))+
  theme(axis.text.x = element_text(angle=45,hjust=1),legend.position="bottom",
        panel.grid.minor = element_blank(),
        legend.box="vertical",
        legend.title.align=0.5)+
  guides(colour=guide_legend(title.position="top"),
         shape=guide_legend(order=2))+
  scale_colour_viridis_d(end=0.9,begin=0.3)+
  scale_fill_viridis_d(end=0.9,begin=0.3)+
  labs(x="Date",y="IgM seroprevalence (%)",col="IgM seroreversion half-life",
       fill="IgM seroreversion half-life",shape=" ",tag="A")+
  scale_shape_manual(values=17)



plot_excess_sero_b <- ggplot(aden_excess_sero_igm %>% filter(igm_scenario %in% c("30 days","35 days","40 days","45 days",
                                                                                 "50 days","55 days","60 days","65 days","70 days")),
                             aes(x=date,y=med*100,ymin=min*100,ymax=max*100))+
  geom_rect(aes(xmin=as.Date("2020-11-28"),xmax=as.Date("2020-12-13"),ymin=-Inf,ymax=Inf),fill="#E8E8E8")+
  geom_line(aes(col=igm_scenario))+
  geom_ribbon(alpha=0.2,aes(fill=igm_scenario))+
  theme_bw()+
  geom_point(aes(x=as.Date("2020-12-06"),y=0.002*100,shape="Reported IgM"),size=3)+
  geom_segment(aes(x=as.Date("2020-12-06"),xend=as.Date("2020-12-06"),y=0.001*100,yend=0.004*100))+
  scale_x_date(date_breaks="1.5 months",date_labels = "%b-%y",limits=c(as.Date("2020-11-01"),as.Date("2020-12-31")))+
  scale_y_continuous(limits=c(0,3))+
  theme(axis.text.x = element_text(angle=45,hjust=1),legend.position="bottom",
        panel.grid.minor = element_blank(),
        legend.box="vertical",
        legend.title.align=0.5)+
  guides(colour=guide_legend(title.position="top"),
         shape=guide_legend(order=2))+
  scale_colour_viridis_d(end=0.9,begin=0.3)+
  scale_fill_viridis_d(end=0.9,begin=0.3)+
  labs(x="Date",y="IgM seroprevalence (%)",col="IgM seroreversion half-life",
       fill="IgM seroreversion half-life",shape=" ",tag="B")+
  scale_shape_manual(values=17)



aden_excess_sero_iggm <- readRDS("analysis/data/derived/seroprevalence/Aden/aden_sero_excess_time_series.RDS") %>%
  filter(antibody=="iggm",
         ifr==0.3,
         igg_scenario=="180 days") %>%
  mutate(igm_scenario=factor(igm_scenario,
                             levels=c("5 days","10 days","15 days","20 days","25 days","30 days","35 days","40 days","45 days",
                                      "50 days","55 days","60 days","65 days","70 days")))




plot_excess_sero_c <- ggplot(aden_excess_sero_iggm %>% filter(igm_scenario %in% c("30 days","35 days","40 days","45 days",
                                                                                  "50 days","55 days","60 days","65 days","70 days")),
                             aes(x=date,y=med*100,ymin=min*100,ymax=max*100))+
  geom_rect(aes(xmin=as.Date("2020-11-28"),xmax=as.Date("2020-12-13"),ymin=-Inf,ymax=Inf),fill="#E8E8E8")+
  geom_ribbon(alpha=0.05,aes(fill=igm_scenario))+
  geom_line(aes(col=igm_scenario))+
  theme_bw()+
  geom_point(aes(x=as.Date("2020-12-06"),y=0.274*100,shape="Reported combined IgG/IgM"),size=3)+
  geom_segment(aes(x=as.Date("2020-12-06"),xend=as.Date("2020-12-06"),y=0.256*100,yend=0.293*100))+
  scale_x_date(date_breaks="1.5 months",date_labels = "%b-%y",limits=c(as.Date("2020-04-01"),as.Date("2020-12-31")))+
  theme(axis.text.x = element_text(angle=45,hjust=1),legend.position="bottom",
        panel.grid.minor = element_blank(),
        legend.box="vertical",
        legend.title.align=0.5)+
  guides(colour=guide_legend(title.position="top"),
         shape=guide_legend(order=2))+
  scale_colour_viridis_d(end=0.9,begin=0.3)+
  scale_fill_viridis_d(end=0.9,begin=0.3)+
  labs(x="Date",y="Combined IgG/IgM seroprevalence (%)",col="IgM seroreversion half-life",
       fill="IgM seroreversion half-life",shape=" ",tag="C")+
  scale_shape_manual(values=15)

# plot_grid(plot_grid(plot_excess_sero_a+guides(colour="none",fill="none"),plot_excess_sero_b+guides(colour="none",fill="none"),nrow=1),
#           get_legend(plot_excess_sero_a+guides(shape="none")),
#           nrow=2,rel_heights=c(1,0.1))


plot_grid(plot_excess_sero_a+guides(colour="none",fill="none"),plot_excess_sero_b+guides(colour="none",fill="none"),
          plot_excess_sero_c+guides(colour="none",fill="none"),
          get_legend(plot_excess_sero_a+guides(shape="none",fill=guide_legend(ncol=2),
                                               colour=guide_legend(ncol=2))+theme(legend.position = "right")))
ggsave("analysis/figures/supplementary/aden_igm_sero.png",height=6,width=5)




# excess fit

aden_fit_ifr_raw <- readRDS("analysis/data/derived/model_fits/Aden/aden_excess_fit_complete.RDS")
aden_fit_ifr_draws <-  generate_draws(aden_fit_ifr_raw,generate_parameters(aden_fit_ifr_raw,draws=100),parallel=TRUE)

aden_deaths_ifr_plot <- plot(aden_fit_ifr_draws,var_select="deaths",particle_fit = TRUE)+
  theme(axis.text.x = element_text(angle=45,hjust=1),legend.position="none",
        panel.grid.minor = element_blank()) +
  geom_vline(xintercept=as.Date("2020-09-16"),linetype="dashed")+
  labs(tag="A") + scale_x_date(date_labels = "%b-%y",breaks="1.5 months")+
  scale_y_continuous(breaks = c(0,5,10,15,20,25))


aden_fit_ifr <- squire::projections(aden_fit_ifr_draws,time_period=125)

aden_deaths_proj_ifr_plot <- squire::projection_plotting(r_list=list(aden_fit_ifr,
                                                                aden_fit_ifr_draws),
                                                      scenarios=c("Model fit","Projection to end of 2020"),
                                                      var_select="deaths",
                                                      date_0 = as.Date("2020-09-16"),x_var="date")+
  scale_x_date(date_labels = "%b-%y",breaks="1.5 months")+
  scale_fill_manual(values=c("grey","red"))+
  scale_colour_manual(values=c("grey","red"))+
  theme(legend.position="none",
        axis.text.x = element_text(angle=45,hjust=1),
        panel.grid.minor = element_blank())+
  geom_vline(xintercept=as.Date("2020-09-16"),linetype="dashed")+
  labs(x="Date",y="Deaths",tag="B")+
  scale_y_continuous(breaks = c(0,5,10,15,20,25),limits=c(0,25))

plot_grid(aden_deaths_ifr_plot,aden_deaths_proj_ifr_plot,nrow=1,align="v")
ggsave("analysis/figures/supplementary/aden_excess.png",height=3,width=7)

# IFR 0.2

aden_fit_ifr02_raw <- readRDS("analysis/data/derived/model_fits/Aden/aden_excess_fit_complete_IFR02.RDS")
aden_fit_ifr02_draws <-  generate_draws(aden_fit_ifr02_raw,generate_parameters(aden_fit_ifr02_raw,draws=100),parallel=TRUE)

aden_deaths_ifr02_plot <- plot(aden_fit_ifr02_draws,var_select="deaths",particle_fit = TRUE)+
  theme(axis.text.x = element_text(angle=45,hjust=1),legend.position="none",
        panel.grid.minor = element_blank()) +
  labs(tag="A") + scale_x_date(date_labels = "%b-%y",breaks="1.5 months")+
  scale_y_continuous(breaks=c(0,5,10,15,20,25),limits=c(0,28))

#
# aden_fit_ifr02 <- squire::projections(aden_fit_ifr02_draws,time_period=125)
#
# aden_deaths_ifr02_plot <- squire::projection_plotting(r_list=list(aden_fit_ifr02,
#                                                                   aden_fit_ifr02_draws),
#                                                       scenarios=c("Model fit","Projection to end of 2020"),
#                                                       var_select="deaths",
#                                                       date_0 = as.Date("2020-09-16"),x_var="date")+
#   scale_x_date(date_labels = "%b-%y",breaks="1.5 months")+
#   scale_fill_manual(values=c("grey","red"))+
#   scale_colour_manual(values=c("grey","red"))+
#   theme(legend.position="none",
#         axis.text.x = element_text(angle=45,hjust=1),
#         panel.grid.minor = element_blank())+
#   geom_vline(xintercept=as.Date("2020-09-16"),linetype="dashed")+
#   labs(x="Date",y="Deaths")


# IFR 0.3

aden_fit_ifr03_raw <- readRDS("analysis/data/derived/model_fits/Aden/aden_excess_fit_complete_IFR03.RDS")
aden_fit_ifr03_draws <-  generate_draws(aden_fit_ifr03_raw,generate_parameters(aden_fit_ifr03_raw,draws=100),parallel=TRUE)

aden_deaths_ifr03_plot <- plot(aden_fit_ifr03_draws,var_select="deaths",particle_fit = TRUE)+
  theme(axis.text.x = element_text(angle=45,hjust=1),legend.position="none",
        panel.grid.minor = element_blank()) +
  labs(tag="B") + scale_x_date(date_labels = "%b-%y",breaks="1.5 months")+
  scale_y_continuous(breaks=c(0,5,10,15,20,25),limits=c(0,28))



# IFR 0.4

aden_fit_ifr04_raw <- readRDS("analysis/data/derived/model_fits/Aden/aden_excess_fit_complete_IFR04.RDS")
aden_fit_ifr04_draws <-  generate_draws(aden_fit_ifr04_raw,generate_parameters(aden_fit_ifr04_raw,draws=100),parallel=TRUE)

aden_deaths_ifr04_plot <- plot(aden_fit_ifr04_draws,var_select="deaths",particle_fit = TRUE)+
  theme(axis.text.x = element_text(angle=45,hjust=1),legend.position="none",
        panel.grid.minor = element_blank()) +
  labs(tag="C") + scale_x_date(date_labels = "%b-%y",breaks="1.5 months")+
  scale_y_continuous(breaks=c(0,5,10,15,20,25),limits=c(0,28))


# IFR 0.5

aden_fit_ifr05_raw <- readRDS("analysis/data/derived/model_fits/Aden/aden_excess_fit_complete_IFR05.RDS")
aden_fit_ifr05_draws <-  generate_draws(aden_fit_ifr05_raw,generate_parameters(aden_fit_ifr05_raw,draws=100),parallel=TRUE)

aden_deaths_ifr05_plot <- plot(aden_fit_ifr05_draws,var_select="deaths",particle_fit = TRUE)+
  theme(axis.text.x = element_text(angle=45,hjust=1),legend.position="none",
        panel.grid.minor = element_blank()) +
  labs(tag="D") + scale_x_date(date_labels = "%b-%y",breaks="1.5 months")+
  scale_y_continuous(breaks=c(0,5,10,15,20,25),limits=c(0,28))

# bring together

plot_grid(aden_deaths_ifr02_plot,aden_deaths_ifr03_plot,aden_deaths_ifr04_plot,aden_deaths_ifr05_plot,nrow=1,align="v")
ggsave("analysis/figures/supplementary/aden_IFR_fits.png",height=3,width=10)




