# plots of model fits and estimated seroprevalence for supplementary material for Addis

source("R/mcmc_utils.R")

# COVID-19 deaths

addis_covid_fit_raw <- readRDS("analysis/data/derived/model_fits/Addis/addis_covid_fit.RDS")
addis_covid_fit <-  generate_draws(addis_covid_fit_raw,generate_parameters(addis_covid_fit_raw,draws=100),parallel=TRUE)

covid_deaths_plot <- plot(addis_covid_fit,var_select="deaths",particle_fit = TRUE)+
  theme(axis.text.x = element_text(angle=45,hjust=1),legend.position="none",
        panel.grid.minor = element_blank()) +
  labs(tag="A") + scale_x_date(date_labels = "%b-%y",breaks="1.5 months")

addis_covid_sero <- readRDS("analysis/data/derived/seroprevalence/Addis/addis_covid_sero.RDS")

covid_sero_plot <- ggplot(addis_covid_sero)+
  theme_bw()+
  geom_ribbon(aes(x=date,ymin=min*100,ymax=max*100,group=antibody,fill=antibody),alpha=0.2)+
  geom_line(aes(x=date,y=med*100,group=antibody,col=antibody),lwd=1)+
  geom_segment(aes(x=as.Date("2020-07-31"),xend=as.Date("2020-07-31"),y=0.004*100,yend=0.037*100))+
  geom_segment(aes(x=as.Date("2020-07-22"),xend=as.Date("2020-08-10"),y=0.037*100,yend=0.037*100))+
  geom_segment(aes(x=as.Date("2020-07-22"),xend=as.Date("2020-08-10"),y=0.004*100,yend=0.004*100))+
  geom_segment(aes(x=as.Date("2020-07-31"),xend=as.Date("2020-07-31"),y=0.017*100,yend=0.054*100))+
  geom_segment(aes(x=as.Date("2020-07-22"),xend=as.Date("2020-08-10"),y=0.054*100,yend=0.054*100))+
  geom_segment(aes(x=as.Date("2020-07-22"),xend=as.Date("2020-08-10"),y=0.017*100,yend=0.017*100))+
  geom_point(aes(x=as.Date("2020-07-31"),y=0.019*100,shape="IgG"),size=5)+
  geom_point(aes(x=as.Date("2020-07-31"),y=0.035*100,shape="IgG + IgM"),size=5)+
  labs(fill="Model predicted",x="Date",y="Seroprevalence (%)",shape="Observed (Abdella et al.)",
       col="Model predicted",tag="B")+
  scale_shape_manual(values=c(18,15))+
  guides(shape=guide_legend(nrow=1,byrow=TRUE,title.position = "top",order=2),
         col=guide_legend(nrow=1,byrow=TRUE,title.position = "top",order=1),
         fill=guide_legend(nrow=1,byrow=TRUE,title.position = "top",order=1))+
  scale_x_date(date_labels = "%b-%y",breaks="1.5 months")+
  scale_y_continuous(breaks=c(0,2,4,6,8,10),limits=c(0,8))+
  scale_colour_viridis_d(begin=0.6,end=0.8)+
  scale_fill_viridis_d(begin=0.6,end=0.8)+
  theme(axis.text.x = element_text(angle=45,hjust=1),
        legend.position="bottom",
        panel.grid.minor = element_blank())

plot_grid(covid_deaths_plot,covid_sero_plot,nrow=2,rel_heights = c(1,1.5),align="v")
ggsave("analysis/figures/supplementary/addis_covid.png",height=6,width=7)


# excess deaths with 2015 - 2019 baseline

addis_excess_all_years_April_fit_raw <- readRDS("analysis/data/derived/model_fits/Addis/addis_excess_fit_all_years_April.RDS")
addis_excess_all_years_April_fit <-  generate_draws(addis_excess_all_years_April_fit_raw,
                                                    generate_parameters(addis_excess_all_years_April_fit_raw,
                                                                        draws=100),parallel=TRUE)

excess_deaths_all_years_April_plot <- plot(addis_excess_all_years_April_fit,var_select="deaths",particle_fit = TRUE)+
  theme(axis.text.x = element_text(angle=45,hjust=1),legend.position="none",
        panel.grid.minor = element_blank()) +
  labs(tag="A") + scale_x_date(date_labels = "%b-%y",breaks="1.5 months")

addis_all_years_April_sero <- readRDS("analysis/data/derived/seroprevalence/Addis/addis_excess_sero_all_years_April.RDS")

all_years_April_sero_plot <- ggplot(addis_all_years_April_sero)+
  theme_bw()+
  geom_ribbon(aes(x=date,ymin=min*100,ymax=max*100,group=antibody,fill=antibody),alpha=0.2)+
  geom_line(aes(x=date,y=med*100,group=antibody,col=antibody),lwd=1)+
  geom_segment(aes(x=as.Date("2020-07-31"),xend=as.Date("2020-07-31"),y=0.004*100,yend=0.037*100))+
  geom_segment(aes(x=as.Date("2020-07-22"),xend=as.Date("2020-08-10"),y=0.037*100,yend=0.037*100))+
  geom_segment(aes(x=as.Date("2020-07-22"),xend=as.Date("2020-08-10"),y=0.004*100,yend=0.004*100))+
  geom_segment(aes(x=as.Date("2020-07-31"),xend=as.Date("2020-07-31"),y=0.017*100,yend=0.054*100))+
  geom_segment(aes(x=as.Date("2020-07-22"),xend=as.Date("2020-08-10"),y=0.054*100,yend=0.054*100))+
  geom_segment(aes(x=as.Date("2020-07-22"),xend=as.Date("2020-08-10"),y=0.017*100,yend=0.017*100))+
  geom_point(aes(x=as.Date("2020-07-31"),y=0.019*100,shape="IgG"),size=5)+
  geom_point(aes(x=as.Date("2020-07-31"),y=0.035*100,shape="IgG + IgM"),size=5)+
  labs(fill="Model predicted",x="Date",y="Seroprevalence (%)",shape="Observed (Abdella et al.)",
       col="Model predicted",tag="C")+
  scale_shape_manual(values=c(18,15))+
  guides(shape=guide_legend(nrow=1,byrow=TRUE,title.position = "top",order=2),
         col=guide_legend(nrow=1,byrow=TRUE,title.position = "top",order=1),
         fill=guide_legend(nrow=1,byrow=TRUE,title.position = "top",order=1))+
  scale_x_date(date_labels = "%b-%y",breaks="1.5 months")+
  scale_y_continuous(breaks=c(0,2,4,6,8,10),limits=c(0,10))+
  scale_colour_viridis_d(begin=0.6,end=0.8)+
  scale_fill_viridis_d(begin=0.6,end=0.8)+
  theme(axis.text.x = element_text(angle=45,hjust=1),
        legend.position="bottom",
        panel.grid.minor = element_blank())

# excess deaths with 2015 - 2019 baseline + no first peak

addis_excess_all_years_June_fit_raw <- readRDS("analysis/data/derived/model_fits/Addis/addis_excess_fit_all_years_June.RDS")
addis_excess_all_years_June_fit <-  generate_draws(addis_excess_all_years_June_fit_raw,
                                                    generate_parameters(addis_excess_all_years_June_fit_raw,
                                                                        draws=100),parallel=TRUE)

excess_deaths_all_years_June_plot <- plot(addis_excess_all_years_June_fit,var_select="deaths",particle_fit = TRUE)+
  theme(axis.text.x = element_text(angle=45,hjust=1),legend.position="none",
        panel.grid.minor = element_blank()) +
  labs(tag="B") + scale_x_date(date_labels = "%b-%y",breaks="1.5 months")

addis_all_years_June_sero <- readRDS("analysis/data/derived/seroprevalence/Addis/addis_excess_sero_all_years_June.RDS")

all_years_June_sero_plot <- ggplot(addis_all_years_June_sero)+
  theme_bw()+
  geom_ribbon(aes(x=date,ymin=min*100,ymax=max*100,group=antibody,fill=antibody),alpha=0.2)+
  geom_line(aes(x=date,y=med*100,group=antibody,col=antibody),lwd=1)+
  geom_segment(aes(x=as.Date("2020-07-31"),xend=as.Date("2020-07-31"),y=0.004*100,yend=0.037*100))+
  geom_segment(aes(x=as.Date("2020-07-22"),xend=as.Date("2020-08-10"),y=0.037*100,yend=0.037*100))+
  geom_segment(aes(x=as.Date("2020-07-22"),xend=as.Date("2020-08-10"),y=0.004*100,yend=0.004*100))+
  geom_segment(aes(x=as.Date("2020-07-31"),xend=as.Date("2020-07-31"),y=0.017*100,yend=0.054*100))+
  geom_segment(aes(x=as.Date("2020-07-22"),xend=as.Date("2020-08-10"),y=0.054*100,yend=0.054*100))+
  geom_segment(aes(x=as.Date("2020-07-22"),xend=as.Date("2020-08-10"),y=0.017*100,yend=0.017*100))+
  geom_point(aes(x=as.Date("2020-07-31"),y=0.019*100,shape="IgG"),size=5)+
  geom_point(aes(x=as.Date("2020-07-31"),y=0.035*100,shape="IgG + IgM"),size=5)+
  labs(fill="Model predicted",x="Date",y="Seroprevalence (%)",shape="Observed (Abdella et al.)",
       col="Model predicted",tag="D")+
  scale_shape_manual(values=c(18,15))+
  guides(shape=guide_legend(nrow=1,byrow=TRUE,title.position = "top",order=2),
         col=guide_legend(nrow=1,byrow=TRUE,title.position = "top",order=1),
         fill=guide_legend(nrow=1,byrow=TRUE,title.position = "top",order=1))+
  scale_x_date(date_labels = "%b-%y",breaks="1.5 months")+
  scale_y_continuous(breaks=c(0,2,4,6,8,10),limits=c(0,10))+
  scale_colour_viridis_d(begin=0.6,end=0.8)+
  scale_fill_viridis_d(begin=0.6,end=0.8)+
  theme(axis.text.x = element_text(angle=45,hjust=1),
        legend.position="bottom",
        panel.grid.minor = element_blank())


plot_grid(plot_grid(excess_deaths_all_years_April_plot,excess_deaths_all_years_June_plot,
          all_years_April_sero_plot+theme(legend.position = "None"),
          all_years_June_sero_plot+theme(legend.position = "None"),
          ncol=2,rel_heights = c(1,1),align="v"),
          get_legend(all_years_April_sero_plot),rel_heights = c(1,0.1),nrow=2)
ggsave("analysis/figures/supplementary/addis_excess_all_years.png",height=6,width=7)


# excess deaths with 2019 baseline

addis_excess_only_2019_April_fit_raw <- readRDS("analysis/data/derived/model_fits/Addis/addis_excess_fit_only_2019_April.RDS")
addis_excess_only_2019_April_fit <-  generate_draws(addis_excess_only_2019_April_fit_raw,
                                                    generate_parameters(addis_excess_only_2019_April_fit_raw,
                                                                        draws=100),parallel=TRUE)

excess_deaths_only_2019_April_plot <- plot(addis_excess_only_2019_April_fit,var_select="deaths",particle_fit = TRUE)+
  theme(axis.text.x = element_text(angle=45,hjust=1),legend.position="none",
        panel.grid.minor = element_blank()) +
  labs(tag="A") + scale_x_date(date_labels = "%b-%y",breaks="1.5 months")

addis_only_2019_April_sero <- readRDS("analysis/data/derived/seroprevalence/Addis/addis_excess_sero_only_2019_April.RDS")

only_2019_April_sero_plot <- ggplot(addis_only_2019_April_sero)+
  theme_bw()+
  geom_ribbon(aes(x=date,ymin=min*100,ymax=max*100,group=antibody,fill=antibody),alpha=0.2)+
  geom_line(aes(x=date,y=med*100,group=antibody,col=antibody),lwd=1)+
  geom_segment(aes(x=as.Date("2020-07-31"),xend=as.Date("2020-07-31"),y=0.004*100,yend=0.037*100))+
  geom_segment(aes(x=as.Date("2020-07-22"),xend=as.Date("2020-08-10"),y=0.037*100,yend=0.037*100))+
  geom_segment(aes(x=as.Date("2020-07-22"),xend=as.Date("2020-08-10"),y=0.004*100,yend=0.004*100))+
  geom_segment(aes(x=as.Date("2020-07-31"),xend=as.Date("2020-07-31"),y=0.017*100,yend=0.054*100))+
  geom_segment(aes(x=as.Date("2020-07-22"),xend=as.Date("2020-08-10"),y=0.054*100,yend=0.054*100))+
  geom_segment(aes(x=as.Date("2020-07-22"),xend=as.Date("2020-08-10"),y=0.017*100,yend=0.017*100))+
  geom_point(aes(x=as.Date("2020-07-31"),y=0.019*100,shape="IgG"),size=5)+
  geom_point(aes(x=as.Date("2020-07-31"),y=0.035*100,shape="IgG + IgM"),size=5)+
  labs(fill="Model predicted",x="Date",y="Seroprevalence (%)",shape="Observed (Abdella et al.)",
       col="Model predicted",tag="C")+
  scale_shape_manual(values=c(18,15))+
  guides(shape=guide_legend(nrow=1,byrow=TRUE,title.position = "top",order=2),
         col=guide_legend(nrow=1,byrow=TRUE,title.position = "top",order=1),
         fill=guide_legend(nrow=1,byrow=TRUE,title.position = "top",order=1))+
  scale_x_date(date_labels = "%b-%y",breaks="1.5 months")+
  scale_y_continuous(breaks=c(0,4,8,12,16),limits=c(0,16))+
  scale_colour_viridis_d(begin=0.6,end=0.8)+
  scale_fill_viridis_d(begin=0.6,end=0.8)+
  theme(axis.text.x = element_text(angle=45,hjust=1),
        legend.position="bottom",
        panel.grid.minor = element_blank())

# excess deaths with 2019 baseline + no first peak

addis_excess_only_2019_June_fit_raw <- readRDS("analysis/data/derived/model_fits/Addis/addis_excess_fit_only_2019_June.RDS")
addis_excess_only_2019_June_fit <-  generate_draws(addis_excess_only_2019_June_fit_raw,
                                                   generate_parameters(addis_excess_only_2019_June_fit_raw,
                                                                       draws=100),parallel=TRUE)

excess_deaths_only_2019_June_plot <- plot(addis_excess_only_2019_June_fit,var_select="deaths",particle_fit = TRUE)+
  theme(axis.text.x = element_text(angle=45,hjust=1),legend.position="none",
        panel.grid.minor = element_blank()) +
  labs(tag="B") + scale_x_date(date_labels = "%b-%y",breaks="1.5 months")

addis_only_2019_June_sero <- readRDS("analysis/data/derived/seroprevalence/Addis/addis_excess_sero_only_2019_June.RDS")

only_2019_June_sero_plot <- ggplot(addis_only_2019_June_sero)+
  theme_bw()+
  geom_ribbon(aes(x=date,ymin=min*100,ymax=max*100,group=antibody,fill=antibody),alpha=0.2)+
  geom_line(aes(x=date,y=med*100,group=antibody,col=antibody),lwd=1)+
  geom_segment(aes(x=as.Date("2020-07-31"),xend=as.Date("2020-07-31"),y=0.004*100,yend=0.037*100))+
  geom_segment(aes(x=as.Date("2020-07-22"),xend=as.Date("2020-08-10"),y=0.037*100,yend=0.037*100))+
  geom_segment(aes(x=as.Date("2020-07-22"),xend=as.Date("2020-08-10"),y=0.004*100,yend=0.004*100))+
  geom_segment(aes(x=as.Date("2020-07-31"),xend=as.Date("2020-07-31"),y=0.017*100,yend=0.054*100))+
  geom_segment(aes(x=as.Date("2020-07-22"),xend=as.Date("2020-08-10"),y=0.054*100,yend=0.054*100))+
  geom_segment(aes(x=as.Date("2020-07-22"),xend=as.Date("2020-08-10"),y=0.017*100,yend=0.017*100))+
  geom_point(aes(x=as.Date("2020-07-31"),y=0.019*100,shape="IgG"),size=5)+
  geom_point(aes(x=as.Date("2020-07-31"),y=0.035*100,shape="IgG + IgM"),size=5)+
  labs(fill="Model predicted",x="Date",y="Seroprevalence (%)",shape="Observed (Abdella et al.)",
       col="Model predicted",tag="D")+
  scale_shape_manual(values=c(18,15))+
  guides(shape=guide_legend(nrow=1,byrow=TRUE,title.position = "top",order=2),
         col=guide_legend(nrow=1,byrow=TRUE,title.position = "top",order=1),
         fill=guide_legend(nrow=1,byrow=TRUE,title.position = "top",order=1))+
  scale_x_date(date_labels = "%b-%y",breaks="1.5 months")+
  scale_y_continuous(breaks=c(0,4,8,12,16),limits=c(0,16))+
  scale_colour_viridis_d(begin=0.6,end=0.8)+
  scale_fill_viridis_d(begin=0.6,end=0.8)+
  theme(axis.text.x = element_text(angle=45,hjust=1),
        legend.position="bottom",
        panel.grid.minor = element_blank())


plot_grid(plot_grid(excess_deaths_only_2019_April_plot,excess_deaths_only_2019_June_plot,
                    only_2019_April_sero_plot+theme(legend.position = "None"),
                    only_2019_June_sero_plot+theme(legend.position = "None"),
                    ncol=2,rel_heights = c(1,1),align="v"),
          get_legend(only_2019_April_sero_plot),rel_heights = c(1,0.1),nrow=2)
ggsave("analysis/figures/supplementary/addis_excess_only_2019.png",height=6,width=7)



