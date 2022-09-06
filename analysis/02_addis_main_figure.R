# Plot results of Addis Ababa analysis

library(tidyverse)
library(cowplot)


# (A) Estimated excess vs reported COVID-19 deaths over time

addis_covid <- readRDS("analysis/data/derived/deaths_time_series/addis_covid_deaths.RDS") %>%
  filter(date >= as.Date("2020-04-05")&date < as.Date("2020-11-01"))

addis_excess <- readRDS("analysis/data/derived/deaths_time_series/addis_excess_deaths.RDS") %>%
  select(dateburialgc,weekly,baseline) %>% rename(date=dateburialgc,deaths=weekly) %>%
  mutate(deaths=round(deaths),date=as.Date(date,format="%Y-%m-%d")) %>%
  filter(date >= as.Date("2020-04-01")&date < as.Date("2020-11-01"))

addis_deaths <- merge(addis_excess %>% rename(excess_deaths=deaths),
                      addis_covid %>% rename(covid_deaths=deaths),
                      by="date") %>%
  mutate(baseline_label = ifelse(baseline=="only_2019","2019 baseline","2015 - 2019 baseline"))


## extract legend from here
addis_plot_a_legend <- ggplot(addis_deaths)+
  #geom_rect(aes(xmin=as.Date("2020-05-09"),xmax=as.Date("2020-06-09"),ymin=-Inf,ymax=Inf),fill="#E8E8E8")+
  geom_line(aes(x=date,y=excess_deaths,group=baseline_label,col="Estimated excess deaths"))+
  facet_wrap(~baseline_label,nrow=1)+
  geom_line(aes(x=date,y=covid_deaths,col="Reported COVID-19 deaths"))+
  theme_bw()+
  theme(legend.position="top",strip.background = element_rect(fill="white"),
        axis.text.x = element_text(angle=45,hjust=1),
        panel.grid.minor = element_blank())+
  scale_colour_manual(values=c("cornflowerblue","#0818A8"))+
  scale_x_date(date_labels = "%b-%y",breaks="1.5 months")+
  geom_vline(xintercept=as.Date("2020-06-05"),linetype="dashed")+
  #geom_hline(yintercept=0,linetype="dashed")+
  labs(x="Date",y="Daily deaths",col="",tag="A")

addis_plot_a <- ggplot(addis_deaths %>% filter(baseline_label=="2015 - 2019 baseline"))+
  #geom_rect(aes(xmin=as.Date("2020-05-09"),xmax=as.Date("2020-06-09"),ymin=-Inf,ymax=Inf),fill="#E8E8E8")+
  geom_line(aes(x=date,y=excess_deaths,group=baseline_label,col="Estimated excess deaths"))+
  facet_wrap(~baseline_label,nrow=1)+
  geom_line(aes(x=date,y=covid_deaths,col="Reported COVID-19 deaths"))+
  theme_bw()+
  theme(legend.position="none",strip.background = element_rect(fill="white"),
        axis.text.x = element_text(angle=45,hjust=1),
        panel.grid.minor = element_blank())+
  scale_colour_manual(values=c("cornflowerblue","#0818A8"))+
  scale_x_date(date_labels = "%b-%y",breaks="1.5 months")+
  scale_y_continuous(limits=c(-12,35))+
  geom_vline(xintercept=as.Date("2020-06-05"),linetype="dashed")+
  #geom_hline(yintercept=0,linetype="dashed")+
  labs(x="Date",y="Daily deaths",col="",tag="A")

addis_plot_b <- ggplot(addis_deaths %>% filter(baseline_label=="2019 baseline"))+
  #geom_rect(aes(xmin=as.Date("2020-05-09"),xmax=as.Date("2020-06-09"),ymin=-Inf,ymax=Inf),fill="#E8E8E8")+
  geom_line(aes(x=date,y=excess_deaths,group=baseline_label,col="Estimated excess deaths"))+
  facet_wrap(~baseline_label,nrow=1)+
  geom_line(aes(x=date,y=covid_deaths,col="Reported COVID-19 deaths"))+
  theme_bw()+
  theme(legend.position="none",strip.background = element_rect(fill="white"),
        axis.text.x = element_text(angle=45,hjust=1),
        panel.grid.minor = element_blank())+
  scale_colour_manual(values=c("cornflowerblue","#0818A8"))+
  scale_x_date(date_labels = "%b-%y",breaks="1.5 months")+
  scale_y_continuous(limits=c(-12,35))+
  geom_vline(xintercept=as.Date("2020-06-05"),linetype="dashed")+
  #geom_hline(yintercept=0,linetype="dashed")+
  labs(x="Date",y="Daily deaths",col="",tag="B")




# (B) Estimated seroprevalence over time from model fit to excess mortality with 2015 - 2019 baseline

excess_all_years_April_sero <- readRDS("analysis/data/derived/seroprevalence/Addis/addis_excess_sero_all_years_April.RDS")

abdella <- data.frame(#"date"=rep(as.Date("2020-07-31"),2),
  "date"= c(as.Date("2020-07-27"),as.Date("2020-08-04")),
                      "med"=c(1.9,3.5),
                      "antibody"=factor(c("IgG","Combined IgG/IgM"),
                                        levels=c("IgG","Combined IgG/IgM")))

# abdella <- data.frame("date"=c(as.Date("2020-07-24"),as.Date("2020-08-03"),
#                                as.Date("2020-07-29"),as.Date("2020-08-08")),
#                       "med"=c(1.9,3.5,3.2,4.7),
#                       "antibody" = factor(c("IgG Weighted","Combined IgG/IgM Weighted",
#                                             "IgG Unweighted","Combined IgG/IgM Unweighted"),
#                                           levels=c("IgG Weighted","Combined IgG/IgM Weighted",
#                                                    "IgG Unweighted","Combined IgG/IgM Unweighted")))

# addis_plot_c <- ggplot(excess_all_years_April_sero %>% filter(date>=as.Date("2020-04-01")))+
#   theme_bw()+
#   geom_rect(aes(xmin=as.Date("2020-07-22"),xmax=as.Date("2020-08-10"),ymin=-Inf,ymax=Inf),fill="#E8E8E8")+
#   geom_ribbon(aes(x=date,ymin=min*100,ymax=max*100,group=antibody,fill=antibody),alpha=0.2)+
#   geom_line(aes(x=date,y=med*100,group=antibody,col=antibody),lwd=1)+
#   geom_segment(aes(x=as.Date("2020-07-24"),xend=as.Date("2020-07-24"),y=0.004*100,yend=0.037*100))+
#   geom_segment(aes(x=as.Date("2020-07-29"),xend=as.Date("2020-07-29"),y=0.022*100,yend=0.046*100))+
#   # geom_segment(aes(x=as.Date("2020-07-22"),xend=as.Date("2020-07-31"),y=0.037*100,yend=0.037*100))+
#   # geom_segment(aes(x=as.Date("2020-07-22"),xend=as.Date("2020-07-31"),y=0.004*100,yend=0.004*100))+
#   geom_segment(aes(x=as.Date("2020-08-03"),xend=as.Date("2020-08-03"),y=0.017*100,yend=0.054*100))+
#   geom_segment(aes(x=as.Date("2020-08-08"),xend=as.Date("2020-08-08"),y=0.035*100,yend=0.062*100))+
#   # geom_segment(aes(x=as.Date("2020-07-31"),xend=as.Date("2020-08-10"),y=0.054*100,yend=0.054*100))+
#   # geom_segment(aes(x=as.Date("2020-07-31"),xend=as.Date("2020-08-10"),y=0.017*100,yend=0.017*100))+
#   geom_point(data=abdella,aes(x=date,y=med,shape=antibody),size=1.5)+
#   labs(fill="Modelled",x="Date",y="Seroprevalence (%)",shape="Observed",
#        col="Modelled",tag="C")+
#   scale_shape_manual(values=c(17,15,2,0))+
#   guides(shape=guide_legend(nrow=2,byrow=TRUE,title.position = "top"),
#          col=guide_legend(nrow=2,byrow=TRUE,title.position = "top"),
#          fill=guide_legend(nrow=2,byrow=TRUE,title.position = "top"))+
#   scale_x_date(date_labels = "%b-%y",breaks="1.5 months")+
#   scale_y_continuous(breaks=c(0,2,4,6,8,10),limits=c(0,10))+
#   scale_colour_viridis_d(begin=0.6,end=0.8)+
#   scale_fill_viridis_d(begin=0.6,end=0.8)+
#   theme(axis.text.x = element_text(angle=45,hjust=1),
#         legend.position="none",
#         panel.grid.minor = element_blank())


addis_plot_c <- ggplot(excess_all_years_April_sero %>% filter(date>=as.Date("2020-04-01")&date<as.Date("2020-11-01")))+
  theme_bw()+
  geom_rect(aes(xmin=as.Date("2020-07-22"),xmax=as.Date("2020-08-10"),ymin=-Inf,ymax=Inf),fill="#E8E8E8")+
  geom_ribbon(aes(x=date,ymin=min*100,ymax=max*100,group=antibody,fill=antibody),alpha=0.2)+
  geom_line(aes(x=date,y=med*100,group=antibody,col=antibody),lwd=1)+
  geom_segment(aes(x=as.Date("2020-07-27"),xend=as.Date("2020-07-27"),y=0.004*100,yend=0.037*100))+
  #geom_segment(aes(x=as.Date("2020-07-22"),xend=as.Date("2020-07-31"),y=0.037*100,yend=0.037*100))+
  #geom_segment(aes(x=as.Date("2020-07-22"),xend=as.Date("2020-07-31"),y=0.004*100,yend=0.004*100))+
  geom_segment(aes(x=as.Date("2020-08-04"),xend=as.Date("2020-08-04"),y=0.017*100,yend=0.054*100))+
  #geom_segment(aes(x=as.Date("2020-08-01"),xend=as.Date("2020-08-10"),y=0.054*100,yend=0.054*100))+
  #geom_segment(aes(x=as.Date("2020-08-01"),xend=as.Date("2020-08-10"),y=0.017*100,yend=0.017*100))+
  geom_point(data=abdella,aes(x=date,y=med,shape=antibody),size=3)+
  labs(fill="Modelled",x="Date",y="Seroprevalence (%)",shape="Observed",
       col="Modelled",tag="C")+
  scale_shape_manual(values=c(18,15))+
  guides(shape=guide_legend(nrow=2,byrow=TRUE,title.position = "top"),
         col=guide_legend(nrow=2,byrow=TRUE,title.position = "top"),
         fill=guide_legend(nrow=2,byrow=TRUE,title.position = "top"))+
  scale_x_date(date_labels = "%b-%y",breaks="1.5 months")+
  scale_y_continuous(breaks=c(0,2,4,6,8,10),limits=c(0,10))+
  scale_colour_viridis_d(begin=0.6,end=0.8)+
  scale_fill_viridis_d(begin=0.6,end=0.8)+
  theme(axis.text.x = element_text(angle=45,hjust=1),
        legend.position="none",
        panel.grid.minor = element_blank())

# (c) Model-estimated seroprevalence under difference scenarios

model_estimates <- readRDS("analysis/data/derived/seroprevalence/Addis/addis_seroprevalence_model_estimates.RDS") %>%
  mutate(model = ifelse(model=="Model predicted (COVID-19)","Modelled (COVID-19)",
                        ifelse(model=="Model predicted (2015 - 19 baseline)","Modelled (2015 - 19 baseline)",
                               ifelse(model=="Model predicted (2015 - 19 baseline without 1st peak)",
                                      "Modelled (2015 - 19 baseline without 1st peak)",
                                      ifelse(model=="Model predicted (2019 baseline)","Modelled (2019 baseline)",
                                             ifelse(model=="Model predicted (2019 baseline without 1st peak)",
                                                    "Modelled (2019 baseline without 1st peak)",
                                                    ifelse(model=="Observed (Abdella et al.)","Reported (weighted)",
                                                           NA))))))) %>%
  mutate(median = median*100,
         lower = lower*100,
         upper = upper*100)

# model_estimates <- rbind(model_estimates,
#                          setNames(data.frame("IgG",3.2,2.2,4.6,"Abdella et al. (unweighted)"),
#                                   names(model_estimates)),
#                          setNames(data.frame("Combined IgG/IgM",4.7,3.5,6.2,"Abdella et al. (unweighted)"),
#                                   names(model_estimates))) #%>%
#  # mutate(median=median*100,lower=lower*100,upper=upper*100)

model_estimates$model <- factor(model_estimates$model,
                                levels=c(#"Abdella et al. (unweighted)",
                                  "Reported (weighted)",
                                         "Modelled (COVID-19)",
                                         "Modelled (2015 - 19 baseline)", "Modelled (2015 - 19 baseline without 1st peak)",
                                         "Modelled (2019 baseline)","Modelled (2019 baseline without 1st peak)"))
model_estimates$antibody <- factor(model_estimates$antibody,levels=c("IgG","Combined IgG/IgM"))

model_estimates <- model_estimates %>%
  mutate(colour_label = ifelse(model=="Reported (weighted)"&antibody=="IgG","Abdella et al. weighted IgG",
                               ifelse(model=="Reported (weighted)"&antibody=="Combined IgG/IgM",
                                      "Abdella et al. weighted combined IgG/IgM",
                                      #ifelse(model=="Abdella et al. (unweighted)"&antibody=="IgG",
                                             #"Abdella et al. unweighted IgG",
                                             #ifelse(model=="Abdella et al. (unweighted)"&antibody=="Combined IgG/IgM",
                                                    #"Abdella et al. unweighted combined IgG/IgM",
                                                    ifelse((model!="Abdella et al. (weighted)"&model!="Abdella et al. (unweighted)")&antibody=="IgG",
                                                           "Modelled IgG",
                                                           ifelse((model!="Abdella et al. (weighted)"&model!="Abdella et al. (unweighted)")&antibody=="Combined IgG/IgM",
                                                                  "Modelled combined IgG/IgM",NA)))))#))
model_estimates$colour_label <- factor(model_estimates$colour_label,
                                       levels=c("Modelled IgG","Modelled combined IgG/IgM",
                                                "Abdella et al. weighted IgG","Abdella et al. weighted combined IgG/IgM"#,
                                                #"Abdella et al. unweighted IgG","Abdella et al. unweighted combined IgG/IgM"
                                                ))




addis_plot_d <- ggplot(model_estimates)+
  #geom_segment(aes(x=model,xend=model,y=lower,yend=upper,group=antibody,col=colour_label),position=position_dodge(width=0.5))+
  #geom_errorbar(aes(x=model,ymin=lower,ymax=upper,group=antibody,col=colour_label),position=position_dodge(width=0.5))+
  geom_linerange(aes(x=model,ymin=lower,ymax=upper,group=antibody,col=colour_label),position=position_dodge(width=0.5))+
  geom_point(aes(x=model,y=median,shape=colour_label,group=antibody,col=colour_label,fill=colour_label),
             position=position_dodge(width=0.5),size=3)+
  theme_bw()+
  scale_shape_manual(values=c(18,15,18,15,2,0))+
  scale_colour_manual(values=c("#22a884","#7ad151","black","black","black","black"))+
  scale_fill_manual(values=c("#22a884","#7ad151","black","black","white","white"))+
  #scale_colour_viridis_d(option="magma",labels = function(x) str_wrap(x, width = 30),end=0.9)+
  labs(x="Model",y="Seroprevalence (%)",tag="D",col="",shape="")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 16))+
  scale_y_continuous(breaks=c(0,2,4,6,8,10),limits=c(0,10))+
  guides(colour=guide_legend(nrow=2))+
  theme(legend.position = "none",
        axis.text.x = element_text(angle=90,vjust=0.5),
        panel.grid.minor = element_blank())


# model_estimates <- model_estimates %>% mutate(legend_label = ifelse(colour_label=="Modelled IgG","IgG",
#                                                                     ifelse(colour_label=="Modelled combined IgG/IgM",
#                                                                            "Combined IgG/IgM",
#                                                                            ifelse(colour_label=="Abdella et al. weighted IgG",
#                                                                                   "Weighted IgG",
#                                                                                   ifelse(colour_label=="Abdella et al. weighted combined IgG/IgM",
#                                                                                          "Weighted combined IgG/IgM",
#                                                                                          ifelse(colour_label=="Abdella et al. unweighted IgG",
#                                                                                                 "Unweighted IgG",
#                                                                                                 ifelse(colour_label=="Abdella et al. unweighted combined IgG/IgM",
#                                                                                                        "Unweighted combined IgG/IgM",
#                                                                                                        NA)))))),
#                                               legend_label = factor(legend_label,
#                                                                     levels=c("IgG","Combined IgG/IgM",
#                                                                              "Weighted IgG","Weighted combined IgG/IgM",
#                                                                              "Unweighted IgG","Unweighted combined IgG/IgM")))


model_estimates <- model_estimates %>% mutate(legend_label2 = ifelse(colour_label=="Modelled IgG","Modelled IgG",
                                                                     ifelse(colour_label=="Modelled combined IgG/IgM",
                                                                            "Modelled combined IgG/IgM",
                                                                            ifelse(colour_label=="Abdella et al. weighted IgG",
                                                                                   "Reported IgG (weighted)",
                                                                                   ifelse(colour_label=="Abdella et al. weighted combined IgG/IgM",
                                                                                          "Reported combined IgG/IgM (weighted)",NA)))),
                                              legend_label2 = factor(legend_label2,
                                                                     levels=c("Modelled IgG",
                                                                              "Modelled combined IgG/IgM",
                                                                              "Reported IgG (weighted)",
                                                                              "Reported combined IgG/IgM (weighted)"
                                                                              )))



addis_plot_d_legend <- ggplot(model_estimates)+
  geom_point(aes(x=model,y=median,shape=legend_label2,group=antibody,col=legend_label2),
             position=position_dodge(width=0.5),size=3)+
  geom_errorbar(aes(x=model,ymin=lower,ymax=upper,group=antibody,col=legend_label2),position=position_dodge(width=0.5))+
  theme_bw()+
  scale_shape_manual(values=c(18,15,18,15,2,0))+
  scale_colour_manual(values=c("#22a884","#7ad151","black","black","black","black"))+
  #scale_colour_viridis_d(option="magma",labels = function(x) str_wrap(x, width = 30),end=0.9)+
  labs(x="Model",y="Seroprevalence (%)",tag="C",
       col="",
       shape="")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 16))+
  scale_y_continuous(breaks=c(0,2,4,6,8,10),limits=c(0,10))+
  guides(colour=guide_legend(nrow=2,title.position = "top"),
         shape=guide_legend(nrow=2,title.position = "top"))+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle=90,vjust=0.5),
        panel.grid.minor = element_blank())


# combine


#plot_grid(addis_plot_a,plot_grid(addis_plot_b,addis_plot_c,ncol=2,align="hv"),nrow=2,rel_heights = c(1,1.5))
#ggsave("analysis/figures/addis_figure_main.png",height=8,width=7)
#ggsave("analysis/figures/addis_figure_main_poster.png",height=7,width=7)

#
# plot_grid(addis_plot_a,
#           addis_plot_b,
#           addis_plot_c,
#           addis_plot_d,
#           nrow=2,
#           align="hv",axis="tb")


plot_grid(cowplot::get_legend(addis_plot_a_legend),
          plot_grid(addis_plot_a,
                    addis_plot_b,
                    addis_plot_c,
                    addis_plot_d,
                    nrow=2,
                    align="vh",axis="tb",
                    rel_heights=c(1,1.5)),
          cowplot::get_legend(addis_plot_d_legend),
          nrow=3,rel_heights=c(0.1,1,0.1))
ggsave("analysis/figures/addis_figure_main_updated.png",height=8,width=7)
### OJ, then in paint reduce the gap between the rows and sort the legend


