## chi-squared tests for stat significant differences between observed and estimated

reported_sero <- data.frame("antibody" = c("igg","iggm","igm"),
                            "unweighted_med"=c(25/100,27.4/100,0.2/100),
                            "unweighted_lower"=c(23.2/100,25.6/100,0.1/100),
                            "unweighted_upper"=c(26.9/100,29.3/100,0.4/100)) %>%
  mutate(unweighted_var = ((unweighted_upper - unweighted_lower)/3.92)^2)


cs_test <- function(est1,est2,var_est1,var_est2){
  test <- ((est1-est2)^2)/(var_est1+var_est2)
  p_val <- 1-pchisq(test,df=1)
  return(c("test"=test,
           "p-val"=p_val))
}

options(scipen=999)

#### excess mortality

excess_summary <- readRDS("analysis/data/derived/seroprevalence/Aden/aden_sero_excess_summary.rds")

chisq_output <- merge(excess_summary,reported_sero,by="antibody",all.x=TRUE) %>%
  group_by(antibody,ifr,igg_scenario,igm_scenario) %>%
  mutate(unweighted_p = cs_test(median,unweighted_med,var,unweighted_var)[2]) %>% data.frame()

chisq_output %>% filter(antibody=="igg") %>% dplyr::select(ifr,igg_scenario,unweighted_p) %>%
  mutate(unweighted_p = round(unweighted_p,3)) %>% arrange(ifr)

chisq_output %>% filter(antibody=="iggm",igm_scenario=="50 days") %>% dplyr::select(ifr,igg_scenario,unweighted_p) %>%
  mutate(unweighted_p = round(unweighted_p,3)) %>% arrange(ifr,igg_scenario)

chisq_output %>% filter(antibody=="igm") %>% dplyr::select(ifr,igm_scenario,unweighted_p) %>%
  mutate(unweighted_p = round(unweighted_p,3)) %>% arrange(ifr,igm_scenario) %>%
  filter(igm_scenario %in% c("30 days","35 days","40 days","45 days",
                             "50 days","55 days","60 days","65 days","70 days")) %>%
  mutate(igm_scenario = factor(igm_scenario,
                               levels=c("30 days","35 days","40 days","45 days",
                                        "50 days","55 days","60 days","65 days","70 days"))) %>% arrange(ifr,igm_scenario)


#### covid mortality (referred to in text only)

covid_summary <- readRDS("analysis/data/derived/seroprevalence/Aden/aden_sero_covid_summary.rds")

chisq_output_covid <- merge(covid_summary,reported_sero,by="antibody",all.x=TRUE) %>%
  group_by(antibody,igg_scenario,igm_scenario) %>%
  mutate(unweighted_p = cs_test(median,unweighted_med,var,unweighted_var)[2]) %>% data.frame()

chisq_output_covid %>% filter(antibody=="igg")
chisq_output_covid %>% filter(antibody=="iggm") #%>% summary()
chisq_output_covid %>% filter(antibody=="igm")



