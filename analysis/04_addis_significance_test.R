## chi-squared tests for stat significant differences between observed and estimated

library(tidyverse)

model_estimates <- readRDS("analysis/data/derived/seroprevalence/Addis/addis_seroprevalence_model_estimates.RDS") %>%
  mutate(model = ifelse(model=="Model predicted (COVID-19)","Modelled (COVID-19)",
                        ifelse(model=="Model predicted (2015 - 19 baseline)","Modelled (2015 - 19 baseline)",
                               ifelse(model=="Model predicted (2015 - 19 baseline without 1st peak)",
                                      "Modelled (2015 - 19 baseline without 1st peak)",
                                      ifelse(model=="Model predicted (2019 baseline)","Modelled (2019 baseline)",
                                             ifelse(model=="Model predicted (2019 baseline without 1st peak)",
                                                    "Modelled (2019 baseline without 1st peak)",
                                                    ifelse(model=="Observed (Abdella et al.)","Observed (Abdella et al.)",
                                                           NA)))))))
cs_test <- function(est1,est2,var_est1,var_est2){
  test <- ((est1-est2)^2)/(var_est1+var_est2)
  p_val <- 1-pchisq(test,df=1)
  return(c("test"=test,
           "p-val"=p_val))
}


model_estimates %>% filter(model!="Observed (Abdella et al.)",antibody=="IgG") %>%
  mutate(#var = ((upper-lower)/3.92)^2,
         weighted_med = 1.9/100,
         weighted_var = ((3.7/100-0.4/100)/3.92)^2,
         unweighted_med = 3.2/100,
         unweighted_var = ((4.6/100-2.2/100)/3.92)^2) %>% group_by(model) %>%
  mutate(weighted_p = cs_test(median,weighted_med,var,weighted_var)[2],
         unweighted_p = cs_test(median,unweighted_med,var,unweighted_var)[2]) %>%
  dplyr::select(antibody,model,weighted_p,unweighted_p) %>%
  mutate(weighted_p=round(weighted_p,3),
         unweighted_p=round(unweighted_p,3)) %>%
  data.frame()


model_estimates %>% filter(model!="Observed (Abdella et al.)",antibody=="Combined IgG/IgM") %>%
  mutate(var = ((upper-lower)/3.92)^2,
         weighted_med = 3.5/100,
         weighted_var = ((5.4/100-1.7/100)/3.92)^2,
         unweighted_med = 4.7/100,
         unweighted_var = ((6.2/100-3.5/100)/3.92)^2) %>% group_by(model) %>%
  mutate(weighted_p = cs_test(median,weighted_med,var,weighted_var)[2],
         unweighted_p = cs_test(median,unweighted_med,var,unweighted_var)[2]) %>%
  dplyr::select(antibody,model,weighted_p,unweighted_p) %>%
  mutate(weighted_p=round(weighted_p,3),
         unweighted_p=round(unweighted_p,3)) %>%
  data.frame()






igg <- model_estimates %>% filter(antibody=="IgG") %>%
  mutate(var = ((upper-lower)/3.92)^2)



cs_test(2.03,1.9,0.121,0.709)
cs_test(6.26,1.9,0.288,0.709)


# ## IgG example
#
# model_estimates %>% filter(antibody=="IgG")
#
# ## take abdella et al. as "observed"
#
# 0.019*956 # observed
# 0.023*956 # expected
#
# 1 - pchisq(((18-22)^2)/22,df=1)
#
# chisq.test(x=c(18,22))
#
#
# ## take modelled as observed
# 0.023*4800000 #observed
# 0.019*4800000 #expected
#
#
# 1 - pchisq(((110400-91200)^2)/91200,df=1)
# chisq.test(x=c(110400,91200))
#
#
# chisq.test(x=c(11,9))
#
#
# ### try reading in time series of one of them
#
# excess_all_years_April_sero <- readRDS("analysis/data/derived/seroprevalence/Addis/addis_excess_sero_all_years_April.RDS")
#
# data <- excess_all_years_April_sero %>% filter(date>=as.Date("2020-07-22"),
#                                                date<=as.Date("2020-08-10")) %>% filter(antibody!="IgG") %>%
#   mutate(obs = round(med*4800000),
#          exp = round(0.035*4800000))
#
# data %>% mutate(test_stat = (obs-exp)^2/exp) %>% dplyr::select(test_stat) %>% sum()
#
# data %>% mutate(test_stat = (med - 3.5)^2/3.5) %>% dplyr::select(test_stat) %>% sum()
#
#
# # model_estimates %>% filter(model!="Observed (Abdella et al.)"&antibody=="IgG") %>%
# #   mutate(observed = round(median/100*956),
# #          expected = round(0.019*956))
# #
# #
# # (97322 - 91200)^2/91200
# #
# # (2.03-1.9)^2/1.9
# # pchisq(410,1)
# #
# # t <- prop.test(x=c(97322,91200),n=c(4800000,4800000),alternative="two.sided")
# # t$statistic
# #
# # prop.test(x=c(19,18),n=c(956,956),alternative="two.sided")
#
#
#
# # ## chi square test stat
# # model_igg <- model_estimates %>% filter(model!="Observed (Abdella et al.)"&antibody=="IgG") %>%
# #   mutate(observed = 1.9)
# #
# # model_igg %>% mutate(test_stat = ((observed - median)^2)/median,
# #                      p_val = pchisq(test_stat,df=1))
# #
# #
# #
# # model_iggm <- model_estimates %>% filter(model!="Observed (Abdella et al.)"&antibody!="IgG") %>%
# #   mutate(observed = 3.5)
# #
# # model_iggm %>% mutate(test_stat = ((observed - median)^2)/median,
# #                      p_val = pchisq(test_stat,df=1))
# #
# # chisq.test(x=c(2.43,3.5))
# #
# #
# #
# #
# # ## try prop test?
# # model_igg %>% mutate(exp_n = 4800000,
# #                      exp_x = round(median * exp_n),
# #                      obs_n = 956,
# #                      obs_x = round(observed/100 * obs_n))
# #
# # prop.test(x=c(9732206,18),
# #           n=c(4800000,956))
# # ### can't do this as have one x greater than one n
# #
# #
# #
# # ### compare to other test
# # # exp <- 12862
# # # obs <- 11256
# #
