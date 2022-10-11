library(tidyverse)
devtools::load_all()

## Symptomatic infection prevalence estimate as of 2nd June which we compare to cum_hosp_inc + cum_icu_inc
# https://www.researchgate.net/publication/342122988_Estimation_of_Coronavirus_COVID-19_Infections_in_Khartoum_State_Sudan
p_e_low <- 0.083
p_e_med <- 0.11
p_e_high <- 0.137

## Omdurman deaths and seroprevalence
# https://www.medrxiv.org/content/10.1101/2021.08.22.21262294v1.full.pdf
omd_deaths_med <- 7113 / (3040604/8877146)
omd_deaths_low <- 5015 / (3040604/8877146)
omd_deaths_high <- 9505 / (3040604/8877146)

omd_sero_med <- 0.546
omd_sero_low <- 0.514
omd_sero_high <- 0.578

# -----------------------------------------------
## Fig a
# -----------------------------------------------

data <- readRDS("analysis/data/derived/deaths_time_series/khartoum_covid_deaths.rds")
gg_a <- ggplot(data, aes(as.Date(date), deaths)) +
  geom_point(aes(fill = "Daily"), shape = 21) +
  geom_line(aes(y = zoo::rollmean(deaths, 7, na.pad = TRUE), color = "Weekly Mean"), size = 1) +
  ylab("Daily Reported COVID-19 Deaths\n") +
  theme_bw() +
  xlab("Date") +
  ggpubr::theme_pubclean() +
  scale_x_date(date_breaks = "3 month", date_labels = "%b-%Y") +
  theme(axis.line = element_line(),
        panel.grid.major.x = element_blank()) +
  scale_fill_manual(name = "", values=c("Black"), label = c("Daily Deaths")) +
  scale_color_manual(name = "", values=c("#21448B"), label = c("Weekly Mean")) +
  theme(legend.position = c(0.85,0.9),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.spacing.y = unit(0.2, "mm")) +
  guides(fill = guide_legend(order = 1),
       color = guide_legend(order = 2))

# -----------------------------------------------
## Fig b
# -----------------------------------------------

ll <- readRDS("analysis/data/derived/seroprevalence/Khartoum/likelihood_rf.rds")

gg_b <- ll %>% group_by(rf) %>% summarise(across(ll_low:p_e_high, .fns = mean)) %>%
  ggplot(aes(rf, p_e_med, ymin = p_e_low, ymax = p_e_high)) +
  scale_x_continuous(labels = scales::percent, limits = c(0.02,0.1), expand = c(-0,0.001)) +
  scale_y_continuous(labels = scales::percent) +
  geom_rect(ymin = p_e_low, ymax = p_e_high,
            xmin = -Inf, xmax = Inf, fill = 'blue', alpha = 0.01) +
  geomtextpath::geom_texthline(yintercept = p_e_med, label = "Observed Cumulative Symptomatic Cases", hjust = 0.9, cex=2.5, linewidth = 0.5) +
  geomtextpath::geom_texthline(yintercept = p_e_low, label = "97.5%", hjust = 0.75, linetype = "dashed", cex=3, linewidth = 0.5) +
  geomtextpath::geom_texthline(yintercept = p_e_high, label = "2.5%", hjust = 0.75, linetype = "dashed", cex=3, linewidth = 0.5) +
  #geom_smooth() +
  geom_pointrange() +
  xlab("Reporting Fraction") +
  ylab("Cumulative proportion of symptomatic \ninfections in 15+ year olds\n") +
  ggpubr::theme_pubclean() +
  theme(axis.line = element_line())


# -----------------------------------------------
## Fig c
# -----------------------------------------------

res <- readRDS("analysis/data/derived/model_fits/Khartoum/khartoum_fit_rf_0.04.rds")
res <- squire::projections(res, 10)
difff <- function(x) {c(0, diff(x))}

# In Moser et al. they have already adjusted for test sensitivity and waning performance
# so just work out the proportion not in S lagged for seroconversion time
sero_det <- sero_det("iggm", igg_sens = 1, igg_conv = 13.3, igg_scale = Inf)

sero_df <- squire::format_output(res, "S", date_0 = max(res$pmcmc_results$inputs$data$date)) %>%
  group_by(date, compartment) %>%
  summarise(med = median(y, na.rm = TRUE),
            min = quantile(y, 0.025, na.rm = TRUE),
            max = quantile(y, 0.975, na.rm = TRUE),
            var = var(y, na.rm = TRUE)) %>%
  mutate(across(med:max, ~1-(.x/sum(res$parameters$population)))) %>%
  ungroup %>%
  mutate(across(med:max, difff)) %>%
  mutate(across(med:max, lag, 5, 0)) %>%
  mutate(across(med:max, roll_func, sero_det))

gg_c <- ggplot() +
  annotate("rect",
           xmin = as.Date("2021-03-01"), xmax = as.Date("2021-04-10"),
           ymin = -Inf, ymax = Inf,  fill = "grey", alpha=0.5) +
  geom_ribbon(aes(date, med, ymin = min, ymax = max), alpha = 0.2, fill = viridis::plasma(1, begin = 0.8, end = 0.8), data = sero_df) +
  geom_point(aes(date, med), size = 0, color = NA, data = sero_df) +
  geom_line(aes(date, med), color = viridis::plasma(1, begin = 0.8, end = 0.8), data = sero_df) +
  geom_point(aes(x = as.Date("2021-03-20"), y = omd_sero_med)) +
  geom_errorbar(aes(x = as.Date("2021-03-20"), y = omd_sero_med, ymin = omd_sero_low, ymax = omd_sero_high)) +
  #geom_errorbarh(aes(xmin = as.Date("2021-03-01"), xmax = as.Date("2021-04-10"), y = omd_sero_med), height = 0) +
  ggrepel::geom_text_repel(aes(x = as.Date("2021-03-20"), y = omd_sero_med, label = "Omdurman \nSeroprevalence"),
                           point.padding = 0.8,
                           nudge_x = 0.4,
                           nudge_y = -.08,
                           segment.curvature = -1e-20,
                           arrow = arrow(length = unit(0.015, "npc"))) +
  xlab("Date") +
  ylab("Adjusted Seroprevalence (%)") +
  scale_y_continuous(labels = scales::percent) +
  ggpubr::theme_pubclean() +
  theme(axis.line = element_line(),
        panel.grid.major.x = element_blank())


# -----------------------------------------------
## Fig d
# -----------------------------------------------

res <- readRDS("analysis/data/derived/model_fits/Khartoum/khartoum_fit_rf_0.04.rds")
res <- squire::projections(res, 10)

# In Moser et al. they have already adjusted for test sensitivity and waning performance
# so just work out the proportion not in S lagged for seroconversion time
sero_det <- sero_det("iggm", igg_sens = 1, igg_conv = 13.3, igg_scale = Inf)

D_df <- squire::format_output(res, "D", date_0 = max(res$pmcmc_results$inputs$data$date)) %>%
  group_by(date, compartment) %>%
  summarise(med = median(y, na.rm = TRUE),
            min = quantile(y, 0.025, na.rm = TRUE),
            max = quantile(y, 0.975, na.rm = TRUE))

# number for manuscript
D_df %>% filter(date == "2021-03-20")

gg_d <- ggplot() +
  annotate("rect",
           xmin = as.Date("2021-03-01"), xmax = as.Date("2021-04-10"),
           ymin = -Inf, ymax = Inf,  fill = "grey", alpha=0.5) +
  geom_ribbon(aes(date, med, ymin = min, ymax = max), alpha = 0.2, fill = viridis::plasma(1, begin = 0.4, end = 0.4), data = D_df) +
  geom_point(aes(date, med), size = 0, color = NA, data = D_df) +
  geom_line(aes(date, med), color = viridis::plasma(1, begin = 0.4, end = 0.4), data = D_df) +
  geom_point(aes(x = as.Date("2021-03-20"), y = omd_deaths_med)) +
  geom_errorbar(aes(x = as.Date("2021-03-20"), y = omd_deaths_med, ymin = omd_deaths_low, ymax = omd_deaths_high)) +
  #geom_errorbarh(aes(xmin = as.Date("2021-03-01"), xmax = as.Date("2021-04-10"), y = omd_deaths_med), height = 0) +
  ggrepel::geom_text_repel(aes(x = as.Date("2021-03-20"), y = omd_deaths_med, label = "Omdurman Excess Mortality"),
                           point.padding = 0.8,
                           nudge_x = -150,
                           nudge_y = 2000,
                           segment.curvature = -1e-20,
                           arrow = arrow(length = unit(0.015, "npc"))) +
  xlab("Date") +
  ylab("Cumulative COVID-19 Deaths") +
  ggpubr::theme_pubclean() +
  theme(axis.line = element_line(),
        panel.grid.major.x = element_blank())

# -----------------------------------------------
## Fig 3 alltogether
# -----------------------------------------------

fig3 <- cowplot::plot_grid(gg_a, NULL, gg_b, gg_c, NULL, gg_d,
                           labels = c("A","","B","C","","D"),
                           rel_widths = c(1,0.1,1)) +
  theme(plot.background = element_rect(fill = "white"))
save_figs("khartoum_figure_main", fig3, width = 10, height = 8, root = "analysis/figures/")


# -----------------------------------------------
## Supplement Best Fit to Moser
# -----------------------------------------------

res <- readRDS("analysis/data/derived/model_fits/Khartoum/khartoum_fit_rf_0.045.rds")
res <- squire::projections(res, 10)

# In Moser et al. they have already adjusted for test sensitivity and waning performance
# so just work out the proportion not in S lagged for seroconversion time
sero_det <- sero_det("iggm", igg_sens = 1, igg_conv = 13.3, igg_scale = Inf)

sero_df <- squire::format_output(res, "S", date_0 = max(res$pmcmc_results$inputs$data$date)) %>%
  group_by(date, compartment) %>%
  summarise(med = median(y, na.rm = TRUE),
            min = quantile(y, 0.025, na.rm = TRUE),
            max = quantile(y, 0.975, na.rm = TRUE)) %>%
  mutate(across(med:max, ~1-(.x/sum(res$parameters$population)))) %>%
  ungroup %>%
  mutate(across(med:max, difff)) %>%
  mutate(across(med:max, lag, 5, 0)) %>%
  mutate(across(med:max, roll_func, sero_det))

supp_2a <- ggplot() +
  annotate("rect",
           xmin = as.Date("2021-03-01"), xmax = as.Date("2021-04-10"),
           ymin = -Inf, ymax = Inf,  fill = "grey", alpha=0.5) +
  geom_ribbon(aes(date, med, ymin = min, ymax = max), alpha = 0.2, fill = viridis::plasma(1, begin = 0.8, end = 0.8), data = sero_df) +
  geom_point(aes(date, med), size = 0, color = NA, data = sero_df) +
  geom_line(aes(date, med), color = viridis::plasma(1, begin = 0.8, end = 0.8), data = sero_df) +
  geom_point(aes(x = as.Date("2021-03-20"), y = omd_sero_med)) +
  geom_errorbar(aes(x = as.Date("2021-03-20"), y = omd_sero_med, ymin = omd_sero_low, ymax = omd_sero_high)) +
  #geom_errorbarh(aes(xmin = as.Date("2021-03-01"), xmax = as.Date("2021-04-10"), y = omd_sero_med), height = 0) +
  ggrepel::geom_text_repel(aes(x = as.Date("2021-03-20"), y = omd_sero_med, label = "Omdurman \nSeroprevalence"),
                           point.padding = 0.8,
                           nudge_x = 0.4,
                           nudge_y = -.08,
                           segment.curvature = -1e-20,
                           arrow = arrow(length = unit(0.015, "npc"))) +
  xlab("Date") +
  ylab("Adjusted Seroprevalence (%)") +
  scale_y_continuous(labels = scales::percent) +
  ggpubr::theme_pubclean() +
  theme(axis.line = element_line(),
        panel.grid.major.x = element_blank())


D_df <- squire::format_output(res, "D", date_0 = max(res$pmcmc_results$inputs$data$date)) %>%
  group_by(date, compartment) %>%
  summarise(med = median(y, na.rm = TRUE),
            min = quantile(y, 0.025, na.rm = TRUE),
            max = quantile(y, 0.975, na.rm = TRUE))

supp_2b <- ggplot() +
  annotate("rect",
           xmin = as.Date("2021-03-01"), xmax = as.Date("2021-04-10"),
           ymin = -Inf, ymax = Inf,  fill = "grey", alpha=0.5) +
  geom_ribbon(aes(date, med, ymin = min, ymax = max), alpha = 0.2, fill = viridis::plasma(1, begin = 0.4, end = 0.4), data = D_df) +
  geom_point(aes(date, med), size = 0, color = NA, data = D_df) +
  geom_line(aes(date, med), color = viridis::plasma(1, begin = 0.4, end = 0.4), data = D_df) +
  geom_point(aes(x = as.Date("2021-03-20"), y = omd_deaths_med)) +
  geom_errorbar(aes(x = as.Date("2021-03-20"), y = omd_deaths_med, ymin = omd_deaths_low, ymax = omd_deaths_high)) +
  #geom_errorbarh(aes(xmin = as.Date("2021-03-01"), xmax = as.Date("2021-04-10"), y = omd_deaths_med), height = 0) +
  ggrepel::geom_text_repel(aes(x = as.Date("2021-03-20"), y = omd_deaths_med, label = "Omdurman Excess Mortality"),
                           point.padding = 0.8,
                           nudge_x = -150,
                           nudge_y = 2000,
                           segment.curvature = -1e-20,
                           arrow = arrow(length = unit(0.015, "npc"))) +
  xlab("Date") +
  ylab("Cumulative COVID-19 Deaths") +
  ggpubr::theme_pubclean() +
  theme(axis.line = element_line(),
        panel.grid.major.x = element_blank())

# Bring it together
supp_2 <- cowplot::plot_grid(supp_2a, NULL, supp_2b,
                           labels = c("A","","B"),
                           rel_widths = c(1,0.1,1), ncol = 3) +
  theme(plot.background = element_rect(fill = "white"))
save_figs("khartoum_best_fit", supp_2, width = 10, height = 4, root = "analysis/figures/supplementary")

# -----------------------------------------------
## Supplement Individual Model Fits
# -----------------------------------------------

fls <- list.files("analysis/data/derived/model_fits/Khartoum/", full.names = TRUE)
fls <- grep("fit_rf", fls, value = TRUE)
pl <- list()
pl2 <- list()
for(i in seq_along(fls)) {
  r <- readRDS(fls[i])
  rf <- r$pmcmc_results$inputs$pars_obs$phi_death
  pl[[i]] <- plot(r, particle_fit = TRUE) +
    geom_line(aes(date, deaths),
              r$pmcmc_results$inputs$data %>% mutate(deaths = zoo::rollmean(deaths, 7, na.pad = TRUE)/rf),
              color = "#21448B") +
    ggtitle(paste0("Reporting Fraction = ", rf)) + theme(legend.position = "none") +
    theme(axis.title.x = element_blank(), plot.margin = margin(10, 10, 10, 10)) +
    ylab("Daily Deaths") + scale_x_date(breaks = "4 months", date_labels = "%b %Y")
  pl[[i]]$layers <- pl[[i]]$layers[c(5,6,1,2,3,4)]

  pl2[[i]] <- plot(r, "D", date_0 = max(r$pmcmc_results$inputs$data$date), x_var = "date") +
    geom_line(aes(date, deaths), r$pmcmc_results$inputs$data %>% mutate(deaths = (cumsum(deaths/rf)))) +
    ggtitle(paste0("Reporting Fraction = ", rf)) + theme(legend.position = "none") +
    theme(axis.title.x = element_blank(), plot.margin = margin(10, 10, 10, 10)) +
    ylab("Cumulative Deaths") + scale_x_date(breaks = "4 months", date_labels = "%b %Y")
}

supp_1a <- cowplot::plot_grid(plotlist = pl, ncol = 3) + theme(plot.background = element_rect(fill="white"))
supp_1b <- cowplot::plot_grid(plotlist = pl2, ncol = 3) + theme(plot.background = element_rect(fill="white"))
save_figs("khartoum_model_fits_daily", supp_1a, width = 10, height = 12, root = "analysis/figures/supplementary")
save_figs("khartoum_model_fits_cumulative", supp_1b, width = 10, height = 12, root = "analysis/figures/supplementary")

# -----------------------------------------------
## Chi-square statistical comparison
# -----------------------------------------------

# explore all models looked at
rfs <- c(seq(0.02,0.06, 0.01), 0.07, 0.08, 0.09, 0.1,0.2,0.4, 0.045)
fits <- paste0("analysis/data/derived/model_fits/Khartoum/khartoum_fit_rf_", rfs, ".rds")


res_list <- list(
  readRDS("analysis/data/derived/model_fits/Khartoum/khartoum_fit_rf_0.04.rds") %>% squire::projections(10),
  readRDS("analysis/data/derived/model_fits/Khartoum/khartoum_fit_rf_0.045.rds") %>% squire::projections(10)
)
difff <- function(x) {c(0, diff(x))}

# In Moser et al. they have already adjusted for test sensitivity and waning performance
# so just work out the proportion not in S lagged for seroconversion time
sero_det <- sero_det("iggm", igg_sens = 1, igg_conv = 13.3, igg_scale = Inf)

sero_df_list <- lapply(
  fits,
  function(y){
    res <- readRDS(y)
    res %>%
      squire::projections(10) %>%
  squire::format_output("S", date_0 = max(res$pmcmc_results$inputs$data$date)) %>%
  group_by(date, compartment, replicate) %>%
  mutate(across(y, ~1-(.x/sum(res$parameters$population)))) %>%
  ungroup %>%
  group_by(replicate) %>%
  mutate(across(y, difff)) %>%
  mutate(across(y, lag, 5, 0)) %>%
  mutate(across(y, roll_func, sero_det)) %>%
  group_by(date, compartment) %>%
  summarise(med = median(y, na.rm = TRUE),
            min = quantile(y, 0.025, na.rm = TRUE),
            max = quantile(y, 0.975, na.rm = TRUE),
            var = var(y, na.rm = TRUE))
})

# number for manuscript of best fitting model at mid point of omdurman sero survey
sero_comp <- sero_df_list %>% lapply(filter, date == "2021-03-21")

# chi_sq_stat
chi_sq <- lapply(sero_comp, function(x) {
  (x$med - omd_sero_med)^2 / (x$var + ((omd_sero_high - omd_sero_low)/3.92)^2)
})

#
data.frame("chi-sq" = unlist(chi_sq),
           "p" = 1-pchisq(unlist(chi_sq), df=1),
           "rf" = rfs) %>%
  arrange(rf)


# And for comparison against FETP

1 - pchisq(
  lapply(sero_df_list[4], function(x) {filter(x, date == "2020-06-13")}) %>%
  lapply(function(x) {
    (x$med - 0.183)^2 / (x$var + ((0.209 - 0.160)/3.92)^2)
  }) %>%
  unlist,
  df = 1)
