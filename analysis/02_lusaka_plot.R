postmortem_fits <- readRDS("analysis/data/derived/seroprevalence/Lusaka/Lusaka_postmortem_lls.rds")
deaths <- readRDS("analysis/data/derived/deaths_time_series/lusaka_postmortem_deaths.rds")
official <- readRDS("analysis/data/derived/deaths_time_series/lusaka_official_deaths.rds")

# suplementary check fits to deaths
for(i in 1:4) {
  postmortem_fits[[i]]$fit$pmcmc_results$inputs$data$deaths <- postmortem_fits[[i]]$fit$pmcmc_results$inputs$data$deaths/7
  }
lusaka_particle_fits <- cowplot::plot_grid(plotlist = lapply(seq_along(postmortem_fits), function(x) {

  title <- match(names(postmortem_fits)[x],
                 c("total_burial_registry_comp",
                   "total_burial_registry_comp_strict",
                   "total_sample_effort",
                   "total_sample_effort_strict"))
  title <- c("Ct < 45, Low", "Ct < 40, Low",
             "Ct < 45, High", "Ct < 40, High")[title]
  ggx <- plot(postmortem_fits[[x]]$fit, particle_fit = TRUE) +
    theme(legend.position = "none") +
    theme(axis.title.x = element_blank()) +
    ggtitle(title)

  ggx <- ggx + geom_step(aes(date, deaths), ggx$layers[[5]]$data)
  ggx$layers[[5]] <- NULL
  ggx
}))
save_figs("lusaka_particle_fits", lusaka_particle_fits, width = 6, height = 6)
root = "analysis/figures"



sero_df <- do.call(rbind,lapply(postmortem_fits, "[[", "sero"))
sero_df$threshold <- "Ct < 45"
sero_df$threshold[grepl("strict", sero_df$death_source)] <- "Ct < 40"
sero_df$sampling <- "Low"
sero_df$sampling[grepl("sample", sero_df$death_source)] <- "High"


sero_fit_gg <- sero_df %>% na.omit %>%
  group_by(date, sampling, threshold) %>%
  summarise(sero_perc = median(sero_perc, na.omit = TRUE)) %>%
  pivot_wider(names_from = sampling, values_from = sero_perc) %>%
  na.omit %>%
  ggplot(aes(date, ymin = Low, ymax = High, color = threshold, fill = threshold)) +
  geom_ribbon(alpha = 0.2) +
  geom_point(aes(x= as.Date("2020-07-15"),y=0.021), color = "black", inherit.aes = FALSE) +
  geom_errorbar(aes(ymin=0.011,ymax=0.031,x=as.Date("2020-07-15"), width=10), color = "black", inherit.aes = FALSE) +
  geom_errorbarh(aes(xmin=as.Date("2020-07-04"),xmax=as.Date("2020-07-27"),y=0.021, height=0),color = "black", inherit.aes = FALSE) +
  geom_point(aes(x= as.Date("2020-07-29"),y=0.106), color="black", inherit.aes = FALSE) +
  annotate("text", label = "Mulenga et al.*", x= as.Date("2020-07-29")+5,y=0.110, hjust= 0) +
  annotate("text", label = "Mulenga et al.", x= as.Date("2020-07-15")+5,y=0.025, hjust= 0) +
  geom_errorbar(aes(ymin=0.073,ymax=0.139,x=as.Date("2020-07-29"), width=10), color="black", inherit.aes = FALSE) +
  geom_errorbarh(aes(xmin=as.Date("2020-07-18"),xmax=as.Date("2020-08-10"),y=0.106, height=0), color="black", inherit.aes = FALSE) +
  ggpubr::theme_pubclean(base_size = 14) +
  theme(axis.line = element_line(), panel.grid.major.x = element_blank(), axis.title.x = element_blank()) +
  scale_color_viridis_d(name = "Ct Threshold:", end = 0.8) +
  scale_fill_viridis_d(name = "Ct Threshold:", end = 0.8) +
  scale_y_continuous(labels = scales::percent) +
  ylab("Seropositivity (%)")

pcr_fit_gg <- sero_df %>% na.omit %>%
  group_by(date, sampling, threshold) %>%
  summarise(pcr_perc = median(pcr_perc, na.omit = TRUE)) %>%
  pivot_wider(names_from = sampling, values_from = pcr_perc) %>%
  na.omit %>%
  ggplot(aes(date, ymin = Low, ymax = High, color = threshold, fill = threshold)) +
  geom_ribbon(alpha = 0.2) +
  geom_point(aes(x= as.Date("2020-07-15"),y=0.076), color = "black", inherit.aes = FALSE) +
  annotate("text", x=as.Date("2020-07-15")+5, y=0.076+0.005, label="Mulenga et al.", hjust =0, color = "black") +
  geom_errorbar(aes(ymin=0.047,ymax=0.106,x=as.Date("2020-07-15"), width=10), color = "black") +
  geom_errorbarh(aes(xmin=as.Date("2020-07-04"),xmax=as.Date("2020-07-27"),y=0.076, height=0), color = "black") +
  ggpubr::theme_pubclean(base_size = 14) +
  theme(axis.line = element_line(), panel.grid.major.x = element_blank(), axis.title.x = element_blank()) +
  scale_color_viridis_d(name = "Ct Threshold:", end = 0.8) +
  scale_fill_viridis_d(name = "Ct Threshold:", end = 0.8) +
  scale_y_continuous(labels = scales::percent) +
  ylab("PCR positivity (%)")

data_gg <- deaths %>%
  select(-(official_covid_deaths:total_sample_effort_mid)) %>%
  pivot_longer(-1) %>%
  mutate(threshold = "Ct < 45") %>%
  mutate(threshold = replace(threshold, grepl("strict", name), "Ct < 40")) %>%
  mutate(sampling = "Low") %>%
  mutate(sampling = replace(sampling, grepl("sample", name), "High")) %>%
  select(-name) %>%
  pivot_wider(names_from = sampling, values_from = value) %>%
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = Low, ymax = High, color = threshold, fill = threshold), linetype = "solid", alpha = 0.2) +
  geom_bar(aes(date, deaths, color = "Govt. Reported Deaths", fill = "Govt. Reported Deaths"), official, stat = "identity", alpha = 0.2) +
  ggpubr::theme_pubclean(base_size = 14) +
  theme(axis.line = element_line(), panel.grid.major.x = element_blank(), axis.title.x = element_blank()) +
  scale_color_manual(name = "Source:", values = c(viridis::viridis(2, end = 0.8), "#000000")) +
  scale_fill_manual(name = "Source:", values = c(viridis::viridis(2, end = 0.8), "#000000")) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.2))) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  ylab("COVID-19 Deaths")

######## ll

ll_df <- do.call(rbind, lapply(postmortem_fits, "[[", "ll"))
ll_df <- ll_df %>% mutate(
  threshold = "Ct < 45",
  threshold = replace(threshold, grepl("strict", death_source), "Ct < 40"),
  sampling = "Low",
  sampling = replace(sampling, grepl("sample", death_source), "High"),
  ll = sero_15 + sero_29 + pcr_15
)

library(scales)
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv,
            log_breaks(base = base),
            domain = c(1e-100, Inf))
}
ll_gg <- ll_df %>% ggplot(aes(fill = threshold, x = sampling, y = -ll)) +
  geom_boxplot(alpha = 0.2) +
  scale_y_continuous(trans=reverselog_trans(base=10),
                     labels=trans_format("identity", function(x) -x)) +
  ggpubr::theme_pubclean(base_size = 14) +
  theme(axis.line = element_line(), panel.grid.major.x = element_blank(),
        legend.key = element_rect(fill = "white")) +
  scale_fill_viridis_d(name = "Ct Threshold:", end = 0.8) +
  ylab("Log Likelihood") +
  xlab("Sampling Effort")



######## rf boxplot lls

rf_fits <- readRDS("analysis/data/derived/seroprevalence/Lusaka/Lusaka_official_lls.rds")

ll_rf_df <- do.call(rbind, lapply(rf_fits, "[[", "ll"))
ll_rf_df <- ll_rf_df %>% mutate(
  ll = sero_15 + sero_29 + pcr_15,
  rf = scales::percent(rf, accuracy = 1)
)

ll_rf_df %>% ggplot(aes(x = rf, y = -ll)) +
  geom_boxplot() +
  scale_y_continuous(trans=reverselog_trans(base=10),
                     labels=trans_format("identity", function(x) -x)) +
  ggpubr::theme_pubclean(base_size = 14) +
  theme(axis.line = element_line(), panel.grid.major.x = element_blank()) +
  ylab("Log Likelihood") +
  xlab("Reporting Fraction")



# real data


rf_fit <- readRDS("analysis/data/derived/seroprevalence/Lusaka/Lusaka_official_rf_est_ll.rds")

# observed data
obs <- data.frame(
  "date" = as.Date(c("2020-07-11","2020-07-11","2020-07-11", "2020-07-25")),
  "date_min"  = as.Date(c("2020-07-03", "2020-07-03", "2020-07-03", "2020-07-17")),
  "date_max"  = as.Date(c("2020-07-19", "2020-07-19", "2020-07-19", "2020-08-02")),
  "ymin"  = c((2.1/(9.7))*9.1, (7.6/(9.7))*9.1, 9.1, 9.1)/100,
  "ymax"  = c((2.1/(9.7))*9.1, (7.6/(9.7))*9.1, 9.1, 9.1)/100,
  "y"  = c((2.1/(9.7))*9.1, (7.6/(9.7))*9.1, 9.1, 9.1)/100,
  "col" = c("Serology", "PCR", "Serology+PCR", "Serology")
)

# plot the observed real data comparison

dat <- rf_fit$fit$pmcmc_results$inputs$data
lusaka_e <- squire::format_output(rf_fit$fit, "deaths",date_0 = max(dat$date)) %>%
  na.omit() %>%
  group_by(date) %>%
  summarise(med = median(y),
            min = quantile(y, 0.025),
            max = quantile(y, 0.975)) %>%
  left_join(dat %>%
              mutate(deaths2 = round(deaths/mean(rf_fit$fit$replicate_parameters$rf))) %>%
              complete(date = seq.Date(min(date), max(date), 1)) %>%
              mutate(deaths = replace_na(deaths, 0)) %>%
              mutate(deaths2 = replace_na(deaths2, 0)) %>%
              rename("Govt. Reported Deaths" = deaths) %>%
              rename("Govt. Reported Deaths Scaled" = deaths2)
            ) %>%
  ggplot(aes(date, med, ymin = min, ymax = max, fill = "Model Fit")) +
  geom_ribbon(alpha = 0.4) +
  geom_line(color = viridis::viridis(1, end = 0.5, option = "C")) +
  geom_bar(aes(date, `Govt. Reported Deaths`,  fill = "Govt. Reported Deaths"),
           stat = "identity", alpha = 0.8, size = 0) +
  geom_point(aes(date, `Govt. Reported Deaths Scaled`,
                 color = "Govt. Reported Deaths \nScaled by \nReporting Fraction")) +
  ggpubr::theme_pubclean(base_size = 14) +
  theme(axis.line = element_line(), panel.grid.major.x = element_blank(), axis.title.x = element_blank()) +
  scale_color_manual(name = "Source:", values = c("#000000")) +
  scale_fill_manual(name = "", values = c(viridis::viridis(2, end = 0.5, option = "C", direction = -1))) +
  guides(fill = guide_legend(override.aes = list(alpha = c(0.5), shape = NA), nrow = 2)) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linetype = NA), keywidth = 0, keyheight = 0)) +
  #guides(fill = guide_legend(override.aes = list(alpha = c(0,1,1)))) +
  #guides(color = guide_legend(override.aes = list(alpha = c(1,1,1)))) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  ylab("COVID-19 Deaths")


# our sero comparison

lusaka_g <- rf_fit$sero %>% group_by(date) %>%
  summarise(Serology_med = median(sero_perc, na.omit = TRUE),
            Serology_min = quantile(sero_perc, 0.025, na.omit = TRUE),
            Serology_max = quantile(sero_perc, 0.975, na.omit = TRUE),
            PCR_med = median(pcr_perc, na.omit = TRUE),
            PCR_min = quantile(pcr_perc, 0.025, na.omit = TRUE),
            PCR_max = quantile(pcr_perc, 0.975, na.omit = TRUE),
            `Serology+PCR_med` = median(combined_perc, na.omit = TRUE),
            `Serology+PCR_min` = quantile(combined_perc, 0.025, na.omit = TRUE),
            `Serology+PCR_max` = quantile(combined_perc, 0.975, na.omit = TRUE)) %>%
  na.omit %>%
  pivot_longer(Serology_med:`Serology+PCR_max`, names_pattern = "(.*)_(.*)", names_to = c("name", "stat")) %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  ggplot(aes(date, y = med, ymin = min, ymax = max, color = name, fill = name)) +
  geom_ribbon(alpha = 0.2) +
  geom_line(lwd = 1.5) +
  geom_point(aes(date, y, color = col), data = obs[3:4,], inherit.aes = FALSE) +
  geom_errorbar(aes(ymin=ymin,ymax=ymax,x=date, width=10, color = col), data = obs[3:4,], inherit.aes = FALSE) +
  geom_errorbarh(aes(y=y,xmin=date_min,xmax=date_max, width=2, height = 0.01, color = col), data = obs[3:4,], inherit.aes = FALSE) +
  ggpubr::theme_pubclean(base_size = 14) +
  theme(axis.line = element_line(), panel.grid.major.x = element_blank(), axis.title.x = element_blank()) +
  scale_color_viridis_d(name = "Test:", end = 0.8) +
  scale_fill_viridis_d(name = "Test:", end = 0.8) +
  scale_y_continuous(labels = scales::percent) +
  ylab("Test Positivity (%)")

# likelihood

mc <- rbind(rf_fit$fit$pmcmc_results$chains$chain1$results[5000:10000,])


lusaka_f <- mc[,c("rf", "log_likelihood")] %>%
  unique %>%
ggplot(aes(rf, log_likelihood)) +
  geom_density_2d(contour_var = "ndensity", bins = 5) +
  geom_point(alpha = 0.4) +
  ggpubr::theme_pubclean(base_size = 14) +
  theme(axis.line = element_line(), panel.grid.major.x = element_blank()) +
  ylab("Log Likelihood") +
  scale_x_continuous(labels = scales::percent) +
  ylim(c(-264, -248)) +
  xlab("Reporting Fraction %")


## bring it all together


lusaka_abcd <- cowplot::plot_grid(
  cowplot::plot_grid(
    data_gg,
    pcr_fit_gg + theme(legend.position = "none"),
    sero_fit_gg + theme(legend.position = "none"),
    labels = "auto", align = "v", ncol = 1)
  ,
  cowplot::plot_grid(
    ll_gg + theme(legend.position = "none"),
    ncol = 1, labels = "d"),
  ncol = 1, rel_heights = c(1, 1/3)
)

lusaka_efg <- cowplot::plot_grid(
    lusaka_e,
    lusaka_f,
    lusaka_g,
    labels = c("e", "f", "g"),
    ncol = 1
)


lusaka_all <- cowplot::plot_grid(lusaka_abcd, ggplot(), lusaka_efg, rel_widths = c(1,0.03,1), ncol = 3 )


save_figs("lusaka_main", lusaka_all, root = "analysis/figures", height = 16, width = 14)
