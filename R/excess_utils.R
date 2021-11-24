# draw from these negative binomials to create an excess mortality trend
create_realisation <- function(rep = 1, deaths_df, baseline) {
  deaths_df$baseline <- vapply(deaths_df$month, function(x) {
    rnbinom(1, size = baseline$nb[[x]]$estimate[1], mu = baseline$nb[[x]]$estimate[2])
  }, numeric(1))

  deaths_df$deaths <-  round((deaths_df$sum - deaths_df$baseline))
  deaths_df$rep <- rep
  return(deaths_df)
}
