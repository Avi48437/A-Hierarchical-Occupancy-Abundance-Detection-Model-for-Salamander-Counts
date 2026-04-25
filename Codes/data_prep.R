library(dplyr)
library(glmmTMB)

data("Salamanders", package = "glmmTMB")
Data <- Salamanders

prepare_salamander_data <- function(Data) {
  Data %>%
    mutate(
      site = factor(site),
      spp = factor(spp),
      mined_fac = factor(mined, levels = c("no", "yes")),
      mined = as.integer(mined == "yes"),
      cover = as.numeric(scale(cover)),
      DOP = as.numeric(scale(DOP)),
      Wtemp = as.numeric(scale(Wtemp)),
      DOY = as.numeric(scale(DOY))
    ) %>%
    arrange(spp, site, sample)
}

make_pair_data <- function(Data_fit) {
  Data_fit %>%
    group_by(spp, site) %>%
    group_split()
}

make_empirical_summaries <- function(Data_fit) {
  occ_summary <- Data_fit %>%
    group_by(site, spp, mined) %>%
    summarise(Z_si = as.integer(sum(count) > 0), .groups = "drop") %>%
    group_by(spp, mined) %>%
    summarise(psi_hat = mean(Z_si), .groups = "drop")

  abund_summary <- Data_fit %>%
    group_by(site, spp, mined) %>%
    summarise(total_count = sum(count), .groups = "drop") %>%
    filter(total_count > 0) %>%
    group_by(spp, mined) %>%
    summarise(lambda_hat = mean(total_count), .groups = "drop")

  list(occ_summary = occ_summary, abund_summary = abund_summary)
}

Data_fit <- prepare_salamander_data(Data)
pair_data <- make_pair_data(Data_fit)
S <- nlevels(Data_fit$spp)

emp_summary <- make_empirical_summaries(Data_fit)
occ_summary <- emp_summary$occ_summary
abund_summary <- emp_summary$abund_summary
