library(dplyr)
library(glmmTMB)

fit_zinb_models <- function(Data) {
  Data_glmm <- Data %>%
    mutate(
      site = factor(site),
      spp = factor(spp),
      mined = factor(mined)
    )

  models <- list(
    zinb_main = glmmTMB(
      count ~ spp + mined + cover + DOP + Wtemp + DOY + (1 | site),
      ziformula = ~ spp + mined,
      family = nbinom2,
      data = Data_glmm
    ),
    zinb_spp_mined = glmmTMB(
      count ~ spp + mined + spp:mined + cover + DOP + Wtemp + DOY + (1 | site),
      ziformula = ~ spp + mined,
      family = nbinom2,
      data = Data_glmm
    ),
    zinb_spp_random = glmmTMB(
      count ~ mined + cover + DOP + Wtemp + DOY + (1 | site) + (1 | spp),
      ziformula = ~ mined + (1 | spp),
      family = nbinom2,
      data = Data_glmm
    )
  )

  comparison <- data.frame(
    Model = names(models),
    AIC = sapply(models, AIC),
    BIC = sapply(models, BIC),
    logLik = sapply(models, function(x) as.numeric(logLik(x))),
    row.names = NULL
  )

  list(models = models, comparison = comparison)
}

zinb_results <- fit_zinb_models(Data)
print(zinb_results$comparison)
