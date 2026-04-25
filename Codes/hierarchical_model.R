library(dplyr)

unpack_spdet_parameters <- function(par, S) {
  list(
    gamma0 = par[1:S],
    gamma1 = par[S + 1],
    beta0 = par[(S + 2):(2 * S + 1)],
    beta1 = par[(2 * S + 2):(3 * S + 1)],
    beta2 = par[3 * S + 2],
    alpha0 = par[(3 * S + 3):(4 * S + 2)],
    alpha1 = par[(4 * S + 3):(5 * S + 2)],
    alpha2 = par[(5 * S + 3):(6 * S + 2)],
    alpha3 = par[(6 * S + 3):(7 * S + 2)],
    alpha4 = par[(7 * S + 3):(8 * S + 2)],
    phi = exp(par[8 * S + 3])
  )
}

nll_nb_spdet <- function(par, pair_data, S, K_extra = 80) {
  theta <- unpack_spdet_parameters(par, S)
  loglik <- 0

  for (df in pair_data) {
    s <- as.integer(df$spp[1])
    y <- df$count

    psi <- plogis(theta$gamma0[s] + theta$gamma1 * df$mined[1])
    mu <- exp(theta$beta0[s] + theta$beta1[s] * df$mined[1] + theta$beta2 * df$cover[1])

    p <- plogis(
      theta$alpha0[s] +
        theta$alpha1[s] * df$DOP +
        theta$alpha2[s] * df$Wtemp +
        theta$alpha3[s] * df$DOY +
        theta$alpha4[s] * df$Wtemp * df$DOY
    )

    y_max <- max(y)
    n_vals <- seq.int(max(1, y_max), y_max + K_extra)

    log_nb <- dnbinom(n_vals, size = theta$phi, mu = mu, log = TRUE)
    log_nb_trunc <- log_nb - log1p(-dnbinom(0, size = theta$phi, mu = mu))

    log_cond_y <- vapply(
      n_vals,
      function(n) sum(dbinom(y, size = n, prob = p, log = TRUE)),
      numeric(1)
    )

    log_terms <- log_nb_trunc + log_cond_y
    m <- max(log_terms)
    if (!is.finite(m)) return(1e12)

    log_sum <- m + log(sum(exp(log_terms - m)))

    if (all(y == 0)) {
      L_si <- (1 - psi) + psi * exp(log_sum)
      if (!is.finite(L_si) || L_si <= 0) return(1e12)
      loglik <- loglik + log(L_si)
    } else {
      if (!is.finite(psi) || psi <= 0) return(1e12)
      loglik <- loglik + log(psi) + log_sum
    }
  }

  if (!is.finite(loglik)) return(1e12)
  -loglik
}

make_initial_values <- function(occ_summary, abund_summary, S) {
  eps <- 1e-4

  occ_no <- occ_summary %>% filter(mined == 0) %>% arrange(spp)
  occ_yes <- occ_summary %>% filter(mined == 1) %>% arrange(spp)

  abund_no <- abund_summary %>% filter(mined == 0) %>% arrange(spp)
  abund_yes <- abund_summary %>% filter(mined == 1) %>% arrange(spp)

  gamma0_init <- qlogis(pmin(pmax(occ_no$psi_hat, eps), 1 - eps))
  gamma1_init <- mean(
    qlogis(pmin(pmax(occ_yes$psi_hat, eps), 1 - eps)) - gamma0_init
  )

  beta0_init <- log(pmax(abund_no$lambda_hat, eps))
  beta1_init <- log(pmax(abund_yes$lambda_hat, eps)) - beta0_init

  c(
    gamma0_init,
    gamma1_init,
    beta0_init,
    beta1_init,
    0,
    rep(qlogis(0.4), S),
    rep(0, 4 * S),
    log(5)
  )
}

fit_spdet_model <- function(pair_data, S, occ_summary, abund_summary, K_extra = 80) {
  par0 <- make_initial_values(occ_summary, abund_summary, S)

  nlminb(
    start = par0,
    objective = nll_nb_spdet,
    pair_data = pair_data,
    S = S,
    K_extra = K_extra,
    control = list(eval.max = 1000, iter.max = 1000)
  )
}

summarise_spdet_fit <- function(fit, Data_fit, pair_data) {
  spp <- levels(Data_fit$spp)

  param_names <- c(
    paste0("gamma0_", spp),
    "gamma1",
    paste0("beta0_", spp),
    paste0("beta1_", spp),
    "beta2",
    paste0("alpha0_", spp),
    paste0("alpha1_DOP_", spp),
    paste0("alpha2_Wtemp_", spp),
    paste0("alpha3_DOY_", spp),
    paste0("alpha4_Wtemp_DOY_", spp),
    "log_phi"
  )

  logLik <- -fit$objective
  k <- length(fit$par)
  n <- length(pair_data)

  list(
    estimates = data.frame(parameter = param_names, estimate = fit$par),
    fit_stats = data.frame(
      Model = "Hierarchical occupancy-abundance-detection model",
      AIC = -2 * logLik + 2 * k,
      BIC = -2 * logLik + log(n) * k,
      logLik = logLik,
      k = k,
      n = n
    )
  )
}

fit_nb_spdet <- fit_spdet_model(pair_data, S, occ_summary, abund_summary)
spdet_summary <- summarise_spdet_fit(fit_nb_spdet, Data_fit, pair_data)

print(spdet_summary$estimates)
print(spdet_summary$fit_stats)
