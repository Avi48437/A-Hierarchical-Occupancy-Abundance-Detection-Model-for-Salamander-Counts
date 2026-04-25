library(DHARMa)

simulate_spdet_model <- function(par, pair_data, S) {
  theta <- unpack_spdet_parameters(par, S)
  sim_list <- vector("list", length(pair_data))

  for (k in seq_along(pair_data)) {
    df <- pair_data[[k]]
    s <- as.integer(df$spp[1])

    psi <- plogis(theta$gamma0[s] + theta$gamma1 * df$mined[1])
    mu <- exp(theta$beta0[s] + theta$beta1[s] * df$mined[1] + theta$beta2 * df$cover[1])

    p <- plogis(
      theta$alpha0[s] +
        theta$alpha1[s] * df$DOP +
        theta$alpha2[s] * df$Wtemp +
        theta$alpha3[s] * df$DOY +
        theta$alpha4[s] * df$Wtemp * df$DOY
    )

    Z <- rbinom(1, 1, psi)
    N <- 0

    if (Z == 1) {
      repeat {
        N <- rnbinom(1, size = theta$phi, mu = mu)
        if (N > 0) break
      }
    }

    sim_list[[k]] <- rbinom(length(p), size = N, prob = p)
  }

  sim_list
}

run_spdet_diagnostics <- function(fit, pair_data, S, nsim = 500, seed = 123) {
  set.seed(seed)

  sim_mat <- replicate(
    nsim,
    unlist(simulate_spdet_model(fit$par, pair_data, S))
  )

  obs_vec <- unlist(lapply(pair_data, function(df) df$count))

  res <- createDHARMa(
    simulatedResponse = sim_mat,
    observedResponse = obs_vec,
    integerResponse = TRUE
  )

  list(
    residuals = res,
    uniformity = testUniformity(res),
    dispersion = testDispersion(res),
    zero_inflation = testZeroInflation(res),
    outliers = testOutliers(res)
  )
}

spdet_diagnostics <- run_spdet_diagnostics(fit_nb_spdet, pair_data, S)
plot(spdet_diagnostics$residuals)
spdet_diagnostics[setdiff(names(spdet_diagnostics), "residuals")]
