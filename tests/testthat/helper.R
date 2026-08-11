# Shared test fixtures.

#' A weekly epidemic curve for tests
#'
#' Simulated from the package's own SIR model with a fixed seed, so the tests
#' exercise the real model rather than a synthetic shape, and are reproducible
#' without an RNG in the test body.
make_epidemic <- function(
  n = 30L,
  population = 1e6,
  beta = 2,
  gamma = 1,
  rho = 0.4,
  seed = 1L
) {
  sys <- dust2::dust_system_create(
    sir,
    list(N = population, beta = beta, gamma = gamma, I0 = 20, rho = rho, phi = 20),
    n_particles = 1L,
    dt = 0.25,
    seed = seed
  )
  dust2::dust_system_set_state_initial(sys)
  y <- dust2::dust_system_simulate(sys, seq_len(n) - 1L)
  round(rho * as.numeric(dust2::dust_unpack_state(sys, y)$incidence))
}


#' A keyed modelling tsibble, as incast's `as_model_ts()` would build it
make_model_ts <- function(locations = "NY", n = 30L, ...) {
  df <- do.call(
    rbind,
    lapply(seq_along(locations), function(i) {
      data.frame(
        location = locations[[i]],
        target_end_date = as.Date("2025-01-04") + 7L * (seq_len(n) - 1L),
        observation = make_epidemic(n = n, seed = i, ...)
      )
    })
  )
  tsibble::as_tsibble(df, index = target_end_date, key = location)
}


#' A short prior over the built-in SIR parameters
test_prior <- function() {
  monty::monty_dsl({
    beta ~ Exponential(mean = 2)
    gamma ~ Exponential(mean = 1)
  })
}

