# The odin2 -> dust2 -> monty pipeline used by the fable adapter.


#' Fit one series
#'
#' @keywords internal
#' @noRd
fit_odin2 <- function(spec, observations, time) {
  data <- data.frame(time = as.integer(time), cases = observations)
  likelihood_engine <- dust2::dust_unfilter_create(
    spec$generator,
    time_start = 0,
    data = data,
    dt = spec$dt
  )

  packer <- monty::monty_packer(spec$estimated, fixed = spec$fixed)
  likelihood <- dust2::dust_likelihood_monty(
    likelihood_engine,
    packer,
    save_state = TRUE
  )
  posterior <- spec$prior + likelihood
  initial <- odin2_initial(posterior, spec)

  # A misspelled data variable is otherwise silently ignored by dust2.
  log_likelihood <- dust2::dust_likelihood_run(
    likelihood_engine,
    packer$unpack(initial)
  )
  if (identical(as.numeric(log_likelihood), 0) && any(observations > 0, na.rm = TRUE)) {
    stop(
      "The odin2 likelihood did not use the observations. ",
      "Declare them as `cases <- data()` in the odin2 model.",
      call. = FALSE
    )
  }

  sampler <- monty::monty_sampler_adaptive(initial_vcv = spec$vcv)
  samples <- monty::monty_sample(
    posterior,
    sampler,
    n_steps = spec$n_steps,
    initial = initial,
    burnin = spec$n_steps %/% 2L
  )
  samples <- monty::monty_flatten_chains(samples)
  keep <- thin_to(ncol(samples$pars), spec$n_draws)

  list(
    parameters = samples$pars[, keep, drop = FALSE],
    state = samples$observations$state[, keep, drop = FALSE],
    density = as.numeric(samples$density)[keep],
    acceptance = mean(diff(as.numeric(samples$density)) != 0)
  )
}


#' Choose a valid MCMC starting point
#'
#' @keywords internal
#' @noRd
odin2_initial <- function(posterior, spec, attempts = 10L) {
  rng <- monty::monty_rng_create()
  candidates <- vapply(
    seq_len(attempts),
    function(i) monty::monty_model_direct_sample(spec$prior, rng),
    numeric(length(spec$estimated))
  )
  density <- apply(candidates, 2L, function(candidate) {
    tryCatch(
      monty::monty_model_density(posterior, candidate),
      error = function(e) -Inf
    )
  })

  if (all(!is.finite(density))) {
    stop(
      "Could not find a valid starting point from the prior. ",
      "Check `prior` and `fixed`.",
      call. = FALSE
    )
  }

  candidates[, which.max(density)]
}


#' Continue posterior states into the forecast period
#'
#' @keywords internal
#' @noRd
project_odin2 <- function(object, horizon) {
  spec <- object$spec
  packer <- monty::monty_packer(spec$estimated, fixed = spec$fixed)
  n_draws <- ncol(object$parameters)
  parameters <- lapply(
    seq_len(n_draws),
    function(i) packer$unpack(object$parameters[, i])
  )

  system <- dust2::dust_system_create(
    spec$generator,
    parameters,
    n_particles = 1L,
    n_groups = n_draws,
    time = object$end_time,
    dt = spec$dt,
    deterministic = TRUE,
    preserve_particle_dimension = TRUE
  )

  state <- object$state
  state[state < 0] <- 0
  dust2::dust_system_set_state(
    system,
    array(state, c(nrow(state), 1L, n_draws))
  )

  incidence_index <- dust2::dust_unpack_index(system)$incidence
  if (length(incidence_index) != 1L) {
    stop("`incidence` must be a scalar state.", call. = FALSE)
  }

  simulation <- dust2::dust_system_simulate(
    system,
    object$end_time + seq_len(horizon)
  )
  incidence <- dust2::dust_unpack_state(system, simulation)$incidence
  paths <- matrix(as.numeric(incidence), nrow = n_draws, ncol = horizon)

  for (i in seq_len(n_draws)) {
    paths[i, ] <- spec$observation_model(paths[i, ], parameters[[i]])
  }
  paths
}


#' Select evenly spaced posterior draws
#'
#' @keywords internal
#' @noRd
thin_to <- function(n, n_draws) {
  unique(round(seq(1, n, length.out = n_draws)))
}


#' Summarise posterior parameters
#'
#' @keywords internal
#' @noRd
odin2_posterior_summary <- function(object) {
  result <- t(apply(object$parameters, 1L, function(x) {
    c(
      median = stats::median(x),
      lower = unname(stats::quantile(x, 0.025)),
      upper = unname(stats::quantile(x, 0.975))
    )
  }))
  colnames(result) <- c("median", "lower", "upper")
  result
}
