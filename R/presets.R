#' Ready-to-use SIR model
#'
#' `odin_sir()` is the shortest way to add a mechanistic model to an
#' incast forecast. It uses a precompiled susceptible-infectious-recovered
#' model with a negative binomial observation model.
#'
#' The model estimates transmission (`beta`), recovery (`gamma`), reporting
#' (`rho`), initial infections (`I0`) and overdispersion (`phi`). Model rates
#' are measured per reporting interval.
#'
#' @param formula The series to forecast, normally `observation`.
#' @param population Population size used by every fitted series.
#' @param prior Optional replacement prior created with
#'   [monty::monty_dsl()]. It must include `rho` and `phi`, which are also used
#'   to simulate future observations. Other omitted parameters use their odin2
#'   defaults. Use [ODIN2()] directly if `rho` or `phi` should be fixed.
#' @param ... Computational settings passed to [ODIN2()]: `n_steps`,
#'   `n_draws` or `dt`.
#'
#' @return A fable model definition.
#'
#' @examplesIf FALSE
#' library(incast)
#'
#' models <- c(
#'   default_models(),
#'   list(SIR = odin_sir(observation, population = 2e7))
#' )
#'
#' check_data(example_data) |>
#'   get_fcast(models = models)
#'
#' @seealso [ODIN2()], [odin_seir()]
#' @export
odin_sir <- function(formula, population = 1e6, prior = NULL, ...) {
  ODIN2(
    {{ formula }},
    generator = sir,
    prior = prior %||% sir_prior(),
    fixed = list(N = population),
    observation_model = odin_negbin(),
    ...
  )
}


#' Ready-to-use SEIR model
#'
#' This is the same interface as [odin_sir()], with an exposed compartment and
#' an additional rate `sigma` for progression from exposed to infectious.
#'
#' @inheritParams odin_sir
#' @return A fable model definition.
#'
#' @seealso [ODIN2()], [odin_sir()]
#' @export
odin_seir <- function(formula, population = 1e6, prior = NULL, ...) {
  ODIN2(
    {{ formula }},
    generator = seir,
    prior = prior %||% seir_prior(),
    fixed = list(N = population),
    observation_model = odin_negbin(),
    ...
  )
}


#' @keywords internal
#' @noRd
sir_prior <- function() {
  monty::monty_dsl({
    beta ~ Exponential(mean = 2)
    gamma ~ Exponential(mean = 1)
    rho ~ Beta(a = 1, b = 1)
    I0 ~ Exponential(mean = 100)
    phi ~ Exponential(mean = 10)
  })
}


#' @keywords internal
#' @noRd
seir_prior <- function() {
  monty::monty_dsl({
    beta ~ Exponential(mean = 2)
    sigma ~ Exponential(mean = 1)
    gamma ~ Exponential(mean = 1)
    rho ~ Beta(a = 1, b = 1)
    I0 ~ Exponential(mean = 100)
    phi ~ Exponential(mean = 10)
  })
}


#' Observation models for odin2 forecasts
#'
#' These helpers turn projected `incidence` into reported cases. The choice
#' should match the likelihood in the odin2 model.
#'
#' @param size Negative binomial size parameter: either its name in the odin2
#'   model or a fixed number.
#' @param scale Reporting parameter: either its name, a fixed number, or `NULL`
#'   for no scaling.
#'
#' @return A function suitable for the `observation_model` argument of
#'   [ODIN2()].
#'
#' @examples
#' simulate_cases <- odin_negbin(size = "phi", scale = "rho")
#' simulate_cases(c(10, 20), list(phi = 10, rho = 0.5))
#'
#' @seealso [ODIN2()]
#' @export
#' @importFrom stats rnbinom rpois
odin_negbin <- function(size = "phi", scale = "rho") {
  model <- function(incidence, parameters) {
    mean <- scale_incidence(incidence, scale, parameters)
    stats::rnbinom(
      length(mean),
      size = parameter_value(size, parameters),
      mu = mean
    )
  }
  attr(model, "odin_parameters") <- observation_parameters(size, scale)
  model
}


#' @rdname odin_negbin
#' @export
odin_poisson <- function(scale = "rho") {
  model <- function(incidence, parameters) {
    mean <- scale_incidence(incidence, scale, parameters)
    stats::rpois(length(mean), mean)
  }
  attr(model, "odin_parameters") <- observation_parameters(scale)
  model
}


#' @keywords internal
#' @noRd
observation_parameters <- function(...) {
  values <- list(...)
  unique(unlist(values[vapply(values, is.character, logical(1L))]))
}


#' @keywords internal
#' @noRd
scale_incidence <- function(incidence, scale, parameters) {
  incidence <- pmax(as.numeric(incidence), 0)
  if (is.null(scale)) {
    return(incidence)
  }
  pmax(incidence * parameter_value(scale, parameters), 0)
}


#' @keywords internal
#' @noRd
parameter_value <- function(value, parameters) {
  if (!is.character(value)) {
    return(value)
  }
  if (!value %in% names(parameters)) {
    stop(
      "Observation model parameter `", value, "` is missing.",
      call. = FALSE
    )
  }
  parameters[[value]]
}
