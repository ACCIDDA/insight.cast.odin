#' Use an odin2 model with incast
#'
#' `ODIN2()` connects a discrete-time [odin2::odin()] model to the ordinary
#' fable interface used by [incast::get_cv()] and
#' [incast::get_fcast()]. The model is fitted with dust2 and monty,
#' then projected from its fitted final state.
#'
#' The odin2 model follows three conventions:
#'
#' - the observed data are declared as `cases <- data()`;
#' - the scalar state to forecast is named `incidence`;
#' - one model-time unit is one reporting interval.
#'
#' The model must also contain a likelihood such as
#' `cases ~ Poisson(incidence + 1e-6)`. See
#' `vignette("incast-odin2")` for a complete SIR example.
#'
#' `ODIN2()` uses dust2's deterministic likelihood and deterministic projection.
#' Forecast uncertainty comes from posterior parameter draws and
#' `observation_model`; process-noise particle filtering is outside this small
#' initial interface.
#'
#' @param formula The series to forecast, normally `observation`.
#' @param generator A system generator returned by [odin2::odin()].
#' @param prior A prior created with [monty::monty_dsl()]. Parameters named in
#'   the prior are estimated. Other parameters use their odin2 defaults or
#'   values supplied in `fixed`. Parameters named by [odin_negbin()] or
#'   [odin_poisson()] must be present in `prior` or `fixed`, because the R
#'   observation function cannot read defaults from the compiled generator.
#' @param fixed Named list of parameters held fixed during fitting.
#' @param observation_model Function that simulates reported cases from
#'   `(incidence, parameters)`. It should match the likelihood in the odin2
#'   model. Use [odin_negbin()] or [odin_poisson()].
#' @param n_steps Number of MCMC steps. The first half is discarded.
#' @param n_draws Number of posterior draws retained for the forecast.
#' @param dt Simulation step as a fraction of one reporting interval. For
#'   example, `0.25` gives four model updates per interval.
#'
#' @return A fable model definition for `incast` or `fabletools`.
#'
#' @examplesIf FALSE
#' library(incast)
#'
#' model <- ODIN2(
#'   observation,
#'   generator = my_sir,
#'   prior = monty::monty_dsl({
#'     beta ~ Exponential(mean = 2)
#'     gamma ~ Exponential(mean = 1)
#'   }),
#'   fixed = list(N = 1e6, I0 = 10),
#'   observation_model = odin_poisson(scale = NULL)
#' )
#'
#' check_data(example_data) |>
#'   get_fcast(models = list(SIR = model))
#'
#' @seealso [odin_sir()], `vignette("incast-odin2")`
#' @export
#' @importFrom fabletools new_model_definition
ODIN2 <- function(
  formula,
  generator,
  prior,
  fixed = list(),
  observation_model,
  n_steps = 2000L,
  n_draws = 500L,
  dt = 0.25
) {
  spec <- validate_odin2_spec(
    generator = generator,
    prior = prior,
    fixed = fixed,
    observation_model = observation_model,
    n_steps = n_steps,
    n_draws = n_draws,
    dt = dt
  )

  fabletools::new_model_definition(model_odin2, {{ formula }}, spec = spec)
}


#' Fit an odin2 model to one time series
#'
#' @keywords internal
#' @noRd
#' @importFrom tsibble measured_vars
train_odin2 <- function(.data, specials, spec, ...) {
  measured <- tsibble::measured_vars(.data)
  if (length(measured) != 1L) {
    stop("`ODIN2()` needs one response series.", call. = FALSE)
  }

  observations <- as.numeric(.data[[measured]])
  n_observed <- sum(!is.na(observations))
  if (n_observed < 5L) {
    stop(
      "`ODIN2()` needs at least 5 observations; got ", n_observed, ".",
      call. = FALSE
    )
  }
  if (any(observations < 0, na.rm = TRUE)) {
    stop("`ODIN2()` observations must be non-negative.", call. = FALSE)
  }

  time <- odin2_times(.data)
  fit <- fit_odin2(spec, observations, time)

  structure(
    c(
      fit,
      list(
        spec = spec,
        end_time = max(time),
        observations = observations
      )
    ),
    class = "model_odin2"
  )
}


#' Convert an index to model time
#'
#' Complete data produce `1, 2, 3, ...`; implicit gaps remain gaps.
#'
#' @keywords internal
#' @noRd
odin2_times <- function(.data) {
  index <- as.numeric(.data[[tsibble::index_var(.data)]])
  if (length(index) < 2L) {
    return(seq_along(index))
  }

  interval <- min(diff(sort(unique(index))))
  as.integer(round((index - min(index)) / interval)) + 1L
}


specials_odin2 <- fabletools::new_specials(
  xreg = function(...) {
    stop("`ODIN2()` does not support formula regressors.", call. = FALSE)
  }
)


model_odin2 <- fabletools::new_model_class(
  "odin2",
  train = train_odin2,
  specials = specials_odin2,
  check = function(.data) {
    if (!tsibble::is_regular(.data)) {
      stop("`ODIN2()` needs a regular time series.", call. = FALSE)
    }
  }
)


# -----------------------------------------------------------------------------
# fable methods
# -----------------------------------------------------------------------------

#' @importFrom fabletools model_sum
#' @export
model_sum.model_odin2 <- function(x) {
  sprintf("ODIN2[%s, %d draws]", x$spec$name, ncol(x$parameters))
}


#' @importFrom fabletools report
#' @export
report.model_odin2 <- function(object, ...) {
  cat("\n--- odin2 model ---\n\n")
  cat(sprintf("  Model        : %s\n", object$spec$name))
  cat(sprintf("  Observations : %d\n", sum(!is.na(object$observations))))
  cat(sprintf("  MCMC steps   : %d\n", object$spec$n_steps))
  cat(sprintf("  Draws kept   : %d\n", ncol(object$parameters)))
  cat(sprintf("  Acceptance   : %.2f\n\n", object$acceptance))
  cat("  Posterior:\n")
  print(round(odin2_posterior_summary(object), 4L))
  invisible(object)
}


#' @importFrom fabletools tidy
#' @export
tidy.model_odin2 <- function(x, ...) {
  summary <- odin2_posterior_summary(x)
  data.frame(
    term = rownames(summary),
    estimate = summary[, "median"],
    lower = summary[, "lower"],
    upper = summary[, "upper"],
    row.names = NULL
  )
}


#' @importFrom fabletools glance
#' @export
glance.model_odin2 <- function(x, ...) {
  data.frame(
    model = x$spec$name,
    n_draws = ncol(x$parameters),
    acceptance = x$acceptance,
    log_posterior = stats::median(x$density),
    n_observations = sum(!is.na(x$observations))
  )
}


#' @importFrom stats fitted
#' @export
fitted.model_odin2 <- function(object, ...) {
  rep(NA_real_, length(object$observations))
}


#' @importFrom stats residuals
#' @export
residuals.model_odin2 <- function(object, ...) {
  rep(NA_real_, length(object$observations))
}


#' @importFrom fabletools forecast
#' @importFrom distributional dist_sample
#' @export
forecast.model_odin2 <- function(object, new_data, specials = NULL, ...) {
  horizon <- NROW(new_data)
  paths <- tryCatch(
    project_odin2(object, horizon),
    error = function(e) {
      stop(
        "The ", object$spec$name, " projection failed: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  distributional::dist_sample(
    lapply(seq_len(horizon), function(i) paths[, i])
  )
}
