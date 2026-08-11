# Validation for the small public ODIN2() contract.


#' Validate an odin2 model specification
#'
#' @keywords internal
#' @noRd
validate_odin2_spec <- function(
  generator,
  prior,
  fixed,
  observation_model,
  n_steps,
  n_draws,
  dt
) {
  if (!inherits(generator, "dust_system_generator")) {
    stop(
      "`generator` must be returned by `odin2::odin()`.",
      call. = FALSE
    )
  }

  properties <- attr(generator, "properties")
  if (!identical(properties$time_type, "discrete")) {
    stop("`ODIN2()` currently supports discrete-time models only.", call. = FALSE)
  }
  if (!isTRUE(properties$has_compare)) {
    stop(
      "The odin2 model needs a likelihood, for example ",
      "`cases <- data(); cases ~ Poisson(incidence + 1e-6)`.",
      call. = FALSE
    )
  }

  if (!inherits(prior, "monty_model")) {
    stop("`prior` must be created with `monty::monty_dsl()`.", call. = FALSE)
  }
  estimated <- prior$parameters
  if (length(estimated) == 0L) {
    stop("`prior` must declare at least one parameter.", call. = FALSE)
  }
  if (!isTRUE(prior$properties$has_direct_sample)) {
    stop("`prior` must support direct sampling.", call. = FALSE)
  }

  validate_fixed(fixed, estimated)
  validate_parameters(generator, estimated, fixed)

  if (!is.function(observation_model)) {
    stop(
      "`observation_model` must be a function. Use `odin_negbin()` or ",
      "`odin_poisson()`.",
      call. = FALSE
    )
  }
  validate_observation_parameters(observation_model, estimated, fixed)

  n_steps <- validate_count(n_steps, "n_steps", 20L)
  n_draws <- validate_count(n_draws, "n_draws", 2L)
  if (n_draws > n_steps %/% 2L) {
    stop("`n_draws` cannot exceed half of `n_steps`.", call. = FALSE)
  }
  if (!is.numeric(dt) || length(dt) != 1L || is.na(dt) || dt <= 0 || dt > 1) {
    stop("`dt` must be one number in (0, 1].", call. = FALSE)
  }
  if (abs(1 / dt - round(1 / dt)) > 1e-8) {
    stop("`dt` must be 1 divided by an integer, such as 1, 0.5 or 0.25.", call. = FALSE)
  }

  validate_incidence_state(generator, fixed, dt)

  list(
    name = attr(generator, "name") %||% "custom",
    generator = generator,
    prior = prior,
    estimated = estimated,
    fixed = fixed,
    observation_model = observation_model,
    n_steps = n_steps,
    n_draws = n_draws,
    dt = dt,
    vcv = prior_vcv(prior, estimated)
  )
}


#' Check named parameters used by the built-in observation helpers
#'
#' odin2 knows the defaults declared in a generator, but a separate R
#' observation function does not. Requiring these few values here prevents a
#' surprising failure only after an otherwise successful model fit.
#'
#' @keywords internal
#' @noRd
validate_observation_parameters <- function(observation_model, estimated, fixed) {
  required <- attr(observation_model, "odin_parameters")
  if (is.null(required)) {
    return(invisible(NULL))
  }

  missing <- setdiff(required, c(estimated, names(fixed)))
  if (length(missing)) {
    stop(
      "Observation model parameter",
      if (length(missing) > 1L) "s `" else " `",
      paste(missing, collapse = "`, `"),
      "` must appear in `prior` or `fixed`.",
      call. = FALSE
    )
  }
  invisible(NULL)
}


#' Validate fixed parameters
#'
#' @keywords internal
#' @noRd
validate_fixed <- function(fixed, estimated) {
  if (!is.list(fixed)) {
    stop("`fixed` must be a named list.", call. = FALSE)
  }
  if (length(fixed) == 0L) {
    return(invisible(NULL))
  }
  if (is.null(names(fixed)) || any(!nzchar(names(fixed)))) {
    stop("`fixed` must be a named list.", call. = FALSE)
  }

  duplicated <- intersect(names(fixed), estimated)
  if (length(duplicated)) {
    stop(
      "Parameter", if (length(duplicated) > 1L) "s " else " ",
      paste(duplicated, collapse = ", "),
      " cannot appear in both `prior` and `fixed`.",
      call. = FALSE
    )
  }
  invisible(NULL)
}


#' Validate parameters against the generator
#'
#' @keywords internal
#' @noRd
validate_parameters <- function(generator, estimated, fixed) {
  table <- attr(generator, "parameters")
  supplied <- c(estimated, names(fixed))
  unknown <- setdiff(supplied, table$name)
  if (length(unknown)) {
    stop(
      "Unknown odin2 parameter", if (length(unknown) > 1L) "s: " else ": ",
      paste(unknown, collapse = ", "), ".",
      call. = FALSE
    )
  }

  missing <- setdiff(table$name[table$required], supplied)
  if (length(missing)) {
    stop(
      "Required odin2 parameter", if (length(missing) > 1L) "s are" else " is",
      " missing: ", paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }
  invisible(NULL)
}


#' Check the required scalar incidence state when possible
#'
#' @keywords internal
#' @noRd
validate_incidence_state <- function(generator, fixed, dt) {
  states <- tryCatch(
    {
      system <- dust2::dust_system_create(generator, fixed, dt = dt)
      dust2::dust_unpack_index(system)
    },
    error = function(e) NULL
  )
  if (is.null(states)) {
    return(invisible(NULL))
  }
  if (is.null(states$incidence)) {
    stop("The odin2 model needs a state named `incidence`.", call. = FALSE)
  }
  if (length(states$incidence) != 1L) {
    stop("The odin2 `incidence` state must be scalar.", call. = FALSE)
  }
  invisible(NULL)
}


#' Initial proposal covariance estimated from the prior
#'
#' @keywords internal
#' @noRd
prior_vcv <- function(prior, estimated) {
  rng <- monty::monty_rng_create()
  draws <- vapply(
    seq_len(100L),
    function(i) monty::monty_model_direct_sample(prior, rng),
    numeric(length(estimated))
  )
  scale <- apply(matrix(draws, nrow = length(estimated)), 1L, stats::sd)
  scale[!is.finite(scale) | scale <= 0] <- 1
  diag((0.1 * scale)^2, nrow = length(estimated))
}


#' Validate an integer control
#'
#' @keywords internal
#' @noRd
validate_count <- function(x, name, minimum) {
  if (
    !is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != round(x) || x < minimum
  ) {
    stop("`", name, "` must be a single integer >= ", minimum, ".", call. = FALSE)
  }
  as.integer(x)
}


#' @keywords internal
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x
