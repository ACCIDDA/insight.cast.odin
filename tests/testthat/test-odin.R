# Keep fits short while exercising the real compiled model.

test_that("SIR fits and forecasts through fable", {
  ts <- make_model_ts("NY", n = 30L)
  fit <- fabletools::model(
    ts,
    SIR = odin_sir(
      observation,
      population = 1e6,
      n_steps = 200L,
      n_draws = 50L
    )
  )

  model <- fit[["SIR"]][[1L]]$fit
  expect_s3_class(model, "model_odin2")
  expect_equal(ncol(model$parameters), 50L)

  forecast <- dplyr::as_tibble(fabletools::forecast(fit, h = 4L))
  expect_equal(nrow(forecast), 4L)
  expect_s3_class(forecast$observation, "distribution")
  expect_true(all(is.finite(forecast$.mean)))
  expect_true(all(forecast$.mean >= 0))
})


test_that("forecast paths are non-negative count draws", {
  fit <- fabletools::model(
    make_model_ts("NY", n = 25L),
    SIR = odin_sir(
      observation,
      n_steps = 200L,
      n_draws = 40L
    )
  )
  paths <- project_odin2(fit[["SIR"]][[1L]]$fit, horizon = 3L)

  expect_equal(dim(paths), c(40L, 3L))
  expect_true(all(is.finite(paths)))
  expect_true(all(paths >= 0))
  expect_true(all(paths == round(paths)))
})


test_that("fable fits each insight.cast series independently", {
  fit <- fabletools::model(
    make_model_ts(c("NY", "CT"), n = 25L),
    SIR = odin_sir(
      observation,
      n_steps = 200L,
      n_draws = 20L
    )
  )
  forecast <- dplyr::as_tibble(fabletools::forecast(fit, h = 2L))

  expect_equal(nrow(fit), 2L)
  expect_setequal(forecast$location, c("NY", "CT"))
})


test_that("missing observations are allowed", {
  ts <- make_model_ts("NY", n = 30L)
  ts$observation[c(5L, 17L)] <- NA_real_
  fit <- fabletools::model(
    ts,
    SIR = odin_sir(
      observation,
      n_steps = 200L,
      n_draws = 20L
    )
  )

  expect_true(all(is.finite(fit[["SIR"]][[1L]]$fit$density)))
})


test_that("model time preserves implicit gaps", {
  ts <- make_model_ts("NY", n = 20L)
  gappy <- ts[-c(8L, 9L), ]

  expect_equal(odin2_times(ts), 1:20)
  expect_equal(odin2_times(gappy), c(1:7, 10:20))

  fit <- fabletools::model(
    gappy,
    SIR = odin_sir(
      observation,
      n_steps = 200L,
      n_draws = 20L
    )
  )
  expect_equal(fit[["SIR"]][[1L]]$fit$end_time, 20L)
})


test_that("cross-validation origins fit through the same fable path", {
  ts <- make_model_ts("NY", n = 30L)
  origins <- max(ts$target_end_date) - c(21L, 14L, 7L)
  cv_ts <- tidyr::expand_grid(.id = 1:3, dplyr::as_tibble(ts)) |>
    dplyr::filter(target_end_date < origins[.id]) |>
    tsibble::as_tsibble(index = target_end_date, key = c(location, .id))

  fit <- fabletools::model(
    cv_ts,
    SIR = odin_sir(
      observation,
      n_steps = 200L,
      n_draws = 20L
    )
  )

  expect_equal(nrow(fit), 3L)
  expect_equal(nrow(dplyr::as_tibble(fabletools::forecast(fit, h = 2L))), 6L)
})


test_that("a custom odin2 SIR model uses the same adapter", {
  skip_if_not_installed("odin2")
  generator <- odin2::odin(
    {
      initial(S) <- N - I0
      initial(I) <- I0
      initial(incidence, zero_every = 1) <- 0
      update(S) <- S - n_SI
      update(I) <- I + n_SI - n_IR
      update(incidence) <- incidence + n_SI
      n_SI <- Binomial(max(S, 0), 1 - exp(-beta * max(I, 0) / N * dt))
      n_IR <- Binomial(max(I, 0), 1 - exp(-gamma * dt))
      beta <- parameter(2)
      gamma <- parameter(1)
      N <- parameter(1e6)
      I0 <- parameter(10)
      cases <- data()
      cases ~ Poisson(incidence + 1e-6)
    },
    quiet = TRUE
  )

  fit <- fabletools::model(
    make_model_ts("NY", n = 25L, rho = 1),
    SIR = ODIN2(
      observation,
      generator,
      prior = test_prior(),
      fixed = list(N = 1e6, I0 = 10),
      observation_model = odin_poisson(scale = NULL),
      n_steps = 200L,
      n_draws = 20L
    )
  )

  expect_equal(nrow(dplyr::as_tibble(fabletools::forecast(fit, h = 3L))), 3L)
})


test_that("the model runs through insight.cast", {
  skip_if_not_installed("insight.cast")
  data <- dplyr::as_tibble(make_model_ts("NY", n = 25L))
  data$target <- "cases"

  forecast <- insight.cast::check_data(data) |>
    insight.cast::get_fcast(
      models = list(
        SIR = odin_sir(
          observation,
          n_steps = 200L,
          n_draws = 20L
        )
      ),
      h = 2L
    )

  expect_s3_class(forecast, "insightcast_fcast")
  expect_contains(unique(forecast$hub$model_out_tbl$model_id), "SIR")
})


test_that("the model runs through the complete insight.cast pipeline", {
  skip_if_not_installed("insight.cast")
  set.seed(1)

  data <- insight.cast::example_data |>
    dplyr::filter(
      location == "NY",
      target_end_date >= as.Date("2024-10-19"),
      as_of <= as.Date("2024-12-29")
    )

  models <- list(
    SIR = odin_sir(
      observation,
      population = 1e6,
      n_steps = 200L,
      n_draws = 20L
    )
  )

  forecast <- data |>
    insight.cast::check_data() |>
    insight.cast::get_ncast(max_delay = 2L, draws = 50L) |>
    insight.cast::get_cv(h = 2L, n_origins = 2L, step = 1L, models = models) |>
    insight.cast::get_fcast(top_n = 1L)

  expect_s3_class(forecast, "insightcast_fcast")
  expect_contains(unique(forecast$hub$model_out_tbl$model_id), "SIR")
})


test_that("too little or negative data fail clearly", {
  short <- make_model_ts("NY", n = 4L)
  expect_error(
    train_odin2(short, specials = NULL, spec = odin_sir(observation)$args$spec),
    "at least 5"
  )

  negative <- make_model_ts("NY", n = 10L)
  negative$observation[3L] <- -1
  expect_error(
    train_odin2(negative, specials = NULL, spec = odin_sir(observation)$args$spec),
    "non-negative"
  )
})


test_that("standard fable summaries are available", {
  fit <- fabletools::model(
    make_model_ts("NY", n = 25L),
    SIR = odin_sir(
      observation,
      n_steps = 200L,
      n_draws = 20L
    )
  )
  model <- fit[["SIR"]][[1L]]$fit

  expect_match(fabletools::model_sum(model), "^ODIN2\\[sir, 20 draws\\]$")
  expect_output(fabletools::report(model), "odin2 model")
  expect_equal(fabletools::glance(fit)$n_draws, 20L)
  expect_setequal(
    fabletools::tidy(fit)$term,
    c("beta", "gamma", "rho", "I0", "phi")
  )
  expect_length(stats::fitted(model), 25L)
  expect_length(stats::residuals(model), 25L)
})
