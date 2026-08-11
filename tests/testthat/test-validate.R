test_that("ODIN2() and the presets build fable definitions", {
  model <- ODIN2(
    observation,
    sir,
    prior = test_prior(),
    fixed = list(N = 1e6),
    observation_model = odin_poisson(scale = NULL)
  )

  expect_s3_class(model, "mdl_defn")
  expect_s3_class(odin_sir(observation), "mdl_defn")
  expect_s3_class(odin_seir(observation), "mdl_defn")
})


test_that("generator and prior errors are direct", {
  expect_error(
    ODIN2(
      observation,
      "sir",
      prior = test_prior(),
      observation_model = odin_poisson()
    ),
    "odin2::odin"
  )
  expect_error(
    ODIN2(
      observation,
      sir,
      prior = list(beta = 1),
      observation_model = odin_poisson()
    ),
    "monty::monty_dsl"
  )
})


test_that("a generator needs an odin2 likelihood", {
  skip_if_not_installed("odin2")
  generator <- odin2::odin({
    initial(incidence) <- 0
    update(incidence) <- incidence + 1
  }, quiet = TRUE)

  expect_error(
    ODIN2(
      observation,
      generator,
      prior = test_prior(),
      observation_model = odin_poisson()
    ),
    "needs a likelihood"
  )
})


test_that("parameters are either estimated or fixed", {
  prior <- monty::monty_dsl({
    beta ~ Exponential(mean = 2)
    unknown ~ Exponential(mean = 1)
  })
  expect_error(
    ODIN2(
      observation,
      sir,
      prior = prior,
      observation_model = odin_poisson()
    ),
    "Unknown odin2 parameter"
  )
  expect_error(
    ODIN2(
      observation,
      sir,
      prior = test_prior(),
      fixed = list(beta = 2),
      observation_model = odin_poisson()
    ),
    "both `prior` and `fixed`"
  )
})


test_that("an observation model is required", {
  expect_error(
    ODIN2(
      observation,
      sir,
      prior = test_prior(),
      observation_model = NULL
    ),
    "must be a function"
  )
})


test_that("observation helper parameters are supplied explicitly", {
  expect_error(
    ODIN2(
      observation,
      sir,
      prior = test_prior(),
      fixed = list(N = 1e6),
      observation_model = odin_negbin()
    ),
    "phi.*rho.*prior.*fixed"
  )

  expect_s3_class(
    ODIN2(
      observation,
      sir,
      prior = test_prior(),
      fixed = list(N = 1e6, rho = 0.4, phi = 20),
      observation_model = odin_negbin()
    ),
    "mdl_defn"
  )
})


test_that("computational settings are validated", {
  build <- function(...) {
    ODIN2(
      observation,
      sir,
      prior = test_prior(),
      observation_model = odin_poisson(scale = NULL),
      ...
    )
  }

  expect_error(build(n_steps = 19L), "n_steps")
  expect_error(build(n_draws = 1L), "n_draws")
  expect_error(build(n_steps = 100L, n_draws = 51L), "half")
  expect_error(build(dt = 0), "in \\(0, 1\\]")
  expect_error(build(dt = 0.3), "1 divided by an integer")
  expect_s3_class(build(dt = 0.125), "mdl_defn")
})


test_that("observation helpers are small and predictable", {
  set.seed(1)
  negbin <- odin_negbin(size = 1e6, scale = 0.5)
  expect_equal(
    negbin(c(100, 200), list()),
    c(50, 100),
    tolerance = 0.1
  )

  poisson <- odin_poisson(scale = NULL)
  expect_length(poisson(c(10, 20, 30), list()), 3L)
  expect_error(odin_negbin(size = "phi", scale = NULL)(10, list()), "phi")
})
