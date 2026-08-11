# incast.odin

`incast.odin` lets an [odin2](https://mrc-ide.github.io/odin2/) transmission
model behave like any other forecasting model in
[incast](https://github.com/ACCIDDA/incast).

## Quick start

The ready-to-use SIR model follows the ordinary `incast` pipeline:

```r
library(incast)
library(incast.odin)

models <- list(
  SIR = odin_sir(observation, population = 19.5e6)
)

set.seed(1)

forecast <- example_data |>
  dplyr::filter(
    location == "NY",
    target_end_date >= as.Date("2024-10-19"),
    as_of <= as.Date("2025-02-23")
  ) |>
  check_data() |>
  get_ncast(max_delay = 2, draws = 200) |>
  get_cv(h = 3, n_origins = 3, step = 1, models = models) |>
  get_fcast(top_n = 1)

autoplot(forecast)
```

That is all that is required to use the built-in model. `odin_seir()` provides
the same interface with an exposed compartment. The date filters isolate the
2024–25 New York influenza wave and recreate the data available on 23 February
2025. The vignette adds influenza-informed priors and compares SIR with the
default statistical models.

## What happens

No one package performs the whole job:

| Package | Why it is needed |
| --- | --- |
| `odin2` | describes the dynamical model and its comparison with data |
| `dust2` | runs the odin2 model and calculates its likelihood |
| `monty` | combines the likelihood with a prior and samples from the posterior |
| `incast` | cross-validates, scores, ensembles and formats the forecast |

`incast.odin::ODIN2()` is the small adapter between these pieces. It returns a
normal `fable` model definition, so `incast` needs no special odin2 workflow.

Users write custom models with `odin2` and priors with `monty`. `ODIN2()` calls
the required dust2 and MCMC code internally.

The current interface uses dust2's deterministic likelihood and projection.
Forecast uncertainty comes from the posterior parameter draws and the
observation model. Particle filtering is intentionally outside this first,
teaching-focused interface.

## Write your own SIR model

The complete workflow has four steps:

1. Write and compile the transmission model with `odin2::odin()`.
2. Describe the parameter prior with `monty::monty_dsl()`.
3. Connect them with `ODIN2()`.
4. Put the result in the usual `models` list.

For the complete model, a prior suited to the example data, the full
`ncast |> cv |> fcast` pipeline and an explanation of each line, read
`vignette("incast-odin2")`.

## Three conventions

A model passed to `ODIN2()` is deliberately simple and well defined:

- it is a discrete-time odin2 model;
- observed counts are declared as `cases <- data()`;
- the scalar state being forecast is named `incidence`.

One model-time unit equals one reporting interval. For weekly input, rates are
per week; for daily input, rates are per day.

## Installation

```r
pak::pak("ACCIDDA/incast.odin")
```

The built-in SIR and SEIR models are precompiled. The `odin2` package and a
compiler are needed only when writing a new model.

## Development

After editing a model in `inst/odin/`, regenerate its compiled sources with:

```r
Rscript data-raw/generate.R
devtools::document()
```

## Licence

MIT
