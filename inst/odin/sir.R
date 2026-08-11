# Stochastic SIR with a negative binomial observation model.
#
# Time is measured in REPORTING INTERVALS, not days, so the same compiled model
# serves weekly, daily or monthly data. `beta` and `gamma` are therefore rates
# per interval, and `incidence` resets every interval (zero_every = 1), which
# `odin` requires to be an integer literal.

# `I0` is estimated on a continuous scale, so it is rounded: a fractional
# compartment lets a binomial draw exceed the compartment it drains, which
# sends the next state negative and then aborts the run on a negative
# probability. `max(I, 0)` guards the same failure for a state set from outside
# the model, as the forecast does.
seed <- max(round(I0), 0)

initial(S) <- N - seed
initial(I) <- seed
initial(R) <- 0
initial(incidence, zero_every = 1) <- 0

update(S) <- S - n_SI
update(I) <- I + n_SI - n_IR
update(R) <- R + n_IR
update(incidence) <- incidence + n_SI

p_SI <- 1 - exp(-beta * max(I, 0) / N * dt)
p_IR <- 1 - exp(-gamma * dt)
n_SI <- Binomial(max(S, 0), p_SI)
n_IR <- Binomial(max(I, 0), p_IR)

beta <- parameter(1)
gamma <- parameter(0.5)
N <- parameter(1e6)
I0 <- parameter(10)
rho <- parameter(1)
phi <- parameter(10)

cases <- data()
cases ~ NegativeBinomial(size = phi, mu = rho * incidence + 1e-6)
