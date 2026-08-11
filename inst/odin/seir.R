# Stochastic SEIR with a negative binomial observation model. Incidence counts
# entries to the infectious compartment, so it aligns with symptom onset rather
# than infection. See inst/odin/sir.R for the time-unit convention.

# See inst/odin/sir.R for why the seeds are rounded and the compartments
# clamped before every binomial draw.
seed_E <- max(round(E0), 0)
seed_I <- max(round(I0), 0)

initial(S) <- N - seed_E - seed_I
initial(E) <- seed_E
initial(I) <- seed_I
initial(R) <- 0
initial(incidence, zero_every = 1) <- 0

update(S) <- S - n_SE
update(E) <- E + n_SE - n_EI
update(I) <- I + n_EI - n_IR
update(R) <- R + n_IR
update(incidence) <- incidence + n_EI

p_SE <- 1 - exp(-beta * max(I, 0) / N * dt)
p_EI <- 1 - exp(-sigma * dt)
p_IR <- 1 - exp(-gamma * dt)
n_SE <- Binomial(max(S, 0), p_SE)
n_EI <- Binomial(max(E, 0), p_EI)
n_IR <- Binomial(max(I, 0), p_IR)

beta <- parameter(1)
sigma <- parameter(1)
gamma <- parameter(0.5)
N <- parameter(1e6)
E0 <- parameter(10)
I0 <- parameter(10)
rho <- parameter(1)
phi <- parameter(10)

cases <- data()
cases ~ NegativeBinomial(size = phi, mu = rho * incidence + 1e-6)
