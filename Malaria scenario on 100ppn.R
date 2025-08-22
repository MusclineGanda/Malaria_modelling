
# collaboration during malaria modelling:

# loading packages
library(tidyverse)
library(deSolve)

### Define the model equations 
vector_human <- function(t, x, parms) {
  with(as.list(c(parms, x)), {

    # Total Populations
    M <- Sm + Im
    H <- S + I

    # Vector Equations
    dSm <- mu * M - beta * I / H * Sm - mu * Sm
    dIm <- beta * I / H * Sm - mu * Im

    # Human Equations
    dS <- -alpha * Im / H * S + gamma * I
    dI <- alpha * Im / H * S - gamma * I

    output <- c(dSm, dIm, dS, dI)
    list(output)
  })
}

# Initialize the compartments
start <- c(
  Sm = 4900,
  Im = 100,
  S = 990,
  I = 10
)

# Define parameters
parms <- c(
  mu = 0.1,
  alpha = (0.3*0.01*0.03),
  beta = 0.5,
  gamma = 1/25
)

# Define time interval
times <- seq(0, 365)

# Solve vector equations
Vector_model <- ode(
  times = times,
  parms = parms,
  func = vector_human,
  y = start
)
View(Vector_model)
# Plotting the results
vmplot <- as_tibble(as.data.frame(Vector_model)) %>%
  pivot_longer(names_to = "variable", cols = !1)

vmplot %>%
  group_by(variable) %>%
  ggplot(aes(x = time, y = value, colour = variable)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Vector Human Compartments", y = "Population", colour = "Species") +
  facet_wrap(~variable)

# Adding seasonality to the model
### Definition of equations
vector_human_seas <- function(t, x, parms) {
  with(as.list(c(parms, x)), {

    # Total Population
    M = Sm + Im
    H = S + I

    # Adding seasonality
    seas <- 1 + amp * cos(2 * pi * (t / 365 - phi))  #

    # Vector Equations
    dSm <- mu * M - seas * beta * I / H * Sm - mu * Sm
    dIm <- seas * beta * I / H * Sm - mu * Im

    # Human equations
    dS <- -seas * alpha * Im / H * S + gamma * I
    dI <- seas * alpha * Im / H * S - gamma * I

    output <- c(dSm, dIm, dS, dI)
    list(output)
  })
}

# Initialization:Initial conditions
start <- c(
  Sm = 1970,
  Im = 30,
  S = 990,
  I = 10
)

# Add new parameter
parms <- c(
  mu = 1/15,
  beta = 0.25,
  alpha = 0.12,
  amp = 0.5,
  gamma = 1/20,
  phi = 1
)

# Define time interval
times <- seq(0, 100, 0.1)

# Solve vector equations
Vector_model_seas <- ode(
  times = times,
  parms = parms,
  func = vector_human_seas,
  y = start
)

# Plotting the results
vmplot <- as_tibble(as.data.frame(Vector_model_seas)) %>%
  pivot_longer(names_to = "variable", cols = !1)

vmplot %>%
  group_by(variable) %>%
  ggplot(aes(x = time, y = value, colour = variable)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Vector Human Compartments with Seasonality", y = "Population", colour = "Species") +
  facet_wrap(~variable)
