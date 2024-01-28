rm(list = ls())

library(tidyverse)

get_beta <- function(t, p, a, theta, omega) { a * p - (t - theta) / (omega * omega) }

get_dt <- function(t, p, Gt, beta) { 0.5 * Gt * beta }

get_dp <- function(t, p, Gtp, beta) { 0.5 * Gtp * beta }

get_nextgen <- function(t, p, Gt, Gtp, a, theta, omega) {

  beta <- get_beta(t, p, a, theta, omega)

  dt <- get_dt(t, p, Gt, beta)
  dp <- get_dp(t, p, Gtp, beta)

  newt <- t + dt
  newp <- p + dp

  return(list(t = newt, p = newp))

}

get_nextgen(0, 0, 1, 1, 1, 1, 1)

sim_model <- function(t0, p0, Gt, Gtp, a, theta, omega, tend) {

  t <- t0
  p <- p0

  data <- tibble(time = 0, t = t, p = p)

  for (i in 1:tend) {

    newgen <- get_nextgen(t, p, Gt, Gtp, a, theta, omega)

    t <- newgen$t
    p <- newgen$p

    data <- data %>% add_row(time = i, t = t, p = p)

  }

  return(data)

}

data <- sim_model(
  t0 = 1,
  p0 = 2,
  Gt = 0.1,
  Gtp = 0.11,
  a = 1,
  theta = 0,
  omega = 2,
  tend = 100
)

data %>%
  pivot_longer(t:p) %>%
  ggplot(aes(x = time, y = value, group = name, color = name)) +
  geom_line()

get_teq <- function(p, a, theta, omega) {

  p * a * omega * omega + theta

}

equil <- tibble(
  p = seq(0, 5, 0.01),
  t = get_teq(p, a = 1, theta = 0, omega = 2)
)

plot <- equil %>%
  ggplot(aes(x = p, y = t)) +
  geom_line()

plot

plot +
  geom_point(data = data, aes(color = time))
