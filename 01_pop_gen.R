rm(list = ls())

library(tidyverse)

theme_set(theme_classic())

get_dt <- function(t, A) { 0.5 * t * (1 - t) * A}

get_dp <- function(D, A) { 0.5 * D * A }

get_dD <- function(t, p, D, A, B, r) {

  D * ((1 - r) * A * (0.5 - t) - 0.25 * A * A * t * (1 - t) - r) +
    0.5 * r * B * (D * D + p * (1 - p) * t * (1 - t))

}

get_A <- function(t, p, tm, U01, U11) { (p * (U11 - U01) - (t - tm)) / (t * (1 - t)) }

get_B <- function(t, U01, U11) { (U11 - U01) / (t * (1 - t)) }

get_tm <- function(t, s) {

  v <- 1 - s
  vbar <- 1 - s * t
  t * v / vbar

}

get_U00 <- function(tm) 1 - tm
get_U01 <- function(tm) tm
get_U10 <- function(tm, a) (1 - tm) / (1 - tm + a * tm)
get_U11 <- function(tm, a) a * tm / (1 - tm + a * tm)

s <- 0.4
a <- 3
r <- 0.5

t <- 0.1
p <- 0.1
D <- 0.5

get_nextgen <- function(t, p, D, s, a, r) {

  tm <- get_tm(t, s)
  U01 <- get_U01(tm)
  U11 <- get_U11(tm, a)

  A <- get_A(t, p, tm, U01, U11)
  B <- get_B(t, U01, U11)

  dt <- get_dt(t, A)
  dp <- get_dp(D, A)
  dD <- get_dD(t, p, D, A, B, r)

  newt <- t + dt
  newp <- p + dp
  newD <- D + dD

  return(list(t = newt, p = newp, D = newD))

}

get_nextgen(0.1, 0.1, 0.5, 0.4, 3, 0.5)

sim_model <- function(t0, p0, D0, s, a, r, tend) {

  t <- t0
  p <- p0
  D <- D0

  data <- tibble(time = 0, t = t, p = p, D = D)

  for (i in 1:tend) {

    newgen <- get_nextgen(t, p, D, s, a, r)

    t <- newgen$t
    p <- newgen$p
    D <- newgen$D

    data <- data %>% add_row(time = i, t = t, p = p, D = D)

  }

  return(data)

}

data <- sim_model(
  t0 = 0.8,
  p0 = 0.5,
  D0 = 0.2,
  s = 0.4,
  a = 3,
  r = 0.2,
  tend = 1000
)

data %>%
  pivot_longer(t:p) %>%
  ggplot(aes(x = time, y = value, group = name, color = name)) +
  geom_line() +
  ylim(c(0, 1))

get_teq <- function(p, s, a) {

  teq <- (p * (1 - s) * (a - 1) / s - 1) / (a * (1 - s) - 1)
  teq[teq < 0] <- 0
  teq[teq > 1] <- 1

  return(teq)

}

equil <- tibble(
  p = seq(0, 1, 0.01),
  t = get_teq(p, s = 0.4, a = 3)
)

plot <- equil %>%
  ggplot(aes(x = p, y = t)) +
  geom_line()

plot

plot +
  geom_point(data = data, aes(color = time))
