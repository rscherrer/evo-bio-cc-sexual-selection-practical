rm(list = ls())

library(tidyverse)

get_gradt <- function(s, vm0, vm1) { (1 - s) * vm1 - vm0 }

get_gradp <- function(bpt, s, vm0, vm1) { bpt * ((1 - s) * vm1 - vm0) }

get_nextgen <- function(t, p, bpt, s, vm0, vm1, z, M, step) {

  gradt <- get_gradt(s, vm0, vm1)
  gradp <- get_gradp(bpt, s, vm0, vm1)

  deriv <- c(z * M %*% c(gradt, gradp))

  dt <- deriv[1]
  dp <- deriv[2]

  a <- dt
  b <-
  newt <- t + step * (a + 2 * b + 2 * c + d) / 6


  newt <- t + dt
  newp <- p + dp

  return(list(t = newt, p = newp))

}

get_nextgen(1, 1, 0.6, 0.3, 1, 1)

sim_model <- function(t0, p0, bpt, s, vm0, vm1, tend) {

  t <- t0
  p <- p0

  data <- tibble(time = 0, t = t, p = p)

  for (i in 1:tend) {

    newgen <- get_nextgen(t, p, bpt, s, vm0, vm1)

    t <- newgen$t
    p <- newgen$p

    data <- data %>% add_row(time = i, t = t, p = p)

  }

  return(data)

}

data <- sim_model(1, 1, 0.6, 0.3, 1, 1, 100)

data %>%
  pivot_longer(t:p) %>%
  ggplot(aes(x = time, y = value, group = name, color = name)) +
  geom_line()

get_teq <- function(p, a, theta, omega) {

  p * a * omega * omega + theta

}
