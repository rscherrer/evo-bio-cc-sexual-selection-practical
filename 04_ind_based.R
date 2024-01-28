rm(list = ls())

library(tidyverse)
library(mnormt)

N0 <- 400
var <- 1
cov <- 0.99
c <- 0.1
tol <- 1
r <- 0.5
mut <- 0.001
mup <- 0.001
sigma <- 0.1
endtime <- 1000

varcov <- matrix(c(var, cov, cov, var), nrow = 2)

varcov <- matrix(c(0.2, 0.3, 0.3, 0.6), nrow = 2)

set.seed(42)

D <- tibble(
  time = 0:endtime,
  data = list(NULL)
)

D$data[[1]] <- rmnorm(N0, mean = c(0, 0), varcov = varcov) %>%
  as_tibble() %>%
  rename(t = "V1", p = "V2") %>%
  mutate(sex = seq(n()) %in% sample(n(), n() / 2, replace = FALSE))

for (curr_time in seq(endtime)) {

  print(curr_time)

  data <- D$data[[which(D$time == curr_time - 1)]]

  data <- data %>% mutate(survival = if_else(sex, exp(-c * t * t), 1))

  data <- data %>%
    mutate(survives = rbinom(n(), 1, prob = survival)) %>%
    filter(survives == 1) %>%
    select(-survives, -survival)

  data_m <- data %>% filter(sex)

  data <- data %>%
    filter(!sex) %>%
    mutate(

      i = map_int(p, function(p) {

        mates <- 0
        i <- NA

        while (mates == 0) {

          i <- sample(nrow(data_m), 1)

          #exp(a * p * data_m$t[i])

          prob <- exp(-(data_m$t[i] - p)^2 / tol^2)
          mates <- rbinom(1, 1, prob = prob)

        }

        return(i)

      }),

      tm = data_m$t[i],
      pm = data_m$p[i]

    ) %>%
    select(-i)

  data <- data %>%
    mutate(noff = 2) %>%
    uncount(noff) %>%
    mutate(

      cross = rbinom(n(), 1, prob = r),
      thap = rbinom(n(), 1, prob = 0.5),
      phap = if_else(cross == 1, as.numeric(!thap), thap),
      toff = if_else(thap == 1, tm, t),
      poff = if_else(phap == 1, pm, p)

    ) %>%
    select(-cross, -thap, -phap)

  data <- data %>%
    mutate(

      tmut = rbinom(n(), 1, prob = mut),
      pmut = rbinom(n(), 1, prob = mup),
      toff = toff + if_else(tmut == 1, rnorm(n(), mean = 0, sd = sigma), 0),
      poff = poff + if_else(pmut == 1, rnorm(n(), mean = 0, sd = sigma), 0)

    ) %>%
    select(-tmut, -pmut)

  data <- data %>%
    select(-t, -p, -sex, -tm, -pm) %>%
    rename(t = "toff", p = "poff") %>%
    mutate(sex = rep(c(TRUE, FALSE), n() / 2))

  D$data[[which(D$time == curr_time)]] <- data

}

D %>%
  unnest(data) %>%
  pivot_longer(t:p) %>%
  group_by(time, name) %>%
  summarize(value = mean(value)) %>%
  ggplot(aes(x = time, y = value, color = name)) +
  geom_point()
