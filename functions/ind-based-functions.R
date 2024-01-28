# Function to run the simulation
sim_model_ibm <- function(

  N, t0, p0, varcov, r, tend, c = 0.1, topt = 1, a = 1,
  mut = 0.01, mup = 0.01, sigma = 0.1, light = TRUE

) {

  # Initialization
  D <- tibble(
    time = 0:tend,
    data = list(NULL)
  )

  # Sample the starting population from a bivariate normal distribution
  D$data[[1]] <- mnormt::rmnorm(N, mean = c(t0, p0), varcov = varcov) %>%
    as_tibble() %>%
    rename(t = "V1", p = "V2") %>%
    mutate(sex = seq(n()) %in% sample(n(), n() / 2, replace = FALSE))

  # For each time step...
  for (curr_time in seq(tend)) {

    print(curr_time)

    # Extract the data from the previous time step
    data <- D$data[[which(D$time == curr_time - 1)]]

    # Compute chances of survival of males based on ornament
    data <- data %>% mutate(survival = if_else(sex, exp(-c * (t - topt) * (t - topt)), 1))

    # Remove all the males that do not survive
    data <- data %>%
      mutate(survives = rbinom(n(), 1, prob = survival)) %>%
      filter(survives == 1) %>%
      select(-survives, -survival)

    # Extra the male data
    data_m <- data %>% filter(sex)

    # For each female...
    data <- data %>%
      filter(!sex) %>%
      mutate(

        # Pick a mate based on the psychophysical model...
        i = map_int(p, function(p) {

          # By weighing males depending on their ornamennt
          probs <- exp(a * p * data_m$t)
          i <- sample(seq(probs), size = 1, prob = probs)

          return(i)

        }),

        # Ornamenent and preference of the mate
        tm = data_m$t[i],
        pm = data_m$p[i]

      ) %>%
      select(-i)

    # For each offspring made (two per mating pair)...
    data <- data %>%
      mutate(noff = 2) %>%
      uncount(noff) %>%
      mutate(

        # Is there recombination?
        cross = rbinom(n(), 1, prob = r),

        # Which haplotypes are inherited from the parents?
        thap = rbinom(n(), 1, prob = 0.5),
        phap = if_else(cross == 1, as.numeric(!thap), thap),

        # What are the trait values in the offspring?
        toff = if_else(thap == 1, tm, t),
        poff = if_else(phap == 1, pm, p)

      ) %>%
      select(-cross, -thap, -phap)

    # Now each offspring mutates with a certain probability
    data <- data %>%
      mutate(

        # Pick which offspring mutates
        tmut = rbinom(n(), 1, prob = mut),
        pmut = rbinom(n(), 1, prob = mup),

        # Then implement a phenotypic deviation upon mutation
        toff = toff + rnorm(n(), mean = 0, sd = sigma) * as.numeric(tmut == 1),
        poff = poff + rnorm(n(), mean = 0, sd = sigma) * as.numeric(pmut == 1)

      ) %>%
      select(-tmut, -pmut)

    # Rearrange the table
    data <- data %>%
      select(-t, -p, -sex, -tm, -pm) %>%
      rename(t = "toff", p = "poff")

    # Assign sexes (balanced sex ratio)
    data <- data %>%
      mutate(sex = rep(c(TRUE, FALSE), n() / 2))

    # Append the data to the master database
    D$data[[which(D$time == curr_time)]] <- data

    # Keep only the mean traits if needed
    if (light) D$data[[which(D$time == curr_time - 1)]] <- data %>%
      summarize(t = mean(t), p = mean(p))

  }

  if (light) D$data[[nrow(D)]] <- D$data[[nrow(D)]] %>%
    summarize(t = mean(t), p = mean(p))

  return(D %>% unnest(data))

}
