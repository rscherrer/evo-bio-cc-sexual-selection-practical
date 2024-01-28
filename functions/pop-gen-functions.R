# Function to get the change in ornament frequency
get_dtfreq <- function(t, A) { 0.5 * t * (1 - t) * A}

# Function to get the change in preference frequency
get_dpfreq <- function(D, A) { 0.5 * D * A }

# Function to compute the change in linkage disequilibrium
get_dD <- function(t, p, D, A, B, r) {

  D * ((1 - r) * A * (0.5 - t) - 0.25 * A * A * t * (1 - t) - r) +
    0.5 * r * B * (D * D + p * (1 - p) * t * (1 - t))

}

# Functions to get some particular quantities
get_A <- function(t, p, tm, U01, U11) { (p * (U11 - U01) - (t - tm)) / (t * (1 - t)) }
get_B <- function(t, U01, U11) { (U11 - U01) / (t * (1 - t)) }

# Function to get ornament frequency after selection
get_tm <- function(t, s) {

  v <- 1 - s
  vbar <- 1 - s * t
  t * v / vbar

}

# Functions to get different assortment quantities
get_U00 <- function(tm) { 1 - tm }
get_U01 <- function(tm) { tm }
get_U10 <- function(tm, a) { (1 - tm) / (1 - tm + a * tm) }
get_U11 <- function(tm, a) { a * tm / (1 - tm + a * tm) }

# Function to compute allele frequencies and linkage at the next generation
get_nextgen_pgn <- function(t, p, D, s, a, r) {

  # Compute ornament frequency after viability selection
  tm <- get_tm(t, s)

  # Get assortment quantities
  U01 <- get_U01(tm)
  U11 <- get_U11(tm, a)

  # Compute some other useful quantities
  A <- get_A(t, p, tm, U01, U11)
  B <- get_B(t, U01, U11)

  # Compute the change in allele frequencies and linkage
  dt <- get_dtfreq(t, A)
  dp <- get_dpfreq(D, A)
  dD <- get_dD(t, p, D, A, B, r)

  # Update allele frequencies and linkage
  newt <- t + dt
  newp <- p + dp
  newD <- D + dD

  return(list(t = newt, p = newp, D = newD))

}

# Function to run a simulation
sim_model_pgn <- function(t0, p0, D0, r, tend, s = 0.1, a = 1.5) {

  # Initialization
  t <- t0
  p <- p0
  D <- D0
  data <- tibble(time = 0, t = t, p = p, D = D)

  # For each time step...
  for (i in 1:tend) {

    # Compute the variables at the next generation
    newgen <- get_nextgen_pgn(t, p, D, s, a, r)

    # Update the system
    t <- newgen$t
    p <- newgen$p
    D <- newgen$D

    t <- if_else(t < 0, 0, if_else(t > 1, 1, t))
    p <- if_else(p < 0, 0, if_else(p > 1, 1, p))

    # Append to a table
    data <- data %>% add_row(time = i, t = t, p = p, D = D)

  }

  return(data)

}

# Function to compute the equilibrium ornament frequency given a preference
get_teq_pgn <- function(p, s = 0.1, a = 1.5) {

  teq <- (p * (1 - s) * (a - 1) / s - 1) / (a * (1 - s) - 1)

  # Clamp between zero and one
  teq[teq < 0] <- 0
  teq[teq > 1] <- 1

  return(teq)

}
