# Function to compute the beta term
get_beta <- function(t, p, a, theta, omega) { a * p - (t - theta) / (omega * omega) }

# Function to compute the selection gradient for ornaments
get_dt <- function(t, p, Gt, beta) { 0.5 * Gt * beta }

# Same for preferences
get_dp <- function(t, p, Gtp, beta) { 0.5 * Gtp * beta }

# Function to compute traits at the next generation
get_nextgen_qgn <- function(t, p, Gt, Gtp, a, theta, omega) {

  # Compute the necessary quantities
  beta <- get_beta(t, p, a, theta, omega)
  dt <- get_dt(t, p, Gt, beta)
  dp <- get_dp(t, p, Gtp, beta)

  # Update trait values
  newt <- t + dt
  newp <- p + dp

  return(list(t = newt, p = newp))

}

# Function to simulate the model over generations
sim_model_qgn <- function(t0, p0, Gt, Gtp, tend, a = 1, theta = 2, omega = 1) {

  # Initialization
  t <- t0
  p <- p0
  data <- tibble(time = 0, t = t, p = p)

  # For each time step...
  for (i in 1:tend) {

    # Compute the next generation
    newgen <- get_nextgen_qgn(t, p, Gt, Gtp, a, theta, omega)

    # Update the trait values
    t <- newgen$t
    p <- newgen$p

    # Store them in a table
    data <- data %>% add_row(time = i, t = t, p = p)

  }

  return(data)

}

# Function to compute the equilibrium ornament value given a preference
get_teq_qgn <- function(p, a = 1, theta = 2, omega = 1) {

  p * a * omega * omega + theta

}
