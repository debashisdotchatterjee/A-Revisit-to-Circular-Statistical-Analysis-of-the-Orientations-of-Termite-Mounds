# Load necessary libraries
library(circular)
library(ggplot2)

# Function to compute the likelihood of the von Mises distribution
von_mises_likelihood <- function(angles, mu, kappa) {
  n <- length(angles)
  # Calculate the Bessel function I_0(kappa)
  I0_kappa <- besselI(kappa, nu = 0)
  
  # Compute the likelihood
  log_likelihood <- n * log(1 / (2 * pi * I0_kappa)) + kappa * sum(cos(angles - mu))
  
  return(log_likelihood)
}

# Function for the log-prior
log_prior <- function(mu, kappa) {
  # Prior for mu (uniform on [0, 2pi])
  if (mu < 0 || mu > 2 * pi) {
    return(-Inf)  # Log-prior is -Inf if mu is outside the range
  }
  
  # Prior for kappa (Gamma prior)
  if (kappa <= 0) {
    return(-Inf)  # Log-prior is -Inf if kappa is non-positive
  }
  
  log_prior_mu <- log(1 / (2 * pi))  # Uniform prior for mu
  log_prior_kappa <- dgamma(kappa, shape = 2, rate = 0.1, log = TRUE)  # Gamma prior for kappa
  
  return(log_prior_mu + log_prior_kappa)
}

# Posterior function (log-posterior)
log_posterior <- function(angles, mu, kappa) {
  log_likelihood <- von_mises_likelihood(angles, mu, kappa)
  log_prior_val <- log_prior(mu, kappa)
  
  return(log_likelihood + log_prior_val)
}

# Metropolis-Hastings MCMC function
mcmc_von_mises <- function(angles, n_iter = 5000, burn_in = 1000, proposal_sd_mu = 0.1, proposal_sd_kappa = 0.5) {
  # Initialize the chain
  chain_mu <- numeric(n_iter)
  chain_kappa <- numeric(n_iter)
  
  # Starting values for mu and kappa
  chain_mu[1] <- runif(1, 0, 2 * pi)  # Initialize mu
  chain_kappa[1] <- rgamma(1, shape = 2, rate = 0.1)  # Initialize kappa
  
  # Metropolis-Hastings algorithm
  for (i in 2:n_iter) {
    # Propose new mu and kappa
    proposed_mu <- rnorm(1, mean = chain_mu[i - 1], sd = proposal_sd_mu)
    proposed_mu <- (proposed_mu %% (2 * pi))  # Ensure mu stays within [0, 2*pi]
    
    proposed_kappa <- abs(rnorm(1, mean = chain_kappa[i - 1], sd = proposal_sd_kappa))  # Ensure kappa is positive
    
    # Compute the log-posterior for the current and proposed parameters
    current_log_posterior <- log_posterior(angles, chain_mu[i - 1], chain_kappa[i - 1])
    proposed_log_posterior <- log_posterior(angles, proposed_mu, proposed_kappa)
    
    # Acceptance probability (Metropolis-Hastings ratio)
    acceptance_prob <- exp(proposed_log_posterior - current_log_posterior)
    acceptance_prob <- min(1, acceptance_prob)
    
    # Accept or reject the proposal
    if (runif(1) < acceptance_prob) {
      chain_mu[i] <- proposed_mu
      chain_kappa[i] <- proposed_kappa
    } else {
      chain_mu[i] <- chain_mu[i - 1]
      chain_kappa[i] <- chain_kappa[i - 1]
    }
  }
  
  # Remove burn-in period
  chain_mu <- chain_mu[(burn_in + 1):n_iter]
  chain_kappa <- chain_kappa[(burn_in + 1):n_iter]
  
  return(list(mu = chain_mu, kappa = chain_kappa))
}

# Function to plot the posterior density of mu as a rose plot
plot_mu_posterior_rose <- function(mu_samples, set_name) {
  # Convert mu samples to a circular object
  mu_circular <- circular(mu_samples, units = "radians", template = "none")
  
  # Plot the rose plot (circular histogram)
  p <- ggplot(data = data.frame(mu = as.numeric(mu_circular)), aes(x = mu)) +
    geom_histogram(binwidth = pi / 18, fill = "lightblue", color = "black") + 
    coord_polar(theta = "x") +
    ggtitle(paste("Posterior Density of Mu for", set_name)) +
    theme_minimal()
  
  # Save the plot
  ggsave(paste0(set_name, "_mu_posterior_rose_plot.png"), plot = p, width = 6, height = 6)
}

# Loop through all sets in fisherB13c
set_names <- names(fisherB13c)
posterior_results <- list()

for (set_name in set_names) {
  # Extract the angles for the current set and convert to radians
  angles <- as.numeric(fisherB13c[[set_name]])
  angles_rad <- circular::rad(angles)
  
  # Run MCMC for the current set
  posterior <- mcmc_von_mises(angles = angles_rad, n_iter = 5000, burn_in = 1000)
  
  # Store the results
  posterior_results[[set_name]] <- posterior
  
  # Plot the posterior density of mu as a rose plot
  plot_mu_posterior_rose(posterior$mu, set_name)
}
