# Load necessary libraries
library(circular)
library(CircStats)
library(ggplot2)
library(gridExtra)
library(tibble)
library(dplyr)

# Create a directory to save all outputs
output_dir <- "termite_orientation_analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Load the dataset
data(fisherB13c)
data(fisherB13c)
plot(fisherB13c$set1, stack=TRUE, shrink=1.5)


# Load necessary library
library(circular)

# Create a summary of the sites
site_summary <- data.frame(
  Site = 1:14,
  Latitude = c("-15°43'", "-15°32'", "-15°38'", "-15°40'", "-15°42'", 
               "-15°36'", "-15°35'", "-15°39'", "-15°33'", "-15°37'", 
               "-15°41'", "-15°31'", "-15°34'", "-15°30'"),
  Longitude = c("144°42'", "144°17'", "144°25'", "144°38'", "144°40'", 
                "144°18'", "144°19'", "144°22'", "144°23'", "144°20'", 
                "144°35'", "144°36'", "144°37'", "144°39'"),
  Observations = sapply(fisherB13c, length) # Assuming length gives the number of observations per set
)

# Print the table to a CSV file
output_file <- "termite_orientation_analysis/Table_1.csv"
write.csv(site_summary, output_file, row.names = FALSE)

# Print the summary to console for verification
print(site_summary)

# Define a helper function to save plots and tables
save_output <- function(plot, filename) {
  ggsave(paste0(output_dir, "/", filename), plot = plot, width = 8, height = 6)
}

# Set up a list to store the outputs (e.g., tables, summaries)
output_list <- list()

# Function to analyze each dataset in fisherB13c
analyze_circular_data <- function(dataset, set_name) {
  
  # Plot the circular data
  plot_title <- paste("Circular Plot for", set_name)
  plot_data <- ggplot() +
    geom_histogram(
      aes(x = as.numeric(dataset)),
      binwidth = 10, fill = "lightblue", color = "black"
    ) +
    coord_polar(theta = "x") +
    ggtitle(plot_title) +
    theme_minimal()
  
  # Save the plot
  save_output(plot_data, paste0(set_name, "_circular_plot.png"))
  
  # Circular summary statistics
  mean_dir <- mean.circular(dataset)
  circ_var <- var.circular(dataset)
  
  # Fit a von Mises distribution to the data
  von_mises_fit <- mle.vonmises(dataset)
  mu <- von_mises_fit$mu
  kappa <- von_mises_fit$kappa
  
  # Rayleigh test for uniformity
  rayleigh_test <- rayleigh.test(dataset)
  
  # Summary output
  summary_data <- tibble(
    Set = set_name,
    Mean_Direction = as.numeric(mean_dir),
    Circular_Variance = as.numeric(circ_var),
    Von_Mises_Mean = as.numeric(mu),
    Von_Mises_Kappa = as.numeric(kappa),
    Rayleigh_Test_Statistic = rayleigh_test$statistic,
    Rayleigh_Test_p_value = rayleigh_test$p.value
  )
  
  # Add to output list
  output_list[[set_name]] <- summary_data
  
  # Print the results
  print(summary_data)
  
  # Save the table to a CSV
  write.csv(summary_data, file = paste0(output_dir, "/", set_name, "_summary.csv"), row.names = FALSE)
  
  # Return the summary data for further use
  return(summary_data)
}

# Loop through all sets in fisherB13c and apply the analysis
all_summaries <- lapply(names(fisherB13c), function(set_name) {
  dataset <- fisherB13c[[set_name]]
  analyze_circular_data(dataset, set_name)
})

# Combine all summaries into a single table
combined_summary <- do.call(rbind, all_summaries)

# Save the combined summary as a CSV file
write.csv(combined_summary, file = paste0(output_dir, "/combined_summary.csv"), row.names = FALSE)

# Print the combined summary
print(combined_summary)

# Generate a table of all the model fits and save it
combined_summary_plot <- ggplot(combined_summary, aes(x = Set, y = Von_Mises_Kappa)) +
  geom_bar(stat = "identity", fill = "lightgreen", color = "black") +
  ggtitle("Von Mises Kappa for each Set") +
  xlab("Set") + ylab("Von Mises Kappa") +
  theme_minimal()

save_output(combined_summary_plot, "von_mises_kappa_comparison.png")

# Print the combined plot
print(combined_summary_plot)

############################


# Load necessary libraries
library(circular)
library(CircStats)
library(ggplot2)
library(gridExtra)
library(dplyr)

# Create a directory to save all outputs
output_dir <- "extended_termite_orientation_analysis_with_degrees"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Define a helper function to save plots and tables
save_output <- function(plot, filename) {
  ggsave(paste0(output_dir, "/", filename), plot = plot, width = 8, height = 6)
}

# Set up a list to store the outputs (e.g., tables, summaries)
output_list <- list()

# Function to analyze each dataset in fisherB13c and generate multiple plots
analyze_circular_data <- function(dataset, set_name) {
  
  # 1. Circular Histogram Plot (Rose Diagram)
  plot_title <- paste("Circular Histogram for", set_name)
  rose_plot <- ggplot() +
    geom_histogram(
      aes(x = as.numeric(dataset)),
      binwidth = 15, fill = "lightblue", color = "black"
    ) +
    coord_polar(theta = "x") +
    ggtitle(plot_title) +
    theme_minimal()
  
  save_output(rose_plot, paste0(set_name, "_circular_histogram.png"))
  
  # 2. Circular Density Plot
  density_plot <- ggplot() +
    geom_density(aes(x = as.numeric(dataset)), fill = "blue", alpha = 0.3) +
    coord_polar(theta = "x") +
    ggtitle(paste("Circular Density Plot for", set_name)) +
    theme_minimal()
  
  save_output(density_plot, paste0(set_name, "_density_plot.png"))
  
  # Circular summary statistics
  mean_dir <- mean.circular(dataset)  # In degrees, no conversion needed
  circ_var <- var.circular(dataset)
  
  # Fit a von Mises distribution to the data
  von_mises_fit <- mle.vonmises(dataset)
  mu <- von_mises_fit$mu
  kappa <- von_mises_fit$kappa
  
  # Rayleigh test for uniformity
  rayleigh_test <- rayleigh.test(dataset)
  
  # Summary output
  summary_data <- tibble(
    Set = set_name,
    Mean_Direction = as.numeric(mean_dir),   # Using mean direction in degrees directly
    Circular_Variance = as.numeric(circ_var),
    Von_Mises_Mean = as.numeric(mu),
    Von_Mises_Kappa = as.numeric(kappa),
    Rayleigh_Test_Statistic = rayleigh_test$statistic,
    Rayleigh_Test_p_value = rayleigh_test$p.value
  )
  
  # Add to output list
  output_list[[set_name]] <- summary_data
  
  # Print the results
  print(summary_data)
  
  # Save the table to a CSV
  write.csv(summary_data, file = paste0(output_dir, "/", set_name, "_summary.csv"), row.names = FALSE)
  
  # Return the summary data for further use
  return(summary_data)
}

# Loop through all sets in fisherB13c and apply the analysis
all_summaries <- lapply(names(fisherB13c), function(set_name) {
  dataset <- fisherB13c[[set_name]]
  analyze_circular_data(dataset, set_name)
})

# Combine all summaries into a single table
combined_summary <- do.call(rbind, all_summaries)

# Save the combined summary as a CSV file
write.csv(combined_summary, file = paste0(output_dir, "/combined_summary.csv"), row.names = FALSE)

# Print the combined summary
print(combined_summary)

# 3. Polar Plot of Mean Directions for all sets with arrows (directly in degrees)
polar_plot_data <- data.frame(
  Set = combined_summary$Set,
  Mean_Direction = combined_summary$Mean_Direction
)

# Set arrow length and x, y positions for each segment (using degrees)
polar_plot_data$arrow_x <- cos(polar_plot_data$Mean_Direction * pi / 180)  # cos of degree
polar_plot_data$arrow_y <- sin(polar_plot_data$Mean_Direction * pi / 180)  # sin of degree

# Plot with arrows denoting mean directions
polar_plot_with_arrows <- ggplot() +
  geom_segment(data = polar_plot_data, aes(x = 0, y = 0, xend = arrow_x, yend = arrow_y),
               arrow = arrow(length = unit(0.2, "inches")), size = 1, color = "blue") +
  coord_polar(theta = "x", start = 0) +
  ggtitle("Mean Directions for All Sets with Arrows") +
  theme_minimal() +
  xlab("Angle (Degrees)") + ylab("Direction")

save_output(polar_plot_with_arrows, "mean_directions_polar_plot_with_arrows.png")
print(polar_plot_with_arrows)

# 4. Bar Plot Comparing Circular Variance for Each Set
variance_bar_plot <- ggplot(combined_summary, aes(x = Set, y = Circular_Variance)) +
  geom_bar(stat = "identity", fill = "lightgreen", color = "black") +
  ggtitle("Comparison of Circular Variance Across All Sets") +
  xlab("Set") + ylab("Circular Variance") +
  theme_minimal()

save_output(variance_bar_plot, "circular_variance_comparison.png")
print(variance_bar_plot)

# 5. Bar Plot Comparing Concentration Parameter (Kappa) for Each Set
kappa_bar_plot <- ggplot(combined_summary, aes(x = Set, y = Von_Mises_Kappa)) +
  geom_bar(stat = "identity", fill = "lightcoral", color = "black") +
  ggtitle("Comparison of Von Mises Kappa Across All Sets") +
  xlab("Set") + ylab("Von Mises Kappa") +
  theme_minimal()

save_output(kappa_bar_plot, "von_mises_kappa_comparison.png")
print(kappa_bar_plot)


######################

# Load required libraries
library(circular)
library(ggplot2)

# Create a directory to save all outputs
output_dir <- "triangular_fit_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Helper function to fit triangular distribution
fit_triangular_distribution <- function(data) {
  # Convert data to radians if needed
  data_rad <- circular(data, units = "degrees")
  
  # Estimate mean direction and circular variance
  mean_direction <- mean.circular(data_rad)
  circ_variance <- var.circular(data_rad)
  
  # Fit triangular distribution: PDF as (1 + 2 * cos(theta - mu)) / 2pi
  # We use a basic estimate: mean direction is the peak
  fitted_values <- (1 + 2 * cos(as.numeric(data_rad) - as.numeric(mean_direction))) / (2 * pi)
  
  # Goodness of fit: Calculate MAD (Mean Absolute Deviation)
  mad_value <- mean(abs(as.numeric(data_rad) - as.numeric(mean_direction)))
  
  return(list(mean_direction = mean_direction, circular_variance = circ_variance, MAD = mad_value, fitted_values = fitted_values))
}

# Function to plot and save the results
plot_and_save_fit <- function(data, set_name, fit_result) {
  # Create a circular plot comparing the actual data and the fit
  plot_data <- ggplot() +
    geom_histogram(aes(x = as.numeric(data)), binwidth = 5, fill = "lightblue", color = "black") +
    coord_polar(theta = "x") +
    ggtitle(paste("Triangular Circular Fit for", set_name)) +
    theme_minimal()
  
  # Save the plot
  ggsave(paste0(output_dir, "/", set_name, "_triangular_fit_plot.png"), plot = plot_data, width = 8, height = 6)
  
  # Create a summary table and save it
  summary_table <- data.frame(
    Set = set_name,
    Mean_Direction = as.numeric(fit_result$mean_direction),
    Circular_Variance = fit_result$circular_variance,
    MAD = fit_result$MAD
  )
  
  write.csv(summary_table, file = paste0(output_dir, "/", set_name, "_summary.csv"), row.names = FALSE)
  
  # Print the summary table to console
  print(summary_table)
}

# Analyze each set in the fisherB13c dataset
for (set_name in names(fisherB13c)) {
  dataset <- fisherB13c[[set_name]]
  
  # Fit the triangular distribution
  fit_result <- fit_triangular_distribution(dataset)
  
  # Plot and save the results
  plot_and_save_fit(dataset, set_name, fit_result)
}

# Combine all summaries into a single table
combined_summary <- do.call(rbind, lapply(names(fisherB13c), function(set_name) {
  read.csv(file.path(output_dir, paste0(set_name, "_summary.csv")))
}))

# Save the combined summary
write.csv(combined_summary, file = paste0(output_dir, "/combined_summary.csv"), row.names = FALSE)

# Print the combined summary
print(combined_summary)

############################

# Load necessary libraries
library(circular)
library(ggplot2)

# Create a directory to save all outputs
output_dir <- "triangular_fit1_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Helper function to fit triangular distribution
fit_triangular_distribution <- function(data) {
  # Convert data to radians
  data_rad <- circular(data, units = "degrees")
  
  # Estimate mean direction and circular variance
  mean_direction <- mean.circular(data_rad)
  circ_variance <- var.circular(data_rad)
  
  # Fit triangular distribution: PDF as (1 + 2 * cos(theta - mu)) / 2pi
  # The peak of the triangular distribution is around the mean direction (estimated mode)
  fitted_values <- (1 + 2 * cos(as.numeric(data_rad) - as.numeric(mean_direction))) / (2 * pi)
  
  # Goodness of fit: Calculate Mean Absolute Deviation (MAD)
  mad_value <- mean(abs(as.numeric(data_rad) - as.numeric(mean_direction)))
  
  return(list(mean_direction = mean_direction, circular_variance = circ_variance, MAD = mad_value, fitted_values = fitted_values))
}

# Function to plot and save the results with triangular density overlay
plot_and_save_triangular_fit <- function(data, set_name, fit_result) {
  # Convert data to radians
  data_rad <- circular(data, units = "degrees")
  
  # Triangular distribution density function
  triangular_density <- function(theta, mu) {
    (1 + 2 * cos(theta - mu)) / (2 * pi)
  }
  
  # Generate values for plotting the triangular density curve
  theta_seq <- seq(0, 2 * pi, length.out = 360)  # Sequence of angles in radians
  density_values <- triangular_density(theta_seq, as.numeric(fit_result$mean_direction))
  
  # Create a data frame for the density plot
  density_df <- data.frame(theta = theta_seq, density = density_values)
  
  # Plot the circular data with triangular density overlaid
  plot_data <- ggplot() +
    geom_histogram(aes(x = as.numeric(data_rad)), binwidth = 10, fill = "lightblue", color = "black") +
    coord_polar(theta = "x") +
    geom_line(data = density_df, aes(x = theta * 180 / pi, y = density * max(table(data_rad)) / max(density)), 
              color = "red", size = 1) +  # Overlay triangular density
    ggtitle(paste("Triangular Circular Fit for", set_name)) +
    theme_minimal() +
    xlab("Angle (degrees)") +
    ylab("Density")
  
  # Save the plot
  ggsave(paste0(output_dir, "/", set_name, "_triangular_fit_plot.png"), plot = plot_data, width = 8, height = 6)
  
  # Create a summary table and save it
  summary_table <- data.frame(
    Set = set_name,
    Mean_Direction = as.numeric(fit_result$mean_direction),
    Circular_Variance = fit_result$circular_variance,
    MAD = fit_result$MAD
  )
  
  write.csv(summary_table, file = paste0(output_dir, "/", set_name, "_summary.csv"), row.names = FALSE)
  
  # Print the summary table to console
  print(summary_table)
}

# Analyze each set in the fisherB13c dataset
for (set_name in names(fisherB13c)) {
  dataset <- fisherB13c[[set_name]]
  
  # Fit the triangular distribution
  fit_result <- fit_triangular_distribution(dataset)
  
  # Plot and save the results
  plot_and_save_triangular_fit(dataset, set_name, fit_result)
}

# Combine all summaries into a single table
combined_summary <- do.call(rbind, lapply(names(fisherB13c), function(set_name) {
  read.csv(file.path(output_dir, paste0(set_name, "_summary.csv")))
}))

# Save the combined summary
write.csv(combined_summary, file = paste0(output_dir, "/combined_summary.csv"), row.names = FALSE)

# Print the combined summary
print(combined_summary)


