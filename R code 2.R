# Load necessary libraries
library(circular)
library(CircStats)
library(ggplot2)
library(tibble)
library(dplyr)

# Create a directory to save all outputs
output_dir <- "termite_orientation1_analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Load the dataset
data(fisherB13c)

# Define a helper function to save plots and tables
save_output <- function(plot, filename) {
  ggsave(paste0(output_dir, "/", filename), plot = plot, width = 8, height = 6)
}

# Function to compute von Mises density
von_mises_density <- function(theta, mu, kappa) {
  dvonmises(theta, mu, kappa)
}

# Function to analyze each dataset in fisherB13c
analyze_circular_data <- function(dataset, set_name) {
  
  # Convert dataset to radians for von Mises distribution
  data_rad <- circular(dataset, units = "degrees")
  
  # Circular summary statistics
  mean_dir <- mean.circular(data_rad)
  circ_var <- var.circular(data_rad)
  
  # Fit a von Mises distribution to the data
  von_mises_fit <- mle.vonmises(data_rad)
  mu <- von_mises_fit$mu
  kappa <- von_mises_fit$kappa
  
  # Rayleigh test for uniformity
  rayleigh_test <- rayleigh.test(data_rad)
  
  # Generate the von Mises density curve based on fitted mu and kappa
  theta_seq <- seq(0, 2 * pi, length.out = 360)  # Sequence of angles in radians
  density_values <- von_mises_density(theta_seq, mu, kappa)  # Compute density
  
  # Create a data frame for the density plot
  density_df <- data.frame(theta = theta_seq, density = density_values)
  
  # Plot the circular data with von Mises density overlaid
  plot_title <- paste("Circular Plot with Von Mises Fit for", set_name)
  plot_data <- ggplot() +
    geom_histogram(aes(x = as.numeric(data_rad)), 
                   binwidth = 10, fill = "lightblue", color = "black", 
                   position = "identity", alpha = 0.6) +
    coord_polar(theta = "x") +
    geom_line(data = density_df, aes(x = theta * 180 / pi, y = density * max(table(data_rad)) / max(density)), 
              color = "red", size = 1) +  # Overlay von Mises density
    ggtitle(plot_title) +
    theme_minimal() +
    xlab("Angle (degrees)") +
    ylab("Density")
  
  # Save the plot
  save_output(plot_data, paste0(set_name, "_circular_plot.png"))
  
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
  
  # Print the results
  print(summary_data)
  
  # Save the summary table to a CSV
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
