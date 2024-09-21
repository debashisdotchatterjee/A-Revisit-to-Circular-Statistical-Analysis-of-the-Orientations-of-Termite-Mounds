# Load necessary libraries
library(circular)
library(ggplot2)
library(tibble)
library(dplyr)

# Create a directory to save all outputs
output_dir <- "wrapped_normal_fit_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Load the dataset
data(fisherB13c)

# Define a helper function to save plots and tables
save_output <- function(plot, filename) {
  ggsave(paste0(output_dir, "/", filename), plot = plot, width = 8, height = 6)
}

# Function to fit Wrapped Normal distribution and calculate density
wrapped_normal_density <- function(theta, mu, rho) {
  # Wrapped normal distribution density function
  # mu: mean direction, rho: concentration parameter
  dwrappednormal(theta, mu, rho)
}

# Function to analyze each dataset in fisherB13c
analyze_circular_data_wrapped_normal <- function(dataset, set_name) {
  
  # Convert dataset to radians for the analysis
  data_rad <- circular(dataset, units = "degrees")
  
  # Estimate mean direction and standard deviation for the wrapped normal distribution
  mean_dir <- mean.circular(data_rad)
  circ_sd <- sd.circular(data_rad)
  
  # Convert standard deviation into rho (concentration parameter for wrapped normal)
  rho <- 1 / circ_sd^2
  
  # Generate the Wrapped Normal density curve based on estimated mean direction (mu) and rho
  theta_seq <- seq(0, 2 * pi, length.out = 360)  # Sequence of angles in radians
  density_values <- wrapped_normal_density(theta_seq, mean_dir, rho)  # Compute density
  
  # Scale the density values to match the height of the histogram
  hist_values <- as.numeric(table(cut(as.numeric(data_rad), breaks = 36)))  # 36 bins for rose plot
  max_hist <- max(hist_values)
  scaled_density <- density_values * max_hist / max(density_values)  # Scale to match histogram
  
  # Create a data frame for the density plot
  density_df <- data.frame(theta = theta_seq * 180 / pi, density = scaled_density)
  
  # Plot the circular data with Wrapped Normal distribution density overlaid
  plot_title <- paste("Circular Plot with Wrapped Normal Fit for", set_name)
  plot_data <- ggplot() +
    geom_histogram(aes(x = as.numeric(data_rad)), 
                   binwidth = 10, fill = "lightblue", color = "black", 
                   position = "identity", alpha = 0.6) +
    coord_polar(theta = "x") +
    geom_line(data = density_df, aes(x = theta, y = density), 
              color = "red", size = 1.2) +  # Overlay scaled Wrapped Normal density
    ggtitle(plot_title) +
    theme_minimal() +
    xlab("Angle (degrees)") +
    ylab("Density")
  
  # Save the plot with Wrapped Normal distribution overlay
  save_output(plot_data, paste0(set_name, "_circular_plot_wrapped_normal.png"))
  
  # Save plot without Wrapped Normal overlay for comparison
  plot_no_density <- ggplot() +
    geom_histogram(aes(x = as.numeric(data_rad)), 
                   binwidth = 10, fill = "lightblue", color = "black", 
                   position = "identity", alpha = 0.6) +
    coord_polar(theta = "x") +
    ggtitle(paste("Circular Plot for", set_name)) +
    theme_minimal() +
    xlab("Angle (degrees)") +
    ylab("Density")
  
  # Save the plot without Wrapped Normal fit
  save_output(plot_no_density, paste0(set_name, "_circular_plot_no_density.png"))
  
  # Summary output
  summary_data <- tibble(
    Set = set_name,
    Mean_Direction = as.numeric(mean_dir),
    Circular_Standard_Deviation = as.numeric(circ_sd),
    Wrapped_Normal_Concentration = rho
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
  analyze_circular_data_wrapped_normal(dataset, set_name)
})

# Combine all summaries into a single table
combined_summary <- do.call(rbind, all_summaries)

# Save the combined summary as a CSV file
write.csv(combined_summary, file = paste0(output_dir, "/combined_summary.csv"), row.names = FALSE)

# Print the combined summary
print(combined_summary)

# Generate a table of all the model fits and save it
combined_summary_plot <- ggplot(combined_summary, aes(x = Set, y = Wrapped_Normal_Concentration)) +
  geom_bar(stat = "identity", fill = "lightgreen", color = "black") +
  ggtitle("Wrapped Normal Concentration for each Set") +
  xlab("Set") + ylab("Concentration (rho)") +
  theme_minimal()

save_output(combined_summary_plot, "wrapped_normal_concentration_comparison.png")

# Print the combined plot
print(combined_summary_plot)
