# Benchmarking script for graphSmooth: R vs Rust backend
library(hySpc.hpc)
library(ggplot2)

# Ensure hyperSpec is available
if (!requireNamespace("hyperSpec", quietly = TRUE)) {
  stop("hyperSpec package is required but not installed.")
}
library(hyperSpec)

# Function to create mock hyperspectral dataset
create_mock_data <- function(width, height, bands) {
  n_pixels <- width * height
  # Create random data: pixels in rows, bands in columns
  spc <- matrix(rnorm(n_pixels * bands), nrow = n_pixels, ncol = bands)
  x <- new("hyperSpec", spc = spc)
  return(x)
}

# Define dataset sizes to benchmark: c(width, height, bands)
sizes <- list(
  c(20, 20, 50),
  c(40, 40, 50),
  c(60, 60, 50),
  c(80, 80, 50),
  c(100, 100, 50),
  c(150, 150, 50),
  c(200, 200, 50),
  c(400, 400, 50),
  c(800, 800, 50),
  c(1000, 1000, 50)
)

results <- data.frame(Width = integer(), Height = integer(), Bands = integer(),
                      TotalPixels = integer(), Backend = character(), Time_Seconds = numeric())

cat("Starting benchmark...\n")

for (s in sizes) {
  w <- s[1]
  h <- s[2]
  b <- s[3]
  n_pixels <- w * h
  
  cat(sprintf("\nBenchmarking size: %dx%d (pixels: %d), bands: %d\n", w, h, n_pixels, b))
  x <- create_mock_data(w, h, b)
  
  # Measure Rust time
  # Run once to warmup (if applicable)
  invisible(graphSmooth(x, width = w, height = h, alpha = 1.0, neighbors = 4, backend = "rust"))
  
  t_rust <- system.time({
    graphSmooth(x, width = w, height = h, alpha = 1.0, neighbors = 4, backend = "rust")
  })["elapsed"]
  cat(sprintf("  Rust time: %.3f seconds\n", t_rust))
  
  # Measure R time
  # Run once to warmup (if applicable)
  invisible(graphSmooth(x, width = w, height = h, alpha = 1.0, neighbors = 4, backend = "r"))
  
  t_r <- system.time({
    graphSmooth(x, width = w, height = h, alpha = 1.0, neighbors = 4, backend = "r")
  })["elapsed"]
  cat(sprintf("  R time:    %.3f seconds\n", t_r))
  
  results <- rbind(results, data.frame(
    Width = w, Height = h, Bands = b, TotalPixels = n_pixels, Backend = "Rust", Time_Seconds = t_rust
  ))
  
  results <- rbind(results, data.frame(
    Width = w, Height = h, Bands = b, TotalPixels = n_pixels, Backend = "R", Time_Seconds = t_r
  ))
}

# Save results
write.csv(results, "benchmarking/benchmark_results.csv", row.names = FALSE)

# Calculate speedup
rust_res <- subset(results, Backend == "Rust")
r_res <- subset(results, Backend == "R")

speedup_df <- data.frame(
  TotalPixels = rust_res$TotalPixels,
  Width = rust_res$Width,
  Height = rust_res$Height,
  Speedup = r_res$Time_Seconds / rust_res$Time_Seconds
)
write.csv(speedup_df, "benchmarking/speedup_results.csv", row.names = FALSE)

cat("\nBenchmark complete. Generating plots...\n")

# Plotting times
p1 <- ggplot(results, aes(x = TotalPixels, y = Time_Seconds, color = Backend)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  labs(title = "Benchmark: Rust vs R Spatial Smoothing",
       x = "Total Pixels (width * height)",
       y = "Time (seconds)",
       subtitle = paste("Bands =", sizes[[1]][3])) +
  theme_minimal() +
  scale_color_manual(values = c("R" = "#F8766D", "Rust" = "#00BFC4"))

# Plotting speedup
p2 <- ggplot(speedup_df, aes(x = TotalPixels, y = Speedup)) +
  geom_line(linewidth = 1, color = "#00BA38") +
  geom_point(size = 3, color = "#00BA38") +
  labs(title = "Speedup Achieved (R time / Rust time)",
       x = "Total Pixels (width * height)",
       y = "Speedup Factor",
       subtitle = paste("Bands =", sizes[[1]][3])) +
  theme_minimal() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red")

# Save plots
ggsave("benchmarking/benchmark_time_plot.png", plot = p1, width = 8, height = 6)
ggsave("benchmarking/benchmark_speedup_plot.png", plot = p2, width = 8, height = 6)

cat("Plots saved in benchmarking/ directory.\n")
