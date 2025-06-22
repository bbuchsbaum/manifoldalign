# Example: Using FPGW for heterogeneous domain alignment
library(tibble)  # Required by multidesign
library(manifoldalign)
library(multidesign)

# Set up the problem: aligning visual and textual representations
set.seed(42)
n <- 50

# Visual features (3D)
visual_features <- matrix(rnorm(n * 3), n, 3)
visual_features[1:25, ] <- visual_features[1:25, ] + 2  # Class 1
visual_features[26:50, ] <- visual_features[26:50, ] - 2  # Class 2

# Textual features (10D) 
text_features <- matrix(rnorm(n * 10), n, 10)
text_features[1:25, ] <- text_features[1:25, ] + 1
text_features[26:50, ] <- text_features[26:50, ] - 1

# Add noise to text domain
noise_samples <- 10
text_features <- rbind(text_features, 
                      matrix(rnorm(noise_samples * 10, sd = 3), noise_samples, 10))

# Create multidesign objects
design_visual <- data.frame(
  id = 1:n,
  class = factor(rep(c("A", "B"), each = 25))
)

design_text <- data.frame(
  id = 1:(n + noise_samples),
  class = factor(c(rep(c("A", "B"), each = 25), rep("noise", noise_samples)))
)

md_visual <- multidesign(visual_features, design_visual)
md_text <- multidesign(text_features, design_text)

# Create hyperdesign
hd <- hyperdesign(list(
  visual = md_visual,
  text = md_text
))

# 1. Classical FGW alignment
cat("1. Classical Fused Gromov-Wasserstein\n")
cat("=====================================\n")
result_classical <- fpgw(hd, omega1 = 0.3, verbose = TRUE)
print(result_classical)

# 2. Mass-constrained to handle noise
cat("\n2. Mass-constrained FPGW (ignore noise)\n")
cat("========================================\n")
# Transport only 80% of mass to avoid noise samples
result_partial <- fpgw(hd, omega1 = 0.3, rho = 0.8, verbose = TRUE)

# Compare transported mass
P_classical <- result_classical$transport_plans[[1]]
P_partial <- result_partial$transport_plans[[1]]
cat("\nTransported mass comparison:\n")
cat("  Classical FGW:", sum(P_classical), "\n")
cat("  Mass-constrained:", sum(P_partial), "\n")

# 3. Analyze transport plan
cat("\n3. Transport Plan Analysis\n")
cat("==========================\n")

# Which text samples are matched to visual?
matched_text <- colSums(P_partial) > 0.01
cat("Text samples with matches:", sum(matched_text), "out of", ncol(P_partial), "\n")

# Are noise samples excluded?
noise_indices <- (n + 1):(n + noise_samples)
noise_mass <- sum(colSums(P_partial)[noise_indices])
cat("Mass assigned to noise samples:", noise_mass, "\n")

# 4. Use S3 methods
cat("\n4. Using S3 Methods\n")
cat("===================\n")

# Summary
summary(result_partial)

# Extract coefficients
params <- coef(result_partial, type = "parameters")
cat("\nExtracted parameters:\n")
print(params)

# Plot distance matrix (in real use, would create actual plot)
# plot(result_partial, which = 1)

# 5. Multi-domain example
cat("\n5. Multi-domain Alignment\n")
cat("=========================\n")

# Add a third domain (audio features, 7D)
audio_features <- matrix(rnorm(n * 7), n, 7)
audio_features[1:25, ] <- audio_features[1:25, ] + 1.5

md_audio <- multidesign(audio_features, design_visual)

# Create a subset of text data without noise
md_text_clean <- multidesign(text_features[1:n, ], design_visual)

hd_multi <- hyperdesign(list(
  visual = md_visual,
  text = md_text_clean,
  audio = md_audio
))

result_multi <- fpgw(hd_multi, omega1 = 0.4, max_iter = 30)
cat("\nPairwise distances between modalities:\n")
print(round(result_multi$distances, 3))

# 6. Different omega values
cat("\n6. Effect of omega1 Parameter\n")
cat("=============================\n")

omega_values <- c(0.01, 0.5, 0.99)
distances <- numeric(length(omega_values))

for (i in seq_along(omega_values)) {
  res <- fpgw(hd, omega1 = omega_values[i], verbose = FALSE)
  distances[i] <- res$distances[1, 2]
}

cat("\nDistance vs omega1:\n")
for (i in seq_along(omega_values)) {
  cat(sprintf("  omega1 = %.2f: distance = %.4f\n", omega_values[i], distances[i]))
}

cat("\nExample complete!\n")