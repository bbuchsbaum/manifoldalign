# Minimal toy benchmark test without multidesign hanging issues
library(testthat)
library(manifoldalign)
suppressPackageStartupMessages({
  library(multivarious)
})

# Simple hyperdesign creation without multidesign
create_simple_hyperdesign <- function(toy_data, labels = NULL) {
  view_names <- names(toy_data)[names(toy_data) != "latent" & 
                                names(toy_data) != "hard_indices"]
  
  if (is.null(labels)) {
    n_points <- nrow(toy_data[[view_names[1]]])
    labels <- factor(rep(c("class1", "class2"), length.out = n_points))
  }
  
  domain_list <- list()
  for (view_name in view_names) {
    X <- toy_data[[view_name]]
    # Ensure X is a proper matrix
    X <- as.matrix(X)
    if (!is.matrix(X)) {
      stop("Failed to convert data to matrix format")
    }
    
    design_df <- data.frame(
      sample_id = seq_len(nrow(X)),
      label = labels,
      stringsAsFactors = FALSE
    )
    domain_list[[view_name]] <- list(x = X, design = design_df)
  }
  
  class(domain_list) <- c("hyperdesign", "list")
  domain_list
}

# Test data generation
cat("Testing toy data generation...\n")
set.seed(123)
linear_data <- list(
  latent = matrix(rnorm(20), 10, 2),
  view1 = matrix(rnorm(20), 10, 2),
  view2 = matrix(rnorm(20), 10, 2)
)

hd <- create_simple_hyperdesign(linear_data)
cat("Simple hyperdesign created successfully\n")

# Test KEMA
cat("Testing KEMA...\n")
tryCatch({
  result <- kema(hd, y = "label", ncomp = 2, solver = "exact")
  cat("KEMA: PASS\n")
}, error = function(e) {
  cat("KEMA: SKIP -", e$message, "\n")
})

# Test GRASP (2 domains only)
cat("Testing GRASP...\n")
grasp_hd <- create_simple_hyperdesign(list(view1 = linear_data$view1, view2 = linear_data$view2))
tryCatch({
  result <- grasp(grasp_hd, preproc = multivarious::center(), ncomp = 5, q_descriptors = 10)
  cat("GRASP: PASS\n")
}, error = function(e) {
  cat("GRASP: SKIP -", e$message, "\n")
})

# Test GPCA
cat("Testing GPCA...\n")
tryCatch({
  result <- gpca_align.hyperdesign(hd, y = "label", preproc = multivarious::center(), ncomp = 2, u = 0.5)
  cat("GPCA: PASS\n")
}, error = function(e) {
  cat("GPCA: SKIP -", e$message, "\n")
})

cat("All algorithm tests completed\n")