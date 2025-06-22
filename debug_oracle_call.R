library(manifoldalign)

# Test the oracle directly with high costs
n <- 4
m <- 5
p <- rep(1/n, n)
q <- rep(1/m, m)

# Low cost
cost_low <- matrix(runif(n*m, 0, 1), n, m)
gamma_low <- manifoldalign:::partial_ot_tv(cost_low, p, q, lambda = 1)
cat("Low cost (0-1): mass =", sum(gamma_low), "\n")

# Medium cost  
cost_med <- matrix(runif(n*m, 10, 20), n, m)
gamma_med <- manifoldalign:::partial_ot_tv(cost_med, p, q, lambda = 1)
cat("Medium cost (10-20): mass =", sum(gamma_med), "\n")

# High cost
cost_high <- matrix(runif(n*m, 100, 200), n, m)
gamma_high <- manifoldalign:::partial_ot_tv(cost_high, p, q, lambda = 1)
cat("High cost (100-200): mass =", sum(gamma_high), "\n")

# Very high cost
cost_vhigh <- matrix(runif(n*m, 1000, 2000), n, m)
gamma_vhigh <- manifoldalign:::partial_ot_tv(cost_vhigh, p, q, lambda = 1)
cat("Very high cost (1000-2000): mass =", sum(gamma_vhigh), "\n")

# Check what happens with uniform high cost
cost_uniform <- matrix(1000, n, m)
gamma_uniform <- manifoldalign:::partial_ot_tv(cost_uniform, p, q, lambda = 1)
cat("Uniform high cost (all 1000): mass =", sum(gamma_uniform), "\n")