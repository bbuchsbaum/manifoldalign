#include "tv_proximal.h"
#include <algorithm>
#include <vector>

// Pool Adjacent Violators Algorithm (PAVA) for isotonic regression
// Solves: argmin_x { sum_i w_i * (x_i - y_i)^2 } subject to x_1 <= x_2 <= ... <= x_n
void TVProximal::pava_increasing(arma::vec& y, const arma::vec& w) {
  int n = y.n_elem;
  if (n <= 1) return;
  
  // Pool-adjacent-violators algorithm
  std::vector<double> pooled_y(n);
  std::vector<double> pooled_w(n);
  std::vector<int> pool_start(n);
  std::vector<int> pool_end(n);
  
  // Initialize pools
  for (int i = 0; i < n; i++) {
    pooled_y[i] = y[i];
    pooled_w[i] = w[i];
    pool_start[i] = i;
    pool_end[i] = i;
  }
  
  int num_pools = n;
  
  // Main PAVA loop
  bool changed = true;
  while (changed) {
    changed = false;
    
    int write_idx = 0;
    int i = 0;
    
    while (i < num_pools) {
      // Check if we should merge with next pool
      if (i < num_pools - 1 && pooled_y[i] > pooled_y[i + 1]) {
        // Merge pools i and i+1
        double total_weight = pooled_w[i] + pooled_w[i + 1];
        double weighted_sum = pooled_w[i] * pooled_y[i] + pooled_w[i + 1] * pooled_y[i + 1];
        
        pooled_y[write_idx] = weighted_sum / total_weight;
        pooled_w[write_idx] = total_weight;
        pool_start[write_idx] = pool_start[i];
        pool_end[write_idx] = pool_end[i + 1];
        
        i += 2;  // Skip the merged pool
        changed = true;
      } else {
        // Keep current pool
        if (write_idx != i) {
          pooled_y[write_idx] = pooled_y[i];
          pooled_w[write_idx] = pooled_w[i];
          pool_start[write_idx] = pool_start[i];
          pool_end[write_idx] = pool_end[i];
        }
        i++;
      }
      write_idx++;
    }
    
    num_pools = write_idx;
  }
  
  // Assign pooled values back to y
  for (int p = 0; p < num_pools; p++) {
    for (int i = pool_start[p]; i <= pool_end[p]; i++) {
      y[i] = pooled_y[p];
    }
  }
}

// 1D TV proximal operator using Condat's algorithm
// Reference: "A Direct Algorithm for 1D Total Variation Denoising", L. Condat, 2013
arma::vec TVProximal::prox_tv_1d(const arma::vec& y, double lambda) {
  int n = y.n_elem;
  if (n <= 1 || lambda <= 0) return y;
  
  arma::vec x = y;  // Copy input
  
  // Condat's algorithm for TV denoising
  // This is more efficient than PAVA for general TV
  int k = 0;  // Current sample
  int k0 = 0; // Beginning of current segment
  int km = 0; // Last position where umax was attained
  int kp = 0; // First position where umin was attained
  
  double vmin = y[0] - lambda;
  double vmax = y[0] + lambda;
  double umin = lambda;
  double umax = -lambda;
  
  for (;;) {
    // Check if we can extend current segment
    while (k < n - 1) {
      if (y[k + 1] + umin < vmin - lambda) {
        // Start new segment
        for (int i = k0; i <= km; i++) {
          x[i] = vmin;
        }
        k = km + 1;
        k0 = k;
        km = k;
        kp = k;
        vmin = y[k];
        vmax = y[k] + 2 * lambda;
        umin = lambda;
        umax = -lambda;
      } else if (y[k + 1] + umax > vmax + lambda) {
        // Start new segment
        for (int i = k0; i <= kp; i++) {
          x[i] = vmax;
        }
        k = kp + 1;
        k0 = k;
        km = k;
        kp = k;
        vmin = y[k] - 2 * lambda;
        vmax = y[k];
        umin = lambda;
        umax = -lambda;
      } else {
        // Continue current segment
        k++;
        umin += y[k] - vmin;
        umax += y[k] - vmax;
        
        if (umin >= lambda) {
          vmin += (umin - lambda) / (k - k0 + 1);
          umin = lambda;
          km = k;
        }
        if (umax <= -lambda) {
          vmax += (umax + lambda) / (k - k0 + 1);
          umax = -lambda;
          kp = k;
        }
      }
    }
    
    // Process last segment
    if (umin < 0) {
      for (int i = k0; i <= km; i++) {
        x[i] = vmin;
      }
      k = km + 1;
      k0 = k;
      km = k;
      vmin = y[k];
      umin = lambda;
      umax = y[k] + lambda - vmax;
    } else if (umax > 0) {
      for (int i = k0; i <= kp; i++) {
        x[i] = vmax;
      }
      k = kp + 1;
      k0 = k;
      kp = k;
      vmax = y[k];
      umax = -lambda;
      umin = y[k] - lambda - vmin;
    } else {
      for (int i = k0; i <= n - 1; i++) {
        x[i] = vmin + umin / (n - k0);
      }
      break;
    }
    
    if (k >= n - 1) {
      if (k0 < n) {
        x[k0] = vmin + umin;
      }
      break;
    }
  }
  
  return x;
}

// Apply TV proximal operator to rows of a matrix
arma::mat TVProximal::prox_tv_rows(const arma::mat& Y, double lambda) {
  arma::mat X = Y;
  for (arma::uword i = 0; i < Y.n_rows; i++) {
    arma::vec row = Y.row(i).t();
    X.row(i) = prox_tv_1d(row, lambda).t();
  }
  return X;
}

// Apply TV proximal operator to columns of a matrix
arma::mat TVProximal::prox_tv_cols(const arma::mat& Y, double lambda) {
  arma::mat X = Y;
  for (arma::uword j = 0; j < Y.n_cols; j++) {
    X.col(j) = prox_tv_1d(Y.col(j), lambda);
  }
  return X;
}

// Apply 2D TV proximal operator (alternating directions)
arma::mat TVProximal::prox_tv_2d(const arma::mat& Y, double lambda, int max_iter) {
  arma::mat X = Y;
  
  for (int iter = 0; iter < max_iter; iter++) {
    // Apply to rows
    X = prox_tv_rows(X, lambda);
    
    // Apply to columns
    X = prox_tv_cols(X, lambda);
    
    // Could add convergence check here
  }
  
  return X;
}

// Compute TV subgradient
arma::vec TVProximal::tv_subgradient(const arma::vec& x) {
  int n = x.n_elem;
  arma::vec g(n, arma::fill::zeros);
  
  // TV subgradient: g[i] = sign(x[i] - x[i-1]) - sign(x[i+1] - x[i])
  for (int i = 0; i < n; i++) {
    if (i > 0) {
      double diff = x[i] - x[i-1];
      if (std::abs(diff) > 1e-10) {
        g[i] += (diff > 0) ? 1.0 : -1.0;
      }
    }
    if (i < n - 1) {
      double diff = x[i+1] - x[i];
      if (std::abs(diff) > 1e-10) {
        g[i] -= (diff > 0) ? 1.0 : -1.0;
      }
    }
  }
  
  return g;
}

// Export to R
// [[Rcpp::export]]
arma::mat apply_tv_proximal_cpp(const arma::mat& Y, double lambda) {
  return TVProximal::prox_tv_2d(Y, lambda, 5);  // 5 iterations of alternating directions
}