#ifndef NETWORK_SIMPLEX_ADAPTER_H
#define NETWORK_SIMPLEX_ADAPTER_H

#include <RcppArmadillo.h>
#include <vector>
#include <algorithm>
#include <cmath>

// Forward declare the LEMON network simplex class
namespace lemon {
  template<typename GR, typename V, typename C>
  class NetworkSimplex;
  
  class FullBipartiteDigraph;
}

// Adapter class for nbonneel/network_simplex library
class NetworkSimplexAdapter {
private:
  // Padding strategy for odd-sized problems
  static bool needsPadding(int n, int m) {
    return (n % 2 == 1) || (m % 2 == 1);
  }
  
  // Pad marginals to even size by adding dummy node with zero mass
  static void padMarginals(arma::vec& a_padded, arma::vec& b_padded,
                          const arma::vec& a, const arma::vec& b) {
    int n = a.n_elem;
    int m = b.n_elem;
    
    // Determine padding needed
    int n_pad = (n % 2 == 1) ? n + 1 : n;
    int m_pad = (m % 2 == 1) ? m + 1 : m;
    
    // Initialize padded vectors
    a_padded = arma::vec(n_pad, arma::fill::zeros);
    b_padded = arma::vec(m_pad, arma::fill::zeros);
    
    // Copy original values
    a_padded.head(n) = a;
    b_padded.head(m) = b;
    
    // Add tiny mass to dummy nodes to maintain feasibility
    double eps = 1e-12;
    if (n_pad > n) a_padded(n) = eps;
    if (m_pad > m) b_padded(m) = eps;
    
    // Adjust marginals to maintain balance
    // Use a simpler approach: ensure total sums match
    double sum_a = arma::sum(a_padded);
    double sum_b = arma::sum(b_padded);
    
    if (std::abs(sum_a - sum_b) > 1e-15) {
      // Add difference to dummy node if available, otherwise to first element
      if (sum_a > sum_b) {
        if (m_pad > m) {
          b_padded(m) += (sum_a - sum_b);  // Add to dummy sink
        } else {
          b_padded(0) += (sum_a - sum_b);  // Add to first element
        }
      } else {
        if (n_pad > n) {
          a_padded(n) += (sum_b - sum_a);  // Add to dummy source
        } else {
          a_padded(0) += (sum_b - sum_a);  // Add to first element
        }
      }
    }
  }
  
  // Pad cost matrix with high costs for dummy nodes
  static arma::mat padCostMatrix(const arma::mat& cost, int n_pad, int m_pad) {
    int n = cost.n_rows;
    int m = cost.n_cols;
    
    // Create padded matrix with high default cost
    double high_cost = cost.max() * 1000;
    arma::mat cost_padded(n_pad, m_pad);
    cost_padded.fill(high_cost);
    
    // Copy original costs
    cost_padded.submat(0, 0, n-1, m-1) = cost;
    
    return cost_padded;
  }
  
  // Extract solution removing padding
  static arma::mat extractSolution(const arma::mat& P_padded, int n, int m) {
    return P_padded.submat(0, 0, n-1, m-1);
  }

public:
  // Main solver interface
  static arma::mat solve(const arma::mat& cost, 
                        const arma::vec& a, 
                        const arma::vec& b,
                        double eps = 1e-9);
  
  // Check if LEMON library is available at compile time
  static bool isAvailable();
};

// Fallback implementation using our existing network simplex
// (to be used if LEMON is not available)
namespace fallback {
  arma::mat network_simplex_ot(const arma::mat& cost, 
                               const arma::vec& a, 
                               const arma::vec& b,
                               double eps = 1e-9);
}

#endif // NETWORK_SIMPLEX_ADAPTER_H