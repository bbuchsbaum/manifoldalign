#include "network_simplex_adapter.h"
#include "ot_common.h"
#include <vector>
#include <cstdint>
#include <cmath>
#include <cassert>
#include <limits>

// Check if we have the LEMON headers available
#ifdef __has_include
  #if __has_include("network_simplex_simple.h") && __has_include("full_bipartitegraph.h")
    #define LEMON_AVAILABLE 1
    #include "network_simplex_simple.h"
    #include "full_bipartitegraph.h"
  #else
    #define LEMON_AVAILABLE 0
  #endif
#else
  // Try to include anyway for older compilers
  #define LEMON_AVAILABLE 1
  #include "network_simplex_simple.h"
  #include "full_bipartitegraph.h"
#endif

// Implementation of the adapter
arma::mat NetworkSimplexAdapter::solve(const arma::mat& cost, 
                                       const arma::vec& a, 
                                       const arma::vec& b,
                                       double eps) {
  
  // Check marginals first
  check_marginals(a, b);
  
  int n_orig = a.n_elem;
  int m_orig = b.n_elem;
  
#if LEMON_AVAILABLE
  // Use LEMON implementation with padding for odd sizes
  
  // Check if we need padding
  if (needsPadding(n_orig, m_orig)) {
    // Pad to even sizes to avoid performance issues
    arma::vec a_padded, b_padded;
    padMarginals(a_padded, b_padded, a, b);
    
    int n_pad = a_padded.n_elem;
    int m_pad = b_padded.n_elem;
    
    // Pad cost matrix
    arma::mat cost_padded = padCostMatrix(cost, n_pad, m_pad);
    
    // Create bipartite graph
    typedef lemon::FullBipartiteDigraph Digraph;
    Digraph graph(n_pad, m_pad);
    
    // Set up types - use int64_t for costs and flows to avoid overflow
    typedef int64_t FlowType;
    typedef int64_t CostType;
    
    // Create the network simplex solver
    typedef lemon::NetworkSimplexSimple<Digraph, FlowType, CostType> NetworkSimplex;
    
    // Scale to integers for precision
    const double scale = 1e6;
    
    // Create cost and supply vectors
    std::vector<CostType> costs(n_pad * m_pad);
    std::vector<FlowType> supplies_src(n_pad);
    std::vector<FlowType> supplies_dst(m_pad);
    
    // Set costs
    size_t arc_id = 0;
    for (int i = 0; i < n_pad; i++) {
      for (int j = 0; j < m_pad; j++) {
        assert(arc_id < costs.size());
        costs[arc_id] = static_cast<CostType>(std::round(cost_padded(i, j) * scale));
        arc_id++;
      }
    }
    
    // Set supplies - sources positive, sinks negative
    FlowType total_supply = 0;
    FlowType total_demand = 0;
    for (int i = 0; i < n_pad; i++) {
      supplies_src[i] = static_cast<FlowType>(std::round(a_padded(i) * scale));
      total_supply += supplies_src[i];
    }
    for (int j = 0; j < m_pad; j++) {
      supplies_dst[j] = -static_cast<FlowType>(std::round(b_padded(j) * scale));
      total_demand -= supplies_dst[j];  // Note: subtracting negative = adding
    }
    
    // Debug: check mass balance
    if (total_supply != total_demand) {
      Rcpp::Rcout << "Mass imbalance detected: supply = " << total_supply 
                  << ", demand = " << total_demand 
                  << ", difference = " << (total_supply - total_demand) << std::endl;
      
      // Fix imbalance by adjusting last supply
      if (n_pad > 0) {
        supplies_src[n_pad - 1] -= (total_supply - total_demand);
      }
    }
    
    // Create and configure solver
    NetworkSimplex ns(graph, true, n_pad + m_pad, n_pad * m_pad);
    
    // Set all arc costs using LEMON's arc indexing
    for (int i = 0; i < n_pad; i++) {
      for (int j = 0; j < m_pad; j++) {
        // LEMON uses arc_id = source * n2 + (target - n1)
        // where source is in [0, n1) and target is in [n1, n1+n2)
        int lemon_arc_id = i * m_pad + j;
        int cost_idx = i * m_pad + j;  // Our sequential index for costs array
        assert(cost_idx < costs.size());
        assert(lemon_arc_id <= std::numeric_limits<int>::max());
        
        // Debug: verify arc is valid before using it
        if (lemon_arc_id >= n_pad * m_pad) {
          Rcpp::Rcout << "ERROR: Arc ID " << lemon_arc_id << " out of bounds (max " 
                      << (n_pad * m_pad - 1) << ")" << std::endl;
          Rcpp::stop("Arc ID out of bounds");
        }
        
        try {
          ns.setCost(graph.arcFromId(lemon_arc_id), costs[cost_idx]);
        } catch (...) {
          Rcpp::Rcout << "ERROR setting cost for arc " << lemon_arc_id 
                      << " (i=" << i << ", j=" << j << ")" << std::endl;
          throw;
        }
      }
    }
    
    // Set supplies
    ns.supplyMap(&supplies_src[0], n_pad, &supplies_dst[0], m_pad);
    
    // Run solver
    int result = ns.run();
    if (result != NetworkSimplex::OPTIMAL) {
      Rcpp::Rcout << "Network simplex error code: " << result << std::endl;
      Rcpp::Rcout << "Problem size: " << n_pad << " x " << m_pad << std::endl;
      if (result == NetworkSimplex::INFEASIBLE) {
        Rcpp::stop("Network simplex: Problem is infeasible");
      } else if (result == NetworkSimplex::UNBOUNDED) {
        Rcpp::stop("Network simplex: Problem is unbounded");
      } else {
        Rcpp::stop("Network simplex: Unknown error");
      }
    }
    
    // Extract solution using LEMON's arc indexing
    arma::mat P_padded(n_pad, m_pad, arma::fill::zeros);
    for (int i = 0; i < n_pad; i++) {
      for (int j = 0; j < m_pad; j++) {
        int lemon_arc_id = i * m_pad + j;
        assert(lemon_arc_id <= std::numeric_limits<int>::max());
        FlowType flow = ns.flow(graph.arcFromId(lemon_arc_id));
        P_padded(i, j) = static_cast<double>(flow) / scale;
      }
    }
    
    // Remove padding and return
    return extractSolution(P_padded, n_orig, m_orig);
    
  } else {
    // Even sizes - no padding needed
    int n = n_orig;
    int m = m_orig;
    
    // Create bipartite graph
    typedef lemon::FullBipartiteDigraph Digraph;
    Digraph graph(n, m);
    
    // Set up types - use int64_t for costs and flows to avoid overflow
    typedef int64_t FlowType;
    typedef int64_t CostType;
    
    // Create the network simplex solver
    typedef lemon::NetworkSimplexSimple<Digraph, FlowType, CostType> NetworkSimplex;
    
    // Scale to integers for precision
    const double scale = 1e6;
    
    // Create cost and supply vectors
    std::vector<CostType> costs(n * m);
    std::vector<FlowType> supplies_src(n);
    std::vector<FlowType> supplies_dst(m);
    
    // Set costs
    int arc_id = 0;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        costs[arc_id] = static_cast<CostType>(std::round(cost(i, j) * scale));
        arc_id++;
      }
    }
    
    // Set supplies - sources positive, sinks negative
    FlowType total_supply = 0;
    FlowType total_demand = 0;
    for (int i = 0; i < n; i++) {
      supplies_src[i] = static_cast<FlowType>(std::round(a(i) * scale));
      total_supply += supplies_src[i];
    }
    for (int j = 0; j < m; j++) {
      supplies_dst[j] = -static_cast<FlowType>(std::round(b(j) * scale));
      total_demand -= supplies_dst[j];  // Note: subtracting negative = adding
    }
    
    // Debug: check mass balance
    if (total_supply != total_demand) {
      Rcpp::Rcout << "Mass imbalance detected: supply = " << total_supply 
                  << ", demand = " << total_demand 
                  << ", difference = " << (total_supply - total_demand) << std::endl;
      
      // Fix imbalance by adjusting last supply
      if (n > 0) {
        supplies_src[n - 1] -= (total_supply - total_demand);
      }
    }
    
    // Create and configure solver
    NetworkSimplex ns(graph, true, n + m, n * m);
    
    // Set all arc costs using LEMON's arc indexing
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        int lemon_arc_id = i * m + j;
        int cost_idx = i * m + j;
        ns.setCost(graph.arcFromId(lemon_arc_id), costs[cost_idx]);
      }
    }
    
    // Set supplies
    ns.supplyMap(&supplies_src[0], n, &supplies_dst[0], m);
    
    // Run solver
    int result = ns.run();
    if (result != NetworkSimplex::OPTIMAL) {
      Rcpp::Rcout << "Network simplex error code: " << result << std::endl;
      Rcpp::Rcout << "Problem size: " << n << " x " << m << std::endl;
      if (result == NetworkSimplex::INFEASIBLE) {
        Rcpp::stop("Network simplex: Problem is infeasible");
      } else if (result == NetworkSimplex::UNBOUNDED) {
        Rcpp::stop("Network simplex: Problem is unbounded");
      } else {
        Rcpp::stop("Network simplex: Unknown error");
      }
    }
    
    // Extract solution using LEMON's arc indexing
    arma::mat P(n, m, arma::fill::zeros);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        int lemon_arc_id = i * m + j;
        FlowType flow = ns.flow(graph.arcFromId(lemon_arc_id));
        P(i, j) = static_cast<double>(flow) / scale;
      }
    }
    
    return P;
  }
  
#else
  // Fallback to existing implementation with warning
  Rcpp::warning("LEMON network simplex not available, using fallback solver (may be less accurate)");
  return fallback::network_simplex_ot(cost, a, b, eps);
#endif
}

bool NetworkSimplexAdapter::isAvailable() {
#if LEMON_AVAILABLE
  return true;
#else
  return false;
#endif
}

// Fallback implementation using lpSolve via the R helper
namespace fallback {

arma::mat network_simplex_ot(const arma::mat& cost,
                              const arma::vec& a,
                              const arma::vec& b,
                              double /*eps*/) {
  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("manifoldalign");
  Rcpp::Function lp_fun = pkg["classical_ot_lp_r"];

  Rcpp::NumericMatrix cost_r = Rcpp::wrap(cost);
  Rcpp::NumericVector a_r = Rcpp::wrap(a);
  Rcpp::NumericVector b_r = Rcpp::wrap(b);

  Rcpp::NumericMatrix sol = lp_fun(cost_r, a_r, b_r);

  return Rcpp::as<arma::mat>(sol);
}

} // namespace fallback
