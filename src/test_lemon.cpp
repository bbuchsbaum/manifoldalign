// Test if LEMON headers are available
#include <RcppArmadillo.h>

// Try to include LEMON headers
#ifdef __has_include
  #if __has_include("network_simplex_simple.h")
    #define HAVE_LEMON_NETWORK_SIMPLEX 1
    #include "network_simplex_simple.h"
    #include "full_bipartitegraph.h"
  #else
    #define HAVE_LEMON_NETWORK_SIMPLEX 0
  #endif
#else
  // Older compilers - try to include anyway
  #include "network_simplex_simple.h"
  #include "full_bipartitegraph.h"
  #define HAVE_LEMON_NETWORK_SIMPLEX 1
#endif

// [[Rcpp::export]]
bool has_lemon_network_simplex() {
#if HAVE_LEMON_NETWORK_SIMPLEX
  return true;
#else
  return false;
#endif
}