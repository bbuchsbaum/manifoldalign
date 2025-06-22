#include <Rcpp.h>
#include "full_bipartitegraph.h"
#include <iostream>

// [[Rcpp::export]]
void test_lemon_arc_indexing() {
  using namespace lemon;
  
  int n1 = 3;  // sources
  int n2 = 2;  // targets
  
  FullBipartiteDigraph graph(n1, n2);
  
  Rcpp::Rcout << "Testing LEMON arc indexing:\n";
  Rcpp::Rcout << "Sources: " << n1 << ", Targets: " << n2 << "\n";
  Rcpp::Rcout << "Total nodes: " << graph.nodeNum() << "\n";
  Rcpp::Rcout << "Total arcs: " << graph.arcNum() << "\n\n";
  
  // Test arc creation and indexing
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      // Create arc from source i to target (n1 + j)
      FullBipartiteDigraph::Arc arc = graph.arc(i, n1 + j);
      int64_t arc_id = FullBipartiteDigraph::id(arc);
      
      // Get arc from ID
      FullBipartiteDigraph::Arc arc_from_id = graph.arcFromId(arc_id);
      
      // Get source and target
      int source = graph.source(arc_from_id);
      int target = graph.target(arc_from_id);
      
      Rcpp::Rcout << "Arc from " << i << " to " << (n1 + j) 
                  << ": ID = " << arc_id 
                  << ", retrieved source = " << source
                  << ", retrieved target = " << target << "\n";
                  
      // Verify the formula
      int64_t expected_id = (int64_t)i * (int64_t)n2 + (int64_t)j;
      if (arc_id != expected_id) {
        Rcpp::Rcout << "  WARNING: Expected ID " << expected_id 
                    << " but got " << arc_id << "\n";
      }
    }
  }
  
  // Test sequential arc IDs
  Rcpp::Rcout << "\nTesting sequential arc IDs:\n";
  for (int64_t id = 0; id < graph.arcNum(); id++) {
    FullBipartiteDigraph::Arc arc = graph.arcFromId(id);
    int source = graph.source(arc);
    int target = graph.target(arc);
    Rcpp::Rcout << "Arc ID " << id << ": " << source << " -> " << target << "\n";
  }
}