// Standalone serial build plus independent exhaustive 2x2 transport oracle.
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include "network_simplex_simple.h"

int main() {
#ifdef _OPENMP
#error This test must compile without OpenMP
#endif
  using Graph = lemon::FullBipartiteDigraph;
  using Solver = lemon::NetworkSimplexSimple<Graph, int64_t, int64_t>;
  int cases = 0;
  for (int code = 0; code < 256; ++code) {
    int64_t cost[4];
    int digits = code;
    for (int j = 0; j < 4; ++j) { cost[j] = digits % 4 - 1; digits /= 4; }
    for (int first = 1; first <= 2; ++first) {
      int64_t a[2] = {first, 3 - first}, b[2] = {-2, -1};
      Graph graph(2, 2);
      Solver solver(graph, true, 4, 4);
      for (int j = 0; j < 4; ++j) solver.setCost(graph.arcFromId(j), cost[j]);
      solver.supplyMap(a, 2, b, 2);
      assert(solver.run() == Solver::OPTIMAL);
      int64_t flow[4], value = 0;
      for (int j = 0; j < 4; ++j) {
        flow[j] = solver.flow(graph.arcFromId(j));
        assert(flow[j] >= 0);
        value += flow[j] * cost[j];
      }
      assert(flow[0] + flow[1] == a[0] && flow[2] + flow[3] == a[1]);
      assert(flow[0] + flow[2] == 2 && flow[1] + flow[3] == 1);
      int64_t oracle = INT64_MAX;
      for (int64_t x = 0; x <= a[0]; ++x) {
        int64_t y = a[0] - x, z = 2 - x, w = a[1] - z;
        if (z >= 0 && w >= 0 && y + w == 1) {
          oracle = std::min(oracle, x*cost[0] + y*cost[1] + z*cost[2] + w*cost[3]);
        }
      }
      assert(value == oracle);
      ++cases;
    }
  }
  std::cout << cases << " serial transport oracles passed\n";
}
