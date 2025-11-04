.onLoad <- function(libname, pkgname) {
  ns <- asNamespace(pkgname)
  cb <- function(name, fun) {
    unlockBinding(name, ns)
    assign(name, fun, envir = ns)
    lockBinding(name, ns)
  }
  cb("compute_edge_gradient_cpp", compute_edge_gradient_cpp_fallback)
  cb("compute_squared_distances_cpp", compute_squared_distances_cpp_fallback)
  cb("compute_rwr_vectorized_cpp", compute_rwr_vectorized_cpp_fallback)
  cb("solve_sinkhorn_stabilized_cpp", solve_sinkhorn_stabilized_cpp_fallback)
  cb("solve_sylvester_rwr_cpp", solve_sylvester_rwr_cpp_fallback)
  cb("compute_parrot_cost_cpp", compute_parrot_cost_cpp_fallback)
}
