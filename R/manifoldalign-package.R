#' manifoldalign package
#'
#' Internal package documentation and dynamic library registration.
#'
#' @return NULL
#' @keywords internal
#' @useDynLib manifoldalign, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom graphics axis hist
#' @importFrom grDevices heat.colors
#' @importFrom stats aggregate ave coef cov IQR lm median model.frame model.matrix model.response prcomp predict runif sd setNames
#' @importFrom utils capture.output modifyList setTxtProgressBar txtProgressBar write.csv
"_PACKAGE"

utils::globalVariables(c("anchors", "benchmark_label", "i", "j", "k", "v", "value"))
