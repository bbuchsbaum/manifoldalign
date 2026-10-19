#' manifoldalign package
#'
#' Internal package documentation and dynamic library registration.
#'
#' @return NULL
#' @keywords internal
#' @useDynLib manifoldalign, .registration = TRUE
#' @importFrom graphics axis hist
#' @importFrom grDevices heat.colors
#' @importFrom stats ave coef sd setNames aggregate cov IQR lm median model.frame model.matrix model.response prcomp predict runif
#' @importFrom utils capture.output modifyList setTxtProgressBar txtProgressBar write.csv
"_PACKAGE"

utils::globalVariables(c("benchmark_label", "anchors", "i", "j", "k", "v", "value"))
