#' manifoldalign package
#'
#' Internal package documentation and dynamic library registration.
#'
#' @return NULL
#' @keywords internal
#' @useDynLib manifoldalign, .registration = TRUE
#' @importFrom graphics axis hist
#' @importFrom grDevices heat.colors
#' @importFrom stats aggregate cov IQR lm median model.frame model.matrix model.response prcomp predict runif
#' @importFrom utils modifyList setTxtProgressBar txtProgressBar write.csv
"_PACKAGE"

utils::globalVariables(c("anchors", "i", "j", "k", "v", "value"))
