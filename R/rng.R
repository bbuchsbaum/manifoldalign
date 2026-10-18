# Seed the RNG for the rest of the calling function and restore the caller's
# `.Random.seed` (or its absence) when that function exits, including on
# error. Use inside package functions instead of set.seed(); calling it more
# than once in the same frame still restores the state from before the first
# call, because exit handlers run in reverse order.
.local_seed <- function(seed, frame = parent.frame()) {
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv()) else NULL
  restore <- function() {
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }
  # after = FALSE prepends, so handlers run last-in-first-out.
  do.call(base::on.exit, list(as.call(list(restore)), TRUE, FALSE),
          envir = frame)
  set.seed(seed)
  invisible(NULL)
}
