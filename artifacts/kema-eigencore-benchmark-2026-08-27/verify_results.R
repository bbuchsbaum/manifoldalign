#!/usr/bin/env Rscript

root <- "artifacts/kema-eigencore-benchmark-2026-08-27"

toy <- utils::read.csv(file.path(root, "toy-feature", "raw_results.csv"))
toy_paired <- utils::read.csv(file.path(root, "comparison", "toy_paired.csv"))
toy_full <- utils::read.csv(file.path(root, "toy-feature-full", "raw_results.csv"))
tournament <- utils::read.csv(file.path(root, "comparison", "tournament_paired.csv"))

stopifnot(
  nrow(toy) == 45L,
  sum(is.na(toy$error)) == 36L,
  all(toy$fidelity_passed[is.na(toy$error)]),
  nrow(toy_paired) == 45L,
  all(toy_paired$old_ok == toy_paired$new_ok),
  nrow(toy_full) == 63L,
  sum(is.na(toy_full$error)) == 54L,
  all(toy_full$fidelity_passed[is.na(toy_full$error)]),
  identical(unique(toy_full$scenario[!is.na(toy_full$error)]), "partial_open_set"),
  nrow(tournament) == 42L,
  all(tournament$status_old == tournament$status_new)
)

tournament_ok <- tournament[
  tournament$status_old == "ok" & tournament$status_new == "ok",
  , drop = FALSE
]
stopifnot(
  nrow(tournament_ok) == 12L,
  all(tournament_ok$delta_top1 == 0),
  all(tournament_ok$delta_classification_oa == 0),
  mean(tournament_ok$peak_rss_mb_new) < mean(tournament_ok$peak_rss_mb_old)
)

toy_ok <- toy_paired[toy_paired$old_ok & toy_paired$new_ok, , drop = FALSE]
speedup <- vapply(split(toy_ok, toy_ok$method), function(x) {
  exp(mean(log(x$speedup_old_over_new)))
}, numeric(1))
stopifnot(
  speedup[["kema_rbf_full_label"]] > 50,
  speedup[["kema_rbf_reduced_label"]] > 50,
  speedup[["kema_linear_label"]] < 1
)

cat("KEMA benchmark evidence verified.\n")
