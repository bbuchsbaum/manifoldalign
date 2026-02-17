rd_files <- list.files("man", pattern = "[.]Rd$", full.names = TRUE)
missing_value <- character(0)
missing_examples <- character(0)

for (f in rd_files) {
  content <- readLines(f, warn = FALSE)
  txt <- paste(content, collapse = "\n")

  if (!grepl("\\\\value", txt)) {
    missing_value <- c(missing_value, basename(f))
  }

  if (grepl("\\\\keyword[{]internal[}]", txt)) next
  if (!grepl("\\\\examples", txt)) {
    missing_examples <- c(missing_examples, basename(f))
  }
}

cat("=== .Rd files missing value ===\n")
if (length(missing_value)) cat(missing_value, sep = "\n") else cat("NONE\n")
cat("\n=== Exported .Rd files missing examples ===\n")
if (length(missing_examples)) cat(missing_examples, sep = "\n") else cat("NONE\n")
cat("\nTotal .Rd files:", length(rd_files), "\n")
cat("Missing value:", length(missing_value), "\n")
cat("Missing examples (exported):", length(missing_examples), "\n")
