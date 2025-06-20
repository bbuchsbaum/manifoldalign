# Development Scripts

This directory contains development and testing scripts that are not part of the main package but are useful for development, testing, and benchmarking.

## Structure

### `test-scripts/`
Contains ad-hoc test scripts used during development:
- `test_*.R` - Various test scripts for different components
- These are separate from the formal testthat tests in `tests/testthat/`

### `benchmarks/`
Contains benchmark scripts and results:
- `benchmark_parrot.R` - Performance benchmarks for PARROT algorithm
- `eigen_*.png` - Benchmark visualization plots
- `eigen_benchmark_results.csv` - Benchmark result data

### Root level scripts
- `fix_parrot_preproc.R` - Script to fix preprocessing issues
- `run_parrot_full_validation.R` - Comprehensive validation suite

## Note

These files are excluded from the R package build via `.Rbuildignore` but are kept in the repository for development purposes. 