# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

manifoldalign is an R package implementing Kernel Manifold Alignment (KEMA) and related algorithms for multi-domain data analysis. The package focuses on domain adaptation methods that project data from multiple domains into a shared latent space while preserving manifold structure and aligning same-class samples across domains.

## Development Commands

### Build and Check
```bash
# Build the package
R CMD build .

# Check the package
R CMD check manifoldalign_*.tar.gz

# Install the package locally
R CMD INSTALL .
```

### Testing
```bash
# Run all tests
Rscript -e "devtools::test()"

# Run specific test file
Rscript -e "testthat::test_file('tests/testthat/test-kema.R')"

# Run tests with more verbose output
Rscript -e "devtools::test(reporter = 'progress')"
```

### Documentation
```bash
# Generate documentation from roxygen2 comments
Rscript -e "devtools::document()"

# Build package website (if pkgdown is set up)
Rscript -e "pkgdown::build_site()"
```

### Development Workflow
```bash
# Load package for interactive development
Rscript -e "devtools::load_all()"

# Check code style (if lintr is available)
Rscript -e "lintr::lint_package()"
```

## Architecture and Key Components

### Core Algorithm Files
- **R/kema.R**: Main KEMA implementation with both exact and regression-based solvers. Contains extensive documentation about implementation choices and bug fixes.
- **R/genprocrustes.R**: Generalized Procrustes Analysis for aligning multiple datasets with partial observations.
- **R/pseudolabel.R**: Semi-supervised learning algorithms for label propagation.
- **R/gpca_align.R**: Generalized PCA alignment methods.

### Important Conventions
1. **Data Format**: All implementations follow "samples x features" convention (rows are samples, columns are features).
2. **Matrix Operations**: Heavy use of sparse matrices from the Matrix package for memory efficiency.
3. **Graph Construction**: Uses neighborweights package for building similarity graphs.
4. **Multi-block Data**: Uses multivarious package structures for handling multi-domain data.

### Current Development Status
- The package is under active development with ongoing bug fixes.
- Current focus is on resolving matrix dimension mismatches in graph construction (see KEMA_TEST_STATUS.md).
- Some tests are temporarily skipped due to these issues.

### Key Dependencies
Critical packages that must be available:
- Matrix (sparse matrix operations)
- PRIMME (eigenvalue computations)
- kernlab (kernel methods)
- glmnet (regression solvers)
- multivarious (multi-block data structures)
- neighborweights (graph construction)

#### External Package Information
For reference when working with this codebase:

**multivarious package**: Multi-block data analysis and dimensionality reduction
- Location: ~/code/multivarious/R or Context7 ID: /bbuchsbaum/multivarious
- Core concepts: projectors (dimensionality reduction), bi_projectors (two-way mappings), hyperdesign objects
- Key classes: projector, bi_projector, cross_projector, multiblock_projector
- Data format: "samples x features" convention (rows=samples, cols=features)
- Preprocessing: Sophisticated pipeline with forward/apply/reverse methods
- Caching: Expensive computations cached in projector .cache environment

**multidesign package**: Complex experimental design and multi-block data management  
- Location: ~/code/multidesign/R or Context7 (no exact match found)
- Purpose: Tools for managing complex experimental designs and multi-block data
- Core concepts: design matrix manipulation, cross-validation, data splitting
- Integrates with tidyverse ecosystem (tibble, dplyr)
- Supports lazy evaluation and memory-efficient processing

### Testing Strategy
- Comprehensive test suite in tests/testthat/
- Validation framework in R/kema-validation.R compares implementations against reference results
- Tests include both unit tests and integration tests with synthetic data

### Test Debugging Guidelines
When test failures occur, follow this systematic debugging approach:

1. **Read the exact error message first** - Don't jump to complex explanations
   - `is.matrix(x) is not TRUE` often means x is a Matrix (capital M) object, not matrix
   - Focus on what the error literally says before investigating deeper

2. **Check object types immediately** - Use simple diagnostics first:
   ```r
   class(result)           # What type of object is it?
   names(result)           # What components does it have?
   str(result, max.level=1) # Structure overview
   ```

3. **Matrix vs Matrix distinction** - This package uses sparse Matrix objects:
   - KEMA functions return sparse Matrix objects for efficiency
   - Tests should check `is.matrix(x) || methods::is(x, "Matrix")`
   - Both behave similarly but have different classes

4. **S3 dispatch verification** - For generic function issues:
   ```r
   methods('function_name')  # Check registered methods
   class(object)            # Confirm object class for dispatch
   ```

5. **Minimal reproduction** - Create the smallest possible test case:
   - Use the same data generation as working tests
   - Strip down to essential parameters only
   - Test one component at a time

6. **Common R package patterns**:
   - Functions may return S4 objects, lists, or Matrix objects
   - Check `$` vs `@` access patterns
   - Sparse matrices are the default for efficiency

### Known Issues and TODOs
1. Matrix dimension mismatches in multi-domain graph construction
2. Package DESCRIPTION needs completion (license, authors)
3. README.md is missing
4. Some files are empty placeholders (R/ssma.R, data-raw/kema_overview.md)