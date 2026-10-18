# Alignment benchmark fixtures

These files are governed inputs to real-data alignment protocol 1.0.0. Exact
byte counts and SHA-256 digests are frozen in
`inst/benchmarks/real_alignment_artifacts.csv`.

## `aucs.mpx` and `monastery.mpx`

The network fixtures are unmodified data files extracted from the CRAN source
archive for `multinet` 4.3.4. They are redistributed under the Apache License
2.0. The versioned extraction command is:

```sh
Rscript inst/benchmarks/prepare_real_alignment_data.R \
  --output-dir=/path/to/cache \
  --extract-network-fixtures-to=inst/extdata/alignment-real
```

Source: <https://cran.r-project.org/package=multinet>

## `so2sat_lcz42_ci_features.csv`

This derived fixture contains summary attributes from actual Sentinel-1 and
Sentinel-2 patches in So2Sat-LCZ42 v4.2, frozen Hugging Face revision
`b5c817486899935e864832b93086cc87f3eee473`. It selects the first four samples
per class from both the validation-west and testing-east files, yielding 68
rows per split. Geographic metadata and labels come directly from the matching
official HDF5 artifacts.

Means, standard deviations, and extrema use every finite patch value. Five
quantiles use a deterministic evenly spaced 512-value grid, which keeps the
full 48,238-patch extraction linear-time and reproducible.

The fixture is redistributed under CC BY 4.0. It exists to test schemas,
identity, blocked split handling, and absence of target-test label leakage. Its
class-balanced accuracy is not representative and must not be used to rank
methods. Rebuild it with:

```sh
export MANIFOLDALIGN_ACCEPT_SO2SAT_LICENSE=yes
Rscript inst/benchmarks/prepare_real_alignment_data.R \
  --output-dir=/path/to/cache \
  --download-so2sat \
  --so2sat-artifacts=validation,validation_geo,testing,testing_geo
Rscript inst/benchmarks/extract_so2sat_ci_fixture.R \
  --data-dir=/path/to/cache \
  --output=inst/extdata/alignment-real/so2sat_lcz42_ci_features.csv \
  --per-class=4
```

Source: <https://huggingface.co/datasets/zhu-xlab/So2Sat-LCZ42>
