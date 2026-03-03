# Evaluate pseudolabeling performance against ground truth

Compares pseudolabeling results with known true cluster assignments.
Uses various clustering evaluation metrics for assessment.

## Usage

``` r
evaluate_pseudolabeling(predicted_labels, true_labels, verbose = TRUE)
```

## Arguments

- predicted_labels:

  Factor vector of predicted pseudolabels

- true_labels:

  Factor vector of true cluster assignments

- verbose:

  Whether to print detailed results (default: TRUE)

## Value

A list containing evaluation metrics:

- n_predicted_clusters:

  Number of predicted clusters found

- n_true_clusters:

  Number of true clusters in ground truth

- coverage:

  Proportion of samples assigned to clusters

- purity:

  Average purity of predicted clusters (proportion of dominant class)

- completeness:

  Average majority-class recall (proportion of true cluster members
  captured by their dominant predicted cluster)

- confusion_matrix:

  Confusion matrix between predicted and true labels

## Examples

``` r
# \donttest{
# Create synthetic data and apply pseudolabeling
synthetic <- create_synthetic_similarity_matrix(n_samples = 500)
result <- assign_pseudolabels(synthetic$sim_matrix)

# Evaluate performance
eval_result <- evaluate_pseudolabeling(result$labels, 
                                       synthetic$true_labels)
#> Pseudolabeling Evaluation Results:
#> ==================================
#> Predicted clusters: 33 
#> True clusters: 20 
#> Coverage: 91.4 %
#> Average purity: 100 %
#> Average completeness: 91.7 %
#> 
#> Confusion Matrix (top 10x10):
#>             true_clean
#> pred_clean   cluster_1 cluster_10 cluster_11 cluster_12 cluster_13 cluster_14
#>   anchor_001        21          0          0          0          0          0
#>   anchor_002         0          0          0          0          0          0
#>   anchor_003         0          0          0          0          0          0
#>   anchor_004         0          0          0          0          0          0
#>   anchor_005         0          0          0          0          0          0
#>   anchor_006         0          0          0          0          0          0
#>   anchor_007         0          0          0          0          0          0
#>   anchor_008         0          0          0          0          0          0
#>   anchor_009         0          0          0          0          0          0
#>   anchor_010         0          0          0          0          0          0
#>             true_clean
#> pred_clean   cluster_15 cluster_16 cluster_17 cluster_18
#>   anchor_001          0          0          0          0
#>   anchor_002          0          0          0          0
#>   anchor_003          0          0          0          0
#>   anchor_004          0          0          0          0
#>   anchor_005          0          0          0          0
#>   anchor_006          0          0          0          0
#>   anchor_007          0          0          0          0
#>   anchor_008          0          0          0          0
#>   anchor_009          0          0          0          0
#>   anchor_010          0          0          0          0
print(eval_result)
#> $n_predicted_clusters
#> [1] 33
#> 
#> $n_true_clusters
#> [1] 20
#> 
#> $coverage
#> [1] 0.914
#> 
#> $purity
#> [1] 1
#> 
#> $completeness
#> [1] 0.91653
#> 
#> $confusion_matrix
#>             true_clean
#> pred_clean   cluster_1 cluster_10 cluster_11 cluster_12 cluster_13 cluster_14
#>   anchor_001        21          0          0          0          0          0
#>   anchor_002         0          0          0          0          0          0
#>   anchor_003         0          0          0          0          0          0
#>   anchor_004         0          0          0          0          0          0
#>   anchor_005         0          0          0          0          0          0
#>   anchor_006         0          0          0          0          0          0
#>   anchor_007         0          0          0          0          0          0
#>   anchor_008         0          0          0          0          0          0
#>   anchor_009         0          0          0          0          0          0
#>   anchor_010         0          0          0          0          0          0
#>   anchor_011         0          0          0          0          0          0
#>   anchor_012         0          0          0          0          0          0
#>   anchor_013         0          0          0          0          0          0
#>   anchor_014         0          0          0          0          0          0
#>   anchor_015         0          0          0          0          0          0
#>   anchor_016         0         18          0          0          0          0
#>   anchor_017         0          2          0          0          0          0
#>   anchor_018         0          3          0          0          0          0
#>   anchor_019         0          2          0          0          0          0
#>   anchor_020         0          0         21          0          0          0
#>   anchor_021         0          0          2          0          0          0
#>   anchor_022         0          0          0         24          0          0
#>   anchor_023         0          0          0          0         12          0
#>   anchor_024         0          0          0          0          2          0
#>   anchor_025         0          0          0          0         10          0
#>   anchor_026         0          0          0          0          0         24
#>   anchor_027         0          0          0          0          0          0
#>   anchor_028         0          0          0          0          0          0
#>   anchor_029         0          0          0          0          0          0
#>   anchor_030         0          0          0          0          0          0
#>   anchor_031         0          0          0          0          0          0
#>   anchor_032         0          0          0          0          0          0
#>   anchor_033         0          0          0          0          0          0
#>             true_clean
#> pred_clean   cluster_15 cluster_16 cluster_17 cluster_18 cluster_19 cluster_2
#>   anchor_001          0          0          0          0          0         0
#>   anchor_002          0          0          0          0          0        23
#>   anchor_003          0          0          0          0          0         0
#>   anchor_004          0          0          0          0          0         0
#>   anchor_005          0          0          0          0          0         0
#>   anchor_006          0          0          0          0          0         0
#>   anchor_007          0          0          0          0          0         0
#>   anchor_008          0          0          0          0          0         0
#>   anchor_009          0          0          0          0          0         0
#>   anchor_010          0          0          0          0          0         0
#>   anchor_011          0          0          0          0          0         0
#>   anchor_012          0          0          0          0          0         0
#>   anchor_013          0          0          0          0          0         0
#>   anchor_014          0          0          0          0          0         0
#>   anchor_015          0          0          0          0          0         0
#>   anchor_016          0          0          0          0          0         0
#>   anchor_017          0          0          0          0          0         0
#>   anchor_018          0          0          0          0          0         0
#>   anchor_019          0          0          0          0          0         0
#>   anchor_020          0          0          0          0          0         0
#>   anchor_021          0          0          0          0          0         0
#>   anchor_022          0          0          0          0          0         0
#>   anchor_023          0          0          0          0          0         0
#>   anchor_024          0          0          0          0          0         0
#>   anchor_025          0          0          0          0          0         0
#>   anchor_026          0          0          0          0          0         0
#>   anchor_027         23          0          0          0          0         0
#>   anchor_028          0         23          0          0          0         0
#>   anchor_029          0          0         24          0          0         0
#>   anchor_030          0          0          0         23          0         0
#>   anchor_031          0          0          0          0         21         0
#>   anchor_032          0          0          0          0          0         0
#>   anchor_033          0          0          0          0          0         0
#>             true_clean
#> pred_clean   cluster_20 cluster_3 cluster_4 cluster_5 cluster_6 cluster_7
#>   anchor_001          0         0         0         0         0         0
#>   anchor_002          0         0         0         0         0         0
#>   anchor_003          0        23         0         0         0         0
#>   anchor_004          0         2         0         0         0         0
#>   anchor_005          0         0        21         0         0         0
#>   anchor_006          0         0         2         0         0         0
#>   anchor_007          0         0         0        25         0         0
#>   anchor_008          0         0         0         0        22         0
#>   anchor_009          0         0         0         0         2         0
#>   anchor_010          0         0         0         0         0        21
#>   anchor_011          0         0         0         0         0         0
#>   anchor_012          0         0         0         0         0         0
#>   anchor_013          0         0         0         0         0         0
#>   anchor_014          0         0         0         0         0         0
#>   anchor_015          0         0         0         0         0         0
#>   anchor_016          0         0         0         0         0         0
#>   anchor_017          0         0         0         0         0         0
#>   anchor_018          0         0         0         0         0         0
#>   anchor_019          0         0         0         0         0         0
#>   anchor_020          0         0         0         0         0         0
#>   anchor_021          0         0         0         0         0         0
#>   anchor_022          0         0         0         0         0         0
#>   anchor_023          0         0         0         0         0         0
#>   anchor_024          0         0         0         0         0         0
#>   anchor_025          0         0         0         0         0         0
#>   anchor_026          0         0         0         0         0         0
#>   anchor_027          0         0         0         0         0         0
#>   anchor_028          0         0         0         0         0         0
#>   anchor_029          0         0         0         0         0         0
#>   anchor_030          0         0         0         0         0         0
#>   anchor_031          0         0         0         0         0         0
#>   anchor_032         19         0         0         0         0         0
#>   anchor_033          3         0         0         0         0         0
#>             true_clean
#> pred_clean   cluster_8 cluster_9
#>   anchor_001         0         0
#>   anchor_002         0         0
#>   anchor_003         0         0
#>   anchor_004         0         0
#>   anchor_005         0         0
#>   anchor_006         0         0
#>   anchor_007         0         0
#>   anchor_008         0         0
#>   anchor_009         0         0
#>   anchor_010         0         0
#>   anchor_011        18         0
#>   anchor_012         2         0
#>   anchor_013         0         3
#>   anchor_014         0         3
#>   anchor_015         0        13
#>   anchor_016         0         0
#>   anchor_017         0         0
#>   anchor_018         0         0
#>   anchor_019         0         0
#>   anchor_020         0         0
#>   anchor_021         0         0
#>   anchor_022         0         0
#>   anchor_023         0         0
#>   anchor_024         0         0
#>   anchor_025         0         0
#>   anchor_026         0         0
#>   anchor_027         0         0
#>   anchor_028         0         0
#>   anchor_029         0         0
#>   anchor_030         0         0
#>   anchor_031         0         0
#>   anchor_032         0         0
#>   anchor_033         0         0
#> 
# }
```
