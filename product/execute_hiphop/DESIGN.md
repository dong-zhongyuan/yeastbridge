# Product step 3: execute the yeast tasks on the public HIP/HOP chemogenomic screen

Registered 2026-08-27.

## Semantics

Heterozygous-deletion hypersensitivity (HIP logic): a compound that
inhibits pathway P preferentially slows the deletion strains of P's genes.
The execution statistic asks: does compound C's response profile
preferentially involve the yeast genes of target T's transferred task?
If yes, C executes the transferred task for T.

## Endpoint (cutoff-free)

FULL-PROFILE Spearman correlation between each target's COMPLETE step-2
ranking (rank vector over all measured strains) and each non-vehicle
condition's |z| sensitivity profile. No gene-set constant of any kind.

- Null: strain-label permutation (1,000 draws, seed 42, target margins
  preserved), one-sided empirical p.
- Aggregation: best dose per InChIKey (dose and dose count reported).
- Multiple testing: Benjamini-Hochberg FDR over all (target, InChIKey)
  pairs (~3.8M).
- Effect size is Spearman rho itself.

## Inputs and outputs

Engine: the processed Lee-2014 HIP/HOP E-MTAB-2391 barcode response
matrix (5,668 strains x 3,850 arrays; vehicle-contrast z-scores; the
registered hiphop_response_record.json processing), read-only reuse from
yeastbridge_vs. Task lists: all 1,177 from product step 2. Outputs:
results/exec_matrix.tsv, results/per_target_top/<TARGET>.tsv,
results/execute_summary.json.

## Registered observation for step 4

Compound specificity is heterogeneous: the median compound hits 4
significant targets, but a tail hits 160-317. Broad-hit compounds carry
global response structure that correlates with many transfer rankings;
the step-4 convergence layer must apply a specificity treatment
(excess-over-breadth or per-compound shrinkage) before treating a pair
as target-specific.

## Claim boundary

Enrichment says the compound's yeast sensitivity profile overlaps the
transferred task genes. It is NOT a human-target claim (back-mapping is
step 4) and NOT compound-target identity (HIP sensitivity marks pathway
involvement, not binding). Doses differ per condition and are reported,
not harmonized.
