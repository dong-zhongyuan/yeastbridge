# Product step 3: execute yeast tasks via the public HIP/HOP chemogenomic screen

## UNIVERSE-WIDE AMENDMENT (2026-08-27, user directive; follows the step-2 amendment)

Step 3 re-executed at universe scope: all 1,177 target task lists from
product/transfer_route_b/results_universe against the same HIP/HOP matrix
and the same v2 endpoint (full-profile Spearman, strain-label permutation
n=1,000 seed 42, best dose per InChIKey, BH-FDR over all
target x InChIKey pairs - now ~3.8M pairs). Config:
configs/product_execute_universe.json; outputs: results_universe/.
The 47-panel execution (results_v2/) is retained as a reference set.
Multiple-testing burden rises accordingly; significance counts will
shrink - that is the honest cost of removing the pre-filter, and final
target selection is deferred to the step-4 convergence layer regardless.

Registered 2026-08-27 before execution. Input: the 47 per-target yeast task
rankings from product step 2 (product/transfer_route_b/results/
yeast_task_<TARGET>.tsv). Engine: the processed Lee-2014 HIP/HOP
E-MTAB-2391 barcode response matrix (5,668 strains x 3,850 arrays;
vehicle-contrast z-scores; the registered hiphop_response_record.json
processing), read-only reuse from yeastbridge_vs.

## AMENDMENT v2 (2026-08-27, superseding the statistic below; see SENSITIVITY.md)

The registered v1 statistic (top-67 task-gene enrichment) proved
K-fragile: individual (target, compound) pairs churn almost completely
across K in {34, 67, 134, 337} (Jaccard 0.01-0.10; cross-K core = 1 pair).
Endpoint v2, executed after the amendment: FULL-PROFILE Spearman
correlation between each target's COMPLETE step-2 ranking (rank vector
over all measured strains) and each non-vehicle condition's |z| profile -
cutoff-free, no arbitrary constant. Null = strain-label permutation
(1,000 draws, seed 42, margins preserved), empirical p -> best dose per
InChIKey -> BH-FDR over all (target, InChIKey) pairs. Effect size is rho
itself. v1 remains in the script for reproducibility only; its 268-pair
list is void as claims. Canonical v2 outputs: results_v2/.

## v1 protocol (superseded; kept for the record)

## Semantics

Heterozygous-deletion hypersensitivity (HIP logic): a compound that
inhibits pathway P preferentially slows the deletion strains of P's genes.
So the registered execution statistic is: does compound C's response
profile preferentially involve the yeast genes of target T's transferred
task? If yes, C "executes" the transferred task for T.

## Registered protocol

- Task set per target: top 67 genes of the step-2 ranking (top 1% of
  6,733; registered constant in configs/product_execute.json).
- Compound conditions: all non-vehicle arrays (InChIKey non-empty),
  grouped later by InChIKey (best dose kept, dose count reported).
- Statistic per (target, condition): mean |z| of the task genes'
  strains, standardized by an EMPIRICAL per-target null from 1,000 random
  67-gene sets drawn from the 5,668 measured strains (seed 42) - the
  reported z_exec is (observed mean - null mean) / null sd.
- Multiple testing: Benjamini-Hochberg FDR over all (target, InChIKey)
  best-dose pairs; q < 0.1 flagged as executors; per-target top-15
  reported regardless.
- Strains not present in the task ranking's gene universe are ignored;
  task genes without a measured strain are dropped with counts reported.

## Claim boundary

Enrichment says the compound's yeast sensitivity profile overlaps the
transferred task genes. It is NOT a human-target claim (back-mapping is
step 4) and NOT target identity for the compound (HIP sensitivity marks
pathway involvement, not binding). Doses differ per condition and are
reported, not harmonized.

## Outputs

results/exec_matrix.tsv (all target x InChIKey best-dose with z_exec, q),
results/per_target_top/<TARGET>.tsv (top-15 with dose, z, q, smiles),
results/execute_summary.json (counts, method, hashes).
