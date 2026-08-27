# Universe-wide execution results (steps 2+3 re-run at 1177-target scope)

Executed 2026-08-27 under the universe-wide amendments (commits 087224e,
bee7279; configs product_transfer_universe.json / product_execute_universe.json).

## Step 2 (transfer, route B verbatim)

1,177/1,177 universe proteins (UniProt reviewed, fetched by gene_exact;
ESM2-650M L33, 0 NaN) through the trained injection projection -> 1,177
yeast task rankings in product/transfer_route_b/results_universe/.

## Step 3 (HIP/HOP execution, v2 endpoint)

Full-profile Spearman, strain-label permutation null (1,000, seed 42),
best-dose per InChIKey, BH-FDR over 3,825,250 pairs.

- **7,100 pairs q<0.1 over 1,038 targets and 484 compounds**; rho up to 0.273.
- The 47-panel's 257 significant pairs: 105 survive universe-level FDR.
- **1,000 targets outside the old panel** carry significant executors - the
  removed pre-screen was hiding most of the executable landscape.
- Top emerging targets (max rho): STING1 0.273, GPR149 0.254, ADGRG4 0.248,
  GPR75 0.246, GPR160 0.244, ADGRV1 0.239, RYR3 0.230, RYR2 0.229, KCNH7
  0.213, PIEZO2 0.187, CATSPER3 0.186. Several carry high scF direction
  scores (STING1 0.41, RYR3 0.41, RYR2 0.40) - the convergence layer will
  use these.

## Known issue for step 4 (registered observation)

Compound specificity is heterogeneous: median compound hits 4 significant
targets, but a tail hits 160-317 (QZESEGHSLFKZIV 317, ZKDHPMYHSXIDJT 255,
...). Broad-hit compounds likely carry global response structure that
correlates with many transfer rankings; the step-4 convergence layer must
include a specificity treatment (excess-over-breadth or per-compound
empirical-Bayes shrinkage) before treating a pair as target-specific.

## Provenance note

No upstream statistical pre-filter is inherited anymore: the only inputs
from outside this project are DATA assets (frozen universe snapshot, HIP/HOP
matrix, KEGG/GO files). The 47-panel outputs are retained as reference sets.
