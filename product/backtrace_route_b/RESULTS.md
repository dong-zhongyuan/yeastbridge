# Step-4 results (first execution, 2026-08-27)

## Calibration (OrthoDB as exam)

Gold set = 8 yeast-human 1:1 pairs whose human ortholog lies inside the
1,177 membrane universe (structural intersection: yeast has almost no
GPCR/channel orthologs). top-1 1/8, top-5 1/8, MRR 0.133, median rank
138.5 vs random 589 (4.3x enrichment). MECHANISM SANITY ONLY - n=8
supports no precision claim; recorded as such.

## Closure (transposed back-map)

- 97 of 7,100 significant pairs at z_closure >= 2 (per-compound permuted-
  weight null); 86 pairs with back-mapped rank in the top 10 of 1,177.
- Top tierB rows: KCNQ5 closure rank 33 z=3.4; TRPC5 rank 18 z=3.3;
  ADGRF4 rank 17 z=3.1; compound SPABMCLDHGNLFQ back-maps coherently
  into the channel family (KCNQ5/SCN8A/SCN1A) - family-level consistency
  as expected for protein-space transfer.

## Knowledge channel - STRUCTURAL GAP (registered finding)

Only 22 of 484 significant compounds appear in moa_labels.tsv (that file
covers the LINCS paired-compound selection, a different population), and
its known targets are mostly outside the membrane universe -> zero
matches. Fix registered: query ChEMBL directly for the 484 significant
compounds (inchikey -> activities -> targets), intersect with the
universe, then recompute tiers. Until then tier A is undefined; tier B
(execution + closure + scF, 97 pairs) is the current top layer.

## Outputs

results/backtrace_matrix.tsv (7,100 rows), calibration.json,
convergence_top.tsv (top 100).
