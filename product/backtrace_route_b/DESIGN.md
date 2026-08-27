# Product step 4: back-trace yeast results to human targets (route B, transposed) + convergence

Registered 2026-08-27 before execution. Config:
configs/product_backtrace.json. Symmetry principle: the back-map is the
transpose of the forward transfer - the SAME trained components (proj +
route-B gene table), zero new parameters.

## Mechanism

- Forward (step 2): score_yeast(g) = cos(proj(esm2(T)), table(g)).
- Backward (this step): for compound C with its condition-averaged
  sensitivity profile |z_C| over strains (target-neutral: no peeking at
  which dose won for which target), q_C = sum_g |z_C(g)| * table(g)
  (normalized); then score_human(h) = cos(proj(esm2(h)), q_C) over the
  1,177-universe protein embeddings. The compound's yeast action profile
  is weighted by its full |z| spectrum - cutoff-free, consistent with the
  step-3 endpoint.

## Part 1 - Precision calibration (runs first; OrthoDB as EXAM, not mechanism)

For every yeast gene with a unique human ortholog inside the universe
(OrthoDB 1:1 pairs, yeast_sgd <-> human_symbol), rank all universe
proteins by cos(proj(esm2(h)), table(g)) and record the true ortholog's
rank. Report top-1 / top-5 accuracy and MRR, plus the descriptive top-1
for yeast genes whose ortholog falls outside the universe. OrthoDB is
used ONLY as an evaluation gold standard here - it participates in no
product computation.

## Part 2 - Closure consistency + specificity + knowledge channel

For each significant (target T, compound C) pair from step 3 (q < 0.1):

- Closure: T's rank in C's back-mapped universe ranking.
- Specificity: per-compound null from n_perm_compound label permutations
  of |z_C| (same weight multiset, shuffled gene assignment); T's closure
  score is reported as a standardized excess (z_closure) over C's own
  null - broad compounds' global structure cancels against their own
  baseline.
- Knowledge channel: C's ChEMBL known targets (moa_labels.tsv, joined by
  inchikey) matched against the universe by normalized protein-name
  comparison (registered limitation: name-based matching, matched and
  unmatched counts reported; no action types in the label file, so
  direction alignment is deferred).

## Convergence rule (pre-declared)

Confidence tiers over the four evidence columns (execution q < 0.1 AND
closure z >= 2 AND known-target hit AND scF direction score present):
tier A = all four; tier B = any three; tier C = execution + closure only.
Within tiers, rank by closure excess, then by execution rho. Outputs:
results/backtrace_matrix.tsv, results/calibration.json,
results/convergence_top.tsv.

## Claim boundary

Back-mapped ranks are mechanism outputs of the trained protein-space
transfer; the calibration bounds their precision. Convergent pairs are
candidate hypotheses, not validated targets; wet-lab remains the arbiter.
