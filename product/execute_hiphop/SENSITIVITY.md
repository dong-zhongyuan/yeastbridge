# Step-3 K-sensitivity check (diagnostic, post-registered-run)

Question: the registered task-set cutoff K=67 (top 1% of 6,733 - a
convention, frozen before execution, not data-derived). Do the significant
(target, compound) pairs survive other K choices?

## Sweep (same protocol, K varied via CLI override; q<0.1)

| K | pairs | targets | overlap with K=67 | Jaccard vs K=67 |
|---|---|---|---|---|
| 34 (0.5%) | 435 | 34 | 62 (23%) | 0.10 |
| **67 (registered)** | **268** | **36** | - | - |
| 134 (2%) | 218 | 38 | 28 (10%) | 0.06 |
| 337 (5%) | 366 | 44 | 5 (2%) | 0.01 |

Robust core (q<0.1 at ALL K): 1 pair. K=67 ∩ K=134: 28 pairs.

## Finding

The aggregate phenomenon is present at every K (hundreds of pairs, 34-44
targets) but the INDIVIDUAL pair list is K-fragile: which compounds top a
target's list depends strongly on the arbitrary cutoff. The per-pair
claims from the registered run must not be treated as individually robust.

## Consequence (registered follow-up)

Step-3 endpoint v2: replace the hard top-K cutoff with a cutoff-free
RANK-WEIGHTED statistic (mean |z| over the full task ranking, weights
decaying with rank, e.g. exp(-rank/tau)), same empirical null and FDR.
Removes the arbitrary constant entirely. Step 4 (back-trace) consumes the
v2 output, not the K=67 pair list. Until v2 runs, no individual pair from
this execution is carried forward as a claim.
