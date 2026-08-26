# Five-route five-task transfer evaluation v3 (corrected endpoint)

Registered 2026-08-27 before any v3 execution. Supersedes the T2' route
selection of REGISTRATION.md (see ERRATUM.md). T1'/T3'/T4' definitions
are unchanged from REGISTRATION.md and are re-executed unchanged.

## Corrections over v2

1. PRIMARY ENDPOINT = T2'' balanced-pool anchor recovery: each held-out
   anchor is ranked against a fixed global pool of the 2,246
   ortholog-mapped, non-anchor, on-graph genes. Every route scores the
   same pool; a pool gene absent from a route's embedding table scores 0.
   This removes the mapped/unmapped stratum separation that produced the
   voided A' result. The v2 full-pool ranks are co-reported for
   continuity.
2. SIX ROUTES: the five of REGISTRATION.md plus B''_original_routeb, the
   representation of the originally selected feasibility winner (route B:
   ESM2 protein-injection layer expanded into the scGPT vocabulary and
   finetuned on the yeast corpus; models/routeB/
   routeb_gene_embeddings.npy, loaded via the original eval.data loader).
   This enters the competition because the original YeastBridge-Eval
   selected route B and that asset had not yet been tested on the
   transfer endpoint. Queries for all dense routes remain signed
   weighted means of the reduced anchors' own embeddings; B' keeps its
   human-ortholog ESM2 query.
3. G2 water line = 20-draw permutation null (clean null: the full
   weight assignment including signs is permuted; seeds 20260827..20260846).
   Pass rule: the route's balanced-pool median normalized rank is better
   (lower) than the permuted median in >= 19/20 draws (one-sided, p<=0.05).
4. G1 unchanged in form: paired anchor-level bootstrap (10,000 resamples,
   seed 20260827) of the balanced-pool median rank, best route vs
   incumbent E'; material rule CI95 upper < 0 AND delta <= -0.02; if no
   route passes with G2, E' is retained and labeled as such.

## Interpretation notes (frozen)

- The balanced pool answers "can the route rank the held-out anchor
  against genes that are equally ortholog-mapped" - i.e. transfer signal
  within the mapped stratum only. It does not measure reach into the
  unmapped 61% of the graph; that reach is exactly what the v2 artifact
  faked and what no current endpoint measures legitimately.
- A' and B'' both carry ortholog-dependent init (A' via the ortholog-init
  table, B'' via the ESM2 injection layer), so their balanced-pool
  results still contain any mapped-stratum-internal advantage those
  inits may confer. The balanced pool removes the stratum artifact, not
  the shared provenance; claims for these routes are bounded accordingly.
- T5' anchor-efficiency curves are recomputed under the balanced-pool
  metric; the v2 curves (full-pool) are void for A' per ERRATUM.md and
  superseded for all routes by this run.

## Claim boundary

Unchanged from REGISTRATION.md: transfer-mechanism selection only; no
therapeutic or efficacy claim; wet-lab remains the final arbiter.
