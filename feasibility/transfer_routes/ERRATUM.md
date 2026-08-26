# ERRATUM: five-route T2' selection of route A' is void (stratum-separation artifact)

Registered 2026-08-27 after read-only verification. No project file was
modified to produce the evidence below; verification scripts live at
dzy/verify_g2.py + dzy/verify_g2.log and the balanced-pool/init-only
checks were run as one-shot read-only computations.

## What is voided

The v2-endpoint T2' selection "A_homolog_repr (passed G1 vs incumbent and
G2 water line)" (results/five_route_results.json, first execution under
REGISTRATION.md). A' median normalized rank 0.012, pooled AUROC 0.837,
hit@1% 0.481.

## Verified mechanism (three checks)

1. INIT-ONLY control: recomputing the T2' endpoint using only
   routeA_init_embeddings.npy (the transplanted human scGPT geometry
   before any yeast finetune) gives median rank 0.535 (random). The
   anchors do NOT cluster in the raw transplanted geometry, so the win is
   not a pass-through of human-side coherence and the label is not
   directly derivable from the feature construction.
2. ORTHOLOG-BALANCED POOL control: ranking each held-out anchor only
   against the 2,246 ortholog-mapped non-anchor genes collapses A' from
   median 0.012 to 0.359, statistically indistinguishable from the
   from-scratch control C' (0.361). The v2 advantage was separation
   between the mapped and unmapped strata of the finetuned space: 61% of
   the full pool is unmapped, and beating that stratum is booked by the
   endpoint as top-1% recovery. Permuting anchor weights does not move
   the strata, which is why the G2 water line could not detect this.
3. Co-reported layers contradicted a transfer claim for A' at execution
   time: T1' Spearman 0.097 CI [-0.128, 0.309] (contains zero, equals
   C'); T3' median rank 0.482 (random level; N=14 descriptive); T4'
   hit@3 = 0.

Consequence: A's T5' efficiency curve from the same execution is void for
the same reason. Routes B' (ESM2 space, orthology-independent), C'
(scratch), D' (scYeast native) and E' (graph propagation) are not
affected by this artifact (C' moves 0.369 -> 0.361 under the balanced
pool).

## G2 finding (registered for v3)

The F1 permuted-anchor arm scored the held-out anchor with its TRUE sign
while permuting the propagated weights (hybrid null). The five-route G2
permutes the full assignment including the sign (clean null). Under the
clean null, single permutation draws at N=79 are unstable: seed 20260826
gives delta -0.048 CI [-0.119, +0.015]; seed 20260827 gives -0.125 CI
[-0.151, +0.044]; neither is significant. The main-arm reproduction is
bit-identical (0.34721), so this is a null-definition and draw-variance
issue, not a statistic implementation error. v3 mandates a 20-draw
permutation null.

## Disposition

The v2 selection is void. Per the registered incumbent-retention rule the
incumbent E_sga_propagation is retained, labeled "retained without
positive evidence" (its own water line is unproven under the clean null;
F1's gates vs no-propagation and vs the rewired graph also failed). The
corrected rerun is registered in REGISTRATION_v3.md and supersedes this
execution for route selection.
