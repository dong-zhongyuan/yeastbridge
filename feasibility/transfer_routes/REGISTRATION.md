# Five-route five-task transfer evaluation (adaptation of YeastBridge-Eval)

Registered 2026-08-27 in yeastbridge_re, before any evaluation below is
executed. This experiment adapts the original feasibility framework
(`/public/home/mengxl/dzy/yeastbridge`, plan v1.0 2026-07-31, four routes x
five tasks) to the functional-target transfer step of the product line.

## Provenance and adaptation rationale

The original YeastBridge-Eval answered "which mechanism best carries human
foundation-model knowledge into yeast REPRESENTATIONS" (route B, protein
bridging, selected). It never evaluated the step the product line needs:
given a signed human intent, produce a ranked list of yeast genes to
perturb. The F1/F2 transfer-method selection (feasibility/transfer/,
v2 endpoint, executed 2026-08-26/27) compared propagation variants on the
SGA graph only: the incumbent signed personalized PageRank (alpha=0.30,
fusion gate on) holds the best median normalized rank (0.347) and every
component of it is ablation-justified, but the pre-registered gates versus
no-propagation and versus the degree-preserving rewired graph remain
unpassed (paired CIs contain zero). Method families other than graph
propagation were never entered. This experiment therefore promotes the
comparison to five ROUTES in fair competition, under five adapted TASKS.

## Design decision: scFoundation replaces scGPT at the signal level

In the original framework the human knowledge source was the scGPT-human
backbone being transferred at the WEIGHT level. Here the human knowledge
source is the scFoundation intent already produced by product step 1
(`product/target_scan/model_evidence_scfoundation_v1/signature_model_scores.tsv`
gating `state_signature.tsv`), transferred at the SIGNAL level. No
yeast-side scFoundation finetune is performed in this round; it is
registered as a possible future route, not one of the five.

## Five routes (transfer mechanisms; all consume identical anchor weights)

Anchor weights: scFoundation-fusion-gated signature mapped through OrthoDB
(verbatim `build_anchor_weights` from feasibility/transfer; 79 eligible
anchors on the 5857-node SGA graph).

- Route A' homolog-mapped representation: yeast gene embeddings of the
  original route A asset (scGPT finetuned with ortholog-init,
  `models/routeA/scgpt_yeast_ft.pt` encoder embeddings). Query = signed
  weighted mean of remaining anchors' embeddings; score = cosine.
- Route B' protein bridging: ESM2 650M layer-33 mean-pooled space (same
  extractor and weights as the yeast table). Query = signed weighted mean
  of ESM2 embeddings of the HUMAN ortholog proteins of remaining anchors
  (50 genes, 51 reviewed UniProt entries, downloaded 2026-08-27 into
  assets/). Score = cosine against the yeast ESM2 table. Cross-species by
  construction; motivated by the original no-ortholog-immunity ablation.
- Route C' from-scratch control: original route C asset (scGPT trained
  from scratch on the same yeast corpus). Same query construction in its
  own space. The space carries no human knowledge; water line for the
  embedding-family routes.
- Route D' native yeast representation: scYeast pos_embedding (5812x200,
  w_knowledge checkpoint). Same query construction in its own space.
- Route E' graph propagation: the incumbent signed personalized PageRank
  alpha=0.30 with fusion gate (v2-selected representative of the SGA
  family; all four components ablation-justified in F2).

Per-route permuted-anchor control (fixed seed 20260827 permutation of
anchor weights) is run for every route as its own water line.

## Five tasks (adapted from T1-T5)

- T1' state-direction geometry (from T1): per route, predicted
  post-perturbation state = sum over Kemmeren mutants of route_score(m) *
  KO_log2FC_profile(m); target = intent state signature mapped to yeast
  ORFs (orthology, anchor sign convention); metric = Spearman on shared
  measured genes, gene-level bootstrap CI. Full anchor set, no LOO.
- T2' anchor recovery (from T2; PRIMARY): leave-one-out over the 79
  anchors; score vector times sign(w_o); rank of the held-out anchor in
  the non-anchor pool. Primary statistic = median normalized rank (the
  corrected v2 endpoint). Co-reported: pooled AUROC = mean over anchors of
  1 - (rank - 1)/pool_size, and hit@1%/hit@5%.
- T3' perturbation matching (from T3): LOO anchors that themselves have a
  Kemmeren deletion profile; alignment of each mutant's KO profile to the
  route vector by cosine; metric = median normalized rank of the anchor's
  own mutant among all mutants with profiles, hit@5%. Descriptive if
  eligible N < 15 (registered downgrade, no gate).
- T4' engineering-record ranking (from T4): reader-generality probe,
  intent-free as in the original T4. For each literature engineering
  record (pathway P, target gene g, direction): query built from P's other
  genes (embedding-mean for A'-D'; uniform positive anchors propagated for
  E'); rank g within P's genes. Metric hit@3/hit@5 over records.
  Inherits the original T4 data-ceiling caveat (few records; essential
  genes structurally absent from the deletion-library-derived sources).
- T5' anchor-efficiency curve (from T5): anchor subsets at 25/50/75/100%
  x seeds {1,2,3}; per route the T2' median normalized rank on each
  subsample; reports the marginal ranking value of additional validated
  human anchors.

## Gates (pre-registered)

- G1 route selection: paired anchor-level bootstrap (10000 resamples, seed
  20260827) on T2' median normalized rank, best route vs incumbent E';
  material rule unchanged from the project contract: CI95 lower bound > 0
  AND delta >= 0.02. If no route beats E' by this rule, E' is retained
  (incumbent-retention rule, same as model selection).
- G2 water line: every route must beat its own permuted-anchor control
  with CI95 lower bound > 0 (no material-delta requirement); a route
  failing G2 is reported as not structure-informed regardless of G1.
- T1'/T3'/T4'/T5' are co-reported evidence layers, not claim-bearing
  gates.

## Claim boundary

This experiment selects the TRANSFER MECHANISM for ranking yeast genes
given a human intent, on LOO anchor recovery with co-reported layers. It
makes no therapeutic or efficacy claim, no cross-domain superiority
claim, and does not reopen any response-regression claim. Final
validation of any transferred target remains wet-lab (sentinel strains /
HIP-HOP), per the pipeline contract.

## Assets

- feasibility/transfer/assets/* (orthology, SGA corrected network,
  signature, scFoundation model scores) - reused verbatim
- yeastbridge/ models/routeA, models/routeC, scYeast checkpoint,
  data/processed/esm2_650m - read-only reuse
- assets/human_anchor_proteins.fasta + assets/esm2_human_anchors/
  (ESM2 650M layer 33, same extractor script, GPU 1) - new, hash to be
  recorded in results
- yeastbridge/data/ko_kemmeren/kemmeren_2014/kemmeren_t3_expr_log2fc.parquet
  + sample metadata - read-only reuse
- yeastbridge/data/pathway/engineering_records.tsv + pathway gene lists -
  read-only reuse

Driver: scripts/five_route_five_task.py (imports the registered
feasibility/transfer driver and the original eval loaders; only route
producers, task endpoints and statistics are new).
