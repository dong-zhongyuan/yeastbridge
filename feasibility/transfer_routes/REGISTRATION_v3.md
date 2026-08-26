# Five-route five-task transfer evaluation v3 (final form: the original framework's own tasks)

Registered 2026-08-27 (amended before any v3 execution; replaces the two
earlier v3 drafts of the same date, which imported the anchor-recovery
endpoint from the yeastbridge_vs SGA-era chain. User directive 2026-08-27:
the evaluation is the ORIGINAL framework's five tasks on the new
scFoundation-backbone routes; the anchor-recovery design — including its
OrthoDB anchor construction — is an SGA-era residue and is dropped.)

## Routes

- A'' scF homolog-mapped: scFoundation MaeAutobin, yeast pos_emb rows
  initialized from human-ortholog scF embeddings (three-source orthology
  union, verbatim rule of scripts/routeA/build_routeA_init.py), finetuned
  on the yeast corpus (training protocol as registered below).
- B'' scF protein bridging: scFoundation backbone with the UCE-style ESM2
  injection layer (verbatim ProteinEmbeddingInjector pattern), finetuned
  identically.
- C'' scF scratch control: identical protocol, random init.
- D' native scYeast: original asset, evaluated as-is (feature
  routed_scyeast).
- E' SGA graph propagation: the incumbent yeastbridge_vs mechanism. It
  produces no gene-embedding table, so it enters the tasks by graph-native
  protocols where definable (below) and is reported N/A elsewhere.

Reference feature: esm2_mean (the original strongest T2 feature) is
co-reported for continuity with the original evaluation table.

## Training protocol (all three scF routes; unchanged from the previous
v3 drafts, already executing)

Yeast vocab = gene_master order (6,733 genes) + 2 resolution + 1 pad row;
input per cell = log1p(CPM-10k) of a uniform random sample of <=1,200
expressed genes (original train_args max_seq_len=1200, trunc_by_sample) +
t4 resolution tokens; masking p=0.30 expressed / 0.03 zero, pure mask;
MSE at masked positions; official finetune granularity (encoder layers
0-9 frozen, last 2 + token_emb + pos_emb/injector + decoder trainable);
AdamW lr 1e-4, clip 1.0, batch 32 x 1, epochs 6, seed 42, bf16 autocast
with full gradient checkpointing.

## Tasks — the original eval harness verbatim (yeastbridge/eval)

The three scF gene tables are registered at runtime into
eval.data._FEATURE_LOADERS (read-only reuse; the old project's files are
not modified) and evaluated by the original frozen protocols, same seeds
(42) and splits as the 2026-08-09 table:

- T2 essential genes: run_t2 logistic, AUROC/AUPRC (primary robustness
  task, as in the original route-B selection).
- T3 perturbation response: run_t3 ridge with val-alpha, Spearman +
  top-100 recall (the transfer-relevant yeast response task).
- T4 pathway engineering ranking: run_t4, hit@3/hit@5.
- T5 data efficiency: run_t5, T3 metric at 1%..100%.
- T1 cell state: DEFERRED (requires cell-level embedding extraction
  through the finetuned scF models; the original harness itself gates T1
  on the cell encoder; registered as a follow-up, not blocking).

Route E task-native protocols (registered):
- E-T2: leave-one-out personalized PageRank (the registered incumbent
  propagation, alpha=0.30) seeded from the other essential genes; the
  held-out essential gene is ranked within the non-essential on-graph
  pool; report median normalized rank and pooled AUROC.
- E-T4: for each engineering record, propagate from the pathway's other
  genes (+1 seeds) and rank the target within the pathway; hit@3/hit@5.
- E-T1/T3/T5: N/A (the mechanism produces no gene-embedding table; no
  ridge/KMeans input). Reported honestly, not imputed.

## Selection rule (pre-declared)

Primary: T3 Spearman. A route is selected if it holds the highest T3
Spearman among A''/B''/C''/D' AND its T2 AUROC is not below the esm2_mean
reference by more than 0.01. If no route satisfies this, no new
representation route is selected and the table is reported as-is (E
remains the incumbent graph mechanism by production continuity, labeled
"retained without new positive evidence"). Ties at T3 within 0.002 go to
the higher T2.

## Claim boundary

Unchanged: mechanism selection only; no therapeutic, efficacy, or
response-regression superiority claim beyond the tested tasks; wet-lab
remains the final arbiter. The membrane-target anchor-population
limitation documented in ERRATUM.md concerns intent-transfer endpoints
and does not apply to these yeast-benchmark tasks; the product branch for
membrane targets (module transfer / humanized yeast) remains separately
registered work.

## Dropped assets (retained on disk, not used by this registration)

signature_proteins.fasta / esm2_signature (full-intent protein set),
module_anchors.tsv (KEGG+GO module sets): inputs of the superseded
anchor-recovery design; kept for the future membrane-target registration.
