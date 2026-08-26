# Five-route five-task transfer evaluation v3 (scFoundation backbone, corrected endpoint)

Registered 2026-08-27 (amended before any v3 execution; the earlier v3 draft
of the same date was registered but never executed and is replaced by this
file). Supersedes the T2' route selection of REGISTRATION.md (see
ERRATUM.md). T1'/T3'/T4' task definitions are unchanged from
REGISTRATION.md; T2'' and G2 are corrected as below.

## Design decision (user directive 2026-08-27)

The framework is the original YeastBridge-Eval reused as-is with ONE
substitution: every scGPT backbone reference is replaced by scFoundation
(the model selected by the product line's retrieval evaluation). The
original routes' construction, corpora and training recipes are reused;
only the backbone changes. Pure-ESM2-space and old-scGPT-table arms from
the voided v2 execution are not part of the route set.

## Five routes

- A'' scF homolog-mapped: scFoundation MaeAutobin with yeast gene rows of
  pos_emb initialized from the mean human-ortholog scF pos_emb rows
  (OrthoDB + OMA + InParanoid union, verbatim rule of
  scripts/routeA/build_routeA_init.py; unmatched genes per-dim Gaussian,
  seed 42), then finetuned on the yeast corpus.
- B'' scF protein bridging: scFoundation backbone with the gene identity
  embedding replaced by a UCE-style injection layer proj(ESM2-650M vector)
  (verbatim ProteinEmbeddingInjector pattern of
  scripts/routeB/protein_inject.py, proj 1280->768 learnable, protein
  matrix frozen, missing genes share the bias row), then finetuned.
- C'' scF scratch control: identical construction and training with
  all-random pos_emb init (per-dim Gaussian, seed 42). Water line for the
  scF family.
- D' native scYeast: scYeast pos_embedding (unchanged original asset).
- E' SGA graph propagation: incumbent signed personalized PageRank
  alpha=0.30 with fusion gate (F1/F2 v2 representative; unchanged).

## Training protocol (registered, all three scF routes)

- Yeast vocab: gene_master systematic order (6,736 genes, identical row
  order to routeA_vocab.json) + 2 resolution positions + 1 pad row;
  scF config patched to seq_len 6,738. Pretrained weights loaded for all
  shared modules (token_emb, encoder, decoder); pos_emb replaced by the
  route init.
- Corpus: GSE125162 raw counts (data/single_cell/
  GSE125162_ALL-fastqTomat0-Counts.tsv, 38,225 cells), aligned to the
  6,736-gene vocab; input per cell = log1p(CPM-10k) of a uniform random
  sample of at most 1,200 expressed genes (the original
  finetune_scgpt.py train_args max_seq_len=1200 with the official
  trunc_by_sample behavior - amended 2026-08-27 after the first launch
  measured 0.56 s/step under the full-transcriptome convention, ~26 h per
  route; the 1,200-gene rule reproduces the original framework's compute
  envelope) + [4.0, log10(total count)] resolution tokens (official 't4'
  recipe). Position ids remain true gene column indices; subset tail
  padding uses the pad row.
- Masking: expressed genes p=0.30, zero-value genes p=0.03, pure mask to
  mask_token_id (the pretrain config's replace/random corruption is not
  reimplemented - registered deviation).
- Loss: MSE at masked positions between decoder reconstruction and true
  continuous values.
- Trainable: encoder layers 10-11 (last two) + token_emb + pos_emb (or
  injection layer) + decoder modules; encoder layers 0-9 frozen. This is
  the official finetune_model.py granularity (LinearProbingClassifier
  unfreezes exactly the last two encoder layers).
- Optimization: AdamW lr 1e-4, grad clip 1.0, batch 32 x accumulation 1
  (= original route A train_args batch_size 32; amended 2026-08-27
  together with the 1,200-gene input rule), epochs 6, seed 42, bf16
  autocast with full gradient checkpointing over encoder and decoder; GPU
  = whichever card has free headroom (user rule 2026-08-27), never
  crowding out the root services.
- Route output for evaluation: the trained gene identity table - pos_emb
  rows for A''/C'', proj(ESM2) rows for B'' (6,736 x 768 each).

## Corrected endpoint and gates

1. PRIMARY ENDPOINT = T2'' balanced-pool anchor recovery: each held-out
   anchor ranked within the fixed pool of 2,246 ortholog-mapped,
   non-anchor, on-graph genes; identical pool for every route; missing
   table rows score 0. Removes the mapped/unmapped stratum artifact
   (ERRATUM.md). v2 full-pool ranks co-reported.
2. Dense-route queries: signed weighted mean of the reduced anchors' own
   rows in the route's own table (uniform construction across A''/B''/
   C''/D'); E' unchanged.
3. G1: paired anchor-level bootstrap (10,000 resamples, seed 20260827) of
   the balanced-pool median rank vs incumbent E'; pass = CI95 upper < 0
   AND delta <= -0.02; if no route passes with G2, E' is retained and
   labeled "retained without positive evidence".
4. G2 water line = 20-draw permutation null (clean null: full weight
   assignment including signs permuted; seeds 20260827..20260846). Pass =
   route median rank beats the permuted median in >= 19/20 draws.
5. T5' anchor-efficiency recomputed on the balanced-pool metric.

## Claim boundary

Unchanged from REGISTRATION.md: transfer-mechanism selection only; no
therapeutic, efficacy, or response-regression claim; wet-lab remains the
final arbiter. The anchor population limitation (housekeeping collapse,
zero GPCR/ion-channel members) is unchanged and is disclosed in
ERRATUM.md; module-level anchors for the membrane-protein universe
require a separate registration.
