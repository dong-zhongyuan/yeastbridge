# Product step 2: human functional targets -> yeast-executable tasks (route B, verbatim)

Registered 2026-08-27.

## Scope

The FULL frozen target universe: 1,177 human membrane proteins (839 GPCR +
338 ion channels; UniProt reviewed 2026-08-13 snapshot, used as DATA).
Every universe member is transferred; target selection is deferred to the
convergence layer (step 4). No upstream pre-filter participates.

## Mechanism (route B, exactly as trained)

- Query: ESM2-650M layer-33 mean-pooled embedding of the target protein,
  passed through the TRAINED injection projection (pos_emb.proj from
  feasibility/transfer_routes/scf_routes/B2/final_model.pt - the selected
  five-route transfer representation).
- Scoring: cosine against the trained route-B gene table (6,733 yeast
  genes x 768).
- Output: per-target full yeast gene ranking = the yeast-executable task
  (results/yeast_task_<TARGET>.tsv, 1,177 files) + transfer_summary.json.

No orthology table participates anywhere in this step; membrane targets
with no yeast ortholog are handled by construction.

## Inputs

- Universe proteins: UniProt reviewed, fetched by gene_exact
  (1,177/1,177; inputs/universe_proteins.fasta).
- ESM2 embeddings: inputs/esm2_universe/ (0 NaN).
- inputs/universe_targets.tsv: family from the snapshot; scFoundation
  model direction score (from the retrieval reader) carried as a
  final-scoring column; 2 unscored genes left blank (= N/A).

## Sanity check (registered, descriptive only)

For GPCR-family targets, the median rank of the yeast pheromone/MAPK
pathway genes (KEGG sce04011, 114 genes) in the transfer ranking is
reported against the random median (~3,367). No gate.

## Claim boundary

Rankings are transfer-mechanism outputs, not perturbation evidence;
wet-lab remains the arbiter. Direction mapping onto yeast perturbation
mode (overexpress/delete) is out of scope here and belongs to the
response layer.
