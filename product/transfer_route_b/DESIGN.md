# Product step 2: human functional target -> yeast-executable task (route B, verbatim)

Registered 2026-08-27. Mechanism: the five-route evaluation's preview (and
pending formal confirmation by the full registered run after route C''
finishes) selects route B'' (scF backbone + ESM2 injection layer) as the
transfer representation. This product step applies route B's mechanism
exactly as trained, with complete components:

- Query: ESM2-650M layer-33 mean-pooled embedding of the human target
  protein (47 registered candidates from product step 1), passed through
  the TRAINED injection projection (pos_emb.proj from
  feasibility/transfer_routes/scf_routes/B2/final_model.pt).
- Scoring: cosine against the trained route-B gene table (6,733 yeast
  genes x 768).
- Output: per-target full yeast gene ranking = the yeast-executable task
  (results/yeast_task_<TARGET>.tsv), top-K report and summary JSON.

No orthology table participates anywhere in this step; membrane targets
with no yeast ortholog are handled by construction. The intended_direction
(activate/inhibit) is carried through from product step 1 as task metadata;
mapping direction onto yeast perturbation mode (overexpress/delete) is
explicitly out of scope for v1 (requires the response layer; registered as
a follow-up).

Sanity check (registered): for GPCR-family targets, the median rank of the
yeast pheromone/MAPK pathway genes (KEGG sce04011, 114 genes) in the
transfer ranking is reported against the random median (~3,367). No gate -
descriptive calibration only.

Claim boundary: rankings are transfer-mechanism outputs, not perturbation
evidence; wet-lab remains the arbiter.

Selection caveat: the formal five-route selection is completed by the
registered run after C'' finishes; the preview (T2 0.8195 / T3 0.168,
only T2-qualified candidate) locks it barring an implausible C'' reversal.
If the formal run contradicts the preview, this product step's config is
re-pointed and re-run (single config change, no code).
