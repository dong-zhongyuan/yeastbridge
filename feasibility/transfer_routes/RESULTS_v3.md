# Five-route five-task results (v3, original-framework tasks)

Executed 2026-08-27 under REGISTRATION_v3.md (original eval/tasks.py
verbatim, seed 42, runtime feature registration; route E by registered
graph-native protocols).

## Selection (registered rule)

**B2_scf_esm2inject selected** - highest T3 Spearman among T2-qualified
routes. Confirms the preview; product step 2 stands as built.

## Task table

| Feature | T2 AUROC | T3 Spearman | T4 hit@3/@5 | notes |
|---|---:|---:|---|---|
| B'' scF+ESM2 inject | 0.8195 | 0.1683 | 0.143/0.143 | selected |
| esm2_mean (reference) | 0.8237 | 0.1721 | 0.143/0.143 | pure-protein baseline |
| D' scYeast | 0.7377 | 0.1487 | 0.143/0.143 | |
| A'' scF+ortholog | 0.4805 | 0.0440 | 0.143/0.143 | random level |
| C'' scF scratch | 0.4687 | 0.0398 | 0.143/0.143 | random level |
| E' SGA (graph-native) | T2: median rank 0.085, AUROC 0.814 (LOO, 972 ess / 4545 pool) | N/A | hit@3 0.857, direction 0.143 | |

## Findings

1. B'' preserves essentially all of the pure-ESM2 signal through the scF
   backbone + finetune (0.819/0.168 vs 0.824/0.172) - the injection layer
   transmits protein-space structure near-losslessly.
2. BOTH plain scF finetunes (ortholog-init and scratch) collapse to random
   on T2/T3: the MaeAutobin value-reconstruction objective with trainable
   pos_emb does not organize a functionally-informative yeast gene table.
   The contrast against the scGPT-era routes (which reached 0.62-0.75 T2
   under full-model finetuning) suggests the scGPT masked-GENE objective
   organized identity embeddings while scF's value reconstruction does
   not. Route-level conclusion: on the scF backbone, protein injection is
   not merely robust - it is the only channel that yields a usable table.
3. E' graph protocols: strong essential-gene network coherence (median
   rank 0.085 among 4,545 non-essential) - the SGA graph remains a strong
   yeast-internal mechanism, complementary to B'' (its product use
   continues where a graph, not embeddings, is the natural substrate).
4. T4 direction consistency remains below the majority prior for every
   feature (unchanged from the original framework's own T4 conclusion).

## Consequences

- Product step 2 (transfer via B verbatim) stands; product step 3
  (HIP/HOP execution) results remain valid as executed.
- E' is retained for graph-native product uses (module/statistics on the
  SGA network), labeled as a separate mechanism, not the transfer winner.
