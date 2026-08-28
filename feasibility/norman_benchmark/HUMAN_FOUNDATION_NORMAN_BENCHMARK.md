# Norman multi-human-foundation benchmark

Status: completed on the A6000 server, 2026-08-17 (Asia/Shanghai).

## Decision first

This run establishes a real, reproducible three-backbone comparison and rejects
the idea that scGPT should be the default human model. It does **not** establish
a fusion advantage.

- **GO**: keep scGPT, Geneformer V1-10M and scFoundation as an explicit human
  model panel. All three checkpoints were freshly loaded and executed in their
  declared runtimes on the same task.
- **NO-GO**: do not claim that the three-model mean is better than the
  development-locked best foundation model. Locked-test delta was `0.00456`,
  paired 95% CI `[-0.00919, 0.02066]`.
- **NO-GO**: do not claim that a foundation model or their mean beats the
  additive linear baseline. Mean-versus-linear delta was `0.00079`, paired 95%
  CI `[-0.03083, 0.02974]`.
- **Operational rule**: model choice and fusion remain task-conditioned. The
  linear baseline is a mandatory comparator, not a cosmetic lower bound.

The full machine-readable result is
`reports/human_foundation_norman/full/summary.json` (identical to
`benchmark_full.json`), SHA-256
`721ae8daa966b62c4a18f42e2f9c5c7b6db8c3c5145a0bcfdf2dd7ea73d87050`.

## Question and task

The narrow first question was: *does a frozen pretrained gene-identity prior
transfer to a held-out genetic-perturbation response task, and does a simple
multi-backbone mean outperform the best single prior?*

The source is the existing Norman A549 Perturb-seq H5AD:

- 91,205 cells and 5,045 measured genes;
- input SHA-256
  `23ffb0fac6a847ff927cf7509d80d85052bfefbfb97610786a2dafaaefa0b6a0`;
- raw labels such as `CEBPE+ctrl` and `ctrl+CEBPE` are merged into one canonical,
  unordered perturbation-gene set;
- all cells in one canonical condition are aggregated to one mean expression
  profile;
- the target is condition mean minus the global control mean.

The response gene set is the Norman/Geneformer-Ensembl/scGPT-symbol/
scFoundation-symbol intersection: 3,066 genes. The shared perturbation input
set contains 100 of the 105 Norman perturbation genes. Conditions containing
`C19orf26`, `C3orf72`, `ELMSAN1`, `KIAA1804` or `RHOXF2BB` are excluded for all
models rather than imputed for selected providers. This leaves 223 eligible
non-control canonical conditions.

## Locked split and adapter

The grouping unit is the canonical unordered perturbation set. Fold assignment
uses only condition identity and perturbation cardinality, never response
values:

- 178 development conditions and 45 sealed final-test conditions;
- four development OOF folds;
- nested inner OOF selection of ridge alpha for every outer fold;
- zero canonical-condition overlap at every boundary;
- test labels are not used for model, alpha, feature, or ensemble selection.

Each checkpoint supplies its own legal identity representation:

| Model | Frozen feature | Width |
|---|---|---:|
| scGPT | `model.encoder.embedding(gene_token_id)` | 512 |
| Geneformer V1-10M | `BertForMaskedLM.get_input_embeddings()(Ensembl token)` | 256 |
| scFoundation | gene-branch `model.pos_emb(19,264-gene index)` | 768 |

Within every condition, each gene vector is L2-normalized, vectors are averaged,
and perturbation cardinality is appended. All three then use the same
multi-output ridge adapter. No post-perturbation expression is fed to a
backbone. The classical additive baseline uses the same 100-gene universe as a
multi-hot input; the mean baseline uses training responses only. The registered
fusion is the unweighted mean of the three foundation predictions.

This is intentionally a **frozen gene-prior transfer benchmark**, not a full
contextual cell-state simulation benchmark.

## Fresh checkpoint evidence

| Model | Weight SHA-256 | Parameters | Finite smoke | Time | Peak GPU |
|---|---|---:|---:|---:|---:|
| scGPT | `6cb5d451ab5c4b33eb673adbe4fddc61d2389df1b89b7651a9fe2e557572b922` | 51,856,898 | `(1, 8)` | 209.52 s | 340,684,288 B |
| Geneformer V1-10M | `a5e33a757431643b3697de7ef6127950cdc49e06e58d4266b3a3ab191b683f14` | 10,288,722 | `(1, 8, 256)` | 6.28 s | 51,739,648 B |
| scFoundation | `9f40bf324d3d0084c4b288d06f5af4fddd12206e2a3f022551d12e89e33a0ea9` | 119,252,285 | `(1, 8, 768)` | 13.15 s | 487,575,040 B |

Each provider report records `loaded_in_current_process: true`, interpreter and
package versions, asset hashes, feature-output hash, runtime, GPU and the exact
feature contract. Direct checkpoint/output hash verification is mandatory
before the scoring process accepts a provider.

## Results

Primary metric: mean across conditions of Pearson correlation across the 3,066
expression deltas.

| Method | Development OOF | Locked test | Test RMSE |
|---|---:|---:|---:|
| additive linear | **0.80183** | 0.76575 | **0.04482** |
| scFoundation | 0.79338 | 0.76198 | 0.05518 |
| three-foundation mean | 0.78947 | **0.76654** | 0.05406 |
| scGPT | 0.77640 | 0.76216 | 0.05516 |
| Geneformer V1-10M | 0.73451 | 0.74758 | 0.05986 |
| training mean | 0.52532 | 0.48922 | 0.08432 |

scFoundation was locked as the best foundation model from development OOF.
scGPT happened to be `0.00018` higher on the final test; switching to it after
seeing test values is prohibited. The three-foundation mean has the largest
test point estimate, but both its confidence interval versus locked
scFoundation and its interval versus the linear baseline include zero. The
observed `0.00456` fusion delta is also below the preregistered material delta
of `0.02`.

Additional paired condition-bootstrap estimates (10,000 resamples):

| Comparison | Test delta | 95% CI |
|---|---:|---:|
| three-foundation mean - locked scFoundation | 0.00456 | [-0.00919, 0.02066] |
| three-foundation mean - additive linear | 0.00079 | [-0.03083, 0.02974] |
| scGPT - additive linear | -0.00359 | [-0.03036, 0.02133] |
| scFoundation - additive linear | -0.00377 | [-0.03250, 0.02132] |
| Geneformer - additive linear | -0.01818 | [-0.06433, 0.02624] |

## GEARS disposition

The installed GEARS API supports `prepare_split(split="custom")`, so GEARS is
not conceptually incompatible with the condition split. It is excluded from
this release, fail-closed, because the only existing harness:

- runs its own legacy `simulation` split rather than this locked canonical
  split;
- internally chooses an epoch using validation and does not export condition-
  aligned raw OOF/test predictions;
- still contains a TODO for extracting its test metric;
- does not slice predictions to this immutable 3,066-gene output contract.

A strict GEARS comparison therefore needs five fresh trainings (four OOF plus
one final fit), a fixed epoch or nested validation policy, canonical source-label
mapping, and raw prediction export. The old GEARS number is not imported. Until
that runner exists, GEARS remains `not_run_fail_closed`, not a missing win or a
missing loss.

## Claim boundary and next experiment

This one A549 task cannot establish model superiority across diseases, cell
types or perturbation modalities. Pseudo-bulk aggregation removes cell-state
heterogeneity; five uncovered perturbation genes are excluded; the final test
has 45 conditions; and identity embeddings do not use baseline cell context.
Fine-tuning or contextual prompting may change rankings, but must receive a new
locked comparison rather than reinterpret this result.

The next human-panel benchmark should add baseline-cell context on a second
cell type/dataset, retain the same classical baseline, and preregister whether
fusion means simple averaging, OOF stacking or a task-conditioned router. No
fusion strategy is enabled globally from this result.

## Durable artifacts

- task manifest: `reports/human_foundation_norman/task/task_manifest.json`,
  SHA-256 `4cb2bd8ab2ce12f6bc56b2579d244662c5592c19db3a65b45347a7d086006908`;
- full raw OOF: `reports/human_foundation_norman/full/raw_oof_full.npz`,
  SHA-256 `e40919f37497658738c6a0919295cf0f907ad562c9ec2c8495c3b7b4adf7c8d9`;
- full raw test: `reports/human_foundation_norman/full/raw_test_full.npz`,
  SHA-256 `2341e4fedf9188139a9f4176841befd7e4885d6cfa851513caeee13c0e835146`;
- per-condition scores:
  `reports/human_foundation_norman/full/per_condition_scores_full.csv`,
  SHA-256 `19553ab9bdbe843924b32c7b0ce9097506194444284f986ab2643a4e5f6720f7`;
- provider reports: `reports/human_foundation_norman/provider_runs/`;
- execution-only pilot: `reports/human_foundation_norman/pilot/`.
