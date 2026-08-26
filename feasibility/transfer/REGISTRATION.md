# Transfer-method selection: F1 method comparison and F2 component ablation

Registered 2026-08-26 in yeastbridge_re, before any variant is executed.
This registration closes a recorded gap in the legacy project
(yeastbridge_vs): the human-intent -> yeast transfer step (OrthoDB anchors
-> SGA propagation) was design-specified and never competed against
alternatives; a propagation ablation was proposed on 2026-08-26 in the
legacy session but never executed.

All core code is copied from the legacy project (`second_round.py`,
`workspace.py`, verbatim); only the experiment driver wiring the registered
arms is new. All inputs are byte-copied frozen assets under
`feasibility/transfer/assets/` (orthology snapshot, rebuilt SGA edge list,
CRC state signature, scFoundation signature scores from product step 1).

## Shared evaluation endpoint (frozen)

Leave-one-out (LOO) anchor recovery on the SGA graph, ranked against
**non-anchor nodes only**:

- Anchor ORFs are built exactly as in `_feature_weights` (legacy
  `second_round.py`): each yeast ORF with >=1 human ortholog in the
  signature receives `anchor_weights[orf] = -mean(human weights)`.
- Eligible held-out anchors: ORFs with nonzero anchor weight and at least
  one SGA edge.
- For each eligible ORF o: the method under test runs with
  `anchor_weights` minus o; o's recovery score is
  `value[o] * sign(original weight of o)`.
- **Ranking pool (design note):** o is ranked only among nodes that are
  not anchors in the reduced run. Restart-type methods place the remaining
  true anchors at the top of a full-graph ranking by construction (the
  restart mass lands on them), so an all-node pool would measure
  self-competition of the anchor set rather than transfer fidelity. The
  question this endpoint answers is: given the remaining intent, does the
  method place o's ORF above ordinary non-anchor yeast genes?
- Ties use average rank (a zero-mass node tied with the bulk is not
  credited a top rank).
- Primary metric: normalized median rank of o within the non-anchor pool
  (rank / pool size; lower is better). Secondary metrics: hit rate at the
  top 1% and top 5% of the pool.
- Statistics: paired bootstrap over held-out ORFs (10,000 resamples, seed
  20260826) on the per-ORF normalized-rank difference (arm minus
  incumbent); percentile 95% CI. Decision rule per the project comparison
  contract: delta >= 0.02 AND CI95 lower bound > 0 means the incumbent
  ranks significantly better (lower).

## F1: method comparison arms

Same anchors, same graph, same endpoint.

1. `incumbent` — signed personalized PageRank, alpha = 0.30, per-iteration
   anchor restart, positive/negative channels propagated separately and
   recombined (the production method).
2. `no_propagation` — direct anchor weights only, no network.
3. `one_hop` — a single random-walk step from the anchors.
4. `full_diffusion` — neighbour-mean diffusion iterated to convergence
   without restart (the registered near-uniformity pathology case).
5. `rewired_graph_control` — incumbent propagation on a degree-preserving
   rewired SGA graph (edge-swap, seed 20260826). Tests whether SGA
   structure carries information beyond its degree sequence.
6. `permuted_anchor_control` — incumbent propagation with anchor weights
   randomly permuted across anchor ORFs (seed 20260826). Tests whether the
   anchor identities carry the information.

F1 gates: the incumbent must beat both controls (arms 5, 6) under the
decision rule to keep its "structure-informed transfer" wording; if
`no_propagation` or `one_hop` is not beaten, the network step is reported
as unjustified. All outcomes published.

## F2: component ablation arms

Component-wise removal from the incumbent:

1. `full` — identical to the F1 incumbent (shared run).
2. `no_fusion_gate` — raw statistical signature weights; the predeclared
   scFoundation dual-track modulation (agree x1.0 / disagree x0.5 /
   no-model x0.75) is not applied.
3. `unsigned` — anchor magnitudes propagated as a single unsigned channel
   (signs dropped).
4. `no_restart` — identical to F1 arm 4 (shared run).
5. `alpha_0.10`, `alpha_0.50`, `alpha_0.70` — incumbent with the
   propagation alpha moved off its design value 0.30.
6. `uniform_weights` — every anchor ORF weight set to sign(w) * 1 (intent
   magnitudes dropped).

F2 gates: an arm that is NOT beaten by `full` under the decision rule flags
its component as unjustified complexity; `full` losing to an arm flags the
incumbent configuration as improvable. All outcomes published.

## Claim boundary

This endpoint measures transfer fidelity (recovery of held-out anchor
identities through the graph), not compound-ranking quality; the
compound-level endpoint (same-target retrieval against ChEMBL labels)
remains queued behind the legacy moa_prep crawl and is unchanged by these
experiments. A method that wins here earns the wording "transfer method
selected by held-out anchor recovery"; no screening or product claim is
unlocked.
