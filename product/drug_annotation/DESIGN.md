# Product step 4 (compound-centric): drug result annotation + physical validation

Registered 2026-08-27.

## Architecture

The product output is the significant compounds from step 3 (execution on
HIP/HOP) together with the targets whose transferred tasks they execute
(forward identity: every task originates from one universe target). The
yeast transfer exists precisely because membrane targets have no homologs;
the compounds are the result, and target identity comes from the
compounds' pharmacology.

- Target annotation: for each significant compound, its known human
  targets are fetched from ChEMBL (inchikey -> activities -> target gene
  symbols; results/chembl_targets.tsv, resumable fetch with backoff).
- Convergence (validation, not gate): a compound's ChEMBL targets
  intersected with the universe and with the targets whose tasks it
  executes. Three cases, all reported:
  1. ChEMBL target IS the executed-task target - convergent pair
     (strongest).
  2. Compound has ChEMBL targets but none match the executed task -
     divergence reported honestly.
  3. Compound has no ChEMBL targets (screening molecules) - the
     executed-task target stands as a novel hypothesis, carried to
     docking/MD.
- Promiscuity: compounds with many ChEMBL targets are annotated with
  their target count; intersection with executed tasks remains the
  focusing rule.

## Physical validation (next step)

Docking + MD on the selected (compound, target) pairs:
- Structure inventory per target: deposited PDB structures vs homology
  models required (LPAR1/KCNQ2 precedents; receptor preparation under
  the certified CONF-01 workflow for deposited structures).
- Docking: Vina workflow (container docking env, read-only use).
- MD: the multi-replica release-ensemble protocol for the top pairs; GPU
  queue subject to the free-card rule.

## Claim boundary

Compounds are candidate effectors of the transferred tasks; convergent
ChEMBL annotation strengthens target identity; docking/MD assess physical
plausibility only. Wet-lab remains the arbiter.
