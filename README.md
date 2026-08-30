# YeastBridge-RE

**Cross-species functional transfer: executing human GPCR / ion-channel drug-target tasks in *Saccharomyces cerevisiae***

Entry for the Global Undergraduate Life Science Challenge (全球大学生生命科学挑战赛).

## Motivation

GPCRs and ion channels are the most important drug-target families, yet they have no orthologs in *Saccharomyces cerevisiae* — so the cheapest, highest-throughput eukaryotic model organism cannot be used directly to screen compounds against these targets. At the same time, predictions made by single-cell pretrained foundation models lack a low-cost, high-throughput experimental validation outlet.

This project connects the two problems. The core idea is to convert a functional objective from a higher organism into an executable task in yeast: the foundation model reads out the functional change required by the disease, and the yeast gene network executes that change and produces a screenable phenotype.

## Approach

1. **Disease-state definition.** Differential-expression signatures from colorectal cancer (CRC) single-cell transcriptomes define the disease state and the desired state, from which candidate GPCR / ion-channel targets and their desired regulation directions are identified.
2. **Foundation-model selection.** Candidate foundation models are competitively evaluated on Norman perturbation data, and the best performer is selected.
3. **Transfer-strategy selection.** Five candidate transfer routes are systematically compared. The final route injects human-protein ESM2 embeddings into the yeast gene table, achieving human-to-yeast cross-species functional transfer.
4. **Screenable phenotype.** After wet-lab validation, the transfer-task rankings are correlated genome-wide against HIP/HOP chemogenomic screening profiles, so that "executing the task of a human target" becomes a phenotype measurable at the strain level. Layered supporting evidence includes pharmacological annotation, direction matching, conformational selection, molecular docking, and pharmacodynamics.

## Key results

Four direction-consistent target–compound pairs were obtained:

- **ADRA2C**
- **ADRA2B**
- **KCNK2 / TREK-1**
- **OPRM1**

For the first time, human membrane-protein targets that lack yeast orthologs can be screened for compounds in a yeast system, providing an engineering implementation path for the combined strategy of "human foundation model + low-cost model-organism validation".

## Repository layout

```
configs/       configuration files
scripts/       pipeline entry points (CRC scan, transfer routes, HIP/HOP execution, docking, MD)
src/           yeastbridge_re python package
data/          datasets (large raw data not tracked; see .gitignore)
feasibility/   feasibility and evidence-chain experiments (incl. MODEL_EVIDENCE_CHAIN.md)
product/       delivery artifacts (target screen, ChEMBL branch, HIP/HOP execution)
result/        final results
```

## Installation

```bash
pip install -r requirements.txt
```

External tools used by the docking / MD branches (not pip-installable): AutoDock Vina, OpenBabel, GROMACS.

## Reproduction

The main pipeline entry points live in `scripts/` (e.g. `crc_scan/`, `five_route_five_task.py`, `product_transfer_route_b.py`, `product_execute_hiphop.py`). The full evidence chain is documented in `feasibility/MODEL_EVIDENCE_CHAIN.md`.
