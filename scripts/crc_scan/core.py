"""Reproducible CRC state-to-target scan with a real CRISPRi baseline.

The module deliberately separates three things that are often conflated:

* observational state evidence from GSE201348/Becker;
* measured HCT116 CRISPRi responses from X-Atlas/Orion; and
* proposed yeast/mammalian validation routes.

No virtual edge is accepted as a label, and no target receives a named-target bonus.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import platform
import sys
import tempfile
from collections import defaultdict
from collections.abc import Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

UTC_COMPAT = timezone.utc  # noqa: UP017 - the data-provider environment is Python 3.10


class ScanError(ValueError):
    """Raised when an input violates the frozen scan contract."""


CALIBRATION_TARGETS = {"LPAR1", "KCNQ2"}
SIGNED_DIRECTIONS = {"inhibit", "activate"}


def build_target_universe(
    gpcr_path: str | Path,
    ion_channel_path: str | Path,
) -> list[dict[str, Any]]:
    """Build a deterministic universe from complete snapshot rows, not seed targets."""

    combined: dict[str, dict[str, Any]] = {}
    for family, path in (("gpcr", gpcr_path), ("ion_channel", ion_channel_path)):
        for record in _read_uniprot_snapshot(Path(path), family):
            target = record["target_id"]
            if target in combined:
                previous = combined[target]
                families = sorted(set(previous["target_family"].split(";")) | {family})
                previous["target_family"] = ";".join(families)
                previous["membership_sources"] += f";{record['membership_sources']}"
                previous["uniprot_entries"] += f";{record['uniprot_entries']}"
            else:
                combined[target] = record
    universe = [combined[target] for target in sorted(combined)]
    if not universe:
        raise ScanError("target universe is empty")
    for record in universe:
        record["is_calibration_control"] = record["target_id"] in CALIBRATION_TARGETS
        record["ranking_bonus"] = 0.0
    return universe


def pareto_fronts(
    rows: Sequence[Mapping[str, Any]],
    axes: Sequence[str],
    *,
    eligible_key: str = "baseline_eligible",
) -> list[int | None]:
    """Assign non-dominated fronts without collapsing evidence into a hidden score."""

    if not axes:
        raise ScanError("at least one Pareto axis is required")
    remaining = {
        index
        for index, row in enumerate(rows)
        if bool(row.get(eligible_key)) and all(_finite_number(row.get(axis)) for axis in axes)
    }
    fronts: list[int | None] = [None] * len(rows)
    front_number = 1
    while remaining:
        current: list[int] = []
        for candidate in sorted(remaining):
            if not any(
                _dominates(rows[other], rows[candidate], axes)
                for other in remaining
                if other != candidate
            ):
                current.append(candidate)
        if not current:
            raise ScanError("Pareto assignment failed to make progress")
        for index in current:
            fronts[index] = front_number
        remaining.difference_update(current)
        front_number += 1
    return fronts


def run_scan(
    *,
    project_root: str | Path,
    config_path: str | Path,
    becker_h5ad: str | Path,
    xatlas_dir: str | Path,
    xatlas_gene_metadata: str | Path,
    norman_h5ad: str | Path,
    ortholog_map: str | Path,
    output_dir: str | Path,
    hash_workers: int = 4,
    hash_inputs: bool = True,
) -> dict[str, Any]:
    """Run the frozen state, coverage, perturbation and route audit."""

    started = datetime.now(UTC_COMPAT)
    root = Path(project_root).resolve()
    config_source = Path(config_path).resolve()
    config = _load_config(config_source)
    gpcr_path = _resolve_project_input(root, config["universe"]["gpcr_source"])
    ion_path = _resolve_project_input(root, config["universe"]["ion_channel_source"])
    metadata_path = _resolve_project_input(root, config["state_baseline"]["metadata_source"])
    becker_path = _required_file(becker_h5ad, "Becker h5ad")
    xatlas_root = Path(xatlas_dir).resolve()
    gene_metadata_path = _required_file(xatlas_gene_metadata, "X-Atlas gene metadata")
    norman_path = _required_file(norman_h5ad, "Norman h5ad")
    ortholog_path = _required_file(ortholog_map, "yeast ortholog map")
    destination = Path(output_dir).resolve()
    destination.mkdir(parents=True, exist_ok=True)
    xatlas_files = sorted(xatlas_root.glob("HCT116_Batch*.parquet"), key=_batch_sort_key)
    if not xatlas_files:
        raise ScanError(f"no HCT116_Batch*.parquet files found in {xatlas_root}")

    universe = build_target_universe(gpcr_path, ion_path)
    target_ids = [str(row["target_id"]) for row in universe]
    state, signature, state_qc = _build_state_baseline(
        becker_path,
        metadata_path,
        gene_metadata_path,
        target_ids,
        config["state_baseline"],
    )
    perturbation, xatlas_qc = _scan_xatlas(
        xatlas_files,
        gene_metadata_path,
        target_ids,
        signature,
        str(config["perturbation_baseline"]["control_label"]),
    )
    norman_coverage, norman_qc = _audit_norman(norman_path, set(target_ids))
    orthologs = _load_orthologs(ortholog_path)
    rows = _assemble_candidates(
        universe,
        state,
        perturbation,
        norman_coverage,
        orthologs,
        config,
    )
    minimum_xatlas_cells = int(config["perturbation_baseline"]["minimum_good_perturbed_cells"])
    xatlas_qc["targets_n_ge_minimum_cells"] = sum(
        int(row["xatlas_good_cells"]) >= minimum_xatlas_cells for row in rows
    )
    axes = [str(axis) for axis in config["ranking"]["axes"]]
    fronts = pareto_fronts(rows, axes)
    for row, front in zip(rows, fronts, strict=True):
        row["pareto_front"] = front if front is not None else ""
    rows.sort(
        key=lambda row: (
            0 if row["baseline_eligible"] else 1,
            int(row["pareto_front"]) if row["pareto_front"] != "" else 10**9,
            str(row["target_id"]),
        )
    )

    universe_counts = _universe_counts(universe)
    eligible_count = sum(bool(row["baseline_eligible"]) for row in rows)
    eligible_by_direction = {
        direction: sum(
            bool(row["baseline_eligible"]) and row["intended_direction"] == direction
            for row in rows
        )
        for direction in sorted(SIGNED_DIRECTIONS)
    }
    eligible_direct_count = sum(
        bool(row["baseline_eligible"]) and bool(row["signed_direction_directly_perturbed"])
        for row in rows
    )
    if eligible_direct_count:
        status = "measured_crispri_retrospective_prioritization_not_target_discovery"
    elif eligible_count:
        status = "counterfactual_prioritization_only_no_direct_signed_ranking"
    else:
        status = "coverage_only_no_credible_ranking"
    calibration_audit = _calibration_audit(rows)
    blockers = _blockers(config, state_qc, xatlas_qc, norman_qc, rows, hash_inputs)
    qc = {
        "schema_version": "1.0",
        "status": status,
        "query": config["query"],
        "universe": universe_counts,
        "state": state_qc,
        "xatlas": xatlas_qc,
        "norman": norman_qc,
        "ranking": {
            "eligible_candidates": eligible_count,
            "eligible_by_signed_direction": eligible_by_direction,
            "eligible_direct_same_direction_perturbations": eligible_direct_count,
            "eligible_opposite_direction_counterfactuals": eligible_count - eligible_direct_count,
            "credible_signed_ranking": bool(eligible_direct_count),
            "pareto_front_1_candidates": sum(row["pareto_front"] == 1 for row in rows),
            "axes": axes,
            "hidden_total_score_present": any("total_score" in row for row in rows),
            "threshold_source": "predeclared_config_before_scan",
            "predictive_test_set_used": False,
        },
        "calibration_control_audit": calibration_audit,
        "blockers": blockers,
        "claim_boundary": (
            "This is measured-CRISPRi coverage plus retrospective/counterfactual prioritization, "
            "not prospective target discovery, drug efficacy, or proof that the human-plus-yeast "
            "system outperforms a human-only system. A credible signed ranking requires direct "
            "same-direction perturbation evidence."
        ),
    }

    universe_fields = _ordered_union_fields(universe)
    candidate_fields = _ordered_union_fields(rows)
    signature_fields = _ordered_union_fields(signature)
    _write_tsv(destination / "target_universe.tsv", universe, universe_fields)
    _write_tsv(destination / "candidate_baseline.tsv", rows, candidate_fields)
    _write_tsv(destination / "state_signature.tsv", signature, signature_fields)
    _write_json(destination / "qc.json", qc)

    input_paths = [
        config_source,
        gpcr_path,
        ion_path,
        metadata_path,
        becker_path,
        gene_metadata_path,
        norman_path,
        ortholog_path,
        *xatlas_files,
    ]
    input_manifest = _input_manifest(
        input_paths,
        xatlas_root=xatlas_root,
        workers=hash_workers,
        hash_inputs=hash_inputs,
    )
    _write_json(destination / "input_manifest.json", input_manifest)
    summary = {
        "schema_version": "1.0",
        "status": status,
        "query": config["query"],
        "universe": universe_counts,
        "eligible_candidates": eligible_count,
        "eligible_by_signed_direction": eligible_by_direction,
        "eligible_direct_same_direction_perturbations": eligible_direct_count,
        "eligible_opposite_direction_counterfactuals": eligible_count - eligible_direct_count,
        "credible_signed_ranking": bool(eligible_direct_count),
        "xatlas_targets_n_ge_20": xatlas_qc["targets_n_ge_minimum_cells"],
        "norman_clean_single_target_coverage": norman_qc["clean_single_target_count"],
        "calibration_control_audit": calibration_audit,
        "input_content_addressed": bool(input_manifest["content_addressed"]),
        "blocker_count": len(blockers),
    }
    _write_json(destination / "summary.json", summary)
    implementation_base = Path(__file__).resolve().parent
    implementation_files = sorted(implementation_base.glob("*.py"))
    implementation_records = [
        _file_record(path, base=implementation_base) for path in implementation_files
    ]
    run_manifest = {
        "schema_version": "1.0",
        "started_at_utc": started.isoformat(),
        "completed_at_utc": datetime.now(UTC_COMPAT).isoformat(),
        "command_contract": "python -m yeastbridge_vs.crc_scan with explicit input paths",
        "python": sys.version,
        "platform": platform.platform(),
        "config_sha256": _sha256_file(config_source),
        "implementation_files": implementation_records,
        "implementation_collection_sha256": _collection_digest(implementation_records),
        "random_seed": int(config["state_baseline"]["seed"]),
        "virtual_edges_used_as_labels": False,
        "test_set_threshold_selection": False,
        "ranking_bonus": {target: 0.0 for target in sorted(CALIBRATION_TARGETS)},
    }
    _write_json(destination / "run_manifest.json", run_manifest)

    output_files = sorted(
        path
        for path in destination.iterdir()
        if path.is_file() and path.name != "result_manifest.json"
    )
    result_manifest = {
        "schema_version": "1.0",
        "files": [_file_record(path, base=destination) for path in output_files],
    }
    result_manifest["collection_sha256"] = _collection_digest(result_manifest["files"])
    _write_json(destination / "result_manifest.json", result_manifest)
    return summary


def _read_uniprot_snapshot(path: Path, family: str) -> list[dict[str, Any]]:
    _required_file(path, f"{family} UniProt snapshot")
    with path.open("r", encoding="utf-8", errors="strict", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "Gene Names" not in reader.fieldnames:
            raise ScanError(f"{path} lacks the Gene Names column")
        rows: list[dict[str, Any]] = []
        seen: set[str] = set()
        for line_number, raw in enumerate(reader, start=2):
            names = str(raw.get("Gene Names", "")).split()
            if not names:
                raise ScanError(f"{path}:{line_number} lacks a primary gene symbol")
            target = names[0]
            if target in seen:
                raise ScanError(f"{path} has duplicate primary symbol {target}")
            seen.add(target)
            entry = str(raw.get("Entry", "") or raw.get("Entry Name", "")).strip()
            rows.append(
                {
                    "target_id": target,
                    "target_family": family,
                    "membership_sources": f"UniProt_reviewed_human_keyword_snapshot:{path.name}",
                    "uniprot_entries": entry,
                    "gene_aliases": ";".join(names[1:]),
                    "protein_name": str(raw.get("Protein names", "")).strip(),
                    "protein_length": str(raw.get("Length", "")).strip(),
                    "keywords": str(raw.get("Keywords", "")).strip(),
                }
            )
    return rows


def _build_state_baseline(
    becker_h5ad: Path,
    metadata_path: Path,
    xatlas_gene_metadata: Path,
    target_ids: Sequence[str],
    config: Mapping[str, Any],
) -> tuple[dict[str, dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import sparse

    metadata = pd.read_csv(metadata_path, dtype=str).fillna("")
    required = {"RNA_SampleName", "DiseaseState", "Donor"}
    missing_columns = sorted(required - set(metadata.columns))
    if missing_columns:
        raise ScanError(f"official Becker metadata lacks columns: {', '.join(missing_columns)}")
    metadata = metadata[metadata["RNA_SampleName"].str.strip() != ""].copy()
    metadata["RNA_SampleName"] = metadata["RNA_SampleName"].str.strip()
    duplicate_names = sorted(
        metadata.loc[metadata["RNA_SampleName"].duplicated(keep=False), "RNA_SampleName"].unique()
    )
    if duplicate_names:
        raise ScanError(
            f"official Becker metadata has duplicate RNA sample names: {duplicate_names}"
        )
    by_sample = metadata.set_index("RNA_SampleName").to_dict("index")

    adata = ad.read_h5ad(becker_h5ad)
    if "sample" not in adata.obs:
        raise ScanError("Becker h5ad lacks obs['sample']")
    samples = adata.obs["sample"].astype(str).to_numpy()
    observed_samples = sorted(set(samples))
    unmapped = sorted(set(observed_samples) - set(by_sample))
    mapped_states = np.array(
        [str(by_sample.get(sample, {}).get("DiseaseState", "")) for sample in samples],
        dtype=object,
    )
    disease_values = {str(value) for value in config["disease_values"]}
    desired_values = {str(value) for value in config["desired_values"]}
    disease_cells = np.isin(mapped_states, list(disease_values))
    desired_cells = np.isin(mapped_states, list(desired_values))
    if not disease_cells.any() or not desired_cells.any():
        raise ScanError("official metadata produced an empty disease or desired cell group")

    included_samples = sorted(
        sample
        for sample in observed_samples
        if str(by_sample.get(sample, {}).get("DiseaseState", "")) in disease_values | desired_values
    )
    sample_means: list[np.ndarray] = []
    sample_states: list[str] = []
    sample_donors: list[str] = []
    for sample in included_samples:
        indices = np.flatnonzero(samples == sample)
        mean = np.asarray(adata.X[indices].mean(axis=0)).ravel()
        sample_means.append(mean.astype(np.float64, copy=False))
        sample_states.append(str(by_sample[sample]["DiseaseState"]))
        sample_donors.append(str(by_sample[sample]["Donor"]))
    sample_matrix = np.vstack(sample_means)
    disease_sample_mask = np.isin(sample_states, list(disease_values))
    desired_sample_mask = np.isin(sample_states, list(desired_values))
    if disease_sample_mask.sum() < 2 or desired_sample_mask.sum() < 2:
        raise ScanError("at least two samples per state are required")

    disease_mean = sample_matrix[disease_sample_mask].mean(axis=0)
    desired_mean = sample_matrix[desired_sample_mask].mean(axis=0)
    state_effect_all = disease_mean - desired_mean
    disease_var = sample_matrix[disease_sample_mask].var(axis=0, ddof=1)
    desired_var = sample_matrix[desired_sample_mask].var(axis=0, ddof=1)
    standard_error = np.sqrt(
        disease_var / disease_sample_mask.sum() + desired_var / desired_sample_mask.sum()
    )
    state_snr_all = np.abs(state_effect_all) / (standard_error + 1e-8)
    var_index = {str(gene): index for index, gene in enumerate(adata.var_names)}

    target_columns = [var_index[target] for target in target_ids if target in var_index]
    target_order = [target for target in target_ids if target in var_index]
    disease_detection = _detection_fraction(adata.X[disease_cells], target_columns, sparse)
    desired_detection = _detection_fraction(adata.X[desired_cells], target_columns, sparse)
    detection_by_target = {
        target: (float(disease_detection[index]), float(desired_detection[index]))
        for index, target in enumerate(target_order)
    }

    seed = int(config["seed"])
    replicates = int(config["bootstrap_replicates"])
    if replicates < 1:
        raise ScanError("bootstrap_replicates must be positive")
    target_present_indices = np.array([var_index[target] for target in target_order], dtype=int)
    observed_target_effect = state_effect_all[target_present_indices]
    sign_hits = np.zeros(len(target_order), dtype=np.int64)
    rng = np.random.default_rng(seed)
    disease_rows = sample_matrix[disease_sample_mask][:, target_present_indices]
    desired_rows = sample_matrix[desired_sample_mask][:, target_present_indices]
    for _ in range(replicates):
        boot_disease = disease_rows[
            rng.integers(0, disease_rows.shape[0], size=disease_rows.shape[0])
        ].mean(axis=0)
        boot_desired = desired_rows[
            rng.integers(0, desired_rows.shape[0], size=desired_rows.shape[0])
        ].mean(axis=0)
        boot_effect = boot_disease - boot_desired
        sign_hits += np.signbit(boot_effect) == np.signbit(observed_target_effect)
    sign_probability = sign_hits / replicates

    state: dict[str, dict[str, Any]] = {}
    present_position = {target: index for index, target in enumerate(target_order)}
    for target in target_ids:
        if target not in var_index:
            state[target] = {
                "state_observed": False,
                "intended_direction": "unresolved",
                "state_effect_disease_minus_desired": "",
                "state_signal_to_noise": "",
                "state_sign_probability": "",
                "disease_expression_fraction": 0.0,
                "desired_expression_fraction": 0.0,
            }
            continue
        column = var_index[target]
        position = present_position[target]
        effect = float(state_effect_all[column])
        direction = "inhibit" if effect > 0 else "activate" if effect < 0 else "unresolved"
        disease_fraction, desired_fraction = detection_by_target[target]
        state[target] = {
            "state_observed": True,
            "intended_direction": direction,
            "state_effect_disease_minus_desired": effect,
            "state_signal_to_noise": float(state_snr_all[column]),
            "state_sign_probability": float(sign_probability[position]),
            "disease_expression_fraction": disease_fraction,
            "desired_expression_fraction": desired_fraction,
        }

    gene_metadata = pd.read_parquet(xatlas_gene_metadata)
    if not {"gene_name", "gene_token_id"}.issubset(gene_metadata.columns):
        raise ScanError("X-Atlas gene metadata lacks gene_name/gene_token_id")
    gene_metadata = gene_metadata.sort_values("gene_token_id").drop_duplicates("gene_name")
    token_by_gene = dict(
        zip(gene_metadata["gene_name"], gene_metadata["gene_token_id"], strict=True)
    )
    excluded_prefixes = tuple(str(value) for value in config["excluded_signature_prefixes"])
    eligible_signature: list[tuple[str, int, float, float, float]] = []
    full_disease_detection = _detection_fraction(
        adata.X[disease_cells], list(range(adata.n_vars)), sparse
    )
    full_desired_detection = _detection_fraction(
        adata.X[desired_cells], list(range(adata.n_vars)), sparse
    )
    for index, gene_value in enumerate(adata.var_names):
        gene = str(gene_value)
        if gene not in token_by_gene or gene.startswith(excluded_prefixes):
            continue
        effect = float(state_effect_all[index])
        if not math.isfinite(effect) or effect == 0.0:
            continue
        eligible_signature.append(
            (
                gene,
                int(token_by_gene[gene]),
                effect,
                float(full_disease_detection[index]),
                float(full_desired_detection[index]),
            )
        )
    per_direction = int(config["signature_genes_per_direction"])
    positive = sorted(
        (row for row in eligible_signature if row[2] > 0),
        key=lambda row: (-row[2], row[0]),
    )[:per_direction]
    negative = sorted(
        (row for row in eligible_signature if row[2] < 0),
        key=lambda row: (row[2], row[0]),
    )[:per_direction]
    selected_signature = positive + negative
    if len(positive) != per_direction or len(negative) != per_direction:
        raise ScanError("insufficient signed genes for the frozen Becker state signature")
    signature = [
        {
            "gene": gene,
            "xatlas_token": token,
            "state_effect_disease_minus_desired": effect,
            "disease_expression_fraction": disease_fraction,
            "desired_expression_fraction": desired_fraction,
            "signature_direction": "disease_up" if effect > 0 else "disease_down",
            "selection_rule": f"top_{per_direction}_per_signed_direction",
        }
        for gene, token, effect, disease_fraction, desired_fraction in selected_signature
    ]
    state_counts = defaultdict(int)
    for sample in included_samples:
        state_counts[str(by_sample[sample]["DiseaseState"])] += 1
    old_heuristic_mismatches = [
        {
            "sample": sample,
            "official_disease_state": str(by_sample[sample]["DiseaseState"]),
            "official_binary_state": _official_binary_state(str(by_sample[sample]["DiseaseState"])),
            "old_filename_heuristic_state": _old_tissue_heuristic(sample),
        }
        for sample in included_samples
        if _old_tissue_heuristic(sample)
        != _official_binary_state(str(by_sample[sample]["DiseaseState"]))
    ]
    qc = {
        "h5ad_cells": int(adata.n_obs),
        "h5ad_genes": int(adata.n_vars),
        "observed_samples": len(observed_samples),
        "included_officially_labeled_samples": len(included_samples),
        "unmapped_samples": unmapped,
        "official_state_sample_counts": dict(sorted(state_counts.items())),
        "disease_cells": int(disease_cells.sum()),
        "desired_cells": int(desired_cells.sum()),
        "disease_samples": int(disease_sample_mask.sum()),
        "desired_samples": int(desired_sample_mask.sum()),
        "disease_donors": len(set(np.array(sample_donors)[disease_sample_mask])),
        "desired_donors": len(set(np.array(sample_donors)[desired_sample_mask])),
        "contrast_unit": "sample_pseudobulk",
        "bootstrap_unit": "sample_not_donor_cluster",
        "universe_targets_in_becker_matrix": len(target_order),
        "signature_genes": len(signature),
        "official_metadata_used": True,
        "old_filename_letter_heuristic_used": False,
        "old_filename_heuristic_sample_mismatches": len(old_heuristic_mismatches),
        "old_filename_heuristic_mismatch_records": old_heuristic_mismatches,
    }
    return state, signature, qc


def _scan_xatlas(
    files: Sequence[Path],
    gene_metadata_path: Path,
    target_ids: Sequence[str],
    signature: Sequence[Mapping[str, Any]],
    control_label: str,
) -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
    import numpy as np
    import pandas as pd
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    from scipy import sparse

    gene_metadata = pd.read_parquet(gene_metadata_path).sort_values("gene_token_id")
    gene_metadata = gene_metadata.drop_duplicates("gene_name")
    expected_tokens = [int(row["xatlas_token"]) for row in signature]
    observed_token_by_gene = dict(
        zip(gene_metadata["gene_name"], gene_metadata["gene_token_id"], strict=True)
    )
    for row in signature:
        gene = str(row["gene"])
        if int(observed_token_by_gene.get(gene, -1)) != int(row["xatlas_token"]):
            raise ScanError(f"X-Atlas token drift for signature gene {gene}")
    signature_vector = np.asarray(
        [float(row["state_effect_disease_minus_desired"]) for row in signature],
        dtype=np.float64,
    )
    signature_norm = np.linalg.norm(signature_vector)
    if not math.isfinite(float(signature_norm)) or signature_norm == 0:
        raise ScanError("Becker signature has zero or non-finite norm")
    disease_unit = signature_vector / signature_norm

    target_to_index = {target: index for index, target in enumerate(target_ids)}
    watch = pa.array([*target_ids, control_label])
    signature_tokens = pa.array(expected_tokens, type=pa.int64())
    max_token = max(expected_tokens)
    token_to_position = np.full(max_token + 1, -1, dtype=np.int32)
    token_to_position[np.asarray(expected_tokens, dtype=int)] = np.arange(len(expected_tokens))
    weighted_delta_sum = np.zeros((len(target_ids), len(signature)), dtype=np.float64)
    total_cells = np.zeros(len(target_ids), dtype=np.int64)
    batch_counts = np.zeros(len(target_ids), dtype=np.int64)
    per_batch_reversal: list[list[float]] = [[] for _ in target_ids]
    total_controls = 0
    used_batches = 0

    for path in files:
        table = pq.read_table(
            path,
            columns=[
                "gene_token_id",
                "gene_expression",
                "gene_target",
                "total_counts",
                "pass_guide_filter",
            ],
        )
        mask = pc.and_(
            pc.equal(table["pass_guide_filter"], 1),
            pc.is_in(table["gene_target"], value_set=watch),
        )
        selected = table.filter(mask)
        if selected.num_rows == 0:
            continue
        targets = np.asarray(selected["gene_target"].to_pylist(), dtype=object)
        controls = np.flatnonzero(targets == control_label)
        if controls.size == 0:
            continue
        totals = np.asarray(
            selected["total_counts"].to_numpy(zero_copy_only=False), dtype=np.float64
        )
        token_lists = selected["gene_token_id"].combine_chunks()
        expression_lists = selected["gene_expression"].combine_chunks()
        flat_tokens_arrow = pc.list_flatten(token_lists)
        token_mask = np.asarray(
            pc.is_in(flat_tokens_arrow, value_set=signature_tokens).to_numpy(zero_copy_only=False),
            dtype=bool,
        )
        parents = np.asarray(
            pc.list_parent_indices(token_lists).to_numpy(zero_copy_only=False), dtype=np.int64
        )[token_mask]
        tokens = np.asarray(flat_tokens_arrow.to_numpy(zero_copy_only=False), dtype=np.int64)[
            token_mask
        ]
        expression = np.asarray(
            pc.list_flatten(expression_lists).to_numpy(zero_copy_only=False), dtype=np.float64
        )[token_mask]
        positions = token_to_position[tokens]
        valid = (positions >= 0) & (totals[parents] > 0) & np.isfinite(expression)
        values = np.log1p(expression[valid] / totals[parents[valid]] * 10000.0)
        cell_matrix = sparse.coo_matrix(
            (values, (parents[valid], positions[valid])),
            shape=(selected.num_rows, len(signature)),
            dtype=np.float64,
        ).tocsr()
        control_mean = np.asarray(cell_matrix[controls].mean(axis=0)).ravel()
        total_controls += int(controls.size)

        candidate_rows = np.flatnonzero(targets != control_label)
        if candidate_rows.size == 0:
            continue
        global_target_indices = np.fromiter(
            (target_to_index[str(targets[row])] for row in candidate_rows),
            dtype=np.int64,
            count=candidate_rows.size,
        )
        group_matrix = sparse.coo_matrix(
            (
                np.ones(candidate_rows.size, dtype=np.float64),
                (global_target_indices, candidate_rows),
            ),
            shape=(len(target_ids), selected.num_rows),
        ).tocsr()
        local_counts = np.bincount(global_target_indices, minlength=len(target_ids))
        active = np.flatnonzero(local_counts > 0)
        sums = group_matrix[active] @ cell_matrix
        means = np.asarray(sums.toarray()) / local_counts[active, None]
        deltas = means - control_mean[None, :]
        reversals = -(deltas @ disease_unit)
        weighted_delta_sum[active] += deltas * local_counts[active, None]
        total_cells[active] += local_counts[active]
        batch_counts[active] += 1
        for target_index, score in zip(active, reversals, strict=True):
            per_batch_reversal[int(target_index)].append(float(score))
        used_batches += 1

    result: dict[str, dict[str, Any]] = {}
    for index, target in enumerate(target_ids):
        count = int(total_cells[index])
        if count == 0:
            result[target] = {
                "xatlas_observed": False,
                "xatlas_good_cells": 0,
                "xatlas_batches": 0,
                "xatlas_knockdown_reversal": "",
                "xatlas_knockdown_reversal_cosine": "",
                "xatlas_batch_positive_fraction": "",
                "xatlas_batch_score_se": "",
            }
            continue
        delta = weighted_delta_sum[index] / count
        reversal = float(-(delta @ disease_unit))
        delta_norm = float(np.linalg.norm(delta))
        cosine = reversal / (delta_norm + 1e-12)
        batch_scores = np.asarray(per_batch_reversal[index], dtype=np.float64)
        score_se = (
            float(batch_scores.std(ddof=1) / math.sqrt(batch_scores.size))
            if batch_scores.size > 1
            else float("nan")
        )
        result[target] = {
            "xatlas_observed": True,
            "xatlas_good_cells": count,
            "xatlas_batches": int(batch_counts[index]),
            "xatlas_knockdown_reversal": reversal,
            "xatlas_knockdown_reversal_cosine": float(cosine),
            "xatlas_batch_positive_fraction": float((batch_scores > 0).mean()),
            "xatlas_batch_score_se": score_se,
        }
    qc = {
        "files_discovered": len(files),
        "files_used": used_batches,
        "signature_genes": len(signature),
        "non_targeting_good_cells_summed_by_batch": total_controls,
        "targets_with_any_good_cell": int((total_cells > 0).sum()),
        "targets_n_ge_minimum_cells": 0,
        "normalization": "log1p(count/total_counts*10000), absent signature genes set to zero",
        "batch_adjustment": "within-batch target-minus-non-targeting pseudobulk",
        "measured_perturbation": True,
        "virtual_edge_labels_used": False,
    }
    return result, qc


def _audit_norman(
    path: Path,
    target_ids: set[str],
) -> tuple[dict[str, int], dict[str, Any]]:
    import anndata as ad

    adata = ad.read_h5ad(path, backed="r")
    if "condition" not in adata.obs:
        raise ScanError("Norman h5ad lacks obs['condition']")
    coverage: defaultdict[str, int] = defaultdict(int)
    combo_coverage: defaultdict[str, int] = defaultdict(int)
    for condition, count_value in adata.obs["condition"].value_counts().items():
        count = int(count_value)
        genes = [gene for gene in str(condition).split("+") if gene != "ctrl"]
        for gene in genes:
            if gene in target_ids:
                combo_coverage[gene] += count
        if len(genes) == 1 and genes[0] in target_ids:
            coverage[genes[0]] += count
    return dict(coverage), {
        "cells": int(adata.n_obs),
        "conditions": int(adata.obs["condition"].nunique()),
        "cell_type_values": sorted(map(str, adata.obs["cell_type"].unique()))
        if "cell_type" in adata.obs
        else [],
        "clean_single_target_count": len(coverage),
        "clean_single_targets": dict(sorted(coverage.items())),
        "any_combo_target_count": len(combo_coverage),
        "used_for_global_ranking": False,
        "reason_not_ranked": (
            "Only clean single-target coverage is audited; A549 CRISPRa is not a CRC-context "
            "ground truth and covers too little of the target universe."
        ),
    }


def _assemble_candidates(
    universe: Sequence[Mapping[str, Any]],
    state: Mapping[str, Mapping[str, Any]],
    perturbation: Mapping[str, Mapping[str, Any]],
    norman_coverage: Mapping[str, int],
    orthologs: Mapping[str, Sequence[str]],
    config: Mapping[str, Any],
) -> list[dict[str, Any]]:
    state_config = config["state_baseline"]
    perturb_config = config["perturbation_baseline"]
    min_expression = float(state_config["minimum_expression_fraction"])
    min_sign_probability = float(state_config["minimum_state_sign_probability"])
    min_cells = int(perturb_config["minimum_good_perturbed_cells"])
    min_batches = int(perturb_config["minimum_batches"])
    min_batch_fraction = float(perturb_config["minimum_batch_direction_fraction"])
    rows: list[dict[str, Any]] = []
    for universe_row in universe:
        target = str(universe_row["target_id"])
        row = dict(universe_row)
        row.update(state[target])
        row.update(perturbation[target])
        direction = str(row["intended_direction"])
        knockdown_score = row["xatlas_knockdown_reversal"]
        positive_fraction = row["xatlas_batch_positive_fraction"]
        if direction == "inhibit" and _finite_number(knockdown_score):
            directional_score = float(knockdown_score)
            directional_batch_fraction = float(positive_fraction)
            xatlas_implied = "inhibit" if float(knockdown_score) > 0 else "activate"
        elif direction == "activate" and _finite_number(knockdown_score):
            directional_score = -float(knockdown_score)
            directional_batch_fraction = 1.0 - float(positive_fraction)
            xatlas_implied = "activate" if float(knockdown_score) < 0 else "inhibit"
        else:
            directional_score = ""
            directional_batch_fraction = ""
            xatlas_implied = "unresolved"
        row["directional_xatlas_reversal"] = directional_score
        row["xatlas_batch_direction_fraction"] = directional_batch_fraction
        row["xatlas_implied_direction"] = xatlas_implied
        row["state_xatlas_direction_agreement"] = (
            direction == xatlas_implied if direction in SIGNED_DIRECTIONS else False
        )
        direct, evidence_class, interpretation = _direction_evidence(direction)
        row["signed_direction_directly_perturbed"] = direct
        row["direction_evidence_class"] = evidence_class
        row["direction_evidence_interpretation"] = interpretation
        row["norman_clean_crispra_cells"] = int(norman_coverage.get(target, 0))
        target_orthologs = sorted(set(map(str, orthologs.get(target, ()))))
        row["yeast_orthologs"] = ";".join(target_orthologs)
        yeast_status, yeast_route, yeast_reason = _yeast_route(row, target_orthologs)
        row["yeast_route_status"] = yeast_status
        row["yeast_validation_route"] = yeast_route
        row["yeast_nonexecutability_reason"] = yeast_reason
        row["mammalian_validation_route"] = _mammalian_route(
            direction, int(row["norman_clean_crispra_cells"])
        )

        expression_ok = (
            max(
                float(row["disease_expression_fraction"]),
                float(row["desired_expression_fraction"]),
            )
            >= min_expression
        )
        state_ok = (
            bool(row["state_observed"])
            and direction in SIGNED_DIRECTIONS
            and _finite_number(row["state_sign_probability"])
            and float(row["state_sign_probability"]) >= min_sign_probability
        )
        xatlas_coverage_ok = (
            bool(row["xatlas_observed"])
            and int(row["xatlas_good_cells"]) >= min_cells
            and int(row["xatlas_batches"]) >= min_batches
        )
        perturbation_direction_ok = (
            _finite_number(directional_score)
            and float(directional_score) > 0
            and _finite_number(directional_batch_fraction)
            and float(directional_batch_fraction) >= min_batch_fraction
        )
        row["expression_gate"] = expression_ok
        row["state_direction_gate"] = state_ok
        row["xatlas_coverage_gate"] = xatlas_coverage_ok
        row["xatlas_direction_gate"] = perturbation_direction_ok
        row["baseline_eligible"] = bool(
            expression_ok and state_ok and xatlas_coverage_ok and perturbation_direction_ok
        )
        uncertainty: list[str] = []
        if not row["state_observed"]:
            uncertainty.append("target_absent_from_becker_matrix")
        if not expression_ok:
            uncertainty.append("low_expression_in_becker_states")
        if not state_ok:
            uncertainty.append("state_direction_not_stable")
        if not xatlas_coverage_ok:
            uncertainty.append("insufficient_xatlas_crispri_coverage")
        if xatlas_coverage_ok and not perturbation_direction_ok:
            uncertainty.append("xatlas_does_not_stably_support_state_implied_direction")
        if direction == "activate" and not norman_coverage.get(target):
            uncertainty.append("activation_not_directly_tested_in_crc_context")
        uncertainty.append("becker_h5ad_is_40000_cell_heuristic_epithelial_subsample")
        uncertainty.append("becker_bootstrap_is_sample_level_not_donor_clustered")
        uncertainty.append("xatlas_is_hct116_transcriptomic_not_drug_efficacy")
        uncertainty.append("yeast_route_is_proposed_not_experimentally_confirmed")
        row["uncertainty"] = ";".join(uncertainty)
        row["human_evidence"] = (
            "official-metadata Becker sample-pseudobulk state contrast; measured within-batch "
            "X-Atlas HCT116 CRISPRi transcriptomic response; "
            + str(row["direction_evidence_class"])
        )
        row["ranking_bonus"] = 0.0
        rows.append(row)
    return rows


def _direction_evidence(direction: str) -> tuple[bool, str, str]:
    if direction == "inhibit":
        return (
            True,
            "direct_same_direction_crispri",
            (
                "X-Atlas directly measured target knockdown; the proposed inhibitory direction "
                "therefore has same-direction perturbation evidence, outside a prospective CRC "
                "validation setting."
            ),
        )
    if direction == "activate":
        return (
            False,
            "opposite_direction_crispri_counterfactual_requires_crispra",
            (
                "X-Atlas measured knockdown, not activation. The proposed activating direction "
                "assumes a reversible response and must be tested by CRISPRa or inducible cDNA."
            ),
        )
    return False, "unresolved", "No signed perturbation interpretation is available."


def _yeast_route(
    row: Mapping[str, Any],
    orthologs: Sequence[str],
) -> tuple[str, str, str]:
    family = str(row["target_family"])
    protein = str(row.get("protein_name", "")).lower()
    keywords = str(row.get("keywords", "")).lower()
    if orthologs:
        return (
            "native_proxy_available_but_not_target_pharmacology",
            (
                "Perturb mapped S. cerevisiae ortholog(s) for conserved-network direction, then "
                "express the human target for molecule-specific confirmation."
            ),
            (
                "Native ortholog phenotypes cannot establish binding or modulation of the "
                "human protein."
            ),
        )
    if "gpcr" in family:
        return (
            "conditional_humanized_gpcr",
            (
                "Heterologously express the human GPCR in S. cerevisiae, match it to a chimeric "
                "G-alpha/reporter circuit, and measure signed agonist/antagonist response."
            ),
            (
                "No mapped native yeast ortholog; surface trafficking, ligand access and G-protein "
                "coupling must pass construct-level QC before this route is executable."
            ),
        )
    multicomponent = any(
        token in protein or token in keywords
        for token in ("subunit", "heteromer", "ligand-gated", "regulatory")
    )
    if multicomponent:
        return (
            "conditional_humanized_multicomponent_channel",
            (
                "Reconstitute the human pore-forming complex and required accessory subunits in "
                "yeast; quantify membrane trafficking and ion-selective or membrane-potential "
                "readout."
            ),
            (
                "No mapped native yeast ortholog and the snapshot indicates a multicomponent or "
                "regulatory channel; a single-gene yeast assay would be non-executable."
            ),
        )
    return (
        "conditional_humanized_channel",
        (
            "Express the human channel in yeast and couple it to membrane-potential, ion-flux or "
            "growth-rescue readouts under a predeclared ionic condition."
        ),
        (
            "No mapped native yeast ortholog; membrane trafficking, orientation and functional "
            "conductance must be demonstrated before compound screening."
        ),
    )


def _mammalian_route(direction: str, norman_cells: int) -> str:
    if direction == "inhibit":
        return (
            "Repeat CRISPRi with independent guides and rescue in HCT116, then validate inhibition "
            "in a second CRC line and patient-derived organoid with target-proximal functional "
            "readout."
        )
    if direction == "activate":
        extra = (
            f" Norman A549 contains {norman_cells} clean single-target CRISPRa cells as "
            "cross-context "
            "support only."
            if norman_cells
            else ""
        )
        return (
            "Use CRISPRa or inducible cDNA in HCT116 and an independent CRC organoid; require "
            "rescue/reversal and target-proximal functional readout because X-Atlas CRISPRi tests "
            "the opposite "
            f"direction.{extra}"
        )
    return "Resolve the signed direction before selecting a mammalian perturbation assay."


def _load_orthologs(path: Path) -> dict[str, list[str]]:
    result: defaultdict[str, set[str]] = defaultdict(set)
    with path.open("r", encoding="utf-8-sig", errors="strict", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"human_symbol", "yeast_sgd"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ScanError("ortholog map lacks human_symbol/yeast_sgd")
        for row in reader:
            human = str(row["human_symbol"]).strip()
            yeast = str(row["yeast_sgd"]).strip()
            if human and yeast:
                result[human].add(yeast)
    return {human: sorted(genes) for human, genes in result.items()}


def _load_config(path: Path) -> dict[str, Any]:
    _required_file(path, "CRC scan config")
    try:
        config = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ScanError(f"cannot parse config {path}: {exc}") from exc
    if not isinstance(config, dict):
        raise ScanError("config root must be an object")
    if config.get("selection_mode") != "target_agnostic":
        raise ScanError("selection_mode must be target_agnostic")
    query = config.get("query")
    required_query = {"disease_context_id", "cell_type_id", "disease_state", "desired_state"}
    if not isinstance(query, dict) or not required_query.issubset(query):
        raise ScanError("query must define disease context, cell type and signed state transition")
    controls = config.get("calibration_controls", {})
    if set(controls.get("target_ids", [])) != CALIBRATION_TARGETS:
        raise ScanError("LPAR1 and KCNQ2 must be calibration controls")
    if float(controls.get("ranking_bonus", math.nan)) != 0.0:
        raise ScanError("calibration-control ranking bonus must be zero")
    if config.get("ranking", {}).get("single_hidden_total_score_forbidden") is not True:
        raise ScanError("hidden total scores must be forbidden")
    return config


def _input_manifest(
    paths: Sequence[Path],
    *,
    xatlas_root: Path,
    workers: int,
    hash_inputs: bool,
) -> dict[str, Any]:
    unique = list(dict.fromkeys(path.resolve() for path in paths))
    if not hash_inputs:
        files = [
            {
                "path": _display_input_path(path, xatlas_root),
                "bytes": path.stat().st_size,
                "sha256": None,
            }
            for path in unique
        ]
        return {
            "schema_version": "1.0",
            "content_addressed": False,
            "files": files,
            "collection_sha256": None,
            "blocker": "input hashes explicitly skipped",
        }
    workers = max(1, int(workers))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        digests = list(pool.map(_sha256_file, unique))
    files = [
        {
            "path": _display_input_path(path, xatlas_root),
            "bytes": path.stat().st_size,
            "sha256": digest,
        }
        for path, digest in zip(unique, digests, strict=True)
    ]
    return {
        "schema_version": "1.0",
        "content_addressed": True,
        "files": files,
        "collection_sha256": _collection_digest(files),
    }


def _blockers(
    config: Mapping[str, Any],
    state_qc: Mapping[str, Any],
    xatlas_qc: Mapping[str, Any],
    norman_qc: Mapping[str, Any],
    rows: Sequence[Mapping[str, Any]],
    hash_inputs: bool,
) -> list[dict[str, str]]:
    blockers = [
        {
            "severity": "high",
            "id": "FAMILY-INDEPENDENT-RECONCILIATION",
            "detail": (
                "The frozen complete universe is a verifiable UniProt reviewed keyword snapshot, "
                "but an independent HGNC/IUPHAR reconciliation snapshot is not yet present."
            ),
        },
        {
            "severity": "high",
            "id": "BECKER-CELL-ANNOTATION",
            "detail": (
                "The reused Becker h5ad is a 40,000-cell subsample whose epithelial gate was the "
                "upper half of a marker score, not the authors' released cell annotations. "
                "Official "
                "sample metadata fixes disease labels but not this cell-level gate."
            ),
        },
        {
            "severity": "high",
            "id": "BECKER-DONOR-CLUSTERING",
            "detail": (
                f"The state contrast uses {state_qc['disease_samples']} disease samples from "
                f"{state_qc['disease_donors']} donors and {state_qc['desired_samples']} desired "
                f"samples from {state_qc['desired_donors']} donors. Its bootstrap resamples "
                "samples rather than donor clusters, so state sign probabilities may be "
                "anti-conservative. Refit a donor-aware pseudobulk model before target claims."
            ),
        },
        {
            "severity": "high",
            "id": "NO-PROSPECTIVE-GROUND-TRUTH",
            "detail": (
                "No held-out or prospective target-validation labels exist. Pareto fronts are "
                "prioritization output and cannot be reported as predictive accuracy."
            ),
        },
        {
            "severity": "high",
            "id": "NO-PAIRED-YEAST-EFFECT",
            "detail": (
                "Yeast routes are construct proposals only; no paired human/yeast perturbation or "
                "compound result is used in this scan."
            ),
        },
        {
            "severity": "medium",
            "id": "XATLAS-CONTEXT",
            "detail": (
                "X-Atlas is measured HCT116 CRISPRi transcriptomics, but state reversal is not "
                "viability, binding, efficacy, selectivity or an organoid endpoint."
            ),
        },
        {
            "severity": "medium",
            "id": "CRISPRA-COVERAGE",
            "detail": (
                f"Norman provides clean single-target CRISPRa coverage for only "
                f"{norman_qc['clean_single_target_count']} universe target(s), in A549 rather "
                "than CRC."
            ),
        },
    ]
    eligible_activation_count = sum(
        bool(row.get("baseline_eligible")) and row.get("intended_direction") == "activate"
        for row in rows
    )
    if eligible_activation_count:
        blockers.append(
            {
                "severity": "high",
                "id": "ACTIVATION-NOT-DIRECTLY-PERTURBED",
                "detail": (
                    f"{eligible_activation_count} eligible activation-direction rows are inferred "
                    "from the opposite CRISPRi response. They are counterfactual hypotheses, not "
                    "direct activation evidence, until tested by CRC-context CRISPRa or inducible "
                    "cDNA."
                ),
            }
        )
    if state_qc.get("unmapped_samples"):
        blockers.append(
            {
                "severity": "high",
                "id": "BECKER-UNMAPPED-SAMPLES",
                "detail": f"Unmapped Becker samples: {state_qc['unmapped_samples']}",
            }
        )
    if not hash_inputs:
        blockers.append(
            {
                "severity": "critical",
                "id": "INPUTS-NOT-CONTENT-ADDRESSED",
                "detail": "Production interpretation is blocked because input hashes were skipped.",
            }
        )
    if config["universe"].get("independent_hgnc_iuphar_reconciliation_status") != "ready":
        pass
    if xatlas_qc.get("virtual_edge_labels_used") is not False:
        raise ScanError("virtual edges must not be used as labels")
    return blockers


def _calibration_audit(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    by_target = {str(row["target_id"]): row for row in rows}
    details: dict[str, Any] = {}
    for target in sorted(CALIBRATION_TARGETS):
        row = by_target.get(target)
        details[target] = {
            "present_in_complete_universe": row is not None,
            "ranking_bonus": float(row["ranking_bonus"]) if row else None,
            "selected_by_general_rules_only": True,
            "special_case_feature_or_score": False,
            "pareto_front": row.get("pareto_front", "") if row else "",
        }
    return {
        "passed": all(
            item["present_in_complete_universe"]
            and item["ranking_bonus"] == 0.0
            and item["special_case_feature_or_score"] is False
            for item in details.values()
        ),
        "targets": details,
    }


def _universe_counts(universe: Sequence[Mapping[str, Any]]) -> dict[str, int]:
    gpcr = sum("gpcr" in str(row["target_family"]).split(";") for row in universe)
    channel = sum("ion_channel" in str(row["target_family"]).split(";") for row in universe)
    return {
        "gpcr_snapshot_rows": gpcr,
        "ion_channel_snapshot_rows": channel,
        "unique_primary_symbols": len(universe),
        "cross_family_symbol_overlap": gpcr + channel - len(universe),
    }


def _detection_fraction(matrix: Any, columns: Sequence[int], sparse_module: Any) -> Any:
    import numpy as np

    if not columns:
        return np.asarray([], dtype=np.float64)
    selected = matrix[:, columns]
    if sparse_module.issparse(selected):
        return np.asarray(selected.getnnz(axis=0), dtype=np.float64) / selected.shape[0]
    return np.asarray((np.asarray(selected) > 0).mean(axis=0), dtype=np.float64).ravel()


def _old_tissue_heuristic(sample: str) -> str:
    if sample.startswith("CRC") or sample.startswith("F"):
        return "disease"
    parts = sample.split("-")
    letter = parts[1] if len(parts) > 1 else ""
    return "desired" if letter == "A" else "disease" if letter == "C" else "excluded"


def _official_binary_state(state: str) -> str:
    if state in {"Normal", "Unaffected"}:
        return "desired"
    if state in {"Polyp", "Adenocarcinoma"}:
        return "disease"
    return "excluded"


def _dominates(left: Mapping[str, Any], right: Mapping[str, Any], axes: Sequence[str]) -> bool:
    left_values = [float(left[axis]) for axis in axes]
    right_values = [float(right[axis]) for axis in axes]
    return all(a >= b for a, b in zip(left_values, right_values, strict=True)) and any(
        a > b for a, b in zip(left_values, right_values, strict=True)
    )


def _finite_number(value: Any) -> bool:
    return (
        not isinstance(value, bool)
        and isinstance(value, int | float)
        and math.isfinite(float(value))
    )


def _resolve_project_input(root: Path, value: Any) -> Path:
    path = Path(str(value))
    resolved = path.resolve() if path.is_absolute() else (root / path).resolve()
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise ScanError(f"configured input escapes project root: {value}") from exc
    return _required_file(resolved, "configured project input")


def _required_file(value: str | Path, label: str) -> Path:
    path = Path(value).resolve()
    if not path.is_file():
        raise ScanError(f"missing {label}: {path}")
    return path


def _batch_sort_key(path: Path) -> tuple[int, str]:
    stem = path.stem
    suffix = stem.removeprefix("HCT116_Batch")
    return (int(suffix) if suffix.isdigit() else 10**9, path.name)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _file_record(path: Path, *, base: Path) -> dict[str, Any]:
    return {
        "path": path.resolve().relative_to(base.resolve()).as_posix(),
        "bytes": path.stat().st_size,
        "sha256": _sha256_file(path),
    }


def _collection_digest(files: Sequence[Mapping[str, Any]]) -> str:
    payload = json.dumps(
        [dict(record) for record in files],
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _display_input_path(path: Path, xatlas_root: Path) -> str:
    try:
        relative = path.resolve().relative_to(xatlas_root.resolve())
    except ValueError:
        return str(path.resolve())
    return f"xatlas://{relative.as_posix()}"


def _ordered_union_fields(rows: Sequence[Mapping[str, Any]]) -> list[str]:
    fields: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for field in row:
            if field not in seen:
                seen.add(field)
                fields.append(field)
    return fields


def _write_tsv(path: Path, rows: Sequence[Mapping[str, Any]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="", delete=False, dir=path.parent, prefix=f".{path.name}."
    ) as handle:
        temporary = Path(handle.name)
        writer = csv.DictWriter(
            handle, fieldnames=list(fields), delimiter="\t", extrasaction="raise"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def _write_json(path: Path, value: Any) -> None:
    payload = (json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n").encode(
        "utf-8"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "wb", delete=False, dir=path.parent, prefix=f".{path.name}."
    ) as handle:
        temporary = Path(handle.name)
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)
