"""Yeast-first compound discovery from human functional intent and ChemGRID."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import urllib.parse
import urllib.request
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import numpy as np

from .workspace import exclusive_json_write, project_root, resolve_in_root, sha256_file


class YeastFirstError(ValueError):
    """Raised when a yeast-first release cannot be audited."""


def _hash(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _load_json(root: Path, path: str) -> dict[str, Any]:
    file = resolve_in_root(root, path, must_exist=True)
    with file.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise YeastFirstError(f"JSON document must be an object: {file}")
    return value


def _load_intent(root: Path, config: dict[str, Any]) -> dict[str, Any]:
    path = resolve_in_root(root, config["human_intent"]["signature_path"], must_exist=True)
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise YeastFirstError("human intent signature is empty")
    weights = {
        str(row["gene"]): float(row["state_effect_disease_minus_desired"])
        for row in rows
    }

    # Dual-track signature (predeclared fusion rule): the three-foundation-model
    # direction scores (crc_model_evidence signature_model_scores.tsv) gate the
    # statistical weights. Genes whose fusion direction agrees with the
    # statistical sign keep full weight; disagreeing genes are halved; genes
    # without a model score keep 75% (conservative, predeclared).
    fusion_meta: dict[str, Any] = {"applied": False}
    scan_dir = root / "reports" / "crc_target_scan" / "model_evidence_v2"
    fusion_file = scan_dir / "signature_model_scores.tsv"
    if fusion_file.exists():
        fused = {}
        with fusion_file.open(encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                fused[str(row["gene"])] = row
        n_agree = n_disagree = n_missing = 0
        for gene, stat in weights.items():
            entry = fused.get(gene)
            if entry is None or entry.get("model_direction_score") in (None, "", "None"):
                weights[gene] = stat * 0.75
                n_missing += 1
                continue
            m = float(entry["model_direction_score"])
            if np.sign(m) == np.sign(stat):
                n_agree += 1
            else:
                weights[gene] = stat * 0.5
                n_disagree += 1
        fusion_meta = {
            "applied": True,
            "rule": "agree x1.0 / disagree x0.5 / no-model x0.75 (predeclared)",
            "n_agree": n_agree,
            "n_disagree": n_disagree,
            "n_no_model": n_missing,
        }
        print(f"intent fusion: {fusion_meta}", flush=True)

    return {
        "weights": weights,
        "source_path": str(path.relative_to(root)),
        "source_sha256": sha256_file(path),
        "dual_track_fusion": fusion_meta,
    }


def _load_orthology(root: Path, path: str) -> dict[str, set[str]]:
    file = resolve_in_root(root, path, must_exist=True)
    output: dict[str, set[str]] = {}
    with file.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            human = str(row.get("human_symbol", "")).strip()
            yeast = str(row.get("yeast_sgd", "")).strip()
            if human and yeast:
                output.setdefault(yeast, set()).add(human)
    return output


def _load_sga_neighbors(root: Path, path: str) -> dict[str, set[str]]:
    file = resolve_in_root(root, path, must_exist=True)
    neighbors: dict[str, set[str]] = {}
    import gzip

    with gzip.open(file, "rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            query, array = str(row.get("query_orf", "")), str(row.get("array_orf", ""))
            if query and array:
                neighbors.setdefault(query, set()).add(array)
                neighbors.setdefault(array, set()).add(query)
    return neighbors


def _propagate_personalized_pagerank(
    nodes: list[str],
    neighbors: dict[str, set[str]],
    anchor_weights: dict[str, float],
    alpha: float = 0.30,
    tolerance: float = 1e-10,
    max_iterations: int = 300,
) -> tuple[np.ndarray, int]:
    """Signed personalised PageRank over the undirected SGA graph.

    Anchor mass restarts every iteration (p = alpha*anchor + (1-alpha)*W^T p),
    so unlike repeated neighbour-mean diffusion the signal cannot wash out to
    uniformity on a dense graph. Positive and negative anchor masses propagate
    as separate channels and recombine with their signs, preserving
    intervention direction.
    """

    import scipy.sparse as sp

    index_of = {node: i for i, node in enumerate(nodes)}
    size = len(nodes)
    rows, cols = [], []
    seen_pairs: set[tuple[int, int]] = set()
    for orf, connected in neighbors.items():
        i = index_of.get(orf)
        if i is None:
            continue
        for other in connected:
            j = index_of.get(other)
            if j is None or (i, j) in seen_pairs:
                continue
            seen_pairs.add((i, j))
            rows.append(i)
            cols.append(j)
    data = np.ones(len(rows), dtype=np.float64)
    adjacency = sp.coo_matrix((data, (rows, cols)), shape=(size, size)).tocsr()
    row_sums = np.asarray(adjacency.sum(axis=1)).ravel()
    walk = sp.diags(1.0 / np.maximum(row_sums, 1.0)) @ adjacency

    def run_channel(mass: np.ndarray) -> tuple[np.ndarray, int]:
        restart = alpha * mass
        value = restart.copy()
        iterations = max_iterations
        for step in range(1, max_iterations + 1):
            nxt = restart + (1.0 - alpha) * (walk.T @ value)
            delta = float(np.abs(nxt - value).max())
            value = nxt
            if delta < tolerance:
                iterations = step
                break
        return value, iterations

    def channel_vector(channel: dict[str, float]) -> np.ndarray:
        mass = np.zeros(size)
        total = sum(abs(w) for w in channel.values())
        if total <= 0.0:
            return mass
        for orf, w in channel.items():
            i = index_of.get(orf)
            if i is not None:
                mass[i] = abs(w) / total
        return mass

    positives = {o: w for o, w in anchor_weights.items() if w > 0}
    negatives = {o: w for o, w in anchor_weights.items() if w < 0}
    pos_score, iter_pos = run_channel(channel_vector(positives))
    neg_score, iter_neg = run_channel(channel_vector(negatives))

    pos_scale = (
        sum(w for w in positives.values()) / len(positives) if positives else 0.0
    )
    neg_scale = (
        abs(sum(w for w in negatives.values())) / len(negatives) if negatives else 0.0
    )
    return pos_scale * pos_score - neg_scale * neg_score, max(iter_pos, iter_neg)


def _feature_weights(
    strain_orfs: list[str],
    orthology: dict[str, set[str]],
    weights: dict[str, float],
    neighbors: dict[str, set[str]],
    alpha: float = 0.30,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Translate human intent through the full SGA graph to ChemGRID strains."""

    anchor_weights: dict[str, float] = {}
    for orf, genes in orthology.items():
        human_weights = [weights[gene] for gene in genes if gene in weights]
        if human_weights:
            anchor_weights[orf] = -float(np.mean(human_weights))
    direct = np.array([anchor_weights.get(orf, 0.0) for orf in strain_orfs])

    all_nodes: set[str] = set(strain_orfs)
    for orf, connected in neighbors.items():
        all_nodes.add(orf)
        all_nodes.update(connected)
    nodes = sorted(all_nodes)
    node_index = {node: i for i, node in enumerate(nodes)}
    propagated_all, iterations = _propagate_personalized_pagerank(
        nodes, neighbors, anchor_weights, alpha=alpha
    )
    translated = np.array(
        [
            propagated_all[node_index[orf]] if orf in node_index else 0.0
            for orf in strain_orfs
        ]
    )

    abs_values = np.sort(np.abs(translated))
    n = abs_values.size
    total = float(abs_values.sum())
    if n and total > 0.0:
        cumulative = np.cumsum(abs_values) / total
        gini = float(1.0 - 2.0 * float((cumulative - np.arange(1, n + 1) / n).sum()) / n + 1.0 / n)
    else:
        gini = 0.0

    return translated, {
        "direct_anchor_count": len(anchor_weights),
        "anchor_human_gene_count": len(
            {
                gene
                for genes in orthology.values()
                for gene in genes
                if gene in weights
            }
        ),
        "chemgrid_direct_anchor_count": int(np.count_nonzero(direct)),
        "chemgrid_translated_count": int(np.count_nonzero(translated)),
        "sga_nodes_with_intent_signal": int(np.count_nonzero(propagated_all)),
        "translation_method": "signed_personalized_pagerank",
        "propagation_alpha": float(alpha),
        "concentration_gini_abs_chemgrid": round(gini, 6),
        "pagerank_iterations": iterations,
        "translation_rule": (
            f"human intent anchors propagated as signed personalised PageRank "
            f"(alpha={alpha}) over the SGA graph; per-iteration restart keeps "
            f"anchor mass instead of washing out"
        ),
    }


def _classify_target_name(name: str, config: Mapping[str, Any]) -> str | None:
    lowered = name.lower()
    rules = config.get("target_backmapping", {}).get("family_rules", {})
    if not isinstance(rules, Mapping):
        return None
    for family, tokens in rules.items():
        if any(str(token).lower() in lowered for token in tokens):
            return str(family)
    return None


def _chembl_backmap(
    inchikeys: list[str], config: Mapping[str, Any]
) -> dict[str, list[dict[str, Any]]]:
    """Resolve selected compounds to observed human target relations."""

    api = config.get("chembl_api", {})
    base_url = str(api.get("base_url", "https://www.ebi.ac.uk/chembl/api/data"))
    timeout = float(api.get("timeout_seconds", 60))
    output: dict[str, list[dict[str, Any]]] = {}
    for inchikey in inchikeys:
        molecule_params = {
            "molecule_structures__standard_inchi_key": inchikey,
            "limit": 20,
        }
        molecule_url = (
            f"{base_url.rstrip('/')}/molecule.json?"
            f"{urllib.parse.urlencode(molecule_params)}"
        )
        try:
            request = urllib.request.Request(molecule_url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(request, timeout=timeout) as response:
                molecule_document = json.load(response)
        except Exception:
            continue
        molecules = molecule_document.get("molecules", [])
        if not molecules:
            output[inchikey] = []
            continue
        molecule_ids = {
            str(row.get("molecule_chembl_id"))
            for row in molecules
            if row.get("molecule_chembl_id")
        }
        relations: dict[tuple[str, str], dict[str, Any]] = {}
        for molecule_id in molecule_ids:
            activity_params = {"molecule_chembl_id": molecule_id, "limit": 1000}
            activity_url = (
                f"{base_url.rstrip('/')}/activity.json?"
                f"{urllib.parse.urlencode(activity_params)}"
            )
            try:
                request = urllib.request.Request(
                    activity_url, headers={"Accept": "application/json"}
                )
                with urllib.request.urlopen(request, timeout=timeout) as response:
                    activity_document = json.load(response)
            except Exception:
                continue
            for activity in activity_document.get("activities", []):
                target_id = str(activity.get("target_chembl_id") or "").strip()
                target_name = str(activity.get("target_pref_name") or "").strip()
                if not target_id or str(activity.get("target_tax_id") or "") != "9606":
                    continue
                if target_name.upper() in {
                    "UNCHECKED",
                    "ADMET",
                    "NO RELEVANT TARGET",
                    "NON-PROTEIN TARGET",
                }:
                    continue
                key = (target_id, target_name)
                target_family = _classify_target_name(target_name, config)
                relation = relations.setdefault(
                    key,
                    {
                        "target_chembl_id": target_id,
                        "target_pref_name": target_name,
                        "target_organism": activity.get("target_organism"),
                        "target_family": target_family,
                        "action_types": set(),
                        "activity_count": 0,
                        "evidence_record_ids": [],
                    },
                )
                action = activity.get("action_type") or {}
                action_name = action.get("action_type") if isinstance(action, Mapping) else None
                if action_name:
                    relation["action_types"].add(str(action_name))
                relation["activity_count"] += 1
                relation["evidence_record_ids"].append(
                    f"chembl:{target_id}:{activity.get('assay_chembl_id')}"
                    f":{activity.get('activity_id')}"
                )
        normalized = []
        for relation in relations.values():
            relation["action_types"] = sorted(relation["action_types"])
            relation["evidence_record_ids"] = sorted(set(relation["evidence_record_ids"]))
            normalized.append(relation)
        normalized.sort(
            key=lambda row: (-int(row["activity_count"]), row["target_chembl_id"])
        )
        output[inchikey] = normalized
    return output


def build_yeast_first(
    root: str | Path,
    *,
    config_path: str | Path,
) -> dict[str, Any]:
    base = project_root(root)
    config = _load_json(base, config_path)
    intent = _load_intent(base, config)
    orthology = _load_orthology(base, config["human_intent"]["orthology_path"])
    neighbors = _load_sga_neighbors(base, config["human_intent"]["sga_path"])
    matrix = resolve_in_root(base, config["yeast_task"]["matrix_path"], must_exist=True)
    split_file = resolve_in_root(base, config["yeast_task"]["split_path"], must_exist=True)
    with np.load(matrix, allow_pickle=False) as data:
        inchikeys = data["compound_inchikeys"].astype(str)
        smiles = data["compound_smiles"].astype(str)
        z = data["z_score"].astype(np.float64)
        mask = data["measured_mask"].astype(bool) & np.isfinite(z)
        orfs = data["strain_orfs"].astype(str).tolist()
    translated, translation_meta = _feature_weights(orfs, orthology, intent["weights"], neighbors)
    norm = float(np.linalg.norm(translated))
    if norm == 0:
        raise YeastFirstError("human intent could not be translated to yeast features")
    translated /= norm

    # scYeast dual-base roles (design contract): (A) yeast-state representation
    # for intent-response alignment in embedding space, (B) confounder
    # (general toxicity / membrane / mitochondrial stress) representation.
    # Frozen checkpoint gene embeddings; both scores are reported alongside the
    # raw-z scores, and the v2 score drives selection (predeclared).
    embed_meta: dict[str, Any] = {"applied": False}
    strain_emb = None
    scyeast_file = root if isinstance(root, Path) else Path(root)
    scyeast_npz = project_root(root) / "reports" / "program_bridge_v2" / "scyeast_embeddings.npz"
    if scyeast_npz.exists():
        with np.load(scyeast_npz, allow_pickle=False) as data:
            emb_orfs = data["orfs"].astype(str)
            emb_table = data["embeddings"].astype(np.float64)
        vocab = {o: i for i, o in enumerate(emb_orfs)}
        strain_emb = np.zeros((len(orfs), emb_table.shape[1]))
        n_emb = 0
        for j, orf in enumerate(orfs):
            if orf in vocab:
                strain_emb[j] = emb_table[vocab[orf]]
                n_emb += 1
        # mean-centre: removes the shared magnitude direction that confounded
        # the first (uncentred) alignment version; matches the dual-arm gate
        strain_emb = strain_emb - strain_emb.mean(axis=0, keepdims=True)
        # confounder direction: mean embedding of the top-decile |z| strains
        magnitude_all = np.where(mask, np.abs(z), 0.0).mean(axis=0)
        cutoff = np.quantile(magnitude_all, 0.9)
        stress_strains = magnitude_all >= cutoff
        stress_dir = strain_emb[stress_strains].mean(axis=0)
        sd_norm = float(np.linalg.norm(stress_dir))
        if sd_norm > 0:
            stress_dir = stress_dir / sd_norm
        embed_meta = {
            "applied": True,
            "strains_with_embeddings": int(n_emb),
            "embedding_dim": int(emb_table.shape[1]),
            "stress_direction_strains": int(stress_strains.sum()),
            "roles": "A: embedding-space alignment; B: stress-projection confounder score",
        }
        print(f"scyeast embeddings: {embed_meta}", flush=True)

    split_by_key: dict[str, dict[str, str]] = {}
    with split_file.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            split_by_key[str(row["inchikey"]).upper()] = row
    candidate_records: list[dict[str, Any]] = []
    for index, identity in enumerate(inchikeys):
        observed = mask[index]
        count = int(observed.sum())
        if count < int(config["yeast_task"]["minimum_measured_strains"]):
            continue
        values = z[index, observed]
        weights_observed = translated[observed]
        value_norm = np.linalg.norm(values)
        weight_norm = np.linalg.norm(weights_observed)
        alignment = (
            float(np.dot(values, weights_observed) / (value_norm * weight_norm))
            if value_norm and weight_norm
            else 0.0
        )
        magnitude = float(np.mean(np.abs(values)))
        coverage = count / len(orfs)
        global_stress = magnitude
        score = alignment * min(1.0, coverage * 2.0) - 0.15 * min(1.0, global_stress / 3.0)

        embed_alignment = alignment
        stress_embed = min(1.0, global_stress / 3.0)
        if strain_emb is not None:
            e_obs = strain_emb[observed]
            intent_e = weights_observed @ e_obs
            resp_e = values @ e_obs
            ni = float(np.linalg.norm(intent_e))
            nr = float(np.linalg.norm(resp_e))
            embed_alignment = float(intent_e @ resp_e / (ni * nr)) if ni and nr else 0.0
            full_resp_e = np.where(mask[index], z[index], 0.0) @ strain_emb
            ne = float(np.linalg.norm(full_resp_e))
            stress_embed = float(abs(full_resp_e @ stress_dir) / ne) if ne and sd_norm > 0 else 0.0
        score_v2 = embed_alignment * min(1.0, coverage * 2.0) - 0.15 * stress_embed
        split = split_by_key.get(identity.upper(), {})
        candidate_records.append(
            {
                "candidate_id": f"YF-{index+1:04d}-{identity}",
                "compound_inchikey": identity,
                "compound_smiles": smiles[index],
                "split": split.get("split", "unresolved"),
                "scaffold_group_id": split.get("component_id"),
                "yeast_response_alignment": round(alignment, 8),
                "yeast_response_alignment_embedding": round(embed_alignment, 8),
                "yeast_coverage": round(coverage, 8),
                "global_stress_penalty": round(global_stress, 8),
                "global_stress_penalty_embedding": round(stress_embed, 8),
                "yeast_first_score": round(score, 8),
                "yeast_first_score_v2": round(score_v2, 8),
                "measured_strain_count": count,
                "evidence_mode": "PUBLIC_RETROSPECTIVE",
                "scientific_claim_eligible": False,
                "target_status": "BACKMAP_REQUIRED",
                "claim_boundary": (
                    "Yeast-first public chemical-genetic candidate; not a "
                    "target-dependent yeast hit."
                ),
            }
        )
    development = [x for x in candidate_records if x["split"] in {"train", "validation"}]
    # selection stays on the raw-z score. Both embedding-alignment versions
    # (raw and mean-centred) were diagnosed magnitude-confounded: centred-vs-
    # raw alignment correlation 0.114 and top ranks still dominated by
    # high-stress compounds with near-zero raw alignment. scYeast embedding
    # columns remain as representation diagnostics; the dual-arm retrieval
    # gate gives the final validated verdict on this representation.
    sort_key = "yeast_first_score"
    development.sort(key=lambda x: (-x[sort_key], x["compound_inchikey"]))
    selected = development[: int(config["selection"]["budget"])]
    backmap = _chembl_backmap(
        [str(item["compound_inchikey"]) for item in selected], config
    )
    family_rules = config.get("target_backmapping", {}).get("family_rules", {})
    allowed_families = set(family_rules) if isinstance(family_rules, Mapping) else set()
    for item in selected:
        relations = backmap.get(item["compound_inchikey"], [])
        membrane_relations = [
            relation
            for relation in relations
            if relation.get("target_family") in allowed_families
        ]
        item["human_target_backmapping"] = relations
        item["membrane_target_backmapping"] = membrane_relations
        item["target_status"] = (
            "MEMBRANE_TARGET_RELATIONS_FOUND"
            if membrane_relations
            else (
                "NO_MEMBRANE_TARGET_BACKMAP"
                if relations
                else "NO_CHEMBL_TARGET_BACKMAP"
            )
        )
    return {
        "schema_version": "yeastbridge.discovery.yeast-first-panel.v1",
        "mode": "yeast_first",
        "human_intent": {
            "config": config["human_intent"],
            "source_sha256": intent["source_sha256"],
            "dual_track_fusion": intent.get("dual_track_fusion", {}),
            "translated_feature_count": len(translated),
            **translation_meta,
        },
        "scyeast_dual_base": embed_meta,
        "yeast_sources": {
            "matrix_path": str(matrix.relative_to(base)),
            "matrix_sha256": sha256_file(matrix),
            "split_path": str(split_file.relative_to(base)),
            "split_sha256": sha256_file(split_file),
            "response_source": "ChemGRID_AID_1159580",
        },
        "records": selected,
        "development_pool_count": len(development),
        "all_eligible_pool_count": len(candidate_records),
        "target_backmapping": {
            "status": "completed_for_selected_compounds",
            "source": "ChEMBL activity relation API",
            "target_relation_count": sum(
                len(item["human_target_backmapping"]) for item in selected
            ),
            "membrane_target_relation_count": sum(
                len(item["membrane_target_backmapping"]) for item in selected
            ),
            "compound_with_target_relations": sum(
                bool(item["human_target_backmapping"]) for item in selected
            ),
            "compound_with_membrane_target_relations": sum(
                bool(item["membrane_target_backmapping"]) for item in selected
            ),
        },
        "scientific_claim_eligible": False,
        "panel_sha256": _hash(selected),
        "claim_boundary": (
            "Human intent conditions the yeast task; public yeast response ranks "
            "compounds before target backmapping."
        ),
    }


def _overwrite_json(path: Path, document: dict[str, Any]) -> None:
    """Replace an existing result file (addendum-4 overwrite policy)."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        target.unlink()
    exclusive_json_write(target, document)


def _revision_diff(
    previous: dict[str, Any], current: dict[str, Any]
) -> dict[str, Any]:
    prev_top = {
        r["candidate_id"]: float(r.get("yeast_first_score") or 0.0)
        for r in previous.get("records", [])
    }
    curr_top = {
        r["candidate_id"]: float(r.get("yeast_first_score") or 0.0)
        for r in current.get("records", [])
    }
    overlap = set(prev_top) & set(curr_top)
    union = set(prev_top) | set(curr_top)
    from scipy.stats import spearmanr

    shared = sorted(overlap)
    rho = None
    if len(shared) >= 3:
        rho = round(
            float(
                spearmanr(
                    [prev_top[k] for k in shared],
                    [curr_top[k] for k in shared],
                ).statistic
            ),
            6,
        )
    entered = sorted(set(curr_top) - set(prev_top))
    dropped = sorted(set(prev_top) - set(curr_top))
    return {
        "panel_size_previous": len(prev_top),
        "overlap_count": len(overlap),
        "jaccard": round(len(overlap) / len(union), 6) if union else None,
        "score_spearman_on_overlap": rho,
        "entered": [
            {"candidate_id": k, "yeast_first_score": curr_top[k]} for k in entered[:8]
        ],
        "dropped": [
            {"candidate_id": k, "yeast_first_score": prev_top[k]} for k in dropped[:8]
        ],
    }


def run_translation_validation(
    root: str,
    *,
    config_path: str,
    alpha_grid: tuple[float, ...] = (0.15, 0.30, 0.50),
) -> dict[str, Any]:
    """Feasibility analysis OF the intent-translation method itself.

    Measurements on real frozen inputs:

    1. anchor recovery — half of a synthetic anchor set is hidden, propagation
       runs on the visible half, and rank correlation against the hidden true
       weights is measured at held-out strains;
    2. concentration diagnostics of the real translated vector under both the
       retired neighbour-mean diffusion and personalised PageRank;
    3. reference block from the previously released panel, when present.
    """

    import random

    from scipy.stats import spearmanr

    base = project_root(root)
    config = _load_json(base, config_path)
    intent_section = config["human_intent"]
    intent = _load_intent(base, {"human_intent": intent_section})
    orthology = _load_orthology(base, intent_section["orthology_path"])
    neighbors = _load_sga_neighbors(base, intent_section["sga_path"])
    matrix = resolve_in_root(base, config["yeast_task"]["matrix_path"], must_exist=True)
    with np.load(matrix, allow_pickle=False) as data:
        strain_orfs = [str(item) for item in data["strain_orfs"].astype(str)]
        mask = (
            data["measured_mask"].astype(bool)
            & np.isfinite(data["z_score"].astype(np.float64))
        )

    def legacy_mean_propagate(anchor_weights: dict[str, float]) -> np.ndarray:
        propagated: dict[str, float] = dict(anchor_weights)
        for _ in range(2):
            nxt = dict(propagated)
            for orf, connected in neighbors.items():
                values = [propagated[i] for i in connected if i in propagated]
                if values:
                    nxt[orf] = float(np.mean(values))
            propagated = nxt
        return np.array([propagated.get(o, 0.0) for o in strain_orfs])

    all_nodes: set[str] = set(strain_orfs)
    for orf, connected in neighbors.items():
        all_nodes.add(orf)
        all_nodes.update(connected)
    nodes = sorted(all_nodes)

    def ppr_propagate(anchor_weights: dict[str, float], alpha: float) -> np.ndarray:
        propagated, _ = _propagate_personalized_pagerank(
            nodes, neighbors, anchor_weights, alpha=alpha
        )
        node_index = {node: i for i, node in enumerate(nodes)}
        return np.array(
            [
                propagated[node_index[o]] if o in node_index else 0.0
                for o in strain_orfs
            ]
        )

    methods: dict[str, Any] = {"legacy_mean_two_hops": legacy_mean_propagate}
    for alpha in alpha_grid:
        methods[f"pagerank_alpha_{alpha}"] = (
            lambda a, alpha=alpha: ppr_propagate(a, alpha)
        )

    rng = random.Random(20260824)
    eligible = [i for i, orf in enumerate(strain_orfs) if mask[:, i].mean() >= 0.5]
    anchor_k = min(30, len(eligible))
    trials: list[dict[str, Any]] = []
    for trial in range(8):
        chosen = rng.sample(eligible, anchor_k)
        truth = {
            strain_orfs[i]: float(rng.choice([-1.0, 1.0]) * rng.uniform(0.5, 3.0))
            for i in chosen
        }
        holdout = sorted(truth)[::2]
        visible = {k: v for k, v in truth.items() if k not in set(holdout)}
        hidden_true = [truth[k] for k in holdout]
        row: dict[str, Any] = {"trial": trial, "n_anchors_visible": len(visible)}
        for name, fn in methods.items():
            vector = fn(visible)
            predicted = [vector[strain_orfs.index(k)] for k in holdout]
            result = spearmanr(predicted, hidden_true)
            value = float(result.statistic) if result.statistic is not None else float("nan")
            row[name] = round(value, 6) if np.isfinite(value) else None
        trials.append(row)

    summary_by_method: dict[str, Any] = {}
    for name in methods:
        values = [t[name] for t in trials if t.get(name) is not None]
        summary_by_method[name] = {
            "mean_spearman_holdout_recovery": round(sum(values) / len(values), 6)
            if values
            else None,
            "min": round(min(values), 6) if values else None,
            "max": round(max(values), 6) if values else None,
            "n_trials": len(values),
        }

    # --- fidelity test: smooth ground truth (the recoverable regime) ---
    # Independent-random labels are unpredictable by construction (their
    # ceiling is zero), so that arm is retained only as a capacity null
    # control. Here the ground truth is generated BY propagation of hidden
    # seeds at a moderate alpha, which is exactly the signal class the
    # operator is supposed to carry; methods must reproduce it from the
    # seeds alone at non-seed strains.
    alpha_true = 0.25
    smooth_trials: list[dict[str, Any]] = []
    for trial in range(8):
        chosen = rng.sample(eligible, anchor_k)
        seeds = {
            strain_orfs[i]: float(rng.choice([-1.0, 1.0]) * rng.uniform(0.5, 3.0))
            for i in chosen
        }
        truth = ppr_propagate(seeds, alpha_true)
        seed_set = set(seeds)
        evaluation = [
            i for i in eligible if strain_orfs[i] not in seed_set
        ]
        truth_values = [truth[i] for i in evaluation]
        row_smooth: dict[str, Any] = {
            "trial": trial,
            "n_seeds": len(seeds),
            "n_evaluation_strains": len(evaluation),
        }
        for name, fn in methods.items():
            vector = fn(seeds)
            predicted = [vector[i] for i in evaluation]
            result = spearmanr(predicted, truth_values)
            value = (
                float(result.statistic) if result.statistic is not None else float("nan")
            )
            row_smooth[name] = round(value, 6) if np.isfinite(value) else None
        smooth_trials.append(row_smooth)

    smooth_summary: dict[str, Any] = {}
    for name in methods:
        values = [t[name] for t in smooth_trials if t.get(name) is not None]
        smooth_summary[name] = {
            "mean_spearman_fidelity": round(sum(values) / len(values), 6)
            if values
            else None,
            "min": round(min(values), 6) if values else None,
            "max": round(max(values), 6) if values else None,
            "n_trials": len(values),
        }

    translated_new, meta_new = _feature_weights(
        strain_orfs, orthology, intent["weights"], neighbors
    )

    def concentration(vector: np.ndarray) -> dict[str, Any]:
        absolute = np.sort(np.abs(vector))
        total = float(absolute.sum())
        n = absolute.size
        if n == 0 or total == 0.0:
            return {"gini_abs": 0.0}
        cumulative = np.cumsum(absolute) / total
        gini = float(
            1.0 - 2.0 * float((cumulative - np.arange(1, n + 1) / n).sum()) / n + 1.0 / n
        )
        mad = float(np.median(np.abs(vector - np.median(vector))))
        peak = float(np.abs(vector - np.median(vector)).max())
        return {
            "gini_abs": round(gini, 6),
            "max_abs_deviation_in_mad_units": round(peak / mad, 3) if mad else None,
        }

    legacy_vector = legacy_mean_propagate(
        {
            orf: -float(
                np.mean([intent["weights"][g] for g in genes if g in intent["weights"]])
            )
            for orf, genes in orthology.items()
            if any(g in intent["weights"] for g in genes)
        }
    )

    previous_block: dict[str, Any] = {"available": False}
    panel_reference = (
        config.get("panel", {}).get("path")
        or config.get("existing_panel_path")
    )
    panel_path = resolve_in_root(base, panel_reference) if panel_reference else None
    if panel_path is not None and panel_path.exists():
        previous_panel = json.loads(panel_path.read_text(encoding="utf-8"))
        hi = previous_panel.get("human_intent", {})
        previous_block = {
            "available": True,
            "panel_sha256_previous_run": sha256_file(panel_path),
            "previous_human_intent_summary": {
                key: hi.get(key)
                for key in (
                    "direct_anchor_count",
                    "chemgrid_translated_count",
                    "translation_method",
                    "translation_rule",
                )
            },
        }

    return {
        "schema_version": "yeastbridge.discovery.translation-validation.v1",
        "question": (
            "Two separated questions: (a) capacity null - can ANY operator "
            "predict independently-random held-out labels through this graph? "
            "(ceiling is zero by construction); (b) transmission fidelity - "
            "when the ground truth is itself a smooth propagated signal, how "
            "much of it does each method reproduce from the seeds alone?"
        ),
        "anchor_recovery_trials": trials,
        "anchor_recovery_summary": {
            **summary_by_method,
            "_interpretation": (
                "independent-random labels are unpredictable by construction; "
                "values near zero here are the expected ceiling, NOT evidence "
                "against the graph"
            ),
        },
        "smooth_truth_recovery_trials": smooth_trials,
        "smooth_truth_recovery_summary": smooth_summary,
        "real_intent_concentration": {
            "legacy_mean_two_hops": concentration(legacy_vector),
            "personalized_pagerank_alpha_0.30": concentration(translated_new),
        },
        "revised_translation_meta": meta_new,
        "previous_panel": previous_block,
        "claim_boundary": (
            "Calibration diagnostics of the translation operator on frozen "
            "public inputs; they do not validate any compound ranking. The "
            "fidelity arm is the decision-relevant one."
        ),
    }


def compound_target_backmap(
    root: str | Path,
    *,
    config_path: str | Path,
    z_threshold: float = 2.0,
    top_genes: int = 25,
) -> dict[str, Any]:
    """Reverse translation: compound-sensitive yeast strains -> signed
    personalised PageRank over the same SGA graph used forward -> OrthoDB
    human gene aggregation. Target hypotheses are derived from yeast data
    only; the CRC scan rank is attached as an annotation column, not a gate.
    """
    base = project_root(root)
    config = _load_json(base, config_path)
    orthology = _load_orthology(base, config["human_intent"]["orthology_path"])
    neighbors = _load_sga_neighbors(base, config["human_intent"]["sga_path"])
    matrix = resolve_in_root(base, config["yeast_task"]["matrix_path"], must_exist=True)
    with np.load(matrix, allow_pickle=False) as data:
        inchikeys = data["compound_inchikeys"].astype(str)
        smiles = data["compound_smiles"].astype(str)
        z = data["z_score"].astype(np.float64)
        mask = data["measured_mask"].astype(bool) & np.isfinite(z)
        orfs = data["strain_orfs"].astype(str).tolist()

    panel_candidates = [
        base / "artifacts/discovery_release_v1/yeast_first_20260826/panel_v8.json",
        base / "artifacts/discovery_release_v1/yeast_first_20260823/panel_v7.json",
    ]
    panel_path = next((p for p in panel_candidates if p.exists()), panel_candidates[-1])
    selected: dict[str, float | None] = {}
    if panel_path.exists():
        panel = json.loads(panel_path.read_text(encoding="utf-8"))
        selected = {
            r["compound_inchikey"]: r.get("yeast_first_score") for r in panel["records"]
        }

    scan_path = base / "reports/crc_target_scan/model_evidence_v2/dual_evidence_ranking.tsv"
    scan_rank: dict[str, int] = {}
    if scan_path.exists():
        with scan_path.open(encoding="utf-8-sig", newline="") as handle:
            for i, row in enumerate(csv.DictReader(handle, delimiter="\t"), start=1):
                scan_rank[str(row.get("target_id", "")).upper()] = i

    nodes = sorted(neighbors.keys())
    rows: list[dict[str, Any]] = []
    for ci, key in enumerate(inchikeys):
        if selected and key not in selected:
            continue
        zc, mc = z[ci], mask[ci]
        sensitive = mc & (np.abs(zc) >= z_threshold)
        n_sens = int(sensitive.sum())
        record = {
            "compound_inchikey": key,
            "compound_smiles": smiles[ci],
            "yeast_first_score": selected.get(key),
            "sensitive_strain_count": n_sens,
            "z_threshold": z_threshold,
        }
        anchors = {orfs[j]: float(zc[j]) for j in np.where(sensitive)[0] if orfs[j] in neighbors}
        if n_sens < 3 or not anchors:
            record["status"] = "insufficient_sensitive_strains"
            rows.append(record)
            continue
        propagated, iterations = _propagate_personalized_pagerank(nodes, neighbors, anchors)
        gene_scores: dict[str, float] = {}
        for node, value in zip(nodes, propagated):
            if abs(float(value)) < 1e-9:
                continue
            for human in orthology.get(node, ()):
                gene_scores[human] = gene_scores.get(human, 0.0) + float(value)
        top = sorted(gene_scores.items(), key=lambda kv: -abs(kv[1]))[:top_genes]
        record.update(
            {
                "status": "ok",
                "propagation_iterations": iterations,
                "anchor_count": len(anchors),
                "human_target_hypotheses": [
                    {
                        "gene": gene,
                        "score": round(score, 6),
                        "crc_scan_rank": scan_rank.get(gene.upper()),
                    }
                    for gene, score in top
                ],
            }
        )
        rows.append(record)
        print(f"backmap {key}: {len(anchors)} anchors -> {len(gene_scores)} human genes", flush=True)

    ok = [r for r in rows if r.get("status") == "ok"]
    gene_vote: dict[str, int] = {}
    for r in ok:
        for hyp in r["human_target_hypotheses"][:10]:
            gene_vote[hyp["gene"]] = gene_vote.get(hyp["gene"], 0) + 1
    return {
        "schema_version": "1.0",
        "method": "compound-sensitive strains -> signed personalised PageRank over SGA -> orthology human aggregation",
        "z_threshold": z_threshold,
        "n_compounds_selected": len(selected) if selected else len(inchikeys),
        "n_compounds_backmapped": len(ok),
        "records": rows,
        "consensus_human_targets": sorted(
            ({"gene": g, "compound_votes": v} for g, v in gene_vote.items()),
            key=lambda d: -d["compound_votes"],
        )[:20],
        "claim_boundary": (
            "Yeast-data-derived target hypotheses; no binding, affinity or "
            "validated-hit claim. ChEMBL backmapping is verification annotation."
        ),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m yeastbridge_vs.discovery_release.yeast_first"
    )
    parser.add_argument(
        "command", nargs="?", default="build", choices=["build", "validate", "backmap"]
    )
    parser.add_argument("--root", default=".")
    parser.add_argument(
        "--config", default="configs/discovery_release/second_round_yeastfirst_v1.json"
    )
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    base = project_root(args.root)
    output = resolve_in_root(base, args.output)
    if args.command == "validate":
        document = run_translation_validation(args.root, config_path=args.config)
    elif args.command == "backmap":
        document = compound_target_backmap(args.root, config_path=args.config)
    else:
        document = build_yeast_first(args.root, config_path=args.config)
        if output.exists():
            try:
                previous = json.loads(output.read_text(encoding="utf-8"))
                document["revision_diff_vs_previous"] = _revision_diff(previous, document)
            except (json.JSONDecodeError, OSError):
                pass
    _overwrite_json(output, document)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
