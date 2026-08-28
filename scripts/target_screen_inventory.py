#!/usr/bin/env python
"""target_screen stage 1: structure inventory for the target universe.

Per target (UniProt accession from the universe fasta): query RCSB for
deposited experimental structures (entity-level UniProt exact match), pick
the highest-resolution entry; fall back to the AlphaFold DB v4 predicted
model when nothing is deposited. mmCIF inputs are converted to PDB with
gemmi (AF2 B-factor = pLDDT preserved). Resumable per target via
results/inventory.tsv; network-failed targets are not written and are
retried on rerun. Registered in product/target_screen/DESIGN.md; all
parameters in configs/target_screen.json.
"""
from __future__ import annotations

import argparse
import json
import time
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

SEARCH_ATTR = (
    "rcsb_polymer_entity_container_identifiers."
    "reference_sequence_identifiers.database_accession"
)


def read_pairs(fasta: Path):
    """(acc, gene) pairs; gene from GN= token, entry-name prefix fallback."""
    out = []
    for line in fasta.read_text().splitlines():
        if line.startswith(">"):
            parts = line[1:].split("|")
            if len(parts) >= 3:
                acc = parts[1]
                toks = parts[2].split()
                gene = next((t[3:] for t in toks if t.startswith("GN=")),
                            toks[0].split("_")[0])
                out.append((acc, gene))
            else:
                first = parts[0].split()[0]
                out.append((first, first))
    return out


def http(url, payload=None, tries=6, timeout=180):
    """GET bytes or POST JSON; returns (status, bytes); (None, b"") on network failure."""
    for t in range(tries):
        try:
            data = json.dumps(payload).encode() if payload is not None else None
            req = urllib.request.Request(
                url, data=data,
                headers={"Content-Type": "application/json"} if data else {},
            )
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.status, r.read()
        except urllib.error.HTTPError as e:
            if e.code in (404, 204):
                return e.code, b""
            err = f"HTTP {e.code}"
        except Exception as e:  # noqa: BLE001
            err = str(e)
        print(f"  retry {t + 1}/{tries} {url[:70]}: {err}", flush=True)
        time.sleep(min(300, 5 * 2 ** t))
    return None, b""


def to_pdb(cif_bytes: Path, pdb_out: Path):
    import gemmi

    st = gemmi.read_structure(str(cif_bytes))
    if len(st) > 1:
        del st[1:]
    # PDB format allows only single-character chain IDs; some entries carry
    # multi-character auth chain names (e.g. AAA)
    used = set()
    for model in st:
        for ch in model:
            if len(ch.name) > 1 or ch.name in used or not ch.name:
                new = next((c for c in
                            "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz"
                            if c not in used), "X")
                used.add(new)
                ch.name = new
            else:
                used.add(ch.name)
    st.setup_entities()
    st.write_pdb(str(pdb_out))


def detect_state(pdb_id, title, method):
    """Heuristic conformational state detection from PDB metadata."""
    t = (title or "").lower()
    if any(k in t for k in ["goa", "gi1", "gi2", "gs", "g protein",
                             "gprotein", "mini-g", "nanobody", "active"]):
        return "active"
    if any(k in t for k in ["antagonist", "inverse agonist", "inactive"]):
        return "inactive"
    if "agonist" in t and "inverse" not in t:
        return "active"
    if "blocker" in t or "inhibitor" in t:
        return "inactive"
    return "unannotated"


def pick_best(entries_with_meta):
    """Pick best entry; prefer experimental PDB over predicted AF2."""
    if not entries_with_meta:
        return None
    return min(entries_with_meta, key=lambda x: x[1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/target_screen.json")
    args = ap.parse_args()
    cfg = json.loads((ROOT / args.config).read_text())

    raw_dir = ROOT / cfg["structures_dir"] / "raw"
    res_dir = ROOT / cfg["results_dir"]
    raw_dir.mkdir(parents=True, exist_ok=True)
    res_dir.mkdir(parents=True, exist_ok=True)
    inv_path = res_dir / "inventory.tsv"

    universe = read_pairs(ROOT / cfg["universe_fasta"])
    done = set()
    if inv_path.exists():
        for line in inv_path.read_text().splitlines()[1:]:
            done.add(line.split("\t")[0])
    print(f"targets: {len(universe)}, already inventoried: {len(done)}",
          flush=True)

    fresh = not inv_path.exists()
    n_pdb = n_af2 = n_fail = 0
    with inv_path.open("a", newline="") as fh:
        if fresh:
            fh.write("acc\tgene\tsource\tpdb_id\tresolution\tmethod\tpath\tnote\tconformational_state\n")
        for acc, gene in universe:
            if acc in done:
                continue
            row = None
            q = {
                "query": {
                    "type": "terminal", "service": "text",
                    "parameters": {
                        "attribute": SEARCH_ATTR,
                        "operator": "exact_match", "value": acc,
                    },
                },
                "return_type": "entry",
                "request_options": {
                    "paginate": {"start": 0, "rows": 25},
                    "results_content_type": ["experimental"],
                },
            }
            st, body = http(cfg["rcsb_search"], payload=q)
            if st is None:
                n_fail += 1
                continue  # network failure: no write, retried on rerun
            entries = []
            if st == 200 and body:
                entries = [h["identifier"]
                           for h in json.loads(body).get("result_set", [])]
            if entries:
                # collect ALL entries with metadata + state
                entry_meta = []
                for pid in sorted(entries):
                    st2, b2 = http(f"{cfg['rcsb_entry']}/{pid}")
                    if st2 != 200:
                        continue
                    ent = json.loads(b2)
                    meth = (ent.get("exptl") or [{}])[0].get("method", "")
                    title = ent.get("struct", {}).get("title", "")
                    res_list = (ent.get("rcsb_entry_info", {})
                                .get("resolution_combined") or [])
                    resv = res_list[0] if res_list else 999.0
                    state = detect_state(pid, title, meth)
                    entry_meta.append((pid, resv, meth, state))
                for pid, resv, meth, state in entry_meta:
                    pdb_path = raw_dir / f"{acc}_{pid}.pdb"
                    if pdb_path.exists():
                        pass
                    else:
                        st3, b3 = http(
                            cfg["rcsb_download"].format(pdb_id=pid))
                        if st3 == 200:
                            pdb_path.write_bytes(b3)
                        else:
                            st4, b4 = http(
                                cfg["rcsb_download_cif"].format(
                                    pdb_id=pid))
                            if st4 != 200:
                                n_fail += 1
                                continue
                            cif = raw_dir / f"{acc}_{pid}.cif"
                            cif.write_bytes(b4)
                            to_pdb(cif, pdb_path)
                            cif.unlink()
                    fh.write("\t".join([acc, gene, "pdb", pid,
                                        "" if resv == 999.0 else f"{resv:.2f}",
                                        meth,
                                        str(pdb_path.relative_to(ROOT)),
                                        "", state]) + "\n")
                    n_pdb += 1
                fh.flush()

            # Always also get AF2 model (inactive-like fallback)
            af2_path = raw_dir / f"{acc}_af2.pdb"
            if not af2_path.exists():
                st5, b5 = http(cfg["af2_api"].format(acc=acc))
                url = None
                if st5 == 200 and b5:
                    try:
                        pred = json.loads(b5)[0]
                        url = pred.get("pdbUrl") or pred.get("cifUrl")
                    except Exception:  # noqa: BLE001
                        url = None
                if url:
                    st6, b6 = http(url)
                    if st6 == 200 and len(b6) > 1000:
                        if url.endswith(".pdb"):
                            af2_path.write_bytes(b6)
                        else:
                            cif = raw_dir / f"{acc}_af2.cif"
                            cif.write_bytes(b6)
                            to_pdb(cif, af2_path)
                            cif.unlink()
            if af2_path.exists():
                fh.write("\t".join([acc, gene, "af2", "", "", "predicted",
                                    str(af2_path.relative_to(ROOT)),
                                    "alphafold_db",
                                    "af2_inactive_like"]) + "\n")
                n_af2 += 1
            fh.flush()
            if (n_pdb + n_af2) % 25 == 0:
                print(f"  {n_pdb + n_af2} done (pdb={n_pdb} af2={n_af2} "
                      f"fail={n_fail})", flush=True)
    print(f"FINISHED new: pdb={n_pdb} af2={n_af2} network_fail={n_fail}",
          flush=True)


if __name__ == "__main__":
    main()
