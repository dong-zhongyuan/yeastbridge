#!/usr/bin/env python3
"""Fetch ChEMBL known targets for the significant HIP/HOP compounds
(knowledge-channel fix; design: product/backtrace_route_b/RESULTS.md).
inchikey -> ChEMBL molecule -> activities -> target gene symbols.
Robust to EBI throttling (exponential backoff, resumable)."""
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CFG = json.loads((ROOT / "configs/product_backtrace.json").read_text())
OUT = ROOT / CFG["results_dir"]
BASE = "https://www.ebi.ac.uk/chembl/api/data"


def get(url, tries=6):
    for a in range(tries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            return json.loads(urllib.request.urlopen(req, timeout=60).read().decode())
        except Exception as e:
            wait = min(300, 10 * (2 ** a))
            print(f"backoff {wait}s ({e})", flush=True)
            time.sleep(wait)
    return None


sig = pd.read_csv(ROOT / CFG["exec_matrix"], sep="\t")
iks = sorted(set(sig.loc[sig["q"] < CFG["fdr_alpha"], "inchikey"]))
print(f"compounds to fetch: {len(iks)}", flush=True)

OUT.mkdir(parents=True, exist_ok=True)
done_path = OUT / "chembl_targets.tsv"
done = set()
if done_path.exists():
    for r in pd.read_csv(done_path, sep="\t").fillna("").to_dict("records"):
        if r["gene_symbols"]:
            done.add(r["inchikey"])
print(f"already fetched: {len(done)}", flush=True)

f = open(done_path, "a", newline="")
f.write("inchikey\tmolecule_chembl_id\ttarget_chembl_id\ttarget_gene\tpchembl\n") if not done_path.exists() else None
import csv
w = csv.writer(f, delimiter="\t")
for i, ik in enumerate(iks):
    if ik in done:
        continue
    m = get(f"{BASE}/molecule.json?molecule_structures__standard_inchi_key={ik}&limit=1")
    if not m or not m.get("molecules"):
        w.writerow([ik, "", "", "", ""])
        f.flush()
        continue
    mid = m["molecules"][0]["molecule_chembl_id"]
    off = 0
    rows = []
    while True:
        a = get(f"{BASE}/activity.json?molecule_chembl_id={mid}&limit=100&offset={off}")
        if not a or not a.get("activities"):
            break
        for act in a["activities"]:
            tc = act.get("target_chembl_id") or ""
            if not tc:
                continue
            rows.append([ik, mid, tc, "", act.get("pchembl_value") or ""])
        total = a.get("page_meta", {}).get("total_count", 0)
        off += 100
        if off >= total or off > 2000:
            break
    # 批量解析 target 基因符号
    tcs = sorted({r[2] for r in rows if r[2]})
    gene_of = {}
    for tc in tcs:
        t = get(f"{BASE}/target/{tc}.json")
        if t and t.get("target_components"):
            genes = [c.get("target_component_symbol", "") for c in t["target_components"]]
            gene_of[tc] = " ".join(g for g in genes if g)
        time.sleep(0.3)
    for r in rows:
        w.writerow([r[0], r[1], r[2], gene_of.get(r[2], ""), r[4]])
    f.flush()
    if (i + 1) % 10 == 0:
        print(f"{i+1}/{len(iks)} compounds done", flush=True)
print("[done]", flush=True)
