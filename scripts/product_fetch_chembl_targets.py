#!/usr/bin/env python3
"""Fetch ChEMBL known targets for the significant HIP/HOP compounds
(target annotation; design: product/drug_annotation/DESIGN.md).
inchikey -> ChEMBL molecule -> activities -> target gene symbols.
Robust to EBI throttling (exponential backoff, resumable, final retry).
Gene symbols resolved in a final backfill over the distinct target set:
ChEMBL deprecated target_components[].target_component_symbol, so symbols
are read from target_component_synonyms[syn_type=GENE_SYMBOL] with the
HGNC xref as fallback."""
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CFG = json.loads((ROOT / "configs/drug_annotation.json").read_text())
OUT = ROOT / CFG["results_dir"]
BASE = "https://www.ebi.ac.uk/chembl/api/data"


def get(url, tries=8):
    for a in range(tries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            return json.loads(urllib.request.urlopen(req, timeout=60).read().decode())
        except Exception as e:
            wait = min(600, 10 * (2 ** a))
            print(f"backoff {wait}s ({e})", flush=True)
            time.sleep(wait)
    return None


sig = pd.read_csv(ROOT / CFG["exec_matrix"], sep="\t")
iks = sorted(set(sig.loc[sig["q"] < CFG["fdr_alpha"], "inchikey"]))
print(f"compounds to fetch: {len(iks)}", flush=True)

OUT.mkdir(parents=True, exist_ok=True)
done_path = OUT / "chembl_targets.tsv"
HEADER = "inchikey\tmolecule_chembl_id\ttarget_chembl_id\ttarget_gene\tpchembl\n"
# 断点续传:已有任意行的化合物视为完成;单行空记录=查过且 ChEMBL 无此分子
done = set()
if done_path.exists():
    for r in pd.read_csv(done_path, sep="\t").fillna("").to_dict("records"):
        done.add(r["inchikey"])
print(f"already fetched: {len(done)}", flush=True)

fresh = not done_path.exists()
f = open(done_path, "a", newline="")
if fresh:
    f.write(HEADER)
import csv
w = csv.writer(f, delimiter="\t")
failed = []
for i, ik in enumerate(iks):
    if ik in done:
        continue
    m = get(f"{BASE}/molecule.json?molecule_structures__standard_inchi_key={ik}&limit=1")
    if m is None:
        failed.append(ik)  # 网络失败不落盘,末轮重试
        continue
    if not m.get("molecules"):
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
    for r in rows:
        w.writerow([r[0], r[1], r[2], "", r[4]])
    f.flush()
    if (i + 1) % 10 == 0:
        print(f"{i+1}/{len(iks)} compounds done", flush=True)
# 末轮:网络失败的化合物整体重试一次
if failed:
    print(f"final retry pass: {len(failed)} network-failed compounds", flush=True)
    done_now = {r for r in pd.read_csv(done_path, sep="\t")["inchikey"]}
    for ik in failed:
        if ik in done_now:
            continue
        m = get(f"{BASE}/molecule.json?molecule_structures__standard_inchi_key={ik}&limit=1")
        if m is None or not m.get("molecules"):
            w.writerow([ik, "", "", "", ""])
            f.flush()
f.close()
# 末段:靶点基因符号统一回填(对去重后的靶点集合各查一次,而非逐化合物)
allrows = [ln.split("\t") for ln in done_path.read_text().splitlines()[1:]]
tcs = sorted({r[2] for r in allrows if len(r) > 2 and r[2]})
print(f"distinct targets to resolve: {len(tcs)}", flush=True)
gene_of = {}
for j, tc in enumerate(tcs, 1):
    t = get(f"{BASE}/target/{tc}.json")
    if t:
        syms = []
        for c in t.get("target_components", []) or []:
            for s in c.get("target_component_synonyms", []) or []:
                if s.get("syn_type") == "GENE_SYMBOL":
                    v = s.get("component_synonym", "")
                    if v and v not in syms:
                        syms.append(v)
        if not syms:
            for c in t.get("target_components", []) or []:
                for x in c.get("target_component_xrefs", []) or []:
                    if x.get("xref_src_db") == "HGNC" and x.get("xref_name"):
                        if x["xref_name"] not in syms:
                            syms.append(x["xref_name"])
        if syms:
            gene_of[tc] = " ".join(syms)
    time.sleep(0.2)
    if j % 50 == 0:
        print(f"  resolved {j}/{len(tcs)} targets", flush=True)
tmp = done_path.with_suffix(".tsv.tmp")
with open(tmp, "w", newline="") as fh:
    fh.write(HEADER)
    for r in allrows:
        if len(r) >= 5 and r[2]:
            r[3] = gene_of.get(r[2], "")
        fh.write("\t".join(r[:5]) + "\n")
tmp.replace(done_path)
print(f"gene backfill complete: {sum(1 for r in allrows if len(r) >= 5 and r[3])} rows annotated", flush=True)
print("[done]", flush=True)
