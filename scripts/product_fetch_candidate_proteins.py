#!/usr/bin/env python3
"""Fetch UniProt reviewed proteins for the registered product candidates
(product step 1 output) into the transfer step's inputs directory.

Gene -> Entry accession comes from the frozen universe snapshots; FASTAs are
fetched by accession in batches. All paths from configs/product_transfer.json.
Pattern copied from the earlier signature fetcher."""
import csv
import json
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CFG = json.loads((ROOT / "configs/product_transfer.json").read_text())
YB = Path(CFG["yeastbridge_root"])

RANKING = ROOT / CFG["candidate_ranking"]
SNAP = YB / CFG["universe_snapshots"]
OUT = ROOT / CFG["inputs_dir"]

cands = list(csv.DictReader(open(RANKING), delimiter="\t"))
gene2acc = {}
# 仅 GPCR 快照带 Entry 列;离子通道快照是另一种导出(无 accession)。
# 统一改用 gene_exact 查询(与签名下载器同法),GN= 标签回对基因名。
gene_set = [c["target_id"] for c in cands]
got = {}
for i in range(0, len(gene_set), 25):
    q = " OR ".join("gene_exact:" + g for g in gene_set[i:i + 25])
    url = ("https://rest.uniprot.org/uniprotkb/search?query=(" + urllib.parse.quote(q)
           + ")+AND+organism_id:9606+AND+reviewed:true&format=fasta&size=500")
    url = url.replace("%28", "(").replace("%29", ")").replace("%20", "+").replace("%3A", ":")
    fasta = urllib.request.urlopen(url, timeout=180).read().decode()
    name, chunks = None, []
    for line in fasta.splitlines():
        if line.startswith(">"):
            if name:
                got.setdefault(name, "".join(chunks))
            name, chunks = line, []
        else:
            chunks.append(line)
    if name:
        got.setdefault(name, "".join(chunks))

accs, missing = [], []
for c in cands:
    hit = None
    for k in got:
        if " GN=" + c["target_id"] + " " in k + " ":
            hit = k
            break
    if hit:
        accs.append((c["target_id"], hit.split("|")[1]))
    else:
        missing.append(c["target_id"])
print(f"candidates={len(cands)} mapped={len(accs)} missing={missing}")

acc2seq = {}
for k, v in got.items():
    parts = k.split("|")
    if len(parts) >= 3:
        acc2seq[parts[1]] = v

OUT.mkdir(parents=True, exist_ok=True)
n = 0
with open(OUT / "candidate_proteins.fasta", "w") as f:
    for g, a in accs:
        if a in acc2seq:
            f.write(f">{g}|{a}\n{acc2seq[a]}\n")
            n += 1
with open(OUT / "candidate_accessions.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["target_id", "uniprot_acc", "fetched"])
    for g, a in accs:
        w.writerow([g, a, int(a in acc2seq)])
print(f"wrote {n} sequences -> {OUT}/candidate_proteins.fasta")
