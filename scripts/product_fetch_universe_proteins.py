#!/usr/bin/env python3
"""Fetch UniProt reviewed proteins for the 1177-gene target universe
(primary gene symbols from the frozen snapshots) for the route-B reverse
mapping. Same pattern as the candidate fetcher; paths from
configs/product_transfer.json."""
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CFG = json.loads((ROOT / "configs/product_transfer.json").read_text())
INPUTS = ROOT / CFG["inputs_dir"]

genes = [g.strip() for g in open(INPUTS / "universe_genes.txt") if g.strip()]
got = {}
for i in range(0, len(genes), 25):
    q = " OR ".join("gene_exact:" + g for g in genes[i:i + 25])
    url = ("https://rest.uniprot.org/uniprotkb/search?query=(" + urllib.parse.quote(q)
           + ")+AND+organism_id:9606+AND+reviewed:true&format=fasta&size=500")
    url = url.replace("%28", "(").replace("%29", ")").replace("%20", "+").replace("%3A", ":")
    fasta = None
    for attempt in range(4):  # UniProt 偶发 SSL 断连,退避重试
        try:
            fasta = urllib.request.urlopen(url, timeout=180).read().decode()
            break
        except Exception as e:
            print(f"batch {i//25}: attempt {attempt+1} failed ({e}); backoff {5*(attempt+1)}s", flush=True)
            time.sleep(5 * (attempt + 1))
    if fasta is None:
        raise RuntimeError(f"batch starting at gene {genes[i]} failed after retries")
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

n, missing = 0, []
out_seq = []
for g in genes:
    hit = next((k for k in got if " GN=" + g + " " in k + " "), None)
    if hit:
        acc = hit.split("|")[1]
        # sp|ACC|GENE 三段式头,提取器 parse_header 直接可解析(accession+GN=)
        out_seq.append((g, f">sp|{acc}|{g}_UMAN GN={g}\n{got[hit]}"))
        n += 1
    else:
        missing.append(g)
with open(INPUTS / "universe_proteins.fasta", "w") as f:
    for g, s in out_seq:
        f.write(s + "\n")
print(f"universe={len(genes)} fetched={n} missing={len(missing)} {missing[:8]}")
