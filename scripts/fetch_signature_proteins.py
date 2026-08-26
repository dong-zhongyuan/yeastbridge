#!/usr/bin/env python3
"""Fetch UniProt reviewed proteins for the full 256-gene disease signature
(extends the 50 anchor-gene set already downloaded) into
feasibility/transfer_routes/assets/signature_proteins.fasta."""
import csv
import urllib.parse
import urllib.request

from pathlib import Path

A = Path("/public/home/mengxl/dzy/yeastbridge_re/feasibility/transfer_routes/assets")

sig = [r["gene"] for r in csv.DictReader(open(A.parent.parent / "transfer/assets/state_signature.tsv"), delimiter="\t")]
done = set(l.strip() for l in open(A / "anchor_human_genes.txt") if l.strip())
todo = [g for g in sig if g not in done]
print(f"signature={len(sig)} already_downloaded={len(done)} need={len(todo)}")

got = {}
for i in range(0, len(todo), 25):
    q = " OR ".join("gene_exact:" + g for g in todo[i:i + 25])
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

have_syms = set()
for k in got:
    if " GN=" in k:
        have_syms.add(k.split(" GN=")[1].split()[0])
missing = [g for g in todo if g not in have_syms]

with open(A / "signature_proteins.fasta", "w") as f:
    f.write(open(A / "human_anchor_proteins.fasta").read())
    if got:
        f.write("\n".join(k + "\n" + got[k] for k in got) + "\n")
with open(A / "signature_genes.txt", "w") as f:
    f.write("\n".join(sig) + "\n")
print(f"downloaded_entries={len(got)} missing_genes={len(missing)} {missing[:12]}")
