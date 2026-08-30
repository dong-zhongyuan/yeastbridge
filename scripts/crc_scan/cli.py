"""Command-line entry point for the CRC target-agnostic baseline scan."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .core import run_scan


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--project-root", type=Path, required=True)
    result.add_argument("--config", type=Path, required=True)
    result.add_argument("--becker-h5ad", type=Path, required=True)
    result.add_argument("--xatlas-dir", type=Path, required=True)
    result.add_argument("--xatlas-gene-metadata", type=Path, required=True)
    result.add_argument("--norman-h5ad", type=Path, required=True)
    result.add_argument("--ortholog-map", type=Path, required=True)
    result.add_argument("--output-dir", type=Path, required=True)
    result.add_argument("--hash-workers", type=int, default=4)
    result.add_argument(
        "--skip-input-hashes",
        action="store_true",
        help="Development-only; production reports remain blocked without content hashes.",
    )
    return result


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    report = run_scan(
        project_root=args.project_root,
        config_path=args.config,
        becker_h5ad=args.becker_h5ad,
        xatlas_dir=args.xatlas_dir,
        xatlas_gene_metadata=args.xatlas_gene_metadata,
        norman_h5ad=args.norman_h5ad,
        ortholog_map=args.ortholog_map,
        output_dir=args.output_dir,
        hash_workers=args.hash_workers,
        hash_inputs=not args.skip_input_hashes,
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if report["status"] != "failed" else 2
