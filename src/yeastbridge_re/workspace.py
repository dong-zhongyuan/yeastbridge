"""Workspace boundary, provenance snapshots, hashing, and exclusive writes."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any


class WorkspaceError(ValueError):
    """Raised when a path or workspace operation is unsafe."""


def project_root(path: str | Path | None = None) -> Path:
    """Resolve a Git project root without relying on newer ``git -C`` support."""

    start = Path(path or Path.cwd()).expanduser().resolve()
    cwd = start if start.is_dir() else start.parent
    result = _git(cwd, "rev-parse", "--show-toplevel")
    root = Path(result.strip()).resolve()
    if not root.is_dir():
        raise WorkspaceError(f"Git project root does not exist: {root}")
    return root


def resolve_in_root(root: str | Path, value: str | Path, *, must_exist: bool = False) -> Path:
    base = Path(root).resolve()
    candidate = Path(value).expanduser()
    candidate = candidate if candidate.is_absolute() else base / candidate
    resolved = candidate.resolve(strict=must_exist)
    try:
        resolved.relative_to(base)
    except ValueError as error:
        raise WorkspaceError(f"path escapes project root: {value}") from error
    return resolved


def relative_to_root(root: str | Path, path: str | Path) -> str:
    return resolve_in_root(root, path).relative_to(Path(root).resolve()).as_posix()


def sha256_file(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def exclusive_json_write(path: str | Path, document: Any) -> Path:
    """Create JSON atomically with O_EXCL; never overwrite an existing file."""

    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    descriptor = os.open(target, flags, 0o644)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(document, handle, indent=2, ensure_ascii=False, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
    except BaseException:
        target.unlink(missing_ok=True)
        raise
    return target


def source_fingerprints(
    root: str | Path, sources: Iterable[str | Path]
) -> tuple[dict[str, Any], ...]:
    base = Path(root).resolve()
    records = []
    for source in sorted({relative_to_root(base, item) for item in sources}):
        path = resolve_in_root(base, source, must_exist=True)
        if not path.is_file():
            raise WorkspaceError(f"source is not a regular file: {source}")
        records.append(
            {
                "path": source,
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    return tuple(records)


def snapshot(root: str | Path, sources: Iterable[str | Path] = ()) -> dict[str, Any]:
    base = Path(root).resolve()
    return {
        "schema_version": "yeastbridge.discovery.workspace.v1",
        "project_root": str(base),
        "git_head": _git(base, "rev-parse", "HEAD").strip(),
        "git_status_porcelain": _git(base, "status", "--porcelain").splitlines(),
        "git_diff": _git(base, "diff", "--no-ext-diff", "--binary"),
        "sources": list(source_fingerprints(base, sources)),
    }


def verify_fingerprints(root: str | Path, records: Iterable[Mapping[str, Any]]) -> list[str]:
    errors: list[str] = []
    for index, record in enumerate(records):
        try:
            path = resolve_in_root(root, str(record["path"]), must_exist=True)
            observed = sha256_file(path)
            if observed != record.get("sha256"):
                errors.append(f"sources[{index}] SHA-256 mismatch")
            if path.stat().st_size != record.get("bytes"):
                errors.append(f"sources[{index}] byte count mismatch")
        except (KeyError, OSError, WorkspaceError) as error:
            errors.append(f"sources[{index}] invalid: {error}")
    return errors


def _git(cwd: Path, *arguments: str) -> str:
    try:
        completed = subprocess.run(
            ["git", *arguments],
            cwd=cwd,
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError) as error:
        detail = getattr(error, "stderr", "") or str(error)
        raise WorkspaceError(f"git {' '.join(arguments)} failed: {detail.strip()}") from error
    return completed.stdout
