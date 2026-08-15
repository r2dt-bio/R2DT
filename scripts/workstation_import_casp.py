#!/usr/bin/env python3
"""Import CASP15/16 compare viewers into the local R2DT workstation catalog.

Reads ``output/site/casp{15,16}/results.json`` and, for each ok/cached pair,
creates a ready compare job under ``~/.r2dt-workstation`` (or ``--workspace``)
whose ``viewer/`` is a symlink into the site tree.

Symlinks use the in-container repo path ``/rna/r2dt/...`` so they resolve when
serving via ``just workstation`` (repo mounted at ``/rna/r2dt``). Pass
``--link-style host`` to point at the absolute host path instead.

Idempotent: re-running updates meta and recreates the viewer link.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

REPO_IN_DOCKER = Path("/rna/r2dt")


def utc_now() -> str:
    """Return an ISO-8601 UTC timestamp."""
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def job_id_for(season: str, entry: Dict[str, Any]) -> str:
    """Stable filesystem-safe id, e.g. casp15-R1107-01-R1107TS232_1."""
    rank = int(entry.get("rank") or 0)
    model = str(entry.get("model") or "model")
    target = str(entry.get("target") or "target")
    return f"{season}-{target}-{rank:02d}-{model}"


def pair_dir(site: Path, entry: Dict[str, Any]) -> Path:
    """Directory that contains ``viewer/`` for this results.json row."""
    page = entry.get("page") or ""
    # page: Target/01-Model/viewer/index.html
    return site / Path(page).parent.parent


def metrics_from_entry(entry: Dict[str, Any]) -> Dict[str, Any]:
    """Map results.json INF/diff fields into workstation meta.metrics shape."""
    inf = entry.get("inf") or {}
    diff = entry.get("diff") or {}
    return {
        "inf": {
            "wc": inf.get("wc"),
            "nwc": inf.get("nwc"),
            "all": inf.get("all"),
        },
        "matched": diff.get("matched"),
        "lost": diff.get("lost"),
        "added": diff.get("added"),
        "superpose_rmsd": entry.get("rmsd"),
    }


def chains_csv(entry: Dict[str, Any]) -> str:
    """Join chain ids from a results row into a comma-separated string."""
    chains = entry.get("chains") or []
    if isinstance(chains, list):
        return ",".join(str(c) for c in chains)
    return str(chains)


def link_target(
    *,
    style: str,
    repo: Path,
    pair: Path,
) -> Path:
    """Return the path stored in the viewer symlink."""
    rel_from_repo = pair.relative_to(repo) / "viewer"
    if style == "docker":
        return REPO_IN_DOCKER / rel_from_repo
    return (repo / rel_from_repo).resolve()


def ensure_viewer_link(job_dir: Path, target: Path) -> None:
    """Create or replace ``job_dir/viewer`` as a symlink to ``target``."""
    link = job_dir / "viewer"
    if link.is_symlink() or link.exists():
        if link.is_symlink() or link.is_file():
            link.unlink()
        else:
            # Refuse to delete a real directory tree of generated assets.
            raise RuntimeError(
                f"{link} exists and is not a symlink — remove it manually first"
            )
    link.symlink_to(target)


def import_season(
    *,
    repo: Path,
    workspace: Path,
    season: str,
    link_style: str,
    force: bool,
) -> Tuple[int, int, int]:
    """Return (created_or_updated, skipped, missing)."""
    site = repo / "output" / "site" / season
    results_path = site / "results.json"
    if not results_path.is_file():
        print(f"[skip] {season}: no {results_path}", file=sys.stderr)
        return 0, 0, 0

    payload = json.loads(results_path.read_text(encoding="utf-8"))
    entries: List[Dict[str, Any]] = list(payload.get("results") or [])
    jobs_dir = workspace / "jobs"
    jobs_dir.mkdir(parents=True, exist_ok=True)

    done = skipped = missing = 0
    for entry in entries:
        status = entry.get("status")
        if status not in ("ok", "cached"):
            skipped += 1
            continue
        pair = pair_dir(site, entry)
        viewer_index = pair / "viewer" / "index.html"
        if not viewer_index.is_file():
            print(f"[missing] {season} {entry.get('page')}")
            missing += 1
            continue

        jid = job_id_for(season, entry)
        job_dir = jobs_dir / jid
        meta_path = job_dir / "meta.json"
        if meta_path.is_file() and not force:
            # Still refresh the symlink in case the site moved.
            pass

        job_dir.mkdir(parents=True, exist_ok=True)
        target = link_target(style=link_style, repo=repo, pair=pair)
        ensure_viewer_link(job_dir, target)

        rank = int(entry.get("rank") or 0)
        label = (
            f"{season.upper()} · {entry.get('target')} · #{rank} · {entry.get('model')}"
        )
        chains = chains_csv(entry)
        created = utc_now()
        if meta_path.is_file():
            try:
                old = json.loads(meta_path.read_text(encoding="utf-8"))
                created = old.get("created") or created
            except (OSError, ValueError):
                pass

        meta: Dict[str, Any] = {
            "id": jid,
            "mode": "compare",
            "label": label,
            "notes": f"Imported from output/site/{season}",
            "created": created,
            "finished": utc_now(),
            "status": "ready",
            "params": {
                "chains": chains,
                "model_chains": chains,
                "mode": "auto",
                "basepairs": "fr3d",
                "casp_season": season,
                "casp_target": entry.get("target"),
                "casp_rank": rank,
                "casp_tm": entry.get("tm"),
            },
            "inputs": {
                "ref_name": Path(str(entry.get("reference") or "")).name,
                "model_name": str(entry.get("model") or ""),
                "content_hash": f"casp-import:{season}:{entry.get('page')}",
            },
            "metrics": metrics_from_entry(entry),
            "viewer_url": f"/jobs/{jid}/viewer/",
            "error": None,
            "source": {
                "kind": "casp-import",
                "season": season,
                "page": entry.get("page"),
                "viewer": str(target),
            },
        }
        meta_path.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
        done += 1
        print(f"[ok] {jid}")

    return done, skipped, missing


def refresh_catalog(workspace: Path) -> int:
    """Rewrite catalog.json from jobs/*/meta.json (newest first)."""
    jobs_dir = workspace / "jobs"
    jobs: List[Dict[str, Any]] = []
    if jobs_dir.is_dir():
        for child in sorted(jobs_dir.iterdir(), reverse=True):
            meta_path = child / "meta.json"
            if not meta_path.is_file():
                continue
            try:
                meta = json.loads(meta_path.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                continue
            if isinstance(meta, dict):
                jobs.append(meta)
    catalog = {"updated": utc_now(), "jobs": jobs}
    (workspace / "catalog.json").write_text(
        json.dumps(catalog, indent=2) + "\n", encoding="utf-8"
    )
    return len(jobs)


def main(argv: Optional[List[str]] = None) -> int:
    """CLI entry: import CASP site pairs into the workstation workspace."""
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--repo",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="R2DT repository root",
    )
    ap.add_argument(
        "--workspace",
        type=Path,
        default=None,
        help="Workstation cache (default: ~/.r2dt-workstation)",
    )
    ap.add_argument(
        "--seasons",
        default="casp15,casp16",
        help="Comma-separated site folder names under output/site/",
    )
    ap.add_argument(
        "--link-style",
        choices=("docker", "host"),
        default="docker",
        help="Symlink target path style (docker=/rna/r2dt/... for just workstation)",
    )
    ap.add_argument(
        "--force",
        action="store_true",
        help="Rewrite meta even when the job already exists",
    )
    args = ap.parse_args(argv)

    repo = args.repo.resolve()
    workspace = (
        args.workspace.expanduser().resolve()
        if args.workspace
        else (Path.home() / ".r2dt-workstation").resolve()
    )
    workspace.mkdir(parents=True, exist_ok=True)

    seasons = [s.strip() for s in args.seasons.split(",") if s.strip()]
    total_done = total_skip = total_miss = 0
    for season in seasons:
        done, skipped, missing = import_season(
            repo=repo,
            workspace=workspace,
            season=season,
            link_style=args.link_style,
            force=args.force,
        )
        total_done += done
        total_skip += skipped
        total_miss += missing
        print(
            f"{season}: imported {done}, skipped status {skipped}, missing viewers {missing}"
        )

    n = refresh_catalog(workspace)
    print(
        f"\nWorkspace {workspace}: {total_done} jobs written, "
        f"catalog has {n} job(s) total "
        f"(link-style={args.link_style})"
    )
    if total_miss:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
