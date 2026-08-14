"""Bake workstation edits into a static shareable viewer zip."""

from __future__ import annotations

import io
import json
import re
import zipfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from utils.multichain import build_inf_report, compute_inf, inf_report_to_csv
from utils.viewer_html import _format_inf_metrics
from utils.workstation.catalog import Catalog, utc_now
from utils.workstation.edits import (
    annotations_to_pairs,
    apply_overrides,
    edit_counts,
    load_overrides,
)
from utils.workstation.package import _has_viewer, _safe_filename

# Relative fr3d.json path inside viewer/ → which override panel to apply.
_FR3D_PANEL = (
    ("fr3d.json", "ref"),
    ("ref/fr3d.json", "ref"),
    ("model/fr3d.json", "model"),
    ("model-own/fr3d.json", "model"),
)

_README = """R2DT shareable viewer
=====================

This zip is a static snapshot of an R2DT interactive viewer with any
workstation base-pair edits already baked into the annotation files
(fr3d.json). No R2DT workstation is required to host it.

1. Unzip this archive.
2. Upload the entire folder to any static HTTPS host (GitHub Pages,
   Cloudflare Pages, S3, nginx, …), keeping relative paths intact.
3. Open index.html, or embed it:

   <iframe src="https://your.cdn/…/index.html" width="1232" height="700"
           style="border:0" allow="fullscreen"></iframe>

   Or load the scripts on your own page and call R2DTViewer.create /
   R2DTViewer.createCompare with baseUrl pointing at this folder
   (see docs/pdb-2d-3d-viewer.md).

The 3D pane loads pdbe-molstar from jsDelivr, so visitors need network
access to that CDN.
"""


def export_shareable_html_bytes(  # pylint: disable=too-many-branches
    catalog: Catalog, job_id: str
) -> Tuple[bytes, str]:
    """Build a baked viewer zip. Returns ``(bytes, download_filename)``."""
    meta = catalog.read_meta(job_id)
    if not meta:
        raise ValueError(f"Unknown job: {job_id}")
    job_dir = catalog.job_dir(job_id)
    if not _has_viewer(job_dir):
        raise ValueError("No interactive viewer for this job")

    viewer = (job_dir / "viewer").resolve()
    if not viewer.is_dir():
        raise ValueError("No interactive viewer for this job")

    overrides = load_overrides(job_dir)
    counts = edit_counts(job_dir)
    baked_fr3d = _bake_fr3d_files(viewer, overrides)
    metrics, inf_report = _recompute_metrics(viewer, baked_fr3d, overrides)

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, mode="w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("README.txt", _README)
        zf.writestr(
            "publish.json",
            json.dumps(
                {
                    "kind": "r2dt-shareable-viewer",
                    "schema": 1,
                    "source_job_id": job_id,
                    "label": meta.get("label") or job_id,
                    "published": utc_now(),
                    "edits": counts,
                },
                indent=2,
            )
            + "\n",
        )
        written_inf = False
        for path in sorted(viewer.rglob("*")):
            if not path.is_file():
                continue
            rel = path.relative_to(viewer).as_posix()
            if rel in baked_fr3d:
                payload = baked_fr3d[rel]
                zf.writestr(rel, json.dumps(payload, indent=2) + "\n")
            elif rel == "metrics.json" and metrics is not None:
                zf.writestr(rel, json.dumps(metrics, indent=2) + "\n")
            elif rel == "inf-pairs.json" and inf_report is not None:
                zf.writestr(rel, json.dumps(inf_report, indent=2) + "\n")
                written_inf = True
            elif rel == "inf-pairs.csv" and inf_report is not None:
                zf.writestr(rel, inf_report_to_csv(inf_report))
            elif rel == "index.html" and metrics is not None:
                html = path.read_text(encoding="utf-8")
                zf.writestr(
                    rel,
                    _rewrite_inf_bar(
                        html,
                        metrics.get("inf"),
                        scopes=metrics.get("scopes"),
                    ),
                )
            elif rel == "manifest.json":
                data = _stamp_manifest(path, job_id, counts)
                zf.writestr(rel, json.dumps(data, indent=2) + "\n")
            else:
                zf.write(path, rel)
        if inf_report is not None and not written_inf:
            zf.writestr("inf-pairs.json", json.dumps(inf_report, indent=2) + "\n")
            zf.writestr("inf-pairs.csv", inf_report_to_csv(inf_report))

    label = str(meta.get("label") or job_id)
    filename = f"{_safe_filename(label) or job_id}.r2dt-viewer.zip"
    return buf.getvalue(), filename


def _bake_fr3d_files(
    viewer: Path, overrides: Dict[str, List[Dict[str, Any]]]
) -> Dict[str, Dict[str, Any]]:
    """Return ``{relative_path: baked_fr3d_dict}`` for files that need rewriting."""
    out: Dict[str, Dict[str, Any]] = {}
    for rel, panel in _FR3D_PANEL:
        path = viewer / rel
        if not path.is_file():
            continue
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        if not isinstance(data, dict):
            continue
        anns = data.get("annotations")
        if not isinstance(anns, list):
            continue
        data = dict(data)
        data["annotations"] = apply_overrides(anns, overrides.get(panel) or [])
        out[rel] = data
    return out


def _recompute_metrics(
    viewer: Path,
    baked_fr3d: Dict[str, Dict[str, Any]],
    overrides: Dict[str, List[Dict[str, Any]]],
) -> Tuple[Optional[Dict[str, Any]], Optional[Dict[str, Any]]]:
    """Update metrics + INF report from baked annotations when both panels exist."""
    metrics_path = viewer / "metrics.json"
    if not metrics_path.is_file():
        return None, None
    try:
        metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None, None
    if not isinstance(metrics, dict):
        return None, None

    ref_data = baked_fr3d.get("ref/fr3d.json")
    model_data = baked_fr3d.get("model/fr3d.json")
    if ref_data is None:
        ref_data = _read_fr3d(viewer / "ref" / "fr3d.json")
        if ref_data and overrides.get("ref"):
            ref_data = dict(ref_data)
            ref_data["annotations"] = apply_overrides(
                ref_data.get("annotations") or [], overrides["ref"]
            )
    if model_data is None:
        model_data = _read_fr3d(viewer / "model" / "fr3d.json")
        if model_data and overrides.get("model"):
            model_data = dict(model_data)
            model_data["annotations"] = apply_overrides(
                model_data.get("annotations") or [], overrides["model"]
            )
    if not ref_data or not model_data:
        return metrics, None

    ref_pairs = annotations_to_pairs(ref_data.get("annotations") or [])
    model_pairs = annotations_to_pairs(model_data.get("annotations") or [])
    metrics = dict(metrics)
    metrics["inf"] = compute_inf(ref_pairs, model_pairs)
    metrics["published"] = utc_now()

    inf_report = None
    if boundaries := _boundaries_from_metrics(metrics):
        inf_report = build_inf_report(
            structure_id=str(metrics.get("structure_id") or ""),
            model_id=str(metrics.get("model_id") or ""),
            chains=[cid for cid, _s, _e in boundaries],
            boundaries=boundaries,
            ref_pairs=ref_pairs,
            model_pairs=model_pairs,
            inf=metrics["inf"],
            one_based=True,
            extra={
                "model_simulated": metrics.get("model_simulated"),
                "model_own_layout": metrics.get("model_own_layout"),
                "published": metrics.get("published"),
            },
        )
        metrics["scopes"] = [
            {
                "id": s.get("id"),
                "label": s.get("label"),
                "type": s.get("type"),
                "chains": s.get("chains"),
                "inf": s.get("inf"),
            }
            for s in (inf_report.get("scopes") or [])
        ]
    return metrics, inf_report


def _boundaries_from_metrics(
    metrics: Dict[str, Any],
) -> List[Tuple[str, int, int]]:
    """Parse 0-based boundaries from metrics.json, if present."""
    raw = metrics.get("boundaries")
    out: List[Tuple[str, int, int]] = []
    if isinstance(raw, list):
        for entry in raw:
            if not isinstance(entry, dict):
                continue
            try:
                out.append(
                    (str(entry["chain"]), int(entry["start"]), int(entry["end"]))
                )
            except (KeyError, TypeError, ValueError):
                continue
    if out:
        return out
    chains = metrics.get("chains")
    if isinstance(chains, list) and len(chains) == 1:
        # Single-chain fallback when older metrics.json lacks boundaries.
        return [(str(chains[0]), 0, 10**9)]
    return out


def _read_fr3d(path: Path) -> Optional[Dict[str, Any]]:
    if not path.is_file():
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return data if isinstance(data, dict) else None


def _rewrite_inf_bar(
    html: str,
    inf: Optional[Dict[str, Any]],
    scopes: Optional[List[Dict[str, Any]]] = None,
) -> str:
    """Replace the compare INF bar in index.html with recomputed values."""
    if not inf:
        return html
    new_bar = _format_inf_metrics(inf, scopes=scopes)
    if not new_bar:
        return html
    updated, n = re.subn(
        r'<div class="mc-inf">.*?</div>',
        new_bar,
        html,
        count=1,
        flags=re.DOTALL,
    )
    return updated if n else html


def _stamp_manifest(path: Path, job_id: str, counts: Dict[str, int]) -> Dict[str, Any]:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        data = {}
    if not isinstance(data, dict):
        data = {}
    data = dict(data)
    data["published"] = {
        "at": utc_now(),
        "source_job_id": job_id,
        "edits": counts,
    }
    return data
