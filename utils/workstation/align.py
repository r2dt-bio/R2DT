"""Helpers for workstation Stockholm / R-scape (align) jobs."""

# Job meta shape mirrors draw/pdb creators in jobs.py.
# pylint: disable=duplicate-code

from __future__ import annotations

import hashlib
import re
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional

from utils.workstation.advanced import hash_align_advanced, normalize_align_advanced

if TYPE_CHECKING:
    from utils.workstation.catalog import Catalog
    from utils.workstation.jobs import JobRunner

COV_MARKERS = ("#=GC cov_SS_cons", "#=GC cov_h_SS_cons")

_GF_ID_RE = re.compile(r"^#=GF\s+ID\s+(\S+)", re.MULTILINE)
_GF_AC_RE = re.compile(r"^#=GF\s+AC\s+(\S+)", re.MULTILINE)

# Repo examples served by GET /api/examples/stockholm/<id>
EXAMPLE_STOCKHOLM = {
    "RF00162": "examples/RF00162.stk",
    "hcv": "examples/hcv-alignment.stk",
}


def has_covariation_annotations(stockholm_text: str) -> bool:
    """True when the Stockholm text already has R-scape cov_* GC lines."""
    return any(marker in stockholm_text for marker in COV_MARKERS)


def peek_stockholm_label(stockholm_text: str) -> str:
    """Best-effort display name from #=GF ID / AC."""
    match = _GF_ID_RE.search(stockholm_text)
    if match:
        return match.group(1)
    match = _GF_AC_RE.search(stockholm_text)
    if match:
        return match.group(1)
    return ""


def find_cacofold_r2r_sto(directory: Path) -> Optional[Path]:
    """Return the preferred R-scape CACOfold R2R Stockholm file, if any."""
    hits = sorted(Path(directory).glob("*.cacofold.R2R.sto"))
    return hits[0] if hits else None


def list_align_result_svgs(job_dir: Path) -> List[Path]:
    """SVGs for the results gallery: stitched first, then covariation, then plain."""
    job_dir = Path(job_dir)
    out: List[Path] = []
    stitched = job_dir / "stitched.svg"
    if stitched.is_file():
        out.append(stitched)
    svg_dir = job_dir / "results" / "svg"
    if svg_dir.is_dir():
        cov = sorted(svg_dir.glob("*.covariation.svg"))
        plain = sorted(
            p for p in svg_dir.glob("*.svg") if not p.name.endswith(".covariation.svg")
        )
        out.extend(cov)
        out.extend(plain)
    return out


def svg_gallery_entry(path: Path, job_dir: Path) -> Dict[str, Any]:
    """Metadata for one gallery tile."""
    name = path.name
    if name == "stitched.svg":
        kind = "stitched"
        caption = "Stitched overview"
    elif name.endswith(".covariation.svg"):
        kind = "covariation"
        caption = name.replace(".covariation.svg", "") + " (covariation)"
    else:
        kind = "standard"
        caption = path.stem
    rel = path.relative_to(job_dir).as_posix()
    return {"path": rel, "name": name, "kind": kind, "caption": caption}


def align_content_hash(
    stockholm_text: str,
    stitch: bool,
    advanced: Optional[Dict[str, Any]] = None,
) -> str:
    """Stable hash for one Stockholm alignment + stitch + Advanced flags."""
    digest = hashlib.sha256()
    digest.update(b"align\0")
    hash_align_advanced(digest, stitch, advanced or {})
    digest.update(stockholm_text.encode("utf-8"))
    return "sha256:" + digest.hexdigest()


def create_align_job_from_stockholm(  # pylint: disable=too-many-arguments,too-many-locals
    catalog: "Catalog",
    runner: "JobRunner",
    *,
    stockholm_text: str,
    stitch: bool = True,
    label: str = "",
    notes: str = "",
    force: bool = False,
    advanced: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Create one align job from Stockholm text and enqueue it."""
    # Deferred to avoid circular import with jobs.JobRunner.
    from utils.workstation.catalog import (  # pylint: disable=import-outside-toplevel
        utc_now,
    )
    from utils.workstation.jobs import (  # pylint: disable=import-outside-toplevel
        new_job_id,
    )

    text = stockholm_text.strip()
    if not text:
        raise ValueError("Stockholm alignment is empty")
    if "STOCKHOLM" not in text.upper()[:400]:
        raise ValueError("Input does not look like a Stockholm alignment")

    adv = normalize_align_advanced(advanced)
    digest = align_content_hash(text, stitch, adv)
    if not force:
        existing = catalog.find_by_content_hash(digest, mode="align")
        if existing:
            return {"dedup": True, "job": existing}

    job_id = new_job_id()
    inputs_dir = catalog.inputs_job_dir(job_id)
    inputs_dir.mkdir(parents=True, exist_ok=True)
    stockholm_name = "alignment.stk"
    (inputs_dir / stockholm_name).write_text(text, encoding="utf-8")

    has_cov = has_covariation_annotations(text)
    display = label.strip() or peek_stockholm_label(text) or job_id

    meta: Dict[str, Any] = {
        "id": job_id,
        "mode": "align",
        "label": display,
        "notes": notes.strip(),
        "created": utc_now(),
        "status": "queued",
        "params": {
            "stitch": bool(stitch),
            "has_covariation": has_cov,
            "will_run_rscape": not has_cov,
            "advanced": adv,
        },
        "inputs": {
            "stockholm_name": stockholm_name,
            "content_hash": digest,
        },
        "metrics": {},
        "results_url": f"/jobs/{job_id}/results/",
    }
    catalog.write_meta(job_id, meta)
    catalog.append_log(job_id, "Queued")
    runner.enqueue(job_id)
    return {"dedup": False, "job": meta}
