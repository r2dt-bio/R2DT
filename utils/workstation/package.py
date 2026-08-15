"""Export / import portable ``.r2dt-job.zip`` packages for the workstation."""

from __future__ import annotations

import io
import json
import shutil
import zipfile
from pathlib import Path
from typing import Any, Dict, Iterable, Tuple

from utils.workstation.catalog import Catalog, utc_now
from utils.workstation.chrome import normalize_job_mode
from utils.workstation.edits import edit_counts
from utils.workstation.jobs import new_job_id

SCHEMA_VERSION = 1
PACKAGE_KIND = "r2dt-workstation-job"


def export_job_bytes(catalog: Catalog, job_id: str) -> Tuple[bytes, str]:
    """Build a zip for ``job_id``. Returns ``(bytes, download_filename)``."""
    meta = catalog.read_meta(job_id)
    if not meta:
        raise ValueError(f"Unknown job: {job_id}")

    job_dir = catalog.job_dir(job_id)
    inputs_dir = catalog.inputs_job_dir(job_id)
    mode = normalize_job_mode(meta.get("mode"))
    label = str(meta.get("label") or job_id)
    counts = (
        edit_counts(job_dir)
        if (job_dir / "edits").exists()
        else {
            "ref_count": 0,
            "model_count": 0,
        }
    )

    has_viewer = _has_viewer(job_dir)
    has_results = (job_dir / "results").is_dir() or bool(meta.get("svg_path"))

    package: Dict[str, Any] = {
        "schema": SCHEMA_VERSION,
        "kind": PACKAGE_KIND,
        "mode": mode,
        "label": label,
        "notes": meta.get("notes") or "",
        "exported": utc_now(),
        "source_job_id": job_id,
        "params": meta.get("params") or {},
        "inputs": meta.get("inputs") or {},
        "metrics": meta.get("metrics") or {},
        "edits": {
            **counts,
            "paths": {
                "ref": "edits/ref.overrides.json",
                "model": "edits/model.overrides.json",
            },
        },
        "has_viewer": has_viewer,
        "has_results": has_results,
        "results_url": meta.get("results_url"),
        "viewer_url": meta.get("viewer_url"),
        "svg_path": meta.get("svg_path"),
        "svg_gallery": meta.get("svg_gallery"),
    }

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, mode="w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr(
            "package.json",
            json.dumps(package, indent=2) + "\n",
        )
        # Full meta for round-trip of mode-specific fields.
        export_meta = dict(meta)
        export_meta["exported_from"] = job_id
        zf.writestr("meta.json", json.dumps(export_meta, indent=2) + "\n")

        _zip_tree(zf, inputs_dir, "inputs")
        _zip_viewer(zf, job_dir)
        _zip_tree(zf, job_dir / "edits", "edits")
        _zip_tree(zf, job_dir / "results", "results")
        for name in ("stitched.svg", "log.txt"):
            path = job_dir / name
            if path.is_file():
                zf.write(path, name)
        # Draw/align may store a single svg_path outside results/.
        if svg_rel := meta.get("svg_path"):
            svg_path = job_dir / str(svg_rel)
            if svg_path.is_file() and "results/" not in str(svg_rel):
                zf.write(svg_path, Path(str(svg_rel)).as_posix())

    safe_label = _safe_filename(label) or job_id
    filename = f"{safe_label}.r2dt-job.zip"
    return buf.getvalue(), filename


def import_job_bytes(  # pylint: disable=too-many-locals,too-many-branches,too-many-statements
    catalog: Catalog,
    data: bytes,
    *,
    label: str = "",
    notes: str = "",
) -> Dict[str, Any]:
    """Import a package into a **new** job id. Returns ``{job: meta}``."""
    if not data:
        raise ValueError("Empty package")

    try:
        zf = zipfile.ZipFile(io.BytesIO(data))
    except zipfile.BadZipFile as exc:
        raise ValueError("Not a valid zip file") from exc

    with zf:
        names = zf.namelist()
        _assert_safe_members(names)
        package = _read_json_member(zf, "package.json")
        if package.get("kind") != PACKAGE_KIND:
            raise ValueError("Not an R2DT workstation job package")
        schema = int(package.get("schema") or 0)
        if schema != SCHEMA_VERSION:
            raise ValueError(
                f"Unsupported package schema {schema} (expected {SCHEMA_VERSION})"
            )

        mode = normalize_job_mode(package.get("mode") or "compare")
        job_id = new_job_id()
        job_dir = catalog.job_dir(job_id)
        inputs_dir = catalog.inputs_job_dir(job_id)
        job_dir.mkdir(parents=True, exist_ok=True)
        inputs_dir.mkdir(parents=True, exist_ok=True)

        _extract_prefix(zf, "inputs/", inputs_dir)
        _extract_prefix(zf, "viewer/", job_dir / "viewer")
        _extract_prefix(zf, "edits/", job_dir / "edits")
        _extract_prefix(zf, "results/", job_dir / "results")
        for name in ("stitched.svg", "log.txt"):
            if name in names:
                _extract_file(zf, name, job_dir / name)

        # Optional meta.json supplies extras; package.json wins for identity.
        old_meta: Dict[str, Any] = {}
        if "meta.json" in names:
            try:
                old_meta = _read_json_member(zf, "meta.json")
            except ValueError:
                old_meta = {}

        has_viewer = (job_dir / "viewer" / "index.html").is_file()
        has_results = (job_dir / "results").is_dir() or bool(
            package.get("svg_path") or old_meta.get("svg_path")
        )
        if has_viewer or (mode in ("draw", "align") and has_results):
            status = "ready"
        elif has_results or _inputs_nonempty(inputs_dir):
            status = "needs_generate"
        else:
            status = "ready" if package.get("has_viewer") else "needs_generate"

        display = (label or "").strip() or str(
            package.get("label") or old_meta.get("label") or job_id
        )
        note_text = (notes or "").strip()
        if not note_text:
            note_text = str(package.get("notes") or old_meta.get("notes") or "")
        imported_note = f"Imported from {package.get('source_job_id') or 'package'}"
        if note_text and imported_note not in note_text:
            note_text = f"{note_text} · {imported_note}"
        elif not note_text:
            note_text = imported_note

        meta: Dict[str, Any] = {
            "id": job_id,
            "mode": mode,
            "label": display,
            "notes": note_text,
            "created": utc_now(),
            "finished": utc_now() if status == "ready" else None,
            "status": status,
            "params": package.get("params") or old_meta.get("params") or {},
            "inputs": package.get("inputs") or old_meta.get("inputs") or {},
            "metrics": package.get("metrics") or old_meta.get("metrics") or {},
            "error": None,
            "source": {
                "kind": "package-import",
                "package_source_job_id": package.get("source_job_id"),
                "imported": utc_now(),
            },
        }
        if has_viewer:
            meta["viewer_url"] = f"/jobs/{job_id}/viewer/"
        if mode in ("draw", "align") or package.get("results_url") or has_results:
            meta["results_url"] = f"/jobs/{job_id}/results/"
        for key in ("svg_path", "svg_gallery"):
            if package.get(key) is not None:
                meta[key] = package[key]
            elif old_meta.get(key) is not None:
                meta[key] = old_meta[key]

        counts = edit_counts(job_dir)
        if counts["ref_count"] or counts["model_count"]:
            meta["edits"] = counts

        catalog.write_meta(job_id, meta)
        catalog.append_log(job_id, f"Imported package ({status})")
        return {"job": meta, "dedup": False}


def _has_viewer(job_dir: Path) -> bool:
    viewer = job_dir / "viewer"
    if viewer.is_symlink() or viewer.is_dir():
        return (viewer / "index.html").is_file() or (
            viewer.resolve() / "index.html"
        ).is_file()
    return False


def _inputs_nonempty(inputs_dir: Path) -> bool:
    if not inputs_dir.is_dir():
        return False
    return any(p.is_file() for p in inputs_dir.iterdir())


def _safe_filename(label: str) -> str:
    cleaned = "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in label)
    cleaned = cleaned.strip("._")
    return cleaned[:80]


def _assert_safe_members(names: Iterable[str]) -> None:
    for name in names:
        if name.startswith("/") or name.startswith("\\") or ".." in Path(name).parts:
            raise ValueError(f"Unsafe path in package: {name}")


def _read_json_member(zf: zipfile.ZipFile, name: str) -> Dict[str, Any]:
    try:
        raw = zf.read(name)
    except KeyError as exc:
        raise ValueError(f"Package missing {name}") from exc
    try:
        data = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, ValueError) as exc:
        raise ValueError(f"Invalid JSON in {name}") from exc
    if not isinstance(data, dict):
        raise ValueError(f"{name} must be a JSON object")
    return data


def _zip_tree(zf: zipfile.ZipFile, src: Path, arc_prefix: str) -> None:
    if not src.is_dir():
        return
    for path in sorted(src.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(src).as_posix()
        zf.write(path, f"{arc_prefix}/{rel}")


def _zip_viewer(zf: zipfile.ZipFile, job_dir: Path) -> None:
    viewer = job_dir / "viewer"
    if not viewer.exists() and not viewer.is_symlink():
        return
    try:
        real = viewer.resolve()
    except OSError:
        return
    if not real.is_dir():
        return
    for path in sorted(real.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(real).as_posix()
        zf.write(path, f"viewer/{rel}")


def _extract_prefix(zf: zipfile.ZipFile, prefix: str, dest: Path) -> None:
    members = [n for n in zf.namelist() if n.startswith(prefix) and not n.endswith("/")]
    if not members:
        return
    dest.mkdir(parents=True, exist_ok=True)
    for name in members:
        rel = name[len(prefix) :]
        if not rel or ".." in Path(rel).parts:
            continue
        out = dest / rel
        out.parent.mkdir(parents=True, exist_ok=True)
        with zf.open(name) as src, out.open("wb") as handle:
            shutil.copyfileobj(src, handle)


def _extract_file(zf: zipfile.ZipFile, name: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with zf.open(name) as src, dest.open("wb") as handle:
        shutil.copyfileobj(src, handle)
