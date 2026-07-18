"""RNA chain auto-detection for uploaded structures."""

from __future__ import annotations

import gzip
import shutil
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional

from utils import fr3d as fr3d_utils

ALLOWED_SUFFIXES = {".cif", ".pdb", ".cif.gz", ".pdb.gz"}


def normalize_suffix(path: Path) -> str:
    """Return a lowercase structure suffix, including ``.gz`` when present."""
    name = path.name.lower()
    for suffix in (".cif.gz", ".pdb.gz", ".cif", ".pdb"):
        if name.endswith(suffix):
            return suffix
    return path.suffix.lower()


def structure_stem(path: Path) -> str:
    """Basename without ``.pdb`` / ``.cif`` / ``.gz`` suffixes."""
    name = path.name
    lower = name.lower()
    for suffix in (".cif.gz", ".pdb.gz", ".cif", ".pdb"):
        if lower.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def is_structure_filename(filename: str) -> bool:
    """True if filename looks like a PDB/mmCIF (optionally gzipped)."""
    return normalize_suffix(Path(filename)) in ALLOWED_SUFFIXES


def ensure_mmcif(path: Path, out_dir: Path, label: Optional[str] = None) -> Path:
    """Return an mmCIF path for ``path``, converting PDB when needed.

    Compare / multi-chain mode reads structures through FR3D's mmCIF reader
    (same constraint as ``r2dt.py pdb --compare`` and CASP ``ensure_cif``).
    Plain ``.cif`` is copied as-is; ``.pdb`` / gzipped inputs are written as
    ``{label}.cif`` under ``out_dir``.
    """
    # pylint: disable=import-outside-toplevel
    path = Path(path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    suffix = normalize_suffix(path)
    stem = label or structure_stem(path)
    out = out_dir / f"{stem}.cif"

    if suffix == ".cif":
        if path.resolve() != out.resolve():
            shutil.copy2(path, out)
        return out

    if suffix == ".cif.gz":
        with gzip.open(path, "rb") as src, out.open("wb") as dest:
            shutil.copyfileobj(src, dest)
        return out

    from Bio.PDB import MMCIFIO, PDBParser

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        parser = PDBParser(QUIET=True)
        if suffix == ".pdb.gz":
            with gzip.open(path, "rt") as handle:
                structure = parser.get_structure(stem, handle)
        else:
            structure = parser.get_structure(stem, str(path))
        io = MMCIFIO()
        io.set_structure(structure)
        io.save(str(out))
    return out


def list_rna_chains(structure_path: Path) -> Dict[str, Any]:
    """Return RNA chain ids for the first model of a structure file.

    Uses ``fr3d.get_structure_info`` (FR3D for mmCIF, BioPython for PDB).
    """
    path = Path(structure_path)
    suffix = normalize_suffix(path)
    info = fr3d_utils.get_structure_info(str(path))
    models = info.get("models") or []
    chains_by_model = info.get("chains") or {}
    model_id = models[0] if models else None
    chains: List[str] = []
    if model_id is not None:
        chains = list(chains_by_model.get(model_id) or [])
    fmt = "cif" if suffix.startswith(".cif") else "pdb"
    return {
        "filename": path.name,
        "format": fmt,
        "suffix": suffix,
        "model": model_id,
        "chains": chains,
        # PDB refs are auto-converted to mmCIF at job create time.
        "compare_ready": True,
        "needs_cif_conversion": fmt == "pdb",
    }
