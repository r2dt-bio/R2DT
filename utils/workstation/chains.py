"""RNA chain auto-detection for uploaded structures."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List

from utils import fr3d as fr3d_utils

ALLOWED_SUFFIXES = {".cif", ".pdb", ".cif.gz", ".pdb.gz"}


def normalize_suffix(path: Path) -> str:
    """Return a lowercase structure suffix, including ``.gz`` when present."""
    name = path.name.lower()
    for suffix in (".cif.gz", ".pdb.gz", ".cif", ".pdb"):
        if name.endswith(suffix):
            return suffix
    return path.suffix.lower()


def is_structure_filename(filename: str) -> bool:
    """True if filename looks like a PDB/mmCIF (optionally gzipped)."""
    return normalize_suffix(Path(filename)) in ALLOWED_SUFFIXES


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
    # Prefer mmCIF for compare (reference must be cif in the pipeline).
    fmt = "cif" if suffix.startswith(".cif") else "pdb"
    return {
        "filename": path.name,
        "format": fmt,
        "suffix": suffix,
        "model": model_id,
        "chains": chains,
        "compare_ready": fmt == "cif",
    }
