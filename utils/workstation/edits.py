"""Durable base-pair override files for workstation jobs."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List

VALID_ACTIONS = frozenset({"add", "delete", "refamily"})


def apply_overrides(
    baseline_anns: List[Dict[str, Any]], overrides: List[Dict[str, Any]]
) -> List[Dict[str, Any]]:
    """Return a deep-copied annotation list with overrides applied.

    Mirrors ``R2DTBpEdit._applyOverrides`` in the workstation edit UI.
    """
    anns = [dict(a) for a in (baseline_anns or [])]
    for op in overrides or []:
        try:
            i = int(op["i"])
            j = int(op["j"])
        except (KeyError, TypeError, ValueError):
            continue
        if i == j:
            continue
        if i > j:
            i, j = j, i
        idx = _find_ann_index(anns, i, j)
        action = op.get("action")
        if action == "delete":
            if idx >= 0:
                anns.pop(idx)
        elif action == "add":
            family = str(op.get("family") or "cWW")
            if idx < 0:
                anns.append(
                    {
                        "seq_id1": str(i),
                        "seq_id2": str(j),
                        "3d_id1": str(i),
                        "3d_id2": str(j),
                        "nt1": "N",
                        "nt2": "N",
                        "unit1": "N",
                        "unit2": "N",
                        "bp": family,
                        "crossing": "0",
                    }
                )
            else:
                anns[idx]["bp"] = family
        elif action == "refamily" and idx >= 0:
            family = op.get("to") or op.get("family") or anns[idx].get("bp")
            anns[idx]["bp"] = str(family)
    return anns


def annotations_to_pairs(anns: List[Dict[str, Any]]) -> List[tuple]:
    """Convert fr3d annotations to ``(i, j, family)`` tuples for INF."""
    pairs = []
    for ann in anns or []:
        try:
            i = int(ann["seq_id1"])
            j = int(ann["seq_id2"])
        except (KeyError, TypeError, ValueError):
            continue
        if i == j:
            continue
        if i > j:
            i, j = j, i
        family = str(ann.get("bp") or "cWW")
        pairs.append((i, j, family))
    return pairs


def _find_ann_index(anns: List[Dict[str, Any]], i: int, j: int) -> int:
    for n, ann in enumerate(anns):
        try:
            a = int(ann["seq_id1"])
            b = int(ann["seq_id2"])
        except (KeyError, TypeError, ValueError):
            continue
        if (a == i and b == j) or (a == j and b == i):
            return n
    return -1


def edits_dir(job_dir: Path) -> Path:
    """Return ``jobs/<id>/edits/``, creating it if needed."""
    path = Path(job_dir) / "edits"
    path.mkdir(parents=True, exist_ok=True)
    return path


def overrides_path(job_dir: Path, panel: str) -> Path:
    """Path to ``ref.overrides.json`` or ``model.overrides.json``."""
    if panel not in ("ref", "model"):
        raise ValueError(f"panel must be 'ref' or 'model', got {panel!r}")
    return edits_dir(job_dir) / f"{panel}.overrides.json"


def load_overrides(job_dir: Path) -> Dict[str, List[Dict[str, Any]]]:
    """Load both panels' override lists (empty if missing)."""
    return {
        "ref": _read_list(overrides_path(job_dir, "ref")),
        "model": _read_list(overrides_path(job_dir, "model")),
    }


def save_overrides(
    job_dir: Path, payload: Dict[str, Any]
) -> Dict[str, List[Dict[str, Any]]]:
    """Validate and write both panels. Returns the normalised payload."""
    ref = _validate_list(payload.get("ref") or [])
    model = _validate_list(payload.get("model") or [])
    _write_list(overrides_path(job_dir, "ref"), ref)
    _write_list(overrides_path(job_dir, "model"), model)
    return {"ref": ref, "model": model}


def edit_counts(job_dir: Path) -> Dict[str, int]:
    """Return ``{ref_count, model_count}`` for catalog badges."""
    data = load_overrides(job_dir)
    return {"ref_count": len(data["ref"]), "model_count": len(data["model"])}


def _read_list(path: Path) -> List[Dict[str, Any]]:
    if not path.is_file():
        return []
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return []
    if not isinstance(raw, list):
        return []
    try:
        return _validate_list(raw)
    except ValueError:
        return []


def _write_list(path: Path, rows: List[Dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")


def _validate_list(rows: Any) -> List[Dict[str, Any]]:
    if not isinstance(rows, list):
        raise ValueError("overrides must be a list")
    return [_validate_one(idx, row) for idx, row in enumerate(rows)]


def _validate_one(idx: int, row: Any) -> Dict[str, Any]:
    if not isinstance(row, dict):
        raise ValueError(f"override[{idx}] must be an object")
    action = row.get("action")
    if action not in VALID_ACTIONS:
        raise ValueError(f"override[{idx}].action invalid: {action!r}")
    try:
        i = int(row["i"])
        j = int(row["j"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"override[{idx}] needs integer i and j") from exc
    if i == j:
        raise ValueError(f"override[{idx}] i and j must differ")
    # Canonicalise order for storage stability.
    if i > j:
        i, j = j, i
    rec: Dict[str, Any] = {"action": action, "i": i, "j": j}
    if action == "add":
        family = str(row.get("family") or "").strip()
        if not family:
            raise ValueError(f"override[{idx}] add needs family")
        rec["family"] = family
    elif action == "delete":
        family = row.get("family")
        if family is not None and str(family).strip():
            rec["family"] = str(family).strip()
    else:  # refamily
        to_fam = str(row.get("to") or row.get("family") or "").strip()
        if not to_fam:
            raise ValueError(f"override[{idx}] refamily needs to/family")
        rec["to"] = to_fam
        from_fam = row.get("from")
        if from_fam is not None and str(from_fam).strip():
            rec["from"] = str(from_fam).strip()
    return rec
