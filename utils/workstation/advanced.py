"""Advanced CLI flag helpers for workstation job create + runners."""

from __future__ import annotations

from typing import Any, Dict, List

TF_ENGINES = frozenset({"auto", "rscape", "rnartist", "rnapuzzler"})
COLOR_BY = frozenset({"none", "structure", "region"})


def normalize_draw_advanced(raw: Any) -> Dict[str, Any]:
    """Validate draw/templatefree Advanced fields from an API payload."""
    data = raw if isinstance(raw, dict) else {}
    force_template = str(data.get("force_template") or "").strip()
    tf_engine = str(data.get("tf_engine") or "auto").strip().lower()
    if tf_engine not in TF_ENGINES:
        raise ValueError("tf_engine must be auto, rscape, rnartist, or rnapuzzler")
    return {
        "force_template": force_template,
        "constraint": bool(data.get("constraint")),
        "skip_ribovore_filters": bool(data.get("skip_ribovore_filters")),
        "tf_engine": tf_engine,
    }


def normalize_structure_advanced(raw: Any) -> Dict[str, Any]:
    """Validate pdb / compare Advanced fields (pseudoknots, rnapuzzler)."""
    data = raw if isinstance(raw, dict) else {}
    # Default on — match CLI --pseudoknots/--no-pseudoknots default=True
    if "pseudoknots" not in data:
        pseudoknots = True
    else:
        pseudoknots = bool(data.get("pseudoknots"))
    return {
        "pseudoknots": pseudoknots,
        "rnapuzzler": bool(data.get("rnapuzzler")),
    }


def normalize_align_advanced(raw: Any) -> Dict[str, Any]:
    """Validate stockholm Advanced fields."""
    data = raw if isinstance(raw, dict) else {}
    color_by = str(data.get("color_by") or "none").strip().lower()
    if color_by not in COLOR_BY:
        raise ValueError("color_by must be none, structure, or region")
    try:
        max_unpaired = int(data.get("max_unpaired", 30))
    except (TypeError, ValueError) as exc:
        raise ValueError("max_unpaired must be an integer") from exc
    if max_unpaired < 0:
        raise ValueError("max_unpaired must be >= 0")
    if "monochrome" not in data:
        monochrome = True
    else:
        monochrome = bool(data.get("monochrome"))
    return {
        "include_novel": bool(data.get("include_novel")),
        "all_structures": bool(data.get("all_structures")),
        "auto_repair": bool(data.get("auto_repair")),
        "max_unpaired": max_unpaired,
        "monochrome": monochrome,
        "color_by": color_by,
    }


def append_draw_flags(cmd: List[str], layout: str, adv: Dict[str, Any]) -> None:
    """Extend a draw/templatefree argv with Advanced flags."""
    if layout != "templatefree":
        if adv.get("force_template"):
            cmd.extend(["--force_template", str(adv["force_template"])])
        if adv.get("constraint"):
            cmd.append("--constraint")
        if adv.get("skip_ribovore_filters"):
            cmd.append("--skip_ribovore_filters")
        return
    engine = adv.get("tf_engine") or "auto"
    if engine == "rscape":
        cmd.append("--rscape")
    elif engine == "rnartist":
        cmd.append("--rnartist")
    elif engine == "rnapuzzler":
        cmd.append("--rnapuzzler")
    # auto: leave defaults (templatefree --auto is default behaviour)


def append_structure_flags(cmd: List[str], adv: Dict[str, Any]) -> None:
    """Extend a pdb / pdb_2d_3d argv with Advanced flags."""
    if adv.get("pseudoknots", True):
        cmd.append("--pseudoknots")
    else:
        cmd.append("--no-pseudoknots")
    if adv.get("rnapuzzler"):
        cmd.append("--rnapuzzler")


def append_stockholm_flags(
    cmd: List[str], *, stitch: bool, adv: Dict[str, Any]
) -> None:
    """Extend a stockholm argv with stitch + Advanced flags."""
    cmd.append("--stitch" if stitch else "--no-stitch")
    if adv.get("include_novel"):
        cmd.append("--include-novel")
    if adv.get("all_structures"):
        cmd.append("--all-structures")
    else:
        cmd.append("--named-only")
    if adv.get("auto_repair"):
        cmd.append("--auto-repair")
    else:
        cmd.append("--no-auto-repair")
    cmd.extend(["--max-unpaired", str(int(adv.get("max_unpaired", 30)))])
    if adv.get("monochrome", True):
        cmd.append("--monochrome")
    else:
        cmd.append("--color")
    if (color_by := adv.get("color_by") or "none") != "none":
        cmd.extend(["--color-by", color_by])


def hash_draw_advanced(digest, adv: Dict[str, Any]) -> None:
    """Fold draw Advanced fields into a content hash digest."""
    digest.update(b"adv-draw\0")
    digest.update((adv.get("force_template") or "").encode("utf-8"))
    digest.update(b"\0")
    digest.update(b"1" if adv.get("constraint") else b"0")
    digest.update(b"1" if adv.get("skip_ribovore_filters") else b"0")
    digest.update((adv.get("tf_engine") or "auto").encode("utf-8"))
    digest.update(b"\0")


def hash_structure_advanced(digest, adv: Dict[str, Any]) -> None:
    """Fold pdb/compare Advanced fields into a content hash digest."""
    digest.update(b"adv-struct\0")
    digest.update(b"1" if adv.get("pseudoknots", True) else b"0")
    digest.update(b"1" if adv.get("rnapuzzler") else b"0")
    digest.update(b"\0")


def hash_align_advanced(digest, stitch: bool, adv: Dict[str, Any]) -> None:
    """Fold align Advanced fields into a content hash digest."""
    digest.update(b"adv-align\0")
    digest.update(b"1" if stitch else b"0")
    digest.update(b"1" if adv.get("include_novel") else b"0")
    digest.update(b"1" if adv.get("all_structures") else b"0")
    digest.update(b"1" if adv.get("auto_repair") else b"0")
    digest.update(b"1" if adv.get("monochrome", True) else b"0")
    digest.update(str(int(adv.get("max_unpaired", 30))).encode("utf-8"))
    digest.update(b"\0")
    digest.update((adv.get("color_by") or "none").encode("utf-8"))
    digest.update(b"\0")
