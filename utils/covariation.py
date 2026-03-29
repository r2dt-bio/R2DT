"""Extract covariation annotations from R-scape Stockholm files.

R-scape annotates covariation in Stockholm files via ``#=GC`` lines:

- ``cov_SS_cons``: ``2`` at paired positions with significant covariation
- ``cov_h_SS_cons``: ``3`` at positions in helices with significant covariation

This module detects those annotations and provides per-basepair covariation
data for rectangle-based visualisation in the SVG output.

Copyright [2009-present] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0
"""

import logging
from pathlib import Path

from .stockholm import StockholmAlignment, compute_rf_consensus, rf_match_columns

logger = logging.getLogger(__name__)


def has_covariation(alignment: StockholmAlignment) -> bool:
    """Return True if *alignment* contains R-scape covariation annotations."""
    return "cov_SS_cons" in (alignment.gc_annotations or {})


def _strip_annotation_for_region(
    annotation_full: str,
    alignment: StockholmAlignment,
    region_start: int,
    region_end: int,
    original_structure: str,
) -> str:
    """Strip a per-column annotation to match a region's gap-free consensus.

    When the alignment has ``#=GC RF``, uses RF match columns (non-dot)
    for filtering — matching :func:`stockholm.extract_named_regions`.
    Otherwise falls back to the standard gap-column removal logic from
    :func:`stockholm.remove_gap_columns`.
    """
    ann_slice = annotation_full[region_start:region_end]

    if alignment.rf:
        rf_slice = alignment.rf[region_start:region_end]
        ss_slice = original_structure
        keep = set(rf_match_columns(rf_slice, ss_slice))
        return "".join(ann for idx, ann in enumerate(ann_slice) if idx in keep)

    raw_consensus = compute_rf_consensus(alignment.sequences, region_start, region_end)

    result = []
    for ann, ss_char, seq_char in zip(ann_slice, original_structure, raw_consensus):
        is_gap_seq = seq_char in "-."
        is_gap_ss = ss_char in ".-"
        if not is_gap_seq or not is_gap_ss:
            result.append(ann)
    return "".join(result)


def _build_partner_map(structure: str) -> dict[int, int]:
    """Build a base-pair partner map from a dot-bracket structure.

    Handles regular brackets ``()``, ``<>``, ``[]``, ``{}`` and
    letter-pair pseudoknots (``Aa``, ``Bb``, …).

    Returns a dict mapping each paired position index to its partner.
    """
    open_chars = "(<[{"
    close_chars = ")>]}"
    close_to_open = dict(zip(close_chars, open_chars))
    stack: dict[str, list[int]] = {c: [] for c in open_chars}
    partner: dict[int, int] = {}

    # Also handle letter-pair pseudoknots
    letter_stacks: dict[str, list[int]] = {}

    for i, char in enumerate(structure):
        if char in open_chars:
            stack[char].append(i)
        elif char in close_chars:
            opener = close_to_open[char]
            if stack[opener]:
                j = stack[opener].pop()
                partner[i] = j
                partner[j] = i
        elif char.isupper() and char.isalpha():
            letter_stacks.setdefault(char, []).append(i)
        elif char.islower() and char.isalpha():
            upper = char.upper()
            if upper in letter_stacks and letter_stacks[upper]:
                j = letter_stacks[upper].pop()
                partner[i] = j
                partner[j] = i
    return partner


_PRIORITY = {"dark_green": 2, "light_green": 1, "none": 0}


# pylint: disable=too-many-arguments,too-many-positional-arguments
def _merge_pk_annotations(base: str, gc: dict[str, str], prefix: str) -> str:
    """Overlay ``_pk_*`` covariation layers onto *base*.

    For each ``{prefix}_pk_1``, ``{prefix}_pk_2``, … annotation,
    non-dot characters replace dots in *base* so pseudoknot covariation
    is included.
    """
    result = list(base)
    for key in sorted(k for k in gc if k.startswith(f"{prefix}_pk_")):
        layer = gc[key]
        for i, c in enumerate(layer):
            if i < len(result) and c != "." and result[i] == ".":
                result[i] = c
    return "".join(result)


# pylint: disable=too-many-arguments,too-many-positional-arguments
def _prepare_annotations(
    alignment: StockholmAlignment,
    region_start: int,
    region_end: int,
    original_structure: str,
    consensus: str,
    structure: str,
) -> tuple[str, str, dict[int, int]] | None:
    """Extract stripped covariation annotations and partner map.

    Returns ``(cov_stripped, hcov_stripped, partner_map)`` or ``None``
    when covariation data is absent or lengths don't match.
    """
    gc = alignment.gc_annotations or {}
    cov_full = gc.get("cov_SS_cons", "")
    hcov_full = gc.get("cov_h_SS_cons", "")

    if not cov_full:
        return None

    # Merge pseudoknot layers (cov_SS_cons_pk_1, …)
    cov_full = _merge_pk_annotations(cov_full, gc, "cov_SS_cons")
    hcov_full = _merge_pk_annotations(
        hcov_full if hcov_full else "." * len(cov_full), gc, "cov_h_SS_cons"
    )

    cov_stripped = _strip_annotation_for_region(
        cov_full, alignment, region_start, region_end, original_structure
    )
    hcov_stripped = _strip_annotation_for_region(
        hcov_full, alignment, region_start, region_end, original_structure
    )

    if len(cov_stripped) != len(consensus):
        logger.warning(
            "cov_SS_cons length %d != consensus length %d; skipping covariation",
            len(cov_stripped),
            len(consensus),
        )
        return None

    if len(hcov_stripped) != len(consensus):
        hcov_stripped = "." * len(consensus)

    partner = _build_partner_map(structure)
    return cov_stripped, hcov_stripped, partner


# pylint: disable=too-many-arguments,too-many-positional-arguments
def _classify_positions(
    alignment: StockholmAlignment,
    region_start: int,
    region_end: int,
    original_structure: str,
    consensus: str,
    structure: str,
) -> tuple[list[str], dict[int, int]] | None:
    """Classify each position as dark_green / light_green / none.

    Returns ``(categories, partner_map)`` or ``None`` on failure.
    """
    result = _prepare_annotations(
        alignment, region_start, region_end, original_structure, consensus, structure
    )
    if result is None:
        return None

    cov_stripped, hcov_stripped, partner = result

    categories: list[str] = []
    for cv, hv in zip(cov_stripped, hcov_stripped):
        if cv == "2":
            categories.append("dark_green")
        elif hv == "3":
            categories.append("light_green")
        else:
            categories.append("none")

    for i, j in partner.items():
        if i < j:
            best = (
                categories[i]
                if _PRIORITY[categories[i]] >= _PRIORITY[categories[j]]
                else categories[j]
            )
            categories[i] = best
            categories[j] = best

    return categories, partner


# pylint: disable=too-many-arguments,too-many-positional-arguments
def get_covariation_pairs(
    alignment: StockholmAlignment,
    region_start: int,
    region_end: int,
    original_structure: str,
    consensus: str,
    structure: str,
) -> list[tuple[int, int, str]]:
    """Return covarying base pairs as ``(i, j, category)`` tuples.

    Positions *i* and *j* are **1-based** residue indices (matching
    Traveler SVG ``Position:`` labels).  A pair that has **both**
    pair-level and helix-level covariation produces two entries
    (one ``light_green``, one ``dark_green``) so that both rectangles
    are drawn — the wider light-green rect behind the narrower
    dark-green one.

    Categories:

    - **dark_green** — significant pair covariation
    - **light_green** — helix-level aggregated covariation

    Returns an empty list when no covariation data is available.
    """
    result = _prepare_annotations(
        alignment, region_start, region_end, original_structure, consensus, structure
    )
    if result is None:
        return []

    cov_stripped, hcov_stripped, partner = result
    pairs: list[tuple[int, int, str]] = []
    for i, j in partner.items():
        if i >= j:
            continue
        has_pair = cov_stripped[i] == "2" or cov_stripped[j] == "2"
        has_helix = hcov_stripped[i] == "3" or hcov_stripped[j] == "3"
        if has_helix:
            pairs.append((i + 1, j + 1, "light_green"))
        if has_pair:
            pairs.append((i + 1, j + 1, "dark_green"))
    return pairs


# pylint: disable=too-many-arguments,too-many-positional-arguments
def generate_covariation_tsv(
    alignment: StockholmAlignment,
    region_start: int,
    region_end: int,
    original_structure: str,
    consensus: str,
    structure: str,
    output_path: Path,
) -> bool:
    """Generate a covariation annotation TSV for one region.

    Categories written to the ``covariation`` column:

    - **dark_green** — ``cov_SS_cons == '2'`` (significant pair covariation)
    - **light_green** — ``cov_h_SS_cons == '3'`` only (helix-level, no
      individual pair significance)
    - **none** — everything else (uncolored)

    Both sides of every base pair receive the same (highest-priority)
    colour so that covariation is always shown symmetrically.

    Args:
        alignment: Parsed Stockholm alignment with ``gc_annotations``.
        region_start: Alignment column start (inclusive, 0-based).
        region_end: Alignment column end (exclusive, 0-based).
        original_structure: Alignment-length ``SS_cons`` for this region.
        consensus: Gap-free consensus sequence for this region.
        structure: Gap-free dot-bracket structure for this region.
        output_path: Where to write the TSV file.

    Returns:
        ``True`` if the TSV was written and contains at least one colored
        position, ``False`` otherwise.
    """
    result = _classify_positions(
        alignment, region_start, region_end, original_structure, consensus, structure
    )
    if result is None:
        return False

    categories, _partner = result

    has_color = False
    lines = ["residue_index\tresidue_name\tcovariation\n"]
    for idx, (nt_char, cat) in enumerate(zip(consensus, categories), 1):
        if cat != "none":
            has_color = True
        lines.append(f"{idx}\t{nt_char}\t{cat}\n")

    output_path = Path(output_path)
    output_path.write_text("".join(lines), encoding="utf-8")
    return has_color
