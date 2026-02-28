"""Extract covariation annotations from R-scape Stockholm files.

R-scape annotates covariation in Stockholm files via ``#=GC`` lines:

- ``cov_SS_cons``: ``2`` at paired positions with significant covariation
- ``cov_h_SS_cons``: ``3`` at positions in helices with significant covariation

This module detects those annotations and generates a TSV file compatible
with the Traveler enrichment pipeline (``enrich_json.py`` + ``json2svg.py``).

Copyright [2009-present] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0
"""

import logging
from pathlib import Path

from .stockholm import StockholmAlignment, compute_rf_consensus

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

    Applies the same column-removal logic as
    :func:`stockholm.remove_gap_columns` so the returned string is 1:1
    with the region's gap-free consensus.
    """
    ann_slice = annotation_full[region_start:region_end]
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

    Returns a dict mapping each paired position index to its partner.
    """
    open_chars = "(<[{"
    close_chars = ")>]}"
    close_to_open = dict(zip(close_chars, open_chars))
    stack: dict[str, list[int]] = {c: [] for c in open_chars}
    partner: dict[int, int] = {}

    for i, char in enumerate(structure):
        if char in open_chars:
            stack[char].append(i)
        elif char in close_chars:
            opener = close_to_open[char]
            if stack[opener]:
                j = stack[opener].pop()
                partner[i] = j
                partner[j] = i
    return partner


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
    cov_full = alignment.gc_annotations.get("cov_SS_cons", "")
    hcov_full = alignment.gc_annotations.get("cov_h_SS_cons", "")

    if not cov_full:
        return False

    # Pad hcov with dots if absent
    if not hcov_full:
        hcov_full = "." * len(cov_full)

    # Strip gaps in sync with the consensus
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
        return False

    if len(hcov_stripped) != len(consensus):
        # Fall back to dots if helix annotation has unexpected length
        hcov_stripped = "." * len(consensus)

    # Initial per-position categorisation (priority: dark > light > none)
    _priority = {"dark_green": 2, "light_green": 1, "none": 0}
    categories: list[str] = []
    for cv, hv in zip(cov_stripped, hcov_stripped):
        if cv == "2":
            categories.append("dark_green")
        elif hv == "3":
            categories.append("light_green")
        else:
            categories.append("none")

    # Symmetrise: both sides of a base pair get the higher-priority colour
    partner = _build_partner_map(structure)
    for i, j in partner.items():
        if i < j:  # process each pair once
            best = (
                categories[i]
                if _priority[categories[i]] >= _priority[categories[j]]
                else categories[j]
            )
            categories[i] = best
            categories[j] = best

    # Build TSV rows
    has_color = False
    lines = ["residue_index\tresidue_name\tcovariation\n"]
    for idx, (nt_char, cat) in enumerate(zip(consensus, categories), 1):
        if cat != "none":
            has_color = True
        lines.append(f"{idx}\t{nt_char}\t{cat}\n")

    output_path = Path(output_path)
    output_path.write_text("".join(lines), encoding="utf-8")
    return has_color
