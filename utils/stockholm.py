# pylint: disable=line-too-long,too-many-lines
"""
Parse Stockholm alignments and extract named secondary structure regions.

This module provides functionality to:
1. Parse Stockholm format multiple sequence alignments
2. Extract named regions from #=GC annotations:
   - New format: #=GC structureID (structure names) + #=GC regionID (parent regions)
   - Legacy format: #=GC knownSS_names / #=GC SS_names
3. Compute Infernal-style RF consensus sequences
4. Validate secondary structures
5. Generate template-free visualizations for each region

Copyright [2009-present] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
     http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import logging
import re
import shutil
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

try:
    from rich import print as rprint
except ImportError:
    rprint = print

from .svg import soften_long_basepair_lines

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class NamedRegion:  # pylint: disable=too-many-instance-attributes
    """A named secondary structure region from a Stockholm alignment."""

    name: str
    start: int  # 0-based column index (inclusive)
    end: int  # 0-based column index (exclusive)
    structure: str  # Secondary structure string (gap-stripped)
    consensus: str  # RF-style consensus sequence
    original_structure: str  # Original structure with gaps
    alignment_start: int  # Position in alignment (for ordering)
    region: Optional[str] = None  # Parent region from #=GC regionID
    is_unnamed: bool = False  # True for auto-generated unnamed regions
    parent_name: Optional[str] = None  # Set on sub-panels after splitting


@dataclass
class StockholmAlignment:  # pylint: disable=too-many-instance-attributes
    """Parsed Stockholm alignment data."""

    sequences: dict[str, str]  # seq_id -> aligned sequence
    ss_cons: str  # Consensus secondary structure
    known_ss_names: str  # Named regions annotation (legacy)
    novel_ss_names: Optional[str] = None  # Optional novel regions (legacy)
    structure_id: Optional[str] = None  # #=GC structureID (new format)
    region_id: Optional[str] = None  # #=GC regionID (new format)
    gc_annotations: dict = None  # Other #=GC annotations
    family_id: str = ""  # #=GF ID (e.g. "SAM")
    family_accession: str = ""  # #=GF AC (e.g. "RF00162")
    rf: str = ""  # #=GC RF annotation

    def __post_init__(self):
        if self.gc_annotations is None:
            self.gc_annotations = {}


# IUPAC ambiguity codes for nucleotides
IUPAC_CODES = {
    frozenset("A"): "A",
    frozenset("C"): "C",
    frozenset("G"): "G",
    frozenset("U"): "U",
    frozenset("T"): "U",  # Treat T as U
    frozenset("AG"): "R",  # puRine
    frozenset("CU"): "Y",  # pYrimidine
    frozenset("CT"): "Y",
    frozenset("AC"): "M",  # aMino
    frozenset("GU"): "K",  # Keto
    frozenset("GT"): "K",
    frozenset("AU"): "W",  # Weak
    frozenset("AT"): "W",
    frozenset("CG"): "S",  # Strong
    frozenset("CGU"): "B",  # not A
    frozenset("CGT"): "B",
    frozenset("AGU"): "D",  # not C
    frozenset("AGT"): "D",
    frozenset("ACU"): "H",  # not G
    frozenset("ACT"): "H",
    frozenset("ACG"): "V",  # not U
    frozenset("ACGU"): "N",
    frozenset("ACGT"): "N",
}


# pylint: disable-next=too-many-branches,too-many-statements
def parse_stockholm(filepath: Path) -> StockholmAlignment:
    """
    Parse a Stockholm format alignment file.

    Args:
        filepath: Path to the Stockholm file

    Returns:
        StockholmAlignment object with parsed data
    """
    sequences = {}
    ss_cons = ""
    known_ss_names = ""
    novel_ss_names = ""
    structure_id = ""
    region_id = ""
    gc_annotations = {}
    family_id = ""
    family_accession = ""
    rf = ""

    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n\r")

            # Skip comments and empty lines
            if not line or line.startswith("# "):
                continue

            # End of alignment
            if line.startswith("//"):
                break

            # Header
            if line.startswith("# STOCKHOLM"):
                continue

            # GF (per-file) annotations
            if line.startswith("#=GF"):
                parts = line.split(None, 2)
                if len(parts) >= 3:
                    tag = parts[1]
                    if tag == "ID":
                        family_id = parts[2].strip()
                    elif tag == "AC":
                        family_accession = parts[2].strip()
                continue

            # GC (per-column) annotations
            if line.startswith("#=GC"):
                parts = line.split(None, 2)
                if len(parts) >= 3:
                    annotation_name = parts[1]
                    annotation_value = parts[2]

                    if annotation_name == "SS_cons":
                        ss_cons += annotation_value
                    elif annotation_name == "structureID":
                        structure_id += annotation_value
                    elif annotation_name == "regionID":
                        region_id += annotation_value
                    elif annotation_name in ("knownSS_names", "SS_names"):
                        known_ss_names += annotation_value
                    elif annotation_name == "novelSS_names":
                        novel_ss_names += annotation_value
                    elif annotation_name == "RF":
                        rf += annotation_value
                    else:
                        gc_annotations[annotation_name] = (
                            gc_annotations.get(annotation_name, "") + annotation_value
                        )
                continue

            # GR (per-residue) annotations - skip
            if line.startswith("#=GR"):
                continue

            # GS (per-sequence) annotations - skip
            if line.startswith("#=GS"):
                continue

            # Sequence lines: "seq_id    sequence"
            if line and not line.startswith("#"):
                parts = line.split(None, 1)
                if len(parts) == 2:
                    seq_id, seq_data = parts
                    if seq_id in sequences:
                        sequences[seq_id] += seq_data
                    else:
                        sequences[seq_id] = seq_data

    return StockholmAlignment(
        sequences=sequences,
        ss_cons=ss_cons,
        known_ss_names=known_ss_names,
        novel_ss_names=novel_ss_names if novel_ss_names else None,
        structure_id=structure_id if structure_id else None,
        region_id=region_id if region_id else None,
        gc_annotations=gc_annotations,
        family_id=family_id,
        family_accession=family_accession,
        rf=rf,
    )


def wuss_to_dotbracket(ss_cons: str) -> str:
    """Convert WUSS notation to dot-bracket with letter-pair pseudoknots.

    In WUSS, ``<>``, ``[]``, ``{}`` are nesting-depth indicators for
    regular (non-crossing) base pairs — they all map to ``()``.  Only
    letter pairs (``Aa``, ``Bb``, …) denote pseudoknots and are kept
    as-is.  Unpaired WUSS characters (``,``, ``-``, ``_``, ``:``,
    ``~``) become dots.

    Args:
        ss_cons: Secondary structure in WUSS notation.

    Returns:
        Structure string using ``()`` and letter-pair pseudoknots.
    """
    result = []
    for c in ss_cons:
        if c in "(<[{":
            result.append("(")
        elif c in ")>]}":
            result.append(")")
        elif c.isalpha():
            result.append(c)  # letter-pair pseudoknot
        else:
            result.append(".")  # '.', ',', '-', '_', ':', '~'
    return "".join(result)


def extract_named_regions(
    alignment: StockholmAlignment,
    include_novel: bool = False,
    include_all: bool = False,
    fallback_name: str = "",
) -> list[NamedRegion]:
    """
    Extract named regions from the Stockholm alignment annotations.

    Supports two annotation formats (new format takes priority):

    **New format** (preferred):
    - ``#=GC structureID`` — pipe-delimited structure names (e.g. SLI, IRES)
    - ``#=GC regionID``   — pipe-delimited parent region names (e.g. 5'UTR)

    **Legacy format** (fallback):
    - ``#=GC knownSS_names`` (or ``SS_names``) — pipe-delimited region names
    - ``#=GC novelSS_names`` — optional additional regions

    Regions are delimited by pipe (|) symbols. The name is extracted
    from between the pipes, and the corresponding SS_cons columns are
    used for the secondary structure.

    Args:
        alignment: Parsed Stockholm alignment
        include_novel: If True, also extract regions from novelSS_names (legacy)
        include_all: If True, also include unnamed segments that contain
            secondary structure (base pairs).  These are assigned synthetic
            names based on their neighbouring named regions.

    Returns:
        List of NamedRegion objects, sorted by alignment position
    """
    regions = []

    # New format: prefer structureID over legacy knownSS_names
    if alignment.structure_id:
        regions.extend(
            _extract_regions_from_annotation(
                alignment.structure_id,
                alignment.ss_cons,
                alignment.sequences,
                include_all=include_all,
            )
        )
        # Assign parent regions from regionID if available
        if alignment.region_id:
            _assign_parent_regions(regions, alignment.region_id)
    else:
        # Legacy format: use knownSS_names / SS_names
        if alignment.known_ss_names:
            regions.extend(
                _extract_regions_from_annotation(
                    alignment.known_ss_names,
                    alignment.ss_cons,
                    alignment.sequences,
                    include_all=include_all,
                )
            )

        # Optionally process novel SS names (legacy only)
        if include_novel and alignment.novel_ss_names:
            regions.extend(
                _extract_regions_from_annotation(
                    alignment.novel_ss_names,
                    alignment.ss_cons,
                    alignment.sequences,
                    include_all=include_all,
                )
            )

    # Whole-alignment fallback: no structureID or knownSS_names
    if not regions:
        name = (
            alignment.family_id
            or alignment.family_accession
            or fallback_name
            or "unknown"
        )
        if alignment.rf:
            # Rfam-style: use RF to identify match columns
            match_cols = [i for i, c in enumerate(alignment.rf) if c != "."]
            ss_match = "".join(alignment.ss_cons[i] for i in match_cols)
            structure = wuss_to_dotbracket(ss_match)
            consensus = "".join(alignment.rf[i].upper() for i in match_cols)
        else:
            # Plain alignment: compute consensus and strip gap columns
            consensus_raw = compute_rf_consensus(
                alignment.sequences, 0, len(alignment.ss_cons)
            )
            structure, consensus = remove_gap_columns(alignment.ss_cons, consensus_raw)
        regions.append(
            NamedRegion(
                name=name,
                start=0,
                end=len(alignment.ss_cons),
                structure=structure,
                consensus=consensus,
                original_structure=alignment.ss_cons,
                alignment_start=0,
            )
        )

    # Sort by alignment position
    regions.sort(key=lambda r: r.alignment_start)

    return regions


def _extract_segment_name(segment: str) -> str:
    """Extract a clean name from a pipe-delimited annotation segment.

    In Stockholm annotations, region names are tiled across columns so
    a segment might look like ``...5'UTR...5'UTR...5'UTR...``.  After
    stripping dot fillers we get ``5'UTR5'UTR5'UTR``.

    This function finds the shortest repeating unit, returning ``5'UTR``.
    It also handles near-repeats where a single-character typo breaks the
    pattern (e.g. ``core_proteinfcore_protein``).
    """
    name = segment.replace(".", "").strip()
    if not name:
        return ""
    # Find shortest exact repeating substring
    for length in range(1, len(name) // 2 + 1):
        candidate = name[:length]
        if len(name) % length == 0 and name == candidate * (len(name) // length):
            return candidate
    # Handle near-repeats: split on the first occurrence of the candidate
    # substring appearing later in the string (e.g. "core_proteinfcore_protein")
    for length in range(2, len(name) // 2 + 1):
        candidate = name[:length]
        second = name.find(candidate, 1)
        if second != -1 and second <= length + 1:
            # Found a near-repeat — use the candidate as the name
            return candidate
    return name


_BRACKET_CHARS = set("([{<ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz)>]}")


def _has_base_pairs(structure: str) -> bool:
    """Return True if *structure* contains at least one bracket character."""
    return any(c in _BRACKET_CHARS for c in structure)


def split_region_at_unpaired(  # pylint: disable=too-many-locals,too-many-branches,too-many-statements,too-many-return-statements
    region: NamedRegion,
    max_unpaired: int,
    sequences: dict[str, str],
    context_nt: int = 5,
) -> list[NamedRegion]:
    """Split a region at long unpaired stretches that occur at nesting depth 0.

    Only runs of consecutive unpaired (dot) characters where **all**
    bracket types have depth 0 are eligible for splitting.  This
    guarantees that no base-pair is broken across sub-panels.

    Each split keeps *context_nt* unpaired positions on either side so
    that neighbouring stems retain some single-stranded context in the
    visualisation.

    Args:
        region: The NamedRegion to potentially split.
        max_unpaired: Minimum dot-run length to trigger a split.
            Use ``0`` to disable splitting.
        sequences: Alignment sequences (needed for rebuilding the column
            map between gap-stripped and alignment coordinates).
        context_nt: Number of context dots to keep on each side of a cut.

    Returns:
        List of NamedRegion objects.  Returns ``[region]`` unchanged when
        no eligible split points are found.
    """
    if max_unpaired <= 0:
        return [region]

    structure = region.structure  # gap-stripped
    num = len(structure)
    if num == 0:
        return [region]

    # ── Bracket types present in this structure ───────────────────────
    bracket_pairs: list[tuple[str, str]] = [
        ("(", ")"),
        ("[", "]"),
        ("{", "}"),
        ("<", ">"),
    ]
    for letter in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
        if letter in structure or letter.lower() in structure:
            bracket_pairs.append((letter, letter.lower()))

    openers = {o for o, _ in bracket_pairs}
    closers = {c for _, c in bracket_pairs}

    # ── Compute depth-before for every position ───────────────────────
    depth_before = [0] * num
    depth = 0
    for i, char in enumerate(structure):
        depth_before[i] = depth
        if char in openers:
            depth += 1
        elif char in closers:
            depth -= 1

    # ── Find runs of unpaired characters at depth 0 ───────────────────
    is_unpaired_d0 = [
        depth_before[i] == 0 and char not in openers and char not in closers
        for i, char in enumerate(structure)
    ]

    runs: list[tuple[int, int]] = []  # (start_inclusive, end_exclusive)
    i = 0
    while i < num:
        if is_unpaired_d0[i]:
            j = i
            while j < num and is_unpaired_d0[j]:
                j += 1
            if j - i >= max_unpaired:
                runs.append((i, j))
            i = j
        else:
            i += 1

    if not runs:
        return [region]

    # ── Build column map: gap-stripped position → absolute alignment col
    original_ss = region.original_structure
    original_consensus = compute_rf_consensus(sequences, region.start, region.end)
    col_map: list[int] = []
    for col_offset, (ss_char, seq_char) in enumerate(
        zip(original_ss, original_consensus)
    ):
        is_gap_seq = seq_char in "-."
        is_gap_ss = ss_char in ".-"
        if not is_gap_seq or not is_gap_ss:
            col_map.append(region.start + col_offset)

    if len(col_map) != num:
        logger.warning(
            "Column map length %d != structure length %d for %s; skipping split",
            len(col_map),
            num,
            region.name,
        )
        return [region]

    # ── Compute cut intervals ─────────────────────────────────────────
    cuts: list[tuple[int, int]] = []
    for run_start, run_end in runs:
        cut_start = run_start + context_nt
        cut_end = run_end - context_nt
        if cut_start < cut_end:
            cuts.append((cut_start, cut_end))

    if not cuts:
        return [region]

    # ── Derive segment boundaries ─────────────────────────────────────
    segments: list[tuple[int, int]] = []
    prev = 0
    for cut_start, cut_end in cuts:
        segments.append((prev, cut_start))
        prev = cut_end
    segments.append((prev, num))

    # ── Build sub-regions ─────────────────────────────────────────────
    sub_regions: list[NamedRegion] = []
    part_num = 0
    for seg_start, seg_end in segments:
        sub_structure = structure[seg_start:seg_end]
        sub_consensus = region.consensus[seg_start:seg_end]

        if not sub_structure or not _has_base_pairs(sub_structure):
            continue

        part_num += 1
        sub_align_start = col_map[seg_start]
        sub_align_end = col_map[seg_end - 1] + 1

        orig_lo = sub_align_start - region.start
        orig_hi = sub_align_end - region.start
        sub_original = original_ss[orig_lo:orig_hi]

        sub_name = f"{region.name}_part{part_num}" if len(segments) > 1 else region.name

        sub_regions.append(
            NamedRegion(
                name=sub_name,
                start=sub_align_start,
                end=sub_align_end,
                structure=sub_structure,
                consensus=sub_consensus,
                original_structure=sub_original,
                alignment_start=sub_align_start,
                region=region.region,
                is_unnamed=region.is_unnamed,
                parent_name=region.name,
            )
        )

    if len(sub_regions) <= 1:
        return [region]

    return sub_regions


def _extract_regions_from_annotation(  # pylint: disable=too-many-locals,too-many-branches
    names_annotation: str,
    ss_cons: str,
    sequences: dict[str, str],
    include_all: bool = False,
) -> list[NamedRegion]:
    """
    Extract regions from a names annotation line.

    Pipe (``|``) characters in the annotation mark region boundaries.
    Because pipes occupy alignment columns, the corresponding SS_cons
    characters at those positions would normally be lost.  To avoid
    dropping structurally important brackets, this function inspects
    the SS_cons character at each boundary:

    * An opening bracket (``(``, ``<``, ``[``, ``{``) at the **left**
      pipe is prepended to the region (it is the region's outermost
      opening bracket).
    * A closing bracket (``)``, ``>``, ``]``, ``}``) at the **right**
      pipe is appended to the region (it is the region's outermost
      closing bracket).

    This ensures that each pipe-position character belongs to exactly
    one region and that balanced structures remain balanced after
    extraction.

    When *include_all* is ``True``, unnamed segments whose SS_cons
    slice contains base-pair brackets are kept as well.  They receive
    a synthetic name derived from their neighbouring named regions
    (e.g. ``_between_SL833_SL1412``) and are marked with
    ``is_unnamed=True``.

    Args:
        names_annotation: The knownSS_names or novelSS_names string
        ss_cons: The SS_cons string
        sequences: Dictionary of aligned sequences
        include_all: If True, keep unnamed segments that have structure

    Returns:
        List of NamedRegion objects
    """
    # ------------------------------------------------------------------
    # Pass 1: extract every segment (named *and* unnamed-but-structured)
    # ------------------------------------------------------------------
    raw_segments: list[tuple[str, bool, int, int]] = []  # (name, unnamed, start, end)

    pipe_positions = [i for i, c in enumerate(names_annotation) if c == "|"]
    boundaries = [-1] + pipe_positions + [len(names_annotation)]

    open_brackets = set("(<[{")
    close_brackets = set(")>]}")

    for seg_idx in range(len(boundaries) - 1):
        left_pipe = boundaries[seg_idx]
        right_pipe = boundaries[seg_idx + 1]

        seg_start = left_pipe + 1
        seg_end = right_pipe
        segment_text = names_annotation[seg_start:seg_end]

        name = _extract_segment_name(segment_text)

        # Determine actual SS_cons column range
        actual_start = seg_start
        actual_end = seg_end

        if 0 <= left_pipe < len(ss_cons):
            if ss_cons[left_pipe] in open_brackets:
                actual_start = left_pipe

        if right_pipe < len(ss_cons):
            if ss_cons[right_pipe] in close_brackets:
                actual_end = right_pipe + 1

        if name:
            raw_segments.append((name, False, actual_start, actual_end))
        elif include_all:
            # Keep the segment only if SS_cons has structure here
            region_ss = ss_cons[actual_start:actual_end]
            if _has_base_pairs(region_ss):
                raw_segments.append(("", True, actual_start, actual_end))

    # ------------------------------------------------------------------
    # Pass 2: assign synthetic names to unnamed segments
    # ------------------------------------------------------------------
    if include_all:
        for idx, (name, is_unnamed, start, end) in enumerate(raw_segments):
            if not is_unnamed:
                continue
            # Find previous named region
            prev_name = "start"
            for j in range(idx - 1, -1, -1):
                if not raw_segments[j][1]:  # not unnamed
                    prev_name = raw_segments[j][0]
                    break
            # Find next named region
            next_name = "end"
            for j in range(idx + 1, len(raw_segments)):
                if not raw_segments[j][1]:
                    next_name = raw_segments[j][0]
                    break
            synthetic = f"_between_{prev_name}_{next_name}"
            raw_segments[idx] = (synthetic, True, start, end)

    # ------------------------------------------------------------------
    # Pass 3: build NamedRegion objects
    # ------------------------------------------------------------------
    regions: list[NamedRegion] = []
    for name, is_unnamed, actual_start, actual_end in raw_segments:
        region_structure = ss_cons[actual_start:actual_end]
        consensus = compute_rf_consensus(sequences, actual_start, actual_end)

        structure_no_gaps, consensus_no_gaps = remove_gap_columns(
            region_structure, consensus
        )

        regions.append(
            NamedRegion(
                name=name,
                start=actual_start,
                end=actual_end,
                structure=structure_no_gaps,
                consensus=consensus_no_gaps,
                original_structure=region_structure,
                alignment_start=actual_start,
                is_unnamed=is_unnamed,
            )
        )

    return regions


def _assign_parent_regions(
    regions: list[NamedRegion], region_id_annotation: str
) -> None:
    """
    Assign parent region names to structure regions using #=GC regionID.

    The regionID annotation uses the same pipe-delimited format as structureID.
    For each structure in ``regions``, we find the regionID segment whose
    column range overlaps with the structure's position in the alignment.

    Args:
        regions: List of NamedRegion objects to update in-place
        region_id_annotation: The #=GC regionID string
    """
    # Parse regionID into (name, start, end) spans
    region_spans = []
    pipe_positions = [i for i, c in enumerate(region_id_annotation) if c == "|"]
    boundaries = [-1] + pipe_positions + [len(region_id_annotation)]

    for seg_idx in range(len(boundaries) - 1):
        left_pipe = boundaries[seg_idx]
        right_pipe = boundaries[seg_idx + 1]
        seg_start = left_pipe + 1
        seg_end = right_pipe
        segment_text = region_id_annotation[seg_start:seg_end]
        name = _extract_segment_name(segment_text)
        if name:
            region_spans.append((name, seg_start, seg_end))

    # For each structure region, find the parent region that overlaps
    for region in regions:
        midpoint = (region.start + region.end) // 2
        for rname, rstart, rend in region_spans:
            if rstart <= midpoint < rend:
                region.region = rname
                break


def compute_rf_consensus(sequences: dict[str, str], start: int, end: int) -> str:
    """
    Compute an Infernal-style RF consensus sequence for alignment columns.

    Uses IUPAC ambiguity codes based on nucleotide frequencies:
    - Single nucleotide if >50% of sequences have that nucleotide
    - IUPAC ambiguity code otherwise
    - Lowercase for poorly conserved positions (<50% non-gap)

    Args:
        sequences: Dictionary of aligned sequences
        start: Start column (0-based, inclusive)
        end: End column (0-based, exclusive)

    Returns:
        Consensus sequence string
    """
    if not sequences:
        return ""

    consensus = []
    seq_list = list(sequences.values())
    num_seqs = len(seq_list)

    for col in range(start, end):
        # Collect nucleotides at this column
        nucleotides = []
        for seq in seq_list:
            if col < len(seq):
                nt = seq[col].upper()
                if nt in "ACGUT":
                    # Normalize T to U
                    if nt == "T":
                        nt = "U"
                    nucleotides.append(nt)

        if not nucleotides:
            # All gaps
            consensus.append("-")
            continue

        # Calculate frequencies
        total = len(nucleotides)
        gap_fraction = 1 - (total / num_seqs)

        counts = {}
        for nt in nucleotides:
            counts[nt] = counts.get(nt, 0) + 1

        # Find the most common nucleotide(s)
        max_count = max(counts.values())
        max_freq = max_count / total

        # Determine consensus character
        if max_freq > 0.5:
            # Single nucleotide dominates
            consensus_nt = [nt for nt, c in counts.items() if c == max_count][0]
        else:
            # Use IUPAC code for observed nucleotides
            observed = frozenset(counts.keys())
            consensus_nt = IUPAC_CODES.get(observed, "N")

        # Use lowercase for poorly conserved positions (>50% gaps)
        if gap_fraction > 0.5:
            consensus_nt = consensus_nt.lower()

        consensus.append(consensus_nt)

    return "".join(consensus)


def remove_gap_columns(structure: str, sequence: str) -> tuple[str, str]:
    """
    Remove columns that are gaps in both structure and sequence.

    Specifically, removes columns where:
    - The sequence is a gap character (-, .)
    - The structure is a gap character (., -)

    Args:
        structure: Secondary structure string
        sequence: Consensus sequence string

    Returns:
        Tuple of (filtered_structure, filtered_sequence)
    """
    if len(structure) != len(sequence):
        # Pad shorter one
        max_len = max(len(structure), len(sequence))
        structure = structure.ljust(max_len, ".")
        sequence = sequence.ljust(max_len, "-")

    filtered_structure = []
    filtered_sequence = []

    for ss, seq in zip(structure, sequence):
        # Keep the column if:
        # - The sequence is a nucleotide (not gap)
        # - OR the structure has a bracket (even if sequence is gap)
        is_gap_seq = seq in "-."
        is_gap_ss = ss in ".-"

        if not is_gap_seq or not is_gap_ss:
            filtered_structure.append(ss)
            filtered_sequence.append(seq if not is_gap_seq else "N")

    return "".join(filtered_structure), "".join(filtered_sequence)


def repair_secondary_structure(structure: str) -> tuple[str, list[str]]:
    """
    Attempt to repair a secondary structure with unbalanced brackets.

    Strategy:
    - Unmatched ')' with no corresponding '(' → convert to '.'
    - Unmatched '(' with no corresponding ')' → convert to '.'
    - Same logic applies to other bracket types and pseudoknot notation

    Args:
        structure: Secondary structure in dot-bracket notation

    Returns:
        Tuple of (repaired_structure, list_of_repairs_made)
    """
    repairs = []
    struct_list = list(structure)

    # Define bracket pairs to check
    bracket_pairs = [
        ("(", ")"),
        ("[", "]"),
        ("{", "}"),
        ("<", ">"),
    ]

    # Add letter pairs for pseudoknot notation
    for letter in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
        if letter in structure or letter.lower() in structure:
            bracket_pairs.append((letter, letter.lower()))

    for open_b, close_b in bracket_pairs:
        # First pass: fix unmatched close brackets
        stack = []
        for i, c in enumerate(struct_list):
            if c == open_b:
                stack.append(i)
            elif c == close_b:
                if not stack:
                    # Unmatched close bracket - convert to dot
                    struct_list[i] = "."
                    repairs.append(
                        f"Position {i+1}: '{close_b}' → '.' (no matching '{open_b}')"
                    )
                else:
                    stack.pop()

        # Second pass: fix unmatched open brackets (remaining in stack)
        for pos in stack:
            struct_list[pos] = "."
            repairs.append(
                f"Position {pos+1}: '{open_b}' → '.' (no matching '{close_b}')"
            )

    return "".join(struct_list), repairs


# pylint: disable-next=too-many-return-statements
def validate_secondary_structure(structure: str) -> tuple[bool, str]:
    """
    Validate a secondary structure string.

    Checks for:
    - Properly nested parentheses (including pseudoknot notation)
    - Valid characters
    - Minimum length

    Args:
        structure: Secondary structure in dot-bracket notation

    Returns:
        Tuple of (is_valid, error_message)
    """
    if not structure:
        return False, "Empty structure"

    if len(structure) < 3:
        return False, f"Structure too short: {len(structure)} characters"

    # Check for valid characters
    valid_chars = set(
        ".()-[]{}<>AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXxYyZz_:,~"
    )
    if invalid := set(structure) - valid_chars:
        return False, f"Invalid characters in structure: {invalid}"

    # Check proper nesting for each bracket type
    bracket_pairs = [
        ("(", ")"),
        ("[", "]"),
        ("{", "}"),
        ("<", ">"),
    ]

    # Also add letter pairs for pseudoknot notation
    for letter in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
        if letter in structure or letter.lower() in structure:
            bracket_pairs.append((letter, letter.lower()))

    for open_b, close_b in bracket_pairs:
        stack = []
        for i, c in enumerate(structure):
            if c == open_b:
                stack.append(i)
            elif c == close_b:
                if not stack:
                    return (
                        False,
                        f"Unmatched '{close_b}' at position {i+1} (no matching '{open_b}')",
                    )
                stack.pop()

        if stack:
            return (
                False,
                f"Unmatched '{open_b}' at positions: {[s+1 for s in stack[:5]]}{'...' if len(stack) > 5 else ''}",
            )

    # Check that there's at least one base pair
    has_pairs = any(
        c in structure
        for c in "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
    )
    if not has_pairs:
        return False, "Structure has no base pairs"

    return True, ""


def write_fasta_with_structure(
    output_path: Path, name: str, sequence: str, structure: str
) -> None:
    """
    Write a FASTA file with secondary structure annotation.

    Format:
        >name
        sequence
        structure

    Args:
        output_path: Path to output file
        name: Sequence name/ID
        sequence: Nucleotide sequence
        structure: Secondary structure in dot-bracket notation
    """
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(f">{name}\n")
        f.write(f"{sequence}\n")
        f.write(f"{structure}\n")


def _normalize_svg_scale(svg_path: Path, target_font_size: float = 9.0) -> None:
    """Normalise font sizes and stroke widths in a Traveler SVG.

    Only CSS values (``font-size``, ``stroke-width``) are modified;
    coordinates and SVG dimensions are left **untouched** so that the
    overall structure layout is preserved.  This keeps large structures
    (e.g. the full IRES domain) visually prominent while ensuring
    consistent text and line weight across panels.

    The function is a no-op when the detected font size is already within
    5 % of the target.

    Args:
        svg_path: Path to the SVG file (modified in place).
        target_font_size: Desired nucleotide font size in SVG user units.
    """
    tree = ET.parse(svg_path)
    root = tree.getroot()

    # ── Detect current nucleotide font size from CSS ──────────────────
    current_font_size = None
    for elem in root.iter():
        tag = elem.tag.split("}")[-1] if "}" in elem.tag else elem.tag
        if tag == "style" and elem.text:
            match = re.search(r"text\s*\{[^}]*font-size:\s*([\d.]+)", elem.text)
            if match:
                current_font_size = float(match.group(1))
                break

    if current_font_size is None or current_font_size <= 0:
        return  # Cannot determine font size – leave SVG unchanged

    scale = target_font_size / current_font_size
    if 0.95 <= scale <= 1.05:
        return  # Already close enough

    # ── Scale only font-size and stroke-width in CSS <style> blocks ───
    # Coordinates and SVG width/height are intentionally left unchanged
    # so the panel retains its natural size in the stitched output.
    _numeric_css = re.compile(r"(font-size:|stroke-width:)\s*([\d.]+)")

    for elem in root.iter():
        tag = elem.tag.split("}")[-1] if "}" in elem.tag else elem.tag
        if tag == "style" and elem.text:
            elem.text = _numeric_css.sub(
                lambda m: f"{m.group(1)} {round(float(m.group(2)) * scale, 6)}",
                elem.text,
            )

    tree.write(svg_path, xml_declaration=False)


def _enrich_covariation(
    alignment: StockholmAlignment,
    sub_region: NamedRegion,
    r2r_folder: Path,
    src_svg: Path,
    quiet: bool = False,
) -> Path:
    """Enrich a region SVG with covariation colours if data is available.

    Looks for the ``.colored.json`` produced by Traveler, generates a
    covariation TSV from the R-scape ``#=GC`` annotations, then calls
    ``enrich_json.py`` + ``json2svg.py`` to produce an ``.enriched.svg``.

    Returns the path to the enriched SVG on success, or *src_svg*
    unchanged on failure.
    """
    # pylint: disable=import-outside-toplevel
    from . import covariation
    from .runner import runner

    json_candidates = list(r2r_folder.glob(f"*{sub_region.name}*.colored.json")) + list(
        r2r_folder.glob("*.colored.json")
    )
    if not json_candidates:
        return src_svg

    colored_json = json_candidates[0]
    tsv_path = r2r_folder / "covariation.tsv"

    wrote_tsv = covariation.generate_covariation_tsv(
        alignment=alignment,
        region_start=sub_region.start,
        region_end=sub_region.end,
        original_structure=sub_region.original_structure,
        consensus=sub_region.consensus,
        structure=sub_region.structure,
        output_path=tsv_path,
    )
    if not wrote_tsv:
        return src_svg

    enriched_json = r2r_folder / f"{sub_region.name}.enriched.json"
    enriched_svg = r2r_folder / f"{sub_region.name}.enriched.svg"

    cmd = (
        f"python3 /rna/traveler/utils/enrich_json.py "
        f"--input-json {colored_json} "
        f"--input-data {tsv_path} "
        f"--output {enriched_json}"
    )
    result = runner.run(cmd)
    if result != 0:
        return src_svg

    cmd = (
        f"python3 /rna/traveler/utils/json2svg.py "
        f"-p /rna/r2dt/utils/covariation_colorscheme.json "
        f"-i {enriched_json} -o {enriched_svg}"
    )
    runner.run(cmd)

    if enriched_svg.exists():
        if not quiet:
            rprint("    [green]✓[/green] Added covariation colouring")
        return enriched_svg
    return src_svg


# pylint: disable-next=too-many-locals,too-many-branches,too-many-statements,too-many-arguments,too-many-positional-arguments
def process_stockholm_alignment(
    stockholm_path: Path,
    output_folder: Path,
    include_novel: bool = False,
    quiet: bool = False,
    auto_repair: bool = False,
    include_all: bool = False,
    max_unpaired: int = 0,
    fallback_name: str = "",
) -> list[dict]:
    """
    Process a Stockholm alignment and generate template-free visualizations.

    This function:
    1. Parses the Stockholm alignment
    2. Extracts named regions from #=GC structureID (preferred) or
       #=GC knownSS_names (legacy fallback)
    3. Assigns parent regions from #=GC regionID when available
    4. Computes RF consensus for each region
    5. Validates secondary structures (with optional auto-repair)
    6. Runs template-free R2DT for valid regions
    7. Returns information for stitching

    Args:
        stockholm_path: Path to the Stockholm alignment file
        output_folder: Output folder for results
        include_novel: If True, also process novelSS_names regions
        quiet: If True, suppress progress output
        auto_repair: If True, attempt to repair unbalanced bracket structures

    Returns:
        List of dicts with region info and SVG paths for stitching
    """
    # pylint: disable=import-outside-toplevel
    # Import here to avoid circular imports
    from . import covariation, r2r, rfam, rnapuzzler
    from .scale_template import scale_coordinates

    if not quiet:
        rprint(f"[blue]Processing Stockholm alignment:[/blue] {stockholm_path}")

    # Parse the alignment
    alignment = parse_stockholm(stockholm_path)

    if not quiet:
        rprint(f"  Found {len(alignment.sequences)} sequences")
        rprint(f"  SS_cons length: {len(alignment.ss_cons)}")
        if alignment.structure_id:
            rprint("  Annotation format: structureID + regionID (new)")
        elif alignment.known_ss_names:
            rprint("  Annotation format: knownSS_names (legacy)")
        else:
            rprint("  Annotation format: simple alignment (whole-alignment fallback)")
            if alignment.rf:
                rprint("  Using #=GC RF for consensus")

    # Detect R-scape covariation annotations
    has_cov = covariation.has_covariation(alignment)
    if has_cov and not quiet:
        rprint("  [green]R-scape covariation annotations detected[/green]")

    # Extract named regions
    regions = extract_named_regions(
        alignment,
        include_novel=include_novel,
        include_all=include_all,
        fallback_name=fallback_name,
    )

    if not quiet:
        named_count = sum(1 for r in regions if not r.is_unnamed)
        unnamed_count = sum(1 for r in regions if r.is_unnamed)
        rprint(f"  Found {named_count} named regions")
        if unnamed_count:
            rprint(f"  Found {unnamed_count} unnamed structured regions")

    # Create output directories
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    results_folder = output_folder / "results"
    results_folder.mkdir(exist_ok=True)

    svg_folder = results_folder / "svg"
    svg_folder.mkdir(exist_ok=True)

    fasta_folder = results_folder / "fasta"
    fasta_folder.mkdir(exist_ok=True)

    regions_folder = output_folder / "regions"
    regions_folder.mkdir(exist_ok=True)

    # Process each region
    processed_regions = []
    skipped_regions = []
    repaired_regions = []

    for region in regions:  # pylint: disable=too-many-nested-blocks
        if not quiet:
            rprint(f"  Processing region: [green]{region.name}[/green]")

        # Validate structure
        is_valid, error_msg = validate_secondary_structure(region.structure)

        if not is_valid:
            if auto_repair:
                # Attempt to repair the structure
                repaired_structure, repairs = repair_secondary_structure(
                    region.structure
                )
                is_valid_now, new_error = validate_secondary_structure(
                    repaired_structure
                )

                if is_valid_now:
                    if not quiet:
                        rprint(
                            f"    [cyan]Auto-repaired:[/cyan] {len(repairs)} fix(es)"
                        )
                        for repair in repairs[:3]:
                            rprint(f"      - {repair}")
                        if len(repairs) > 3:
                            rprint(f"      ... and {len(repairs) - 3} more")

                    region = NamedRegion(
                        name=region.name,
                        start=region.start,
                        end=region.end,
                        structure=repaired_structure,
                        consensus=region.consensus,
                        original_structure=region.original_structure,
                        alignment_start=region.alignment_start,
                        region=region.region,
                        is_unnamed=region.is_unnamed,
                        parent_name=region.parent_name,
                    )
                    repaired_regions.append((region.name, repairs))
                else:
                    skipped_regions.append((region.name, f"Repair failed: {new_error}"))
                    if not quiet:
                        rprint(
                            f"    [yellow]Skipping:[/yellow] Repair failed - {new_error}"
                        )
                    logger.warning(
                        "Skipping region %s: Repair failed - %s", region.name, new_error
                    )
                    continue
            else:
                skipped_regions.append((region.name, error_msg))
                if not quiet:
                    rprint(f"    [yellow]Skipping:[/yellow] {error_msg}")
                logger.warning("Skipping region %s: %s", region.name, error_msg)
                continue

        # Check sequence/structure length match
        if len(region.consensus) != len(region.structure):
            error_msg = (
                f"Length mismatch: sequence={len(region.consensus)}, "
                f"structure={len(region.structure)}"
            )
            skipped_regions.append((region.name, error_msg))
            if not quiet:
                rprint(f"    [yellow]Skipping:[/yellow] {error_msg}")
            logger.warning("Skipping region %s: %s", region.name, error_msg)
            continue

        # Split at long unpaired stretches (depth-0 only)
        sub_regions = split_region_at_unpaired(
            region, max_unpaired, alignment.sequences
        )
        if not quiet and len(sub_regions) > 1:
            rprint(
                f"    Split into {len(sub_regions)} sub-panels "
                f"(max_unpaired={max_unpaired})"
            )

        for sub_region in sub_regions:
            # Safety: re-validate sub-region structure
            sub_valid, sub_err = validate_secondary_structure(sub_region.structure)
            if not sub_valid:
                skipped_regions.append((sub_region.name, sub_err))
                if not quiet:
                    rprint(
                        f"    [yellow]Skipping sub-panel "
                        f"{sub_region.name}:[/yellow] {sub_err}"
                    )
                continue

            if len(sub_region.consensus) != len(sub_region.structure):
                err = (
                    f"Length mismatch: sequence={len(sub_region.consensus)}, "
                    f"structure={len(sub_region.structure)}"
                )
                skipped_regions.append((sub_region.name, err))
                continue

            # Create safe filename
            safe_name = re.sub(r"[^\w\-]", "_", sub_region.name)

            # Write FASTA file
            region_folder = regions_folder / safe_name
            region_folder.mkdir(exist_ok=True)

            fasta_path = region_folder / f"{safe_name}.fasta"
            write_fasta_with_structure(
                fasta_path,
                sub_region.name,
                sub_region.consensus,
                sub_region.structure,
            )

            # Copy to results fasta folder
            shutil.copy(fasta_path, fasta_folder / f"{safe_name}.fasta")

            # Run template-free visualization (RNApuzzler first, R2R fallback)
            r2r_folder = region_folder / "r2r"
            r2r_folder.mkdir(exist_ok=True)

            try:
                # Try RNApuzzler first (overlap-free)
                try:
                    rnapuzzler.run_puzzler_pipeline(
                        str(fasta_path),
                        r2r_folder,
                        sub_region.name,
                        sub_region.consensus,
                        sub_region.structure,
                    )
                except Exception:  # pylint: disable=broad-except
                    # Fall back to R2R
                    r2r.generate_r2r_input_file(
                        sub_region.consensus,
                        sub_region.structure,
                        r2r_folder,
                    )
                    r2r_svg = r2r.run_r2r(r2r_folder)
                    rscape_one_line_svg = rfam.convert_rscape_svg_to_one_line(
                        r2r_svg, r2r_folder
                    )
                    rfam.convert_rscape_svg_to_traveler(rscape_one_line_svg, r2r_folder)
                    scale_coordinates(
                        r2r_folder / "traveler-template.xml",
                        scaling_factor=3,
                    )
                    r2r.run_traveler(fasta_path, r2r_folder, sub_region.name)

                # Find the output SVG
                svg_candidates = list(
                    r2r_folder.glob(f"*{sub_region.name}*.colored.svg")
                ) + list(r2r_folder.glob("*.colored.svg"))

                if svg_candidates:
                    src_svg = svg_candidates[0]
                    base_stem = f"{safe_name}_{sub_region.start}-{sub_region.end}"

                    # Always copy the normal (uncolored) SVG
                    dest_svg = svg_folder / f"{base_stem}.svg"
                    shutil.copy(src_svg, dest_svg)

                    # Covariation enrichment (R-scape annotations)
                    if has_cov:
                        enriched = _enrich_covariation(
                            alignment,
                            sub_region,
                            r2r_folder,
                            src_svg,
                            quiet,
                        )
                        if enriched != src_svg:
                            dest_enriched = svg_folder / f"{base_stem}.covariation.svg"
                            shutil.copy(enriched, dest_enriched)
                            _normalize_svg_scale(dest_enriched, target_font_size=9.0)
                            soften_long_basepair_lines(dest_enriched)

                    # Normalise the SVG so that its nucleotide font size
                    # matches the typical Traveler output (~9 px).  Large
                    # regions (e.g. the full IRES) can be generated at a
                    # much bigger scale; rescaling here keeps all panels
                    # visually consistent when stitched together.
                    _normalize_svg_scale(dest_svg, target_font_size=9.0)
                    soften_long_basepair_lines(dest_svg)

                    processed_regions.append(
                        {
                            "name": sub_region.name,
                            "safe_name": safe_name,
                            "start": sub_region.start,
                            "end": sub_region.end,
                            "svg_path": str(dest_svg),
                            "sequence": sub_region.consensus,
                            "structure": sub_region.structure,
                            "region": sub_region.region,
                            "is_unnamed": sub_region.is_unnamed,
                            "parent_name": sub_region.parent_name,
                        }
                    )

                    if not quiet:
                        rprint(f"    [green]✓[/green] Generated SVG: {dest_svg.name}")
                else:
                    error_msg = "No SVG output generated"
                    skipped_regions.append((sub_region.name, error_msg))
                    if not quiet:
                        rprint(f"    [red]✗[/red] {error_msg}")
                    logger.error("No SVG generated for region %s", sub_region.name)

            except Exception as e:  # pylint: disable=broad-exception-caught
                error_msg = str(e)
                skipped_regions.append((sub_region.name, error_msg))
                if not quiet:
                    rprint(f"    [red]✗[/red] Error: {error_msg}")
                logger.error(
                    "Error processing region %s: %s", sub_region.name, error_msg
                )

    # Write summary
    summary_path = output_folder / "processing_summary.txt"
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write(f"Stockholm alignment: {stockholm_path}\n")
        f.write(f"Total regions found: {len(regions)}\n")
        f.write(f"Successfully processed: {len(processed_regions)}\n")
        f.write(f"Auto-repaired: {len(repaired_regions)}\n")
        f.write(f"Skipped: {len(skipped_regions)}\n\n")

        if processed_regions:
            f.write("Processed regions:\n")
            for r in processed_regions:
                f.write(f"  - {r['name']} (columns {r['start']}-{r['end']})\n")

        if repaired_regions:
            f.write("\nAuto-repaired regions:\n")
            for name, repairs in repaired_regions:
                f.write(f"  - {name}:\n")
                for repair in repairs:
                    f.write(f"      {repair}\n")

        if skipped_regions:
            f.write("\nSkipped regions:\n")
            for name, reason in skipped_regions:
                f.write(f"  - {name}: {reason}\n")

    if not quiet:
        rprint("\n[blue]Summary:[/blue]")
        rprint(f"  Processed: {len(processed_regions)} regions")
        if repaired_regions:
            rprint(f"  Auto-repaired: {len(repaired_regions)} regions")
        rprint(f"  Skipped: {len(skipped_regions)} regions")
        if skipped_regions:
            rprint(f"  See {summary_path} for details")

    return processed_regions
