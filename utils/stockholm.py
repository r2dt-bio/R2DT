# pylint: disable=line-too-long
"""
Parse Stockholm alignments and extract named secondary structure regions.

This module provides functionality to:
1. Parse Stockholm format multiple sequence alignments
2. Extract named regions from #=GC knownSS_names annotation
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
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

try:
    from rich import print as rprint
except ImportError:
    rprint = print

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class NamedRegion:
    """A named secondary structure region from a Stockholm alignment."""

    name: str
    start: int  # 0-based column index (inclusive)
    end: int  # 0-based column index (exclusive)
    structure: str  # Secondary structure string (gap-stripped)
    consensus: str  # RF-style consensus sequence
    original_structure: str  # Original structure with gaps
    alignment_start: int  # Position in alignment (for ordering)


@dataclass
class StockholmAlignment:
    """Parsed Stockholm alignment data."""

    sequences: dict[str, str]  # seq_id -> aligned sequence
    ss_cons: str  # Consensus secondary structure
    known_ss_names: str  # Named regions annotation
    novel_ss_names: Optional[str] = None  # Optional novel regions
    gc_annotations: dict = None  # Other #=GC annotations

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


# pylint: disable-next=too-many-branches
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
    gc_annotations = {}

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

            # GF (per-file) annotations - skip
            if line.startswith("#=GF"):
                continue

            # GC (per-column) annotations
            if line.startswith("#=GC"):
                parts = line.split(None, 2)
                if len(parts) >= 3:
                    annotation_name = parts[1]
                    annotation_value = parts[2]

                    if annotation_name == "SS_cons":
                        ss_cons += annotation_value
                    elif annotation_name in ("knownSS_names", "SS_names"):
                        known_ss_names += annotation_value
                    elif annotation_name == "novelSS_names":
                        novel_ss_names += annotation_value
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
        gc_annotations=gc_annotations,
    )


def extract_named_regions(
    alignment: StockholmAlignment, include_novel: bool = False
) -> list[NamedRegion]:
    """
    Extract named regions from the knownSS_names annotation.

    Regions are delimited by pipe (|) symbols. The name is extracted
    from between the pipes, and the corresponding SS_cons columns are
    used for the secondary structure.

    Args:
        alignment: Parsed Stockholm alignment
        include_novel: If True, also extract regions from novelSS_names

    Returns:
        List of NamedRegion objects, sorted by alignment position
    """
    regions = []

    # Process known SS names
    if alignment.known_ss_names:
        regions.extend(
            _extract_regions_from_annotation(
                alignment.known_ss_names, alignment.ss_cons, alignment.sequences
            )
        )

    # Optionally process novel SS names
    if include_novel and alignment.novel_ss_names:
        regions.extend(
            _extract_regions_from_annotation(
                alignment.novel_ss_names, alignment.ss_cons, alignment.sequences
            )
        )

    # Sort by alignment position
    regions.sort(key=lambda r: r.alignment_start)

    return regions


def _extract_regions_from_annotation(
    names_annotation: str, ss_cons: str, sequences: dict[str, str]
) -> list[NamedRegion]:
    """
    Extract regions from a names annotation line.

    Args:
        names_annotation: The knownSS_names or novelSS_names string
        ss_cons: The SS_cons string
        sequences: Dictionary of aligned sequences

    Returns:
        List of NamedRegion objects
    """
    regions = []

    # Split by pipe delimiter to get segments
    # Pipes act as delimiters, so we split by them and track positions
    segments = names_annotation.split("|")

    # Track the current position in the original string
    current_pos = 0

    for segment in segments:
        segment_len = len(segment)
        segment_end = current_pos + segment_len

        # Extract the name from the segment (strip dots and whitespace)
        name = segment.replace(".", "").strip()

        if name:
            # Get the structure for this segment
            region_structure = ss_cons[current_pos:segment_end]

            # Get aligned sequences for this region to compute consensus
            consensus = compute_rf_consensus(sequences, current_pos, segment_end)

            # Remove gap columns from structure and consensus
            structure_no_gaps, consensus_no_gaps = remove_gap_columns(
                region_structure, consensus
            )

            regions.append(
                NamedRegion(
                    name=name,
                    start=current_pos,
                    end=segment_end,
                    structure=structure_no_gaps,
                    consensus=consensus_no_gaps,
                    original_structure=region_structure,
                    alignment_start=current_pos,
                )
            )

        # Move past this segment and the pipe delimiter
        current_pos = segment_end + 1  # +1 for the pipe character

    return regions


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


# pylint: disable-next=too-many-locals,too-many-branches,too-many-statements
def process_stockholm_alignment(
    stockholm_path: Path,
    output_folder: Path,
    include_novel: bool = False,
    quiet: bool = False,
    auto_repair: bool = False,
) -> list[dict]:
    """
    Process a Stockholm alignment and generate template-free visualizations.

    This function:
    1. Parses the Stockholm alignment
    2. Extracts named regions from knownSS_names
    3. Computes RF consensus for each region
    4. Validates secondary structures (with optional auto-repair)
    5. Runs template-free R2DT for valid regions
    6. Returns information for stitching

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
    from . import r2r, rfam
    from .scale_template import scale_coordinates

    if not quiet:
        rprint(f"[blue]Processing Stockholm alignment:[/blue] {stockholm_path}")

    # Parse the alignment
    alignment = parse_stockholm(stockholm_path)

    if not quiet:
        rprint(f"  Found {len(alignment.sequences)} sequences")
        rprint(f"  SS_cons length: {len(alignment.ss_cons)}")

    # Extract named regions
    regions = extract_named_regions(alignment, include_novel=include_novel)

    if not quiet:
        rprint(f"  Found {len(regions)} named regions")

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

        # Create safe filename
        safe_name = re.sub(r"[^\w\-]", "_", region.name)

        # Write FASTA file
        region_folder = regions_folder / safe_name
        region_folder.mkdir(exist_ok=True)

        fasta_path = region_folder / f"{safe_name}.fasta"
        write_fasta_with_structure(
            fasta_path, region.name, region.consensus, region.structure
        )

        # Copy to results fasta folder
        shutil.copy(fasta_path, fasta_folder / f"{safe_name}.fasta")

        # Run template-free visualization using R2R
        r2r_folder = region_folder / "r2r"
        r2r_folder.mkdir(exist_ok=True)

        try:
            r2r.generate_r2r_input_file(region.consensus, region.structure, r2r_folder)
            r2r_svg = r2r.run_r2r(r2r_folder)
            rscape_one_line_svg = rfam.convert_rscape_svg_to_one_line(
                r2r_svg, r2r_folder
            )
            rfam.convert_rscape_svg_to_traveler(rscape_one_line_svg, r2r_folder)
            scale_coordinates(r2r_folder / "traveler-template.xml", scaling_factor=3)
            r2r.run_traveler(fasta_path, r2r_folder, region.name)

            # Find the output SVG
            svg_candidates = list(
                r2r_folder.glob(f"*{region.name}*.colored.svg")
            ) + list(r2r_folder.glob("*.colored.svg"))

            if svg_candidates:
                src_svg = svg_candidates[0]
                # Create filename with position info for sorting
                dest_svg = svg_folder / f"{safe_name}_{region.start}-{region.end}.svg"
                shutil.copy(src_svg, dest_svg)

                processed_regions.append(
                    {
                        "name": region.name,
                        "safe_name": safe_name,
                        "start": region.start,
                        "end": region.end,
                        "svg_path": str(dest_svg),
                        "sequence": region.consensus,
                        "structure": region.structure,
                    }
                )

                if not quiet:
                    rprint(f"    [green]✓[/green] Generated SVG: {dest_svg.name}")
            else:
                error_msg = "No SVG output generated"
                skipped_regions.append((region.name, error_msg))
                if not quiet:
                    rprint(f"    [red]✗[/red] {error_msg}")
                logger.error("No SVG generated for region %s", region.name)

        except Exception as e:  # pylint: disable=broad-exception-caught
            error_msg = str(e)
            skipped_regions.append((region.name, error_msg))
            if not quiet:
                rprint(f"    [red]✗[/red] Error: {error_msg}")
            logger.error("Error processing region %s: %s", region.name, error_msg)

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
