"""
RNApuzzler layout engine for overlap-free RNA 2D structure drawings.

Uses ViennaRNA's RNAplot with layout type 4 (RNApuzzler) via the Python
bindings that are already installed in the Docker image.  The algorithm
is guaranteed to produce intersection-free (outerplanar) drawings.

Reference:
    Wiegreffe et al. (2018) "RNApuzzler: efficient outerplanar drawing
    of RNA secondary structures", Bioinformatics.
"""

import re
from pathlib import Path

import RNA  # ViennaRNA Python bindings  # pylint: disable=import-error

from .runner import runner

# ViennaRNA layout type constants
_LAYOUT_PUZZLER = 4


def generate_layout_svg(sequence, structure, output_dir):
    """Generate an overlap-free SVG layout using ViennaRNA's RNApuzzler.

    Calls ``RNA.svg_rna_plot`` with ``rna_plot_type = 4`` (RNApuzzler)
    to produce an SVG file containing ``<text>`` elements with (x, y)
    coordinates for every nucleotide.

    Args:
        sequence:   RNA sequence string (e.g. ``"GGGAAACCC"``).
        structure:  Dot-bracket secondary structure.
        output_dir: Directory to write the SVG into.

    Returns:
        Path to the generated SVG file.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    svg_path = output_dir / "rnapuzzler.svg"

    # ViennaRNA only understands () and . — strip pseudoknot characters
    clean_structure = re.sub(r"[^().]", ".", structure)

    # Set layout type to RNApuzzler (type 4)
    RNA.cvar.rna_plot_type = _LAYOUT_PUZZLER
    result = RNA.svg_rna_plot(sequence, clean_structure, str(svg_path))

    if result != 1 or not svg_path.exists():
        raise RuntimeError(
            f"RNApuzzler (RNA.svg_rna_plot) failed for "
            f"sequence of length {len(sequence)}"
        )
    return svg_path


def parse_svg_coordinates(svg_path):
    """Extract nucleotide (x, y, letter) tuples from an RNApuzzler SVG.

    ViennaRNA's ``svg_rna_plot`` writes nucleotides as::

        <text x="123.456" y="789.012">A</text>

    inside a ``<g>`` element with ``id="seq"``.

    Args:
        svg_path: Path to the SVG produced by :func:`generate_layout_svg`.

    Returns:
        List of ``(x, y, nucleotide)`` tuples.
    """
    with open(svg_path, "r", encoding="utf-8") as fh:
        content = fh.read()

    # ViennaRNA wraps nucleotides in <g id="seq" ...> ... </g>
    # Each nucleotide is: <text x="..." y="...">N</text>
    seq_block_match = re.search(r'<g[^>]*id="seq"[^>]*>(.*?)</g>', content, re.DOTALL)
    if seq_block_match:
        block = seq_block_match.group(1)
    else:
        # Fallback: scan the whole file
        block = content

    pattern = re.compile(
        r'<text\s[^>]*?x="([^"]+)"\s*y="([^"]+)"[^>]*>([A-Za-z])</text>'
    )
    coords = []
    for match in pattern.finditer(block):
        x_val = float(match.group(1))
        y_val = float(match.group(2))
        nucleotide = match.group(3)
        coords.append((x_val, y_val, nucleotide))

    if not coords:
        raise ValueError(f"No nucleotide coordinates found in {svg_path}")

    return coords


def svg_to_traveler_template(svg_path, output_dir):
    """Convert an RNApuzzler SVG into a Traveler template XML file.

    Produces ``traveler-template.xml`` with the format::

        <structure>
          <point x="..." y="..." b="N" />
          ...
        </structure>

    Args:
        svg_path:   Path to the RNApuzzler SVG.
        output_dir: Directory to write the template into.

    Returns:
        Path to the generated ``traveler-template.xml``.
    """
    output_dir = Path(output_dir)
    coords = parse_svg_coordinates(svg_path)

    xml_path = output_dir / "traveler-template.xml"
    with open(xml_path, "w", encoding="utf-8") as fh:
        fh.write("<structure>\n")
        for x_val, y_val, nucleotide in coords:
            fh.write(f'<point x="{x_val:.3f}" y="{y_val:.3f}" b="{nucleotide}" />\n')
        fh.write("</structure>\n")

    return xml_path


def run_puzzler_pipeline(_fasta_input, output_dir, seq_id, sequence, structure):
    """Full RNApuzzler pipeline: layout → Traveler template → Traveler render.

    This is the high-level entry point that mirrors the R2R and RNArtist
    pipelines.  It:

    1. Generates an overlap-free SVG via ViennaRNA's RNApuzzler
    2. Converts the SVG to a Traveler template XML
    3. Scales coordinates (3×) to avoid Traveler precision issues
    4. Runs Traveler to produce the final colored SVG

    Args:
        fasta_input: Path to 3-line FASTA (header, sequence, structure).
        output_dir:  Working directory for intermediate files.
        seq_id:      Sequence identifier.
        sequence:    RNA sequence string.
        structure:   Dot-bracket secondary structure.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Normalise the structure for ViennaRNA layout generation.
    # ViennaRNA only accepts () and . — pseudoknot characters (Aa, Bb,
    # [], {} etc.) must be stripped for the coordinate calculation.
    # However, the *original* structure is passed to Traveler so that
    # pseudoknot arcs are still rendered in the final SVG.
    _open = set("([{<")
    _close = set(")]}>")
    norm_structure = "".join(
        "(" if c in _open else ")" if c in _close else "." for c in structure
    )

    # Template FASTA: normalised structure matching the RNApuzzler layout.
    template_fasta = output_dir / "rnapuzzler-template.fasta"
    with open(template_fasta, "w", encoding="utf-8") as fh:
        fh.write(f">{seq_id}\n{sequence}\n{norm_structure}\n")

    # Target FASTA: original structure preserving pseudoknots for Traveler.
    target_fasta = output_dir / "rnapuzzler-target.fasta"
    with open(target_fasta, "w", encoding="utf-8") as fh:
        fh.write(f">{seq_id}\n{sequence}\n{structure}\n")

    # Step 1: Generate overlap-free layout
    svg_path = generate_layout_svg(sequence, norm_structure, output_dir)

    # Step 2: Convert to Traveler template XML
    xml_path = svg_to_traveler_template(svg_path, output_dir)

    # ViennaRNA already produces well-spaced coordinates (~25 units
    # between adjacent nucleotides), unlike R2R which needs 3× scaling.
    # No additional scaling is required here.

    # Step 3: Run Traveler — target keeps pseudoknots, template matches layout
    has_structure = any(c != "." for c in structure)
    if has_structure:
        cmd = (
            f"traveler --verbose "
            f"--target-structure {target_fasta} "
            f"--template-structure --file-format traveler "
            f"{xml_path} {template_fasta} "
            f"--all {output_dir}/{seq_id} "
            f"> {output_dir}/traveler.log"
        )
        runner.run(cmd)
    else:
        # No base pairs — Traveler crashes on all-dots structures.
        # Fall back to the R2R direct-SVG converter.
        from . import r2r  # pylint: disable=import-outside-toplevel

        r2r.r2r_svg_to_colored(
            str(output_dir / "rnapuzzler.svg"),
            str(output_dir / f"{seq_id}.colored.svg"),
        )
