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

# Minimum number of unpaired bases inside a hairpin loop for ViennaRNA.
# RNA.svg_rna_plot produces nan coordinates for hairpins smaller than this.
_MIN_HAIRPIN_SIZE = 3


def _remove_small_hairpins(structure, min_loop=_MIN_HAIRPIN_SIZE):
    """Remove base pairs that form hairpins with fewer than *min_loop* bases.

    ViennaRNA's layout engine cannot handle empty pairs ``()`` or very
    small hairpin loops like ``(.)`` — it produces ``nan`` coordinates.
    This function iteratively strips such pairs until all remaining
    hairpins have at least *min_loop* unpaired bases.

    Only ``(``, ``)``, and ``.`` characters are expected in *structure*.
    """
    chars = list(structure)
    changed = True
    while changed:
        changed = False
        # Build a pairing map using a stack
        stack = []
        pair_map = {}
        for i, char in enumerate(chars):
            if char == "(":
                stack.append(i)
            elif char == ")":
                if stack:
                    j = stack.pop()
                    pair_map[j] = i
                    pair_map[i] = j
        # Find hairpin-forming pairs with < min_loop unpaired bases
        for open_pos, close_pos in list(pair_map.items()):
            if open_pos >= close_pos:
                continue
            # Check if this pair forms a hairpin (no nested pairs inside)
            inner = chars[open_pos + 1 : close_pos]
            if all(c == "." for c in inner) and len(inner) < min_loop:
                chars[open_pos] = "."
                chars[close_pos] = "."
                changed = True
    return "".join(chars)


def _write_linear_colored_svg(rnapuzzler_svg, output_path):
    """Write a Traveler-compatible SVG for an all-dots (no pairs) structure.

    When hairpin removal strips every base pair, Traveler cannot run.
    This function reads nucleotide coordinates from the ViennaRNA SVG
    and writes them in Traveler's coloured SVG format so the rest of
    the pipeline (stitching, font adjustment, etc.) works unchanged.
    """
    coords = parse_svg_coordinates(rnapuzzler_svg)
    if not coords:
        return

    padding = 15
    font_half = 6
    label_offset = 18
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    label_y_first = coords[0][1] + label_offset
    label_y_last = coords[-1][1] + label_offset
    min_x = min(xs) - padding
    min_y = min(ys + [label_y_first, label_y_last]) - font_half - padding
    max_x = max(xs) + padding
    max_y = max(ys + [label_y_first, label_y_last]) + font_half + padding
    vb_w = max_x - min_x
    vb_h = max_y - min_y

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{vb_w:.2f}" height="{vb_h:.2f}" '
        f'viewBox="{min_x:.2f} {min_y:.2f} {vb_w:.2f} {vb_h:.2f}">',
        "<style>.black{font-size:12px;}</style>",
    ]
    # Backbone
    for i in range(len(coords) - 1):
        x1, y1, _ = coords[i]
        x2, y2, _ = coords[i + 1]
        lines.append(
            f'<line x1="{x1:.3f}" y1="{y1:.3f}" '
            f'x2="{x2:.3f}" y2="{y2:.3f}" '
            f'stroke="gray" stroke-width="1.5" />'
        )
    # Nucleotides
    for idx, (x, y, letter) in enumerate(coords):
        lines.append(
            f"<g><title>{idx} (position.label in template: {idx})</title>"
            f'<text x="{x:.3f}" y="{y:.3f}" class="black">{letter}</text></g>'
        )
    # 5'/3' labels
    lines.append(
        f"<g><title>0 (position.label in template: 0.5')</title>"
        f'<text x="{coords[0][0]:.3f}" y="{coords[0][1] + label_offset:.3f}" '
        f'class="green">5\'</text></g>'
    )
    lines.append(
        f"<g><title>{len(coords) - 1} "
        f"(position.label in template: {len(coords) - 1}.3')</title>"
        f'<text x="{coords[-1][0]:.3f}" y="{coords[-1][1] + label_offset:.3f}" '
        f'class="green">3\'</text></g>'
    )
    lines.append("</svg>")

    output_path = Path(output_path)
    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))


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

    # Remove pairs that form hairpins too small for ViennaRNA to lay out.
    clean_structure = _remove_small_hairpins(clean_structure)

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

    # Build two normalised structures:
    #
    # 1. *layout_structure* — for ViennaRNA coordinate generation.
    #    ViennaRNA only accepts () and ., so ALL bracket types are
    #    converted to () here.
    #
    # 2. *template_structure* — for the Traveler template FASTA.
    #    Traveler interprets () as nested pair-tree nodes and []{}
    #    as pseudoknot arcs.  The template FASTA must have the SAME
    #    pair-tree topology as the target, so pseudoknot characters
    #    must become dots (not parentheses) in the template.
    #
    # The target FASTA keeps the original structure unchanged so
    # pseudoknot arcs are still rendered in the final SVG.
    _open = set("([{<")
    _close = set(")]}>")

    layout_structure = "".join(
        "(" if c in _open else ")" if c in _close else "." for c in structure
    )
    template_structure = "".join(c if c in "()" else "." for c in structure)

    # Remove pairs that form hairpins too small for ViennaRNA to lay out.
    # This must happen before writing the template FASTA so that the
    # structure in the FASTA matches the coordinates in the XML template.
    layout_structure = _remove_small_hairpins(layout_structure)
    template_structure = _remove_small_hairpins(template_structure)

    # Also remove the same pairs from the target structure so Traveler
    # does not try to draw pairs that have no coordinates.
    # Compare against the original ()(). only template to detect removals.
    orig_template = "".join(c if c in "()" else "." for c in structure)
    target_structure = list(structure)
    for i, (old, new) in enumerate(zip(orig_template, template_structure)):
        if old != new:  # pair was removed
            target_structure[i] = "."
    target_structure = "".join(target_structure)

    # Template FASTA: cleaned structure matching the RNApuzzler layout.
    # Uses only () pairs — same tree topology as the target.
    template_fasta = output_dir / "rnapuzzler-template.fasta"
    with open(template_fasta, "w", encoding="utf-8") as fh:
        fh.write(f">{seq_id}\n{sequence}\n{template_structure}\n")

    # Target FASTA: structure preserving pseudoknots for Traveler,
    # but with small hairpins removed to match the template layout.
    target_fasta = output_dir / "rnapuzzler-target.fasta"
    with open(target_fasta, "w", encoding="utf-8") as fh:
        fh.write(f">{seq_id}\n{sequence}\n{target_structure}\n")

    # Step 1: Generate overlap-free layout
    svg_path = generate_layout_svg(sequence, layout_structure, output_dir)

    # Step 2: Convert to Traveler template XML
    xml_path = svg_to_traveler_template(svg_path, output_dir)

    # ViennaRNA already produces well-spaced coordinates (~25 units
    # between adjacent nucleotides), unlike R2R which needs 3× scaling.
    # No additional scaling is required here.

    # Step 3: Run Traveler — target keeps pseudoknots, template matches layout
    has_structure = any(c != "." for c in target_structure)
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
        # Build a minimal colored SVG directly from the ViennaRNA
        # coordinates (the R2R converter expects a different SVG format).
        _write_linear_colored_svg(
            svg_path,
            output_dir / f"{seq_id}.colored.svg",
        )
