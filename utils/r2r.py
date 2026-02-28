"""
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

import os
import re

from .runner import runner

R2R_META = "r2r.meta.txt"
R2R_INPUT = "r2r-input.sto.txt"
R2R_GSC = "r2r-input.gsc.sto.txt"
R2R_SVG = "r2r.svg"
R2R_SEQ_ID = "input"


def run_traveler(fasta_input, output_folder, seq_id):
    """Run traveler."""
    traveler_params = (
        f"--template-structure --file-format traveler "
        f"{output_folder}/traveler-template.xml "
        f"{fasta_input}"
    )
    cmd = (
        "traveler --verbose "
        f"--target-structure {fasta_input} {traveler_params} "
        f"--all {output_folder}/{seq_id} > {output_folder}/traveler.log"
    )
    runner.run(cmd)


def parse_fasta(fasta_input):
    """Parse a fasta file and return the sequence and structure."""
    with open(fasta_input, "r", encoding="utf-8") as f_in:
        lines = f_in.readlines()
    seq_id = ""
    sequence = ""
    structure = ""
    for i, line in enumerate(lines):
        if i == 0 and line.startswith(">"):
            seq_id = line.replace(">", "").split(" ")[0].strip()
        if i == 1:
            sequence = line.strip()
        if i == 2:
            structure = line.strip()
    if not sequence and not structure:
        raise ValueError("Invalid fasta file.")
    return seq_id, sequence, structure


def generate_r2r_input_file(sequence, structure, output_folder):
    """Generate an input file for R2R."""
    with open(os.path.join(output_folder, R2R_INPUT), "w", encoding="utf-8") as f_out:
        f_out.write("# STOCKHOLM 1.0\n")
        f_out.write(f"{R2R_SEQ_ID}          {sequence}\n")
        f_out.write(f"#=GC SS_cons   {structure}\n")
        f_out.write("#=GF R2R tick_label_disable_default_numbering\n")
        f_out.write("//\n")
    with open(os.path.join(output_folder, R2R_META), "w", encoding="utf-8") as f_out:
        f_out.write(f"{output_folder}/{R2R_GSC}\toneseq\t{R2R_SEQ_ID}\n")


def run_r2r(output_folder):
    """Run R2R: generate GSC alignment and SVG output."""
    cmd = (
        f"r2r --GSC-weighted-consensus {output_folder}/{R2R_INPUT} "
        f"{output_folder}/{R2R_GSC} 3 0.97 0.9 0.75 4 0.97 0.9 0.75 0.5 0.1"
    )
    runner.run(cmd)
    cmd = (
        f"r2r --disable-usage-warning {output_folder}/{R2R_META} "
        f"{output_folder}/{R2R_SVG} > {output_folder}/r2r.log"
    )
    runner.run(cmd)
    return f"{output_folder}/{R2R_SVG}"


def clean_r2r_output(output_folder):
    """Remove filename from R2R output."""
    r2r_svg_file = os.path.join(output_folder, R2R_SVG)
    with open(r2r_svg_file, "r", encoding="utf-8") as f_in:
        lines = f_in.readlines()
    with open(r2r_svg_file, "w", encoding="utf-8") as f_out:
        for line in lines:
            if R2R_GSC in line:
                continue
            f_out.write(line)


# Traveler-compatible CSS style block (matches Traveler output format).
_TRAVELER_CSS = """\
<g><style type="text/css" >
<!-- create color definitions -->
<![CDATA[
circle.red {stroke: rgb(255, 0, 255); fill: none; }
circle.green {stroke: rgb(0, 255, 0); fill: none; }
circle.blue {stroke: rgb(0, 0, 255); fill: none; }
circle.black {stroke: rgb(0, 0, 0); fill: none; }
circle.gray {stroke: rgb(204, 204, 204); fill: none; }
circle.brown {stroke: rgb(211.65, 104.55, 30.6); fill: none; }
line.red {stroke: rgb(255, 0, 255); stroke-width: 1.5; }
line.green {stroke: rgb(0, 255, 0); stroke-width: 1.5; }
line.blue {stroke: rgb(0, 0, 255); stroke-width: 1.5; }
line.black {stroke: rgb(0, 0, 0); stroke-width: 1.5; }
line.gray {stroke: rgb(204, 204, 204); stroke-width: 1.5; }
line.brown {stroke: rgb(211.65, 104.55, 30.6); stroke-width: 1.5; }
line {stroke: rgb(0, 0, 0); }
text.red {fill: rgb(255, 0, 255); font-size: 12px; font-weight: bold; font-family: Helvetica; }
text.green {fill: rgb(0, 255, 0); font-size: 12px; font-weight: bold; font-family: Helvetica; }
text.blue {fill: rgb(0, 0, 255); font-size: 12px; font-weight: bold; font-family: Helvetica; }
text.black {fill: rgb(0, 0, 0); font-size: 12px; font-weight: bold; font-family: Helvetica; }
text.gray {fill: rgb(204, 204, 204); font-size: 12px; font-weight: bold; font-family: Helvetica; }
text.brown {fill: rgb(211.65, 104.55, 30.6); font-size: 12px; font-weight: bold; font-family: Helvetica; }
text {fill: rgb(0, 0, 0); font-size: 12px; font-weight: bold; font-family: Helvetica; text-anchor: middle; alignment-baseline: middle; }
line.numbering-line{stroke: rgb(204, 204, 204); stroke-width: 0.5; }
polyline{fill: none; stroke-linejoin: round; }
template{visibility: hidden; }
text.numbering-label{fill: rgb(204, 204, 204); }
]]>
</style></g>"""


def r2r_svg_to_colored(r2r_svg_path, output_path):  # pylint: disable=too-many-locals
    """Convert an R2R SVG into a Traveler-compatible colored SVG.

    R2R produces a simple SVG with ``<text>``/``<tspan>`` elements for each
    nucleotide.  When Traveler cannot run (e.g. the dot-bracket is all dots),
    this function extracts nucleotide positions from the R2R SVG and writes
    them in Traveler's format so the rest of the pipeline works unchanged.

    Args:
        r2r_svg_path: Path to the R2R SVG (``traveler-template.svg``).
        output_path:  Destination path (e.g. ``{r2r_folder}/{seq_id}.colored.svg``).
    """
    with open(r2r_svg_path, "r", encoding="utf-8") as fh:
        svg_text = fh.read()

    # Extract nucleotide coordinates and letters from <text>/<tspan> elements.
    # R2R format:  <text x="19.245" y="23.985" ...><tspan ...>G</tspan></text>
    nucleotides = re.findall(
        r'<text\s+x="([^"]+)"\s+y="([^"]+)"[^>]*>'
        r"\s*<tspan[^>]*>([A-Za-z])</tspan>\s*</text>",
        svg_text,
    )
    if not nucleotides:
        return

    # Scale coordinates up (R2R uses ~7.5px font; Traveler uses ~12px)
    scale = 1.8
    padding = 15
    label_offset = 18  # vertical offset for 5'/3' labels
    font_half = 6  # half of font-size (12px) for text extent
    coords = []
    for x_str, y_str, letter in nucleotides:
        coords.append((float(x_str) * scale + padding, float(y_str) * scale, letter))

    # Compute bounding box including 5'/3' labels below first/last nucleotide
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    label_y_first = coords[0][1] + label_offset
    label_y_last = coords[-1][1] + label_offset
    min_x = min(xs) - padding
    min_y = min(ys + [label_y_first, label_y_last]) - font_half - padding
    max_x = max(xs) + padding
    max_y = max(ys + [label_y_first, label_y_last]) + font_half + padding
    vb_width = max_x - min_x
    vb_height = max_y - min_y

    lines = []
    lines.append(
        f'<svg\n        xmlns="http://www.w3.org/2000/svg"\n'
        f'        width="{vb_width:.2f}"\n        height="{vb_height:.2f}"\n'
        f'        viewBox="{min_x:.2f} {min_y:.2f} '
        f'{vb_width:.2f} {vb_height:.2f}">'
    )
    lines.append(_TRAVELER_CSS)

    # 5' / 3' labels
    first = coords[0]
    last = coords[-1]
    lines.append(
        f"<g><title>0 (position.label in template: 0.5')</title>"
        f'<text x="{first[0]:.3f}" y="{first[1] + 18:.3f}" '
        f'class="green" >5\'</text></g>'
    )

    # Backbone lines between consecutive nucleotides
    for i in range(len(coords) - 1):
        x1, y1, _ = coords[i]
        x2, y2, _ = coords[i + 1]
        lines.append(
            f'<g><line x1="{x1:.3f}" y1="{y1:.3f}" '
            f'x2="{x2:.3f}" y2="{y2:.3f}" class="gray" /></g>'
        )

    # Nucleotide text elements (Traveler format with <title>)
    for idx, (nx, ny, letter) in enumerate(coords, start=1):
        lines.append(
            f"<g><title>{idx} (position.label in template: {idx}.{letter})"
            f"</title>"
            f'<text x="{nx:.3f}" y="{ny:.3f}" class="black" '
            f">{letter}</text></g>"
        )

    # 3' label
    lines.append(
        f"<g><title>{len(coords) + 1} (position.label in template: "
        f"{len(coords) + 1}.3')</title>"
        f'<text x="{last[0]:.3f}" y="{last[1] + 18:.3f}" '
        f'class="green" >3\'</text></g>'
    )

    lines.append("</svg>")

    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
