"""
SVG namespacing utilities.

Provides functions to scope CSS rules in Traveler-generated SVGs so that
multiple diagrams can be safely embedded in the same document (e.g. via
stitching or inline embedding on a web page) without style leaks.

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

import math
import re
from pathlib import Path

_NON_CSS_IDENT_RE = re.compile(r"[^a-zA-Z0-9_-]")


def soften_long_basepair_lines(svg_path: str | Path, multiplier: float = 2.0) -> None:
    """Recolour abnormally long base-pair lines to light grey.

    Traveler draws base-pair lines between paired nucleotides.  In
    well-laid-out stems these are short and uniform, but in multi-way
    junctions or irregularly spaced layouts some lines can be many
    times longer than normal.  These long dark lines visually dominate
    the diagram.

    This function computes the **median** base-pair line length in the
    SVG, then recolours any line longer than ``multiplier × median``
    from ``class="black"`` to ``class="gray"``.  The threshold is
    self-calibrating: panels where all lines are similar length are
    left untouched.

    The function is a no-op when there are fewer than two base-pair
    lines or no outliers exceed the threshold.

    Args:
        svg_path: Path to the SVG file (modified in place).
        multiplier: Lines longer than ``multiplier * median`` are
            recoloured.  Default 2.0.
    """
    svg_path = Path(svg_path)
    content = svg_path.read_text(encoding="utf-8")

    # Match base-pair lines: <line x1=".." y1=".." x2=".." y2=".." class="black" />
    bp_pattern = re.compile(
        r"(<line\s+"
        r'x1="([\d.]+)"\s+y1="([\d.]+)"\s+'
        r'x2="([\d.]+)"\s+y2="([\d.]+)"\s+)'
        r'(class="black"\s*/>)'
    )

    matches = list(bp_pattern.finditer(content))
    if len(matches) < 2:
        return

    # Compute lengths
    lengths = []
    for m in matches:
        x1, y1 = float(m.group(2)), float(m.group(3))
        x2, y2 = float(m.group(4)), float(m.group(5))
        lengths.append(math.hypot(x2 - x1, y2 - y1))

    sorted_lengths = sorted(lengths)
    median = sorted_lengths[len(sorted_lengths) // 2]
    threshold = median * multiplier

    # Only proceed if there are actual outliers
    if sorted_lengths[-1] <= threshold:
        return

    # Rebuild content, replacing class on long lines.
    # Mark them "gray softened" so downstream tools (thumbnail,
    # outline) can identify and strip former base-pair lines.
    parts = []
    prev_end = 0
    for m, length in zip(matches, lengths):
        parts.append(content[prev_end : m.start()])
        if length > threshold:
            parts.append(m.group(1))
            parts.append('class="gray softened" />')
        else:
            parts.append(m.group(0))
        prev_end = m.end()
    parts.append(content[prev_end:])

    svg_path.write_text("".join(parts), encoding="utf-8")


def _collect_arc_points(chain):
    """Return polyline points for an insertion arc from *chain* nucleotides."""
    return [(nt["x"], nt["y"]) for nt in chain]


def _line_inside_bbox(line, bboxes, padding=15):
    """Return True if both endpoints of *line* fall inside any padded bbox."""
    for x_min, y_min, x_max, y_max in bboxes:
        x_lo = x_min - padding
        y_lo = y_min - padding
        x_hi = x_max + padding
        y_hi = y_max + padding
        if (
            x_lo <= line["x1"] <= x_hi
            and y_lo <= line["y1"] <= y_hi
            and x_lo <= line["x2"] <= x_hi
            and y_lo <= line["y2"] <= y_hi
        ):
            return True
    return False


_NT_RE = re.compile(
    r"<g><title>(\d+)\s[^<]*</title>"
    r'<text\s+x="([\d.]+)"\s+y="([\d.]+)"\s+'
    r'class="\w+"\s*>([^<]+)</text></g>'
)

_GRAY_LINE_RE = re.compile(
    r'(?:<g>)?<line\s+x1="([-\d.]+)"\s+y1="([-\d.]+)"\s+'
    r'x2="([-\d.]+)"\s+y2="([-\d.]+)"\s+'
    r'class="gray(?:\s+softened)?"\s*/>(?:</g>)?'
)


def _parse_nucleotides(content):
    """Return list of nucleotide dicts parsed from SVG *content*."""
    return [
        {
            "idx": int(m.group(1)),
            "x": float(m.group(2)),
            "y": float(m.group(3)),
            "letter": m.group(4).strip(),
            "start": m.start(),
            "end": m.end(),
        }
        for m in _NT_RE.finditer(content)
    ]


def _parse_backbone_lines(content):
    """Return list of gray-line dicts parsed from SVG *content*."""
    return [
        {
            "x1": float(m.group(1)),
            "y1": float(m.group(2)),
            "x2": float(m.group(3)),
            "y2": float(m.group(4)),
            "text": m.group(0),
        }
        for m in _GRAY_LINE_RE.finditer(content)
    ]


def _group_consecutive_x(all_nts):
    """Return groups of consecutive X nucleotides."""
    x_nts = [nt for nt in all_nts if nt["letter"] == "X"]
    if not x_nts:
        return []
    groups: list[list[dict]] = []
    current: list[dict] = [x_nts[0]]
    for nt in x_nts[1:]:
        if nt["idx"] == current[-1]["idx"] + 1:
            current.append(nt)
        else:
            groups.append(current)
            current = [nt]
    groups.append(current)
    return groups


def _max_deviation(points):
    """Return the max perpendicular distance of intermediate points from the chord.

    The chord is the straight line connecting the first and last point.
    Returns 0 when there are fewer than 3 points.
    """
    if len(points) < 3:
        return 0.0
    x0, y0 = points[0]
    x1, y1 = points[-1]
    chord = math.hypot(x1 - x0, y1 - y0)
    if chord < 1e-6:
        return 0.0
    max_d = 0.0
    for p_x, p_y in points[1:-1]:
        # Perpendicular distance from point to line (x0,y0)→(x1,y1)
        dist = abs((y1 - y0) * p_x - (x1 - x0) * p_y + x1 * y0 - y1 * x0) / chord
        max_d = max(max_d, dist)
    return max_d


def _build_arc_element(points, tip):
    """Build an SVG ``<g class="insertion-arc">`` element from *points*.

    When the backbone points naturally curve, the path traces through
    them as line segments.  When they are approximately collinear (a
    straight line), a quadratic Bézier with a perpendicular bulge is
    used so the arc remains visually distinct from the backbone.
    """
    if _max_deviation(points) < 3.0:
        # Points are approximately collinear — use a Bézier bulge
        x0, y0 = points[0]
        x1, y1 = points[-1]
        mid_x, mid_y = (x0 + x1) / 2, (y0 + y1) / 2
        chord = math.hypot(x1 - x0, y1 - y0)
        bulge = max(chord * 0.5, 20)
        # Perpendicular direction (rotate chord 90°)
        if chord > 1e-6:
            ctrl_x = mid_x - ((y1 - y0) / chord) * bulge
            ctrl_y = mid_y + ((x1 - x0) / chord) * bulge
        else:
            ctrl_x = mid_x + bulge
            ctrl_y = mid_y
        path_d = (
            f"M {x0:.1f},{y0:.1f} " f"Q {ctrl_x:.1f},{ctrl_y:.1f} {x1:.1f},{y1:.1f}"
        )
    else:
        # Natural curve — trace through backbone waypoints
        path_d = f"M {points[0][0]:.1f},{points[0][1]:.1f}"
        for p_x, p_y in points[1:]:
            path_d += f" L {p_x:.1f},{p_y:.1f}"
    return (
        f'<g class="insertion-arc">'
        f"<title>{tip}</title>"
        f'<path d="{path_d}" '
        f'stroke="black" stroke-width="2" '
        f'stroke-linecap="round" stroke-linejoin="round" fill="none"/>'
        f"</g>"
    )


def replace_xxxx_with_arc(  # pylint: disable=too-many-branches
    svg_path: str | Path,
    insertion_sizes: list[int] | None = None,
) -> None:
    """Replace XXXX placeholder nucleotides with a smooth insertion arc.

    When large insertions are removed from the input sequence they are
    replaced with four ``X`` characters.  Traveler renders these as
    literal ``X`` text elements in the SVG.  This function finds each
    run of consecutive X nucleotides, removes the ``X`` text elements,
    and replaces them with a thick arc that follows the nucleotide
    positions.  Gray lines whose both endpoints fall inside the arc
    bounding box are also removed so no stray lines show behind the arc.

    A native SVG ``<title>`` element is added so that hovering over
    the arc shows how many nucleotides were removed.

    Args:
        svg_path: Path to the SVG file (modified in place).
        insertion_sizes: Number of nucleotides removed for each XXXX
            group, in left-to-right order.  When *None* the tooltip
            says "nucleotides not shown" without a count.
    """
    svg_path = Path(svg_path)
    if not svg_path.exists() or svg_path.stat().st_size == 0:
        return

    content = svg_path.read_text(encoding="utf-8")
    all_nts = _parse_nucleotides(content)
    if not all_nts:
        return

    by_index = {nt["idx"]: nt for nt in all_nts}
    groups = _group_consecutive_x(all_nts)
    if not groups:
        return

    # Collect bounding boxes for every arc so we can remove gray lines later.
    arc_bboxes: list[tuple[float, float, float, float]] = []

    # Process in reverse so earlier string offsets stay valid.
    for rev_i, group in enumerate(reversed(groups)):
        fwd_i = len(groups) - 1 - rev_i
        before = by_index.get(group[0]["idx"] - 1)
        after = by_index.get(group[-1]["idx"] + 1)

        if before is None and after is None:
            continue
        if before is None:
            before = after
        if after is None:
            after = before

        chain = [before] + list(group) + [after]
        points = _collect_arc_points(chain)
        if len(points) < 2:
            continue

        # Bounding box of this arc's points.
        xs = [p[0] for p in points]
        ys = [p[1] for p in points]
        arc_bboxes.append((min(xs), min(ys), max(xs), max(ys)))

        if insertion_sizes and fwd_i < len(insertion_sizes):
            tip = f"{insertion_sizes[fwd_i]:,} nucleotides not shown"
        else:
            tip = "Nucleotides not shown"

        arc = _build_arc_element(points, tip)
        content = content[: group[0]["start"]] + arc + content[group[-1]["end"] :]

    # Remove gray lines that fall inside an arc bounding box.
    if arc_bboxes:
        gray_lines = _parse_backbone_lines(content)
        for gline in reversed(gray_lines):
            if _line_inside_bbox(gline, arc_bboxes):
                content = content.replace(gline["text"], "", 1)

    svg_path.write_text(content, encoding="utf-8")


def _sanitize_css_id(scope_id: str) -> str:
    """Sanitize a string so it is a valid CSS identifier.

    CSS identifiers cannot start with a digit and must not contain dots,
    spaces, or other special characters.  We replace all non-alphanumeric
    characters (except hyphens and underscores) with hyphens.
    """
    sanitized = _NON_CSS_IDENT_RE.sub("-", scope_id)
    # Ensure it doesn't start with a digit or hyphen-digit
    if sanitized and (sanitized[0].isdigit() or sanitized[0] == "-"):
        sanitized = "r" + sanitized
    return sanitized


def namespace_svg_file(svg_path: str, scope_id: str) -> None:
    """Namespace an SVG file on disk, rewriting it in place.

    Wraps the SVG content in a scoped ``<g>`` element and prefixes all CSS
    selectors in ``<style>`` blocks so they only apply within that scope.

    Args:
        svg_path: Path to the SVG file.
        scope_id: Unique identifier used as the wrapper ``id`` attribute
            (e.g. ``"r2dt-URS123"``).
    """
    scope_id = _sanitize_css_id(scope_id)
    path = Path(svg_path)
    if not path.exists() or path.stat().st_size == 0:
        return
    content = path.read_text(encoding="utf-8")
    content = namespace_svg_text(content, scope_id)
    path.write_text(content, encoding="utf-8")


def namespace_svg_text(svg_text: str, scope_id: str) -> str:
    """Namespace raw SVG text.

    Scopes all CSS selectors inside ``<style>`` blocks so they are prefixed
    with ``#scope_id``.  A wrapper ``<g id="scope_id">`` is inserted around
    the SVG body content (everything between the opening and closing ``<svg>``
    tags) so the scoped selectors match.

    Args:
        svg_text: Raw SVG markup.
        scope_id: Unique scope identifier.

    Returns:
        The rewritten SVG text with scoped CSS.
    """
    scope_id = _sanitize_css_id(scope_id)

    # 1. Scope CSS selectors inside every <style> block
    svg_text = _scope_style_blocks(svg_text, scope_id)

    # 2. Wrap body content in a <g id="scope_id">
    svg_text = _wrap_body(svg_text, scope_id)

    return svg_text


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

_STYLE_RE = re.compile(r"(<style[^>]*>)(.*?)(</style>)", re.DOTALL | re.IGNORECASE)

_CDATA_RE = re.compile(r"(<!\[CDATA\[)(.*?)(\]\]>)", re.DOTALL)


def _scope_style_blocks(svg_text: str, scope_id: str) -> str:
    """Prefix every CSS selector inside <style> blocks with ``#scope_id``.

    Traveler wraps CSS inside ``<![CDATA[...]]>`` sections and may include
    XML comments like ``<!-- create color definitions -->``.  We must scope
    only the actual CSS text inside the CDATA block and preserve the XML
    structure around it so that downstream XML parsers (e.g. cairosvg) can
    still read the SVG correctly.
    """

    def _rewrite(match: re.Match) -> str:
        open_tag = match.group(1)
        raw = match.group(2)
        close_tag = match.group(3)

        # If there is a CDATA block, scope only the CSS inside it
        if cdata_match := _CDATA_RE.search(raw):
            before = raw[: cdata_match.start()]
            cdata_open = cdata_match.group(1)  # <![CDATA[
            css = cdata_match.group(2)
            cdata_close = cdata_match.group(3)  # ]]>
            after = raw[cdata_match.end() :]
            scoped_css = _scope_css(css, scope_id)
            return (
                f"{open_tag}{before}{cdata_open}"
                f"{scoped_css}{cdata_close}{after}{close_tag}"
            )

        # No CDATA — scope the raw text directly
        scoped_css = _scope_css(raw, scope_id)
        return f"{open_tag}{scoped_css}{close_tag}"

    return _STYLE_RE.sub(_rewrite, svg_text)


# Regex that matches a CSS rule: selector(s) { declarations }
# Handles multi-line selectors and declarations.
_RULE_RE = re.compile(
    r"([^{}]+?)\{([^{}]*)\}",
    re.DOTALL,
)


def _scope_css(css: str, scope_id: str) -> str:
    """Prefix each selector in a CSS block with ``#scope_id``."""

    def _rewrite_rule(match: re.Match) -> str:
        selectors_raw = match.group(1)
        declarations = match.group(2)
        # Split comma-separated selectors and prefix each one
        selectors = selectors_raw.split(",")
        scoped = []
        for sel in selectors:
            sel = sel.strip()
            if not sel:
                continue
            scoped.append(f"#{scope_id} {sel}")
        return ",".join(scoped) + "{" + declarations + "}"

    return _RULE_RE.sub(_rewrite_rule, css)


def _wrap_body(svg_text: str, scope_id: str) -> str:
    """Insert ``<g id="scope_id">`` right after the opening ``<svg>`` tag
    and ``</g>`` right before the closing ``</svg>`` tag."""
    # Find the end of the opening <svg ...> tag
    open_match = re.search(r"<svg\b[^>]*>", svg_text, re.IGNORECASE)
    if not open_match:
        return svg_text

    # Find the last </svg>
    close_idx = svg_text.rfind("</svg>")
    if close_idx == -1:
        return svg_text

    insert_after = open_match.end()
    wrapper_open = f'<g id="{scope_id}">'
    wrapper_close = "</g>"

    return (
        svg_text[:insert_after]
        + wrapper_open
        + svg_text[insert_after:close_idx]
        + wrapper_close
        + svg_text[close_idx:]
    )
