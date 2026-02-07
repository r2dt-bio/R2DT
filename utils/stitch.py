# pylint: disable=too-many-lines

"""
Stitch multiple R2DT SVG diagrams into one combined SVG.

This module provides functions to place multiple R2DT secondary structure SVGs
into one combined SVG, arranged left-to-right, with panel i's 3′ visually
joining panel i+1's 5′.

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
import sys
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path
from typing import NamedTuple, Optional

# SVG namespace
SVG_NS = "http://www.w3.org/2000/svg"
NSMAP = {"svg": SVG_NS}

# Register default namespace to avoid ns0: prefixes
ET.register_namespace("", SVG_NS)


class Anchor(NamedTuple):
    """Represents a 5' or 3' anchor position."""

    x: float
    y: float


class ViewBox(NamedTuple):
    """SVG viewBox parameters."""

    minx: float
    miny: float
    width: float
    height: float


# Font size multipliers relative to detected nucleotide font size
# These create a visual hierarchy where labels are appropriately sized
FONT_MULTIPLIERS = {
    "anchor_label": 4.0,  # 5'/3' labels: prominent, ~4x nucleotide size
    "caption": 6.0,  # Captions: largest, ~6x nucleotide size
    "gap_label": 2.5,  # Gap distance labels: moderate, ~2.5x nucleotide size
}

# Default nucleotide font size if detection fails
DEFAULT_NUCLEOTIDE_FONT_SIZE = 10.0


class PanelData(NamedTuple):
    """Data extracted from a single SVG panel."""

    filepath: Path
    root: ET.Element
    viewbox: ViewBox
    anchor_5: Anchor
    anchor_3: Anchor
    anchor_5_elem: Optional[ET.Element]  # Reference to 5' text element
    anchor_3_elem: Optional[ET.Element]  # Reference to 3' text element
    children: list
    sort_key: int  # Genomic start position for sorting
    genomic_start: int  # Genomic start position
    genomic_end: int  # Genomic end position
    visual_bbox: Optional["BBox"] = None  # Visual bounding box (set after parsing)
    font_size: float = DEFAULT_NUCLEOTIDE_FONT_SIZE  # Detected nucleotide font size


def parse_viewbox(svg_root: ET.Element) -> ViewBox:
    """
    Extract viewBox from SVG root element.
    Falls back to width/height attributes if viewBox is missing.
    """
    if viewbox_attr := svg_root.get("viewBox"):
        parts = viewbox_attr.split()
        if len(parts) == 4:
            return ViewBox(
                minx=float(parts[0]),
                miny=float(parts[1]),
                width=float(parts[2]),
                height=float(parts[3]),
            )

    # Fallback to width/height
    width = svg_root.get("width")
    height = svg_root.get("height")

    if width and height:
        # Strip units if present (e.g., "100px" -> "100")
        width_val = float(re.sub(r"[a-zA-Z%]+", "", width))
        height_val = float(re.sub(r"[a-zA-Z%]+", "", height))
        return ViewBox(minx=0, miny=0, width=width_val, height=height_val)

    raise ValueError("SVG has no viewBox and no parseable width/height attributes")


def get_text_content(elem: ET.Element) -> str:
    """
    Get all text content from an element, including nested tspan elements.
    """
    texts = []
    if elem.text:
        texts.append(elem.text)
    for child in elem:
        texts.append(get_text_content(child))
        if child.tail:
            texts.append(child.tail)
    return "".join(texts)


def find_anchor(
    svg_root: ET.Element, pattern: str, anchor_name: str, filepath: Path
) -> tuple[Anchor, Optional[ET.Element]]:
    """
    Find 5' or 3' anchor position in SVG.

    Args:
        svg_root: Root SVG element
        pattern: Regex pattern to match (e.g., r"5[′']" for 5')
        anchor_name: Name for error messages (e.g., "5′")
        filepath: Source file path for error messages

    Returns:
        Tuple of (Anchor position, text element reference)

    Raises:
        ValueError: If anchor not found or missing coordinates
    """
    candidates = []

    # Find all text elements (with or without namespace)
    for elem in svg_root.iter():
        tag_local = elem.tag.split("}")[-1] if "}" in elem.tag else elem.tag
        if tag_local == "text":
            content = get_text_content(elem)
            if re.search(pattern, content):
                x_attr = elem.get("x")
                y_attr = elem.get("y")

                if x_attr and y_attr:
                    # Handle multiple values (use first)
                    x = float(x_attr.split()[0])
                    y = float(y_attr.split()[0])
                    candidates.append((Anchor(x=x, y=y), elem))

    if not candidates:
        raise ValueError(
            f"Could not find {anchor_name} anchor in {filepath}; "
            "ensure R2DT outputs include end labels with x/y."
        )

    if len(candidates) > 1:
        print(
            f"Warning: Multiple {anchor_name} labels found in {filepath}; using the first",
            file=sys.stderr,
        )

    return candidates[0]


def style_anchor_label(
    elem: ET.Element, font_size: float, font_weight: str = "bold"
) -> None:
    """
    Style a 5' or 3' anchor label element to make it more distinct.
    Uses inline style attribute to override CSS rules in <style> blocks.

    Args:
        elem: The text element to style
        font_size: Font size to apply
        font_weight: Font weight (default: bold)
    """
    # Get existing style and strip any existing font-size to avoid conflicts
    existing_style = elem.get("style", "")
    # Remove any existing font-size declarations from the style (handles optional semicolon)
    existing_style = re.sub(r"font-size:\s*[\d.]+px;?\s*", "", existing_style)
    # Build new style with our font-size
    # pylint: disable-next=line-too-long
    new_style = f"font-size: {font_size}px; font-weight: {font_weight}; font-family: Helvetica, Arial, sans-serif;"
    if existing_style.strip():
        elem.set("style", f"{new_style} {existing_style.strip()}")
    else:
        elem.set("style", new_style)


def extract_genomic_coords(filepath: Path) -> tuple[int, int]:
    """
    Extract genomic start and end positions from filename.
    Expects format like: name_START-END.svg or name_START-END_...

    Returns:
        Tuple of (start, end) positions
    """
    filename = filepath.stem
    # Try pattern with underscores on both sides: _START-END_
    match = re.search(r"_(\d+)-(\d+)_", filename)
    if match:
        return int(match.group(1)), int(match.group(2))
    # Try pattern with underscore before and end of string: _START-END$
    match = re.search(r"_(\d+)-(\d+)$", filename)
    if match:
        return int(match.group(1)), int(match.group(2))
    # Fallback: try to find any START-END pattern
    match = re.search(r"(\d+)-(\d+)", filename)
    if match:
        return int(match.group(1)), int(match.group(2))
    return 0, 0  # Default if no number found


# pylint: disable-next=too-many-branches,too-many-statements
def detect_nucleotide_font_size(svg_root: ET.Element) -> float:
    """
    Detect the font size used for nucleotide letters in an SVG.

    Looks for font-size in CSS style blocks and on text elements.
    Returns the most common font size found, or a default.

    Args:
        svg_root: Root SVG element

    Returns:
        Detected font size in SVG units, or DEFAULT_NUCLEOTIDE_FONT_SIZE if detection fails
    """
    font_sizes = []
    nucleotide_pattern = re.compile(r"^[ACGUT]$", re.IGNORECASE)

    # First, try to extract font-size from CSS style blocks
    css_font_size = None
    for elem in svg_root.iter():
        tag_local = elem.tag.split("}")[-1] if "}" in elem.tag else elem.tag
        if tag_local == "style" and elem.text:
            # Look for "text { ... font-size: Xpx; ... }" pattern
            text_rule = re.search(r"text\s*\{[^}]*font-size:\s*([\d.]+)", elem.text)
            if text_rule:
                try:
                    css_font_size = float(text_rule.group(1))
                    break
                except ValueError:
                    pass
            # Also try ".font { ... font-size: Xpx; ... }" pattern (R2DT output)
            if css_font_size is None:
                font_rule = re.search(
                    r"\.font\s*\{[^}]*font-size:\s*([\d.]+)", elem.text
                )
                if font_rule:
                    try:
                        css_font_size = float(font_rule.group(1))
                        break
                    except ValueError:
                        pass

    # If found in CSS, return it directly
    if css_font_size is not None and css_font_size > 0:
        return css_font_size

    # Otherwise, scan individual text elements
    for elem in svg_root.iter():
        tag_local = elem.tag.split("}")[-1] if "}" in elem.tag else elem.tag
        if tag_local != "text":
            continue

        # Skip numbering labels and other special elements
        elem_class = elem.get("class", "")
        if "numbering-label" in elem_class:
            continue

        # Get text content
        content = get_text_content(elem).strip()
        if not nucleotide_pattern.match(content):
            continue

        # Try to get font-size from various sources
        font_size = None

        # Check direct font-size attribute
        fs_attr = elem.get("font-size")
        if fs_attr:
            try:
                font_size = float(re.sub(r"[a-zA-Z%]+", "", fs_attr))
            except ValueError:
                pass

        # Check style attribute
        if font_size is None:
            style = elem.get("style", "")
            fs_match = re.search(r"font-size:\s*([\d.]+)", style)
            if fs_match:
                try:
                    font_size = float(fs_match.group(1))
                except ValueError:
                    pass

        if font_size is not None and font_size > 0:
            font_sizes.append(font_size)

    if not font_sizes:
        return DEFAULT_NUCLEOTIDE_FONT_SIZE

    # Return the most common font size (mode)
    counter = Counter(font_sizes)
    return counter.most_common(1)[0][0]


# pylint: disable-next=too-many-branches
def normalize_text_font_sizes(svg_root: ET.Element, default_font_size: float) -> None:
    """
    Add inline style font-size to all text elements and remove font-size from CSS.

    In SVG, CSS rules in <style> blocks override presentation attributes (font-size="X"),
    but inline style attributes (style="font-size: Xpx") have the highest priority.
    This function:
    1. Extracts font-size values from CSS rules
    2. Applies them as inline style attributes to text elements
    3. Removes font-size from CSS rules to prevent override issues

    Args:
        svg_root: Root SVG element
        default_font_size: Font size to use for text elements without explicit font-size
    """
    # First, extract font-sizes from CSS rules in the SVG
    css_font_sizes = {}  # Maps class name to font-size

    for elem in svg_root.iter():
        tag_local = elem.tag.split("}")[-1] if "}" in elem.tag else elem.tag
        if tag_local == "style" and elem.text:
            # Extract font-size for each CSS class
            # Pattern: .classname { ... font-size: Xpx; ... }
            for match in re.finditer(
                r"\.?([\w-]+)\s*\{[^}]*font-size:\s*([\d.]+)", elem.text
            ):
                class_name = match.group(1)
                try:
                    font_size = float(match.group(2))
                    css_font_sizes[class_name] = font_size
                except ValueError:
                    pass

            # Also match "text.classname" or just "text"
            text_match = re.search(r"text\s*\{[^}]*font-size:\s*([\d.]+)", elem.text)
            if text_match:
                try:
                    css_font_sizes["text"] = float(text_match.group(1))
                except ValueError:
                    pass

            # Remove font-size from all CSS rules to prevent them from overriding inline styles
            elem.text = re.sub(r"font-size:\s*[\d.]+px?;\s*", "", elem.text)

    # Now apply font-size as inline style to all text elements
    for elem in svg_root.iter():
        tag_local = elem.tag.split("}")[-1] if "}" in elem.tag else elem.tag
        if tag_local != "text":
            continue

        # Check if element already has font-size in style attribute
        style = elem.get("style", "")
        if "font-size" in style:
            continue

        # Get font-size from presentation attribute if present
        font_size = None
        fs_attr = elem.get("font-size")
        if fs_attr:
            try:
                font_size = float(re.sub(r"[a-zA-Z%]+", "", fs_attr))
            except ValueError:
                pass

        # Otherwise, find font-size from CSS based on element's class
        if font_size is None:
            elem_classes = elem.get("class", "").split()
            for class_name in elem_classes:
                if class_name in css_font_sizes:
                    font_size = css_font_sizes[class_name]
                    break

        # Fall back to 'text' rule if no class-specific rule found
        if font_size is None and "text" in css_font_sizes:
            font_size = css_font_sizes["text"]

        # Fall back to default
        if font_size is None:
            font_size = default_font_size

        # Add font-size to inline style (highest CSS priority)
        if style:
            new_style = f"font-size: {font_size}px; {style}"
        else:
            new_style = f"font-size: {font_size}px"
        elem.set("style", new_style)

        # Remove the presentation attribute since we're using style now
        if "font-size" in elem.attrib:
            del elem.attrib["font-size"]


def extract_nucleotide_positions(svg_root: ET.Element) -> list[tuple[float, float]]:
    """
    Extract nucleotide (text element) positions from SVG.
    Returns list of (x, y) coordinates in document order.
    Skips numbering labels.
    """
    points = []
    for elem in svg_root.iter():
        tag_local = elem.tag.split("}")[-1] if "}" in elem.tag else elem.tag
        if tag_local == "text":
            # Skip numbering labels
            elem_class = elem.get("class", "")
            if "numbering-label" in elem_class:
                continue

            x_attr = elem.get("x")
            y_attr = elem.get("y")
            if x_attr and y_attr:
                try:
                    x = float(x_attr.split()[0])
                    y = float(y_attr.split()[0])
                    points.append((x, y))
                except ValueError:
                    continue
    return points


class BBox(NamedTuple):
    """Visual bounding box of an SVG element or panel."""

    min_x: float
    min_y: float
    max_x: float
    max_y: float

    @property
    def width(self) -> float:
        """Return the width of the bounding box."""
        return self.max_x - self.min_x

    @property
    def height(self) -> float:
        """Return the height of the bounding box."""
        return self.max_y - self.min_y


# pylint: disable-next=too-many-branches,too-many-statements
def calculate_visual_bbox(svg_root: ET.Element, viewbox: ViewBox) -> BBox:
    """
    Calculate the actual visual bounding box of SVG content.

    This scans all visual elements (text, circles, lines, paths) to find
    the true extent of the drawing, which may differ from the viewBox.

    Args:
        svg_root: Root SVG element
        viewbox: The SVG's declared viewBox (used as fallback)

    Returns:
        BBox with the actual visual extent
    """
    min_x = float("inf")
    min_y = float("inf")
    max_x = float("-inf")
    max_y = float("-inf")

    # Estimate text element size (half of typical font size as radius)
    text_padding = 6.0

    for elem in svg_root.iter():
        tag_local = elem.tag.split("}")[-1] if "}" in elem.tag else elem.tag

        if tag_local == "text":
            x_attr = elem.get("x")
            y_attr = elem.get("y")
            if x_attr and y_attr:
                try:
                    x = float(x_attr.split()[0])
                    y = float(y_attr.split()[0])
                    min_x = min(min_x, x - text_padding)
                    max_x = max(max_x, x + text_padding)
                    min_y = min(min_y, y - text_padding)
                    max_y = max(max_y, y + text_padding)
                except ValueError:
                    pass

        elif tag_local == "circle":
            cx = elem.get("cx")
            cy = elem.get("cy")
            r = elem.get("r", "0")
            if cx and cy:
                try:
                    cx_val = float(cx)
                    cy_val = float(cy)
                    r_val = float(r)
                    min_x = min(min_x, cx_val - r_val)
                    max_x = max(max_x, cx_val + r_val)
                    min_y = min(min_y, cy_val - r_val)
                    max_y = max(max_y, cy_val + r_val)
                except ValueError:
                    pass

        elif tag_local == "line":
            x1 = elem.get("x1")
            y1 = elem.get("y1")
            x2 = elem.get("x2")
            y2 = elem.get("y2")
            if x1 and y1 and x2 and y2:
                try:
                    min_x = min(min_x, float(x1), float(x2))
                    max_x = max(max_x, float(x1), float(x2))
                    min_y = min(min_y, float(y1), float(y2))
                    max_y = max(max_y, float(y1), float(y2))
                except ValueError:
                    pass

        elif tag_local == "rect":
            x = elem.get("x", "0")
            y = elem.get("y", "0")
            w = elem.get("width", "0")
            h = elem.get("height", "0")
            try:
                x_val = float(x)
                y_val = float(y)
                w_val = float(w)
                h_val = float(h)
                min_x = min(min_x, x_val)
                max_x = max(max_x, x_val + w_val)
                min_y = min(min_y, y_val)
                max_y = max(max_y, y_val + h_val)
            except ValueError:
                pass

    # Fallback to viewBox if no elements found
    if min_x == float("inf"):
        return BBox(
            min_x=viewbox.minx,
            min_y=viewbox.miny,
            max_x=viewbox.minx + viewbox.width,
            max_y=viewbox.miny + viewbox.height,
        )

    return BBox(min_x=min_x, min_y=min_y, max_x=max_x, max_y=max_y)


# pylint: disable-next=too-many-arguments,too-many-positional-arguments
def create_outline_path(
    panels: list,
    translations: list[tuple[float, float]],
    scale_factors: list[float] = None,
    stroke_color: str = "#cccccc",
    stroke_width: float = 3.0,
    stroke_opacity: float = 0.6,
) -> Optional[ET.Element]:
    """
    Create an SVG path element that traces through all nucleotide positions
    across all panels, creating a connecting outline.

    Args:
        panels: List of PanelData objects
        translations: List of (tx, ty) translation offsets for each panel
        scale_factors: Optional list of scale factors for each panel
        stroke_color: Color of the outline stroke
        stroke_width: Width of the outline stroke
        stroke_opacity: Opacity of the outline (0-1)

    Returns:
        SVG path element, or None if insufficient points
    """
    all_points = []

    for i, (panel, (tx, ty)) in enumerate(zip(panels, translations)):
        # Get scale factor for this panel
        scale = scale_factors[i] if scale_factors else 1.0

        # Extract nucleotide positions from this panel
        local_points = extract_nucleotide_positions(panel.root)

        # Transform to global coordinates (apply scale then translate)
        for lx, ly in local_points:
            all_points.append((lx * scale + tx, ly * scale + ty))

    if len(all_points) < 2:
        return None

    # Build path data
    path_data = f"M{all_points[0][0]:.1f} {all_points[0][1]:.1f}"
    for x, y in all_points[1:]:
        path_data += f" L{x:.1f} {y:.1f}"

    # Create path element
    path = ET.Element("path")
    path.set("d", path_data)
    path.set("fill", "none")
    path.set("stroke", stroke_color)
    path.set("stroke-width", str(stroke_width))
    path.set("stroke-opacity", str(stroke_opacity))
    path.set("stroke-linecap", "round")
    path.set("stroke-linejoin", "round")

    return path


def parse_svg(filepath: Path) -> PanelData:
    """
    Parse an R2DT SVG file and extract necessary data.

    Args:
        filepath: Path to SVG file

    Returns:
        PanelData with all extracted information

    Raises:
        FileNotFoundError: If file doesn't exist
        ET.ParseError: If XML is invalid
        ValueError: If required data is missing
    """
    if not filepath.exists():
        raise FileNotFoundError(f"Input file not found: {filepath}")

    tree = ET.parse(filepath)
    root = tree.getroot()

    # Parse viewBox
    viewbox = parse_viewbox(root)

    # Find anchors (5' and 3')
    # Match: 5′ (prime symbol U+2032) or 5' (apostrophe)
    anchor_5, anchor_5_elem = find_anchor(root, r"5[′']", "5′", filepath)
    anchor_3, anchor_3_elem = find_anchor(root, r"3[′']", "3′", filepath)

    # Get all children (everything except the outer svg wrapper)
    children = list(root)

    # Extract genomic coordinates for sorting and gap calculation
    genomic_start, genomic_end = extract_genomic_coords(filepath)

    # Calculate visual bounding box
    visual_bbox = calculate_visual_bbox(root, viewbox)

    # Detect nucleotide font size
    font_size = detect_nucleotide_font_size(root)

    return PanelData(
        filepath=filepath,
        root=root,
        viewbox=viewbox,
        anchor_5=anchor_5,
        anchor_3=anchor_3,
        anchor_5_elem=anchor_5_elem,
        anchor_3_elem=anchor_3_elem,
        children=children,
        sort_key=genomic_start,
        genomic_start=genomic_start,
        genomic_end=genomic_end,
        visual_bbox=visual_bbox,
        font_size=font_size,
    )


# pylint: disable-next=too-many-statements
def create_glyph(
    glyph_type: str, center_x: float, center_y: float, size: float = 12
) -> Optional[ET.Element]:
    """
    Create a join glyph element.

    Args:
        glyph_type: "bead", "bar", "break", or "none"
        center_x, center_y: Center position
        size: Size of the glyph

    Returns:
        SVG element for the glyph, or None if glyph_type is "none"
    """
    if glyph_type == "none":
        return None

    if glyph_type == "bead":
        circle = ET.Element("circle")
        circle.set("cx", str(center_x))
        circle.set("cy", str(center_y))
        circle.set("r", str(size))
        circle.set("fill", "black")
        circle.set("stroke", "black")
        circle.set("stroke-width", "2")
        return circle

    if glyph_type == "bar":
        line = ET.Element("line")
        line.set("x1", str(center_x))
        line.set("y1", str(center_y - size))
        line.set("x2", str(center_x))
        line.set("y2", str(center_y + size))
        line.set("stroke", "black")
        line.set("stroke-width", "3")
        return line

    if glyph_type == "break":
        # Break symbol: horizontal line with two diagonal slashes
        g = ET.Element("g")
        g.set("class", "break-glyph")

        stroke_width = "4"
        stroke_color = "#444"

        # Scale up the glyph
        s = size * 1.5  # Make it larger

        # Horizontal line (left segment)
        line_left = ET.SubElement(g, "line")
        line_left.set("x1", str(center_x - s * 3))
        line_left.set("y1", str(center_y))
        line_left.set("x2", str(center_x - s * 0.9))
        line_left.set("y2", str(center_y))
        line_left.set("stroke", stroke_color)
        line_left.set("stroke-width", stroke_width)

        # Horizontal line (right segment)
        line_right = ET.SubElement(g, "line")
        line_right.set("x1", str(center_x + s * 0.9))
        line_right.set("y1", str(center_y))
        line_right.set("x2", str(center_x + s * 3))
        line_right.set("y2", str(center_y))
        line_right.set("stroke", stroke_color)
        line_right.set("stroke-width", stroke_width)

        # First diagonal slash (/)
        slash1 = ET.SubElement(g, "line")
        slash1.set("x1", str(center_x - s * 0.6))
        slash1.set("y1", str(center_y + s * 0.9))
        slash1.set("x2", str(center_x - s * 0.1))
        slash1.set("y2", str(center_y - s * 0.9))
        slash1.set("stroke", stroke_color)
        slash1.set("stroke-width", stroke_width)

        # Second diagonal slash (/)
        slash2 = ET.SubElement(g, "line")
        slash2.set("x1", str(center_x + s * 0.1))
        slash2.set("y1", str(center_y + s * 0.9))
        slash2.set("x2", str(center_x + s * 0.6))
        slash2.set("y2", str(center_y - s * 0.9))
        slash2.set("stroke", stroke_color)
        slash2.set("stroke-width", stroke_width)

        return g

    return None


def create_caption_element(
    text: str, x: float, y: float, font_size: float
) -> ET.Element:
    """Create a caption text element with inline style to override CSS rules."""
    caption = ET.Element("text")
    caption.set("x", str(x))
    caption.set("y", str(y))
    caption.set("text-anchor", "middle")
    # Use inline style attribute (highest CSS priority) to override stylesheet rules
    caption.set(
        "style", f"font-size: {font_size}px; font-family: Helvetica, Arial, sans-serif;"
    )
    caption.text = text
    return caption


def remove_text_element(root: ET.Element, target_elem: ET.Element) -> bool:
    """
    Remove a text element from the SVG tree.
    Handles the case where text is wrapped in a <g> element.
    """
    # Find all text elements matching the target's position and text
    target_text = get_text_content(target_elem)

    removed = False

    # Iterate and find matching elements to remove
    for parent in root.iter():  # pylint: disable=too-many-nested-blocks
        for child in list(parent):
            tag_local = child.tag.split("}")[-1] if "}" in child.tag else child.tag

            # Check if it's a <g> element that might wrap our text
            if tag_local == "g":
                for grandchild in list(child):
                    gc_tag = (
                        grandchild.tag.split("}")[-1]
                        if "}" in grandchild.tag
                        else grandchild.tag
                    )
                    if gc_tag == "text":
                        gc_text = get_text_content(grandchild)
                        if gc_text == target_text:
                            parent.remove(child)  # Remove the whole <g>
                            removed = True
                            break

            # Check if it's a direct text element
            elif tag_local == "text":
                child_text = get_text_content(child)
                if child_text == target_text:
                    parent.remove(child)
                    removed = True

    return removed


def calculate_visual_gap(
    nt_distance: int,
    min_gap: float = 100,
    max_gap: float = 400,
    scale_factor: float = 50,
) -> float:
    """
    Calculate visual gap width based on nucleotide distance.
    Uses a logarithmic scale to compress large distances while keeping small
    distances visible.

    Args:
        nt_distance: Number of nucleotides between regions
        min_gap: Minimum visual gap width
        max_gap: Maximum visual gap width
        scale_factor: Scaling factor for logarithmic compression

    Returns:
        Visual gap width in SVG units
    """
    if nt_distance <= 0:
        return min_gap

    # Logarithmic scaling: gap = min_gap + scale_factor * log10(1 + nt_distance)
    visual_gap = min_gap + scale_factor * math.log10(1 + nt_distance)
    return min(visual_gap, max_gap)


# pylint: disable-next=too-many-branches
def strip_svg_styling(root: ET.Element) -> None:
    """
    Strip color styling from SVG elements to create a monochrome version.
    - Removes residue circles (colored backgrounds behind nucleotides)
    - Makes all text black
    - Modifies the CSS style block to override colors
    """
    elements_to_remove = []

    for elem in root.iter():
        tag_local = elem.tag.split("}")[-1] if "}" in elem.tag else elem.tag

        # Find and remove residue circles
        if tag_local == "circle":
            class_attr = elem.get("class", "")
            if "residue-circle" in class_attr or "color-posterior" in class_attr:
                # Mark for removal (can't remove while iterating)
                elements_to_remove.append(elem)

        # Make all text black
        if tag_local == "text":
            class_attr = elem.get("class", "")
            # Remove color classes and set fill to black
            if any(
                c in class_attr
                for c in ["text-green", "text-red", "text-blue", "text-brown"]
            ):
                # Replace color class with text-black
                new_class = class_attr
                for color_class in [
                    "text-green",
                    "text-red",
                    "text-blue",
                    "text-brown",
                ]:
                    new_class = new_class.replace(color_class, "text-black")
                elem.set("class", new_class)
            # Also set explicit fill attribute to ensure black
            elem.set("fill", "black")

        # Modify the CSS style block
        if tag_local == "style":
            if elem.text:
                css = elem.text
                # Override all text colors to black
                css = css.replace(
                    "fill: rgb(0, 255, 0)", "fill: rgb(0, 0, 0)"
                )  # green -> black
                css = css.replace(
                    "fill: rgb(255, 0, 255)", "fill: rgb(0, 0, 0)"
                )  # magenta -> black
                css = css.replace(
                    "fill: rgb(0, 0, 255)", "fill: rgb(0, 0, 0)"
                )  # blue -> black
                css = css.replace(
                    "fill: rgb(211.65, 104.55, 30.6)", "fill: rgb(0, 0, 0)"
                )  # brown -> black
                # Make residue circles transparent/invisible
                css = css.replace(
                    ".residue-circle { fill: rgb(255,255,255);  }",
                    ".residue-circle { fill: none; stroke: none;  }",
                )
                # Make all posterior probability colors transparent
                for i in range(11):
                    old_pattern = f".color-posterior_probability-{i} {{"
                    if old_pattern in css:
                        # Find and replace the whole rule
                        css = re.sub(
                            rf"\.color-posterior_probability-{i} \{{ fill: rgb\([^)]+\);  \}}",
                            f".color-posterior_probability-{i} {{ fill: none;  }}",
                            css,
                        )
                css = re.sub(
                    r"\.color-posterior_probability-_ \{ fill: rgb\([^)]+\);  \}",
                    ".color-posterior_probability-_ { fill: none;  }",
                    css,
                )
                elem.text = css

    # Remove marked elements
    for elem in elements_to_remove:
        for parent in root.iter():
            if elem in list(parent):
                parent.remove(elem)
                break


# pylint: disable-next=too-many-branches
def create_outline_svg(root: ET.Element, stroke_width: float = 3.0) -> None:
    """
    Transform an SVG into an outline-only version for high-level overview.

    This function:
    - Removes nucleotide text (single letters A, U, G, C, N, R, Y, etc.)
    - Keeps structural labels (5', 3', numbering)
    - Makes backbone lines much thicker with solid strokes
    - Removes colored circles/backgrounds
    - Creates a simplified silhouette view

    Args:
        root: SVG root element (modified in place)
        stroke_width: Width for backbone strokes (default 3.0)
    """
    elements_to_remove = []
    nucleotide_pattern = re.compile(
        r"^[AUGCNRYWSMKBDHVaugcnrywsmkbdhv]$|^[AUGC]-[AUGC]$"
    )

    for elem in root.iter():
        tag_local = elem.tag.split("}")[-1] if "}" in elem.tag else elem.tag

        # Remove nucleotide text elements (single letters)
        if tag_local == "text":
            text_content = elem.text.strip() if elem.text else ""
            # Keep 5', 3', numbering labels, and captions
            if nucleotide_pattern.match(text_content):
                elements_to_remove.append(elem)
                continue

        # Remove residue circles
        if tag_local == "circle":
            class_attr = elem.get("class", "")
            if "residue" in class_attr or "posterior" in class_attr:
                elements_to_remove.append(elem)

        # Make lines thicker and solid
        if tag_local == "line":
            class_attr = elem.get("class", "")
            # Skip numbering lines
            if "numbering" not in class_attr:
                elem.set("stroke-width", str(stroke_width))
                elem.set("stroke", "black")
                # Remove any dash pattern
                if elem.get("stroke-dasharray"):
                    del elem.attrib["stroke-dasharray"]

        # Make polylines thicker and solid
        if tag_local == "polyline":
            elem.set("stroke-width", str(stroke_width))
            elem.set("stroke", "black")
            if elem.get("stroke-dasharray"):
                del elem.attrib["stroke-dasharray"]

        # Modify CSS style block
        if tag_local == "style":
            if elem.text:
                css = elem.text
                # Increase stroke width for all line classes
                css = re.sub(
                    r"stroke-width:\s*[\d.]+;",
                    f"stroke-width: {stroke_width};",
                    css,
                )
                # Remove all dash patterns
                css = re.sub(
                    r"stroke-dasharray:\s*[^;]+;",
                    "",
                    css,
                )
                # Make all fills black
                css = re.sub(
                    r"fill:\s*rgb\([^)]+\)",
                    "fill: rgb(0, 0, 0)",
                    css,
                )
                # Make all strokes black
                css = re.sub(
                    r"stroke:\s*rgb\([^)]+\)",
                    "stroke: rgb(0, 0, 0)",
                    css,
                )
                # Hide residue circles
                css = re.sub(
                    r"\.residue-circle\s*\{[^}]+\}",
                    ".residue-circle { fill: none; stroke: none; }",
                    css,
                )
                elem.text = css

    # Remove marked elements
    for elem in elements_to_remove:
        for parent in root.iter():
            if elem in list(parent):
                parent.remove(elem)
                break


# pylint: disable-next=too-many-arguments,too-many-positional-arguments,too-many-locals,too-many-branches,too-many-statements
def stitch_svgs(
    panels: list[PanelData],
    gap: float = 100,
    glyph_type: str = "break",
    captions: Optional[list[str]] = None,
    caption_font_size: Optional[float] = None,
    caption_pad: float = 8,
    keep_intermediate_labels: bool = False,
    show_gap_labels: bool = True,
    gap_label_font_size: Optional[float] = None,
    monochrome: bool = False,
    outline: bool = True,
    outline_color: str = "#cccccc",
    outline_width: float = 3.0,
    outline_opacity: float = 0.6,
    anchor_label_font_size: Optional[float] = None,
    normalize_font_size: bool = False,
) -> ET.Element:
    """
    Stitch multiple SVG panels into one combined SVG.

    Args:
        panels: List of parsed panel data
        gap: Base horizontal gap between join points (used if genomic coords unavailable)
        glyph_type: Type of join glyph ("none", "bead", "bar", "break")
        captions: List of caption strings (or None for no captions)
        caption_font_size: Font size for captions (None = auto-detect from nucleotide size)
        caption_pad: Padding between panel bottom and caption
        keep_intermediate_labels: If False, remove intermediate 5'/3' labels
        show_gap_labels: If True, show nucleotide distance labels above gaps
        gap_label_font_size: Font size for gap labels (None = auto-detect from nucleotide size)
        monochrome: If True, strip colors and make all text black
        outline: If True, add a connecting outline stroke through all nucleotides
        outline_color: Color of the outline stroke
        outline_width: Width of the outline stroke
        outline_opacity: Opacity of the outline (0-1)
        anchor_label_font_size: Font size for 5'/3' labels (None = auto-detect from nucleotide size)
        normalize_font_size: If True, scale panels to match the first panel's nucleotide font size

    Returns:
        Root element of combined SVG
    """
    # Auto-detect font sizes from the first panel if not explicitly provided
    if panels:
        detected_nt_size = detect_nucleotide_font_size(panels[0].root)
    else:
        detected_nt_size = DEFAULT_NUCLEOTIDE_FONT_SIZE

    # Calculate scale factors for each panel if normalizing font size
    scale_factors = []
    if normalize_font_size and panels:
        reference_font_size = panels[0].font_size
        for panel in panels:
            if panel.font_size > 0:
                scale = reference_font_size / panel.font_size
            else:
                scale = 1.0
            scale_factors.append(scale)
    else:
        scale_factors = [1.0] * len(panels)

    # Apply defaults using multipliers if not explicitly set
    # Use a minimum caption size for readability in wide diagrams
    min_caption_font_size = 24.0
    _ = (  # Computed but used later via final_caption_font_size
        caption_font_size
        if caption_font_size is not None
        else max(min_caption_font_size, detected_nt_size * FONT_MULTIPLIERS["caption"])
    )
    effective_gap_label_font_size: float = (
        gap_label_font_size
        if gap_label_font_size is not None
        else detected_nt_size * FONT_MULTIPLIERS["gap_label"]
    )
    effective_anchor_label_font_size: float = (
        anchor_label_font_size
        if anchor_label_font_size is not None
        else detected_nt_size * FONT_MULTIPLIERS["anchor_label"]
    )

    # Apply monochrome styling to each panel if requested
    if monochrome:
        for panel in panels:
            strip_svg_styling(panel.root)

    # Normalize font sizes: add inline font-size to all text elements
    # This prevents CSS cascade issues when panels with different stylesheets are combined
    for panel in panels:
        normalize_text_font_sizes(panel.root, detected_nt_size)

    # Calculate genomic gaps between consecutive panels
    genomic_gaps = []  # List of nucleotide distances between panels
    for i in range(len(panels) - 1):
        # Gap is from end of panel i to start of panel i+1
        nt_gap = panels[i + 1].genomic_start - panels[i].genomic_end
        genomic_gaps.append(max(0, nt_gap))

    # Calculate visual gaps based on genomic distances
    visual_gaps = []
    for nt_gap in genomic_gaps:
        if nt_gap > 0:
            visual_gap = calculate_visual_gap(nt_gap)
        else:
            visual_gap = gap  # Use default gap if no genomic info
        visual_gaps.append(visual_gap)

    # Calculate translations for each panel
    translations = []  # List of (tx, ty) tuples

    # Panel 0: place so viewBox top-left maps near origin
    vb0 = panels[0].viewbox
    s0 = scale_factors[0] if scale_factors else 1.0
    t0 = (-vb0.minx * s0, -vb0.miny * s0)
    translations.append(t0)

    # Calculate positions for remaining panels using calculated visual gaps
    # and bounding box overlap prevention
    for i in range(len(panels) - 1):
        # Get scale factors
        si = scale_factors[i] if scale_factors else 1.0
        si_next = scale_factors[i + 1] if scale_factors else 1.0

        # Global position of panel i's 3' (scaled)
        ti = translations[i]
        p3_i = panels[i].anchor_3
        g3 = (ti[0] + p3_i.x * si, ti[1] + p3_i.y * si)

        # Target join point (visual_gap units to the right)
        current_gap = visual_gaps[i]
        j = (g3[0] + current_gap, g3[1])

        # Local position of panel i+1's 5' (scaled)
        p5_next = panels[i + 1].anchor_5

        # Initial translation for panel i+1
        t_next = (j[0] - p5_next.x * si_next, j[1] - p5_next.y * si_next)

        # Check for bounding box overlap and adjust if needed
        bbox_i = panels[i].visual_bbox
        bbox_next = panels[i + 1].visual_bbox

        if bbox_i is not None and bbox_next is not None:
            # Calculate global bounding box of panel i (scaled)
            global_max_x_i = ti[0] + bbox_i.max_x * si

            # Calculate global bounding box of panel i+1 with current translation (scaled)
            global_min_x_next = t_next[0] + bbox_next.min_x * si_next

            # Check for overlap
            overlap = global_max_x_i - global_min_x_next
            if overlap > -gap:  # -gap ensures minimum spacing
                # Shift panel i+1 to the right to eliminate overlap
                additional_shift = overlap + gap
                t_next = (t_next[0] + additional_shift, t_next[1])

        translations.append(t_next)

    # Create output SVG structure
    root = ET.Element(f"{{{SVG_NS}}}svg")

    # Create outline path first (so it renders behind everything)
    if outline:
        outline_path = create_outline_path(
            panels,
            translations,
            scale_factors=scale_factors if scale_factors else None,
            stroke_color=outline_color,
            stroke_width=outline_width,
            stroke_opacity=outline_opacity,
        )
        if outline_path is not None:
            outline_group = ET.SubElement(root, "g")
            outline_group.set("id", "outline")
            outline_path.set("class", "connecting-outline")
            outline_group.append(outline_path)

    # Create panel group
    panels_group = ET.SubElement(root, "g")
    panels_group.set("id", "panels")

    # Track global bounds for final viewBox
    global_min_x = float("inf")
    global_min_y = float("inf")
    global_max_x = float("-inf")
    global_max_y = float("-inf")

    # Panel bounding boxes (for caption placement)
    panel_bboxes = []  # List of (minX, maxX, maxY) in global coords

    for i, (panel, (tx, ty)) in enumerate(zip(panels, translations)):
        # Get scale factor for this panel
        si = scale_factors[i] if scale_factors else 1.0

        # Remove intermediate 5'/3' labels if requested
        if not keep_intermediate_labels:
            # Remove 3' from all but the last panel
            if i < len(panels) - 1 and panel.anchor_3_elem is not None:
                remove_text_element(panel.root, panel.anchor_3_elem)
            # Remove 5' from all but the first panel
            if i > 0 and panel.anchor_5_elem is not None:
                remove_text_element(panel.root, panel.anchor_5_elem)

        # Style kept anchor labels with larger font
        if i == 0 and panel.anchor_5_elem is not None:
            style_anchor_label(panel.anchor_5_elem, effective_anchor_label_font_size)
        if i == len(panels) - 1 and panel.anchor_3_elem is not None:
            style_anchor_label(panel.anchor_3_elem, effective_anchor_label_font_size)
        if keep_intermediate_labels:
            # Style all labels if keeping intermediates
            if panel.anchor_5_elem is not None:
                style_anchor_label(
                    panel.anchor_5_elem, effective_anchor_label_font_size
                )
            if panel.anchor_3_elem is not None:
                style_anchor_label(
                    panel.anchor_3_elem, effective_anchor_label_font_size
                )

        # Create panel group with translation and optional scaling
        panel_group = ET.SubElement(panels_group, "g")
        panel_group.set("id", f"panel_{i}")
        if si != 1.0:
            panel_group.set("transform", f"translate({tx},{ty}) scale({si})")
        else:
            panel_group.set("transform", f"translate({tx},{ty})")

        # Copy all children from source SVG (re-fetch after potential removals)
        for child in list(panel.root):
            panel_group.append(child)

        # Calculate panel bbox in global coordinates (accounting for scale)
        vb = panel.viewbox
        min_x = tx + vb.minx * si
        max_x = tx + (vb.minx + vb.width) * si
        min_y = ty + vb.miny * si
        max_y = ty + (vb.miny + vb.height) * si

        panel_bboxes.append((min_x, max_x, max_y))

        # Update global bounds
        global_min_x = min(global_min_x, min_x)
        global_min_y = min(global_min_y, min_y)
        global_max_x = max(global_max_x, max_x)
        global_max_y = max(global_max_y, max_y)

    # Create join glyphs and gap labels
    if glyph_type != "none" or show_gap_labels:
        joins_group = ET.SubElement(root, "g")
        joins_group.set("id", "joins")

        for i in range(len(panels) - 1):
            ti = translations[i]
            ti_next = translations[i + 1]
            si = scale_factors[i] if scale_factors else 1.0
            si_next = scale_factors[i + 1] if scale_factors else 1.0

            # Get global positions of 3' of panel i and 5' of panel i+1 (scaled)
            p3_i = panels[i].anchor_3
            p5_next = panels[i + 1].anchor_5
            g3 = (ti[0] + p3_i.x * si, ti[1] + p3_i.y * si)
            g5_next = (
                ti_next[0] + p5_next.x * si_next,
                ti_next[1] + p5_next.y * si_next,
            )

            # Glyph center is halfway between actual 3' and 5' positions
            gc_x = (g3[0] + g5_next[0]) / 2
            gc_y = (g3[1] + g5_next[1]) / 2

            # Create glyph
            if glyph_type != "none":
                glyph = create_glyph(glyph_type, gc_x, gc_y, size=12)
                if glyph is not None:
                    joins_group.append(glyph)

            # Add gap distance label
            if show_gap_labels and genomic_gaps[i] > 0:
                gap_label = ET.Element("text")
                gap_label.set("x", str(gc_x))
                gap_label.set("y", str(gc_y - 30))  # Above the glyph
                gap_label.set("text-anchor", "middle")
                gap_label.set("font-size", str(effective_gap_label_font_size))
                gap_label.set("font-family", "Helvetica, Arial, sans-serif")
                gap_label.set("fill", "#666")
                gap_label.text = f"{genomic_gaps[i]:,} nt"
                joins_group.append(gap_label)

                # Update global bounds to include gap label
                global_min_y = min(
                    global_min_y, gc_y - 30 - effective_gap_label_font_size
                )

    # Create captions
    if captions:
        captions_group = ET.SubElement(root, "g")
        captions_group.set("id", "captions")

        # Calculate caption font size for consistent display
        # Assume typical display: 800px wide, 400px tall container
        # Use the more restrictive constraint (fit-to-container scaling)
        current_width = global_max_x - global_min_x
        current_height = global_max_y - global_min_y

        scale_by_width = 800.0 / current_width  # scale if width-constrained
        scale_by_height = 400.0 / current_height  # scale if height-constrained
        effective_scale = min(
            scale_by_width, scale_by_height
        )  # browser uses the smaller one

        if caption_font_size is not None:
            # caption_font_size = desired screen pixels
            # Convert to SVG units: SVG_size = screen_pixels / scale
            final_caption_font_size = caption_font_size / effective_scale
        else:
            # Use nucleotide-relative sizing for consistent visual proportion
            # Caption size = 6x nucleotide size (from FONT_MULTIPLIERS["caption"])
            # This ensures captions scale properly with the RNA structure
            final_caption_font_size = detected_nt_size * FONT_MULTIPLIERS["caption"]

            # Cap caption size at 3.5% of diagram height to avoid oversized captions
            # on shorter diagrams (like dengue with only 4 panels)
            max_caption_size = current_height * 0.035
            final_caption_font_size = min(final_caption_font_size, max_caption_size)

        for caption_text, bbox in zip(captions, panel_bboxes):
            min_x, max_x, max_y = bbox

            # Caption position: centered horizontally, below panel
            caption_x = (min_x + max_x) / 2
            caption_y = max_y + caption_pad + final_caption_font_size

            caption_elem = create_caption_element(
                caption_text, caption_x, caption_y, final_caption_font_size
            )
            captions_group.append(caption_elem)

            # Update global bounds to include caption
            global_max_y = max(global_max_y, caption_y + final_caption_font_size)

    # Set final viewBox with padding around the image
    margin = 50
    final_width = global_max_x - global_min_x + 2 * margin
    final_height = global_max_y - global_min_y + 2 * margin
    root.set(
        "viewBox",
        f"{global_min_x - margin} {global_min_y - margin} {final_width} {final_height}",
    )
    root.set("width", str(final_width))
    root.set("height", str(final_height))

    return root


def write_svg(root: ET.Element, output_path: Path) -> None:
    """Write SVG element to file."""
    tree = ET.ElementTree(root)
    ET.indent(tree, space="  ")

    with open(output_path, "wb") as f:
        tree.write(f, encoding="utf-8", xml_declaration=True)
