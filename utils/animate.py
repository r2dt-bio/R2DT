"""Generate an SVG animation."""

import argparse
import re
from math import ceil, floor

DURATION = 2

SVG_HEADER = """
<svg
	xmlns="http://www.w3.org/2000/svg"
    viewBox="MIN_WIDTH MIN_HEIGHT WIDTH HEIGHT"
    width="WIDTH"
    height="HEIGHT">
<g><style type="text/css">
<![CDATA[
circle.red {stroke: rgb(255, 0, 255); fill: none; }
circle.green {stroke: rgb(0, 255, 0); fill: none; }
circle.blue {stroke: rgb(0, 0, 255); fill: none; }
circle.black {stroke: rgb(0, 0, 0); fill: none; }
circle.gray {stroke: rgb(204, 204, 204); fill: none; }
circle.brown {stroke: rgb(211.65, 104.55, 30.6); fill: none; }
line.red {stroke: rgb(255, 0, 255); stroke-width: 0.382365; }
line.green {stroke: rgb(0, 255, 0); stroke-width: 0.382365; }
line.blue {stroke: rgb(0, 0, 255); stroke-width: 0.382365; }
line.black {stroke: rgb(0, 0, 0); stroke-width: 0.382365; }
line.gray {stroke: rgb(204, 204, 204); stroke-width: 0.382365; }
line.brown {stroke: rgb(211.65, 104.55, 30.6); stroke-width: 0.382365; }
line {stroke: rgb(0, 0, 0); }
text.red {fill: rgb(255, 0, 255); font-size: FONT_SIZEpx; font-weight: bold; font-family: Helvetica; }
text.green {fill: rgb(0, 255, 0); font-size: FONT_SIZEpx; font-weight: bold; font-family: Helvetica; }
text.blue {fill: rgb(0, 0, 255); font-size: FONT_SIZEpx; font-weight: bold; font-family: Helvetica; }
text.black {fill: rgb(0, 0, 0); font-size: FONT_SIZEpx; font-weight: bold; font-family: Helvetica; }
text.gray {fill: rgb(204, 204, 204); font-size: FONT_SIZEpx; font-weight: bold; font-family: Helvetica; }
text.brown {fill: rgb(211.65, 104.55, 30.6); font-size: FONT_SIZEpx; font-weight: bold; font-family: Helvetica; }
text {fill: rgb(0, 0, 0); font-size: FONT_SIZEpx; font-weight: bold; font-family: Helvetica; text-anchor: middle; alignment-baseline: middle; }
.pseudoknot_connection{stroke-linecap: round; stroke-opacity: 0.2; stroke-width: 1.5; }
.pseudoknot_segment1{stroke-linecap: round; stroke-opacity: 0.4; stroke-width: 3.05892; }
.pseudoknot_segment2{stroke-linecap: round; stroke-opacity: 0.4; stroke-width: 3.05892; }
line.numbering-line{stroke: rgb(204, 204, 204); stroke-width: 0.191182; }
line.predicted{stroke: rgb(0, 0, 0); stroke-dasharray: 2; }
polyline{fill: none; stroke-linejoin: round; }
template{visibility: hidden; }
text.numbering-label{fill: rgb(204, 204, 204); }
]]>
</style></g>
"""


# pylint: disable=too-many-statements,too-many-branches
def move_to_same_start(lines1, lines2):
    """Move the second SVG to the same start position as the first SVG."""
    for line in lines1:
        if "5'" in line:
            match = re.search(r'x="(-?\d+(\.\d+)?)" y="(-?\d+(\.\d+)?)"', line)
            if match:
                x1_coord, _, y1_coord, _ = match.groups()
                break
    for line in lines2:
        if "5'" in line:
            match = re.search(r'x="(-?\d+(\.\d+)?)" y="(-?\d+(\.\d+)?)"', line)
            if match:
                x2_coord, _, y2_coord, _ = match.groups()
                break
    delta_x = float(x1_coord) - float(x2_coord)
    delta_y = float(y1_coord) - float(y2_coord)
    for i, line in enumerate(lines2):
        if 'x="' in line:
            match = re.search(r'x="(-?\d+(\.\d+)?)"', line)
            if match:
                x_coord = match.group(1)
                x_coord = str(float(x_coord) + delta_x)
                line = line.replace('x="' + match.group(1), 'x="' + x_coord)
                lines2[i] = line
        if 'y="' in line:
            match = re.search(r'y="(-?\d+(\.\d+)?)"', line)
            if match:
                y_coord = match.group(1)
                y_coord = str(float(y_coord) + delta_y)
                line = line.replace('y="' + match.group(1), 'y="' + y_coord)
                lines2[i] = line
        if 'x1="' in line:
            match = re.search(r'x1="(-?\d+(\.\d+)?)"', line)
            if match:
                x_coord = match.group(1)
                x_coord = str(float(x_coord) + delta_x)
                line = line.replace('x1="' + match.group(1), 'x1="' + x_coord)
                lines2[i] = line
            match = re.search(r'x2="(-?\d+(\.\d+)?)"', line)
            if match:
                x_coord = match.group(1)
                x_coord = str(float(x_coord) + delta_x)
                line = line.replace('x2="' + match.group(1), 'x2="' + x_coord)
                lines2[i] = line
            match = re.search(r'y1="(-?\d+(\.\d+)?)"', line)
            if match:
                y_coord = match.group(1)
                y_coord = str(float(y_coord) + delta_y)
                line = line.replace('y1="' + match.group(1), 'y1="' + y_coord)
                lines2[i] = line
            match = re.search(r'y2="(-?\d+(\.\d+)?)"', line)
            if match:
                y_coord = match.group(1)
                y_coord = str(float(y_coord) + delta_y)
                line = line.replace('y2="' + match.group(1), 'y2="' + y_coord)
                lines2[i] = line
    return lines2


def add_transform(svg):
    """Add a transform to the SVG to move the text."""

    max_height, min_height, max_width, min_width = 0, 0, 0, 0
    structure1, structure2 = [], []

    for line in svg:
        if "green" in line or "numbering-label" in line:
            continue  # skip basepair lines and basepair numbers
        if "<text" in line and "<g" in line and "second" not in line:
            structure1.append(line)
        if "<text" in line and "<g" in line and "second" in line:
            structure2.append(line)
    print(len(structure1), len(structure2))
    animated_lines = []
    for index, line in enumerate(structure1):
        x1_coord, _, y1_coord, _ = re.search(
            r'x="(-?\d+(\.\d+)?)" y="(-?\d+(\.\d+)?)"', line
        ).groups()
        try:
            x2_coord, _, y2_coord, _ = re.search(
                r'x="(-?\d+(\.\d+)?)" y="(-?\d+(\.\d+)?)"', structure2[index]
            ).groups()
        except IndexError:
            continue

        x1_coord, y1_coord, x2_coord, y2_coord = (
            float(x1_coord),
            float(y1_coord),
            float(x2_coord),
            float(y2_coord),
        )

        delta_x = x2_coord - x1_coord
        delta_y = y2_coord - y1_coord
        transform_line = (
            f'<animateTransform id="anim{index}" '
            f'attributeName="transform" type="translate" from="0 0" '
            f'to="{delta_x} {delta_y}" restart="always" repeatCount="1" '
            f'dur="{DURATION}s" fill="freeze" />'
        )
        line = line.replace("</text>", f"{transform_line}</text>")
        animated_lines.append(line)

        max_height = max(max_height, y2_coord + delta_y, y1_coord)
        max_width = max(max_width, x2_coord + delta_x, x1_coord)
        min_height = min(min_height, y2_coord + delta_y, y1_coord)
        min_width = min(min_width, x2_coord + delta_x, x1_coord)

    return animated_lines, max_height, max_width, min_height, min_width


def read_in_svg(filename):
    """Read in an SVG file."""
    with open(filename, "r", encoding="utf-8") as f_svg:
        return f_svg.readlines()


def generate_combined_svg(svg1, svg2):
    """Generate the combined SVG."""
    combined_svg = svg1[:-1]  # first SVG without closing svg tag
    for line in svg2[1:]:  # second SVG without opening svg tag
        combined_svg.append(
            line.replace('class="black"', 'class="black second"').replace(
                'class="gray"', 'class="gray second"'
            )
        )
    return combined_svg


def get_lines(svg, include_class=None, add_class=None):
    """Get line tags from SVG."""
    lines = []
    for line in svg:
        if "<line" in line:
            if add_class:
                line = line.replace('class="', f' class="{add_class} ')
            if include_class is None:
                lines.append(line)
            elif include_class in line:
                lines.append(line)
    return lines


def get_font_size(svg1, svg2):
    """Get maximum font size of the SVGs."""
    font_size1 = 0
    font_size2 = 0
    for line in svg1:
        if "text.black" in line and "font-size" in line:
            font_size1 = float(re.search(r"font-size: (\d+(\.\d+)?)px;", line).group(1))
            break
    for line in svg2:
        if "text.black" in line and "font-size" in line:
            font_size2 = float(re.search(r"font-size: (\d+(\.\d+)?)px;", line).group(1))
            break
    return max(font_size1, font_size2)


def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="Generate an SVG animation.")
    parser.add_argument("svg1", help="First SVG file")
    parser.add_argument("svg2", help="Second SVG file")
    parser.add_argument("animated", help="Animated SVG file")
    args = parser.parse_args()

    svg1 = read_in_svg(args.svg1)
    svg2 = read_in_svg(args.svg2)
    lines_moved = move_to_same_start(svg1, svg2)

    moved_svg_lines = generate_combined_svg(svg1, lines_moved)
    animated_text, max_height, max_width, min_height, min_width = add_transform(
        moved_svg_lines
    )

    svg1_lines = get_lines(svg1, include_class=None, add_class="first")
    svg2_moved_lines = get_lines(
        moved_svg_lines, include_class="gray second", add_class=None
    ) + get_lines(moved_svg_lines, include_class="black second", add_class=None)

    height = max_height - min_height
    width = max_width - min_width

    height = ceil(height) + 5
    width = ceil(width) + 5
    min_height = floor(min_height)
    min_width = floor(min_width)

    font_size = get_font_size(svg1, svg2)

    with open(args.animated, "w", encoding="utf-8") as f_svg:
        f_svg.write(
            SVG_HEADER.replace("MIN_HEIGHT", str(min_height))
            .replace("MIN_WIDTH", str(min_width))
            .replace("WIDTH", str(width))
            .replace("HEIGHT", str(height))
            .replace("FONT_SIZE", str(font_size))
        )
        f_svg.write("".join(animated_text))
        f_svg.write("".join(svg1_lines))
        f_svg.write("".join(svg2_moved_lines))
        f_svg.write("</svg>")


if __name__ == "__main__":
    main()
