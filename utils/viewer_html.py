"""
Emit the ``index.html`` shell that pairs ``pdb-rna-viewer`` (2D) with
``pdbe-molstar`` (3D).

The page is intentionally thin: it loads the vendored plugin + CSS, the
pinned pdbe-molstar bundle, and ``viewer.js`` (the interaction glue, also
vendored under ``data/viewer/`` and copied alongside). Per-structure
configuration is injected as ``window.R2DT_CONFIG``; the data files
``api.json`` / ``fr3d.json`` and the structure file are fetched from the
same folder, so the viewer works under any same-origin HTTP server
(Python ``http.server``, GitHub Pages, Cloudflare Pages, etc.).
"""

import json
from pathlib import Path
from typing import Optional

# Pin pdbe-molstar to a known-good version so an upstream release can't
# silently break the viewer. Bump intentionally after testing.
_MOLSTAR_VERSION = "3.12.0"
_MOLSTAR_CDN_JS = (
    f"https://cdn.jsdelivr.net/npm/pdbe-molstar@{_MOLSTAR_VERSION}"
    "/build/pdbe-molstar-plugin.js"
)
_MOLSTAR_CDN_CSS = (
    f"https://cdn.jsdelivr.net/npm/pdbe-molstar@{_MOLSTAR_VERSION}"
    "/build/pdbe-molstar.css"
)

# pdb-rna-viewer's compiled bundle isn't committed to its master branch,
# isn't on npm, and the GitHub release downloads are served with
# ``application/octet-stream`` -- which browsers refuse to load as a
# stylesheet. So we vendor the build files in ``data/viewer/`` and copy
# them next to ``index.html`` at viewer-generation time.
VIEWER_PLUGIN_FILENAME = "pdb-rna-viewer-plugin-0.3.0.js"
VIEWER_CSS_FILENAME = "pdb-rna-viewer-0.3.0.css"
# Bump when r2dt-viewer.css / viewer.js change materially (cache-bust query).
_R2DT_ASSETS_VERSION = "20"
# R2DT-owned overrides (toolbar chrome, toggles, floating buttons).
R2DT_CSS_FILENAME = "r2dt-viewer.css"
# The interaction glue (plain JS, reads window.R2DT_CONFIG).
VIEWER_JS_FILENAME = "viewer.js"

# Colour of the base-pair symbols, matching how pdb-rna-viewer draws them
# in the 2D diagram (light grey rather than stark black).
_BP_SYMBOL_COLOR = "#ccc"

# Credit shown in the header for where the base-pair annotations came from.
# Keyed by the --basepairs source; defaults to FR3D.
ANNOTATION_SOURCE_HTML = {
    "fr3d": (
        '<a href="https://github.com/BGSU-RNA/fr3d-python" target="_blank" '
        'rel="noopener">FR3D (Python)</a>'
    ),
    "rnaview": (
        '<a href="https://github.com/rcsb/RNAView" target="_blank" '
        'rel="noopener">RNAView</a>'
    ),
    "cif": "the mmCIF's own DNATCO/NDB base-pair annotation",
}
_DEFAULT_ANNOTATION_SOURCE = ANNOTATION_SOURCE_HTML["fr3d"]


_TEMPLATE = """<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>R2DT viewer — {structure_id}</title>
<link rel="stylesheet" type="text/css" href="{viewer_css}">
<link rel="stylesheet" type="text/css" href="{r2dt_css}?v={assets_version}">
<link rel="stylesheet" type="text/css" href="{molstar_css}">
<style>
  body {{ font-family: -apple-system, sans-serif; margin: 0; padding: 16px; }}
  h1 {{ font-size: 16px; margin: 0 0 12px; }}
  .meta {{ color: #666; font-size: 12px; margin-bottom: 12px; }}
  .vis {{ display: flex; gap: 16px; align-items: flex-start; }}
  #pdb-rna-viewer {{ width: 600px; height: 600px; flex: none; }}
  #pdb-molstar {{ width: 600px; height: 600px; position: relative; flex: none; }}
  /* --- LBN panel --- */
  #lbn-panel {{ margin-top: 20px; max-width: 1232px; }}
  #lbn-panel h2 {{ font-size: 14px; font-weight: 600; margin: 0 0 6px; }}
  .lbn-body {{
    overflow-x: auto; background: #fafafa;
    border: 1px solid #e0e0e0; border-radius: 4px; padding: 8px 12px;
  }}
  .lbn-block + .lbn-block {{ margin-top: 10px; border-top: 1px solid #eee; padding-top: 8px; }}
  .lbn-block-header {{ color: #aaa; font-size: 11px; margin-bottom: 1px; }}
  .lbn-row {{
    font-family: 'Courier New', Courier, monospace;
    font-size: 12px; white-space: pre; line-height: 1.6;
  }}
  .lbn-label {{ color: #888; }}
  .lbn-nt {{ cursor: pointer; border-radius: 2px; }}
  .lbn-nt:hover {{ background: #d0eaf8; }}
  .lbn-bp {{ cursor: pointer; font-weight: bold; border-radius: 2px; }}
  .lbn-bp:hover {{ background: #d0eaf8; }}
  .lbn-selected {{ background: #ffe066 !important; outline: 1px solid #e6b800; }}
  .lbn-partner {{ background: #ffc0cb !important; outline: 1px solid #d0607a; }}
</style>
</head>
<body>
<h1>{structure_id} — chain {chain_id}</h1>
<div class="meta">
  2D diagram by <a href="https://r2dt.bio" target="_blank" rel="noopener">R2DT</a>,
  3D by <a href="https://github.com/molstar/pdbe-molstar" target="_blank" rel="noopener">pdbe-molstar</a>,
  base-pair annotations from {annotation_source}.
  Unresolved residues are dimmed.
</div>
<div class="vis">
  <div id="pdb-rna-viewer"></div>
  <div id="pdb-molstar"></div>
</div>

{legend}

<div id="lbn-panel" style="display:none">
  <h2>Layered dot-bracket notation (Leontis–Westhof base pairs)</h2>
  <div class="lbn-body" id="lbn-body">Loading…</div>
</div>

<script>
window.R2DT_CONFIG = {config_json};
</script>
<script src="{viewer_plugin_js}"></script>
<script src="{molstar_js}"></script>
<script src="{viewer_js}?v={assets_version}"></script>
</body>
</html>
"""


def _bp_symbol(shapes, filled: bool) -> str:
    """Return an inline SVG glyph for a Leontis-Westhof base-pair symbol.

    ``shapes`` is a list of one or two edge symbols drawn on a short line:
    ``circle`` (Watson-Crick edge), ``square`` (Hoogsteen), ``tri-r`` /
    ``tri-l`` (Sugar edge, pointing right / left). ``filled`` selects the
    orientation: solid (cis) vs open (trans).
    """
    # Match compact in-panel glyph proportions (viewer.js buildBpSymbolSvg).
    width = 28 if len(shapes) == 1 else 40
    height, cy, r, stroke_w = 16, 8, 4, 1.2
    color = _BP_SYMBOL_COLOR
    fill = color if filled else "#fff"
    parts = [
        f'<line x1="2" y1="{cy}" x2="{width - 2}" y2="{cy}" '
        f'stroke="{color}" stroke-width="{stroke_w}"/>'
    ]
    positions = [width / 2] if len(shapes) == 1 else [14, 26]
    for kind, cx in zip(shapes, positions):
        if kind == "circle":
            parts.append(
                f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="{fill}" '
                f'stroke="{color}" stroke-width="{stroke_w}"/>'
            )
        elif kind == "square":
            parts.append(
                f'<rect x="{cx - r}" y="{cy - r}" width="{2 * r}" height="{2 * r}" '
                f'fill="{fill}" stroke="{color}" stroke-width="{stroke_w}"/>'
            )
        elif kind == "tri-r":
            parts.append(
                f'<polygon points="{cx - r},{cy - r} {cx - r},{cy + r} {cx + r},{cy}" '
                f'fill="{fill}" stroke="{color}" stroke-width="{stroke_w}"/>'
            )
        elif kind == "tri-l":
            parts.append(
                f'<polygon points="{cx + r},{cy - r} {cx + r},{cy + r} {cx - r},{cy}" '
                f'fill="{fill}" stroke="{color}" stroke-width="{stroke_w}"/>'
            )
    return (
        f'<svg class="r2dt-bp-glyph" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">{"".join(parts)}</svg>'
    )


def _bp_legend_html() -> str:
    """Build the Leontis-Westhof legend (mounted in-panel by viewer.js)."""
    # Each row: cis label, cis glyph specs, trans label, trans glyph specs.
    # A spec is a list of one or two edge symbols drawn on a short line:
    # ``circle`` (Watson-Crick edge), ``square`` (Hoogsteen), ``tri-r`` /
    # ``tri-l`` (Sugar edge, pointing right / left). ``filled`` selects the
    # orientation: solid (cis) vs open (trans).
    rows = [
        ("cWW", [["circle"]], "tWW", [["circle"]]),
        (
            "cWH, cHW",
            [["circle", "square"], ["square", "circle"]],
            "tWH, tHW",
            [["circle", "square"], ["square", "circle"]],
        ),
        (
            "cWS, cSW",
            [["circle", "tri-r"], ["tri-l", "circle"]],
            "tWS, tSW",
            [["circle", "tri-r"], ["tri-l", "circle"]],
        ),
        ("cHH", [["square"]], "tHH", [["square"]]),
        (
            "cHS, cSH",
            [["square", "tri-r"], ["tri-l", "square"]],
            "tHS, tSH",
            [["square", "tri-r"], ["tri-l", "square"]],
        ),
        ("cSS", [["tri-r"]], "tSS", [["tri-r"]]),
    ]
    body = []
    for cis_label, cis_glyphs, trans_label, trans_glyphs in rows:
        cis_svg = "".join(_bp_symbol(s, True) for s in cis_glyphs)
        trans_svg = "".join(_bp_symbol(s, False) for s in trans_glyphs)
        body.append(
            f"<tr><th>{cis_label}</th><td>{cis_svg}</td>"
            f"<th>{trans_label}</th><td>{trans_svg}</td></tr>"
        )
    rows_html = "\n".join(body)
    return f"""<div id="r2dt-bp-legend-source" hidden>
<div class="r2dt-bp-legend-panel">
<p class="r2dt-bp-legend-intro">Each glyph encodes the two interacting edges
and glycosidic-bond orientation. Edge shape: circle = Watson–Crick,
square = Hoogsteen, triangle = Sugar. Fill: solid = <em>cis</em>,
open = <em>trans</em>.</p>
<table class="r2dt-bp-legend-table">
<thead><tr><th colspan="2">cis</th><th colspan="2">trans</th></tr></thead>
<tbody>
{rows_html}
</tbody>
</table>
<p class="r2dt-bp-legend-cite">Symbols follow
<a href="https://pubmed.ncbi.nlm.nih.gov/11345429/" target="_blank" rel="noopener">Leontis &amp; Westhof, 2001</a>.</p>
</div>
</div>"""


_LEGEND_HTML = _bp_legend_html()


# pylint: disable=too-many-arguments,too-many-positional-arguments
def render(
    out_dir: Path,
    structure_id: str,
    chain_id: Optional[str],
    structure_filename: str,
    structure_format: str,
    annotation_source: Optional[str] = None,
) -> Path:
    """Write ``index.html`` into ``out_dir`` and return its path.

    The caller is responsible for placing ``api.json``, ``fr3d.json``,
    the structure file, and the vendored viewer assets (plugin JS/CSS and
    ``viewer.js``) next to it. The viewer fetches the data via relative
    URLs so it works under any same-origin HTTP server without
    configuration.
    """
    config = {
        "structureId": structure_id,
        "chainId": chain_id or "",
        "structureUrl": f"./{structure_filename}",
        "structureFormat": structure_format,
    }
    # Escape "</" so the JSON can't prematurely close the <script> block.
    config_json = json.dumps(config).replace("</", "<\\/")
    html = _TEMPLATE.format(
        structure_id=structure_id,
        chain_id=chain_id or "",
        config_json=config_json,
        legend=_LEGEND_HTML,
        annotation_source=annotation_source or _DEFAULT_ANNOTATION_SOURCE,
        viewer_plugin_js=VIEWER_PLUGIN_FILENAME,
        viewer_css=VIEWER_CSS_FILENAME,
        r2dt_css=R2DT_CSS_FILENAME,
        viewer_js=VIEWER_JS_FILENAME,
        assets_version=_R2DT_ASSETS_VERSION,
        molstar_js=_MOLSTAR_CDN_JS,
        molstar_css=_MOLSTAR_CDN_CSS,
    )
    target = out_dir / "index.html"
    target.write_text(html)
    return target
