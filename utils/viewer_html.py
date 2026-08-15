"""
Emit the ``index.html`` shell that pairs ``pdb-rna-viewer`` (2D) with
``pdbe-molstar`` (3D).

The page loads the vendored plugin + CSS, the pinned pdbe-molstar bundle,
and ``r2dt-2d-3d-viewer.js``, then bootstraps the viewer with ``R2DTViewer.create()``.
Data files (``api.json``, ``fr3d.json``, structure, ``manifest.json``) live
in the same folder so the viewer works under any same-origin HTTP server.
"""

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

# Pin pdbe-molstar to a known-good version so an upstream release can't
# silently break the viewer. Bump intentionally after testing.
_MOLSTAR_VERSION = "3.12.0"
_MOLSTAR_CDN_JS = (
    f"https://cdn.jsdelivr.net/npm/pdbe-molstar@{_MOLSTAR_VERSION}"
    "/build/pdbe-molstar-plugin.js"
)
_MOLSTAR_CDN_CSS = (
    f"https://cdn.jsdelivr.net/npm/pdbe-molstar@{_MOLSTAR_VERSION}"
    "/build/pdbe-molstar-light.css"
)

# pdb-rna-viewer's compiled bundle isn't committed to its master branch,
# isn't on npm, and the GitHub release downloads are served with
# ``application/octet-stream`` -- which browsers refuse to load as a
# stylesheet. So we vendor the build files in ``data/viewer/`` and copy
# them next to ``index.html`` at viewer-generation time.
VIEWER_PLUGIN_FILENAME = "pdb-rna-viewer-plugin-0.3.0.js"
VIEWER_CSS_FILENAME = "pdb-rna-viewer-0.3.0.css"
# Bump when r2dt-2d-3d-viewer.css / r2dt-2d-3d-viewer.js change materially (cache-bust query).
_R2DT_ASSETS_VERSION = "144"
# R2DT-owned overrides (toolbar chrome, toggles, floating buttons).
R2DT_CSS_FILENAME = "r2dt-2d-3d-viewer.css"
# The interaction glue (``R2DTViewer.create`` API).
VIEWER_JS_FILENAME = "r2dt-2d-3d-viewer.js"

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
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; margin: 0; padding: 16px; }}
  h1 {{ font-size: 16px; margin: 0 0 12px; }}
  .meta {{ color: #666; font-size: 12px; margin-bottom: 12px; }}
</style>
</head>
<body>
<h1>{structure_id} — chain {chain_id}</h1>
<div class="meta">
  2D diagram by <a href="https://r2dt.bio" target="_blank" rel="noopener">R2DT</a>
  rendered with <a href="https://github.com/PDBeurope/pdb-rna-viewer" target="_blank" rel="noopener">pdb-rna-viewer</a>,
  3D by <a href="https://github.com/molstar/pdbe-molstar" target="_blank" rel="noopener">pdbe-molstar</a>,
  base-pair annotations from {annotation_source}.
  Unresolved residues are dimmed.
</div>
<div id="r2dt-viewer-mount"></div>
<script src="{viewer_plugin_js}"></script>
<script src="{molstar_js}"></script>
<script src="{viewer_js}?v={assets_version}"></script>
<script>
R2DTViewer.create({{
  mount: '#r2dt-viewer-mount',
  baseUrl: '.',
  structureId: {structure_id_json},
  chainId: {chain_id_json},
  structureUrl: {structure_url_json},
  structureFormat: {structure_format_json},
}}).catch(function (err) {{
  console.error(err);
  var mount = document.getElementById('r2dt-viewer-mount');
  if (mount) mount.textContent = 'Failed to load viewer: ' + err.message;
}});
</script>
</body>
</html>
"""


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
    ``r2dt-2d-3d-viewer.js``) next to it. The viewer fetches the data via relative
    URLs so it works under any same-origin HTTP server without
    configuration.
    """
    structure_url = f"./{structure_filename}"
    html = _TEMPLATE.format(
        structure_id=structure_id,
        chain_id=chain_id or "",
        structure_id_json=json.dumps(structure_id),
        chain_id_json=json.dumps(chain_id or ""),
        structure_url_json=json.dumps(structure_url),
        structure_format_json=json.dumps(structure_format),
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
    stale_embed = out_dir / "embed-example.html"
    if stale_embed.is_file():
        stale_embed.unlink()
    stale_viewer_js = out_dir / "viewer.js"
    if stale_viewer_js.is_file():
        stale_viewer_js.unlink()
    stale_viewer_css = out_dir / "r2dt-viewer.css"
    if stale_viewer_css.is_file():
        stale_viewer_css.unlink()
    write_manifest(
        out_dir,
        structure_id=structure_id,
        chain_id=chain_id,
        structure_filename=structure_filename,
        structure_format=structure_format,
    )
    return target


_COMPARE_TEMPLATE = """<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>{page_title}</title>
<link rel="stylesheet" type="text/css" href="{viewer_css}">
<link rel="stylesheet" type="text/css" href="{r2dt_css}?v={assets_version}">
<link rel="stylesheet" type="text/css" href="{molstar_css}">
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
         margin: 0; padding: 16px; color: #222; }}
  h1 {{ font-size: 18px; margin: 0 0 4px; }}
  .meta {{ color: #666; font-size: 12px; margin: 0 0 12px; }}
  .mc-inf {{ display: flex; align-items: center; gap: 18px; flex-wrap: wrap;
             margin: 0 0 16px; padding: 8px 14px; background: #f8fafc;
             border: 1px solid #e5e7eb; border-radius: 8px; font-size: 13px; }}
  .mc-inf-title {{ color: #374151; font-weight: 600; }}
  .mc-inf-item {{ display: inline-flex; align-items: baseline; gap: 6px; cursor: help; }}
  .mc-inf-key {{ color: #6b7280; font-size: 11px; font-weight: 600; letter-spacing: 0.02em; }}
  .mc-inf-val {{ font-weight: 700; font-variant-numeric: tabular-nums; font-size: 15px; }}
  .mc-inf-downloads {{ position: relative; margin-left: auto; }}
  .mc-inf-download-summary {{ list-style: none; cursor: pointer; user-select: none;
                              font-size: 12px; font-weight: 600; color: #1d4ed8;
                              white-space: nowrap; }}
  .mc-inf-download-summary::-webkit-details-marker {{ display: none; }}
  .mc-inf-download-summary::marker {{ content: ""; }}
  .mc-inf-download-summary::after {{ content: " ▾"; font-size: 10px; }}
  .mc-inf-download-summary:hover {{ text-decoration: underline; }}
  .mc-inf-download-panel {{ display: none; position: absolute; right: 0; top: calc(100% + 4px);
                            z-index: 20; min-width: 9em; padding: 4px 0;
                            background: #fff; border: 1px solid #e5e7eb; border-radius: 8px;
                            box-shadow: 0 8px 20px rgba(15, 23, 42, 0.12); }}
  .mc-inf-downloads[open] > .mc-inf-download-panel {{ display: flex; flex-direction: column; }}
  .mc-inf-download {{ display: block; padding: 8px 12px; font-size: 12px; font-weight: 600;
                      color: #111827; text-decoration: none; white-space: nowrap; }}
  .mc-inf-download:hover {{ background: #f3f4f6; color: #1d4ed8; }}
  .mc-inf-scopes {{ flex: 1 0 100%; margin: 4px 0 0; }}
  .mc-inf-scopes > summary {{ cursor: pointer; color: #4b5563; font-size: 12px;
                              font-weight: 600; user-select: none; }}
  .mc-inf-scope {{ display: flex; flex-wrap: wrap; align-items: baseline; gap: 12px;
                   margin: 6px 0 0; padding: 6px 8px; background: #fff;
                   border: 1px solid #e5e7eb; border-radius: 6px; }}
  .mc-inf-scope-label {{ min-width: 9em; color: #111827; font-weight: 600; font-size: 12px; }}
</style>
</head>
<body>
<h1>{heading}</h1>
<div class="meta">{subtitle}</div>
{metrics_html}
<div id="r2dt-compare-mount"></div>
<script src="{viewer_plugin_js}"></script>
<script src="{molstar_js}"></script>
<script src="{viewer_js}?v={assets_version}"></script>
<script>
R2DTViewer.createCompare({{
  mount: '#r2dt-compare-mount',
  panels: {panels_json},
  molstar: {molstar_json},
  fetchShim: {fetch_shim_json},
}}).catch(function (err) {{
  console.error(err);
  var mount = document.getElementById('r2dt-compare-mount');
  if (mount) mount.textContent = 'Failed to load comparison viewer: ' + err.message;
}});
</script>
</body>
</html>
"""


def _inf_colour(value: Optional[float]) -> str:
    """Traffic-light colour for an INF value (green good → red poor)."""
    if value is None:
        return "#9ca3af"
    if value >= 0.95:
        return "#16a34a"
    if value >= 0.85:
        return "#65a30d"
    if value >= 0.60:
        return "#d97706"
    return "#dc2626"


def _inf_block_html(m: Dict[str, Any], label: str, desc: str) -> str:
    value = m.get("inf")
    text = "n/a" if value is None else f"{value:.3f}"
    colour = _inf_colour(value)
    title = (
        f"{desc} — TP {m.get('tp', 0)}, "
        f"FP {m.get('fp', 0)} (model-only), "
        f"FN {m.get('fn', 0)} (missing in model)"
    )
    return (
        f'<span class="mc-inf-item" title="{title}">'
        f'<span class="mc-inf-key">{label}</span>'
        f'<span class="mc-inf-val" style="color:{colour}">{text}</span></span>'
    )


def _format_scope_inf(scope: Dict[str, Any]) -> str:
    """One by-chain / inter-chain INF row (spans only — no nested divs)."""
    inf = scope.get("inf") or {}
    label = scope.get("label") or scope.get("id") or ""
    bits = [
        _inf_block_html(inf.get(k) or {}, lab, desc)
        for k, lab, desc in (
            ("wc", "WC", "canonical (cis Watson–Crick) base pairs"),
            ("nwc", "non-WC", "non-canonical base pairs"),
            ("all", "all", "all base pairs"),
        )
    ]
    return (
        f'<span class="mc-inf-scope" data-scope="{scope.get("id", "")}">'
        f'<span class="mc-inf-scope-label">{label}</span>'
        f'{"".join(bits)}</span>'
    )


def _format_inf_metrics(
    metrics: Optional[Dict[str, Any]],
    *,
    scopes: Optional[List[Dict[str, Any]]] = None,
    download_href: str = "./inf-pairs.json",
) -> str:
    """Render the INF bar (with optional by-chain scopes + download), or ''."""
    if not metrics:
        return ""
    rows = [
        ("wc", "INF WC", "canonical (cis Watson–Crick) base pairs"),
        ("nwc", "INF non-WC", "non-canonical base pairs"),
        ("all", "INF all", "all base pairs (canonical + non-canonical)"),
    ]
    items = [
        _inf_block_html(metrics.get(key) or {}, label, desc)
        for key, label, desc in rows
    ]
    download = (
        '<details class="mc-inf-downloads">'
        '<summary class="mc-inf-download-summary" '
        'title="Download INF scores and the reference/model base-pair lists">'
        "Download scores</summary>"
        '<span class="mc-inf-download-panel" role="menu">'
        f'<a class="mc-inf-download" role="menuitem" href="{download_href}" '
        'download="inf-pairs.json" title="Structured JSON report">JSON</a>'
        '<a class="mc-inf-download" role="menuitem" href="./inf-pairs.csv" '
        'download="inf-pairs.csv" title="Spreadsheet-friendly CSV">CSV</a>'
        "</span></details>"
    )
    # By-chain / inter-chain breakdown: skip the aggregate "all" entry and any
    # empty-looking scopes. Use <details>/<span> only (no nested <div>) so the
    # workstation export-menu injector can still match the outer .mc-inf.
    scope_html = ""
    subset = [s for s in (scopes or []) if s.get("type") in ("intra", "inter")]
    if subset:
        scope_rows = "".join(_format_scope_inf(s) for s in subset)
        scope_html = (
            '<details class="mc-inf-scopes">'
            "<summary>By chain</summary>"
            f"{scope_rows}"
            "</details>"
        )
    return (
        '<div class="mc-inf">'
        '<span class="mc-inf-title">Interaction Network Fidelity '
        "(base pairs, model vs reference)</span>"
        + "".join(items)
        + download
        + scope_html
        + "</div>"
    )


def render_compare(
    out_dir: Path,
    *,
    page_title: str,
    heading: str,
    subtitle: str,
    panels: List[Dict[str, Any]],
    molstar: Optional[Dict[str, Any]] = None,
    fetch_shim: bool = True,
    metrics: Optional[Dict[str, Any]] = None,
    scopes: Optional[List[Dict[str, Any]]] = None,
) -> Path:
    """Write a multi-panel ``index.html`` that calls ``R2DTViewer.createCompare()``.

    ``fetch_shim=False`` when each panel carries its data inline (``apiData`` /
    ``fr3dData``): there is nothing to fetch, so the PDB-id fetch router — which
    cannot tell apart panels that share a structure id — is disabled.
    ``metrics`` is the INF dict from ``utils.multichain.compute_inf``.
    ``scopes`` is the optional by-chain / inter-chain INF breakdown.
    """
    html = _COMPARE_TEMPLATE.format(
        page_title=page_title,
        heading=heading,
        subtitle=subtitle,
        metrics_html=_format_inf_metrics(metrics, scopes=scopes),
        panels_json=json.dumps(panels, indent=2),
        molstar_json=json.dumps(molstar or {}, indent=2),
        fetch_shim_json=json.dumps(fetch_shim),
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
    stale_compare = out_dir / "compare.js"
    if stale_compare.is_file():
        stale_compare.unlink()
    return target


def write_manifest(
    out_dir: Path,
    structure_id: str,
    chain_id: Optional[str],
    structure_filename: str,
    structure_format: str,
) -> Path:
    """Write ``manifest.json`` for ``R2DTViewer.create({ baseUrl })`` embeds."""
    manifest = {
        "structureId": structure_id,
        "chainId": chain_id or "",
        "structureUrl": f"./{structure_filename}",
        "structureFormat": structure_format,
        "assetsVersion": _R2DT_ASSETS_VERSION,
    }
    target = out_dir / "manifest.json"
    target.write_text(json.dumps(manifest, indent=2) + "\n")
    return target
