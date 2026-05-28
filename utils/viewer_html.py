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
# The interaction glue (plain JS, reads window.R2DT_CONFIG).
VIEWER_JS_FILENAME = "viewer.js"


_TEMPLATE = """<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>R2DT viewer — {structure_id}</title>
<link rel="stylesheet" type="text/css" href="{viewer_css}">
<link rel="stylesheet" type="text/css" href="{molstar_css}">
<style>
  body {{ font-family: -apple-system, sans-serif; margin: 0; padding: 16px; }}
  h1 {{ font-size: 16px; margin: 0 0 12px; }}
  .meta {{ color: #666; font-size: 12px; margin-bottom: 12px; }}
  .vis {{ display: flex; gap: 16px; align-items: flex-start; }}
  #pdb-rna-viewer {{ width: 600px; height: 600px; flex: none; }}
  #pdb-molstar {{ width: 600px; height: 600px; position: relative; flex: none; }}
</style>
</head>
<body>
<h1>{structure_id} — chain {chain_id}</h1>
<div class="meta">
  2D diagram by <a href="https://r2dt.bio" target="_blank" rel="noopener">R2DT</a>,
  3D by <a href="https://github.com/molstar/pdbe-molstar" target="_blank" rel="noopener">pdbe-molstar</a>,
  base-pair annotations from <a href="https://github.com/BGSU-RNA/fr3d-python" target="_blank" rel="noopener">FR3D (Python)</a>.
  Unresolved residues are dimmed.
</div>
<div class="vis">
  <div id="pdb-rna-viewer"></div>
  <div id="pdb-molstar"></div>
</div>

<script>
window.R2DT_CONFIG = {config_json};
</script>
<script src="{viewer_plugin_js}"></script>
<script src="{molstar_js}"></script>
<script src="{viewer_js}"></script>
</body>
</html>
"""


def render(
    out_dir: Path,
    structure_id: str,
    chain_id: Optional[str],
    structure_filename: str,
    structure_format: str,
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
        viewer_plugin_js=VIEWER_PLUGIN_FILENAME,
        viewer_css=VIEWER_CSS_FILENAME,
        viewer_js=VIEWER_JS_FILENAME,
        molstar_js=_MOLSTAR_CDN_JS,
        molstar_css=_MOLSTAR_CDN_CSS,
    )
    target = out_dir / "index.html"
    target.write_text(html)
    return target
