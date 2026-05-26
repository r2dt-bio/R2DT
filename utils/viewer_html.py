"""
Emit a self-contained ``viewer.html`` that pairs ``pdb-rna-viewer`` (2D)
with ``pdbe-molstar`` (3D) and bridges click/hover between them.

The HTML loads two local JSON files (``api.json`` and ``fr3d.json``) and
a local structure file (``<id>.cif`` / ``<id>.pdb``) -- all in the same
folder -- so the viewer renders without any network calls to RCSB/PDBe.
"""

import json
from pathlib import Path
from typing import Optional

_MOLSTAR_CDN_JS = (
    "https://cdn.jsdelivr.net/npm/pdbe-molstar@latest/build/pdbe-molstar-plugin.js"
)
_MOLSTAR_CDN_CSS = (
    "https://cdn.jsdelivr.net/npm/pdbe-molstar@latest/build/pdbe-molstar.css"
)

# pdb-rna-viewer ships its compiled bundle only under its build/ folder
# (not on a public CDN, and not committed to its master branch).
# These are the filenames the user already has from
# https://github.com/PDBeurope/pdb-rna-viewer .
VIEWER_PLUGIN_FILENAME = "pdb-rna-viewer-plugin-0.3.0.js"
VIEWER_CSS_FILENAME = "pdb-rna-viewer-0.3.0.css"


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
// pdb-rna-viewer 0.2.0 ignores options.apiData / options.FR3DData and
// always fetches from www.ebi.ac.uk.  Intercept those two URLs and
// serve the local files we wrote alongside this page.
(function() {{
  const origFetch = window.fetch.bind(window);
  window.fetch = function(input, init) {{
    const url = typeof input === 'string' ? input : (input && input.url) || '';
    if (url.includes('ebi.ac.uk/pdbe/static/entry/') && url.endsWith('_basepair.json')) {{
      return origFetch('fr3d.json', init);
    }}
    if (url.includes('ebi.ac.uk/pdbe/static/entry/') && url.endsWith('.json')) {{
      return origFetch('api.json', init);
    }}
    return origFetch(input, init);
  }};
}})();
</script>
<script src="{viewer_js}"></script>
<script src="{molstar_js}"></script>
<script>
(async function() {{
  const apiData = await (await fetch('api.json')).json();
  const fr3dData = await (await fetch('fr3d.json')).json();

  const STRUCTURE_ID = {structure_id_json};
  const CHAIN_ID = {chain_id_json};
  const STRUCTURE_URL = {structure_url_json};
  const STRUCTURE_FORMAT = {structure_format_json};

  // 1..N label  ->  PDB author residue number
  const labelToAuth = {{}};
  apiData.label_seq_ids.forEach((label, i) => {{
    if (label !== null && label !== undefined) {{
      labelToAuth[label] = apiData.auth_seq_ids[i];
    }}
  }});
  // PDB author residue number  ->  1..N label
  const authToLabel = {{}};
  Object.entries(labelToAuth).forEach(([label, auth]) => {{
    authToLabel[auth] = parseInt(label);
  }});

  // --- 2D viewer ---
  const rnaPlugin = new PdbRnaViewerPlugin();
  await rnaPlugin.render(
    document.getElementById('pdb-rna-viewer'),
    {{
      pdbId: STRUCTURE_ID.toLowerCase(),
      entityId: '1',
      chainId: CHAIN_ID || 'A',
      subscribeEvents: true,
      apiData: apiData,
      FR3DData: fr3dData,
      theme: {{ unobservedColor: '#bbbbbb' }},
    }}
  );

  // Enable every base-pair family by default (the plugin ships with
  // only cWW visible).  Poll for the "All" checkbox since it's added
  // by the plugin after render() resolves.
  (function enableAllBPs() {{
    let attempts = 0;
    const tick = () => {{
      const all = document.getElementById('Checkbox_All');
      if (all && !all.checked) {{
        all.click();
        return;
      }}
      if (attempts++ > 40) return;
      setTimeout(tick, 100);
    }};
    tick();
  }})();

  // Post-render fixups: (a) dim nucleotide letters for unobserved
  // residues -- the plugin's unobservedColor theme only colors backbone,
  // not text; (b) lighten long-range Watson-Crick pairs (pseudoknots) so
  // they don't dominate the nested cWW ladder.  The minified plugin's
  // async render() doesn't reliably wait for the DOM, so poll briefly.
  const PDB_LOWER = STRUCTURE_ID.toLowerCase();
  const unobserved = apiData.unobserved_label_seq_ids || [];
  const crossingWCPairs = new Set();
  (fr3dData.annotations || []).forEach((a) => {{
    if (a.bp === 'cWW' && a.crossing && String(a.crossing) !== '0') {{
      crossingWCPairs.add(`${{a.seq_id1}}_${{a.seq_id2}}`);
      crossingWCPairs.add(`${{a.seq_id2}}_${{a.seq_id1}}`);
    }}
  }});

  function applyFixups() {{
    let any = false;
    unobserved.forEach((seqId) => {{
      document
        .querySelectorAll(`text.rnaview_${{PDB_LOWER}}_${{seqId}}`)
        .forEach((el) => {{ el.setAttribute('fill', '#bbbbbb'); any = true; }});
    }});
    if (crossingWCPairs.size) {{
      document
        .querySelectorAll('path[class*="cWW_"]')
        .forEach((el) => {{
          const cls = el.getAttribute('class') || '';
          const m = cls.match(/cWW_(\\d+)_(\\d+)/);
          if (!m) return;
          if (!crossingWCPairs.has(`${{m[1]}}_${{m[2]}}`)) return;
          el.setAttribute('stroke', '#cccccc');
          el.setAttribute('fill', 'none');
          any = true;
        }});
    }}
    return any;
  }}
  if (unobserved.length || crossingWCPairs.size) {{
    let attempts = 0;
    const tick = () => {{
      if (applyFixups() || attempts++ > 40) return;
      setTimeout(tick, 100);
    }};
    tick();
    // Re-apply whenever the plugin re-renders nucleotides or base pairs
    // (e.g. user changes the "Filter Base Pairings" selection).
    const container = document.getElementById('pdb-rna-viewer');
    if (container && 'MutationObserver' in window) {{
      const mo = new MutationObserver(() => {{
        // Debounce: schedule a single fixup pass on next tick.
        if (mo._pending) return;
        mo._pending = true;
        setTimeout(() => {{ mo._pending = false; applyFixups(); }}, 0);
      }});
      mo.observe(container, {{ childList: true, subtree: true }});
    }}
  }}

  // --- 3D viewer ---
  const molstar = new PDBeMolstarPlugin();
  await new Promise((resolve) => {{
    molstar.render(
      document.getElementById('pdb-molstar'),
      {{
        customData: {{ url: STRUCTURE_URL, format: STRUCTURE_FORMAT, binary: false }},
        subscribeEvents: true,
        bgColor: {{ r: 255, g: 255, b: 255 }},
        hideControls: true,
        sequencePanel: false,
        loadingOverlay: true,
      }}
    );
    // pdbe-molstar fires loadComplete via events; fall back to a short delay.
    if (molstar.events && molstar.events.loadComplete) {{
      const sub = molstar.events.loadComplete.subscribe((loaded) => {{
        if (loaded) {{ sub.unsubscribe(); resolve(); }}
      }});
    }} else {{
      setTimeout(resolve, 1500);
    }}
  }});

  // Expose handles for debugging in the browser console.
  window.__r2dt = {{ molstar, rnaPlugin, apiData, fr3dData, labelToAuth, authToLabel }};

  function labelsToAuthData(labels) {{
    return labels.map(l => ({{
      auth_asym_id: CHAIN_ID,
      auth_seq_id: labelToAuth[l],
    }})).filter(d => d.auth_seq_id !== undefined && d.auth_seq_id !== null);
  }}

  // 2D -> 3D.  pdb-rna-viewer emits {{label_seq_id}} for single clicks
  // and {{label_seq_ids: [...]}} for range selections; accept both.
  function labelsFromEvent(e) {{
    const d = e.eventData || e.detail || {{}};
    if (Array.isArray(d.label_seq_ids)) return d.label_seq_ids;
    if (d.label_seq_id !== undefined && d.label_seq_id !== null) return [d.label_seq_id];
    return [];
  }}
  async function selectInMolstar(labels) {{
    const data = labelsToAuthData(labels);
    if (data.length === 0) return;
    await molstar.visual.select({{
      data: data.map(d => ({{ ...d, color: {{ r: 255, g: 112, b: 67 }}, focus: false }})),
      keepRepresentations: true,
    }});
    molstar.visual.focus(data);
  }}

  document.addEventListener('PDB.RNA.viewer.click', (e) => {{
    selectInMolstar(labelsFromEvent(e));
  }});

  // Base-pair lines in 2D don't emit PDB.RNA.viewer.click; wire our own
  // click handler that picks both partners out of the class name
  // (e.g. "cWW_5_27" -> labels [5, 27]).
  function attachBPClicks() {{
    document.querySelectorAll('path[class*="rnaviewBP"]').forEach((el) => {{
      if (el.dataset.r2dtBpClickBound) return;
      el.dataset.r2dtBpClickBound = '1';
      el.style.cursor = 'pointer';
      el.addEventListener('click', (ev) => {{
        const cls = el.getAttribute('class') || '';
        const m = cls.match(/[a-zA-Z]{{3}}_(\\d+)_(\\d+)/);
        if (!m) return;
        ev.stopPropagation();
        selectInMolstar([parseInt(m[1]), parseInt(m[2])]);
      }});
    }});
  }}
  attachBPClicks();
  // The plugin re-renders bp paths when the filter changes, so reattach.
  const bpObserver = new MutationObserver(() => {{
    if (bpObserver._pending) return;
    bpObserver._pending = true;
    setTimeout(() => {{ bpObserver._pending = false; attachBPClicks(); }}, 0);
  }});
  const bpContainer = document.getElementById('pdb-rna-viewer');
  if (bpContainer) bpObserver.observe(bpContainer, {{ childList: true, subtree: true }});
  document.addEventListener('PDB.RNA.viewer.mouseover', (e) => {{
    const labels = labelsFromEvent(e);
    // Hovering a base-pair line dispatches both partners in this event;
    // we don't want hover to flash the 3D, only the nucleotide hover.
    if (labels.length !== 1) return;
    const data = labelsToAuthData(labels);
    if (data.length === 0) return;
    molstar.visual.highlight({{ data: data, color: {{ r: 255, g: 213, b: 0 }} }});
  }});
  document.addEventListener('PDB.RNA.viewer.mouseout', () => {{
    molstar.visual.clearHighlight();
  }});

  // 3D -> 2D (already wired inside pdb-rna-viewer when subscribeEvents)
  // We just translate auth -> label in the listener that fires here.
  document.addEventListener('PDB.molstar.click', (ev) => {{
    if (!ev.eventData || ev.eventData.auth_asym_id !== CHAIN_ID) return;
    const label = authToLabel[ev.eventData.auth_seq_id];
    if (!label) return;
    document.dispatchEvent(new CustomEvent('protvista-click', {{
      detail: {{ start: label, end: label }}
    }}));
  }});
  document.addEventListener('PDB.molstar.mouseover', (ev) => {{
    if (!ev.eventData || ev.eventData.auth_asym_id !== CHAIN_ID) return;
    const label = authToLabel[ev.eventData.auth_seq_id];
    if (!label) return;
    document.dispatchEvent(new CustomEvent('protvista-mouseover', {{
      detail: {{ start: label, end: label }}
    }}));
  }});
  document.addEventListener('PDB.molstar.mouseout', () => {{
    document.dispatchEvent(new CustomEvent('protvista-mouseout', {{ detail: {{}} }}));
  }});
}})();
</script>
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
    """Write ``viewer.html`` into ``out_dir`` and return its path.

    The caller is responsible for placing ``api.json``, ``fr3d.json``,
    the structure file, and the pdb-rna-viewer assets
    (``pdb-rna-viewer.css`` and ``pdb-rna-viewer-plugin.js``) next to it.
    """
    html = _TEMPLATE.format(
        structure_id=structure_id,
        chain_id=chain_id or "",
        structure_id_json=json.dumps(structure_id),
        chain_id_json=json.dumps(chain_id or ""),
        structure_url_json=json.dumps(f"./{structure_filename}"),
        structure_format_json=json.dumps(structure_format),
        viewer_js=VIEWER_PLUGIN_FILENAME,
        viewer_css=VIEWER_CSS_FILENAME,
        molstar_js=_MOLSTAR_CDN_JS,
        molstar_css=_MOLSTAR_CDN_CSS,
    )
    target = out_dir / "index.html"
    target.write_text(html)
    return target
