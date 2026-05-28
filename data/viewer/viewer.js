/*
 * R2DT 2D+3D viewer glue.
 *
 * Pairs pdb-rna-viewer (2D) with pdbe-molstar (3D) and bridges
 * click selection between them. Configuration (structure id, chain,
 * structure file URL/format) is injected by the generated index.html as
 * ``window.R2DT_CONFIG``; the data files ``api.json`` and ``fr3d.json``
 * plus the structure file are fetched from the same folder.
 *
 * This file is vendored under ``data/viewer/`` and copied next to each
 * generated ``index.html`` by the ``pdb_2d_3d`` command.
 */

// pdb-rna-viewer 0.3.0 ignores options.apiData / options.FR3DData and
// always issues fetches for these EBI URLs. Redirect them to the local
// files written alongside this page so the viewer works same-origin
// (HTTP server, GitHub Pages, Cloudflare Pages, etc.). Installed before
// rnaPlugin.render() runs below.
(function () {
  const origFetch = window.fetch.bind(window);
  window.fetch = function (input, init) {
    const url = typeof input === 'string' ? input : (input && input.url) || '';
    if (url.includes('ebi.ac.uk/pdbe/static/entry/') && url.endsWith('_basepair.json')) {
      return origFetch('fr3d.json', init);
    }
    if (url.includes('ebi.ac.uk/pdbe/static/entry/') && url.endsWith('.json')) {
      return origFetch('api.json', init);
    }
    return origFetch(input, init);
  };
})();

(async function () {
  const CONFIG = window.R2DT_CONFIG || {};
  const STRUCTURE_ID = CONFIG.structureId;
  const CHAIN_ID = CONFIG.chainId || '';
  const STRUCTURE_URL = CONFIG.structureUrl;
  const STRUCTURE_FORMAT = CONFIG.structureFormat;

  const apiData = await (await fetch('api.json')).json();
  const fr3dData = await (await fetch('fr3d.json')).json();

  // The plugin auto-zooms the 2D view to the clicked nucleotide. That
  // creates a jarring "zoom out + zoom back in" jolt because the focus
  // change is animated -- disable it so clicks just highlight.
  if (window.UiActionsService && window.UiActionsService.zoomToNucleotides) {
    window.UiActionsService.zoomToNucleotides = function () {};
  }

  // 1..N label  ->  PDB author residue number
  const labelToAuth = {};
  apiData.label_seq_ids.forEach((label, i) => {
    if (label !== null && label !== undefined) {
      labelToAuth[label] = apiData.auth_seq_ids[i];
    }
  });
  // PDB author residue number  ->  1..N label
  const authToLabel = {};
  Object.entries(labelToAuth).forEach(([label, auth]) => {
    authToLabel[auth] = parseInt(label);
  });

  // --- 2D viewer ---
  const rnaPlugin = new PdbRnaViewerPlugin();
  await rnaPlugin.render(
    document.getElementById('pdb-rna-viewer'),
    {
      pdbId: STRUCTURE_ID.toLowerCase(),
      entityId: '1',
      chainId: CHAIN_ID || 'A',
      subscribeEvents: true,
      apiData: apiData,
      FR3DData: fr3dData,
      theme: { unobservedColor: '#bbbbbb' },
    }
  );

  // Render a faint backbone path overlay UNDER the nucleotide letters,
  // always on by default. Repurposes the plugin's "View as Path" dropdown
  // into a Show/Hide toggle for that overlay (the plugin's letter view
  // remains the only main view -- we no longer expose path-only mode).
  let backboneVisible = true;
  // Derive nucleotide centre coordinates from apiData.svg_paths.
  // Each entry past the two dummy heads encodes prevX,prevY,x,y after
  // splitting on "M" and ","; the "current" point sits at tokens[2..3].
  function buildBackbonePoints() {
    const pts = [];
    const paths = apiData.svg_paths || [];
    for (let i = 2; i < paths.length; i++) {
      const tokens = paths[i].split(/[M,]/).filter((s) => s !== '');
      if (tokens.length < 4) continue;
      const x = parseFloat(tokens[2]);
      const y = parseFloat(tokens[3]);
      if (!isNaN(x) && !isNaN(y)) pts.push(x + ',' + y);
    }
    return pts.join(' ');
  }
  function injectBackboneOverlay() {
    const svg = document.querySelector('svg.rnaTopoSvg');
    if (!svg) return false;
    const existing = svg.querySelector('.r2dt-backbone-overlay');
    if (!backboneVisible) {
      if (existing) existing.remove();
      return true;
    }
    // Idempotent: if the overlay is already present, do nothing. This
    // matters because the re-inject runs from a MutationObserver -- adding
    // the polyline is itself a mutation, so without this guard we'd loop.
    if (existing) return true;
    const points = buildBackbonePoints();
    if (!points) return false;
    const inner = svg.querySelector('[class^="rnaTopoSvg_"]');
    if (!inner) return false;
    const polyline = document.createElementNS('http://www.w3.org/2000/svg', 'polyline');
    polyline.setAttribute('class', 'r2dt-backbone-overlay');
    polyline.setAttribute('points', points);
    polyline.setAttribute('fill', 'none');
    // Scale the overlay with structure size: small RNAs (~100 nt) need a
    // thicker, less transparent line to be legible against the letters;
    // large rRNAs (1000+ nt) look noisy at the same settings.
    const nPts = points.split(' ').length;
    const isSmall = nPts < 500;
    polyline.setAttribute('stroke', '#666');
    polyline.setAttribute('stroke-width', isSmall ? '2' : '1.2');
    polyline.setAttribute('opacity', isSmall ? '0.45' : '0.25');
    polyline.style.pointerEvents = 'none';
    // Insert as the first child of the plugin's inner group so the same
    // zoom/pan transform applies to the overlay as to the nucleotides.
    inner.insertBefore(polyline, inner.firstChild);
    return true;
  }
  (function setupBackboneToggle() {
    let attempts = 0;
    const tick = () => {
      const sel = document.querySelector('.menuSelectbox');
      if (!sel) {
        if (attempts++ > 40) return;
        setTimeout(tick, 100);
        return;
      }
      // Hide the plugin's "View as Path" dropdown entirely. Touching its
      // value/listeners breaks the BP filter (pathOrNucleotide reads the
      // dropdown value as an int when redrawing BP lines). Replace it
      // with our own checkbox that only governs the backbone overlay.
      sel.style.display = 'none';
      const toggle = document.createElement('label');
      toggle.style.cssText = 'font-size: 12px; cursor: pointer; user-select: none;';
      toggle.innerHTML =
        '<input type="checkbox" id="r2dt-backbone-toggle" checked> Show backbone path';
      sel.parentNode.insertBefore(toggle, sel.nextSibling);
      toggle.querySelector('input').addEventListener('change', (ev) => {
        backboneVisible = ev.target.checked;
        injectBackboneOverlay();
      });
      injectBackboneOverlay();
    };
    tick();
  })();

  // Hide base-pair family checkboxes that have no annotations in this
  // structure's fr3d.json, so the dropdown only lists families the user
  // can actually toggle. Then enable the remaining families ("All").
  (function tidyBPFilter() {
    const presentBPs = new Set(
      (fr3dData.annotations || []).map((a) => a.bp).filter(Boolean)
    );
    let attempts = 0;
    const tick = () => {
      const boxes = document.querySelectorAll('input[id^="Checkbox_"]');
      if (boxes.length === 0) {
        if (attempts++ > 40) return;
        setTimeout(tick, 100);
        return;
      }
      boxes.forEach((cb) => {
        if (cb.id === 'Checkbox_All') return;
        const family = cb.id.slice('Checkbox_'.length);
        if (!presentBPs.has(family)) {
          const row = cb.parentElement;
          if (row) row.style.display = 'none';
        }
      });
      const all = document.getElementById('Checkbox_All');
      if (all && !all.checked) all.click();
    };
    tick();
  })();

  // Rename the plugin's "Filter Base Pairings" button to "Filter Base
  // Pairs", preserving its help icon and caret (only the text node).
  (function renameBpFilterBtn() {
    let attempts = 0;
    const tick = () => {
      const btn = document.getElementById('bpFilterBtn');
      if (!btn) {
        if (attempts++ > 40) return;
        setTimeout(tick, 100);
        return;
      }
      btn.childNodes.forEach((node) => {
        if (node.nodeType === 3 && node.nodeValue.includes('Filter Base Pairings')) {
          node.nodeValue = node.nodeValue.replace('Filter Base Pairings', 'Filter Base Pairs');
        }
      });
    };
    tick();
  })();

  // Declared up front so the post-render fixup pass (which skips the
  // currently-selected BP) can reference it before the click handler
  // block initializes -- temporal-dead-zone errors otherwise.
  let lastBPSelected = null;

  // Post-render fixups: (a) dim nucleotide letters for unobserved
  // residues -- the plugin's unobservedColor theme only colors backbone,
  // not text; (b) lighten long-range Watson-Crick pairs (pseudoknots) so
  // they don't dominate the nested cWW ladder. The minified plugin's
  // async render() doesn't reliably wait for the DOM, so poll briefly.
  const PDB_LOWER = STRUCTURE_ID.toLowerCase();
  const unobserved = apiData.unobserved_label_seq_ids || [];
  const crossingWCPairs = new Set();
  (fr3dData.annotations || []).forEach((a) => {
    if (a.bp === 'cWW' && a.crossing && String(a.crossing) !== '0') {
      crossingWCPairs.add(`${a.seq_id1}_${a.seq_id2}`);
      crossingWCPairs.add(`${a.seq_id2}_${a.seq_id1}`);
    }
  });

  function applyFixups() {
    let any = false;
    unobserved.forEach((seqId) => {
      document
        .querySelectorAll(`text.rnaview_${PDB_LOWER}_${seqId}`)
        .forEach((el) => { el.setAttribute('fill', '#bbbbbb'); any = true; });
    });
    if (crossingWCPairs.size) {
      document
        .querySelectorAll('path[class*="cWW_"]')
        .forEach((el) => {
          // Don't recolour the currently-selected BP; the click handler
          // already painted it orange and would otherwise be overwritten.
          if (el === lastBPSelected) return;
          const cls = el.getAttribute('class') || '';
          const m = cls.match(/cWW_(\d+)_(\d+)/);
          if (!m) return;
          if (!crossingWCPairs.has(`${m[1]}_${m[2]}`)) return;
          el.setAttribute('stroke', '#cccccc');
          el.setAttribute('fill', 'none');
          // Remember the greyed-out value so the click handler's
          // restore-on-next-click brings it back to grey, not the
          // plugin's default colour.
          el.dataset.r2dtOrigStroke = '#cccccc';
          any = true;
        });
    }
    return any;
  }
  // Re-apply the dimming fixups and re-inject the backbone overlay
  // whenever the plugin re-renders the SVG -- e.g. when the BP-family
  // filter changes, or right after tidyBPFilter clicks "All" on load,
  // which rewrites the inner group and wipes the overlay. Always active
  // (not gated on there being something to dim) so the backbone is
  // restored even for structures with no unobserved/crossing residues.
  if (unobserved.length || crossingWCPairs.size) {
    let attempts = 0;
    const tick = () => {
      if (applyFixups() || attempts++ > 40) return;
      setTimeout(tick, 100);
    };
    tick();
  }
  const container = document.getElementById('pdb-rna-viewer');
  if (container && 'MutationObserver' in window) {
    const mo = new MutationObserver(() => {
      // Debounce: schedule a single pass on next tick.
      if (mo._pending) return;
      mo._pending = true;
      setTimeout(() => {
        mo._pending = false;
        applyFixups();
        injectBackboneOverlay();
      }, 0);
    });
    mo.observe(container, { childList: true, subtree: true });
  }

  // --- 3D viewer ---
  const molstar = new PDBeMolstarPlugin();
  await new Promise((resolve) => {
    molstar.render(
      document.getElementById('pdb-molstar'),
      {
        customData: { url: STRUCTURE_URL, format: STRUCTURE_FORMAT, binary: false },
        subscribeEvents: true,
        bgColor: { r: 255, g: 255, b: 255 },
        hideControls: true,
        sequencePanel: false,
        loadingOverlay: true,
      }
    );
    // pdbe-molstar fires loadComplete via events; fall back to a short delay.
    if (molstar.events && molstar.events.loadComplete) {
      const sub = molstar.events.loadComplete.subscribe((loaded) => {
        if (loaded) { sub.unsubscribe(); resolve(); }
      });
    } else {
      setTimeout(resolve, 1500);
    }
  });

  // Expose handles for debugging in the browser console.
  window.__r2dt = { molstar, rnaPlugin, apiData, fr3dData, labelToAuth, authToLabel };

  function labelsToAuthData(labels) {
    // pdbe-molstar's visual.select / visual.highlight expect a residue
    // range via start_auth_residue_number / end_auth_residue_number.
    // Passing auth_seq_id alone is silently ignored, which selects the
    // whole chain and makes focus() zoom to nothing.
    return labels.map((l) => {
      const auth = labelToAuth[l];
      if (auth === undefined || auth === null) return null;
      return {
        auth_asym_id: CHAIN_ID,
        start_auth_residue_number: auth,
        end_auth_residue_number: auth,
      };
    }).filter((d) => d !== null);
  }

  // 2D -> 3D. pdb-rna-viewer emits {label_seq_id} for single clicks and
  // {label_seq_ids: [...]} for range selections; accept both.
  function labelsFromEvent(e) {
    const d = e.eventData || e.detail || {};
    if (Array.isArray(d.label_seq_ids)) return d.label_seq_ids;
    if (d.label_seq_id !== undefined && d.label_seq_id !== null) return [d.label_seq_id];
    return [];
  }
  async function selectInMolstar(labels) {
    const data = labelsToAuthData(labels);
    if (data.length === 0) return;
    await molstar.visual.select({
      data: data.map((d) => ({ ...d, color: { r: 255, g: 112, b: 67 }, focus: false })),
      keepRepresentations: true,
    });
    molstar.visual.focus(data);
  }

  document.addEventListener('PDB.RNA.viewer.click', (e) => {
    selectInMolstar(labelsFromEvent(e));
  });

  // Select a base pair: colour its 2D line orange (if its path is in the
  // DOM) and select both partner residues in 3D. `pathEl` may be null --
  // e.g. when triggered from the base-pair list while that pair's line
  // isn't currently rendered.
  const BP_SELECTED_COLOR = 'orange'; // matches the plugin's nucleotide click colour
  function selectBasePair(a, b, pathEl) {
    if (pathEl) {
      // Restore the previously selected BP's stroke, then mark this one so
      // only one base pair appears highlighted at a time.
      if (lastBPSelected && lastBPSelected !== pathEl) {
        const prevOrig = lastBPSelected.dataset.r2dtOrigStroke;
        if (prevOrig !== undefined) lastBPSelected.setAttribute('stroke', prevOrig);
      }
      if (pathEl.dataset.r2dtOrigStroke === undefined) {
        pathEl.dataset.r2dtOrigStroke = pathEl.getAttribute('stroke') || '';
      }
      pathEl.setAttribute('stroke', BP_SELECTED_COLOR);
      lastBPSelected = pathEl;
    }
    selectInMolstar([a, b]);
  }

  // Find a rendered base-pair path by its two seq ids (either order).
  // The bp class is the last token, "<bp>_<a>_<b>", so an attribute
  // "ends-with" match avoids false hits like _5_27 inside _5_271.
  function findBPPath(a, b) {
    return document.querySelector(
      `.rnaviewBP[class$="_${a}_${b}"], .rnaviewBP[class$="_${b}_${a}"]`
    );
  }

  // Base-pair lines in 2D don't emit PDB.RNA.viewer.click; wire our own
  // click handler that picks both partners out of the class name
  // (e.g. "cWW_5_27" -> labels [5, 27]).
  function attachBPClicks() {
    document.querySelectorAll('path[class*="rnaviewBP"]').forEach((el) => {
      // Strip the plugin's inline tooltip handlers so hovering a base pair
      // does nothing -- only clicks should fire any UI response.
      el.removeAttribute('onmouseover');
      el.removeAttribute('onmouseout');
      if (el.dataset.r2dtBpClickBound) return;
      el.dataset.r2dtBpClickBound = '1';
      el.style.cursor = 'pointer';
      el.addEventListener('click', (ev) => {
        const cls = el.getAttribute('class') || '';
        const m = cls.match(/[a-zA-Z]{3}_(\d+)_(\d+)/);
        if (!m) return;
        ev.stopPropagation();
        selectBasePair(parseInt(m[1]), parseInt(m[2]), el);
      });
    });
  }
  attachBPClicks();

  // Base Pairings List rows: the plugin handles a row click by dispatching
  // a click onto the matching SVG path, which silently does nothing when
  // that path isn't currently rendered (e.g. nested pairs in the non-nested
  // view). Handle row clicks directly so every listed pair updates 3D.
  const bpListId = 'bpListDialog-' + STRUCTURE_ID.toLowerCase();
  document.addEventListener('click', (ev) => {
    const li = ev.target.closest && ev.target.closest('#' + bpListId + ' li');
    if (!li) return;
    // Row text looks like "G5 - C27; cWW" -- pull out the two seq ids.
    const m = (li.textContent || '').match(/(\d+)\D*?-\D*?(\d+)/);
    if (!m) return;
    const a = parseInt(m[1]);
    const b = parseInt(m[2]);
    selectBasePair(a, b, findBPPath(a, b));
  });
  // The plugin re-renders bp paths when the filter changes, so reattach.
  const bpObserver = new MutationObserver(() => {
    if (bpObserver._pending) return;
    bpObserver._pending = true;
    setTimeout(() => { bpObserver._pending = false; attachBPClicks(); }, 0);
  });
  const bpContainer = document.getElementById('pdb-rna-viewer');
  if (bpContainer) bpObserver.observe(bpContainer, { childList: true, subtree: true });
  // 2D hover intentionally does not update the 3D view -- only clicks
  // on nucleotides or base-pair lines cause a 3D selection/focus.

  // 3D -> 2D (already wired inside pdb-rna-viewer when subscribeEvents).
  // We just translate auth -> label in the listener that fires here.
  document.addEventListener('PDB.molstar.click', (ev) => {
    if (!ev.eventData || ev.eventData.auth_asym_id !== CHAIN_ID) return;
    const label = authToLabel[ev.eventData.auth_seq_id];
    if (!label) return;
    document.dispatchEvent(new CustomEvent('protvista-click', {
      detail: { start: label, end: label },
    }));
  });
  document.addEventListener('PDB.molstar.mouseover', (ev) => {
    if (!ev.eventData || ev.eventData.auth_asym_id !== CHAIN_ID) return;
    const label = authToLabel[ev.eventData.auth_seq_id];
    if (!label) return;
    document.dispatchEvent(new CustomEvent('protvista-mouseover', {
      detail: { start: label, end: label },
    }));
  });
  document.addEventListener('PDB.molstar.mouseout', () => {
    document.dispatchEvent(new CustomEvent('protvista-mouseout', { detail: {} }));
  });
})();
