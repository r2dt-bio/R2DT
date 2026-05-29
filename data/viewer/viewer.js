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

  // Attach the MutationObserver BEFORE the polling IIFEs below run their
  // first synchronous tick. setupBackboneToggle inserts the overlay
  // polyline, then tidyBPFilter clicks "All" which triggers the plugin
  // to rebuild the inner group's innerHTML -- wiping the polyline. The
  // observer is what restores it; if it isn't attached yet, the wipe
  // goes unnoticed and the backbone only reappears after the next
  // unrelated DOM mutation (e.g. a hover).
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

  // Relabel the plugin's buttons: "Filter Base Pairings" -> "Filter Base
  // Pairs" (text node only, to keep its help icon and caret) and
  // "Base Pairings List" -> "Base Pair List".
  (function relabelButtons() {
    let attempts = 0;
    const tick = () => {
      const filterBtn = document.getElementById('bpFilterBtn');
      const listBtn = document.querySelector('.bp-list-btn');
      if (!filterBtn || !listBtn) {
        if (attempts++ > 40) return;
        setTimeout(tick, 100);
        return;
      }
      filterBtn.childNodes.forEach((node) => {
        if (node.nodeType === 3 && node.nodeValue.includes('Filter Base Pairings')) {
          node.nodeValue = node.nodeValue.replace('Filter Base Pairings', 'Filter Base Pairs');
        }
      });
      listBtn.textContent = 'Base Pair List';
    };
    tick();
  })();

  // Shorten the Filter Base Pairs "?" help tooltip -- the Leontis-Westhof
  // legend below the viewer now explains the symbols, so the long
  // description is redundant. The plugin builds one <div class="help-tooltip">
  // per help icon at bind time; identify the filter one by its content.
  (function shortenFilterHelp() {
    let attempts = 0;
    const tick = () => {
      const tips = document.querySelectorAll('div.help-tooltip');
      if (tips.length === 0) {
        if (attempts++ > 40) return;
        setTimeout(tick, 100);
        return;
      }
      tips.forEach((tip) => {
        if (/Displays checkboxes/i.test(tip.textContent)) {
          tip.innerHTML = 'Show or hide individual base-pair families.';
        }
      });
    };
    tick();
  })();

  // Best-effort early sweep for the dim/recolour fixups in case the
  // MutationObserver above misses the plugin's initial paint (e.g. when
  // the inner group is populated entirely in the same task as the
  // observer attachment).
  if (unobserved.length || crossingWCPairs.size) {
    let attempts = 0;
    const tick = () => {
      if (applyFixups() || attempts++ > 40) return;
      setTimeout(tick, 100);
    };
    tick();
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
    // Colour the residue (chain-scoped) without molstar's "focus", which
    // routes through the structure-focus manager: that renders the
    // residue plus its surroundings and labels the group with a
    // neighbouring residue -- in protein-RNA complexes a contacting amino
    // acid (e.g. "SER 60"). Instead move only the camera onto the
    // selection loci.
    await molstar.visual.select({
      data: data.map((d) => ({ ...d, color: { r: 255, g: 112, b: 67 }, focus: false })),
      keepRepresentations: true,
    });
    try {
      const loci = molstar.getLociForParams(data);
      const camera = molstar.plugin
        && molstar.plugin.managers
        && molstar.plugin.managers.camera;
      if (loci && camera && camera.focusLoci) {
        camera.focusLoci(loci);
      }
    } catch (err) {
      /* camera focus is best-effort; selection colour already applied */
    }
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

  // ── LBN (layered dot-bracket notation) panel ────────────────────────────
  // Load lbn.json (generated alongside api.json / fr3d.json).  The file may
  // be absent for viewers built before this feature; skip gracefully.
  let lbnData = null;
  try {
    const resp = await fetch('lbn.json');
    if (resp.ok) lbnData = await resp.json();
  } catch (_) { /* ignore */ }

  if (lbnData && lbnData.rows && lbnData.rows.length > 0) {
    const lbnPanel = document.getElementById('lbn-panel');
    const lbnBody  = document.getElementById('lbn-body');
    if (lbnPanel && lbnBody) {
      lbnPanel.style.display = '';
      _renderLBN(lbnData, lbnBody);
    }
  }

  function _renderLBN(data, container) {
    const SEQ = data.sequence;
    const N   = SEQ.length;
    const BLOCK = 100;

    // Color per LW family (label strip of overflow suffix before lookup).
    const COLORS = {
      WC:  '#2ca02c', cWW: '#1f77b4',
      tWW: '#ff7f0e',
      cWH: '#c5573a', tWH: '#c5573a',
      cWS: '#9467bd', tWS: '#9467bd',
      cHW: '#8c564b', tHW: '#8c564b',
      cHH: '#e377c2', tHH: '#e377c2',
      cHS: '#7f7f7f', tHS: '#7f7f7f',
      cSW: '#bcbd22', tSW: '#bcbd22',
      cSH: '#17becf', tSH: '#17becf',
      cSS: '#d62728', tSS: '#d62728',
    };
    function colorOf(label) {
      return COLORS[label.replace(/\(\d+\)$/, '')] || '#555';
    }

    // Build all blocks into a document fragment.
    const frag = document.createDocumentFragment();

    for (let s = 0; s < N; s += BLOCK) {
      const e     = Math.min(s + BLOCK, N);
      const block = document.createElement('div');
      block.className = 'lbn-block';
      block.dataset.blockStart = s + 1;

      let html = '';
      if (N > BLOCK) {
        html += `<div class="lbn-block-header">${s + 1}–${e}</div>`;
      }

      // seq row — every character is a clickable span.
      html += '<div class="lbn-row"><span class="lbn-label">seq         : </span>';
      for (let i = s; i < e; i++) {
        html += `<span data-pos="${i + 1}" class="lbn-nt">${SEQ[i]}</span>`;
      }
      html += '</div>';

      // One row per LW layer.
      for (const row of data.rows) {
        const col    = colorOf(row.label);
        const label  = (row.label + '            ').slice(0, 12);
        html += `<div class="lbn-row"><span class="lbn-label" style="color:${col}">${label}: </span>`;
        for (let i = s; i < e; i++) {
          const pos = i + 1;
          const ch  = row.chars[i];
          if (ch === '.') {
            html += '.';
          } else {
            const partner = row.partners[String(pos)];
            const pAttr   = partner != null ? ` data-partner="${partner}"` : '';
            html += `<span data-pos="${pos}"${pAttr} class="lbn-bp" style="color:${col}">${ch}</span>`;
          }
        }
        html += '</div>';
      }

      block.innerHTML = html;
      frag.appendChild(block);
    }

    container.innerHTML = '';
    container.appendChild(frag);

    // Pre-build position → [span, …] index for O(1) highlight lookups.
    const posSpans = {};
    container.querySelectorAll('[data-pos]').forEach((sp) => {
      const p = +sp.dataset.pos;
      (posSpans[p] = posSpans[p] || []).push(sp);
    });

    // Track what is currently highlighted so we can clear it quickly.
    let highlighted = [];

    function _lbnHighlight(positions) {
      highlighted.forEach((sp) => sp.classList.remove('lbn-selected', 'lbn-partner'));
      highlighted = [];
      if (!positions || positions.length === 0) return;

      // First position is the "clicked" one; rest are partners.
      positions.forEach((pos, idx) => {
        const cls = idx === 0 ? 'lbn-selected' : 'lbn-partner';
        (posSpans[pos] || []).forEach((sp) => {
          sp.classList.add(cls);
          highlighted.push(sp);
        });
      });

      // Scroll the block that contains the first span into view.
      if (highlighted.length > 0) {
        const blk = highlighted[0].closest('.lbn-block');
        if (blk) blk.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
      }
    }

    // LBN click → 2D + 3D.
    container.addEventListener('click', (ev) => {
      const sp = ev.target.closest('[data-pos]');
      if (!sp) return;
      const pos     = +sp.dataset.pos;
      const partner = sp.dataset.partner ? +sp.dataset.partner : null;
      const labels  = partner != null ? [pos, partner] : [pos];

      selectInMolstar(labels);
      document.dispatchEvent(new CustomEvent('protvista-click', {
        detail: { start: pos, end: pos },
      }));
      _lbnHighlight(labels);
    });

    // 2D → LBN.
    document.addEventListener('PDB.RNA.viewer.click', (e) => {
      const labels = labelsFromEvent(e);
      if (labels.length > 0) _lbnHighlight([labels[0]]);
    });

    // 3D → LBN.
    document.addEventListener('PDB.molstar.click', (ev) => {
      if (!ev.eventData || ev.eventData.auth_asym_id !== CHAIN_ID) return;
      const lbl = authToLabel[ev.eventData.auth_seq_id];
      if (lbl) _lbnHighlight([lbl]);
    });

    // Expose for console debugging.
    window.__r2dt.lbnData      = data;
    window.__r2dt.lbnHighlight = _lbnHighlight;
  }
  // ── end LBN ─────────────────────────────────────────────────────────────
})();
