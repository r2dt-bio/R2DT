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

  // Make selected 2D nucleotide labels easier to spot: the plugin only
  // recolours them orange, which is easy to miss on busy diagrams.
  (function enhanceSelectedNucleotideStyle() {
    const svc = window.UiActionsService;
    if (!svc || !svc.colorNucleotide) return;

    const FONT_BUMP_PX = 1.33; // ~1 pt
    const HALO_STROKE = '#ffffff';
    const HALO_WIDTH = '2';
    const BG_FILL = '#fff3b0'; // soft yellow pill behind selected letters
    const BG_PAD = 2;
    const origStyle = new Map(); // label -> saved SVG text attrs
    const selectionBgs = new Map(); // label -> background <rect>

    function nucleotideEl(pdbId, label) {
      return document
        .querySelector('svg.rnaTopoSvg')
        .getElementsByClassName(
          `rnaviewEle rnaviewEle_${pdbId} rnaview_${pdbId}_${label}`
        )[0];
    }

    function setOrRemove(el, attr, val) {
      if (val != null) el.setAttribute(attr, val);
      else el.removeAttribute(attr);
    }

    function removeSelectionBg(label) {
      const existing = selectionBgs.get(label);
      if (existing && existing.parentNode) existing.parentNode.removeChild(existing);
      selectionBgs.delete(label);
    }

    function applySelectionBg(el, label) {
      if (!el || el.nodeName !== 'text' || !el.parentNode) return;
      removeSelectionBg(label);
      let bbox;
      try {
        bbox = el.getBBox();
      } catch (_) {
        return;
      }
      const rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
      rect.setAttribute('class', 'r2dt-nt-selection-bg');
      rect.setAttribute('x', String(bbox.x - BG_PAD));
      rect.setAttribute('y', String(bbox.y - BG_PAD));
      rect.setAttribute('width', String(bbox.width + 2 * BG_PAD));
      rect.setAttribute('height', String(bbox.height + 2 * BG_PAD));
      rect.setAttribute('rx', '3');
      rect.setAttribute('fill', BG_FILL);
      rect.setAttribute('stroke', 'none');
      rect.style.pointerEvents = 'none';
      el.parentNode.insertBefore(rect, el);
      selectionBgs.set(label, rect);
    }

    function applySelectionTypography(el, label) {
      if (!el || el.nodeName !== 'text') return;
      if (!origStyle.has(label)) {
        origStyle.set(label, {
          fw: el.getAttribute('font-weight'),
          fs: el.getAttribute('font-size'),
          stroke: el.getAttribute('stroke'),
          sw: el.getAttribute('stroke-width'),
          slj: el.getAttribute('stroke-linejoin'),
          po: el.getAttribute('paint-order'),
        });
      }
      el.setAttribute('font-weight', 'bold');
      const fs = parseFloat(el.getAttribute('font-size'));
      if (!isNaN(fs)) {
        el.setAttribute('font-size', (fs + FONT_BUMP_PX) + 'px');
      }
      el.setAttribute('stroke', HALO_STROKE);
      el.setAttribute('stroke-width', HALO_WIDTH);
      el.setAttribute('stroke-linejoin', 'round');
      el.setAttribute('paint-order', 'stroke fill');
      applySelectionBg(el, label);
    }

    function restoreTypography(pdbId, label) {
      const stored = origStyle.get(label);
      removeSelectionBg(label);
      if (!stored) return;
      const el = nucleotideEl(pdbId, label);
      if (el && el.nodeName === 'text') {
        setOrRemove(el, 'font-weight', stored.fw);
        setOrRemove(el, 'font-size', stored.fs);
        setOrRemove(el, 'stroke', stored.stroke);
        setOrRemove(el, 'stroke-width', stored.sw);
        setOrRemove(el, 'stroke-linejoin', stored.slj);
        setOrRemove(el, 'paint-order', stored.po);
      }
      origStyle.delete(label);
    }

    const origColor = svc.colorNucleotide.bind(svc);
    svc.colorNucleotide = function (pdbId, label, color, mode) {
      origColor(pdbId, label, color, mode);
      if (mode === 'selection') {
        applySelectionTypography(nucleotideEl(pdbId, label), label);
      }
    };

    const origClear = svc.clearNucleotides.bind(svc);
    svc.clearNucleotides = function (pdbId, mode, labels) {
      if (mode === 'selection') {
        const keys = labels != null ? labels : Array.from(svc.selected.keys());
        keys.forEach((label) => restoreTypography(pdbId, label));
      }
      origClear(pdbId, mode, labels);
    };
  })();

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
  // Cool blue-gray — distinct from letter fill (#000 / #bbb) and from
  // lightened base-pair strokes (#ccc for crossing cWW, coloured families).
  const BACKBONE_STROKE = '#8fa3b3';
  const BACKBONE_OPACITY_SMALL = 0.55;
  const BACKBONE_OPACITY_LARGE = 0.4;
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
    polyline.setAttribute('stroke', BACKBONE_STROKE);
    polyline.setAttribute('stroke-width', isSmall ? '2' : '1.2');
    polyline.setAttribute('stroke-linecap', 'round');
    polyline.setAttribute('stroke-linejoin', 'round');
    polyline.setAttribute('opacity', isSmall ? String(BACKBONE_OPACITY_SMALL) : String(BACKBONE_OPACITY_LARGE));
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
  const VIEWPORT_PADDING = 0.10; // keep ~10% margin so diagrams never hug the edge
  let userAdjustedView = false;

  function computeFitTransform() {
    const svg = document.querySelector('svg.rnaTopoSvg');
    if (!svg) return null;
    const inner = svg.querySelector(`.rnaTopoSvg_${PDB_LOWER}`);
    if (!inner) return null;
    let bbox;
    try {
      bbox = inner.getBBox();
    } catch (_) {
      return null;
    }
    if (!(bbox.width > 0 && bbox.height > 0)) return null;

    const vb = svg.viewBox.baseVal;
    const vw = vb.width || svg.clientWidth;
    const vh = vb.height || svg.clientHeight;
    if (!(vw > 0 && vh > 0)) return null;

    const pad = VIEWPORT_PADDING;
    const k = Math.min(
      (vw * (1 - 2 * pad)) / bbox.width,
      (vh * (1 - 2 * pad)) / bbox.height
    );
    const cx = bbox.x + bbox.width / 2;
    const cy = bbox.y + bbox.height / 2;
    return { k, x: vw / 2 - k * cx, y: vh / 2 - k * cy };
  }

  function applyFitTransform() {
    const tr = computeFitTransform();
    if (!tr) return false;
    const tStr = `translate(${tr.x},${tr.y}) scale(${tr.k})`;
    [
      `.rnaTopoSvg_${PDB_LOWER}`,
      `.rnaTopoSvgHighlight_${PDB_LOWER}`,
      `.rnaTopoSvgSelection_${PDB_LOWER}`,
    ].forEach((sel) => {
      const g = document.querySelector(sel);
      if (g) g.setAttribute('transform', tStr);
    });
    const svg = document.querySelector('svg.rnaTopoSvg');
    if (svg) svg.__zoom = tr;
    const svc = window.UiActionsService;
    if (svc) svc.zoomed = false;
    userAdjustedView = false;
    return true;
  }

  function maybeRefit2D() {
    if (userAdjustedView) return;
    applyFitTransform();
  }

  function bindFitControls() {
    const svg = document.querySelector('svg.rnaTopoSvg');
    if (svg && !svg.dataset.r2dtFitListener) {
      svg.dataset.r2dtFitListener = '1';
      svg.addEventListener('mousedown', () => { userAdjustedView = true; });
    }
    [`#rnaTopologyZoomIn-${PDB_LOWER}`, `#rnaTopologyZoomOut-${PDB_LOWER}`].forEach(
      (sel) => {
        const btn = document.querySelector(sel);
        if (btn && !btn.dataset.r2dtFitBound) {
          btn.dataset.r2dtFitBound = '1';
          btn.addEventListener('click', () => { userAdjustedView = true; });
        }
      }
    );
    const resetBtn = document.querySelector(`#rnaTopologyReset-${PDB_LOWER}`);
    if (resetBtn && !resetBtn.dataset.r2dtFitBound) {
      resetBtn.dataset.r2dtFitBound = '1';
      resetBtn.addEventListener(
        'click',
        (ev) => {
          ev.stopImmediatePropagation();
          ev.preventDefault();
          applyFitTransform();
        },
        true
      );
    }
    const svc = window.UiActionsService;
    if (svc && !svc._r2dtZoomResetPatched) {
      svc._r2dtZoomResetPatched = true;
      svc.zoomReset = function () {
        applyFitTransform();
      };
    }
    return !!resetBtn;
  }

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
        maybeRefit2D();
        bindFitControls();
        setupToolbar();
      }, 0);
    });
    mo.observe(container, { childList: true, subtree: true });
  }

  function isFilterCheckboxVisible(cb) {
    for (let el = cb; el && el.id !== 'checkboxes'; el = el.parentElement) {
      if (el.style.display === 'none') return false;
    }
    return true;
  }

  function updateFilterBadge() {
    const badge = document.getElementById('r2dt-filter-badge');
    if (!badge) return;
    const boxes = document.querySelectorAll('input[id^="Checkbox_"]:not(#Checkbox_All)');
    let total = 0;
    let checked = 0;
    boxes.forEach((cb) => {
      if (!isFilterCheckboxVisible(cb)) return;
      total += 1;
      if (cb.checked) checked += 1;
    });
    badge.textContent = total > 0 ? `${checked}/${total}` : '';
  }

  function ensureFilterPanelTitle() {
    const checkboxes = document.getElementById('checkboxes');
    if (!checkboxes || checkboxes.querySelector('.r2dt-filter-panel-title')) {
      return;
    }
    const title = document.createElement('div');
    title.className = 'r2dt-filter-panel-title';
    title.textContent = 'Leontis-Westhof Base Pairs';
    checkboxes.insertBefore(title, checkboxes.firstChild);
  }

  function removeFilterHelp() {
    const filterBtn = document.getElementById('bpFilterBtn');
    if (filterBtn) {
      filterBtn.querySelector('#bpFilterBtnHelp, .help-icon')?.remove();
      filterBtn.setAttribute('aria-label', 'Base pairs filter and symbol legend');
    }
    document.querySelectorAll('div.help-tooltip').forEach((tip) => {
      if (/Displays checkboxes|base-pair famil/i.test(tip.textContent)) {
        tip.remove();
      }
    });
  }

  const BP_GLYPH_COLOR = '#ccc';
  const BP_FAMILY_GLYPH = {
    cWW: { shapes: ['circle'], filled: true },
    tWW: { shapes: ['circle'], filled: false },
    cWH: { shapes: ['circle', 'square'], filled: true },
    tWH: { shapes: ['circle', 'square'], filled: false },
    cWS: { shapes: ['circle', 'tri-r'], filled: true },
    tWS: { shapes: ['circle', 'tri-r'], filled: false },
    cHH: { shapes: ['square'], filled: true },
    tHH: { shapes: ['square'], filled: false },
    cHS: { shapes: ['square', 'tri-r'], filled: true },
    tHS: { shapes: ['square', 'tri-r'], filled: false },
    cSS: { shapes: ['tri-r'], filled: true },
    tSS: { shapes: ['tri-r'], filled: false },
  };

  function buildBpSymbolSvg(shapes, filled) {
    const width = shapes.length === 1 ? 28 : 40;
    const cy = 8;
    const r = 4;
    const color = BP_GLYPH_COLOR;
    const fill = filled ? color : '#fff';
    const positions = shapes.length === 1 ? [width / 2] : [14, 26];
    let markup =
      `<line x1="2" y1="${cy}" x2="${width - 2}" y2="${cy}" ` +
      `stroke="${color}" stroke-width="1.2"/>`;
    shapes.forEach((kind, idx) => {
      const cx = positions[idx];
      if (kind === 'circle') {
        markup +=
          `<circle cx="${cx}" cy="${cy}" r="${r}" fill="${fill}" ` +
          `stroke="${color}" stroke-width="1.2"/>`;
      } else if (kind === 'square') {
        markup +=
          `<rect x="${cx - r}" y="${cy - r}" width="${2 * r}" height="${2 * r}" ` +
          `fill="${fill}" stroke="${color}" stroke-width="1.2"/>`;
      } else if (kind === 'tri-r') {
        markup +=
          `<polygon points="${cx - r},${cy - r} ${cx - r},${cy + r} ${cx + r},${cy}" ` +
          `fill="${fill}" stroke="${color}" stroke-width="1.2"/>`;
      } else if (kind === 'tri-l') {
        markup +=
          `<polygon points="${cx + r},${cy - r} ${cx + r},${cy + r} ${cx - r},${cy}" ` +
          `fill="${fill}" stroke="${color}" stroke-width="1.2"/>`;
      }
    });
    return (
      `<svg class="r2dt-bp-glyph" width="${width}" height="16" ` +
      `viewBox="0 0 ${width} 16" aria-hidden="true">${markup}</svg>`
    );
  }

  function wrapFilterLabelText(cb, fallbackText) {
    const label = cb.closest('label');
    if (!label) return null;
    let textEl = label.querySelector('.r2dt-bp-family-label');
    if (textEl) return textEl;
    for (const node of [...label.childNodes]) {
      if (node.nodeType === 3 && node.textContent.trim()) {
        textEl = document.createElement('span');
        textEl.className = 'r2dt-bp-family-label';
        textEl.textContent = node.textContent.trim();
        node.replaceWith(textEl);
        return textEl;
      }
    }
    if (!fallbackText) return null;
    textEl = document.createElement('span');
    textEl.className = 'r2dt-bp-family-label';
    textEl.textContent = fallbackText;
    cb.insertAdjacentElement('afterend', textEl);
    return textEl;
  }

  function injectFilterGlyphs() {
    document.querySelectorAll('#checkboxes .r2dt-bp-glyph-wrap').forEach((wrap) => {
      const cb = wrap.closest('label')?.querySelector('input[id^="Checkbox_"]');
      if (!cb || !isFilterCheckboxVisible(cb)) wrap.remove();
    });
    document
      .querySelectorAll('input[id^="Checkbox_"]:not(#Checkbox_All)')
      .forEach((cb) => {
        const label = cb.closest('label');
        if (!label || !isFilterCheckboxVisible(cb)) return;
        const family = cb.id.slice('Checkbox_'.length);
        const spec = BP_FAMILY_GLYPH[family];
        if (!spec) return;
        const textEl = wrapFilterLabelText(cb, family);
        let glyphWrap = label.querySelector('.r2dt-bp-glyph-wrap');
        if (!glyphWrap) {
          glyphWrap = document.createElement('span');
          glyphWrap.className = 'r2dt-bp-glyph-wrap';
        }
        glyphWrap.innerHTML = buildBpSymbolSvg(spec.shapes, spec.filled);
        if (textEl) {
          textEl.insertAdjacentElement('afterend', glyphWrap);
        }
      });
    document.querySelectorAll('input[id^="Checkbox_"]').forEach((cb) => {
      if (cb.id === 'Checkbox_All') wrapFilterLabelText(cb, 'All');
    });
  }

  function isBasePairsPanelOpen() {
    const dropdown = document.querySelector('#mainMenu .menu-dropdown');
    const checkboxes = document.getElementById('checkboxes');
    if (dropdown?.classList.contains('show')) return true;
    if (!checkboxes) return false;
    return getComputedStyle(checkboxes).display !== 'none';
  }

  function closeBasePairsPanel() {
    if (!isBasePairsPanelOpen()) return;
    document.getElementById('bpFilterBtn')?.click();
  }

  function bindToolbarDropdowns() {
    if (document.body.dataset.r2dtDropdownBound) return;
    document.body.dataset.r2dtDropdownBound = '1';

    // Capture phase: the plugin stops propagation on diagram clicks, so a
    // bubble-phase listener on document never runs for clicks on the 2D canvas.
    document.addEventListener(
      'click',
      (ev) => {
        if (!isBasePairsPanelOpen()) return;
        if (ev.target.closest('#mainMenu .menu-dropdown')) return;
        closeBasePairsPanel();
      },
      true
    );
  }

  function mountFilterLegend() {
    document.querySelector('.r2dt-symbols-dropdown')?.remove();

    const checkboxes = document.getElementById('checkboxes');
    if (!checkboxes || checkboxes.querySelector('.r2dt-bp-legend-more')) return true;

    const source = document.getElementById('r2dt-bp-legend-source');
    let panelContent = source?.querySelector('.r2dt-bp-legend-panel');
    if (!panelContent) {
      panelContent = document.querySelector('.r2dt-bp-legend-panel');
    }
    if (!panelContent) return false;

    const details = document.createElement('details');
    details.className = 'r2dt-bp-legend-more';
    const summary = document.createElement('summary');
    summary.textContent = 'More info';
    details.append(summary, panelContent);
    checkboxes.appendChild(details);
    source?.remove();
    return true;
  }

  function bindFilterBadgeListener() {
    const checkboxes = document.getElementById('checkboxes');
    if (!checkboxes || checkboxes.dataset.r2dtBadgeBound) return;
    checkboxes.dataset.r2dtBadgeBound = '1';
    checkboxes.addEventListener('change', updateFilterBadge);
  }

  function setupToolbar() {
    const mainMenu = document.getElementById('mainMenu');
    if (!mainMenu) return false;
    if (mainMenu.dataset.r2dtToolbarBound) {
      bindFilterBadgeListener();
      updateFilterBadge();
      removeFilterHelp();
      ensureFilterPanelTitle();
      injectFilterGlyphs();
      mountFilterLegend();
      bindToolbarDropdowns();
      return true;
    }

    const form = mainMenu.querySelector('form');
    if (!form) return false;

    const menuSelect = form.querySelector('.menuSelectbox');
    const filterDropdown = form.querySelector('.menu-dropdown');
    const nestedWrap = document.getElementById('nestedBP')?.parentElement;
    if (!menuSelect || !filterDropdown) return false;

    mainMenu.classList.add('r2dt-toolbar');
    mainMenu.dataset.r2dtToolbarBound = '1';
    menuSelect.style.display = 'none';

    const titleOverlay = document.querySelector('.pdb-rna-view-title');
    if (titleOverlay) titleOverlay.style.display = 'none';

    const titleEl = document.createElement('span');
    titleEl.className = 'r2dt-toolbar-title';
    titleEl.textContent = CHAIN_ID
      ? `${STRUCTURE_ID} · chain ${CHAIN_ID}`
      : STRUCTURE_ID;

    const viewGroup = document.createElement('div');
    viewGroup.className = 'r2dt-toolbar-group r2dt-toolbar-group--view';
    const viewLabel = document.createElement('span');
    viewLabel.className = 'r2dt-toolbar-group-label';
    viewLabel.textContent = 'View';

    const backboneLabel = document.createElement('label');
    backboneLabel.className = 'r2dt-toggle';
    backboneLabel.htmlFor = 'r2dt-backbone-toggle';
    backboneLabel.innerHTML =
      '<input type="checkbox" id="r2dt-backbone-toggle" checked ' +
      'aria-label="Show backbone path">' +
      '<span class="r2dt-toggle-track" aria-hidden="true"></span>' +
      '<span class="r2dt-toggle-label">Backbone</span>';
    backboneLabel.querySelector('input').addEventListener('change', (ev) => {
      backboneVisible = ev.target.checked;
      injectBackboneOverlay();
    });

    viewGroup.append(viewLabel, backboneLabel);

    const pairsGroup = document.createElement('div');
    pairsGroup.className = 'r2dt-toolbar-group r2dt-toolbar-group--pairs';
    const pairsLabel = document.createElement('span');
    pairsLabel.className = 'r2dt-toolbar-group-label';
    pairsLabel.textContent = 'Base pairs';

    pairsGroup.append(pairsLabel, filterDropdown);
    if (nestedWrap) {
      const nestedText = nestedWrap.querySelector('label[for="nestedBP"] span');
      if (nestedText) nestedText.textContent = 'Nested only';
      const nestedInput = nestedWrap.querySelector('#nestedBP');
      if (nestedInput) {
        nestedInput.setAttribute('aria-label', 'Only nested base pairs');
      }
      pairsGroup.append(nestedWrap);
    }

    form.insertBefore(titleEl, form.firstChild);
    form.insertBefore(viewGroup, menuSelect);
    form.insertBefore(pairsGroup, menuSelect);

    const filterBtn = document.getElementById('bpFilterBtn');
    if (filterBtn) {
      filterBtn.childNodes.forEach((node) => {
        if (node.nodeType === 3 && /Filter/i.test(node.nodeValue)) {
          node.nodeValue = 'Base Pairs';
        }
      });
      if (!document.getElementById('r2dt-filter-badge')) {
        const badge = document.createElement('span');
        badge.id = 'r2dt-filter-badge';
        badge.className = 'r2dt-filter-badge';
        badge.setAttribute('aria-hidden', 'true');
        const icon = document.getElementById('bpFilterBtnIcon');
        if (icon) filterBtn.insertBefore(badge, icon);
        else filterBtn.appendChild(badge);
      }
      bindFilterBadgeListener();
      updateFilterBadge();
      removeFilterHelp();
      ensureFilterPanelTitle();
      injectFilterGlyphs();
      mountFilterLegend();
    }

    bindToolbarDropdowns();
    injectBackboneOverlay();
    return true;
  }

  (function initToolbar() {
    let attempts = 0;
    const tick = () => {
      if (setupToolbar()) return;
      if (attempts++ > 40) return;
      setTimeout(tick, 100);
    };
    tick();
  })();

  // Centre the 2D diagram with a comfortable margin. The plugin's default
  // viewBox starts at (0,0) and often places residues flush against an
  // edge; re-fit once content (and the backbone overlay) is in the DOM.
  (function setup2DFit() {
    let attempts = 0;
    const tick = () => {
      bindFitControls();
      if (applyFitTransform()) return;
      if (attempts++ > 40) return;
      setTimeout(tick, 100);
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
      ensureFilterPanelTitle();
      injectFilterGlyphs();
      mountFilterLegend();
      updateFilterBadge();
    };
    tick();
  })();

  // Relabel the plugin's "Base Pairings List" -> "Base Pair List".
  // (Filter button label is handled in setupToolbar.)
  (function relabelButtons() {
    let attempts = 0;
    const tick = () => {
      const listBtn = document.querySelector('.bp-list-btn');
      if (!listBtn) {
        if (attempts++ > 40) return;
        setTimeout(tick, 100);
        return;
      }
      listBtn.textContent = 'Base Pair List';
    };
    tick();
  })();

  // The plugin injects a "?" help icon and floating tooltip on the filter
  // button; remove them once the base-pairs panel is wired up.
  (function stripFilterHelp() {
    let attempts = 0;
    const tick = () => {
      removeFilterHelp();
      const helpLeft = document.getElementById('bpFilterBtnHelp');
      const tipLeft = [...document.querySelectorAll('div.help-tooltip')].some(
        (tip) => /Displays checkboxes|base-pair famil/i.test(tip.textContent)
      );
      if (!helpLeft && !tipLeft) return;
      if (attempts++ > 40) return;
      setTimeout(tick, 100);
    };
    tick();
  })();

  (function initFilterLegend() {
    let attempts = 0;
    const tick = () => {
      if (mountFilterLegend()) return;
      if (attempts++ > 40) return;
      setTimeout(tick, 100);
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

      if (partner != null) {
        // Clicking a bracket selects the whole pair: highlight its glyph/line
        // in 2D and both partners in 3D, exactly like clicking the pair in the
        // 2D diagram. (No protvista-click, which would re-render and wipe the
        // path highlight.)
        selectBasePair(pos, partner, findBPPath(pos, partner));
        _lbnHighlight([pos, partner]);
      } else {
        // A lone position: highlight just that nucleotide in 2D + 3D.
        selectInMolstar([pos]);
        document.dispatchEvent(new CustomEvent('protvista-click', {
          detail: { start: pos, end: pos },
        }));
        _lbnHighlight([pos]);
      }
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
