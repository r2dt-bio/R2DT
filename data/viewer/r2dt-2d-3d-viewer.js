/*
 * R2DT 2D+3D viewer glue (r2dt-2d-3d-viewer.js) — standalone page and embed API.
 */

(function (global) {
  'use strict';

  let activeViewer = null;

  function normalizeBaseUrl(url) {
    if (!url) return './';
    return url.endsWith('/') ? url : url + '/';
  }

  function resolveUrl(baseUrl, path) {
    const base = normalizeBaseUrl(baseUrl);
    if (/^https?:\/\//i.test(base)) {
      return new URL(path.replace(/^\.\//, ''), base).href;
    }
    return base + path.replace(/^\.\//, '');
  }

  async function fetchJson(url) {
    const resp = await fetch(url);
    if (!resp.ok) throw new Error(`Failed to fetch ${url}: ${resp.status}`);
    return resp.json();
  }

  async function loadManifest(baseUrl) {
    try {
      const resp = await fetch(resolveUrl(baseUrl, 'manifest.json'));
      if (resp.ok) return resp.json();
    } catch (_) { /* optional */ }
    return null;
  }

  function resolveMount(mount) {
    if (!mount) throw new Error('R2DTViewer.create: mount is required');
    if (typeof mount === 'string') {
      const el = document.querySelector(mount);
      if (!el) throw new Error(`R2DTViewer.create: mount not found: ${mount}`);
      return el;
    }
    if (mount.nodeType === 1) return mount;
    throw new Error('R2DTViewer.create: mount must be a selector or Element');
  }

  function disableAutoZoomOnce() {
    if (global.__r2dtZoomDisabled) return;
    global.__r2dtZoomDisabled = true;
    if (window.UiActionsService && window.UiActionsService.zoomToNucleotides) {
      window.UiActionsService.zoomToNucleotides = function () {};
    }
  }

  function enhanceSelectedNucleotideStyleOnce() {
    const svc = window.UiActionsService;
    if (!svc || !svc.colorNucleotide || svc.__r2dtStyleEnhanced) return;
    svc.__r2dtStyleEnhanced = true;

    const FONT_BUMP_PX = 1.33;
    const HALO_STROKE = '#ffffff';
    const HALO_WIDTH = '2';
    const BG_FILL = '#fff3b0';
    const BG_PAD = 2;
    const origStyle = new Map();
    const selectionBgs = new Map();

    function nucleotideEl(pdbId, label) {
      return document.getElementsByClassName(
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
  }

  function labelsFromEvent(e) {
    const d = e.eventData || e.detail || {};
    if (Array.isArray(d.label_seq_ids)) return d.label_seq_ids;
    if (d.label_seq_id !== undefined && d.label_seq_id !== null) return [d.label_seq_id];
    return [];
  }

  function buildLabelMaps(apiData) {
    const labelToAuth = {};
    apiData.label_seq_ids.forEach((label, i) => {
      if (label !== null && label !== undefined) {
        labelToAuth[label] = apiData.auth_seq_ids[i];
      }
    });
    const authToLabel = {};
    Object.entries(labelToAuth).forEach(([label, auth]) => {
      authToLabel[auth] = parseInt(label);
    });
    return { labelToAuth, authToLabel };
  }

  async function resolvePanelData(baseUrl, opts) {
    const normalized = normalizeBaseUrl(baseUrl);
    const manifest = await loadManifest(normalized);
    const structureId = opts.structureId || manifest?.structureId;
    const chainId = opts.chainId ?? manifest?.chainId ?? '';
    const structureFormat = opts.structureFormat || manifest?.structureFormat || 'cif';
    let structureUrl = opts.structureUrl || manifest?.structureUrl;
    if (!structureId) {
      throw new Error('structureId required (or manifest.json)');
    }
    if (!structureUrl) {
      const ext = structureFormat === 'pdb' ? 'pdb' : 'cif';
      structureUrl = `${structureId}.${ext}`;
    }
    const resolvedStructureUrl = resolveUrl(normalized, structureUrl);
    const apiData = opts.apiData || await fetchJson(resolveUrl(normalized, 'api.json'));
    const fr3dData = opts.fr3dData || await fetchJson(resolveUrl(normalized, 'fr3d.json'));
    return {
      baseUrl: normalized,
      resolveUrl,
      structureId,
      chainId,
      structureUrl: resolvedStructureUrl,
      structureFormat,
      apiData,
      fr3dData,
      PDB_LOWER: structureId.toLowerCase(),
    };
  }

  function installMultiPanelDomShim() {
    if (global.__r2dtMultiPanelShim) return;
    global.__r2dtMultiPanelShim = true;

    let activeRoot = null;
    document.addEventListener(
      'pointerdown',
      (ev) => {
        const root = ev.target.closest('.r2dt-viewer-root');
        if (root) activeRoot = root;
      },
      true
    );

    const origGetById = document.getElementById.bind(document);
    document.getElementById = function (id) {
      if (/-rnaTopology/.test(id)) return origGetById(id);
      if (activeRoot) {
        const inRoot = activeRoot.querySelector('#' + CSS.escape(id));
        if (inRoot) return inRoot;
      }
      return origGetById(id);
    };

    const origQuery = document.querySelector.bind(document);
    document.querySelector = function (selector) {
      if (activeRoot && selector !== 'svg.rnaTopoSvg') {
        const inRoot = activeRoot.querySelector(selector);
        if (inRoot) return inRoot;
      }
      return origQuery(selector);
    };
  }

  function installFetchShim(routes) {
    if (global.__r2dtFetchShim || !routes || routes.length === 0) return;
    global.__r2dtFetchShim = true;
    const orig = window.fetch.bind(window);
    window.fetch = function (input, init) {
      const url = typeof input === 'string' ? input : (input && input.url) || '';
      if (url.includes('ebi.ac.uk/pdbe/static/entry/')) {
        const isBp = url.endsWith('_basepair.json');
        const file = isBp ? 'fr3d.json' : 'api.json';
        const lower = url.toLowerCase();
        for (let i = 0; i < routes.length; i++) {
          const route = routes[i];
          if (lower.includes(String(route.pdbIdMatch).toLowerCase())) {
            return orig(resolveUrl(route.baseUrl, file), init);
          }
        }
      }
      return orig(input, init);
    };
  }

  async function renderMolstarPlugin(panel3d, structureUrl, structureFormat) {
    const molstar = new PDBeMolstarPlugin();
    await new Promise((resolve) => {
      molstar.render(
        panel3d,
        {
          customData: { url: structureUrl, format: structureFormat, binary: false },
          subscribeEvents: true,
          bgColor: { r: 255, g: 255, b: 255 },
          hideControls: true,
          hideCanvasControls: ['expand'],
          sequencePanel: false,
          loadingOverlay: true,
        }
      );
      if (molstar.events && molstar.events.loadComplete) {
        const sub = molstar.events.loadComplete.subscribe((loaded) => {
          if (loaded) { sub.unsubscribe(); resolve(); }
        });
      } else {
        setTimeout(resolve, 1500);
      }
    });
    return molstar;
  }

  function createMolstarSelector(molstar, labelToAuth, chainId) {
    return async function selectInMolstar(labels) {
      if (!molstar) return;
      const data = labels.map((l) => {
        const auth = labelToAuth[l];
        if (auth === undefined || auth === null) return null;
        return {
          auth_asym_id: chainId,
          start_auth_residue_number: auth,
          end_auth_residue_number: auth,
        };
      }).filter((d) => d !== null);
      if (data.length === 0) return;
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
      } catch (_) { /* best-effort */ }
    };
  }

  const BP_GLYPH_COLOR = '#909090';
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

  const BP_LEGEND_ROWS = [
    { cisLabel: 'cWW', cisGlyphs: [['circle']], transLabel: 'tWW', transGlyphs: [['circle']] },
    {
      cisLabel: 'cWH, cHW',
      cisGlyphs: [['circle', 'square'], ['square', 'circle']],
      transLabel: 'tWH, tHW',
      transGlyphs: [['circle', 'square'], ['square', 'circle']],
    },
    {
      cisLabel: 'cWS, cSW',
      cisGlyphs: [['circle', 'tri-r'], ['tri-l', 'circle']],
      transLabel: 'tWS, tSW',
      transGlyphs: [['circle', 'tri-r'], ['tri-l', 'circle']],
    },
    { cisLabel: 'cHH', cisGlyphs: [['square']], transLabel: 'tHH', transGlyphs: [['square']] },
    {
      cisLabel: 'cHS, cSH',
      cisGlyphs: [['square', 'tri-r'], ['tri-l', 'square']],
      transLabel: 'tHS, tSH',
      transGlyphs: [['square', 'tri-r'], ['tri-l', 'square']],
    },
    { cisLabel: 'cSS', cisGlyphs: [['tri-r']], transLabel: 'tSS', transGlyphs: [['tri-r']] },
  ];

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

  function createBpLegendPanel() {
    const panel = document.createElement('div');
    panel.className = 'r2dt-bp-legend-panel';

    const intro = document.createElement('p');
    intro.className = 'r2dt-bp-legend-intro';
    intro.innerHTML =
      'Each glyph encodes the two interacting edges and glycosidic-bond orientation. ' +
      'Edge shape: circle = Watson–Crick, square = Hoogsteen, triangle = Sugar. ' +
      'Fill: solid = <em>cis</em>, open = <em>trans</em>.';

    const table = document.createElement('table');
    table.className = 'r2dt-bp-legend-table';
    const thead = document.createElement('thead');
    thead.innerHTML = '<tr><th colspan="2">cis</th><th colspan="2">trans</th></tr>';
    const tbody = document.createElement('tbody');
    BP_LEGEND_ROWS.forEach((row) => {
      const tr = document.createElement('tr');
      const cisTh = document.createElement('th');
      cisTh.textContent = row.cisLabel;
      const cisTd = document.createElement('td');
      cisTd.innerHTML = row.cisGlyphs
        .map((shapes) => buildBpSymbolSvg(shapes, true))
        .join('');
      const transTh = document.createElement('th');
      transTh.textContent = row.transLabel;
      const transTd = document.createElement('td');
      transTd.innerHTML = row.transGlyphs
        .map((shapes) => buildBpSymbolSvg(shapes, false))
        .join('');
      tr.append(cisTh, cisTd, transTh, transTd);
      tbody.appendChild(tr);
    });
    table.append(thead, tbody);

    const cite = document.createElement('p');
    cite.className = 'r2dt-bp-legend-cite';
    cite.innerHTML =
      'Symbols follow <a href="https://pubmed.ncbi.nlm.nih.gov/11345429/" ' +
      'target="_blank" rel="noopener">Leontis &amp; Westhof, 2001</a>.';

    panel.append(intro, table, cite);
    return panel;
  }

  function buildViewerDom(mountEl, opts) {
    mountEl.innerHTML = '';
    const root = document.createElement('div');
    root.className = 'r2dt-viewer-root';
    root.dataset.structure = opts.structureId || '';
    if (opts.layout) root.dataset.layout = opts.layout;
    if (opts.height != null) {
      root.style.setProperty('--r2dt-viewer-height', typeof opts.height === 'number' ? `${opts.height}px` : String(opts.height));
    }
    if (opts.panelWidth != null) {
      root.style.setProperty('--r2dt-panel-size', typeof opts.panelWidth === 'number' ? `${opts.panelWidth}px` : String(opts.panelWidth));
    }

    const vis = document.createElement('div');
    vis.className = 'r2dt-viewer-vis';

    const panel2d = document.createElement('div');
    panel2d.id = 'pdb-rna-viewer';
    panel2d.className = 'r2dt-panel r2dt-panel--2d';

    const panel3d = document.createElement('div');
    panel3d.id = 'pdb-molstar';
    panel3d.className = 'r2dt-panel r2dt-panel--3d';
    if (opts.layout === '2d-only') panel3d.hidden = true;
    if (opts.layout === '3d-only') panel2d.hidden = true;

    vis.append(panel2d, panel3d);
    root.appendChild(vis);

    if (opts.showLbn !== false) {
      const lbnPanel = document.createElement('div');
      lbnPanel.id = 'lbn-panel';
      lbnPanel.className = 'r2dt-viewer-lbn';
      lbnPanel.hidden = true;
      const lbnTitle = document.createElement('h2');
      lbnTitle.className = 'r2dt-viewer-lbn-title';
      lbnTitle.textContent = 'Layered dot-bracket notation (Leontis–Westhof base pairs)';
      const lbnBody = document.createElement('div');
      lbnBody.id = 'lbn-body';
      lbnBody.className = 'lbn-body';
      lbnBody.textContent = 'Loading…';
      lbnPanel.append(lbnTitle, lbnBody);
      root.appendChild(lbnPanel);
    }

    if (opts.showLegend !== false) {
      root.appendChild(createBpLegendPanel());
    }

    mountEl.appendChild(root);
    return { root, panel2d, panel3d };
  }

  async function initViewer(ctx) {
    const {
      root, panel2d, panel3d, baseUrl, resolveUrl,
      STRUCTURE_ID, CHAIN_ID, STRUCTURE_URL, STRUCTURE_FORMAT, PDB_LOWER,
      apiData, fr3dData, showLbn,
    } = ctx;
    const link3d = ctx.link3d !== false;

    disableAutoZoomOnce();
    enhanceSelectedNucleotideStyleOnce();

    const { labelToAuth, authToLabel } = buildLabelMaps(apiData);

  // --- 2D viewer ---
  const rnaPlugin = new PdbRnaViewerPlugin();
  await rnaPlugin.render(
    panel2d,
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

  // The plugin's pathOrNucleotide() rebuilds the entire SVG inner group on
  // every filter checkbox change, which makes all base-pair lines blink.
  // Toggle path visibility instead once the diagram is painted.
  let bpPathsMaterialized = false;
  let onBpFilterUpdated = () => {};

  const BP_PATH_ID_RE = /\b([ct](?:WW|WH|WS|HH|HS|SS)_\d+_\d+)\b/;

  function pathIdsInDisplayHtml(html) {
    const ids = new Set();
    if (!html) return ids;
    let m;
    const re = new RegExp(BP_PATH_ID_RE.source, 'g');
    while ((m = re.exec(html))) ids.add(m[1]);
    return ids;
  }

  function bpPathIdFromElement(el) {
    const m = (el.getAttribute('class') || '').match(BP_PATH_ID_RE);
    return m ? m[1] : null;
  }

  function materializeAllBpPaths() {
    const ui = rnaPlugin.uiTemplateService;
    if (bpPathsMaterialized || !ui?.baseStrs) return;
    const inner = root.querySelector(`.rnaTopoSvg_${PDB_LOWER}`);
    if (!inner) return;
    const existing = new Set();
    inner.querySelectorAll('path.rnaviewBP').forEach((p) => {
      const id = bpPathIdFromElement(p);
      if (id) existing.add(id);
    });
    const toAdd = [];
    ui.baseStrs.forEach((entry) => {
      entry[1].forEach((html) => {
        const m = html.match(BP_PATH_ID_RE);
        if (m && !existing.has(m[1])) toAdd.push(html);
      });
    });
    if (toAdd.length) inner.insertAdjacentHTML('beforeend', toAdd.join(''));
    bpPathsMaterialized = true;
  }

  function applyBpVisibility() {
    const ui = rnaPlugin.uiTemplateService;
    if (!ui) return;
    materializeAllBpPaths();
    const nested = root.querySelector('#nestedBP')?.checked;
    const visible = pathIdsInDisplayHtml(
      nested ? ui.displayNestedBaseStrs : ui.displayBaseStrs
    );
    const inner = root.querySelector(`.rnaTopoSvg_${PDB_LOWER}`);
    if (!inner) return;
    inner.querySelectorAll('path.rnaviewBP').forEach((p) => {
      const id = bpPathIdFromElement(p);
      if (!id) return;
      p.style.display = visible.has(id) ? '' : 'none';
    });
    onBpFilterUpdated();
  }

  const uiSvc = rnaPlugin.uiTemplateService;
  if (uiSvc?.pathOrNucleotide) {
    const origPathOrNucleotide = uiSvc.pathOrNucleotide.bind(uiSvc);
    uiSvc.pathOrNucleotide = function r2dtPathOrNucleotide() {
      const menuSelect = uiSvc.containerElement?.querySelector('.menuSelectbox');
      const mode = menuSelect ? parseInt(menuSelect.value, 10) : 0;
      const inner = root.querySelector(`.rnaTopoSvg_${PDB_LOWER}`);
      if (mode === 0 && inner && inner.querySelector('text')) {
        applyBpVisibility();
        uiSvc.renderBpListDialog(false);
        return;
      }
      bpPathsMaterialized = false;
      origPathOrNucleotide();
    };
  }

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
  // Low tension keeps the curve close to the residue centres — just softens corners.
  const BACKBONE_SMOOTH_TENSION = 0.18;
  // Derive nucleotide centre coordinates from apiData.svg_paths.
  // Each entry past the two dummy heads encodes prevX,prevY,x,y after
  // splitting on "M" and ","; the "current" point sits at tokens[2..3].
  function buildBackbonePointList() {
    const pts = [];
    const paths = apiData.svg_paths || [];
    for (let i = 2; i < paths.length; i++) {
      const tokens = paths[i].split(/[M,]/).filter((s) => s !== '');
      if (tokens.length < 4) continue;
      const x = parseFloat(tokens[2]);
      const y = parseFloat(tokens[3]);
      if (!isNaN(x) && !isNaN(y)) pts.push({ x, y });
    }
    return pts;
  }

  function backboneSmoothPath(pts, tension) {
    if (pts.length < 2) return '';
    if (pts.length === 2) {
      return `M ${pts[0].x},${pts[0].y} L ${pts[1].x},${pts[1].y}`;
    }
    let d = `M ${pts[0].x},${pts[0].y}`;
    for (let i = 0; i < pts.length - 1; i++) {
      const p0 = pts[i - 1] || pts[i];
      const p1 = pts[i];
      const p2 = pts[i + 1];
      const p3 = pts[i + 2] || p2;
      const cp1x = p1.x + (p2.x - p0.x) * tension;
      const cp1y = p1.y + (p2.y - p0.y) * tension;
      const cp2x = p2.x - (p3.x - p1.x) * tension;
      const cp2y = p2.y - (p3.y - p1.y) * tension;
      d += ` C ${cp1x},${cp1y} ${cp2x},${cp2y} ${p2.x},${p2.y}`;
    }
    return d;
  }

  function injectBackboneOverlay() {
    const svg = root.querySelector('svg.rnaTopoSvg');
    if (!svg) return false;
    const existing = svg.querySelector('.r2dt-backbone-overlay');
    if (!backboneVisible) {
      if (existing) existing.remove();
      return true;
    }
    // Idempotent: if the overlay is already present, do nothing. This
    // matters because the re-inject runs from a MutationObserver -- adding
    // the path is itself a mutation, so without this guard we'd loop.
    if (existing) return true;
    const ptList = buildBackbonePointList();
    if (!ptList.length) return false;
    const pathD = backboneSmoothPath(ptList, BACKBONE_SMOOTH_TENSION);
    if (!pathD) return false;
    const inner = svg.querySelector('[class^="rnaTopoSvg_"]');
    if (!inner) return false;
    const pathEl = document.createElementNS('http://www.w3.org/2000/svg', 'path');
    pathEl.setAttribute('class', 'r2dt-backbone-overlay');
    pathEl.setAttribute('d', pathD);
    pathEl.setAttribute('fill', 'none');
    // Scale the overlay with structure size: small RNAs (~100 nt) need a
    // thicker, less transparent line to be legible against the letters;
    // large rRNAs (1000+ nt) look noisy at the same settings.
    const nPts = ptList.length;
    const isSmall = nPts < 500;
    pathEl.setAttribute('stroke', BACKBONE_STROKE);
    pathEl.setAttribute('stroke-width', isSmall ? '2' : '1.2');
    pathEl.setAttribute('stroke-linecap', 'round');
    pathEl.setAttribute('stroke-linejoin', 'round');
    pathEl.setAttribute('opacity', isSmall ? String(BACKBONE_OPACITY_SMALL) : String(BACKBONE_OPACITY_LARGE));
    pathEl.style.pointerEvents = 'none';
    // Insert as the first child of the plugin's inner group so the same
    // zoom/pan transform applies to the overlay as to the nucleotides.
    inner.insertBefore(pathEl, inner.firstChild);
    return true;
  }
  // Declared up front so the post-render fixup pass (which skips the
  // currently-selected BP) can reference it before the click handler
  // block initializes -- temporal-dead-zone errors otherwise.
  let lastBPSelected = null;

  const VIEWPORT_PADDING = 0.05; // keep ~5% margin so diagrams never hug the edge
  let userAdjustedView = false;

  // Post-render fixups: (a) dim nucleotide letters for unobserved
  // residues -- the plugin's unobservedColor theme only colors backbone,
  // not text; (b) lighten/thin long-range non-nested pairs (crossing>0) so
  // they don't dominate the nested cWW ladder. The minified plugin's
  // async render() doesn't reliably wait for the DOM, so poll briefly.
  function computeFitTransform() {
    const svg = root.querySelector('svg.rnaTopoSvg');
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

  function parseViewTransform(str) {
    if (!str) return null;
    const m = str.match(
      /translate\(([-\d.e+]+),([-\d.e+]+)\)\s*scale\(([-\d.e+]+)\)/
    );
    if (!m) return null;
    return { x: parseFloat(m[1]), y: parseFloat(m[2]), k: parseFloat(m[3]) };
  }

  // Minimal d3 ZoomTransform stand-in so svg.__zoom stays compatible with
  // the plugin's pan/zoom behaviour (d3 initially stores an accessor fn).
  function createZoomTransform(k, x, y) {
    const t = { k, x, y };
    t.scale = function scaleBy(f) {
      return f === 1 ? t : createZoomTransform(t.k * f, t.x, t.y);
    };
    t.translate = function translateBy(dx, dy) {
      return dx === 0 && dy === 0
        ? t
        : createZoomTransform(t.k, t.x + t.k * dx, t.y + t.k * dy);
    };
    t.invert = function invert(point) {
      return [(point[0] - t.x) / t.k, (point[1] - t.y) / t.k];
    };
    t.apply = function apply(point) {
      return [t.k * point[0] + t.x, t.k * point[1] + t.y];
    };
    t.toString = function toString() {
      return `translate(${t.x},${t.y}) scale(${t.k})`;
    };
    return t;
  }

  function getCurrentViewTransform() {
    const inner = root.querySelector(`.rnaTopoSvg_${PDB_LOWER}`);
    const parsed = parseViewTransform(inner && inner.getAttribute('transform'));
    if (parsed) return parsed;
    return computeFitTransform() || { k: 1, x: 0, y: 0 };
  }

  function applyViewTransform(tr) {
    const tStr = `translate(${tr.x},${tr.y}) scale(${tr.k})`;
    [
      `.rnaTopoSvg_${PDB_LOWER}`,
      `.rnaTopoSvgHighlight_${PDB_LOWER}`,
      `.rnaTopoSvgSelection_${PDB_LOWER}`,
    ].forEach((sel) => {
      const g = root.querySelector(sel);
      if (g) g.setAttribute('transform', tStr);
    });
    const svg = root.querySelector('svg.rnaTopoSvg');
    if (svg) svg.__zoom = createZoomTransform(tr.k, tr.x, tr.y);
  }

  function scaleViewTransform(factor) {
    const svg = root.querySelector('svg.rnaTopoSvg');
    if (!svg) return;
    const vb = svg.viewBox.baseVal;
    const cx = vb.width / 2;
    const cy = vb.height / 2;
    const tr = getCurrentViewTransform();
    applyViewTransform({
      k: tr.k * factor,
      x: cx - factor * (cx - tr.x),
      y: cy - factor * (cy - tr.y),
    });
  }

  function applyFitTransform() {
    const tr = computeFitTransform();
    if (!tr) return false;
    applyViewTransform(tr);
    const svc = window.UiActionsService;
    if (svc) svc.zoomed = false;
    return true;
  }

  function maybeRefit2D() {
    if (userAdjustedView) return;
    applyFitTransform();
  }

  function bindFitControls() {
    const svg = root.querySelector('svg.rnaTopoSvg');
    if (svg && !svg.dataset.r2dtFitListener) {
      svg.dataset.r2dtFitListener = '1';
      svg.addEventListener('mousedown', () => { userAdjustedView = true; });
    }
    const zoomIn = root.querySelector(`#rnaTopologyZoomIn-${PDB_LOWER}`);
    if (zoomIn && !zoomIn.dataset.r2dtZoomBound) {
      zoomIn.dataset.r2dtZoomBound = '1';
      zoomIn.addEventListener(
        'click',
        (ev) => {
          ev.stopImmediatePropagation();
          ev.preventDefault();
          userAdjustedView = true;
          scaleViewTransform(1.2);
        },
        true
      );
    }
    const zoomOut = root.querySelector(`#rnaTopologyZoomOut-${PDB_LOWER}`);
    if (zoomOut && !zoomOut.dataset.r2dtZoomBound) {
      zoomOut.dataset.r2dtZoomBound = '1';
      zoomOut.addEventListener(
        'click',
        (ev) => {
          ev.stopImmediatePropagation();
          ev.preventDefault();
          const tr = getCurrentViewTransform();
          const fit = computeFitTransform();
          const nextK = tr.k * 0.8;
          if (fit && (tr.k <= fit.k * 1.01 || nextK <= fit.k)) {
            userAdjustedView = false;
            applyFitTransform();
          } else {
            userAdjustedView = true;
            scaleViewTransform(0.8);
          }
        },
        true
      );
    }
    const resetBtn = root.querySelector(`#rnaTopologyReset-${PDB_LOWER}`);
    if (resetBtn && !resetBtn.dataset.r2dtFitBound) {
      resetBtn.dataset.r2dtFitBound = '1';
      resetBtn.addEventListener(
        'click',
        (ev) => {
          ev.stopImmediatePropagation();
          ev.preventDefault();
          userAdjustedView = false;
          applyFitTransform();
        },
        true
      );
    }
    const svc = window.UiActionsService;
    if (svc && !svc._r2dtZoomResetPatched) {
      svc._r2dtZoomResetPatched = true;
      svc.zoomReset = function () {
        userAdjustedView = false;
        applyFitTransform();
      };
    }
    return !!resetBtn;
  }

  const unobserved = apiData.unobserved_label_seq_ids || [];
  const CWW_STROKE = '#888888';
  const CWW_CROSSING_STROKE = '#aaaaaa';
  const CROSSING_BP_STROKE_WIDTH = '1.4';
  const crossingPairs = new Set();
  (fr3dData.annotations || []).forEach((a) => {
    if (a.crossing && String(a.crossing) !== '0') {
      crossingPairs.add(`${a.seq_id1}_${a.seq_id2}`);
      crossingPairs.add(`${a.seq_id2}_${a.seq_id1}`);
    }
  });

  function parseBpPathClass(cls) {
    const m = cls.match(/([a-zA-Z]{3})_(\d+)_(\d+)/);
    if (!m) return null;
    return { bp: m[1], a: m[2], b: m[3] };
  }

  function applyFixups() {
    let any = false;
    unobserved.forEach((seqId) => {
      root.querySelectorAll(`text.rnaview_${PDB_LOWER}_${seqId}`)
        .forEach((el) => { el.setAttribute('fill', '#bbbbbb'); any = true; });
    });
    root.querySelectorAll('path.rnaviewBP').forEach((el) => {
      const parsed = parseBpPathClass(el.getAttribute('class') || '');
      if (!parsed) return;
      const isCrossing = crossingPairs.has(`${parsed.a}_${parsed.b}`);
      if (parsed.bp === 'cWW' && el !== lastBPSelected) {
        const stroke = isCrossing ? CWW_CROSSING_STROKE : CWW_STROKE;
        el.setAttribute('stroke', stroke);
        // Remember the greyed-out value so the click handler's
        // restore-on-next-click brings it back to grey, not the
        // plugin's default black.
        el.dataset.r2dtOrigStroke = stroke;
      }
      if (isCrossing) {
        el.setAttribute('stroke-width', CROSSING_BP_STROKE_WIDTH);
      }
      any = true;
    });
    return any;
  }

  // Attach the MutationObserver BEFORE the polling IIFEs below run their
  // first synchronous tick. setupBackboneToggle inserts the overlay
  // polyline, then tidyBPFilter clicks "All" which triggers the plugin
  // to rebuild the inner group's innerHTML -- wiping the polyline. The
  // observer is what restores it; if it isn't attached yet, the wipe
  // goes unnoticed and the backbone only reappears after the next
  // unrelated DOM mutation (e.g. a hover).
  const container = panel2d;
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
    const td = cb.closest('td');
    if (td && td.style.display === 'none') return false;
    for (let el = cb; el && el.id !== 'checkboxes'; el = el.parentElement) {
      if (el.style.display === 'none') return false;
    }
    return true;
  }

  function hideEmptyFilterRows() {
    root.querySelectorAll('#checkboxes tr').forEach((tr) => {
      const anyVisible = [...tr.querySelectorAll('td')].some(
        (td) => td.style.display !== 'none'
      );
      tr.style.display = anyVisible ? '' : 'none';
    });
  }

  function updateFilterBadge() {
    const badge = root.querySelector('#r2dt-filter-badge');
    if (!badge) return;
    const boxes = root.querySelectorAll('input[id^="Checkbox_"]:not(#Checkbox_All)');
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
    const checkboxes = root.querySelector('#checkboxes');
    if (!checkboxes || checkboxes.querySelector('.r2dt-filter-panel-title')) {
      return;
    }
    const title = document.createElement('div');
    title.className = 'r2dt-filter-panel-title';
    title.textContent = 'Leontis-Westhof Base Pairs';
    checkboxes.insertBefore(title, checkboxes.firstChild);
  }

  function removeFilterHelp() {
    const filterBtn = root.querySelector('#bpFilterBtn');
    if (filterBtn) {
      filterBtn.querySelector('#bpFilterBtnHelp, .help-icon')?.remove();
      filterBtn.setAttribute('aria-label', 'Base pairs filter and symbol legend');
    }
    root.querySelectorAll('div.help-tooltip').forEach((tip) => {
      if (/Displays checkboxes|base-pair famil/i.test(tip.textContent)) {
        tip.remove();
      }
    });
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
    root.querySelectorAll('#checkboxes .r2dt-bp-glyph-wrap').forEach((wrap) => {
      const cb = wrap.closest('label')?.querySelector('input[id^="Checkbox_"]');
      if (!cb || !isFilterCheckboxVisible(cb)) wrap.remove();
    });
    root
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
    root.querySelectorAll('input[id^="Checkbox_"]').forEach((cb) => {
      if (cb.id === 'Checkbox_All') wrapFilterLabelText(cb, 'All');
    });
  }

  function isBasePairsPanelOpen() {
    const checkboxes = root.querySelector('#checkboxes');
    if (!checkboxes) return false;
    return getComputedStyle(checkboxes).display !== 'none';
  }

  function closeBasePairsPanel() {
    if (!isBasePairsPanelOpen()) return;
    root.querySelector('#bpFilterBtn')?.click();
  }

  function isBpListPanelOpen() {
    const dialog = root.querySelector(`#bpListDialog-${PDB_LOWER}`);
    return dialog && getComputedStyle(dialog).display !== 'none';
  }

  function closeBpListPanel() {
    if (!isBpListPanelOpen()) return;
    root.querySelector(`#rnaTopologyBPList-${PDB_LOWER}`)?.click();
  }

  function ensureBpListPanelTitle() {
    const dialog = root.querySelector(`#bpListDialog-${PDB_LOWER}`);
    if (!dialog || !dialog.querySelector('ul')) return;
    if (dialog.querySelector('.r2dt-bp-list-panel-title')) return;
    const title = document.createElement('div');
    title.className = 'r2dt-bp-list-panel-title';
    title.textContent = 'Base Pair List';
    dialog.insertBefore(title, dialog.firstChild);
  }

  function applyBpListItemLabels() {
    const dialog = root.querySelector(`#bpListDialog-${PDB_LOWER}`);
    if (!dialog) return;
    dialog.querySelectorAll('ul > li').forEach((li) => {
      if (li.querySelector('.r2dt-bp-list-pair')) return;
      const raw = (li.textContent || '').trim();
      const m = raw.match(/^(.+?)\s*;\s*(\S+)\s*$/);
      if (!m) return;
      li.textContent = '';
      const pair = document.createElement('span');
      pair.className = 'r2dt-bp-list-pair';
      pair.textContent = m[1].trim();
      const family = document.createElement('span');
      family.className = 'r2dt-bp-list-family';
      family.textContent = m[2];
      li.append(pair, family);
    });
  }

  function normalizeBpListScroll() {
    const dialog = root.querySelector(`#bpListDialog-${PDB_LOWER}`);
    const ul = dialog?.querySelector('ul');
    if (!dialog || !ul) return;
    // Plugin toggles display:block, which breaks flex and prevents ul scroll.
    if (dialog.style.display !== 'none') {
      dialog.style.display = 'flex';
    }
    ul.style.maxHeight = '';
    ul.style.overflowY = '';
    ul.style.padding = '';
    ul.style.listStyle = '';
    if (!dialog.dataset.r2dtWheelBound) {
      dialog.dataset.r2dtWheelBound = '1';
      dialog.addEventListener('wheel', (ev) => { ev.stopPropagation(); }, { passive: true });
    }
  }

  function setupBpListToolbar() {
    const btn = root.querySelector(`#rnaTopologyBPList-${PDB_LOWER}`);
    const dialog = root.querySelector(`#bpListDialog-${PDB_LOWER}`);
    const pairsGroup = root.querySelector('.r2dt-toolbar-group--pairs');
    const viewer = panel2d;
    if (!btn || !dialog || !pairsGroup || !viewer || dialog.dataset.r2dtToolbarMoved) {
      return !!dialog?.dataset.r2dtToolbarMoved;
    }

    btn.textContent = 'List';
    btn.classList.add('r2dt-btn');
    btn.setAttribute('aria-label', 'Base pair list');
    btn.setAttribute('aria-haspopup', 'true');

    const wrap = document.createElement('div');
    wrap.className = 'r2dt-bp-list-dropdown';
    pairsGroup.append(wrap);
    wrap.append(btn);
    // Dock the list inside the 2D panel (not over the 3D viewer).
    viewer.appendChild(dialog);
    dialog.classList.add('r2dt-bp-list-panel');
    dialog.dataset.r2dtToolbarMoved = '1';

    const btnGroup = root.querySelector('.pdb-rna-view-btn-group.left');
    if (btnGroup && !btnGroup.querySelector('button, .pdb-rna-view-btn')) {
      btnGroup.remove();
    }

    if (!btn.dataset.r2dtListBound) {
      btn.dataset.r2dtListBound = '1';
      btn.addEventListener(
        'click',
        () => {
          setTimeout(() => {
            closeBasePairsPanel();
            btn.setAttribute('aria-expanded', String(isBpListPanelOpen()));
          }, 0);
        },
        true
      );
    }

    if (!root.querySelector('#bpFilterBtn')?.dataset.r2dtListBound) {
      const filterBtn = root.querySelector('#bpFilterBtn');
      if (filterBtn) {
        filterBtn.dataset.r2dtListBound = '1';
        filterBtn.addEventListener(
          'click',
          () => { setTimeout(closeBpListPanel, 0); },
          true
        );
      }
    }

    const ui = rnaPlugin.uiTemplateService;
    if (ui?.renderBpListDialog && !ui._r2dtBpListPatched) {
      ui._r2dtBpListPatched = true;
      const orig = ui.renderBpListDialog.bind(ui);
      ui.renderBpListDialog = function r2dtRenderBpListDialog(toggle) {
        orig(toggle);
        ensureBpListPanelTitle();
        applyBpListItemLabels();
        normalizeBpListScroll();
        const listBtn = root.querySelector(`#rnaTopologyBPList-${PDB_LOWER}`);
        if (listBtn) {
          listBtn.setAttribute('aria-expanded', String(isBpListPanelOpen()));
        }
      };
    }

    return true;
  }

  function bindToolbarDropdowns() {
    if (root.dataset.r2dtDropdownBound) return;
    root.dataset.r2dtDropdownBound = '1';

    // Capture phase: the plugin stops propagation on diagram clicks, so a
    // bubble-phase listener on document never runs for clicks on the 2D canvas.
    // Defer closing so we do not run a synthetic filter click before other
    // capture handlers (e.g. zoom buttons) on the same event.
    document.addEventListener(
      'click',
      (ev) => {
        if (isBasePairsPanelOpen()) {
          if (!ev.target.closest('#mainMenu .menu-dropdown')) {
            setTimeout(closeBasePairsPanel, 0);
          }
        }
        if (isBpListPanelOpen()) {
          if (!ev.target.closest('.r2dt-bp-list-dropdown, .r2dt-bp-list-panel')) {
            setTimeout(closeBpListPanel, 0);
          }
        }
      },
      true
    );
  }

  function mountFilterLegend() {
    root.querySelector('.r2dt-symbols-dropdown')?.remove();

    const checkboxes = root.querySelector('#checkboxes');
    if (!checkboxes || checkboxes.querySelector('.r2dt-bp-legend-more')) return true;

    const panelContent = root.querySelector('.r2dt-bp-legend-panel');
    if (!panelContent) return false;

    const details = document.createElement('details');
    details.className = 'r2dt-bp-legend-more';
    const summary = document.createElement('summary');
    summary.textContent = 'More info';
    details.append(summary, panelContent);
    checkboxes.appendChild(details);
    return true;
  }

  function bindFilterBadgeListener() {
    const checkboxes = root.querySelector('#checkboxes');
    if (!checkboxes || checkboxes.dataset.r2dtBadgeBound) return;
    checkboxes.dataset.r2dtBadgeBound = '1';
    checkboxes.addEventListener('change', updateFilterBadge);
  }

  function setupToolbar() {
    const mainMenu = root.querySelector('#mainMenu');
    if (!mainMenu) return false;
    if (mainMenu.dataset.r2dtToolbarBound) {
      bindFilterBadgeListener();
      updateFilterBadge();
      removeFilterHelp();
      ensureFilterPanelTitle();
      injectFilterGlyphs();
      mountFilterLegend();
      setupBpListToolbar();
      bindToolbarDropdowns();
      return true;
    }

    const form = mainMenu.querySelector('form');
    if (!form) return false;

    const menuSelect = form.querySelector('.menuSelectbox');
    const filterDropdown = form.querySelector('.menu-dropdown');
    const nestedWrap = root.querySelector('#nestedBP')?.parentElement;
    if (!menuSelect || !filterDropdown) return false;

    mainMenu.classList.add('r2dt-toolbar');
    mainMenu.dataset.r2dtToolbarBound = '1';
    menuSelect.style.display = 'none';

    const titleOverlay = root.querySelector('.pdb-rna-view-title');
    if (titleOverlay) titleOverlay.style.display = 'none';

    const titleText = CHAIN_ID
      ? `${STRUCTURE_ID} · chain ${CHAIN_ID}`
      : STRUCTURE_ID;
    const titleEl = document.createElement('span');
    titleEl.className = 'r2dt-toolbar-title';
    titleEl.textContent = titleText;
    titleEl.title = titleText;

    const viewGroup = document.createElement('div');
    viewGroup.className = 'r2dt-toolbar-group r2dt-toolbar-group--view';
    const viewLabel = document.createElement('span');
    viewLabel.className = 'r2dt-toolbar-group-label';
    viewLabel.textContent = 'View';

    const backboneLabel = document.createElement('label');
    backboneLabel.className = 'r2dt-toggle';
    const backboneId = `r2dt-backbone-toggle-${PDB_LOWER}`;
    backboneLabel.htmlFor = backboneId;
    backboneLabel.innerHTML =
      `<input type="checkbox" id="${backboneId}" checked ` +
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
      nestedWrap.classList.add('r2dt-nested-wrap');
      const nestedText = nestedWrap.querySelector('label span');
      if (nestedText) nestedText.textContent = 'Nested only';
      const nestedInput = nestedWrap.querySelector('#nestedBP');
      const nestedLabel = nestedWrap.querySelector('label');
      if (nestedInput) {
        nestedInput.setAttribute('aria-label', 'Only nested base pairs');
      }
      if (nestedLabel) {
        nestedLabel.classList.add('r2dt-nested-toggle');
      }
      // Duplicate #nestedBP ids on compare pages break label[for] (browser
      // toggles the first checkbox in the document). Scope clicks locally.
      if (nestedLabel && nestedInput) {
        nestedLabel.removeAttribute('for');
        nestedLabel.addEventListener('click', (ev) => {
          if (ev.target === nestedInput) return;
          ev.preventDefault();
          nestedInput.click();
        });
      }
      pairsGroup.append(nestedWrap);
    }

    form.insertBefore(titleEl, form.firstChild);
    form.insertBefore(viewGroup, menuSelect);
    form.insertBefore(pairsGroup, menuSelect);

    const filterBtn = root.querySelector('#bpFilterBtn');
    if (filterBtn) {
      filterBtn.childNodes.forEach((node) => {
        if (node.nodeType === 3 && /Filter/i.test(node.nodeValue)) {
          node.nodeValue = 'Base Pairs';
        }
      });
      if (!root.querySelector('#r2dt-filter-badge')) {
        const badge = document.createElement('span');
        badge.id = 'r2dt-filter-badge';
        badge.className = 'r2dt-filter-badge';
        badge.setAttribute('aria-hidden', 'true');
        const icon = root.querySelector('#bpFilterBtnIcon');
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

    setupBpListToolbar();
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
      const boxes = root.querySelectorAll('input[id^="Checkbox_"]');
      if (boxes.length === 0) {
        if (attempts++ > 40) return;
        setTimeout(tick, 100);
        return;
      }
      boxes.forEach((cb) => {
        if (cb.id === 'Checkbox_All') return;
        const family = cb.id.slice('Checkbox_'.length);
        const td = cb.closest('td');
        if (!presentBPs.has(family)) {
          if (td) td.style.display = 'none';
        } else if (td) {
          td.style.display = '';
        }
      });
      hideEmptyFilterRows();
      const all = root.querySelector('#Checkbox_All');
      if (all && !all.checked) all.click();
      ensureFilterPanelTitle();
      injectFilterGlyphs();
      mountFilterLegend();
      updateFilterBadge();
    };
    tick();
  })();

  // Relabel handled in setupBpListToolbar (short "List" label in toolbar).

  // The plugin injects a "?" help icon and floating tooltip on the filter
  // button; remove them once the base-pairs panel is wired up.
  (function stripFilterHelp() {
    let attempts = 0;
    const tick = () => {
      removeFilterHelp();
      const helpLeft = root.querySelector('#bpFilterBtnHelp');
      const tipLeft = [...root.querySelectorAll('div.help-tooltip')].some(
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
  if (unobserved.length || crossingPairs.size) {
    let attempts = 0;
    const tick = () => {
      if (applyFixups() || attempts++ > 40) return;
      setTimeout(tick, 100);
    };
    tick();
  }

  // --- 3D viewer (optional; compare pages use a shared molstar pane) ---
  let molstar = null;
  if (link3d && panel3d && !panel3d.hidden) {
    molstar = await renderMolstarPlugin(panel3d, STRUCTURE_URL, STRUCTURE_FORMAT);
  }

  ctx.handles = {
    molstar,
    rnaPlugin,
    apiData,
    fr3dData,
    labelToAuth,
    authToLabel,
    root,
  };
  if (link3d) window.__r2dt = ctx.handles;

  const selectInMolstar = createMolstarSelector(molstar, labelToAuth, CHAIN_ID);

  if (link3d && molstar) {
    document.addEventListener('PDB.RNA.viewer.click', (e) => {
      selectInMolstar(labelsFromEvent(e));
    });
  }

  // Set by _renderLBN when the panel is shown; keeps LBN in sync with every
  // base-pair selection path (list, 2D line click, API).
  let lbnHighlightFn = null;

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
    if (lbnHighlightFn) lbnHighlightFn([a, b]);
  }

  // Find a rendered base-pair path by its two seq ids (either order).
  // The bp class is the last token, "<bp>_<a>_<b>", so an attribute
  // "ends-with" match avoids false hits like _5_27 inside _5_271.
  function findBPPath(a, b) {
    return root.querySelector(
      `.rnaviewBP[class$="_${a}_${b}"], .rnaviewBP[class$="_${b}_${a}"]`
    );
  }

  // Base-pair lines in 2D don't emit PDB.RNA.viewer.click; wire our own
  // click handler that picks both partners out of the class name
  // (e.g. "cWW_5_27" -> labels [5, 27]).
  function attachBPClicks() {
    root.querySelectorAll('path[class*="rnaviewBP"]').forEach((el) => {
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
  onBpFilterUpdated = () => {
    attachBPClicks();
    applyFixups();
  };

  // Base Pairings List rows: the plugin handles a row click by dispatching
  // a click onto the matching SVG path, which silently does nothing when
  // that path isn't currently rendered (e.g. nested pairs in the non-nested
  // view). Handle row clicks directly so every listed pair updates 3D.
  const bpListId = 'bpListDialog-' + STRUCTURE_ID.toLowerCase();
  document.addEventListener('click', (ev) => {
    const li = ev.target.closest?.('#' + bpListId + ' li');
    if (!li || !root.contains(li)) return;
    const pairText =
      li.querySelector('.r2dt-bp-list-pair')?.textContent || li.textContent || '';
    const m = pairText.match(/(\d+)\D*?-\D*?(\d+)/);
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
  const bpContainer = panel2d;
  if (bpContainer) bpObserver.observe(bpContainer, { childList: true, subtree: true });
  // 2D hover intentionally does not update the 3D view -- only clicks
  // on nucleotides or base-pair lines cause a 3D selection/focus.

  if (link3d && molstar) {
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
  }

  // ── LBN (layered dot-bracket notation) panel ────────────────────────────
  if (showLbn) {
    let lbnData = null;
    try {
      const resp = await fetch(resolveUrl(baseUrl, 'lbn.json'));
      if (resp.ok) lbnData = await resp.json();
    } catch (_) { /* ignore */ }

    if (lbnData && lbnData.rows && lbnData.rows.length > 0) {
      const lbnPanel = root.querySelector('#lbn-panel');
      const lbnBody = root.querySelector('#lbn-body');
      if (lbnPanel && lbnBody) {
        lbnPanel.hidden = false;
        _renderLBN(lbnData, lbnBody);
      }
    }
  }

  function _renderLBN(data, container) {
    const SEQ = data.sequence;
    const N   = SEQ.length;

    let html = '';
    html += '<div class="lbn-row"><span class="lbn-label">seq</span>: ';
    for (let i = 0; i < N; i++) {
      html += `<span data-pos="${i + 1}" class="lbn-nt">${SEQ[i]}</span>`;
    }
    html += '</div>';

    for (const row of data.rows) {
      html += `<div class="lbn-row"><span class="lbn-label">${row.label}</span>: `;
      for (let i = 0; i < N; i++) {
        const pos = i + 1;
        const ch  = row.chars[i];
        if (ch === '.') {
          html += '<span class="lbn-dot">.</span>';
        } else {
          const partner = row.partners[String(pos)];
          const pAttr   = partner != null ? ` data-partner="${partner}"` : '';
          html += `<span data-pos="${pos}"${pAttr} class="lbn-bp">${ch}</span>`;
        }
      }
      html += '</div>';
    }

    container.innerHTML = html;

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

      if (highlighted.length > 0) {
        highlighted[0].scrollIntoView({
          behavior: 'smooth',
          block: 'nearest',
          inline: 'nearest',
        });
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

    lbnHighlightFn = _lbnHighlight;

    // Expose for console debugging.
    window.__r2dt.lbnData      = data;
    window.__r2dt.lbnHighlight = _lbnHighlight;
  }
  // ── end LBN ─────────────────────────────────────────────────────────────

  ctx.handles.selectResidue = (l) => selectInMolstar([l]);
  ctx.handles.selectBasePair = (a, b) => selectBasePair(a, b, findBPPath(a, b));

  }

  function buildCompare2dDom(slotEl, panelOpts) {
    const root = document.createElement('div');
    root.className = 'r2dt-viewer-root';
    root.dataset.layout = '2d-only';
    if (panelOpts.height != null) {
      root.style.setProperty(
        '--r2dt-viewer-height',
        typeof panelOpts.height === 'number' ? `${panelOpts.height}px` : String(panelOpts.height)
      );
    }
    const vis = document.createElement('div');
    vis.className = 'r2dt-viewer-vis';
    const panel2d = document.createElement('div');
    panel2d.className = 'r2dt-panel r2dt-panel--2d';
    vis.append(panel2d);
    root.append(vis);
    slotEl.appendChild(root);
    return { root, panel2d, panel3d: null };
  }

  function buildCompareSlot(title, subtitle) {
    const slot = document.createElement('div');
    slot.className = 'r2dt-compare-slot';
    const heading = document.createElement('h2');
    heading.className = 'r2dt-compare-slot-title';
    heading.textContent = title;
    if (subtitle) {
      const tag = document.createElement('span');
      tag.className = 'r2dt-compare-tag';
      tag.textContent = subtitle;
      heading.appendChild(document.createTextNode(' '));
      heading.appendChild(tag);
    }
    slot.appendChild(heading);
    return slot;
  }

  async function create(userOptions) {
    if (activeViewer) {
      throw new Error('R2DTViewer.create: only one viewer or compare widget per page (v1)');
    }

    const opts = userOptions || {};
    const mountEl = resolveMount(opts.mount);
    const panelData = await resolvePanelData(opts.baseUrl || '.', opts);

    const dom = buildViewerDom(mountEl, {
      structureId: panelData.structureId,
      layout: opts.layout || 'side-by-side',
      height: opts.height,
      panelWidth: opts.panelWidth,
      showLbn: opts.showLbn !== false,
      showLegend: opts.showLegend !== false,
    });

    const ctx = {
      root: dom.root,
      panel2d: dom.panel2d,
      panel3d: dom.panel3d,
      baseUrl: panelData.baseUrl,
      resolveUrl: panelData.resolveUrl,
      STRUCTURE_ID: panelData.structureId,
      CHAIN_ID: panelData.chainId,
      STRUCTURE_URL: panelData.structureUrl,
      STRUCTURE_FORMAT: panelData.structureFormat,
      PDB_LOWER: panelData.PDB_LOWER,
      apiData: panelData.apiData,
      fr3dData: panelData.fr3dData,
      showLbn: opts.showLbn !== false,
      link3d: true,
      cleanup: [],
      handles: null,
    };

    await initViewer(ctx);

    const handle = {
      root: ctx.root,
      destroy() {
        ctx.cleanup.forEach((fn) => { try { fn(); } catch (_) {} });
        mountEl.innerHTML = '';
        if (activeViewer === handle) activeViewer = null;
        if (window.__r2dt === ctx.handles) window.__r2dt = undefined;
      },
      selectResidue(label) {
        return ctx.handles?.selectResidue?.(label);
      },
      selectBasePair(a, b) {
        return ctx.handles?.selectBasePair?.(a, b);
      },
    };
    activeViewer = handle;
    if (typeof opts.onReady === 'function') opts.onReady(handle);
    return handle;
  }

  async function createCompare(userOptions) {
    if (activeViewer) {
      throw new Error('R2DTViewer.createCompare: only one viewer or compare widget per page (v1)');
    }

    const opts = userOptions || {};
    const mountEl = resolveMount(opts.mount);
    const panels = opts.panels;
    if (!Array.isArray(panels) || panels.length < 2) {
      throw new Error('R2DTViewer.createCompare: panels array required (min 2 entries)');
    }

    const molstarOpts = opts.molstar || {};
    const linkIndex = molstarOpts.panelIndex != null ? molstarOpts.panelIndex : 0;
    if (linkIndex < 0 || linkIndex >= panels.length) {
      throw new Error('R2DTViewer.createCompare: molstar.panelIndex out of range');
    }

    installMultiPanelDomShim();

    if (opts.fetchShim !== false) {
      const routes = [];
      for (let i = 0; i < panels.length; i++) {
        const cfg = await resolvePanelData(panels[i].baseUrl || '.', panels[i]);
        routes.push({ pdbIdMatch: cfg.PDB_LOWER, baseUrl: cfg.baseUrl });
      }
      installFetchShim(routes);
    }

    mountEl.innerHTML = '';
    const compareRoot = document.createElement('div');
    compareRoot.className = 'r2dt-compare-root';
    if (opts.panelHeight != null) {
      compareRoot.style.setProperty(
        '--r2dt-compare-panel-height',
        typeof opts.panelHeight === 'number' ? `${opts.panelHeight}px` : String(opts.panelHeight)
      );
    }
    const grid = document.createElement('div');
    grid.className = 'r2dt-compare-grid';
    compareRoot.appendChild(grid);
    mountEl.appendChild(compareRoot);

    const cleanup = [];
    const panelCtxs = [];
    const panelHandles = [];

    for (let i = 0; i < panels.length; i++) {
      const pOpts = panels[i];
      const panelData = await resolvePanelData(pOpts.baseUrl || '.', pOpts);
      const slot = buildCompareSlot(pOpts.title || panelData.structureId, pOpts.subtitle || '');
      grid.appendChild(slot);
      const dom = buildCompare2dDom(slot, { height: opts.panelHeight });
      const ctx = {
        root: dom.root,
        panel2d: dom.panel2d,
        panel3d: dom.panel3d,
        baseUrl: panelData.baseUrl,
        resolveUrl: panelData.resolveUrl,
        STRUCTURE_ID: panelData.structureId,
        CHAIN_ID: panelData.chainId,
        STRUCTURE_URL: panelData.structureUrl,
        STRUCTURE_FORMAT: panelData.structureFormat,
        PDB_LOWER: panelData.PDB_LOWER,
        apiData: panelData.apiData,
        fr3dData: panelData.fr3dData,
        showLbn: false,
        link3d: false,
        cleanup: [],
        handles: null,
      };
      await initViewer(ctx);
      panelCtxs.push(ctx);
      panelHandles.push({
        root: ctx.root,
        structureId: panelData.structureId,
        selectResidue(label) {
          return ctx.handles?.selectResidue?.(label);
        },
        selectBasePair(a, b) {
          return ctx.handles?.selectBasePair?.(a, b);
        },
      });
    }

    const linkedCtx = panelCtxs[linkIndex];
    const linkedPanel = panels[linkIndex];
    const molBaseUrl = molstarOpts.baseUrl || linkedPanel.baseUrl || '.';
    const molData = await resolvePanelData(molBaseUrl, {
      ...linkedPanel,
      ...molstarOpts,
      structureId: molstarOpts.structureId || linkedPanel.structureId,
      chainId: molstarOpts.chainId ?? linkedPanel.chainId,
    });

    const molSlot = buildCompareSlot(
      molstarOpts.title || molData.structureId,
      molstarOpts.subtitle || '— 3D'
    );
    grid.appendChild(molSlot);
    const molRoot = document.createElement('div');
    molRoot.className = 'r2dt-viewer-root r2dt-compare-molstar-root';
    molRoot.dataset.layout = '3d-only';
    const molVis = document.createElement('div');
    molVis.className = 'r2dt-viewer-vis';
    const panel3d = document.createElement('div');
    panel3d.className = 'r2dt-panel r2dt-panel--3d';
    molVis.append(panel3d);
    molRoot.append(molVis);
    molSlot.appendChild(molRoot);

    const molstar = await renderMolstarPlugin(
      panel3d,
      molData.structureUrl,
      molData.structureFormat
    );
    const selectInMolstar = createMolstarSelector(
      molstar,
      linkedCtx.handles.labelToAuth,
      molData.chainId
    );
    const linkedPdbLower = linkedCtx.PDB_LOWER;

    function onLinked2dClick(e) {
      const d = e.eventData || e.detail || {};
      if ((d.pdbId || '').toLowerCase() !== linkedPdbLower) return;
      selectInMolstar(labelsFromEvent(e));
    }
    document.addEventListener('PDB.RNA.viewer.click', onLinked2dClick);
    cleanup.push(() => document.removeEventListener('PDB.RNA.viewer.click', onLinked2dClick));

    const handle = {
      root: compareRoot,
      panels: panelHandles,
      molstar,
      destroy() {
        cleanup.forEach((fn) => { try { fn(); } catch (_) {} });
        mountEl.innerHTML = '';
        if (activeViewer === handle) activeViewer = null;
      },
      selectResidue(label) {
        return selectInMolstar([label]);
      },
      selectBasePair(a, b) {
        return linkedCtx.handles?.selectBasePair?.(a, b);
      },
    };
    activeViewer = handle;
    if (typeof opts.onReady === 'function') opts.onReady(handle);
    return handle;
  }

  global.R2DTViewer = { create, createCompare };
})(typeof window !== 'undefined' ? window : global);
