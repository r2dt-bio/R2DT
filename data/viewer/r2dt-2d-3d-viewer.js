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
      if (!resp.ok) return null;
      const ct = resp.headers.get('content-type') || '';
      if (!ct.includes('json')) return null;
      return await resp.json();
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

  function patchAutoZoomOnSelect() {
    const svc = window.UiActionsService;
    if (!svc) return;
    // Selection must not pan/zoom the diagram (plugin default zooms to nt).
    svc.zoomToNucleotides = function () {};
  }

  function enhanceSelectedNucleotideStyleOnce() {
    const svc = window.UiActionsService;
    if (!svc || !svc.colorNucleotide || svc.__r2dtStyleEnhanced) return;
    svc.__r2dtStyleEnhanced = true;

    const FONT_BUMP_PX = 1.33;
    const PARTNER_FONT_BUMP_PX = 0.67;
    const HALO_STROKE = '#ffffff';
    const HALO_WIDTH = '2';
    const PARTNER_HALO_WIDTH = '1.25';
    const BG_FILL = '#fff3b0';
    const PARTNER_BG_FILL = '#ffe8c2';
    const BG_PAD = 2;
    const PARTNER_BG_PAD = 1.5;
    const origStyle = new Map();
    const selectionBgs = new Map();
    svc.__r2dtPartnerLabels = svc.__r2dtPartnerLabels || new Set();

    function nucleotideEl(pdbId, label) {
      return document.getElementsByClassName(
        `rnaviewEle rnaviewEle_${pdbId} rnaview_${pdbId}_${label}`
      )[0];
    }

    function setOrRemove(el, attr, val) {
      if (val != null) el.setAttribute(attr, val);
      else el.removeAttribute(attr);
    }

    // Compare mode shares this whole module (and these Maps) across every
    // panel, but residue labels repeat across structures (co-indexed
    // panels both have a "residue 1") -- key every per-nucleotide record by
    // pdbId+label, not label alone, or one panel's badge/typography gets
    // silently deleted as a side effect of the other panel's mirrored click.
    function mapKey(pdbId, label) {
      return `${pdbId}::${label}`;
    }

    function removeSelectionBg(key) {
      const existing = selectionBgs.get(key);
      if (existing && existing.parentNode) existing.parentNode.removeChild(existing);
      selectionBgs.delete(key);
    }

    function applySelectionBg(el, key, fill, pad) {
      if (!el || el.nodeName !== 'text' || !el.parentNode) return;
      removeSelectionBg(key);
      let bbox;
      try {
        bbox = el.getBBox();
      } catch (_) {
        return;
      }
      const rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
      rect.setAttribute('class', 'r2dt-nt-selection-bg');
      rect.setAttribute('x', String(bbox.x - pad));
      rect.setAttribute('y', String(bbox.y - pad));
      rect.setAttribute('width', String(bbox.width + 2 * pad));
      rect.setAttribute('height', String(bbox.height + 2 * pad));
      rect.setAttribute('rx', '3');
      rect.setAttribute('fill', fill);
      rect.setAttribute('stroke', 'none');
      rect.style.pointerEvents = 'none';
      el.parentNode.insertBefore(rect, el);
      selectionBgs.set(key, rect);
    }

    function applySelectionTypography(el, pdbId, label, color) {
      if (!el || el.nodeName !== 'text') return;
      const isPartner = svc.__r2dtPartnerLabels && svc.__r2dtPartnerLabels.has(label);
      const key = mapKey(pdbId, label);
      if (!origStyle.has(key)) {
        origStyle.set(key, {
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
        el.setAttribute(
          'font-size',
          (fs + (isPartner ? PARTNER_FONT_BUMP_PX : FONT_BUMP_PX)) + 'px'
        );
      }
      el.setAttribute('stroke', HALO_STROKE);
      el.setAttribute('stroke-width', isPartner ? PARTNER_HALO_WIDTH : HALO_WIDTH);
      el.setAttribute('stroke-linejoin', 'round');
      el.setAttribute('paint-order', 'stroke fill');
      // Compare mode passes an `rgb(...)` selection colour (per-structure
      // green/blue) -- tint the badge to match. The standalone viewer's
      // 'orange'/'#ffab40' fallbacks keep their original neutral badge.
      const tint = /^rgb\(/.test(color || '')
        ? paleTintFor(color, isPartner ? 0.88 : 0.85)
        : null;
      applySelectionBg(
        el,
        key,
        tint || (isPartner ? PARTNER_BG_FILL : BG_FILL),
        isPartner ? PARTNER_BG_PAD : BG_PAD
      );
    }

    function restoreTypography(pdbId, label) {
      const key = mapKey(pdbId, label);
      const stored = origStyle.get(key);
      removeSelectionBg(key);
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
      origStyle.delete(key);
    }

    const origColor = svc.colorNucleotide.bind(svc);
    svc.colorNucleotide = function (pdbId, label, color, mode) {
      origColor(pdbId, label, color, mode);
      if (mode === 'selection') {
        applySelectionTypography(nucleotideEl(pdbId, label), pdbId, label, color);
      }
    };

    const origClear = svc.clearNucleotides.bind(svc);
    svc.clearNucleotides = function (pdbId, mode, labels) {
      if (mode === 'selection') {
        // Don't enumerate "everything selected" from svc.selected -- in
        // compare mode every panel shares that one Map, keyed by label only
        // (not pdbId). When the *other* panel clears its own selection over
        // the same co-indexed labels, the vendor unconditionally deletes
        // each entry even if it found no matching DOM node under its own
        // pdbId, silently dropping this panel's still-badged labels from
        // the map before this panel ever gets a chance to restore them
        // (orphaning their selection highlight). Our own origStyle records
        // are correctly scoped by pdbId, so enumerate from those instead.
        const prefix = `${pdbId}::`;
        const keys = labels != null
          ? labels
          : Array.from(origStyle.keys())
            .filter((k) => k.startsWith(prefix))
            .map((k) => k.slice(prefix.length));
        keys.forEach((label) => restoreTypography(pdbId, label));
      }
      origClear(pdbId, mode, labels);
    };

    // Hover tooltip reads "Residue N" by default; relabel it to
    // "Highlighted Residue N" so it reads clearly next to the persistent
    // "Selected Residue N" selection tooltip.
    if (svc.highlightNucleotide && !svc.__r2dtHighlightLabelPatched) {
      svc.__r2dtHighlightLabelPatched = true;
      const origHighlight = svc.highlightNucleotide.bind(svc);
      svc.highlightNucleotide = function (pdbId) {
        origHighlight.apply(svc, arguments);
        const tip = document.getElementById(`${pdbId}-rnaTopologyTooltipHighlight`);
        const strong = tip && tip.querySelector('strong');
        if (strong && !/Highlighted/.test(strong.textContent)) {
          strong.textContent = strong.textContent.replace(
            /\bResidue/, 'Highlighted Residue'
          );
        }
      };
    }
  }

  function labelsFromEvent(e) {
    const d = e.eventData || e.detail || {};
    if (Array.isArray(d.label_seq_ids)) return d.label_seq_ids;
    if (d.label_seq_id !== undefined && d.label_seq_id !== null) return [d.label_seq_id];
    return [];
  }

  function buildLabelMaps(apiData) {
    const labelToAuth = {};
    const labelToChain = {};
    const chainIds = apiData.chain_ids;
    apiData.label_seq_ids.forEach((label, i) => {
      if (label !== null && label !== undefined) {
        labelToAuth[label] = apiData.auth_seq_ids[i];
        if (chainIds && chainIds[i] != null) labelToChain[label] = chainIds[i];
      }
    });
    const authToLabel = {};
    Object.entries(labelToAuth).forEach(([label, auth]) => {
      authToLabel[auth] = parseInt(label);
    });
    return { labelToAuth, authToLabel, labelToChain };
  }

  async function resolvePanelData(baseUrl, opts) {
    const normalized = normalizeBaseUrl(baseUrl);
    const manifest = opts.structureId ? null : await loadManifest(normalized);
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
    // otherPairKeys is only present in compare mode (bpCompare set at all);
    // fetch it lazily from this panel's own folder when not passed inline.
    let bpCompare = opts.bpCompare || null;
    if (bpCompare && !bpCompare.otherPairKeys) {
      try {
        const otherPairKeys = await fetchJson(resolveUrl(normalized, 'bp-compare.json'));
        bpCompare = { ...bpCompare, otherPairKeys };
      } catch (_) { /* optional */ }
    }
    return {
      baseUrl: normalized,
      resolveUrl,
      structureId,
      chainId,
      structureUrl: resolvedStructureUrl,
      structureFormat,
      apiData,
      fr3dData,
      lbnData: opts.lbnData || null,
      bpCompare,
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

    // Init and programmatic cross-panel ops (mirroring) happen without a
    // pointer event, so let callers pin the active panel explicitly.  Returns
    // the previous root so callers can restore it.
    global.__r2dtSetActiveRoot = (root) => {
      const prev = activeRoot;
      activeRoot = root || null;
      return prev;
    };

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
      // Scope every selector (including svg.rnaTopoSvg) to the active panel so
      // two panels' identically-classed SVGs don't collide on the first match.
      if (activeRoot) {
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
    keepFogAndOcclusionOff(molstar);
    return molstar;
  }

  // Mol*'s auto camera-clipping ties fog density to the current view radius,
  // and ambient occlusion's shadowing gets a lot more visible too, so both read
  // as "switching on" whenever focusLoci zooms tight on a single nucleotide.
  // cameraFog has a real discriminated "off" state and setting it once here
  // reliably sticks. postprocessing.occlusion and cameraClipping.radius do not
  // -- occlusion gets reset back to Mol*'s default at some later,
  // unpredictable point (confirmed at least once, right after the compare
  // view's overlay structure loads, but re-asserting right after that specific
  // load isn't enough either, so that isn't the only reset), and
  // cameraClipping.radius looks to be recomputed continuously from the camera
  // view every frame rather than a static setting, so nothing short of a
  // per-frame override would pin it. Re-asserting all three on an interval is
  // a best-effort, not a confirmed fix, for the latter two -- see
  // .ai/molstar-fog-occlusion-clipping.md for the full investigation before
  // spending more time on this.
  function keepFogAndOcclusionOff(molstar) {
    const apply = () => {
      if (molstar && molstar.plugin && molstar.plugin.canvas3d) {
        molstar.plugin.canvas3d.setProps({
          cameraFog: { name: 'off', params: {} },
          postprocessing: { occlusion: { name: 'off', params: {} } },
          cameraClipping: { radius: 0 },
        });
      }
    };
    apply();
    setInterval(apply, 1000);
  }

  function createMolstarSelector(molstar, labelToAuth, chainId, labelToChain) {
    return async function selectInMolstar(labels) {
      if (!molstar) return;
      if (!labels || labels.length === 0) {
        try {
          await molstar.visual.select({ data: [], keepRepresentations: true });
        } catch (_) { /* best-effort */ }
        return;
      }
      const data = labels.map((l) => {
        const auth = labelToAuth[l];
        if (auth === undefined || auth === null) return null;
        // Per-nucleotide chain (multi-chain): fall back to the single chainId.
        const chain = (labelToChain && labelToChain[l]) || chainId;
        return {
          auth_asym_id: chain,
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

  // Fallback highlight (deep pink) for any structure without its own base colour.
  const SELECTION_COLOR = { r: 233, g: 30, b: 99 };

  // Per-structure highlight colour: a vivid, brightened version of the
  // structure's own base colour, so a highlighted reference nt reads as bright
  // green (base is CASP green) and a highlighted model nt as bright blue (base
  // is CASP navy). Same hue as the structure it belongs to → the colour itself
  // says "this is the reference" vs "this is the model", while staying clearly
  // distinct from the (darker, less saturated) base colour around it.
  function rgbToHsl(r, g, b) {
    r /= 255; g /= 255; b /= 255;
    const max = Math.max(r, g, b), min = Math.min(r, g, b);
    const l = (max + min) / 2;
    let h = 0, s = 0;
    const d = max - min;
    if (d !== 0) {
      s = d / (1 - Math.abs(2 * l - 1));
      if (max === r) h = ((g - b) / d) % 6;
      else if (max === g) h = (b - r) / d + 2;
      else h = (r - g) / d + 4;
      h *= 60;
      if (h < 0) h += 360;
    }
    return { h, s, l };
  }
  function hslToRgb(h, s, l) {
    const c = (1 - Math.abs(2 * l - 1)) * s;
    const x = c * (1 - Math.abs(((h / 60) % 2) - 1));
    const m = l - c / 2;
    let r = 0, g = 0, b = 0;
    if (h < 60) [r, g, b] = [c, x, 0];
    else if (h < 120) [r, g, b] = [x, c, 0];
    else if (h < 180) [r, g, b] = [0, c, x];
    else if (h < 240) [r, g, b] = [0, x, c];
    else if (h < 300) [r, g, b] = [x, 0, c];
    else [r, g, b] = [c, 0, x];
    return {
      r: Math.round((r + m) * 255),
      g: Math.round((g + m) * 255),
      b: Math.round((b + m) * 255),
    };
  }
  function nudgedHueFor(base) {
    let { h } = rgbToHsl(base.r, base.g, base.b);
    // The reference base is a spring-green (~150°); a same-hue highlight reads
    // too close to it, so nudge green-family hues toward lime (~90°) — clearly
    // brighter and more yellow, unmistakable against the base green.
    if (h > 100 && h < 175) h = 90;
    return h;
  }
  function highlightColorFor(base) {
    if (!base) return SELECTION_COLOR;
    return hslToRgb(nudgedHueFor(base), 0.95, 0.6); // vivid + bright
  }
  // Companion shade for a paired/partner residue: same hue, lighter and
  // softer — mirrors how the old orange/`#ffab40` pair related to each other.
  function partnerHighlightColorFor(base) {
    if (!base) return SELECTION_COLOR;
    return hslToRgb(nudgedHueFor(base), 0.85, 0.72);
  }
  function rgbToCss(c) {
    return `rgb(${c.r},${c.g},${c.b})`;
  }
  // Parses the small set of colour formats this file ever hands to
  // colorNucleotide: the literal 'orange'/'#ffab40' fallbacks (standalone
  // viewer) or an `rgb(r,g,b)` string (compare-mode per-structure colour).
  function parseCssColor(css) {
    if (!css) return null;
    if (css === 'orange') return { r: 255, g: 165, b: 0 };
    let m = /^rgb\((\d+),\s*(\d+),\s*(\d+)\)$/.exec(css);
    if (m) return { r: +m[1], g: +m[2], b: +m[3] };
    m = /^#([0-9a-fA-F]{6})$/.exec(css);
    if (m) {
      const n = parseInt(m[1], 16);
      return { r: (n >> 16) & 255, g: (n >> 8) & 255, b: n & 255 };
    }
    return null;
  }
  // A very pale, same-hue tint for the selection badge behind a clicked
  // letter, so the badge reads as "this colour, softly" rather than a
  // colour-neutral yellow that clashes with a green/blue selection.
  function paleTintFor(css, lightness) {
    const rgb = parseCssColor(css);
    if (!rgb) return null;
    const { h } = rgbToHsl(rgb.r, rgb.g, rgb.b);
    return rgbToCss(hslToRgb(h, 0.9, lightness));
  }

  // Selection fan-out for the compare view's shared Mol*, which holds two
  // superimposed structures (reference + predicted model). Each target scopes to
  // one loaded structure via `structureNumber` and carries that structure's own
  // label→(auth, chain) map, so one 2D click highlights the matching residue in
  // both. `nonSelectedColor` restores each structure's base colour on every call,
  // so the same call both highlights the click and clears the previous one.
  function createMultiMolstarSelector(molstar, targets) {
    function residueData(target, labels) {
      return labels.map((l) => {
        const auth = target.labelToAuth ? target.labelToAuth[l] : undefined;
        if (auth === undefined || auth === null) return null;
        const chain = (target.labelToChain && target.labelToChain[l]) || target.chainId;
        return {
          auth_asym_id: chain,
          start_auth_residue_number: auth,
          end_auth_residue_number: auth,
        };
      }).filter((d) => d !== null);
    }
    // Highlight the clicked residue in EVERY structure at once, each in its own
    // structure-coloured highlight (bright green on the reference, bright blue on
    // the model), so you can see both where the residue sits and which strand is
    // which. `focusStructureNumber` (optional): which structure's residue the
    // camera frames — the clicked panel's own, so the view never drifts to (or,
    // for a reference-only residue, gets stuck on) the wrong structure.
    return async function selectInMolstar(labels, focusStructureNumber) {
      if (!molstar) return;
      const cleared = !labels || labels.length === 0;
      for (let i = 0; i < targets.length; i++) {
        const t = targets[i];
        const data = cleared ? [] : residueData(t, labels);
        const hi = t.highlightColor || SELECTION_COLOR;
        const params = {
          data: data.map((d) => ({ ...d, color: hi, focus: false })),
          structureNumber: t.structureNumber,
          keepRepresentations: true,
        };
        if (t.baseColor) params.nonSelectedColor = t.baseColor;
        try {
          await molstar.visual.select(params);
        } catch (_) { /* best-effort per structure */ }
      }
      // Zoom into the clicked residue. getLociForParams resolves auth params
      // against a *specific* structure; without the structureNumber it falls
      // back to the last-loaded one (the model), so a reference-only residue
      // (its chain/number don't exist in the model) yields an empty loci and
      // the camera flies to blank space. Frame the clicked structure's own
      // residue; superimposition keeps the partner highlight nearby in view.
      // But the clicked *panel* isn't always the clicked *structure*: a
      // widened multi-chain reference can show a residue in 2D (both panels
      // share its layout) that only the reference actually has 3D atoms for
      // (e.g. a real dimer's second chain, never scored against the model).
      // If the preferred structure has nothing here, fall back to whichever
      // target does, rather than focusing an empty loci and going blank.
      if (!cleared) {
        try {
          const fnum = focusStructureNumber != null
            ? focusStructureNumber : targets[0] && targets[0].structureNumber;
          let ft = targets.find((t) => t.structureNumber === fnum) || targets[0];
          let focusData = ft ? residueData(ft, labels) : [];
          if (focusData.length === 0) {
            for (const t of targets) {
              const d = residueData(t, labels);
              if (d.length > 0) { ft = t; focusData = d; break; }
            }
          }
          const loci = ft && ft.structureNumber != null
            ? molstar.getLociForParams(focusData, ft.structureNumber)
            : molstar.getLociForParams(focusData);
          const camera = molstar.plugin
            && molstar.plugin.managers
            && molstar.plugin.managers.camera;
          if (loci && camera && camera.focusLoci) camera.focusLoci(loci);
        } catch (_) { /* best-effort */ }
      }
    };
  }

  // Reference/model visibility toggles for the shared 3D pane. Each checkbox
  // drives `structureVisibility(structureNumber, on)`; a colour swatch matches
  // the structure's base colour so the control doubles as a legend.
  function addStructureToggles(slotEl, molstar, entries) {
    const bar = document.createElement('div');
    bar.className = 'r2dt-3d-toggles';
    entries.forEach((entry) => {
      const label = document.createElement('label');
      label.className = 'r2dt-3d-toggle';
      const cb = document.createElement('input');
      cb.type = 'checkbox';
      cb.checked = true;
      cb.addEventListener('change', () => {
        try { molstar.visual.structureVisibility(entry.structureNumber, cb.checked); } catch (_) {}
      });
      const swatch = document.createElement('span');
      swatch.className = 'r2dt-3d-swatch';
      if (entry.color) {
        swatch.style.background = `rgb(${entry.color.r},${entry.color.g},${entry.color.b})`;
      }
      const text = document.createElement('span');
      text.className = 'r2dt-3d-toggle-label';
      text.textContent = entry.label;
      label.append(cb, swatch, text);
      bar.appendChild(label);
    });
    // Placed at the end of the slot so the legend/toggles sit *below* the 3D
    // pane (the 3D root has already been appended when this runs).
    bar.classList.add('r2dt-3d-toggles--below');
    slotEl.appendChild(bar);
    return bar;
  }

  const BP_GLYPH_COLOR = '#909090';

  // pdb-rna-viewer stores these families under the flipped LW code (cHW -> cWH
  // bucket, tSW -> tWS, …). Filter checkboxes and path visibility use the
  // flipped name; fr3d.json keeps FR3D's original direction.
  const FLIP_LW_FAMILIES = new Set(['cSH', 'tSH', 'cSW', 'tSW', 'cHW', 'tHW']);

  function flipLw(lw) {
    if (lw && lw.length === 3 && (lw[0] === 'c' || lw[0] === 't')) {
      return lw[0] + lw[2] + lw[1];
    }
    return lw;
  }

  function pluginLwFamily(lw) {
    return FLIP_LW_FAMILIES.has(lw) ? flipLw(lw) : lw;
  }

  function presentLwFamilies(annotations) {
    const families = new Set();
    (annotations || []).forEach((a) => {
      if (!a.bp) return;
      families.add(a.bp);
      families.add(pluginLwFamily(a.bp));
      // Dropdown families follow drawable LW codes; ncWW/ntWW map onto WW.
      if (a.bp === 'ncWW') families.add('cWW');
      if (a.bp === 'ntWW') families.add('tWW');
    });
    return families;
  }

  function bpPairKeyFromPathId(pathId) {
    const m = (pathId || '').match(/([a-zA-Z]{3})_(\d+)_(\d+)/);
    if (!m) return pathId || '';
    const a = +m[2];
    const b = +m[3];
    return `${Math.min(a, b)}_${Math.max(a, b)}`;
  }

  // [min, max] residue numbers from a "<family>_<a>_<b>" pathID, for sorting
  // the base-pair list by nucleotide position. Unparseable ids sort last.
  function bpPositionsFromPathId(pathId) {
    const m = (pathId || '').match(/([a-zA-Z]{3})_(\d+)_(\d+)/);
    if (!m) return [Infinity, Infinity];
    const a = +m[2];
    const b = +m[3];
    return a < b ? [a, b] : [b, a];
  }

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
      const lbnCaption = document.createElement('p');
      lbnCaption.className = 'r2dt-viewer-lbn-caption';
      lbnCaption.textContent =
        'Greyed-out pairs are hidden in the current 2D view '
        + '(via "Nested only" or the base-pair family filters).';
      const lbnBody = document.createElement('div');
      lbnBody.id = 'lbn-body';
      lbnBody.className = 'lbn-body';
      lbnBody.textContent = 'Loading…';
      lbnPanel.append(lbnTitle, lbnCaption, lbnBody);
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

    // Compare mode: classify each listed base pair against the other structure's
    // pairs. Reference list → TP (matched) / FN (missed); model list → TP / FP.
    const bpCompare = ctx.bpCompare || null;
    let bpCompareKeys = bpCompare ? new Set(bpCompare.otherPairKeys || []) : null;
    function classifyBpPair(a, b) {
      if (!bpCompareKeys) return null;
      const key = `${Math.min(a, b)}_${Math.max(a, b)}`;
      if (bpCompareKeys.has(key)) return 'TP';
      return bpCompare.role === 'reference' ? 'FN' : 'FP';
    }
    function setOtherPairKeys(keys) {
      if (!bpCompareKeys) return;
      bpCompareKeys.clear();
      (keys || []).forEach((k) => bpCompareKeys.add(String(k)));
      if (bpCompare) bpCompare.otherPairKeys = Array.from(bpCompareKeys);
    }

    // Workstation editing hooks (no-ops unless R2DTBpEdit attaches listeners).
    const pairSelectListeners = [];
    const residueSelectListeners = [];

    enhanceSelectedNucleotideStyleOnce();

    const { labelToAuth, authToLabel, labelToChain } = buildLabelMaps(apiData);
    // Only worth naming the chain in the selection tooltip when there's more
    // than one on this panel (a widened multi-chain reference, e.g. a real
    // dimer) -- for an ordinary single-chain structure it'd just be noise.
    const hasMultipleChains = new Set(Object.values(labelToChain)).size > 1;

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
  patchAutoZoomOnSelect();

  // The plugin's pathOrNucleotide() rebuilds the entire SVG inner group on
  // every filter checkbox change, which makes all base-pair lines blink.
  // Toggle path visibility instead once the diagram is painted.
  let bpPathsMaterialized = false;
  let onBpFilterUpdated = () => {};

  // All 12 LW families, including SH/HW/SW (plugin stores some pairs under
  // canonical swapped codes like tHS while the SVG class stays tSH_…).
  const BP_PATH_ID_RE =
    /\b([ct](?:WW|WH|HW|WS|SW|HH|HS|SH|SS)_\d+_\d+)\b/;

  function nestedBpInput() {
    return root.querySelector(`#nestedBP-${PDB_LOWER}`)
      || root.querySelector('#nestedBP');
  }

  function pathIdsInDisplayHtml(html) {
    const ids = new Set();
    if (!html) return ids;
    let m;
    const re = new RegExp(BP_PATH_ID_RE.source, 'g');
    while ((m = re.exec(html))) ids.add(m[1]);
    return ids;
  }

  function pairKeysInDisplayHtml(html) {
    const keys = new Set();
    pathIdsInDisplayHtml(html).forEach((id) => keys.add(bpPairKeyFromPathId(id)));
    return keys;
  }

  function bpPathIdFromElement(el) {
    const m = (el.getAttribute('class') || '').match(BP_PATH_ID_RE);
    return m ? m[1] : null;
  }

  function familyFromPathId(pathId) {
    const m = (pathId || '').match(/^([a-zA-Z]{3})_/);
    return m ? m[1] : '';
  }

  function isBpFamilyFilterOn(family) {
    const all = root.querySelector('#Checkbox_All');
    if (all?.checked) return true;
    const f = String(family || '');
    if (!f) return true;
    const ids = new Set([f, pluginLwFamily(f)]);
    for (const id of ids) {
      const cb = root.querySelector(`#Checkbox_${id}`);
      if (cb?.checked) return true;
    }
    // Families the filter UI doesn't list (ncWW etc.) stay visible when All
    // is off only if no family boxes exist yet.
    return root.querySelectorAll('input[id^="Checkbox_"]').length === 0;
  }

  function annPairKey(a) {
    return bpPairKeyFromPathId(`cWW_${a.seq_id1}_${a.seq_id2}`);
  }

  function isCrossingAnn(a) {
    return !!(a && a.crossing && String(a.crossing) !== '0');
  }

  /** Whether an edited/live path should be shown given current filters. */
  function shouldShowBpPath(pathId, annByKey) {
    if (!pathId) return false;
    const ui = rnaPlugin.uiTemplateService;
    if (!ui) return true;
    const nested = nestedBpInput()?.checked;
    const displayHtml = nested ? ui.displayNestedBaseStrs : ui.displayBaseStrs;
    const visibleIds = pathIdsInDisplayHtml(displayHtml);
    if (visibleIds.has(pathId)) return true;
    // Refamily: same pair, new LW class — keep if the original pair was shown.
    const key = bpPairKeyFromPathId(pathId);
    if (pairKeysInDisplayHtml(displayHtml).has(key)) return true;
    // Newly added pairs aren't in the plugin display strings.
    const family = familyFromPathId(pathId);
    if (!isBpFamilyFilterOn(family)) return false;
    const ann = annByKey.get(key);
    if (nested && isCrossingAnn(ann)) return false;
    return true;
  }

  function isCanonicalWatsonCrick(family, nt1, nt2) {
    if (family !== 'cWW') return false;
    const a = nt1.toUpperCase();
    const b = nt2.toUpperCase();
    return (a === 'A' && (b === 'U' || b === 'T'))
      || (b === 'A' && (a === 'U' || a === 'T'))
      || (a === 'G' && b === 'C')
      || (a === 'C' && b === 'G');
  }

  let fr3dPairLookup = null;
  function getFr3dPairLookup() {
    if (!fr3dPairLookup) {
      fr3dPairLookup = new Map();
      (fr3dData.annotations || []).forEach((a) => {
        fr3dPairLookup.set(`${a.seq_id1}_${a.seq_id2}`, a);
        fr3dPairLookup.set(`${a.seq_id2}_${a.seq_id1}`, a);
      });
    }
    return fr3dPairLookup;
  }

  function buildCwwNonCanonicalD(m, _, b, w, i) {
    const cx = (m + b) / 2;
    const cy = (_ + w) / 2;
    const halfLen = Math.hypot(b - m, w - _) / 2 || 1;
    const frac = Math.min(i / 4 / halfLen, 0.45);
    const lerp = (a, b, t) => a + (b - a) * t;
    const bx = lerp(m, cx, 1 - frac);
    const by = lerp(_, cy, 1 - frac);
    const nx = lerp(cx, b, frac);
    const ny = lerp(cy, w, frac);
    return (
      `M ${m} ${_} ${bx} ${by} ` +
      `M ${cx - i / 4} ${cy} ` +
      `a ${i / 4},${i / 4} 0 1,0 ${i / 2},0 ` +
      `a ${i / 4},${i / 4} 0 1,0 ${-i / 2},0 ` +
      `M ${nx} ${ny} ${b} ${w}`
    );
  }

  const CWW_LINE_D_RE =
    /M\s*([\d.+-]+)\s+([\d.+-]+)\s+([\d.+-]+)\s+([\d.+-]+)/;

  let labelPositionMap = null;
  function getLabelPositionMap() {
    if (labelPositionMap) return labelPositionMap;
    labelPositionMap = new Map();
    const paths = apiData.svg_paths || [];
    const labels = apiData.label_seq_ids || [];
    for (let n = 2; n < paths.length; n++) {
      const label = labels[n - 1];
      if (label == null) continue;
      const tokens = paths[n].split(/[M,]/).filter((s) => s !== '');
      if (tokens.length < 4) continue;
      const x = parseFloat(tokens[2]);
      const y = parseFloat(tokens[3]);
      if (!isNaN(x) && !isNaN(y)) labelPositionMap.set(+label, [x, y]);
    }
    return labelPositionMap;
  }

  function getResidueLocation(seqId) {
    const ui = rnaPlugin.uiTemplateService;
    const id = +seqId;
    if (ui?.locations?.get(id)) return ui.locations.get(id);
    return getLabelPositionMap().get(id) || null;
  }

  function getCwwBpEndpoints(seq1, seq2, strokeWidth) {
    const p1 = getResidueLocation(seq1);
    const p2 = getResidueLocation(seq2);
    if (!p1 || !p2) return null;
    const i = parseFloat(strokeWidth || '2') * 6;
    const x1 = p1[0] + i / 2.5;
    const x2 = p2[0] + i / 2.5;
    const y1 = p1[1] - i / 2.5;
    const y2 = p2[1] - i / 2.5;
    const len = Math.hypot(x1 - x2, y1 - y2) || 1;
    const frac = i / len;
    const lerp = (a, b, t) => a + (b - a) * t;
    return {
      m: lerp(x1, x2, frac),
      _: lerp(y1, y2, frac),
      b: lerp(x1, x2, 1 - frac),
      w: lerp(y1, y2, 1 - frac),
    };
  }

  function parseCwwLineCoords(text) {
    const coords = (text || '').match(CWW_LINE_D_RE);
    if (!coords) return null;
    return {
      m: parseFloat(coords[1]),
      _: parseFloat(coords[2]),
      b: parseFloat(coords[3]),
      w: parseFloat(coords[4]),
    };
  }

  function getCwwLineCoords(pathId) {
    const ui = rnaPlugin.uiTemplateService;
    if (!ui?.baseStrs) return null;
    let found = null;
    ui.baseStrs.forEach((entry) => {
      if (found || !Array.isArray(entry?.[1])) return;
      entry[1].forEach((html) => {
        if (found || typeof html !== 'string' || !html.includes(pathId)) return;
        found = parseCwwLineCoords(html);
      });
    });
    return found;
  }

  function resolveCwwEndpoints(pathId, seq1, seq2, strokeWidth, dText) {
    return (
      parseCwwLineCoords(dText)
      || getCwwLineCoords(pathId)
      || getCwwBpEndpoints(seq1, seq2, strokeWidth)
    );
  }

  function isCwwLineGlyphD(d) {
    return (d.match(/M/g) || []).length >= 3 && /\ba\s/i.test(d);
  }

  function repairCwwPathHtml(html) {
    if (typeof html !== 'string') return html;
    const idMatch = html.match(BP_PATH_ID_RE);
    if (!idMatch || !idMatch[1].startsWith('cWW_')) return html;
    const nums = idMatch[1].match(/cWW_(\d+)_(\d+)/);
    if (!nums) return html;
    const ann = getFr3dPairLookup().get(`${nums[1]}_${nums[2]}`);
    if (!ann || isCanonicalWatsonCrick('cWW', ann.nt1, ann.nt2)) return html;
    const swMatch = html.match(/stroke-width="([^"]+)"/);
    const strokeWidth = swMatch?.[1] || '2';
    const dMatch = html.match(/d="([^"]*)"/);
    const ends = resolveCwwEndpoints(
      idMatch[1], nums[1], nums[2], strokeWidth, dMatch?.[1]
    );
    if (!ends) return html;
    const i = parseFloat(strokeWidth) * 6;
    const glyphD = buildCwwNonCanonicalD(ends.m, ends._, ends.b, ends.w, i);
    return html
      .replace(/d="[^"]*"/, `d="${glyphD}"`)
      .replace(/\s*data-stroke-color="[^"]*"/, '')
      .replace(/stroke="#000"/, 'stroke="#ccc" fill="#ccc"');
  }

  function repairNonCanonicalCwwGlyphs() {
    const lookup = getFr3dPairLookup();
    root.querySelectorAll('path.rnaviewBP').forEach((p) => {
      const id = bpPathIdFromElement(p);
      if (!id || !id.startsWith('cWW_')) return;
      const d = p.getAttribute('d') || '';
      if (p.dataset.r2dtCwwGlyphFixed === '2' || isCwwLineGlyphD(d)) {
        p.dataset.r2dtCwwGlyphFixed = '2';
        return;
      }
      const nums = id.match(/cWW_(\d+)_(\d+)/);
      if (!nums) return;
      const ann = lookup.get(`${nums[1]}_${nums[2]}`);
      if (!ann || isCanonicalWatsonCrick('cWW', ann.nt1, ann.nt2)) return;
      const strokeWidth = p.getAttribute('stroke-width') || '2';
      const ends = resolveCwwEndpoints(
        id, nums[1], nums[2], strokeWidth, d
      );
      if (!ends) return;
      const i = parseFloat(strokeWidth) * 6;
      p.setAttribute('d', buildCwwNonCanonicalD(ends.m, ends._, ends.b, ends.w, i));
      p.setAttribute('stroke', '#ccc');
      p.setAttribute('fill', '#ccc');
      p.removeAttribute('data-stroke-color');
      p.dataset.r2dtCwwGlyphFixed = '2';
    });
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
        if (m && !existing.has(m[1])) {
          toAdd.push(repairCwwPathHtml(html));
        }
      });
    });
    if (toAdd.length) inner.insertAdjacentHTML('beforeend', toAdd.join(''));
    repairNonCanonicalCwwGlyphs();
    bpPathsMaterialized = true;
  }

  function applyBpVisibility() {
    const ui = rnaPlugin.uiTemplateService;
    if (!ui) return;
    materializeAllBpPaths();
    const annByKey = new Map();
    (fr3dData.annotations || []).forEach((a) => {
      annByKey.set(annPairKey(a), a);
    });
    const inner = root.querySelector(`.rnaTopoSvg_${PDB_LOWER}`);
    if (!inner) return;
    inner.querySelectorAll('path.rnaviewBP').forEach((p) => {
      const id = bpPathIdFromElement(p);
      if (!id) return;
      p.style.display = shouldShowBpPath(id, annByKey) ? '' : 'none';
    });
    repairNonCanonicalCwwGlyphs();
    onBpFilterUpdated();
  }

  // Pair keys ("min_max") for the base pairs currently shown in the 2D
  // diagram. Reads the rendered paths directly so it reflects whatever the
  // "Nested only" toggle and family checkboxes produced, without
  // re-deriving the plugin's filter logic.
  function getVisible2dPairKeys() {
    const keys = new Set();
    root.querySelectorAll('path.rnaviewBP').forEach((p) => {
      if (p.style.display === 'none') return;
      const id = bpPathIdFromElement(p);
      if (id) keys.add(bpPairKeyFromPathId(id));
    });
    return keys;
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
  // Declared up front so the post-render fixup pass (which skips
  // currently-selected BPs) can reference them before the click handler
  // block initializes -- temporal-dead-zone errors otherwise.
  const bpHighlightedPaths = new Set();

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
      patchAutoZoomOnSelect();
      svc.zoomReset = function () {
        // Plugin calls zoomReset on deselect; keep the user's zoom if they panned.
        if (userAdjustedView) return;
        userAdjustedView = false;
        applyFitTransform();
      };
    }
    return !!resetBtn;
  }

  const unobserved = apiData.unobserved_label_seq_ids || [];
  const CWW_STROKE = '#888888';
  const CWW_CROSSING_STROKE = '#aaaaaa';
  // Compare-mode diff highlight: recolour base pairs that disagree with the
  // other structure so they stand out at a glance. In the reference panel that
  // means FN (blue) — pairs the model missed; in the model panel FP (red) —
  // pairs the model added that aren't in the reference. TP pairs keep their
  // neutral colour. Same palette as the LBN rows / base-pair-list badges.
  const DIFF_FN_COLOR = '#4363d8';
  const DIFF_FP_COLOR = '#e6194b';
  // Crossing pairs are drawn thinner than nested cWW; width scales with the
  // same layout metric the plugin uses (calculateFontSize / 6), not a fixed px.
  const CROSSING_BP_WIDTH_RATIO = 0.5;
  const CROSSING_BP_STROKE_WIDTH_MIN = 0.25;
  const CROSSING_BP_STROKE_WIDTH_MAX = 1.15;
  let crossingBpStrokeWidth = null;

  function calculateLayoutScale(paths) {
    const cx = [];
    const cy = [];
    const dist = [];
    const last = paths.length - 1;
    paths.forEach((path, o) => {
      if (o === 0 || o === last) return;
      const s = path.split('M').join(',').split(',');
      cx[o] = (Number(s[1]) + Number(s[3])) / 2;
      cy[o] = (Number(s[2]) + Number(s[4])) / 2;
      if (o > 1) {
        dist[o] = Math.hypot(cx[o] - cx[o - 1], cy[o] - cy[o - 1]);
      }
    });
    const sorted = dist.filter((d) => d > 0).sort((a, b) => a - b);
    if (sorted.length === 0) return 12;
    return 0.9 * sorted[Math.floor(0.05 * sorted.length)];
  }

  function getBaseBpStrokeWidth() {
    return calculateLayoutScale(apiData.svg_paths || []) / 6;
  }

  function getCrossingBpStrokeWidth() {
    if (crossingBpStrokeWidth === null) {
      const scaled = getBaseBpStrokeWidth() * CROSSING_BP_WIDTH_RATIO;
      crossingBpStrokeWidth = String(Math.max(
        CROSSING_BP_STROKE_WIDTH_MIN,
        Math.min(CROSSING_BP_STROKE_WIDTH_MAX, scaled)
      ));
    }
    return crossingBpStrokeWidth;
  }

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

  // The plugin draws each cWW base-pair line between the two nucleotide
  // centres, inset from each end by the layout scale `i` (= one nucleotide
  // radius, so the line clears the letters). Line length is therefore
  // `y - 2i`, where `y` is the distance between the paired bases. When a
  // helix is laid out narrow (paired strands only ~2i apart, common in
  // tightly-wound R2DT templates like 9bz1), the two insets nearly meet and
  // the line collapses to a dot at the midpoint — the base pair looks
  // unmarked. Cap the inset fraction so a legible portion of the line always
  // survives; a no-op for helices wide enough that the plugin's own line is
  // already fine.
  const MAX_BP_INSET_FRAC = 0.32; // keep >=36% of each base-pair line visible

  function fixCollapsedCwwLine(el, aId, bId) {
    const d = el.getAttribute('d') || '';
    if (isCwwLineGlyphD(d)) return;         // non-canonical special glyph
    const nums = d.match(/-?\d[\d.eE+-]*/g);
    if (!nums || nums.length !== 4) return; // only simple 2-point lines
    const p1 = getResidueLocation(aId);
    const p2 = getResidueLocation(bId);
    if (!p1 || !p2) return;
    const i = getBaseBpStrokeWidth() * 6;
    // Same endpoint offsets the plugin uses (see pathOrNucleotide).
    const x1 = p1[0] + i / 2.5, y1 = p1[1] - i / 2.5;
    const x2 = p2[0] + i / 2.5, y2 = p2[1] - i / 2.5;
    const len = Math.hypot(x1 - x2, y1 - y2) || 1;
    if (i / len <= MAX_BP_INSET_FRAC) return; // plugin's line is already fine
    const frac = MAX_BP_INSET_FRAC;
    const lerp = (a, b, t) => a + (b - a) * t;
    const mx = lerp(x1, x2, frac), my = lerp(y1, y2, frac);
    const bx = lerp(x1, x2, 1 - frac), by = lerp(y1, y2, 1 - frac);
    el.setAttribute('d', `M${mx} ${my} ${bx} ${by}`);
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
      if (parsed.bp === 'cWW' && !bpHighlightedPaths.has(el)) {
        const stroke = isCrossing ? CWW_CROSSING_STROKE : CWW_STROKE;
        el.setAttribute('stroke', stroke);
        // Remember the greyed-out value so the click handler's
        // restore-on-next-click brings it back to grey, not the
        // plugin's default black.
        el.dataset.r2dtOrigStroke = stroke;
        // Rescue base-pair lines the narrow-helix layout collapsed to a dot.
        fixCollapsedCwwLine(el, parsed.a, parsed.b);
      }
      if (isCrossing) {
        el.setAttribute('stroke-width', getCrossingBpStrokeWidth());
      }
      // Compare-mode: recolour disagreements (FN in the reference panel, FP in
      // the model panel). Runs after the cWW/glyph default colouring above, for
      // every family. Preserve open (fill:none) vs filled glyphs — only recolour
      // a fill that was actually painted. Skip selected pairs so the orange
      // selection wins; store the diff colour as the restore value.
      if (!bpHighlightedPaths.has(el)) {
        const kind = classifyBpPair(+parsed.a, +parsed.b);
        if (kind === 'FN' || kind === 'FP') {
          const color = kind === 'FN' ? DIFF_FN_COLOR : DIFF_FP_COLOR;
          el.setAttribute('stroke', color);
          const fill = el.getAttribute('fill');
          if (fill && fill !== 'none') el.setAttribute('fill', color);
          el.dataset.r2dtOrigStroke = color;
        }
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
        const spec = BP_FAMILY_GLYPH[family] || BP_FAMILY_GLYPH[flipLw(family)];
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
    const dropdown = root.querySelector('#mainMenu .menu-dropdown');
    if (dropdown?.classList.contains('show')) return true;
    const checkboxes = root.querySelector('#checkboxes');
    if (!checkboxes) return false;
    return getComputedStyle(checkboxes).display !== 'none';
  }

  function closeBasePairsPanel() {
    if (!isBasePairsPanelOpen()) return;
    root.querySelector('#bpFilterBtn')?.click();
  }

  function mountFilterPanel() {
    const checkboxes = root.querySelector('#checkboxes');
    if (!checkboxes || !panel2d || checkboxes.dataset.r2dtFilterMoved) {
      return !!checkboxes?.dataset.r2dtFilterMoved;
    }
    panel2d.appendChild(checkboxes);
    checkboxes.classList.add('r2dt-bp-filter-panel');
    checkboxes.dataset.r2dtFilterMoved = '1';
    return true;
  }

  function installBpFilterPatch() {
    const ui = rnaPlugin.uiTemplateService;
    const filterBtn = root.querySelector('#bpFilterBtn');
    const dropdown = root.querySelector('#mainMenu .menu-dropdown');
    const checkboxes = root.querySelector('#checkboxes');
    if (!ui || !filterBtn || !dropdown || !checkboxes || filterBtn.dataset.r2dtFilterBound) {
      return;
    }
    const fresh = filterBtn.cloneNode(true);
    fresh.setAttribute('aria-haspopup', 'true');
    fresh.setAttribute('aria-expanded', 'false');
    filterBtn.replaceWith(fresh);
    fresh.dataset.r2dtFilterBound = '1';
    fresh.addEventListener('click', () => {
      ui.checkboxesExpanded = !ui.checkboxesExpanded;
      dropdown.classList.toggle('show', ui.checkboxesExpanded);
      root.querySelector('#bpFilterBtnIcon')
        ?.classList.toggle('active', ui.checkboxesExpanded);
      checkboxes.style.display = ui.checkboxesExpanded ? 'block' : 'none';
      fresh.setAttribute('aria-expanded', String(ui.checkboxesExpanded));
      if (ui.checkboxesExpanded) {
        setTimeout(closeBpListPanel, 0);
      }
    });
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
      const pairText = m[1].trim();
      const familyText = pluginLwFamily(m[2]);
      const ntMatch = pairText.match(/([A-Z])(\d+)\s*-\s*([A-Z])(\d+)/i);
      const canonical = ntMatch
        ? isCanonicalWatsonCrick(familyText, ntMatch[1], ntMatch[3])
        : familyText === 'cWW';

      li.textContent = '';
      const pair = document.createElement('span');
      pair.className = 'r2dt-bp-list-pair';
      pair.textContent = pairText;
      const family = document.createElement('span');
      family.className = 'r2dt-bp-list-family';
      family.textContent = familyText;

      // TP/FP/FN badge (compare mode): does this pair exist in the other structure?
      let cmp = null;
      if (bpCompareKeys && ntMatch) {
        const kind = classifyBpPair(+ntMatch[2], +ntMatch[4]);
        if (kind) {
          cmp = document.createElement('span');
          cmp.className = `r2dt-bp-list-cmp r2dt-bp-list-cmp--${kind.toLowerCase()}`;
          cmp.textContent = kind;
          cmp.title = {
            TP: 'True positive — this pair is in both structures',
            FP: 'False positive — this pair is only in the model',
            FN: 'False negative — this pair is in the reference but missing from the model',
          }[kind];
          li.classList.add(`r2dt-bp-list-item--${kind.toLowerCase()}`);
        }
      }

      // non-WC marker as a fixed-width slot (empty but still reserving width on
      // canonical rows) so the family column stays aligned down the list.
      const tag = document.createElement('span');
      tag.className = 'r2dt-bp-list-tag';
      if (!canonical) {
        li.classList.add('r2dt-bp-list-item--nonwc');
        tag.textContent = 'non-WC';
        tag.setAttribute('aria-label', 'Non-Watson–Crick base pair');
      } else {
        tag.classList.add('r2dt-bp-list-tag--empty');
      }
      // Order: nt1-nt2 pair | non-WC (if any) | base-pair family | TP/FP/FN.
      // pair flex-grows so non-WC, family and the TP/FP/FN badge pin to the
      // right edge and line up in their own aligned columns.
      li.append(pair, tag, family, ...(cmp ? [cmp] : []));
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

  function findBpPathByPathId(pathID) {
    const m = (pathID || '').match(BP_PATH_ID_RE);
    if (!m) return null;
    return root.querySelector(`.rnaviewBP[class*="${m[1]}"]`);
  }

  function installBpListDialogPatch() {
    const ui = rnaPlugin.uiTemplateService;
    if (!ui?.renderBpListDialog || ui._r2dtBpListPatched) return;
    ui._r2dtBpListPatched = true;
    ui.renderBpListDialog = function r2dtRenderBpListDialog(toggle) {
      const dialog = root.querySelector(`#bpListDialog-${PDB_LOWER}`);
      const nestedInput = nestedBpInput();
      if (!dialog || !nestedInput) return;
      if (toggle) {
        const open = getComputedStyle(dialog).display !== 'none';
        dialog.style.display = open ? 'none' : 'flex';
      }
      dialog.innerHTML = '';
      const ul = document.createElement('ul');
      // Prefer live annotations (includes workstation edits) over the plugin's
      // frozen bpLabels, which still point at pre-edit path IDs / families.
      const annByKey = new Map();
      (fr3dData.annotations || []).forEach((a) => {
        const key = annPairKey(a);
        if (!annByKey.has(key)) annByKey.set(key, a);
      });
      const visibleKeys = new Set();
      root.querySelectorAll('path.rnaviewBP').forEach((p) => {
        if (p.style.display === 'none') return;
        const id = bpPathIdFromElement(p);
        if (id) visibleKeys.add(bpPairKeyFromPathId(id));
      });
      // If visibility hasn't been applied yet, fall back to display-html keys
      // plus every annotation (first paint).
      const useVisibleFilter = visibleKeys.size > 0
        || root.querySelector('path.rnaviewBP');
      const rows = [...annByKey.values()]
        .filter((a) => {
          const key = annPairKey(a);
          if (!useVisibleFilter) return true;
          if (visibleKeys.has(key)) return true;
          // Path may not exist yet for a just-added pair before geometry sync;
          // still list it when its family filter allows.
          return isBpFamilyFilterOn(a.bp || 'cWW')
            && !(nestedInput.checked && isCrossingAnn(a));
        })
        .sort((x, y) => {
          const ax = +x.seq_id1;
          const ay = +y.seq_id1;
          const bx = +x.seq_id2;
          const by = +y.seq_id2;
          return ax - ay || bx - by;
        });
      rows.forEach((a) => {
        const li = document.createElement('li');
        const nt1 = a.nt1 || a.unit1 || 'N';
        const nt2 = a.nt2 || a.unit2 || 'N';
        const fam = a.bp || 'cWW';
        li.textContent = `${nt1}${a.seq_id1} - ${nt2}${a.seq_id2} ; ${fam}`;
        li.style.cursor = 'pointer';
        const i = +a.seq_id1;
        const j = +a.seq_id2;
        li.addEventListener('mouseenter', () => {
          findBPPath(i, j)?.dispatchEvent(
            new Event('mouseover', { bubbles: true })
          );
          const tip = document.getElementById('tooltip');
          if (tip) tip.style.display = 'none';
        });
        li.addEventListener('mouseleave', () => {
          findBPPath(i, j)?.dispatchEvent(
            new Event('mouseout', { bubbles: true })
          );
        });
        ul.appendChild(li);
      });
      dialog.appendChild(ul);
      ensureBpListPanelTitle();
      applyBpListItemLabels();
      normalizeBpListScroll();
      const listBtn = root.querySelector(`#rnaTopologyBPList-${PDB_LOWER}`);
      if (listBtn) {
        listBtn.setAttribute('aria-expanded', String(isBpListPanelOpen()));
      }
    };
  }

  function setupBpListToolbar() {
    installBpListDialogPatch();
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
          if (!ev.target.closest(
            '#mainMenu .menu-dropdown, #checkboxes, .r2dt-bp-filter-panel'
          )) {
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
      mountFilterPanel();
      installBpFilterPatch();
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
    const nestedWrap = nestedBpInput()?.parentElement;
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
      if (nestedText) {
        // Full + short variants; CSS shows one based on toolbar width.
        nestedText.innerHTML =
          '<span class="r2dt-lbl-full">Nested only</span>'
          + '<span class="r2dt-lbl-short">Nested</span>';
      }
      let nestedInput = nestedWrap.querySelector(
        `#nestedBP-${PDB_LOWER}, #nestedBP`
      );
      const nestedLabel = nestedWrap.querySelector('label');
      const nestedId = `nestedBP-${PDB_LOWER}`;
      if (nestedInput) {
        nestedInput.setAttribute('aria-label', 'Only nested base pairs');
      }
      if (nestedLabel) {
        nestedLabel.classList.add('r2dt-nested-toggle');
      }
      // Plugin binds #nestedBP to the original pathOrNucleotide (full SVG
      // rebuild). Rebind to our patched version so non-canonical cWW glyphs
      // stay repaired after toggling nested view.
      if (nestedInput && !nestedInput.dataset.r2dtNestedBound) {
        const nestedFresh = nestedInput.cloneNode(true);
        nestedFresh.checked = nestedInput.checked;
        nestedFresh.id = nestedId;
        nestedInput.replaceWith(nestedFresh);
        nestedFresh.dataset.r2dtNestedBound = '1';
        nestedInput = nestedFresh;
        nestedInput.addEventListener('change', () => {
          if (uiSvc?.pathOrNucleotide) uiSvc.pathOrNucleotide();
        });
      }
      if (nestedLabel && nestedInput) {
        nestedLabel.setAttribute('for', nestedId);
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
          const label = document.createElement('span');
          label.className = 'r2dt-bp-filter-label';
          label.innerHTML =
            '<span class="r2dt-lbl-full">Base Pairs</span>'
            + '<span class="r2dt-lbl-short">BPs</span>';
          node.replaceWith(label);
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

    mountFilterPanel();
    installBpFilterPatch();
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

  // Hide family checkboxes with no pairs in the current annotations, show
  // newly present ones, and refresh the BPs n/m badge. Safe to call after
  // workstation edits (add / delete / refamily).
  function syncBpFilterToAnnotations(opts) {
    const ensureAllOn = !!(opts && opts.ensureAllOn);
    const presentBPs = presentLwFamilies(fr3dData.annotations);
    const boxes = root.querySelectorAll('input[id^="Checkbox_"]');
    if (boxes.length === 0) return false;
    boxes.forEach((cb) => {
      if (cb.id === 'Checkbox_All') return;
      const family = cb.id.slice('Checkbox_'.length);
      const td = cb.closest('td');
      const present = presentBPs.has(family);
      const wasHidden = !td || td.style.display === 'none'
        || !isFilterCheckboxVisible(cb);
      if (!present) {
        if (td) td.style.display = 'none';
        return;
      }
      if (td) td.style.display = '';
      // A newly appeared family (e.g. user added a tHS) should be checked
      // so the new pairs are visible immediately.
      if (wasHidden && !cb.checked) cb.checked = true;
    });
    hideEmptyFilterRows();
    const all = root.querySelector('#Checkbox_All');
    if (all) {
      if (ensureAllOn && !all.checked) {
        all.click();
      } else {
        // Keep "All" in sync with whether every visible family is checked.
        const famBoxes = [...boxes].filter(
          (cb) => cb.id !== 'Checkbox_All' && isFilterCheckboxVisible(cb)
        );
        const allOn = famBoxes.length > 0 && famBoxes.every((cb) => cb.checked);
        if (all.checked !== allOn) all.checked = allOn;
      }
    }
    ensureFilterPanelTitle();
    injectFilterGlyphs();
    mountFilterLegend();
    updateFilterBadge();
    return true;
  }

  // Hide base-pair family checkboxes that have no annotations in this
  // structure's fr3d.json, so the dropdown only lists families the user
  // can actually toggle. Then enable the remaining families ("All").
  (function tidyBPFilter() {
    let attempts = 0;
    const tick = () => {
      if (!syncBpFilterToAnnotations({ ensureAllOn: true })) {
        if (attempts++ > 40) return;
        setTimeout(tick, 100);
        return;
      }
      applyBpVisibility();
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
    labelToChain,
    root,
  };
  if (link3d) window.__r2dt = ctx.handles;

  const selectInMolstar = createMolstarSelector(
    molstar, labelToAuth, CHAIN_ID, labelToChain
  );

  // Set by _renderLBN when the panel is shown; keeps LBN in sync with every
  // base-pair selection path (list, 2D line click, API).
  let lbnHighlightFn = null;
  // Set by _renderLBN; greys out LBN pairs that are currently hidden in the
  // 2D diagram (via the "Nested only" toggle or the family checkboxes).
  let lbnVisibilityFn = null;

  // Compare mode: match the 2D click highlight to this panel's own 3D colour
  // (reference green / model blue) instead of a colour-neutral orange, so the
  // 2D and 3D highlights read as the same selection. Falls back to the
  // original orange for the standalone viewer, which has no reference/model
  // structure colour to match.
  const baseColor = ctx.baseColor || null;
  const BP_SELECTED_COLOR = baseColor ? rgbToCss(highlightColorFor(baseColor)) : 'orange';
  const RESIDUE_PARTNER_COLOR = baseColor
    ? rgbToCss(partnerHighlightColorFor(baseColor)) : '#ffab40';
  const RNA_ENTITY_ID = rnaPlugin.options?.entityId || '1';
  let lastResidueSelection = null;

  function getPartnersForResidue(label) {
    const partners = new Set();
    const id = +label;
    (fr3dData.annotations || []).forEach((a) => {
      const s1 = +a.seq_id1;
      const s2 = +a.seq_id2;
      if (s1 === id) partners.add(s2);
      else if (s2 === id) partners.add(s1);
    });
    return Array.from(partners);
  }

  function clearBasePairHighlights() {
    bpHighlightedPaths.forEach((el) => {
      const prevOrig = el.dataset.r2dtOrigStroke;
      if (prevOrig !== undefined) el.setAttribute('stroke', prevOrig);
    });
    bpHighlightedPaths.clear();
  }

  function highlightBPPath(pathEl) {
    if (!pathEl) return;
    if (pathEl.dataset.r2dtOrigStroke === undefined) {
      pathEl.dataset.r2dtOrigStroke = pathEl.getAttribute('stroke') || '';
    }
    pathEl.setAttribute('stroke', BP_SELECTED_COLOR);
    bpHighlightedPaths.add(pathEl);
  }

  function setPartnerLabels(partners) {
    const svc = window.UiActionsService;
    if (!svc) return;
    svc.__r2dtPartnerLabels = new Set(partners);
  }

  // selectResidueRange's 3rd argument is a sequence letter for the tooltip,
  // not a colour — pass colours via UiActionsService.selectNucleotide instead.
  function updateSelectionTooltip(labels) {
    const svc = window.UiActionsService;
    const tip = document.getElementById(`${PDB_LOWER}-rnaTopologyTooltip`);
    if (!tip || !svc?.buildResidueLabel) return;
    if (!labels || labels.length === 0) {
      tip.style.display = 'none';
      return;
    }
    tip.style.display = 'inline';
    let text = svc.buildResidueLabel(labels, undefined);
    if (hasMultipleChains) {
      // Grouped by chain (not a flat "(chains 1, 0)" list) so which residue
      // is on which chain is unambiguous -- a flat list's order isn't
      // guaranteed to match buildResidueLabel's own (possibly re-sorted)
      // residue-number order above, which could misattribute a residue to
      // the wrong chain at a glance.
      const byChain = new Map();
      labels.forEach((l) => {
        const c = labelToChain[l];
        if (c == null) return;
        if (!byChain.has(c)) byChain.set(c, []);
        byChain.get(c).push(l);
      });
      if (byChain.size === 1) {
        text += ` (chain ${byChain.keys().next().value})`;
      } else if (byChain.size > 1) {
        const parts = Array.from(byChain, ([c, ls]) => `chain ${c}: ${ls.join(', ')}`);
        text += ` (${parts.join('; ')})`;
      }
    }
    tip.innerHTML =
      `<strong>Selected ${text}</strong>`
      + ' <button type="button" class="r2dt-deselect-btn"'
      + ' aria-label="Clear selection" title="Clear selection (Esc)">&#10005;</button>';
    const clearBtn = tip.querySelector('.r2dt-deselect-btn');
    if (clearBtn) {
      clearBtn.addEventListener('click', (ev) => {
        ev.preventDefault();
        ev.stopPropagation();
        clearResidueSelection();
      });
    }
  }

  function selectNucleotidesIn2d(primary, partnerList) {
    const svc = window.UiActionsService;
    if (!svc) return;
    setPartnerLabels(partnerList);
    rnaPlugin.clearSelection(undefined, true);
    rnaPlugin.clearHighlight(true);
    svc.selectNucleotide(
      PDB_LOWER, RNA_ENTITY_ID, [primary], 'click', false,
      undefined, undefined, BP_SELECTED_COLOR, false, true
    );
    partnerList.forEach((p) => {
      svc.selectNucleotide(
        PDB_LOWER, RNA_ENTITY_ID, [p], 'click', false,
        undefined, undefined, RESIDUE_PARTNER_COLOR, true, true
      );
    });
    updateSelectionTooltip([primary, ...partnerList]);
  }

  function emitResidueSelectEvent(detail) {
    document.dispatchEvent(new CustomEvent('r2dt-residue-select', { detail }));
  }

  function nucleotideLabelFromClickTarget(target) {
    const el = target.nodeName === 'text'
      ? target
      : target.closest?.('text.rnaviewEle');
    if (!el) return null;
    const m = (el.getAttribute('class') || '').match(
      new RegExp(`rnaview_${PDB_LOWER}_(\\d+)\\b`)
    );
    return m ? +m[1] : null;
  }

  function bindNucleotideCaptureClicks(container) {
    if (!container || container.dataset.r2dtNtCaptureBound) return;
    container.dataset.r2dtNtCaptureBound = '1';
    // Capture before the plugin's inline onclick toggles a single nt in 2D only.
    container.addEventListener('click', (ev) => {
      const label = nucleotideLabelFromClickTarget(ev.target);
      if (label == null) return;
      ev.preventDefault();
      ev.stopImmediatePropagation();
      selectResidue(label);
    }, true);
  }

  function hasActiveSelection() {
    if (lastResidueSelection || bpHighlightedPaths.size > 0) return true;
    const svc = window.UiActionsService;
    return !!(svc?.selected && svc.selected.size > 0);
  }

  function isDiagramChromeClick(target) {
    return !!target.closest?.(
      'button, input, select, label, a, .help-tooltip, .pdb-rna-view-btn, ' +
      '[id^="bpListDialog"], [id^="rnaTopologyZoom"], [id^="rnaTopologyReset"], ' +
      '[id$="-rnaTopologyTooltip"]'
    );
  }

  function bindDiagramBackgroundClick(container) {
    if (!container || container.dataset.r2dtBgClickBound) return;
    container.dataset.r2dtBgClickBound = '1';
    // Capture phase: the plugin's container handler calls stopPropagation(),
    // so a bubble listener on panel2d would never run.
    container.addEventListener('click', (ev) => {
      if (nucleotideLabelFromClickTarget(ev.target) != null) return;
      if (ev.target.closest?.('path[class*="rnaviewBP"]')) return;
      if (isDiagramChromeClick(ev.target)) return;
      if (!hasActiveSelection()) return;
      clearResidueSelection();
    }, true);
  }

  function clearResidueSelection() {
    lastResidueSelection = null;
    clearBasePairHighlights();
    setPartnerLabels([]);
    rnaPlugin.clearSelection(undefined, true);
    rnaPlugin.clearHighlight(true);
    updateSelectionTooltip([]);
    if (link3d && molstar) selectInMolstar([]);
    if (lbnHighlightFn) lbnHighlightFn(null, []);
    emitResidueSelectEvent({ pdbId: PDB_LOWER, cleared: true });
  }

  function selectResidue(label) {
    const id = +label;
    if (lastResidueSelection && lastResidueSelection.label === id) {
      clearResidueSelection();
      return;
    }
    const partners = getPartnersForResidue(id);
    lastResidueSelection = { label: id, partners };
    clearBasePairHighlights();
    partners.forEach((p) => highlightBPPath(findBPPath(id, p)));
    selectNucleotidesIn2d(id, partners);
    if (link3d && molstar) selectInMolstar([id, ...partners]);
    if (lbnHighlightFn) lbnHighlightFn(id, partners.map((p) => [id, p]));
    emitResidueSelectEvent({ pdbId: PDB_LOWER, label: id, partners, cleared: false });
    residueSelectListeners.forEach((fn) => {
      try { fn(id); } catch (_) { /* ignore listener errors */ }
    });
  }

  bindNucleotideCaptureClicks(panel2d);
  bindDiagramBackgroundClick(panel2d);

  // Pressing Escape clears the current selection (in addition to clicking
  // the same residue again or clicking empty diagram space).
  function onDeselectKey(ev) {
    if (ev.key !== 'Escape' && ev.key !== 'Esc') return;
    if (!hasActiveSelection()) return;
    clearResidueSelection();
  }
  document.addEventListener('keydown', onDeselectKey);
  ctx.cleanup.push(() => document.removeEventListener('keydown', onDeselectKey));

  // Select a base pair: colour its 2D line orange (if its path is in the
  // DOM) and select both partner residues in 3D. `pathEl` may be null --
  // e.g. when triggered from the base-pair list while that pair's line
  // isn't currently rendered.
  function selectBasePair(a, b, pathEl) {
    lastResidueSelection = null;
    clearBasePairHighlights();
    highlightBPPath(pathEl);
    const svc = window.UiActionsService;
    setPartnerLabels([]);
    rnaPlugin.clearSelection(undefined, true);
    rnaPlugin.clearHighlight(true);
    if (svc) {
      svc.selectNucleotide(
        PDB_LOWER, RNA_ENTITY_ID, [a], 'click', false,
        undefined, undefined, BP_SELECTED_COLOR, false, true
      );
      svc.selectNucleotide(
        PDB_LOWER, RNA_ENTITY_ID, [b], 'click', false,
        undefined, undefined, BP_SELECTED_COLOR, true, true
      );
    }
    updateSelectionTooltip([a, b]);
    if (link3d && molstar) selectInMolstar([a, b]);
    if (lbnHighlightFn) lbnHighlightFn(a, [[a, b]]);
    // Drive the shared 3D in compare mode (panels run link3d=false): emit the
    // same event as a residue click so onAnyResidueSelect selects+zooms both
    // partners. `pair: true` tells the mirror to highlight this exact pair
    // (both residues at the same positions) in the other panel, rather than
    // re-deriving a partner for `a` from that panel's own base pairs.
    emitResidueSelectEvent({
      pdbId: PDB_LOWER, label: a, partners: [b], cleared: false, pair: true,
    });
    const famMatch = (pathEl?.getAttribute('class') || '')
      .match(/([a-zA-Z]{3})_\d+_\d+/);
    const family = famMatch ? famMatch[1] : (
      (fr3dData.annotations || []).find((ann) => {
        const x = +ann.seq_id1;
        const y = +ann.seq_id2;
        return (x === a && y === b) || (x === b && y === a);
      }) || {}
    ).bp || 'cWW';
    pairSelectListeners.forEach((fn) => {
      try { fn(a, b, family); } catch (_) { /* ignore listener errors */ }
    });
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
    if (lbnVisibilityFn) lbnVisibilityFn(getVisible2dPairKeys());
    // Compare mode: this panel has no LBN widget of its own (createCompare
    // renders one shared widget for both structures), so push visibility
    // changes ("Nested only" etc.) to it via this optional hook instead.
    ctx.onBpVisibilityChange?.(getVisible2dPairKeys());
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
      selectResidue(label);
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
    // Prefer inline LBN data (compare panels embed it); else fetch lbn.json.
    let lbnData = ctx.lbnData || null;
    if (!lbnData) {
      try {
        const resp = await fetch(resolveUrl(baseUrl, 'lbn.json'));
        if (resp.ok) lbnData = await resp.json();
      } catch (_) { /* ignore */ }
    }

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
    html += '<div class="lbn-row"><span class="lbn-label">seq:</span> ';
    for (let i = 0; i < N; i++) {
      html += `<span data-pos="${i + 1}" class="lbn-nt">${SEQ[i]}</span>`;
    }
    html += '</div>';

    for (const row of data.rows) {
      html += `<div class="lbn-row"><span class="lbn-label">${row.label}:</span> `;
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

    // Highlight by base pair (edge), not by column. `primary` is the
    // clicked residue (or null); `pairs` is the list of selected [a, b]
    // edges. The clicked residue's nucleotide and its bracket(s) get the
    // "selected" colour; partner nucleotides and the partner side of each
    // selected pair get the "partner" colour. Brackets are highlighted
    // only when their own {pos, partner} edge is selected, so a partner's
    // unrelated pairs in other layers are never partially lit.
    function _lbnHighlight(primary, pairs) {
      highlighted.forEach((sp) => sp.classList.remove('lbn-selected', 'lbn-partner'));
      highlighted = [];
      const prim = primary == null ? null : +primary;
      const edges = pairs || [];
      if (prim == null && edges.length === 0) return;

      const pairKeys = new Set(
        edges.map(([a, b]) => `${Math.min(a, b)}_${Math.max(a, b)}`)
      );
      const partnerPositions = new Set();
      edges.forEach(([a, b]) => {
        if (a !== prim) partnerPositions.add(a);
        if (b !== prim) partnerPositions.add(b);
      });

      // Sequence-row nucleotides: show every residue involved.
      container.querySelectorAll('.lbn-nt').forEach((sp) => {
        const p = +sp.dataset.pos;
        let cls = null;
        if (p === prim) cls = 'lbn-selected';
        else if (partnerPositions.has(p)) cls = 'lbn-partner';
        if (cls) {
          sp.classList.add(cls);
          highlighted.push(sp);
        }
      });

      // Bracket spans: only light the two brackets of a selected pair.
      container.querySelectorAll('.lbn-bp').forEach((sp) => {
        const p = +sp.dataset.pos;
        const partner = sp.dataset.partner ? +sp.dataset.partner : null;
        if (partner == null) return;
        const key = `${Math.min(p, partner)}_${Math.max(p, partner)}`;
        if (!pairKeys.has(key)) return;
        sp.classList.add(p === prim ? 'lbn-selected' : 'lbn-partner');
        highlighted.push(sp);
      });

      const focusSpan = (prim != null && posSpans[prim] && posSpans[prim][0])
        || highlighted[0];
      if (focusSpan) {
        focusSpan.scrollIntoView({
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
      } else {
        selectResidue(pos);
      }
    });

    // Grey out the pairs that are hidden in the 2D diagram so the list
    // mirrors what the "Nested only" toggle / family checkboxes show.
    // `visibleKeys` is a Set of "min_max" pair keys (null => show all).
    function _lbnVisibility(visibleKeys) {
      container.querySelectorAll('.lbn-bp').forEach((sp) => {
        const pos     = +sp.dataset.pos;
        const partner = sp.dataset.partner ? +sp.dataset.partner : null;
        if (partner == null) return;
        const key    = `${Math.min(pos, partner)}_${Math.max(pos, partner)}`;
        const hidden = visibleKeys ? !visibleKeys.has(key) : false;
        sp.classList.toggle('lbn-bp-hidden', hidden);
      });
      // Dim a layer's label when none of its pairs are visible in 2D.
      container.querySelectorAll('.lbn-row').forEach((rowEl) => {
        const bps = rowEl.querySelectorAll('.lbn-bp');
        if (bps.length === 0) return;
        const anyVisible = Array.from(bps).some(
          (sp) => !sp.classList.contains('lbn-bp-hidden')
        );
        rowEl.classList.toggle('lbn-row--hidden', !anyVisible);
      });
    }

    // 2D/3D → LBN sync is handled by selectResidue / selectBasePair via
    // lbnHighlightFn; no separate listeners needed here.

    lbnHighlightFn = _lbnHighlight;
    lbnVisibilityFn = _lbnVisibility;
    // Reflect the current 2D filter state immediately on first render.
    _lbnVisibility(getVisible2dPairKeys());

    // Expose for console debugging (only the single-viewer sets window.__r2dt;
    // compare panels run with link3d=false and leave it undefined).
    if (window.__r2dt) {
      window.__r2dt.lbnData      = data;
      window.__r2dt.lbnHighlight = _lbnHighlight;
    }
  }
  // ── end LBN ─────────────────────────────────────────────────────────────

  ctx.handles.selectResidue = selectResidue;
  // Non-toggling select for cross-panel mirroring (avoids clearing when the
  // mirrored panel already has that residue selected).
  ctx.handles.selectResidueForce = (label) => {
    if (lastResidueSelection && lastResidueSelection.label === +label) return;
    selectResidue(label);
  };
  ctx.handles.clearSelection = clearResidueSelection;
  ctx.handles.selectBasePair = (a, b) => selectBasePair(a, b, findBPPath(a, b));
  // Compare mode: lets createCompare seed the shared LBN widget's initial
  // per-panel visibility state ("Nested only" etc.) without waiting for a
  // filter change.
  ctx.handles.getVisiblePairKeys = getVisible2dPairKeys;

  // Thin adapter for workstation base-pair editing (R2DTBpEdit). Keeps edit
  // logic out of this file while still allowing geometry / list updates.
  ctx.handles.getEditSurface = () => ({
    root,
    panel2d,
    pdbLower: PDB_LOWER,
    getAnnotations: () => fr3dData.annotations || [],
    setAnnotations: (anns) => {
      fr3dData.annotations = anns;
      fr3dPairLookup = null;
    },
    findBPPath,
    getResidueLocation,
    bindBpPathClick: (el) => {
      if (!el || el.dataset.r2dtBpClickBound) return;
      el.dataset.r2dtBpClickBound = '1';
      el.style.cursor = 'pointer';
      el.addEventListener('click', (ev) => {
        const cls = el.getAttribute('class') || '';
        const m = cls.match(/[a-zA-Z]{3}_(\d+)_(\d+)/);
        if (!m) return;
        ev.stopPropagation();
        selectBasePair(parseInt(m[1], 10), parseInt(m[2], 10), el);
      });
    },
    refreshAfterGeometryChange: () => {
      attachBPClicks();
      // Refresh BPs dropdown (hide empty families / show new ones) before
      // re-applying visibility so the badge and checkboxes match edits.
      syncBpFilterToAnnotations();
      // Re-apply family/nested filters using pair-key matching so refamily /
      // add paths are not force-hidden (their SVG class no longer matches the
      // plugin's original path IDs in displayBaseStrs).
      applyBpVisibility();
    },
    refreshBpListLabels: () => {
      applyBpListItemLabels();
      normalizeBpListScroll();
    },
    rebuildBpList: () => {
      const ui = rnaPlugin.uiTemplateService;
      if (ui?.renderBpListDialog) ui.renderBpListDialog(false);
      else {
        applyBpListItemLabels();
        normalizeBpListScroll();
      }
    },
    setOtherPairKeys,
    onPairSelect: (fn) => { pairSelectListeners.push(fn); },
    onResidueSelect: (fn) => { residueSelectListeners.push(fn); },
    getLayoutScale: () => getBaseBpStrokeWidth() * 6,
  });

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

    // No per-panel LBN panel here: createCompare renders one shared LBN
    // widget (both structures' rows, TP/FP/FN coloured) below the grid —
    // see buildCompareLbnDom / renderCompareLbn.

    slotEl.appendChild(root);
    return { root, panel2d, panel3d: null };
  }

  // ── Shared compare-mode LBN widget ──────────────────────────────────────
  // One dot-bracket widget for the whole compare page instead of one per
  // panel: every layer (WC, cWW, …) shows the reference row directly above
  // (or below) the model row for that same layer, brackets coloured by
  // TP (both structures) / FP (model only) / FN (reference only). Each row
  // still belongs to exactly one source structure, so a source's own
  // "Nested only" / family filters only greyscale that source's own rows —
  // see setVisibility below.

  function buildCompareLbnDom(compareRoot) {
    const lbnPanel = document.createElement('div');
    lbnPanel.className = 'r2dt-viewer-lbn r2dt-compare-lbn';
    lbnPanel.hidden = true;
    const lbnTitle = document.createElement('h2');
    lbnTitle.className = 'r2dt-viewer-lbn-title';
    lbnTitle.textContent = 'Layered dot-bracket notation (Leontis–Westhof base pairs)';
    const lbnCaption = document.createElement('p');
    lbnCaption.className = 'r2dt-viewer-lbn-caption';
    lbnCaption.innerHTML =
      'Brackets coloured by agreement: '
      + '<b class="lbn-legend lbn-bp--tp">TP both</b> '
      + '<b class="lbn-legend lbn-bp--fp">FP model only</b> '
      + '<b class="lbn-legend lbn-bp--fn">FN reference only</b>'
      + '. Each row belongs to one structure (labelled, left edge coloured) and '
      + 'follows that structure’s own "Nested only" / family filters.';
    const lbnBody = document.createElement('div');
    lbnBody.className = 'lbn-body';
    lbnPanel.append(lbnTitle, lbnCaption, lbnBody);
    compareRoot.appendChild(lbnPanel);
    return { panel: lbnPanel, body: lbnBody };
  }

  function classifierFor(bpCompare) {
    if (!bpCompare) return null;
    const otherPairKeys = new Set(bpCompare.otherPairKeys || []);
    return (a, b) => {
      const key = `${Math.min(a, b)}_${Math.max(a, b)}`;
      if (otherPairKeys.has(key)) return 'TP';
      return bpCompare.role === 'reference' ? 'FN' : 'FP';
    };
  }

  // sources: [{ panelIdx, label, lbnData, bpCompare }, …] one entry per
  // compare panel that has LBN data. selectInSource(panelIdx, pos, partner)
  // drives that panel's own selection (partner == null => plain residue).
  function renderCompareLbn(dom, sources, selectInSource) {
    const withRows = sources.filter(
      (s) => s.lbnData && s.lbnData.rows && s.lbnData.rows.length > 0
    );
    if (withRows.length === 0) return null;

    dom.panel.hidden = false;
    const SEQ = withRows[0].lbnData.sequence;
    const N = SEQ.length;

    // Row order: labels in first-seen order across sources, so same-layer
    // rows from different sources stay adjacent (ref row directly followed
    // by the model row for that layer).
    const labelOrder = [];
    const seenLabels = new Set();
    withRows.forEach((s) => {
      s.lbnData.rows.forEach((r) => {
        if (!seenLabels.has(r.label)) { seenLabels.add(r.label); labelOrder.push(r.label); }
      });
    });

    const renderRows = [];
    labelOrder.forEach((label) => {
      withRows.forEach((s) => {
        const row = s.lbnData.rows.find((r) => r.label === label);
        if (row) renderRows.push({ source: s, row, classify: classifierFor(s.bpCompare) });
      });
    });

    let html = '';
    html += '<div class="lbn-row"><span class="lbn-label">seq:</span> ';
    for (let i = 0; i < N; i++) {
      html += `<span data-pos="${i + 1}" class="lbn-nt">${SEQ[i]}</span>`;
    }
    html += '</div>';

    renderRows.forEach(({ source, row, classify }) => {
      const role = source.bpCompare?.role || null;
      const marker = role === 'reference' ? 'ref' : role === 'model' ? 'mdl' : '';
      const labelText = marker ? `${row.label} · ${marker}` : row.label;
      html += `<div class="lbn-row${role ? ` lbn-row--${role}` : ''}" data-panel-idx="${source.panelIdx}">`
        + `<span class="lbn-label${role ? ` lbn-label--${role}` : ''}" title="${source.label}">${labelText}:</span> `;
      for (let i = 0; i < N; i++) {
        const pos = i + 1;
        const ch = row.chars[i];
        if (ch === '.') {
          html += '<span class="lbn-dot">.</span>';
        } else {
          const partner = row.partners[String(pos)];
          const pAttr = partner != null ? ` data-partner="${partner}"` : '';
          const cls = partner != null && classify ? classify(pos, partner) : null;
          const clsAttr = cls ? ` lbn-bp--${cls.toLowerCase()}` : '';
          html += `<span data-pos="${pos}"${pAttr} class="lbn-bp${clsAttr}">${ch}</span>`;
        }
      }
      html += '</div>';
    });

    dom.body.innerHTML = html;

    const posSpans = {};
    dom.body.querySelectorAll('[data-pos]').forEach((sp) => {
      const p = +sp.dataset.pos;
      (posSpans[p] = posSpans[p] || []).push(sp);
    });

    let highlighted = [];
    function highlight(primary, pairs) {
      highlighted.forEach((sp) => sp.classList.remove('lbn-selected', 'lbn-partner'));
      highlighted = [];
      const prim = primary == null ? null : +primary;
      const edges = pairs || [];
      if (prim == null && edges.length === 0) return;

      const pairKeys = new Set(
        edges.map(([a, b]) => `${Math.min(a, b)}_${Math.max(a, b)}`)
      );
      const partnerPositions = new Set();
      edges.forEach(([a, b]) => {
        if (a !== prim) partnerPositions.add(a);
        if (b !== prim) partnerPositions.add(b);
      });

      dom.body.querySelectorAll('.lbn-nt').forEach((sp) => {
        const p = +sp.dataset.pos;
        let cls = null;
        if (p === prim) cls = 'lbn-selected';
        else if (partnerPositions.has(p)) cls = 'lbn-partner';
        if (cls) { sp.classList.add(cls); highlighted.push(sp); }
      });

      dom.body.querySelectorAll('.lbn-bp').forEach((sp) => {
        const p = +sp.dataset.pos;
        const partner = sp.dataset.partner ? +sp.dataset.partner : null;
        if (partner == null) return;
        const key = `${Math.min(p, partner)}_${Math.max(p, partner)}`;
        if (!pairKeys.has(key)) return;
        sp.classList.add(p === prim ? 'lbn-selected' : 'lbn-partner');
        highlighted.push(sp);
      });

      const focusSpan = (prim != null && posSpans[prim] && posSpans[prim][0]) || highlighted[0];
      if (focusSpan) {
        focusSpan.scrollIntoView({ behavior: 'smooth', block: 'nearest', inline: 'nearest' });
      }
    }

    dom.body.addEventListener('click', (ev) => {
      const sp = ev.target.closest('[data-pos]');
      if (!sp) return;
      const rowEl = sp.closest('[data-panel-idx]');
      // The shared seq row has no data-panel-idx; default to the first
      // source (mirroring propagates the selection to the other panels).
      const panelIdx = rowEl ? +rowEl.dataset.panelIdx : withRows[0].panelIdx;
      const pos = +sp.dataset.pos;
      const partner = sp.dataset.partner ? +sp.dataset.partner : null;
      selectInSource(panelIdx, pos, partner);
    });

    // Grey out only the rows belonging to the source whose 2D visibility
    // changed (own "Nested only" / family filters), leaving every other
    // source's rows untouched.
    function setVisibility(panelIdx, visibleKeys) {
      dom.body.querySelectorAll(`.lbn-row[data-panel-idx="${panelIdx}"]`).forEach((rowEl) => {
        const bps = rowEl.querySelectorAll('.lbn-bp');
        bps.forEach((sp) => {
          const pos = +sp.dataset.pos;
          const partner = sp.dataset.partner ? +sp.dataset.partner : null;
          if (partner == null) return;
          const key = `${Math.min(pos, partner)}_${Math.max(pos, partner)}`;
          const hidden = visibleKeys ? !visibleKeys.has(key) : false;
          sp.classList.toggle('lbn-bp-hidden', hidden);
        });
        const anyVisible = bps.length === 0
          || Array.from(bps).some((sp) => !sp.classList.contains('lbn-bp-hidden'));
        rowEl.classList.toggle('lbn-row--hidden', !anyVisible);
      });
    }

    return { highlight, setVisibility };
  }
  // ── end shared compare-mode LBN widget ──────────────────────────────────

  // The pdb-rna-viewer selection state (window.UiActionsService.selected) is a
  // single global Map shared by every panel, so when one panel clears it the
  // other panels' already-painted nucleotides are orphaned and never repainted
  // to their deselected colour.  Before mirroring a selection, hard-reset the
  // target panel's nucleotide colours so stale highlights can't accumulate.
  function resetPanelNucleotideColors(root) {
    if (!root) return;
    root.querySelectorAll('text.rnaviewEle').forEach((t) => {
      t.removeAttribute('fill');
      if (t.style) t.style.fill = '';
    });
    root.querySelectorAll('.r2dt-nt-selection-bg').forEach((r) => r.remove());
  }

  // `chainViews`: optional [{label, url, current}], for a reference structure
  // whose RNA chains don't all interact with each other — e.g. two unrelated
  // copies in the asymmetric unit. Picking an option navigates to that
  // chain's own reference-only page (a separate pre-rendered static page, not
  // a live swap — see r2dt.py's `_run_multichain_pdb`/partition_components).
  // Chains that *do* interact are always shown together in one view instead,
  // with no selector.
  function buildCompareSlot(title, subtitle, chainViews) {
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
    if (Array.isArray(chainViews) && chainViews.length > 1) {
      const select = document.createElement('select');
      select.className = 'r2dt-chain-select';
      select.title = 'This structure has independent (non-interacting) chains — pick one to view';
      chainViews.forEach((view) => {
        const option = document.createElement('option');
        option.value = view.url;
        option.textContent = view.label;
        if (view.current) option.selected = true;
        select.appendChild(option);
      });
      select.addEventListener('change', () => {
        if (select.value) window.location.href = select.value;
      });
      slot.appendChild(select);
    }
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
    // Workstation-only editing: probe + dynamic script load. Static hosts
    // never expose /__edit-api, so this is a no-op there.
    maybeEnableWorkstationEditing([ctx]).catch(() => {});
    if (typeof opts.onReady === 'function') opts.onReady(handle);
    return handle;
  }
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
    const lbnSources = [];
    // Mutable indirection so each panel's ctx.onBpVisibilityChange can be
    // wired before the shared LBN widget exists (it's built after every
    // panel has rendered, since row order needs every source's labels).
    const mergedLbn = { highlight() {}, setVisibility() {} };

    for (let i = 0; i < panels.length; i++) {
      const pOpts = panels[i];
      const panelData = await resolvePanelData(pOpts.baseUrl || '.', pOpts);
      const slot = buildCompareSlot(pOpts.title || panelData.structureId, pOpts.subtitle || '', pOpts.chainViews);
      grid.appendChild(slot);
      const dom = buildCompare2dDom(slot, { height: opts.panelHeight });
      const panelIdx = i;
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
        lbnData: panelData.lbnData,
        bpCompare: panelData.bpCompare,
        baseColor: pOpts.baseColor || null,
        showLbn: false,
        link3d: false,
        cleanup: [],
        handles: null,
        onBpVisibilityChange: (keys) => mergedLbn.setVisibility(panelIdx, keys),
      };
      // Pin this panel as active so the plugin's global DOM lookups during
      // init (svg.rnaTopoSvg, zoom setup, …) resolve inside this panel.
      const prevRoot = global.__r2dtSetActiveRoot?.(dom.root);
      try {
        await initViewer(ctx);
      } finally {
        global.__r2dtSetActiveRoot?.(prevRoot);
      }
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
      lbnSources.push({
        panelIdx,
        label: pOpts.title || panelData.structureId,
        lbnData: panelData.lbnData,
        bpCompare: panelData.bpCompare,
      });
    }

    // One shared LBN widget below the grid (both structures' rows, TP/FP/FN
    // coloured) instead of a widget per panel.
    const lbnDom = buildCompareLbnDom(compareRoot);
    const compareLbnApi = renderCompareLbn(lbnDom, lbnSources, (panelIdx, pos, partner) => {
      const handles = panelCtxs[panelIdx]?.handles;
      if (!handles) return;
      if (partner != null) handles.selectBasePair?.(pos, partner);
      else handles.selectResidue?.(pos);
    });
    if (compareLbnApi) {
      mergedLbn.highlight = compareLbnApi.highlight;
      mergedLbn.setVisibility = compareLbnApi.setVisibility;
      // Seed each source's initial visibility (current "Nested only" state).
      panelCtxs.forEach((c, i) => {
        const keys = c.handles?.getVisiblePairKeys?.();
        if (keys) compareLbnApi.setVisibility(i, keys);
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
      // Not used for the 3D pane itself; unset so resolvePanelData doesn't
      // try (and fail) to fetch a bp-compare.json that may not exist at
      // this baseUrl (e.g. when molstarOpts.baseUrl differs from the
      // linked panel's own baseUrl).
      bpCompare: null,
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

    // Structure 1 is the reference (loaded above). Its label→auth / label→chain
    // maps come from the linked 2D panel.
    const refBaseColor = molstarOpts.baseColor || null;
    const molTargets = [{
      structureNumber: 1,
      structureId: molData.structureId,
      labelToAuth: linkedCtx.handles.labelToAuth,
      labelToChain: linkedCtx.handles.labelToChain,
      chainId: molData.chainId,
      baseColor: refBaseColor,
      highlightColor: highlightColorFor(refBaseColor),
    }];
    const toggleEntries = [{
      structureNumber: 1,
      label: molstarOpts.title || molData.structureId,
      color: refBaseColor,
    }];

    // Load each overlay (e.g. the pre-aligned predicted model) as an additional
    // superimposed structure in the same canvas.
    const overlays = molstarOpts.overlays || [];
    for (let k = 0; k < overlays.length; k++) {
      const ov = overlays[k];
      const structureNumber = k + 2;
      const ovUrl = resolveUrl(molBaseUrl, ov.structureUrl);
      const ovFormat = ov.structureFormat === 'pdb' ? 'pdb' : 'mmcif';
      try {
        await molstar.visual.update(
          { customData: { url: ovUrl, format: ovFormat, binary: false } },
          false
        );
      } catch (err) {
        console.error('R2DTViewer.createCompare: overlay load failed', err);
      }
      let ovLabelToAuth = ov.labelToAuth || null;
      let ovLabelToChain = ov.labelToChain || null;
      if (!ovLabelToAuth || !ovLabelToChain) {
        try {
          const maps = await fetchJson(resolveUrl(ov.baseUrl || molBaseUrl, 'label-maps.json'));
          ovLabelToAuth = ovLabelToAuth || maps.labelToAuth || {};
          ovLabelToChain = ovLabelToChain || maps.labelToChain || {};
        } catch (_) { /* optional */ }
      }
      molTargets.push({
        structureNumber,
        structureId: ov.structureId,
        labelToAuth: ovLabelToAuth || {},
        labelToChain: ovLabelToChain || {},
        chainId: ov.chainId || '',
        baseColor: ov.baseColor || null,
        highlightColor: highlightColorFor(ov.baseColor || null),
      });
      toggleEntries.push({
        structureNumber,
        label: ov.title || ov.structureId,
        color: ov.baseColor || null,
      });
    }

    // Give each structure its distinct base colour (empty selection +
    // nonSelectedColor colours the whole structure) so the overlay reads apart.
    for (let i = 0; i < molTargets.length; i++) {
      const t = molTargets[i];
      if (!t.baseColor) continue;
      try {
        await molstar.visual.select({
          data: [],
          nonSelectedColor: t.baseColor,
          structureNumber: t.structureNumber,
          keepRepresentations: true,
        });
      } catch (_) { /* best-effort */ }
    }

    // Loading the overlay via visual.update resets the canvas background to the
    // Mol* default (black); re-assert white to match the init-time bgColor.
    if (overlays.length) {
      try { molstar.canvas.setBgColor({ r: 255, g: 255, b: 255 }); } catch (_) {}
    }

    if (overlays.length) addStructureToggles(molSlot, molstar, toggleEntries);
    // Loading a second structure via visual.update re-shows Mol*'s side controls
    // panel (the init-time hideControls no longer applies). It can't be hidden via
    // the layout state post-init (setProps doesn't re-render the embedded root), so
    // the compare CSS hides .msp-layout-right and expands .msp-layout-main instead.

    const selectInMolstar = createMultiMolstarSelector(molstar, molTargets);

    // Index panels by their (lower-cased) structure id so a selection event
    // can be attributed to its source panel and mirrored to the others.
    const idxByPdb = {};
    panelCtxs.forEach((c, i) => { idxByPdb[String(c.PDB_LOWER).toLowerCase()] = i; });
    // Map each panel's structure id to the 3D structureNumber it owns, so a
    // click in a 2D panel frames that panel's own structure in the shared 3D
    // (both structures are highlighted, but the camera follows the clicked one).
    const structNumById = {};
    molTargets.forEach((t) => {
      if (t.structureId != null) structNumById[String(t.structureId).toLowerCase()] = t.structureNumber;
    });
    let syncing = false;

    function onAnyResidueSelect(e) {
      if (syncing) return;
      const d = e.detail || {};
      const srcIdx = idxByPdb[(d.pdbId || '').toLowerCase()];
      if (srcIdx == null) return;

      // Drive the shared 3D from whichever panel changed: both structures'
      // residues light up (each in its own colour); the camera frames the
      // clicked panel's own structure.
      const srcStructNum = structNumById[(d.pdbId || '').toLowerCase()];
      if (d.cleared) selectInMolstar([]);
      else selectInMolstar([d.label, ...(d.partners || [])], srcStructNum);

      // Mirror the selection into the other 2D panel(s). For a base-pair click,
      // highlight the *same two positions* (co-indexed structures share the
      // label space) so clicking pair (a,b) in one panel selects (a,b) in the
      // other — not a re-derived partner. For a single-nucleotide click, mirror
      // the residue and let each panel recompute its own partners, so
      // agreements/disagreements show against its own pairs.
      // Pin the target panel active around each op so the plugin's global DOM
      // lookups (clearSelection/select) hit that panel, not the source panel.
      syncing = true;
      const prevRoot = global.__r2dtSetActiveRoot?.(null);
      try {
        panelCtxs.forEach((c, i) => {
          if (i === srcIdx || !c.handles) return;
          global.__r2dtSetActiveRoot?.(c.root);
          // Clear prior state, then hard-reset colours to defeat the shared-Map
          // orphaning, then paint the mirrored selection afresh.
          c.handles.clearSelection?.();
          resetPanelNucleotideColors(c.root);
          if (d.cleared) return;
          if (d.pair && d.partners && d.partners.length) {
            c.handles.selectBasePair?.(d.label, d.partners[0]);
          } else {
            c.handles.selectResidue?.(d.label);
          }
        });
      } finally {
        global.__r2dtSetActiveRoot?.(prevRoot);
        syncing = false;
      }

      // Shared LBN widget: light up this position/pair regardless of which
      // source row(s) it appears in (position is shared across all rows).
      const edges = d.cleared ? [] : (d.partners || []).map((p) => [d.label, p]);
      mergedLbn.highlight(d.cleared ? null : d.label, edges);
    }
    document.addEventListener('r2dt-residue-select', onAnyResidueSelect);
    cleanup.push(() => {
      document.removeEventListener('r2dt-residue-select', onAnyResidueSelect);
    });

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
        return linkedCtx.handles?.selectResidue?.(label);
      },
      selectBasePair(a, b) {
        return linkedCtx.handles?.selectBasePair?.(a, b);
      },
    };
    activeViewer = handle;
    // Workstation-only editing: probe + dynamic script load. Static hosts
    // never expose /__edit-api, so this is a no-op there.
    maybeEnableWorkstationEditing(panelCtxs).catch(() => {});
    if (typeof opts.onReady === 'function') opts.onReady(handle);
    return handle;
  }

  function loadScriptOnce(src) {
    return new Promise((resolve, reject) => {
      if (document.querySelector(`script[data-r2dt-src="${src}"]`)) {
        resolve();
        return;
      }
      const s = document.createElement('script');
      s.src = src;
      s.async = true;
      s.dataset.r2dtSrc = src;
      s.onload = () => resolve();
      s.onerror = () => reject(new Error('Failed to load ' + src));
      document.head.appendChild(s);
    });
  }

  async function maybeEnableWorkstationEditing(panelCtxs) {
    let ping;
    try {
      const res = await fetch('/__edit-api/ping', { method: 'GET' });
      if (!res.ok) return;
      ping = await res.json().catch(() => ({ ok: true }));
    } catch (_) {
      return;
    }
    if (!ping || ping.ok === false) return;
    const jobMatch = typeof location !== 'undefined'
      ? location.pathname.match(/\/jobs\/([^/]+)\/viewer\/?/)
      : null;
    const jobId = jobMatch && jobMatch[1];
    if (!jobId) return;
    await loadScriptOnce('/static/r2dt-bp-edit.js');
    if (!global.R2DTBpEdit) return;
    if (panelCtxs.length >= 2
        && typeof global.R2DTBpEdit.attachCompare === 'function') {
      await global.R2DTBpEdit.attachCompare({ panelCtxs, jobId });
      return;
    }
    if (panelCtxs.length === 1
        && typeof global.R2DTBpEdit.attachSingle === 'function') {
      await global.R2DTBpEdit.attachSingle({ panelCtx: panelCtxs[0], jobId });
    }
  }

  global.R2DTViewer = { create, createCompare };
})(typeof window !== 'undefined' ? window : global);
