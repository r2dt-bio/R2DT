/**
 * R2DT base-pair editing — workstation only.
 *
 * Loaded dynamically after GET /__edit-api/ping succeeds. Never referenced by
 * static Cloudflare pages, so edit chrome cannot appear there.
 *
 * Public API: R2DTBpEdit.attachCompare(options)
 */
(function (global) {
  'use strict';

  var LW_FAMILIES = [
    'cWW', 'tWW', 'cWH', 'tWH', 'cHW', 'tHW',
    'cWS', 'tWS', 'cSW', 'tSW', 'cHS', 'tHS',
    'cSH', 'tSH', 'cSS', 'tSS',
    'ncWW', 'ntWW',
  ];

  // ----- INF / diff (mirrors utils/multichain.py) -----

  function pairKey(a, b) {
    var i = Math.min(+a, +b);
    var j = Math.max(+a, +b);
    return i + '_' + j;
  }

  function isCanonical(family) {
    return String(family || '').toUpperCase() === 'CWW';
  }

  function splitFamilies(pairs) {
    var wc = new Set();
    var nwc = new Set();
    var all = new Set();
    (pairs || []).forEach(function (p) {
      var key = pairKey(p.i, p.j);
      all.add(key);
      if (isCanonical(p.family)) wc.add(key);
      else nwc.add(key);
    });
    return { wc: wc, nwc: nwc, all: all };
  }

  function infScore(ref, model) {
    var tp = 0;
    ref.forEach(function (k) { if (model.has(k)) tp += 1; });
    var fp = 0;
    model.forEach(function (k) { if (!ref.has(k)) fp += 1; });
    var fn = 0;
    ref.forEach(function (k) { if (!model.has(k)) fn += 1; });
    if (tp === 0) {
      var empty = ref.size === 0 && model.size === 0;
      return {
        inf: empty ? null : 0,
        ppv: empty ? null : 0,
        sty: empty ? null : 0,
        tp: tp, fp: fp, fn: fn,
      };
    }
    var ppv = tp / (tp + fp);
    var sty = tp / (tp + fn);
    return { inf: Math.sqrt(ppv * sty), ppv: ppv, sty: sty, tp: tp, fp: fp, fn: fn };
  }

  function computeInf(refPairs, modelPairs) {
    var r = splitFamilies(refPairs);
    var m = splitFamilies(modelPairs);
    return {
      wc: infScore(r.wc, m.wc),
      nwc: infScore(r.nwc, m.nwc),
      all: infScore(r.all, m.all),
    };
  }

  function diffCounts(refPairs, modelPairs) {
    var ref = new Set((refPairs || []).map(function (p) { return pairKey(p.i, p.j); }));
    var model = new Set((modelPairs || []).map(function (p) { return pairKey(p.i, p.j); }));
    var matched = 0;
    var lost = 0;
    var added = 0;
    ref.forEach(function (k) {
      if (model.has(k)) matched += 1;
      else lost += 1;
    });
    model.forEach(function (k) {
      if (!ref.has(k)) added += 1;
    });
    return { matched: matched, lost: lost, added: added };
  }

  // ----- Override apply -----

  function annToPair(a) {
    return {
      i: +a.seq_id1,
      j: +a.seq_id2,
      family: a.bp || 'cWW',
      nt1: a.nt1 || a.unit1 || 'N',
      nt2: a.nt2 || a.unit2 || 'N',
      crossing: a.crossing,
    };
  }

  function cloneAnns(anns) {
    return JSON.parse(JSON.stringify(anns || []));
  }

  function findAnnIndex(anns, i, j) {
    for (var n = 0; n < anns.length; n++) {
      var a = +anns[n].seq_id1;
      var b = +anns[n].seq_id2;
      if ((a === i && b === j) || (a === j && b === i)) return n;
    }
    return -1;
  }

  function applyOverrides(baselineAnns, overrides) {
    var anns = cloneAnns(baselineAnns);
    (overrides || []).forEach(function (op) {
      var i = +op.i;
      var j = +op.j;
      if (i > j) { var t = i; i = j; j = t; }
      var idx = findAnnIndex(anns, i, j);
      if (op.action === 'delete') {
        if (idx >= 0) anns.splice(idx, 1);
      } else if (op.action === 'add') {
        if (idx < 0) {
          anns.push({
            seq_id1: String(i),
            seq_id2: String(j),
            '3d_id1': String(i),
            '3d_id2': String(j),
            nt1: 'N',
            nt2: 'N',
            unit1: 'N',
            unit2: 'N',
            bp: op.family || 'cWW',
            crossing: '0',
          });
        } else {
          anns[idx].bp = op.family || anns[idx].bp;
        }
      } else if (op.action === 'refamily') {
        if (idx >= 0) anns[idx].bp = op.to || op.family || anns[idx].bp;
      }
    });
    return anns;
  }

  function annsToPairs(anns) {
    return (anns || []).map(annToPair);
  }

  function pairKeysFromAnns(anns) {
    return (anns || []).map(function (a) {
      return pairKey(a.seq_id1, a.seq_id2);
    });
  }

  // ----- Panel DOM helpers (via edit surface) -----

  function lerp(a, b, t) { return (1 - t) * a + t * b; }

  function residueLetter(surface, seqId) {
    var el = surface.panel2d.querySelector(
      'text.rnaview_' + surface.pdbLower + '_' + seqId
    );
    if (!el) return 'N';
    var t = (el.textContent || 'N').trim();
    return t.charAt(0).toUpperCase() || 'N';
  }

  function isCanonicalWC(family, nt1, nt2) {
    if (String(family) !== 'cWW') return false;
    var a = String(nt1 || '').toUpperCase();
    var b = String(nt2 || '').toUpperCase();
    return (a === 'A' && (b === 'U' || b === 'T'))
      || (b === 'A' && (a === 'U' || a === 'T'))
      || (a === 'G' && b === 'C')
      || (a === 'C' && b === 'G');
  }

  // Normalise dropdown families the plugin doesn't know about.
  function drawFamily(family) {
    var f = String(family || 'cWW');
    if (f === 'ncWW') return 'cWW';
    if (f === 'ntWW') return 'tWW';
    return f;
  }

  /**
   * Build SVG path attrs matching pdb-rna-viewer's Leontis–Westhof glyphs.
   * Port of UiTemplateService.calcBaseStrs geometry (filled cis / open trans).
   */
  function lwGlyphSpec(family, p1, p2, nt1, nt2, scale) {
    var original = String(family || 'cWW');
    var s = drawFamily(original);

    // Plugin flips SH/SW/HW so the first edge is W/H/S in canonical order.
    var flip = { cSH: 1, tSH: 1, cSW: 1, tSW: 1, cHW: 1, tHW: 1 };
    var loc1 = p1;
    var loc2 = p2;
    var n1 = nt1;
    var n2 = nt2;
    if (flip[s]) {
      loc1 = p2;
      loc2 = p1;
      n1 = nt2;
      n2 = nt1;
      s = s.charAt(0) + s.charAt(2) + s.charAt(1);
    }

    var i = scale || 12;
    var d0 = loc1[0] + i / 2.5;
    var f0 = loc2[0] + i / 2.5;
    var v0 = loc1[1] - i / 2.5;
    var g0 = loc2[1] - i / 2.5;
    var y = Math.hypot(d0 - f0, v0 - g0) || 1;
    var m = lerp(d0, f0, i / y);
    var _ = lerp(v0, g0, i / y);
    var b = lerp(d0, f0, 1 - i / y);
    var w = lerp(v0, g0, 1 - i / y);
    var fill = s.charAt(0) === 't' ? 'none' : '#ccc';
    var stroke = '#555';
    var sw = i / 6;
    var x = (m + b) / 2;
    var E = (_ + w) / 2;
    var k = y - 2 * i;
    var I = i / 1.5;
    var A = 270 + 180 * Math.atan2(v0 - g0, d0 - f0) / Math.PI;
    var transform = 'rotate(' + A + ' ' + x + ' ' + E + ')';

    if (s === 'cWW') {
      // Plain WC line for canonical pairs; line+circle+line for ncWW /
      // non-canonical cWW (matches viewer buildCwwNonCanonicalD).
      // Newly added pairs often have unknown letters (N) — treat as WC line
      // unless the user explicitly picked ncWW.
      var useCircle = original === 'ncWW' ||
        (!(n1 === 'N' && n2 === 'N') && !isCanonicalWC('cWW', n1, n2));
      if (!useCircle) {
        return {
          d: 'M' + m + ' ' + _ + ' ' + b + ' ' + w,
          fill: 'none',
          stroke: stroke,
          strokeWidth: sw,
          transform: null,
        };
      }
      var halfLen = Math.hypot(b - m, w - _) / 2 || 1;
      var frac = Math.min((i / 4) / halfLen, 0.45);
      var bx = lerp(m, x, 1 - frac);
      var by = lerp(_, E, 1 - frac);
      var nx = lerp(x, b, frac);
      var ny = lerp(E, w, frac);
      return {
        d: 'M ' + m + ' ' + _ + ' ' + bx + ' ' + by +
          ' M ' + (x - i / 4) + ' ' + E +
          ' a ' + (i / 4) + ',' + (i / 4) + ' 0 1,0 ' + (i / 2) + ',0' +
          ' a ' + (i / 4) + ',' + (i / 4) + ' 0 1,0 ' + (-i / 2) + ',0' +
          ' M ' + nx + ' ' + ny + ' ' + b + ' ' + w,
        fill: '#555',
        stroke: stroke,
        strokeWidth: sw,
        transform: null,
      };
    }

    if (s === 'tWW') {
      var B = lerp(m, x, 1 - (i / 3) / (y / 2));
      var T = lerp(_, E, 1 - (i / 3) / (y / 2));
      var N = lerp(x, b, (i / 3) / (y / 2));
      var C = lerp(E, w, (i / 3) / (y / 2));
      return {
        d: 'M ' + m + ' ' + _ + ' ' + B + ' ' + T +
          ' M ' + (x - i / 3) + ' ' + E +
          ' a ' + (i / 3) + ',' + (i / 3) + ' 0 1,0 ' + (i / 1.5) + ',0' +
          ' a ' + (i / 3) + ',' + (i / 3) + ' 0 1,0 ' + (-i / 1.5) + ',0' +
          ' M ' + N + ' ' + C + ' ' + b + ' ' + w,
        fill: fill,
        stroke: stroke,
        strokeWidth: sw,
        transform: null,
      };
    }

    var d;
    if (s === 'cSS' || s === 'tSS') {
      d = 'M ' + x + ' ' + (E + k / 2) + ' ' + x + ' ' + (E + I / 2) +
        ' l ' + (I / 2) + ' 0' +
        ' l -' + (I / 2) + ' -' + I +
        ' l -' + (I / 2) + ' ' + I +
        ' l ' + (I / 2) + ' 0' +
        ' M ' + x + ' ' + (E - I / 2) + ' ' + x + ' ' + (E - k / 2);
    } else if (s === 'cHS' || s === 'tHS') {
      d = 'M ' + x + ' ' + (E + k / 2) + ' ' + x + ' ' + (E + I + I / 4) +
        ' h -' + (I / 2) +
        ' v -' + I +
        ' h ' + I +
        ' v ' + I +
        ' h -' + (I / 2) +
        ' M ' + x + ' ' + (E + I / 4) + ' ' + x + ' ' + (E - I / 4) +
        ' l ' + (I / 2) + ' 0' +
        ' l -' + (I / 2) + ' -' + I +
        ' l -' + (I / 2) + ' ' + I +
        ' l ' + (I / 2) + ' 0' +
        ' M ' + x + ' ' + (E - I - I / 4) + ' ' + x + ' ' + (E - k / 2);
    } else if (s === 'cWS' || s === 'tWS') {
      d = 'M ' + x + ' ' + (E + k / 2) + ' ' + x + ' ' + (E + I + I / 4) +
        ' M ' + (x - I / 2) + ' ' + (E + 3 * I / 4) +
        ' a ' + (I / 2) + ',' + (I / 2) + ' 0 1,0 ' + I + ',0' +
        ' a ' + (I / 2) + ',' + (I / 2) + ' 0 1,0 ' + (-I) + ',0' +
        ' M ' + x + ' ' + (E + I / 4) + ' ' + x + ' ' + (E - I / 4) +
        ' l ' + (I / 2) + ' 0' +
        ' l -' + (I / 2) + ' -' + I +
        ' l -' + (I / 2) + ' ' + I +
        ' l ' + (I / 2) + ' 0' +
        ' M ' + x + ' ' + (E - I - I / 4) + ' ' + x + ' ' + (E - k / 2);
    } else if (s === 'cWH' || s === 'tWH') {
      d = 'M ' + x + ' ' + (E + k / 2) + ' ' + x + ' ' + (E + I + I / 4) +
        ' M ' + (x - I / 2) + ' ' + (E + 3 * I / 4) +
        ' a ' + (I / 2) + ',' + (I / 2) + ' 0 1,0 ' + I + ',0' +
        ' a ' + (I / 2) + ',' + (I / 2) + ' 0 1,0 ' + (-I) + ',0' +
        ' M ' + x + ' ' + (E + I / 4) + ' ' + x + ' ' + (E - I / 4) +
        ' h -' + (I / 2) +
        ' v -' + I +
        ' h ' + I +
        ' v ' + I +
        ' h -' + (I / 2) +
        ' M ' + x + ' ' + (E - I - I / 4) + ' ' + x + ' ' + (E - k / 2);
    } else if (s === 'cHH' || s === 'tHH') {
      d = 'M ' + x + ' ' + (E + k / 2) + ' ' + x + ' ' + (E + I / 2) +
        ' h -' + (I / 2) +
        ' v -' + I +
        ' h ' + I +
        ' v ' + I +
        ' h -' + (I / 2) +
        ' M ' + x + ' ' + (E - I / 2) + ' ' + x + ' ' + (E - k / 2);
    } else {
      // Fallback: dashed chord.
      return {
        d: 'M ' + m + ' ' + _ + ' ' + b + ' ' + w,
        fill: 'none',
        stroke: stroke,
        strokeWidth: sw,
        transform: null,
        dash: '4,3',
      };
    }

    return {
      d: d,
      fill: fill,
      stroke: stroke,
      strokeWidth: sw,
      transform: transform,
    };
  }

  function removePath(surface, i, j) {
    var el = surface.findBPPath(i, j);
    if (el && el.parentNode) el.parentNode.removeChild(el);
  }

  function drawPairPath(surface, i, j, family) {
    removePath(surface, i, j);
    var p1 = surface.getResidueLocation(i);
    var p2 = surface.getResidueLocation(j);
    if (!p1 || !p2) return null;
    var svg = surface.panel2d.querySelector('svg.rnaTopoSvg') ||
      surface.panel2d.querySelector('svg');
    if (!svg) return null;
    var g = svg.querySelector('g') || svg;
    var nt1 = residueLetter(surface, i);
    var nt2 = residueLetter(surface, j);
    var scale = (surface.getLayoutScale && surface.getLayoutScale()) || 12;
    var spec = lwGlyphSpec(family, p1, p2, nt1, nt2, scale);
    var path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
    // Keep original seq order in the class (not the flipped drawing order).
    var cls = 'rnaviewBP ' + drawFamily(family) + '_' + i + '_' + j;
    path.setAttribute('class', cls);
    path.setAttribute('d', spec.d);
    path.setAttribute('fill', spec.fill);
    path.setAttribute('stroke', spec.stroke);
    path.setAttribute('stroke-width', String(spec.strokeWidth));
    if (spec.transform) path.setAttribute('transform', spec.transform);
    if (spec.dash) path.setAttribute('stroke-dasharray', spec.dash);
    path.style.cursor = 'pointer';
    path.dataset.r2dtEdited = '1';
    path.dataset.r2dtFamily = String(family || 'cWW');
    g.appendChild(path);
    surface.bindBpPathClick(path);
    return path;
  }

  function syncPanelGeometry(surface, anns) {
    // Always rebuild edited geometry from annotations so family changes
    // redraw the correct Leontis–Westhof glyph (not just the class name).
    var want = {};
    (anns || []).forEach(function (a) {
      want[pairKey(a.seq_id1, a.seq_id2)] = a;
    });
    surface.panel2d.querySelectorAll('path.rnaviewBP').forEach(function (el) {
      var cls = el.getAttribute('class') || '';
      var m = cls.match(/([a-zA-Z]{3})_(\d+)_(\d+)/);
      if (!m) return;
      var key = pairKey(m[2], m[3]);
      if (!want[key]) {
        el.parentNode.removeChild(el);
        return;
      }
      // Redraw if family changed or this was an edited path.
      var fam = want[key].bp || m[1];
      var drawnFam = el.dataset.r2dtFamily || m[1];
      if (el.dataset.r2dtEdited || fam !== drawnFam) {
        el.parentNode.removeChild(el);
        return; // will be redrawn below
      }
      delete want[key];
    });
    Object.keys(want).forEach(function (key) {
      var a = want[key];
      drawPairPath(surface, +a.seq_id1, +a.seq_id2, a.bp || 'cWW');
    });
    surface.refreshAfterGeometryChange();
  }

  function rebuildBpList(surface) {
    if (surface.rebuildBpList) {
      surface.rebuildBpList();
      return;
    }
    surface.refreshBpListLabels();
  }

  // ----- UI chrome (per-panel, top-left of each 2D) -----

  var ICON_TRASH =
    '<svg class="r2dt-bp-edit-ico" viewBox="0 0 16 16" aria-hidden="true">' +
    '<path fill="currentColor" d="M6 2h4l.5 1H14v1.5H2V3h3.5L6 2zm1 4v6H5.5V6H7zm3.5 ' +
    '0v6H9V6h1.5zM3.5 4.5h9l-.7 9.2A1.5 1.5 0 0 1 10.3 15H5.7a1.5 1.5 0 0 1-1.5-1.3L3.5 ' +
    '4.5z"/></svg>';

  function infColour(value) {
    if (value == null) return '#9ca3af';
    if (value >= 0.95) return '#16a34a';
    if (value >= 0.85) return '#65a30d';
    if (value >= 0.60) return '#d97706';
    return '#dc2626';
  }

  function updateInfBar(metrics) {
    var bar = document.querySelector('.mc-inf');
    if (!bar || !metrics) return;
    ['wc', 'nwc', 'all'].forEach(function (key) {
      var label = key === 'wc' ? 'INF WC' : key === 'nwc' ? 'INF non-WC' : 'INF all';
      var item = null;
      bar.querySelectorAll('.mc-inf-item').forEach(function (el) {
        var k = el.querySelector('.mc-inf-key');
        if (k && k.textContent === label) item = el;
      });
      if (!item) return;
      var m = metrics[key] || {};
      var val = m.inf;
      var text = val == null ? 'n/a' : Number(val).toFixed(3);
      var v = item.querySelector('.mc-inf-val');
      if (v) {
        v.textContent = text;
        v.style.color = infColour(val);
      }
      item.title = (item.title || '').replace(/TP \d+.*/, '') +
        'TP ' + (m.tp || 0) + ', FP ' + (m.fp || 0) + ' (model-only), FN ' +
        (m.fn || 0) + ' (missing in model)';
    });
  }

  function fillFamilySelect(sel) {
    LW_FAMILIES.forEach(function (f) {
      var opt = document.createElement('option');
      opt.value = f;
      opt.textContent = f;
      sel.appendChild(opt);
    });
  }

  function buildPanelChrome(panelIdx) {
    var wrap = document.createElement('div');
    wrap.className = 'r2dt-bp-edit';
    wrap.dataset.panel = String(panelIdx);
    wrap.innerHTML =
      '<button type="button" class="pdb-rna-view-btn r2dt-bp-edit-toggle" ' +
        'data-act="toggle" title="Edit base pairs in this diagram">Edit</button>' +
      '<div class="r2dt-bp-edit-tools" role="toolbar" aria-label="Base-pair editing">' +
        '<button type="button" class="pdb-rna-view-btn r2dt-bp-edit-iconbtn" ' +
          'data-act="undo" title="Undo" aria-label="Undo" disabled>↶</button>' +
        '<button type="button" class="pdb-rna-view-btn r2dt-bp-edit-iconbtn" ' +
          'data-act="redo" title="Redo" aria-label="Redo" disabled>↷</button>' +
        '<button type="button" class="pdb-rna-view-btn" data-act="add" ' +
          'title="Add pair by clicking two nucleotides">Add</button>' +
        '<button type="button" class="pdb-rna-view-btn r2dt-bp-edit-iconbtn" ' +
          'data-act="delete" title="Delete selected pair" aria-label="Delete" disabled>' +
          ICON_TRASH + '</button>' +
        '<select class="r2dt-bp-edit-select" data-act="family" disabled ' +
          'title="Leontis–Westhof type" aria-label="Pair type"></select>' +
        '<button type="button" class="pdb-rna-view-btn" data-act="done" ' +
          'title="Finish editing this diagram">Done</button>' +
      '</div>';
    fillFamilySelect(wrap.querySelector('[data-act="family"]'));
    return wrap;
  }

  // ----- Controller -----

  function attachCompare(options) {
    var panelCtxs = options.panelCtxs || [];
    if (panelCtxs.length < 2) {
      return Promise.reject(new Error('R2DTBpEdit needs two compare panels'));
    }
    var surfaces = panelCtxs.map(function (ctx) {
      var surface = ctx.handles && ctx.handles.getEditSurface && ctx.handles.getEditSurface();
      if (!surface) throw new Error('Panel missing getEditSurface()');
      return surface;
    });
    var jobId = options.jobId;
    var baselines = surfaces.map(function (s) { return cloneAnns(s.getAnnotations()); });
    var overrides = { ref: [], model: [] };
    var undoStack = [];
    var redoStack = [];
    // Which 2D has its edit tools open (−1 = none). Only that panel accepts edits.
    var activePanel = -1;
    var addMode = false;
    var addFirst = null;
    var selected = null; // { panel, i, j, family }
    var saveTimer = null;
    var chromes = surfaces.map(function (surface, idx) {
      var chrome = buildPanelChrome(idx);
      surface.panel2d.appendChild(chrome);
      return chrome;
    });

    function panelName(idx) { return idx === 0 ? 'ref' : 'model'; }

    function workingAnns(idx) {
      var anns = applyOverrides(baselines[idx], overrides[panelName(idx)]);
      anns.forEach(function (a) {
        var n1 = residueLetter(surfaces[idx], a.seq_id1);
        var n2 = residueLetter(surfaces[idx], a.seq_id2);
        a.nt1 = n1;
        a.nt2 = n2;
        a.unit1 = n1;
        a.unit2 = n2;
      });
      return anns;
    }

    function activeChrome() {
      return activePanel >= 0 ? chromes[activePanel] : null;
    }

    function setStatus(text, kind) {
      void kind;
      chromes.forEach(function (chrome) {
        var tools = chrome.querySelector('.r2dt-bp-edit-tools');
        if (!tools) return;
        var base = tools.getAttribute('data-hint') || '';
        tools.title = [base, text].filter(Boolean).join(' · ');
      });
    }

    function setHint(text) {
      chromes.forEach(function (chrome) {
        var tools = chrome.querySelector('.r2dt-bp-edit-tools');
        if (!tools) return;
        tools.setAttribute('data-hint', text || '');
        var statusPart = (tools.title || '').split(' · ').slice(1).join(' · ');
        tools.title = [text || '', statusPart].filter(Boolean).join(' · ');
      });
    }

    function refreshChrome() {
      chromes.forEach(function (chrome, idx) {
        var open = activePanel === idx;
        chrome.classList.toggle('is-open', open);
        chrome.querySelector('[data-act="undo"]').disabled = !open || undoStack.length === 0;
        chrome.querySelector('[data-act="redo"]').disabled = !open || redoStack.length === 0;
        var selHere = !!(selected && selected.panel === idx);
        chrome.querySelector('[data-act="delete"]').disabled = !open || !selHere;
        var fam = chrome.querySelector('[data-act="family"]');
        fam.disabled = !open || !(addMode || selHere);
        if (selHere && selected.family && !addMode) fam.value = selected.family;
        chrome.querySelector('[data-act="add"]').classList.toggle('active', open && addMode);
      });
    }

    function closeEdit() {
      activePanel = -1;
      addMode = false;
      addFirst = null;
      selected = null;
      setHint('');
      refreshChrome();
    }

    function openEdit(panelIdx) {
      if (activePanel !== panelIdx) {
        addMode = false;
        addFirst = null;
        if (selected && selected.panel !== panelIdx) selected = null;
      }
      activePanel = panelIdx;
      setHint('Editing this 2D — Add, or select a pair to delete / change type');
      refreshChrome();
    }

    function pushHistory() {
      undoStack.push({
        ref: JSON.parse(JSON.stringify(overrides.ref)),
        model: JSON.parse(JSON.stringify(overrides.model)),
      });
      if (undoStack.length > 100) undoStack.shift();
      redoStack.length = 0;
    }

    function commit(mutator) {
      pushHistory();
      mutator();
      selected = null;
      addMode = false;
      addFirst = null;
      setHint('');
      refreshAll();
      scheduleSave();
      refreshChrome();
    }

    function doUndo() {
      if (!undoStack.length) return;
      redoStack.push({
        ref: JSON.parse(JSON.stringify(overrides.ref)),
        model: JSON.parse(JSON.stringify(overrides.model)),
      });
      var prev = undoStack.pop();
      overrides.ref = prev.ref;
      overrides.model = prev.model;
      selected = null;
      addMode = false;
      addFirst = null;
      setHint('');
      refreshAll();
      scheduleSave();
      refreshChrome();
    }

    function doRedo() {
      if (!redoStack.length) return;
      undoStack.push({
        ref: JSON.parse(JSON.stringify(overrides.ref)),
        model: JSON.parse(JSON.stringify(overrides.model)),
      });
      var next = redoStack.pop();
      overrides.ref = next.ref;
      overrides.model = next.model;
      selected = null;
      addMode = false;
      addFirst = null;
      setHint('');
      refreshAll();
      scheduleSave();
      refreshChrome();
    }

    function refreshAll() {
      var anns0 = workingAnns(0);
      var anns1 = workingAnns(1);
      surfaces[0].setAnnotations(anns0);
      surfaces[1].setAnnotations(anns1);
      surfaces[0].setOtherPairKeys(pairKeysFromAnns(anns1));
      surfaces[1].setOtherPairKeys(pairKeysFromAnns(anns0));
      syncPanelGeometry(surfaces[0], anns0);
      syncPanelGeometry(surfaces[1], anns1);
      rebuildBpList(surfaces[0]);
      rebuildBpList(surfaces[1]);
      var infMetrics = computeInf(annsToPairs(anns0), annsToPairs(anns1));
      updateInfBar(infMetrics);
      var diff = diffCounts(annsToPairs(anns0), annsToPairs(anns1));
      setStatus(
        'edits ' + overrides.ref.length + '/' + overrides.model.length +
        ' · BP ' + diff.matched + '/' + diff.lost + '/' + diff.added,
        'ok'
      );
    }

    function scheduleSave() {
      setStatus('Saving…', 'busy');
      if (saveTimer) clearTimeout(saveTimer);
      saveTimer = setTimeout(doSave, 120);
    }

    function doSave() {
      return fetch('/__edit-api/jobs/' + encodeURIComponent(jobId) + '/basepairs', {
        method: 'PUT',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(overrides),
      }).then(function (r) {
        return r.json().then(function (body) {
          if (!r.ok) throw new Error(body.error || ('save failed ' + r.status));
          setStatus('Saved', 'ok');
        });
      }).catch(function (err) {
        setStatus(err.message || String(err), 'err');
      });
    }

    function setOverridesForPair(panelIdx, i, j, nextAction) {
      var panel = panelName(panelIdx);
      var baseline = baselines[panelIdx];
      var i0 = Math.min(+i, +j);
      var j0 = Math.max(+i, +j);
      var baseIdx = findAnnIndex(baseline, i0, j0);
      var baseFam = baseIdx >= 0 ? (baseline[baseIdx].bp || 'cWW') : null;
      overrides[panel] = overrides[panel].filter(function (op) {
        var a = Math.min(+op.i, +op.j);
        var b = Math.max(+op.i, +op.j);
        return !(a === i0 && b === j0);
      });
      if (nextAction.action === 'delete') {
        if (baseIdx >= 0) {
          overrides[panel].push({
            action: 'delete', i: i0, j: j0, family: baseFam,
          });
        }
      } else if (nextAction.action === 'add') {
        if (baseIdx >= 0 && baseFam === nextAction.family) {
          // equals baseline — no override
        } else if (baseIdx >= 0) {
          overrides[panel].push({
            action: 'refamily', i: i0, j: j0, from: baseFam, to: nextAction.family,
          });
        } else {
          overrides[panel].push({
            action: 'add', i: i0, j: j0, family: nextAction.family,
          });
        }
      } else if (nextAction.action === 'refamily') {
        if (baseIdx < 0) {
          overrides[panel].push({
            action: 'add', i: i0, j: j0, family: nextAction.to,
          });
        } else if (baseFam === nextAction.to) {
          // back to baseline
        } else {
          overrides[panel].push({
            action: 'refamily', i: i0, j: j0, from: baseFam, to: nextAction.to,
          });
        }
      }
    }

    function onPairSelected(panelIdx, i, j, family) {
      if (activePanel < 0 || panelIdx !== activePanel) return;
      if (addMode) return;
      selected = { panel: panelIdx, i: +i, j: +j, family: family || 'cWW' };
      refreshChrome();
    }

    function onResidueClicked(panelIdx, label) {
      if (activePanel < 0 || !addMode || panelIdx !== activePanel) return;
      var nt = +label;
      if (!addFirst) {
        addFirst = nt;
        setHint('Selected ' + nt + ' — click partner nucleotide');
        return;
      }
      if (addFirst === nt) {
        setHint('Pick a different nucleotide');
        return;
      }
      var i = Math.min(addFirst, nt);
      var j = Math.max(addFirst, nt);
      var chrome = activeChrome();
      var fam = (chrome && chrome.querySelector('[data-act="family"]').value) || 'cWW';
      commit(function () {
        setOverridesForPair(activePanel, i, j, { action: 'add', family: fam });
      });
      selected = { panel: activePanel, i: i, j: j, family: fam };
      setHint('Added ' + i + '–' + j + ' as ' + fam);
      refreshChrome();
    }

    surfaces.forEach(function (surface, idx) {
      surface.onPairSelect(function (i, j, family) {
        onPairSelected(idx, i, j, family);
      });
      surface.onResidueSelect(function (label) {
        onResidueClicked(idx, label);
      });
    });

    chromes.forEach(function (chrome, idx) {
      chrome.addEventListener('click', function (ev) {
        var btn = ev.target.closest('[data-act]');
        if (!btn || btn.tagName === 'SELECT' || btn.tagName === 'LABEL') return;
        var act = btn.getAttribute('data-act');
        if (act === 'toggle') {
          openEdit(idx);
          return;
        }
        if (act === 'done') {
          closeEdit();
          return;
        }
        if (activePanel !== idx) return;
        if (act === 'undo') {
          doUndo();
        } else if (act === 'redo') {
          doRedo();
        } else if (act === 'add') {
          addMode = !addMode;
          addFirst = null;
          selected = null;
          setHint(addMode
            ? 'Pick Type, then click two nucleotides in this 2D'
            : 'Editing this 2D — Add, or select a pair to delete / change type');
          refreshChrome();
        } else if (act === 'delete') {
          if (!selected || selected.panel !== activePanel) return;
          var si = selected.i;
          var sj = selected.j;
          commit(function () {
            setOverridesForPair(activePanel, si, sj, { action: 'delete' });
          });
        }
      });

      chrome.querySelector('[data-act="family"]').addEventListener('change', function (ev) {
        if (activePanel !== idx) return;
        if (!selected || selected.panel !== activePanel) return;
        var to = ev.target.value;
        var si = selected.i;
        var sj = selected.j;
        commit(function () {
          setOverridesForPair(activePanel, si, sj, { action: 'refamily', to: to });
        });
        selected = { panel: activePanel, i: si, j: sj, family: to };
        refreshChrome();
      });
    });

    document.addEventListener('keydown', function (ev) {
      if (activePanel < 0) return;
      var mod = ev.metaKey || ev.ctrlKey;
      if (!mod) return;
      var key = ev.key.toLowerCase();
      if (key === 'z' && !ev.shiftKey) {
        ev.preventDefault();
        doUndo();
      } else if ((key === 'z' && ev.shiftKey) || key === 'y') {
        ev.preventDefault();
        doRedo();
      }
    });

    return fetch('/__edit-api/jobs/' + encodeURIComponent(jobId) + '/basepairs')
      .then(function (r) { return r.json(); })
      .then(function (data) {
        overrides.ref = data.ref || [];
        overrides.model = data.model || [];
        refreshAll();
        refreshChrome();
        setStatus('Click Edit on a 2D to curate pairs', 'ok');
      });
  }

  function jobIdFromLocation() {
    var m = (global.location && location.pathname || '')
      .match(/\/jobs\/([^/]+)\/viewer\/?/);
    return m ? m[1] : null;
  }

  function probe() {
    return fetch('/__edit-api/ping', { method: 'GET' })
      .then(function (r) {
        if (!r.ok) return null;
        return r.json().catch(function () { return { ok: true }; });
      })
      .catch(function () { return null; });
  }

  global.R2DTBpEdit = {
    probe: probe,
    jobIdFromLocation: jobIdFromLocation,
    attachCompare: attachCompare,
    // exposed for tests
    _computeInf: computeInf,
    _applyOverrides: applyOverrides,
  };
}(typeof window !== 'undefined' ? window : globalThis));
