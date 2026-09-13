(function () {
  var state = {
    jobs: [],
    ref: null,
    model: null,
    suggestion: null,
    highlightId: null,
    pollTimer: null,
    sortKey: "created",
    sortDir: "desc",
  };
  var IDENT_THRESHOLD = 0.9;

  function $(id) { return document.getElementById(id); }

  function fmt(n) {
    if (n === null || n === undefined || n === "") return "—";
    if (typeof n === "number") return n.toFixed(3);
    return String(n);
  }

  function fmtDate(iso) {
    if (!iso) return "—";
    var d = new Date(iso);
    if (isNaN(d.getTime())) return String(iso);
    try {
      return new Intl.DateTimeFormat(undefined, {
        dateStyle: "medium",
        timeStyle: "short",
      }).format(d);
    } catch (_) {
      return d.toLocaleString();
    }
  }

  function showPanel(name) {
    document.querySelectorAll(".panel").forEach(function (p) {
      p.classList.toggle("active", p.id === "panel-" + name);
    });
  }

  function friendlyWorkspace(path) {
    if (!path) return "";
    // Inside the container the mount is /workspace; show the host-facing idea.
    if (path === "/workspace" || path.indexOf("/workspace") === 0) {
      return "~/.r2dt-workstation (mounted at /workspace in Docker)";
    }
    return path;
  }

  function loadRuntime() {
    return fetch("/api/runtime").then(function (r) { return r.json(); }).then(function (data) {
      var el = $("runtime");
      if (!data.ok) {
        el.textContent = "Runtime error: " + data.message;
        el.style.color = "var(--danger)";
      } else {
        var parts = [
          "Server running in Docker",
          "cache: " + friendlyWorkspace(data.workspace),
        ];
        if (data.current_job) parts.push("job running: " + data.current_job);
        el.textContent = parts.join(" · ");
      }
    });
  }

  function setDashStatus(text, cls) {
    var el = $("import-status");
    if (!el) return;
    el.className = "form-status" + (cls ? " " + cls : "");
    el.textContent = text || "";
  }

  function loadJobs() {
    var tbody = $("rows");
    if (tbody && !tbody.children.length) {
      tbody.innerHTML = '<tr><td colspan="12">Loading jobs…</td></tr>';
    }
    return fetch("/api/jobs?mode=compare").then(function (r) { return r.json(); }).then(function (data) {
      state.jobs = data.jobs || [];
      renderJobs();
      maybePoll();
    }).catch(function (err) {
      if (!tbody) tbody = $("rows");
      if (tbody) {
        tbody.innerHTML = '<tr><td colspan="12">Could not load jobs: ' +
          esc(err && err.message ? err.message : err) + "</td></tr>";
      }
    });
  }

  function seqNumbers(jobs) {
    // Stable # by created ascending (oldest = 1), independent of current sort.
    var ordered = jobs.slice().sort(function (a, b) {
      return String(a.created || "").localeCompare(String(b.created || "")) ||
        String(a.id || "").localeCompare(String(b.id || ""));
    });
    var map = {};
    ordered.forEach(function (j, i) { map[j.id] = i + 1; });
    return map;
  }

  function rowValue(j, key, seqMap) {
    var inf = (j.metrics && j.metrics.inf) || {};
    var params = j.params || {};
    var inputs = j.inputs || {};
    var m = j.metrics || {};
    switch (key) {
      case "seq": return seqMap[j.id] || 0;
      case "label": return (j.label || j.id || "").toLowerCase();
      case "ref": return (inputs.ref_name || "").toLowerCase();
      case "model": return (inputs.model_name || "").toLowerCase();
      case "chains": return String(params.chains || "");
      case "created": return j.created || "";
      case "inf_wc": return inf.wc == null ? -1 : +inf.wc;
      case "inf_nwc": return inf.nwc == null ? -1 : +inf.nwc;
      case "inf_all": return inf.all == null ? -1 : +inf.all;
      case "bp": return (m.matched || 0) + (m.lost || 0) + (m.added || 0);
      case "status": return j.status || "";
      default: return "";
    }
  }

  function sortedRows(rows, seqMap) {
    var key = state.sortKey;
    var dir = state.sortDir === "asc" ? 1 : -1;
    return rows.slice().sort(function (a, b) {
      var va = rowValue(a, key, seqMap);
      var vb = rowValue(b, key, seqMap);
      if (typeof va === "number" && typeof vb === "number") {
        if (va === vb) return 0;
        return va < vb ? -dir : dir;
      }
      va = String(va);
      vb = String(vb);
      return va < vb ? -dir : va > vb ? dir : 0;
    });
  }

  function updateSortHeaders() {
    document.querySelectorAll("#tbl thead th[data-sort]").forEach(function (th) {
      th.classList.remove("sort-asc", "sort-desc");
      if (th.getAttribute("data-sort") === state.sortKey) {
        th.classList.add(state.sortDir === "asc" ? "sort-asc" : "sort-desc");
      }
    });
  }

  function renderJobs() {
    var q = ($("filter").value || "").trim().toLowerCase();
    var filtered = state.jobs.filter(function (j) {
      if (!q) return true;
      var hay = [j.label, j.notes, j.id,
        (j.inputs || {}).ref_name, (j.inputs || {}).model_name,
        (j.params || {}).chains].join(" ").toLowerCase();
      return hay.indexOf(q) !== -1;
    });
    var seqMap = seqNumbers(state.jobs);
    var rows = sortedRows(filtered, seqMap);
    $("count").textContent = rows.length + " job" + (rows.length === 1 ? "" : "s");
    updateSortHeaders();
    var tbody = $("rows");
    tbody.innerHTML = "";
    if (!rows.length) {
      tbody.innerHTML = '<tr><td colspan="12">No comparisons yet. Use New comparison to start one.</td></tr>';
      return;
    }
    rows.forEach(function (j) {
      var inf = (j.metrics && j.metrics.inf) || {};
      var params = j.params || {};
      var inputs = j.inputs || {};
      var m = j.metrics || {};
      var bp = [m.matched, m.lost, m.added].every(function (x) { return x !== undefined && x !== null; })
        ? (m.matched + "/" + m.lost + "/" + m.added) : "—";
      var tr = document.createElement("tr");
      tr.innerHTML =
        '<td class="num">' + (seqMap[j.id] || "—") + "</td>" +
        labelCell(j) +
        '<td><span class="badge ' + esc(j.status || "") + '">' + esc(j.status || "") + "</span></td>" +
        '<td class="actions-col"><div class="row-actions"></div></td>' +
        "<td>" + esc(inputs.ref_name || "") + "</td>" +
        "<td>" + esc(inputs.model_name || "") + "</td>" +
        "<td>" + chainsCell(params, m) + "</td>" +
        '<td class="num">' + fmt(inf.wc) + "</td>" +
        '<td class="num">' + fmt(inf.nwc) + "</td>" +
        '<td class="num">' + fmt(inf.all) + "</td>" +
        "<td>" + bp + "</td>" +
        '<td title="' + esc(j.created || "") + '">' + esc(fmtDate(j.created)) + "</td>";
      var actions = tr.querySelector(".row-actions");
      if (window.R2DTTransfer) window.R2DTTransfer.addExportButton(actions, j);
      if (j.status === "running" || j.status === "queued" || j.status === "failed") {
        var logBtn = document.createElement("button");
        logBtn.type = "button";
        logBtn.textContent = "Log";
        logBtn.addEventListener("click", function () { showLog(j.id); });
        actions.appendChild(logBtn);
      }
      var del = document.createElement("button");
      del.type = "button";
      del.className = "danger";
      del.textContent = "Delete";
      del.addEventListener("click", function () { deleteJob(j.id); });
      actions.appendChild(del);
      if (state.highlightId && j.id === state.highlightId) {
        tr.className = "is-highlight";
      }
      tbody.appendChild(tr);
    });
    if (state.highlightId) {
      var hit = tbody.querySelector("tr.is-highlight");
      if (hit && hit.scrollIntoView) {
        hit.scrollIntoView({ block: "center" });
      }
    }
  }

  function chainsCell(params, metrics) {
    var map = esc(params.chains || "") +
      (params.model_chains ? " → " + esc(params.model_chains) : "");
    if (metrics && metrics.display_widened && metrics.display_chains) {
      var disp = Array.isArray(metrics.display_chains)
        ? metrics.display_chains.join(",")
        : String(metrics.display_chains);
      map +=
        ' <span class="hint" title="Reference display includes interacting ' +
        'chain(s) not scored against the model">(display ' + esc(disp) + ")</span>";
    }
    return map;
  }

  function esc(s) {
    return String(s == null ? "" : s)
      .replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;");
  }

  function labelCell(j) {
    var text = esc(j.label || j.id);
    if (j.status === "ready" && j.viewer_url) {
      return (
        '<td><a class="ws-label-link" href="' + esc(j.viewer_url) +
        '" target="_blank" rel="noopener">' + text + "</a></td>"
      );
    }
    return "<td>" + text + "</td>";
  }

  function maybePoll() {
    var busy = state.jobs.some(function (j) {
      return j.status === "running" || j.status === "queued";
    });
    if (state.pollTimer) {
      clearInterval(state.pollTimer);
      state.pollTimer = null;
    }
    if (busy) {
      state.pollTimer = setInterval(function () {
        loadJobs();
        loadRuntime();
      }, 2500);
    }
  }

  function showLog(id) {
    fetch("/api/jobs/" + encodeURIComponent(id) + "/log?tail=300")
      .then(function (r) { return r.json(); })
      .then(function (data) {
        var el = $("log");
        el.classList.remove("hidden");
        el.textContent = "— " + id + " —\n" + (data.log || "(empty)");
      });
  }

  function deleteJob(id) {
    if (!window.confirm("Delete job " + id + " from the local cache?")) return;
    fetch("/api/jobs/" + encodeURIComponent(id), { method: "DELETE" })
      .then(function (r) { return r.json(); })
      .then(function () { return loadJobs(); });
  }

  function selectedChains(containerId) {
    var boxes = $(containerId).querySelectorAll('input[type=checkbox]:checked');
    return Array.prototype.map.call(boxes, function (b) { return b.value; });
  }

  function chainDetails(info) {
    if (!info) return [];
    if (info.chain_details && info.chain_details.length) return info.chain_details;
    return (info.chains || []).map(function (id) {
      return { id: id, length: 0, preview: "", sequence: "" };
    });
  }

  function sumLen(details, ids) {
    var by = {};
    details.forEach(function (d) { by[d.id] = d.length || 0; });
    return ids.reduce(function (n, id) { return n + (by[id] || 0); }, 0);
  }

  function chainIdent(a, b) {
    if (!a || !b || a.length !== b.length) return 0;
    var n = 0;
    for (var i = 0; i < a.length; i++) {
      if (a[i] === b[i] || a[i] === "N" || b[i] === "N") n++;
    }
    return n / a.length;
  }

  function concatSeq(details, ids) {
    var by = {};
    details.forEach(function (d) { by[d.id] = d.sequence || ""; });
    return ids.map(function (id) { return by[id] || ""; }).join("");
  }

  function sameSet(a, b) {
    if (!a || !b || a.length !== b.length) return false;
    var left = a.slice().sort();
    var right = b.slice().sort();
    return left.every(function (id, i) { return id === right[i]; });
  }

  function mappingOrder(suggested, selected) {
    if (suggested && suggested.length && sameSet(suggested, selected)) {
      return suggested.slice();
    }
    return selected;
  }

  function packMapping(refIds, modelIds, reason, extra) {
    extra = extra || {};
    return {
      ref_chains: refIds,
      model_chains: modelIds,
      reason: reason,
      summary: refIds.join("+") + " → " + modelIds.join("+"),
      alternatives: extra.alternatives || [],
      identity: extra.identity,
      by: extra.by || "sequence",
    };
  }

  // Keep in sync with utils.workstation.chains.suggest_chain_mapping.
  function suggestChainMapping(refD, modelD) {
    if (!refD || !modelD || !refD.length || !modelD.length) return null;
    var refIds = refD.map(function (d) { return d.id; });
    var modIds = modelD.map(function (d) { return d.id; });
    var modelSeq = concatSeq(modelD, modIds);
    var modelLen = modelSeq.length || sumLen(modelD, modIds);

    if (refD.length === modelD.length) {
      var ok = true;
      var worst = 1;
      for (var i = 0; i < refD.length; i++) {
        var pair = chainIdent(refD[i].sequence, modelD[i].sequence);
        if (pair < IDENT_THRESHOLD) { ok = false; break; }
        if (pair < worst) worst = pair;
      }
      if (ok) {
        return packMapping(
          refIds, modIds,
          "Each reference chain matches the model chain in the same order.",
          { identity: worst }
        );
      }
    }

    var singles = [];
    refD.forEach(function (d) {
      var score = chainIdent(d.sequence, modelSeq);
      if (score >= IDENT_THRESHOLD) {
        singles.push({ id: d.id, identity: score, length: d.length || (d.sequence || "").length });
      }
    });
    if (singles.length) {
      var bestId = singles[0].identity;
      singles.forEach(function (hit) {
        if (hit.identity > bestId) bestId = hit.identity;
      });
      singles = singles.filter(function (hit) { return hit.identity === bestId; });
      var best = singles[0];
      var alts = singles.slice(1).map(function (hit) {
        return packMapping(
          [hit.id], modIds, "Same sequence as chain " + best.id + ".",
          { identity: hit.identity }
        );
      });
      var reason = "Model " + (modIds.length === 1 ? "chain " + modIds[0] : "chains " + modIds.join("+")) +
        " (" + (best.length || modelLen) + " nt) matches reference chain " + best.id + ".";
      if (alts.length) {
        var altIds = alts.map(function (item) { return item.ref_chains[0]; });
        if (altIds.length === 1) {
          reason += " Chain " + altIds[0] +
            " has the same sequence — pick that chain to score the other monomer. " +
            "Interacting partners can still appear in the display.";
        } else {
          reason += " Chains " + altIds.join(", ") +
            " have the same sequence — pick one of those to score the other monomer. " +
            "Interacting partners can still appear in the display.";
        }
      }
      return packMapping([best.id], modIds, reason, { identity: best.identity, alternatives: alts });
    }

    var width, start, ids, score;
    for (width = 2; width <= refD.length; width++) {
      for (start = 0; start + width <= refD.length; start++) {
        ids = refIds.slice(start, start + width);
        score = chainIdent(concatSeq(refD, ids), modelSeq);
        if (score >= IDENT_THRESHOLD) {
          return packMapping(
            ids, modIds,
            "Model " + (modIds.length === 1 ? "chain " + modIds[0] : "chains " + modIds.join("+")) +
            " (" + modelLen + " nt) matches reference " + ids.join("+") + " concatenated.",
            { identity: score }
          );
        }
      }
    }

    if (modelD.length > 1) {
      var chosen = [];
      var used = {};
      var okG = true;
      var worstG = 1;
      for (var mi = 0; mi < modelD.length; mi++) {
        var pickId = null;
        var pickS = -1;
        for (var ri = 0; ri < refD.length; ri++) {
          if (used[refD[ri].id]) continue;
          var sc = chainIdent(refD[ri].sequence, modelD[mi].sequence);
          if (sc >= IDENT_THRESHOLD && sc > pickS) {
            pickId = refD[ri].id;
            pickS = sc;
          }
        }
        if (!pickId) { okG = false; break; }
        used[pickId] = true;
        chosen.push(pickId);
        if (pickS < worstG) worstG = pickS;
      }
      if (okG) {
        return packMapping(
          chosen, modIds,
          "Each model chain matches a distinct reference chain.",
          { identity: worstG }
        );
      }
    }

    var lenHits = refD.filter(function (d) {
      return (d.length || 0) === modelLen && modelLen > 0;
    });
    if (lenHits.length) {
      var first = lenHits[0];
      var altsL = lenHits.slice(1).map(function (d) {
        return packMapping(
          [d.id], modIds, "Same length as chain " + first.id + ".", { by: "length" }
        );
      });
      var rsn = "Reference chain " + first.id + " and the model are both " + modelLen +
        " nt (sequences differ). Confirm this is the intended pair.";
      if (altsL.length) {
        rsn += " Chain" + (altsL.length > 1 ? "s " : " ") +
          altsL.map(function (item) { return item.ref_chains[0]; }).join(", ") +
          " have the same length.";
      }
      return packMapping(
        [first.id], modIds, rsn, { by: "length", alternatives: altsL }
      );
    }

    for (width = 2; width <= refD.length; width++) {
      for (start = 0; start + width <= refD.length; start++) {
        ids = refIds.slice(start, start + width);
        if (sumLen(refD, ids) === modelLen && modelLen > 0) {
          return packMapping(
            ids, modIds,
            "Reference " + ids.join("+") + " and the model are both " + modelLen +
            " nt concatenated. Sequences were not an automatic match — check the order.",
            { by: "length" }
          );
        }
      }
    }
    return null;
  }

  function chainLabelInner(d) {
    var meta = [];
    if (d.length) meta.push(d.length + " nt");
    if (d.auth_start && d.auth_end) {
      meta.push("residues " + d.auth_start + "–" + d.auth_end);
    }
    if (d.same_sequence_as && d.same_sequence_as.length) {
      meta.push("same sequence as " + d.same_sequence_as.map(function (id) {
        return "chain " + id;
      }).join(", "));
    }
    var html = '<span class="chain-meta"><span class="chain-id">chain ' +
      esc(d.id) + "</span>";
    if (meta.length) {
      html += '<span class="chain-sub">' + esc(meta.join(" · ")) + "</span>";
    }
    if (d.preview) {
      html += '<span class="chain-seq">' + esc(d.preview) + "</span>";
    }
    return html + "</span>";
  }

  function renderChainPicker(containerId, info, role, suggestedIds) {
    var el = $(containerId);
    if (!info) {
      el.classList.add("hidden");
      el.innerHTML = "";
      return;
    }
    el.classList.remove("hidden");
    var details = chainDetails(info);
    var html = "";
    if (role === "ref" && info.needs_cif_conversion) {
      html += '<p class="hint">PDB reference will be converted to mmCIF for compare.</p>';
    }
    if (!details.length) {
      html += '<p class="err">No RNA chains detected.</p>';
      el.innerHTML = html;
      return;
    }
    html += '<p class="hint">' + details.length +
      " RNA chain" + (details.length === 1 ? "" : "s") +
      " — pick by length and sequence; order is diagram order</p>";
    var suggestedSet = {};
    (suggestedIds || []).forEach(function (id) { suggestedSet[id] = true; });
    var applySuggested = suggestedIds && suggestedIds.length;
    details.forEach(function (d) {
      var checked = "";
      if (applySuggested) {
        checked = suggestedSet[d.id] ? " checked" : "";
      } else if (details.length === 1) {
        checked = " checked";
      }
      var cls = "chain" + (suggestedSet[d.id] ? " suggested" : "");
      html += '<label class="' + cls + '"><input type="checkbox" value="' +
        esc(d.id) + '"' + checked + "> " + chainLabelInner(d) + "</label>";
    });
    el.innerHTML = html;
  }

  function refreshChainPickers(applySuggestion) {
    var suggestion = null;
    if (state.ref && state.model) {
      suggestion = suggestChainMapping(chainDetails(state.ref), chainDetails(state.model));
    }
    state.suggestion = suggestion;
    var refSuggested = applySuggestion && suggestion ? suggestion.ref_chains : null;
    var modelSuggested = applySuggestion && suggestion ? suggestion.model_chains : null;
    renderChainPicker("ref-chains", state.ref, "ref", refSuggested);
    renderChainPicker("model-chains", state.model, "model", modelSuggested);
    updateChainCountHint();
  }

  function updateChainCountHint() {
    var el = $("chain-map-hint");
    if (!el) return;
    var refDetails = chainDetails(state.ref);
    var modelDetails = chainDetails(state.model);
    var refSel = selectedChains("ref-chains");
    var modelSel = selectedChains("model-chains");
    var suggestion = state.suggestion;
    var paras = [];
    var kind = "";

    if ((state.ref || state.model) && !(state.ref && state.model)) {
      var one = chainDetails(state.ref || state.model);
      if (one.length > 1) {
        paras.push("Upload the other structure to match chains by sequence and length.");
        kind = "warn";
      }
    }

    if (state.ref && state.model && suggestion) {
      paras.push("Suggested mapping: " + suggestion.summary + ". " + suggestion.reason);
    } else if (state.ref && state.model && !suggestion && (refDetails.length > 1 || modelDetails.length > 1)) {
      paras.push(
        "Could not match sequences automatically. Select chains so both sides " +
        "add up to the same number of nucleotides, in the same order."
      );
      kind = kind || "warn";
    }

    if (refSel.length && modelSel.length) {
      var refN = sumLen(refDetails, refSel);
      var modelN = sumLen(modelDetails, modelSel);
      if (refN && modelN && refN !== modelN) {
        paras.push(
          "These lengths do not match: reference " + refSel.join("+") + " is " +
          refN + " nt, model " + modelSel.join("+") + " is " + modelN + " nt."
        );
        var modelAll = sumLen(modelDetails, modelDetails.map(function (d) { return d.id; }));
        var equalRefs = refDetails.filter(function (d) { return d.length === modelN; });
        if (modelSel.length === 1 && equalRefs.length && refSel.length > 1 && refN > modelN) {
          paras.push(
            "This model is one " + modelN + "-nt chain — pick a single reference chain (" +
            equalRefs.map(function (d) { return d.id; }).join(" or ") + "), not both."
          );
        } else if (modelAll && refN === modelAll && modelSel.length < modelDetails.length) {
          paras.push("The full model is " + modelAll + " nt; include the other model chain(s).");
        }
        kind = "err";
      } else if (refN && modelN && refN === modelN) {
        paras.push(
          "Selection: " + refSel.join("+") + " (" + refN + " nt) → " +
          modelSel.join("+") + " (" + modelN + " nt)."
        );
        if (kind !== "err") kind = "ok";
      }
    }

    if (state.ref && state.model && suggestion && refSel.length && modelSel.length) {
      var usingSuggested = sameSet(suggestion.ref_chains, refSel) &&
        sameSet(suggestion.model_chains, modelSel);
      if (!usingSuggested && kind !== "err") {
        paras.push("Current selection differs from the suggestion.");
        if (kind !== "ok") kind = "warn";
      }
    }

    var refAll = (state.ref && state.ref.chains) || [];
    if (refAll.length > 1 && refSel.length && refSel.length < refAll.length) {
      paras.push(
        "Unselected reference chains may still appear in the results if they " +
        "base-pair with your selection; they are not scored."
      );
    }

    if (!paras.length) {
      el.innerHTML = "";
      el.className = "chain-map-hint hidden";
      return;
    }
    el.className = "chain-map-hint" + (kind ? " " + kind : "");
    el.innerHTML = paras.map(function (text) {
      return "<p>" + esc(text) + "</p>";
    }).join("");
  }

  function uploadFile(file, role) {
    var status = $("form-status");
    status.className = "form-status";
    status.textContent = "Uploading " + file.name + "…";
    return fetch("/api/uploads", {
      method: "POST",
      headers: {
        "Content-Type": "application/octet-stream",
        "X-Filename": file.name,
      },
      body: file,
    }).then(function (r) { return r.json().then(function (body) {
      if (!r.ok) throw new Error(body.error || ("upload failed (" + r.status + ")"));
      return body;
    }); }).then(function (info) {
      if (role === "ref") state.ref = info;
      else state.model = info;
      refreshChainPickers(!!(state.ref && state.model));
      status.textContent = "";
    }).catch(function (err) {
      status.className = "form-status err";
      status.textContent = err.message || String(err);
    });
  }

  function onSubmit(ev) {
    ev.preventDefault();
    var status = $("form-status");
    status.className = "form-status";
    if (!state.ref || !state.model) {
      status.className = "form-status err";
      status.textContent = "Upload both a reference and a model first.";
      return;
    }
    var chains = mappingOrder(
      state.suggestion && state.suggestion.ref_chains,
      selectedChains("ref-chains")
    );
    var modelChains = mappingOrder(
      state.suggestion && state.suggestion.model_chains,
      selectedChains("model-chains")
    );
    if (!chains.length || !modelChains.length) {
      status.className = "form-status err";
      status.textContent = "Select at least one chain on each side.";
      return;
    }
    var refN = sumLen(chainDetails(state.ref), chains);
    var modelN = sumLen(chainDetails(state.model), modelChains);
    if (refN && modelN && refN !== modelN) {
      status.className = "form-status err";
      status.textContent = "Selected chains differ in length (" + refN +
        " vs " + modelN + " nt). Pick chains that add up to the same total.";
      return;
    }
    status.textContent = "Starting job…";
    $("generate").disabled = true;
    fetch("/api/jobs", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        ref_upload_id: state.ref.upload_id,
        model_upload_id: state.model.upload_id,
        chains: chains,
        model_chains: modelChains,
        mode: $("layout-mode").value,
        basepairs: $("basepairs").value,
        label: $("label").value,
        notes: $("notes").value,
        force: $("force").checked,
        advanced: {
          pseudoknots: $("adv-pseudoknots").checked,
          rnapuzzler: $("adv-rnapuzzler").checked,
        },
      }),
    }).then(function (r) {
      return r.json().then(function (body) {
        if (!r.ok) throw new Error(body.error || ("create failed (" + r.status + ")"));
        return body;
      });
    }).then(function (body) {
      $("generate").disabled = false;
      var job = (body && body.job) || {};
      state.highlightId = job.id || null;
      var dashText = body.dedup
        ? "Already cached — showing " + (job.label || job.id) +
          ". Check Force re-run to compute it again."
        : "Job queued: " + (job.label || job.id);
      status.className = "form-status ok";
      status.textContent = dashText;
      setDashStatus(dashText, "ok");
      if (window.location.pathname !== "/compare") {
        window.history.replaceState({}, "", "/compare");
      }
      showPanel("dashboard");
      return loadJobs();
    }).catch(function (err) {
      $("generate").disabled = false;
      status.className = "form-status err";
      status.textContent = err.message || String(err);
    });
  }

  function init() {
    var onNew = window.location.pathname.indexOf("/compare/new") === 0;
    showPanel(onNew ? "new" : "dashboard");
    refreshChainPickers(false);
    var form = $("new-form");
    if (form) {
      form.addEventListener("change", function (ev) {
        if (ev.target && ev.target.closest && ev.target.closest(".chains")) {
          updateChainCountHint();
        }
      });
    }
    $("filter").addEventListener("input", renderJobs);
    $("refresh").addEventListener("click", function () { loadJobs(); loadRuntime(); });
    if (window.R2DTTransfer) {
      window.R2DTTransfer.wireImportControls({
        onDone: function () { loadJobs(); loadRuntime(); },
      });
    }
    if (window.R2DTDropzone) {
      window.R2DTDropzone.wire({
        zone: $("ref-drop"),
        input: $("ref-file"),
        nameEl: $("ref-drop-name"),
        onFile: function (file, err) {
          if (err) {
            var status = $("form-status");
            status.className = "form-status err";
            status.textContent = err.message || String(err);
            return;
          }
          if (file) uploadFile(file, "ref");
        },
      });
      window.R2DTDropzone.wire({
        zone: $("model-drop"),
        input: $("model-file"),
        nameEl: $("model-drop-name"),
        onFile: function (file, err) {
          if (err) {
            var status = $("form-status");
            status.className = "form-status err";
            status.textContent = err.message || String(err);
            return;
          }
          if (file) uploadFile(file, "model");
        },
      });
    }
    $("new-form").addEventListener("submit", onSubmit);
    document.querySelectorAll("#tbl thead th[data-sort]").forEach(function (th) {
      th.addEventListener("click", function () {
        var key = th.getAttribute("data-sort");
        if (state.sortKey === key) {
          state.sortDir = state.sortDir === "asc" ? "desc" : "asc";
        } else {
          state.sortKey = key;
          state.sortDir = key === "created" || key.indexOf("inf_") === 0 ? "desc" : "asc";
        }
        renderJobs();
      });
    });
    loadRuntime();
    loadJobs();
  }

  init();
})();
