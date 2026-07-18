(function () {
  var state = {
    jobs: [],
    pollTimer: null,
    sortKey: "created",
    sortDir: "desc",
  };

  function $(id) { return document.getElementById(id); }

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

  function loadJobs() {
    return fetch("/api/jobs?mode=align").then(function (r) { return r.json(); }).then(function (data) {
      state.jobs = data.jobs || [];
      renderJobs();
      maybePoll();
    });
  }

  function seqNumbers(jobs) {
    var ordered = jobs.slice().sort(function (a, b) {
      return String(a.created || "").localeCompare(String(b.created || "")) ||
        String(a.id || "").localeCompare(String(b.id || ""));
    });
    var map = {};
    ordered.forEach(function (j, i) { map[j.id] = i + 1; });
    return map;
  }

  function rscapeLabel(params) {
    if (params.ran_rscape) return "ran";
    if (params.will_run_rscape) return "will run";
    if (params.has_covariation) return "present";
    return "—";
  }

  function rowValue(j, key, seqMap) {
    var params = j.params || {};
    switch (key) {
      case "seq": return seqMap[j.id] || 0;
      case "label": return (j.label || j.id || "").toLowerCase();
      case "rscape": return rscapeLabel(params);
      case "stitch": return params.stitch ? 1 : 0;
      case "created": return j.created || "";
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

  function esc(s) {
    return String(s == null ? "" : s)
      .replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;");
  }

  function labelCell(j) {
    var text = esc(j.label || j.id);
    if (j.status === "ready" && j.results_url) {
      return (
        '<td><a class="ws-label-link" href="' + esc(j.results_url) +
        '" target="_blank" rel="noopener">' + text + "</a></td>"
      );
    }
    return "<td>" + text + "</td>";
  }

  function renderJobs() {
    var q = ($("filter").value || "").trim().toLowerCase();
    var filtered = state.jobs.filter(function (j) {
      if (!q) return true;
      var hay = [j.label, j.notes, j.id].join(" ").toLowerCase();
      return hay.indexOf(q) !== -1;
    });
    var seqMap = seqNumbers(state.jobs);
    var rows = sortedRows(filtered, seqMap);
    $("count").textContent = rows.length + " job" + (rows.length === 1 ? "" : "s");
    updateSortHeaders();
    var tbody = $("rows");
    tbody.innerHTML = "";
    if (!rows.length) {
      tbody.innerHTML =
        '<tr><td colspan="7">No alignment jobs yet. Use New alignment job to start one.</td></tr>';
      return;
    }
    rows.forEach(function (j) {
      var params = j.params || {};
      var tr = document.createElement("tr");
      tr.innerHTML =
        '<td class="num">' + (seqMap[j.id] || "—") + "</td>" +
        labelCell(j) +
        '<td><span class="badge ' + esc(j.status || "") + '">' + esc(j.status || "") + "</span></td>" +
        '<td class="actions-col"><div class="row-actions"></div></td>' +
        "<td>" + esc(rscapeLabel(params)) + "</td>" +
        "<td>" + (params.stitch ? "yes" : "no") + "</td>" +
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
      tbody.appendChild(tr);
    });
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

  function onSubmit(ev) {
    ev.preventDefault();
    var status = $("form-status");
    status.className = "form-status";
    var stockholm = ($("stockholm").value || "").trim();
    if (!stockholm) {
      status.className = "form-status err";
      status.textContent = "Paste or upload a Stockholm alignment first.";
      return;
    }
    status.textContent = "Starting job…";
    $("generate").disabled = true;
    fetch("/api/jobs", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        job_mode: "align",
        stockholm: stockholm,
        stitch: $("stitch").checked,
        label: $("label").value,
        notes: $("notes").value,
        force: $("force").checked,
        advanced: {
          include_novel: $("adv-include-novel").checked,
          all_structures: $("adv-all-structures").checked,
          auto_repair: $("adv-auto-repair").checked,
          monochrome: $("adv-monochrome").checked,
          max_unpaired: parseInt($("adv-max-unpaired").value, 10) || 0,
          color_by: $("adv-color-by").value || "none",
        },
      }),
    }).then(function (r) {
      return r.json().then(function (body) {
        if (!r.ok) throw new Error(body.error || ("create failed (" + r.status + ")"));
        return body;
      });
    }).then(function (body) {
      $("generate").disabled = false;
      status.className = "form-status ok";
      var job = body.job || {};
      status.textContent = body.dedup
        ? "Identical job already cached: " + (job.label || job.id)
        : "Queued " + (job.label || job.id);
      if (window.location.pathname !== "/align") {
        window.history.replaceState({}, "", "/align");
      }
      showPanel("dashboard");
      return loadJobs();
    }).catch(function (err) {
      $("generate").disabled = false;
      status.className = "form-status err";
      status.textContent = err.message || String(err);
    });
  }

  function initExamples() {
    var host = $("examples");
    var list = (window.R2DT_WS_EXAMPLES && window.R2DT_WS_EXAMPLES.align) || [];
    if (!host || !list.length) return;
    list.forEach(function (ex) {
      var btn = document.createElement("button");
      btn.type = "button";
      btn.className = "ws-example";
      btn.textContent = ex.label;
      if (ex.note) btn.title = ex.note;
      btn.addEventListener("click", function () {
        statusLoad(ex);
        host.querySelectorAll(".ws-example").forEach(function (b) {
          b.classList.toggle("is-active", b === btn);
        });
      });
      host.appendChild(btn);
    });
  }

  function statusLoad(ex) {
    var status = $("form-status");
    status.className = "form-status";
    status.textContent = "Loading " + ex.label + "…";
    fetch("/api/examples/stockholm/" + encodeURIComponent(ex.id))
      .then(function (r) {
        return r.json().then(function (body) {
          if (!r.ok) throw new Error(body.error || "load failed");
          return body;
        });
      })
      .then(function (body) {
        $("stockholm").value = body.stockholm || "";
        if (!$("label").value) $("label").value = ex.label;
        status.textContent = ex.note
          ? "Loaded " + ex.label + " (" + ex.note + ")"
          : "Loaded " + ex.label;
      })
      .catch(function (err) {
        status.className = "form-status err";
        status.textContent = err.message || String(err);
      });
  }

  function init() {
    var onNew = window.location.pathname.indexOf("/align/new") === 0;
    showPanel(onNew ? "new" : "dashboard");
    initExamples();
    $("filter").addEventListener("input", renderJobs);
    $("refresh").addEventListener("click", function () { loadJobs(); loadRuntime(); });
    if (window.R2DTTransfer) {
      window.R2DTTransfer.wireImportControls({
        onDone: function () { loadJobs(); loadRuntime(); },
      });
    }
    $("stk-file").addEventListener("change", function (ev) {
      var f = ev.target.files && ev.target.files[0];
      if (!f) return;
      var reader = new FileReader();
      reader.onload = function () {
        $("stockholm").value = String(reader.result || "");
        if (!$("label").value) {
          $("label").value = f.name.replace(/\.(stk|sto|stockholm|txt)$/i, "");
        }
      };
      reader.readAsText(f);
    });
    $("new-form").addEventListener("submit", onSubmit);
    document.querySelectorAll("#tbl thead th[data-sort]").forEach(function (th) {
      th.addEventListener("click", function () {
        var key = th.getAttribute("data-sort");
        if (state.sortKey === key) {
          state.sortDir = state.sortDir === "asc" ? "desc" : "asc";
        } else {
          state.sortKey = key;
          state.sortDir = key === "created" ? "desc" : "asc";
        }
        renderJobs();
      });
    });
    loadRuntime();
    loadJobs();
  }

  init();
})();
