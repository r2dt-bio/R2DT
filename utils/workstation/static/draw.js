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
    return fetch("/api/jobs?mode=draw").then(function (r) { return r.json(); }).then(function (data) {
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

  function rowValue(j, key, seqMap) {
    var params = j.params || {};
    switch (key) {
      case "seq": return seqMap[j.id] || 0;
      case "label": return (j.label || j.id || "").toLowerCase();
      case "seq_id": return String(params.seq_id || "").toLowerCase();
      case "length": return params.length == null ? -1 : +params.length;
      case "layout": return String(params.layout || "");
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
      var hay = [j.label, j.notes, j.id, (j.params || {}).seq_id].join(" ").toLowerCase();
      return hay.indexOf(q) !== -1;
    });
    var seqMap = seqNumbers(state.jobs);
    var rows = sortedRows(filtered, seqMap);
    $("count").textContent = rows.length + " job" + (rows.length === 1 ? "" : "s");
    updateSortHeaders();
    var tbody = $("rows");
    tbody.innerHTML = "";
    if (!rows.length) {
      tbody.innerHTML = '<tr><td colspan="8">No 2D jobs yet. Use New 2D job to start one.</td></tr>';
      return;
    }
    rows.forEach(function (j) {
      var params = j.params || {};
      var tr = document.createElement("tr");
      tr.innerHTML =
        '<td class="num">' + (seqMap[j.id] || "—") + "</td>" +
        labelCell(j) +
        "<td>" + esc(params.seq_id || "") + "</td>" +
        '<td class="num">' + (params.length != null ? params.length : "—") + "</td>" +
        "<td>" + esc(params.layout || "") + "</td>" +
        '<td title="' + esc(j.created || "") + '">' + esc(fmtDate(j.created)) + "</td>" +
        '<td><span class="badge ' + esc(j.status || "") + '">' + esc(j.status || "") + "</span></td>" +
        '<td class="actions-col"><div class="row-actions"></div></td>';
      var actions = tr.querySelector(".row-actions");
      if (j.status === "ready" && j.results_url) {
        var open = document.createElement("a");
        open.href = j.results_url;
        open.textContent = "Open";
        open.target = "_blank";
        open.rel = "noopener";
        actions.appendChild(open);
      }
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
    var fasta = ($("fasta").value || "").trim();
    if (!fasta) {
      status.className = "form-status err";
      status.textContent = "Paste or upload a FASTA first.";
      return;
    }
    status.textContent = "Starting job(s)…";
    $("generate").disabled = true;
    fetch("/api/jobs", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        job_mode: "draw",
        fasta: fasta,
        layout: $("layout").value,
        label: $("label").value,
        notes: $("notes").value,
        force: $("force").checked,
      }),
    }).then(function (r) {
      return r.json().then(function (body) {
        if (!r.ok) throw new Error(body.error || ("create failed (" + r.status + ")"));
        return body;
      });
    }).then(function (body) {
      $("generate").disabled = false;
      status.className = "form-status ok";
      status.textContent = body.message || ("Queued " + (body.jobs || []).length + " job(s)");
      if (window.location.pathname !== "/2d") {
        window.history.replaceState({}, "", "/2d");
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
    var list = (window.R2DT_WS_EXAMPLES && window.R2DT_WS_EXAMPLES.draw) || [];
    if (!host || !list.length) return;
    list.forEach(function (ex) {
      var btn = document.createElement("button");
      btn.type = "button";
      btn.className = "ws-example";
      btn.textContent = ex.label;
      if (ex.note) btn.title = ex.note;
      btn.addEventListener("click", function () {
        $("fasta").value = ex.fasta || "";
        if (ex.layout && $("layout")) $("layout").value = ex.layout;
        if (!$("label").value) $("label").value = ex.label;
        host.querySelectorAll(".ws-example").forEach(function (b) {
          b.classList.toggle("is-active", b === btn);
        });
        $("form-status").className = "form-status";
        $("form-status").textContent = ex.note
          ? "Loaded " + ex.label + " (" + ex.note + ")"
          : "Loaded " + ex.label;
      });
      host.appendChild(btn);
    });
  }

  function init() {
    var onNew = window.location.pathname.indexOf("/2d/new") === 0;
    showPanel(onNew ? "new" : "dashboard");
    initExamples();
    $("filter").addEventListener("input", renderJobs);
    $("refresh").addEventListener("click", function () { loadJobs(); loadRuntime(); });
    if (window.R2DTTransfer) {
      window.R2DTTransfer.wireImportControls({
        onDone: function () { loadJobs(); loadRuntime(); },
      });
    }
    $("fasta-file").addEventListener("change", function (ev) {
      var f = ev.target.files && ev.target.files[0];
      if (!f) return;
      var reader = new FileReader();
      reader.onload = function () {
        $("fasta").value = String(reader.result || "");
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
          state.sortDir = key === "created" || key === "length" ? "desc" : "asc";
        }
        renderJobs();
      });
    });
    loadRuntime();
    loadJobs();
  }

  init();
})();
