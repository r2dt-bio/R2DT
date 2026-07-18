(function () {
  var state = {
    jobs: [],
    ref: null,
    model: null,
    pollTimer: null,
  };

  function $(id) { return document.getElementById(id); }

  function fmt(n) {
    if (n === null || n === undefined || n === "") return "—";
    if (typeof n === "number") return n.toFixed(3);
    return String(n);
  }

  function setTabs() {
    document.querySelectorAll(".tab").forEach(function (btn) {
      btn.addEventListener("click", function () {
        document.querySelectorAll(".tab").forEach(function (b) { b.classList.remove("active"); });
        document.querySelectorAll(".panel").forEach(function (p) { p.classList.remove("active"); });
        btn.classList.add("active");
        $("panel-" + btn.getAttribute("data-panel")).classList.add("active");
      });
    });
  }

  function loadRuntime() {
    return fetch("/api/runtime").then(function (r) { return r.json(); }).then(function (data) {
      var el = $("runtime");
      if (!data.ok) {
        el.textContent = "Runtime error: " + data.message;
        el.style.color = "var(--danger)";
      } else {
        el.textContent = data.message + " · workspace " + data.workspace +
          (data.current_job ? " · running " + data.current_job : "");
      }
    });
  }

  function loadJobs() {
    return fetch("/api/jobs").then(function (r) { return r.json(); }).then(function (data) {
      state.jobs = data.jobs || [];
      renderJobs();
      maybePoll();
    });
  }

  function renderJobs() {
    var q = ($("filter").value || "").trim().toLowerCase();
    var rows = state.jobs.filter(function (j) {
      if (!q) return true;
      var hay = [j.label, j.notes,
        (j.inputs || {}).ref_name, (j.inputs || {}).model_name,
        (j.params || {}).chains].join(" ").toLowerCase();
      return hay.indexOf(q) !== -1;
    });
    $("count").textContent = rows.length + " job" + (rows.length === 1 ? "" : "s");
    var tbody = $("rows");
    tbody.innerHTML = "";
    if (!rows.length) {
      tbody.innerHTML = '<tr><td colspan="11">No comparisons yet. Start one under New comparison.</td></tr>';
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
        "<td>" + esc(j.label || j.id) + "</td>" +
        "<td>" + esc(inputs.ref_name || "") + "</td>" +
        "<td>" + esc(inputs.model_name || "") + "</td>" +
        "<td>" + esc(params.chains || "") +
          (params.model_chains ? " → " + esc(params.model_chains) : "") + "</td>" +
        "<td>" + esc(j.created || "") + "</td>" +
        '<td class="num">' + fmt(inf.wc) + "</td>" +
        '<td class="num">' + fmt(inf.nwc) + "</td>" +
        '<td class="num">' + fmt(inf.all) + "</td>" +
        "<td>" + bp + "</td>" +
        '<td><span class="badge ' + esc(j.status || "") + '">' + esc(j.status || "") + "</span></td>" +
        '<td class="row-actions"></td>';
      var actions = tr.querySelector(".row-actions");
      if (j.status === "ready" && j.viewer_url) {
        var open = document.createElement("a");
        open.href = j.viewer_url;
        open.textContent = "Open";
        open.target = "_blank";
        open.rel = "noopener";
        actions.appendChild(open);
      }
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

  function esc(s) {
    return String(s == null ? "" : s)
      .replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;");
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

  function renderChainPicker(containerId, info, role) {
    var el = $(containerId);
    if (!info) {
      el.innerHTML = "<p class=\"hint\">Upload a structure to detect RNA chains.</p>";
      return;
    }
    var chains = info.chains || [];
    var html = "";
    if (role === "ref" && !info.compare_ready) {
      html += '<p class="err">Compare needs an mmCIF reference (.cif). This file is ' +
        esc(info.format) + ".</p>";
    }
    if (!chains.length) {
      html += '<p class="err">No RNA chains detected.</p>';
      el.innerHTML = html;
      return;
    }
    html += '<p class="hint">' + esc(info.filename) + " · " + chains.length +
      " RNA chain" + (chains.length === 1 ? "" : "s") +
      (info.compare_ready ? "" : "") + "</p>";
    if (chains.length > 1) {
      html += '<p class="warn">Select the chain(s) to include (order = diagram order).</p>';
    }
    chains.forEach(function (c) {
      var checked = chains.length === 1 ? " checked" : "";
      html += '<label class="chain"><input type="checkbox" value="' + esc(c) + '"' +
        checked + "> " + esc(c) + "</label>";
    });
    el.innerHTML = html;
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
      if (role === "ref") {
        state.ref = info;
        renderChainPicker("ref-chains", info, "ref");
      } else {
        state.model = info;
        renderChainPicker("model-chains", info, "model");
      }
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
    var chains = selectedChains("ref-chains");
    var modelChains = selectedChains("model-chains");
    if (!chains.length || !modelChains.length) {
      status.className = "form-status err";
      status.textContent = "Select at least one chain on each side.";
      return;
    }
    if (chains.length !== modelChains.length) {
      status.className = "form-status err";
      status.textContent = "Select the same number of reference and model chains.";
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
        mode: $("mode").value,
        basepairs: $("basepairs").value,
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
      if (body.dedup) {
        status.className = "form-status ok";
        status.textContent = "Already cached — opening existing job.";
      } else {
        status.className = "form-status ok";
        status.textContent = "Job queued: " + body.job.id;
      }
      document.querySelector('.tab[data-panel="dashboard"]').click();
      return loadJobs();
    }).catch(function (err) {
      $("generate").disabled = false;
      status.className = "form-status err";
      status.textContent = err.message || String(err);
    });
  }

  function init() {
    setTabs();
    renderChainPicker("ref-chains", null, "ref");
    renderChainPicker("model-chains", null, "model");
    $("filter").addEventListener("input", renderJobs);
    $("refresh").addEventListener("click", function () { loadJobs(); loadRuntime(); });
    $("ref-file").addEventListener("change", function (ev) {
      var f = ev.target.files && ev.target.files[0];
      if (f) uploadFile(f, "ref");
    });
    $("model-file").addEventListener("change", function (ev) {
      var f = ev.target.files && ev.target.files[0];
      if (f) uploadFile(f, "model");
    });
    $("new-form").addEventListener("submit", onSubmit);
    loadRuntime();
    loadJobs();
  }

  init();
})();
