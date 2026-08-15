/* Shared export / import helpers for workstation dashboards. */
(function (global) {
  function exportJob(jobId) {
    if (!jobId) return;
    global.location.href = "/api/jobs/" + encodeURIComponent(jobId) + "/export";
  }

  function exportJobHtml(jobId) {
    if (!jobId) return;
    global.location.href =
      "/api/jobs/" + encodeURIComponent(jobId) + "/export/html";
  }

  function addExportButton(actionsEl, job) {
    if (!actionsEl || !job || job.status !== "ready") return;
    var includeHtml = Boolean(job.viewer_url);
    var wrap = document.createElement("details");
    wrap.className = "ws-export-menu ws-export-menu--row";
    var summary = document.createElement("summary");
    summary.textContent = "Export";
    summary.title = "Export this job";
    wrap.appendChild(summary);
    var panel = document.createElement("div");
    panel.className = "ws-export-menu-panel";
    panel.setAttribute("role", "menu");

    var pkg = document.createElement("button");
    pkg.type = "button";
    pkg.textContent = "R2DT work package";
    pkg.title = "Download .r2dt-job.zip for another workstation";
    pkg.addEventListener("click", function () {
      wrap.open = false;
      exportJob(job.id);
    });
    panel.appendChild(pkg);

    if (includeHtml) {
      var htmlBtn = document.createElement("button");
      htmlBtn.type = "button";
      htmlBtn.textContent = "Shareable HTML";
      htmlBtn.title = "Static viewer with edits baked in";
      htmlBtn.addEventListener("click", function () {
        wrap.open = false;
        exportJobHtml(job.id);
      });
      panel.appendChild(htmlBtn);
    }

    wrap.appendChild(panel);
    actionsEl.appendChild(wrap);
  }

  function wireImportControls(opts) {
    var input = document.getElementById(opts.inputId || "import-file");
    var button = document.getElementById(opts.buttonId || "import-btn");
    var statusEl = document.getElementById(opts.statusId || "import-status");
    if (!input || !button) return;

    function setStatus(text, cls) {
      if (!statusEl) return;
      statusEl.className = "form-status" + (cls ? " " + cls : "");
      statusEl.textContent = text || "";
    }

    button.addEventListener("click", function () { input.click(); });
    input.addEventListener("change", function () {
      var file = input.files && input.files[0];
      input.value = "";
      if (!file) return;
      setStatus("Importing " + file.name + "…");
      button.disabled = true;
      file.arrayBuffer().then(function (buf) {
        return fetch("/api/jobs/import", {
          method: "POST",
          headers: {
            "Content-Type": "application/zip",
            "X-Filename": file.name,
          },
          body: buf,
        }).then(function (r) {
          return r.json().then(function (body) {
            if (!r.ok) throw new Error(body.error || ("import failed (" + r.status + ")"));
            return body;
          });
        });
      }).then(function (body) {
        button.disabled = false;
        var job = body.job || {};
        var mode = job.mode || "";
        var path = ({ draw: "/2d", pdb: "/pdb", compare: "/compare", align: "/align" })[mode] || "/";
        setStatus(
          "Imported " + (job.label || job.id) + " — see " + path,
          "ok"
        );
        if (typeof opts.onDone === "function") opts.onDone(body);
      }).catch(function (err) {
        button.disabled = false;
        setStatus(err.message || String(err), "err");
      });
    });
  }

  global.R2DTTransfer = {
    exportJob: exportJob,
    exportJobHtml: exportJobHtml,
    addExportButton: addExportButton,
    wireImportControls: wireImportControls,
  };
})(window);
