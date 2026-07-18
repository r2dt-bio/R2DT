/* Shared export / import helpers for workstation dashboards. */
(function (global) {
  function exportJob(jobId) {
    if (!jobId) return;
    global.location.href = "/api/jobs/" + encodeURIComponent(jobId) + "/export";
  }

  function addExportButton(actionsEl, job) {
    if (!actionsEl || !job || job.status !== "ready") return;
    var btn = document.createElement("button");
    btn.type = "button";
    btn.textContent = "Export";
    btn.title = "Download .r2dt-job.zip";
    btn.addEventListener("click", function () { exportJob(job.id); });
    actionsEl.appendChild(btn);
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
    addExportButton: addExportButton,
    wireImportControls: wireImportControls,
  };
})(window);
