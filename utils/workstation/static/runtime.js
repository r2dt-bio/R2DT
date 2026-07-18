(function () {
  function friendlyWorkspace(path) {
    if (!path) return "";
    if (path === "/workspace" || path.indexOf("/workspace") === 0) {
      return "~/.r2dt-workstation (mounted at /workspace in Docker)";
    }
    return path;
  }

  function loadRuntime() {
    var el = document.getElementById("runtime");
    if (!el) return;
    return fetch("/api/runtime").then(function (r) { return r.json(); }).then(function (data) {
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

  loadRuntime();
})();
