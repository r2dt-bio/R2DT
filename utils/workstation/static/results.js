(function () {
  function activate(tabId) {
    var stage = document.querySelector(".ws-svg-stage");
    if (!stage) return;
    stage.querySelectorAll(".ws-svg-tab").forEach(function (btn) {
      var on = btn.getAttribute("data-tab") === tabId;
      btn.classList.toggle("is-active", on);
      btn.setAttribute("aria-selected", on ? "true" : "false");
    });
    var activeSrc = null;
    stage.querySelectorAll(".ws-svg-panel").forEach(function (panel) {
      var on = panel.getAttribute("data-panel") === tabId;
      panel.classList.toggle("is-active", on);
      if (on) {
        panel.removeAttribute("hidden");
        var img = panel.querySelector(".ws-svg-preview");
        if (img) activeSrc = img.getAttribute("src");
      } else {
        panel.setAttribute("hidden", "");
      }
    });
    var dl = document.getElementById("ws-svg-download");
    if (dl && activeSrc) dl.setAttribute("href", activeSrc);
  }

  function init() {
    var stage = document.querySelector(".ws-svg-stage");
    if (!stage) return;
    stage.querySelectorAll(".ws-svg-tab").forEach(function (btn) {
      btn.addEventListener("click", function () {
        activate(btn.getAttribute("data-tab"));
      });
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
