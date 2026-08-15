/* Shared drag-and-drop upload zone for structure files. */
(function (global) {
  var STRUCTURE_RE = /\.(cif|pdb)(\.gz)?$/i;

  function isStructureFile(file) {
    if (!file || !file.name) return false;
    return STRUCTURE_RE.test(file.name);
  }

  function setDragging(zone, on) {
    zone.classList.toggle("is-dragging", !!on);
  }

  function setHasFile(zone, nameEl, filename) {
    zone.classList.toggle("has-file", !!filename);
    if (nameEl) {
      if (filename) {
        nameEl.textContent = filename;
        nameEl.classList.remove("hidden");
      } else {
        nameEl.textContent = "";
        nameEl.classList.add("hidden");
      }
    }
  }

  /**
   * Wire a drop zone.
   * @param {object} opts
   * @param {HTMLElement} opts.zone
   * @param {HTMLInputElement} opts.input
   * @param {HTMLElement} [opts.nameEl]
   * @param {function(File): void} opts.onFile
   */
  function wire(opts) {
    var zone = opts.zone;
    var input = opts.input;
    var nameEl = opts.nameEl || null;
    var onFile = opts.onFile;
    if (!zone || !input || typeof onFile !== "function") return;

    function takeFile(file) {
      if (!file) return;
      if (!isStructureFile(file)) {
        onFile(null, new Error("Expected a .cif or .pdb file (optionally .gz)"));
        return;
      }
      setHasFile(zone, nameEl, file.name);
      onFile(file, null);
    }

    zone.addEventListener("click", function (ev) {
      if (ev.target === input) return;
      input.click();
    });
    zone.addEventListener("keydown", function (ev) {
      if (ev.key === "Enter" || ev.key === " ") {
        ev.preventDefault();
        input.click();
      }
    });
    input.addEventListener("change", function () {
      var f = input.files && input.files[0];
      if (f) takeFile(f);
    });
    zone.addEventListener("dragenter", function (ev) {
      ev.preventDefault();
      setDragging(zone, true);
    });
    zone.addEventListener("dragover", function (ev) {
      ev.preventDefault();
      if (ev.dataTransfer) ev.dataTransfer.dropEffect = "copy";
      setDragging(zone, true);
    });
    zone.addEventListener("dragleave", function (ev) {
      if (ev.target === zone || !zone.contains(ev.relatedTarget)) {
        setDragging(zone, false);
      }
    });
    zone.addEventListener("drop", function (ev) {
      ev.preventDefault();
      setDragging(zone, false);
      var file = ev.dataTransfer && ev.dataTransfer.files && ev.dataTransfer.files[0];
      takeFile(file);
    });

    return {
      setFilename: function (name) {
        setHasFile(zone, nameEl, name || "");
      },
      clear: function () {
        input.value = "";
        setHasFile(zone, nameEl, "");
      },
    };
  }

  global.R2DTDropzone = {
    wire: wire,
    isStructureFile: isStructureFile,
  };
})(window);
