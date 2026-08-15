# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# pylint: skip-file

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "R2DT"
copyright = "2026"
author = "R2DT Team"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

import os

import sphinx_rtd_theme

extensions = [
    "sphinx_rtd_theme",
    "myst_parser",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

html_theme_options = {
    # Toc options
    "collapse_navigation": False,
    "logo_only": True,
}

# -- Options for MyST parser -------------------------------------------------
myst_heading_anchors = 3
myst_enable_extensions = [
    "colon_fence",  # enable admonitions
    "attrs_inline",  # enable inline attributes for images
]

# -- Read the Docs Canonical URL -----------------------------------------------

# Set canonical URL from the Read the Docs Domain
html_baseurl = os.environ.get("READTHEDOCS_CANONICAL_URL", "")

# Tell Jinja2 templates the build is running on Read the Docs
if os.environ.get("READTHEDOCS", "") == "True":
    if "html_context" not in globals():
        html_context = {}
    html_context["READTHEDOCS"] = True

# -- Logo and favicon -------------------------------------------------
html_logo = "images/r2dt-banner-logo.png"

# -- Workstation starter-pack download ----------------------------------------
# docs/files/r2dt-workstation-start.zip is generated here, at build time, from
# scripts/workstation/ -- it is deliberately not committed, so it can never go
# stale. workstation.md links it with {download}. `just workstation-pack`
# builds the same zip locally.


def _pack_workstation_zip(app):
    import zipfile
    from pathlib import Path

    docs_dir = Path(app.srcdir).resolve()
    src_dir = docs_dir.parent / "scripts" / "workstation"
    out = docs_dir / "files" / "r2dt-workstation-start.zip"
    out.parent.mkdir(parents=True, exist_ok=True)
    members = [
        "README.txt",
        "Start-macOS.command",
        "Start-Windows.bat",
        "Start-Linux.sh",
    ]
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as zf:
        for name in members:
            src = src_dir / name
            info = zipfile.ZipInfo(f"r2dt-workstation-start/{name}")
            # Fixed timestamp and mode keep the zip reproducible; launchers
            # need the executable bit.
            info.date_time = (2026, 1, 1, 0, 0, 0)
            mode = 0o755 if name != "README.txt" else 0o644
            info.external_attr = mode << 16
            info.compress_type = zipfile.ZIP_DEFLATED
            zf.writestr(info, src.read_bytes())


def setup(app):
    app.connect("builder-inited", _pack_workstation_zip)
