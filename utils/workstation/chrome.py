"""Shared workstation chrome (header) and mode registry."""

from __future__ import annotations

import html as html_lib
from typing import Optional

# Short URL path segment → job mode id + UI copy.
MODES = (
    {
        "id": "draw",
        "path": "2d",
        "nav": "2D",
        "title": "2D",
        "blurb": "Draw secondary structure diagrams from FASTA sequences.",
        "cta": "New 2D job",
        "ready": True,
    },
    {
        "id": "pdb",
        "path": "pdb",
        "nav": "2D + 3D",
        "title": "2D + 3D",
        "blurb": "Linked 2D topology and 3D structure for one PDB or mmCIF.",
        "cta": "New 2D + 3D job",
        "ready": True,
    },
    {
        "id": "compare",
        "path": "compare",
        "nav": "2D + 2D + 3D",
        "title": "2D + 2D + 3D",
        "blurb": "CASP-style reference vs model comparison with INF scores and editing.",
        "cta": "New comparison",
        "ready": True,
    },
    {
        "id": "align",
        "path": "align",
        "nav": "Alignments",
        "title": "Alignments",
        "blurb": "Stockholm / R-scape alignments with covariation-annotated diagrams.",
        "cta": "New alignment job",
        "ready": True,
    },
)

MODE_BY_PATH = {m["path"]: m for m in MODES}
MODE_BY_ID = {m["id"]: m for m in MODES}

VALID_JOB_MODES = frozenset(m["id"] for m in MODES)


def normalize_job_mode(mode: Optional[str]) -> str:
    """Return a valid job mode; legacy jobs default to compare."""
    if mode in VALID_JOB_MODES:
        return str(mode)
    return "compare"


def export_menu_html(
    job_id: str,
    *,
    include_html: bool = True,
    variant: str = "chrome",
) -> str:
    """Dropdown: R2DT work package and optional shareable HTML viewer zip."""
    safe_id = html_lib.escape(job_id, quote=True)
    package_href = f"/api/jobs/{safe_id}/export"
    html_href = f"/api/jobs/{safe_id}/export/html"
    if variant == "inf":
        summary_cls = "ws-inf-export"
        menu_cls = "ws-export-menu ws-export-menu--inf"
    else:
        summary_cls = "ws-chrome-export"
        menu_cls = "ws-export-menu"
    html_item = ""
    if include_html:
        html_item = (
            f'<a href="{html_href}" '
            f'title="Static viewer with edits baked in">'
            f"Shareable HTML</a>"
        )
    return (
        f'<details class="{menu_cls}">'
        f'<summary class="{summary_cls}" title="Export this job">Export</summary>'
        f'<div class="ws-export-menu-panel" role="menu">'
        f'<a href="{package_href}" '
        f'title="Download .r2dt-job.zip for another workstation">'
        f"R2DT work package</a>"
        f"{html_item}"
        f"</div></details>"
    )


def chrome_header(
    active_path: Optional[str] = None,
    *,
    job_label: str = "",
    job_id: str = "",
    show_export: bool = False,
    export_html: bool = True,
) -> str:
    """Return the sticky header HTML used on every workstation page."""
    brand = (
        '<a class="ws-chrome-brand" href="/" title="R2DT workstation home">'
        '<img class="ws-chrome-logo" src="/static/r2dt-logo-blue.svg" alt="">'
        "<span>R2DT workstation</span>"
        "</a>"
    )
    links = []
    for mode in MODES:
        path = mode["path"]
        href = f"/{path}"
        cls = "ws-chrome-link"
        if active_path == path:
            cls += " is-active"
        label = html_lib.escape(mode["nav"])
        links.append(f'<a class="{cls}" href="{href}">{label}</a>')
    nav = '<nav class="ws-chrome-nav" aria-label="Modes">' + "".join(links) + "</nav>"
    ext = (
        '<nav class="ws-chrome-ext" aria-label="R2DT links">'
        '<a href="https://docs.r2dt.bio" target="_blank" rel="noopener">Docs</a>'
        '<a href="https://r2dt.bio" target="_blank" rel="noopener">r2dt.bio</a>'
        '<a href="https://github.com/r2dt-bio/r2dt" target="_blank" rel="noopener">'
        "GitHub</a>"
        "</nav>"
    )
    job = ""
    if job_label or job_id:
        display = (job_label or "").strip() or job_id
        safe_label = html_lib.escape(display, quote=True)
        safe_id = html_lib.escape(job_id or display, quote=True)
        export = ""
        if show_export and job_id:
            export = export_menu_html(job_id, include_html=export_html)
        job = (
            f'<span class="ws-chrome-job-wrap">'
            f'<span class="ws-chrome-job" title="{safe_id}">{safe_label}</span>'
            f"{export}"
            f"</span>"
        )
    return f'<header class="ws-chrome">{brand}{nav}{ext}{job}</header>\n'


FAVICON_LINKS = (
    '<link rel="icon" type="image/png" href="/static/favicon/favicon-96x96.png" '
    'sizes="96x96">\n'
    '<link rel="icon" type="image/svg+xml" href="/static/favicon/favicon.svg">\n'
    '<link rel="shortcut icon" href="/static/favicon/favicon.ico">\n'
    '<link rel="apple-touch-icon" sizes="180x180" '
    'href="/static/favicon/apple-touch-icon.png">\n'
)
