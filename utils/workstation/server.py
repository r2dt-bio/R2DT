"""Stdlib HTTP server for the local R2DT workstation."""

from __future__ import annotations

import html as html_lib
import json
import mimetypes
import re
import uuid
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Dict, Optional, Tuple
from urllib.parse import parse_qs, urlparse

from rich import print as rprint

from utils.workstation import edits as edits_mod
from utils.workstation.align import (
    EXAMPLE_STOCKHOLM,
    create_align_job_from_stockholm,
    list_align_result_svgs,
    svg_gallery_entry,
)
from utils.workstation.catalog import Catalog
from utils.workstation.chains import is_structure_filename, list_rna_chains
from utils.workstation.chrome import (
    FAVICON_LINKS,
    MODE_BY_PATH,
    MODES,
    chrome_header,
    export_menu_html,
    normalize_job_mode,
)
from utils.workstation.jobs import (
    JobRunner,
    create_draw_jobs_from_fasta,
    create_job_from_uploads,
    create_pdb_job_from_upload,
    find_draw_svg_variants,
    require_runtime,
)
from utils.workstation.package import export_job_bytes, import_job_bytes
from utils.workstation.publish import export_shareable_html_bytes
from utils.workstation.security import local_mutating_request_error, path_is_within

STATIC_DIR = Path(__file__).resolve().parent / "static"

_HEAD_LINKS = (
    FAVICON_LINKS
    + '<link rel="stylesheet" href="/static/chrome.css">\n'
    + '<link rel="stylesheet" href="/static/style.css">\n'
)


def _fill_page(template: str, *, active_path: Optional[str] = None, **repl) -> str:
    """Substitute <!--HEADER--> / <!--HEAD--> and optional string placeholders."""
    html = template.replace("<!--HEAD-->", _HEAD_LINKS)
    html = html.replace("<!--HEADER-->", chrome_header(active_path=active_path))
    for key, value in repl.items():
        html = html.replace(f"<!--{key}-->", str(value))
    return html


def _mode_cards_html() -> str:
    cards = []
    for mode in MODES:
        path = mode["path"]
        title = mode["title"]
        blurb = mode["blurb"]
        if mode["ready"]:
            href = f"/{path}"
            cta = mode["cta"]
            extra = f'<a class="ws-card-cta" href="{href}/new">{cta}</a>'
            dash = f'<a class="ws-card-dash" href="{href}">View jobs →</a>'
        else:
            href = f"/{path}"
            extra = '<span class="ws-card-soon">Coming soon</span>'
            dash = f'<a class="ws-card-dash" href="{href}">Open mode →</a>'
        cards.append(
            f'<article class="ws-mode-card">'
            f'<h2><a href="{href}">{title}</a></h2>'
            f"<p>{blurb}</p>"
            f'<div class="ws-card-actions">{extra}{dash}</div>'
            f"</article>"
        )
    return "\n".join(cards)


_VIEWER_FIT_SCRIPT = """
<script>
(function () {
  function fitWorkstationViewer() {
    var compareRoot = document.querySelector('.r2dt-compare-root');
    var viewerRoot = document.querySelector('.r2dt-viewer-root');
    var root = compareRoot || viewerRoot;
    if (!root) return;
    var chrome = document.querySelector('.ws-chrome');
    var h1 = document.querySelector('body > h1');
    var meta = document.querySelector('body > .meta');
    var inf = document.querySelector('body > .mc-inf');
    var reserved = 24;
    [chrome, h1, meta, inf].forEach(function (el) {
      if (el) reserved += el.getBoundingClientRect().height;
    });
    var lbn = root.querySelector('.r2dt-compare-lbn, .r2dt-viewer-lbn');
    var lbnBudget = lbn ? Math.min(200, Math.max(96, window.innerHeight * 0.22)) : 0;
    var slotChrome = compareRoot ? 36 : 8;
    var avail = window.innerHeight - reserved - lbnBudget - slotChrome;
    var panelH = Math.max(240, Math.min(520, avail));
    if (compareRoot) {
      compareRoot.style.setProperty('--r2dt-compare-panel-height', panelH + 'px');
    }
    if (viewerRoot) {
      viewerRoot.style.setProperty('--r2dt-viewer-height', panelH + 'px');
      viewerRoot.style.setProperty('--r2dt-panel-size', Math.min(600, panelH) + 'px');
    }
    if (lbn) {
      var used = reserved + panelH + slotChrome + 12;
      lbn.style.maxHeight = Math.max(80, window.innerHeight - used) + 'px';
    }
  }
  window.addEventListener('resize', fitWorkstationViewer);
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', fitWorkstationViewer);
  } else {
    fitWorkstationViewer();
  }
  // createCompare is async; refit as the mount fills in.
  var mount = document.getElementById('r2dt-compare-mount')
    || document.getElementById('r2dt-viewer-mount');
  if (mount && typeof MutationObserver !== 'undefined') {
    var obs = new MutationObserver(function () { fitWorkstationViewer(); });
    obs.observe(mount, { childList: true, subtree: true });
    setTimeout(function () { obs.disconnect(); }, 8000);
  }
  setTimeout(fitWorkstationViewer, 300);
  setTimeout(fitWorkstationViewer, 1200);
})();
</script>
"""


def _inject_viewer_chrome(
    html: str,
    job_id: str,
    label: str = "",
    active_path: str = "compare",
) -> str:
    """Add shared workstation chrome to a viewer index.html (idempotent)."""
    if "ws-chrome" in html:
        return html
    head = FAVICON_LINKS + '<link rel="stylesheet" href="/static/chrome.css">\n'
    if "</head>" in html:
        html = html.replace("</head>", head + "</head>", 1)
    else:
        html = head + html
    # Compare pages have an INF metrics bar — put Export there instead of the
    # header so it sits with the scoring summary.
    has_inf = 'class="mc-inf"' in html
    header = chrome_header(
        active_path=active_path,
        job_label=label,
        job_id=job_id,
        show_export=bool(job_id) and not has_inf,
    )
    match = re.search(r"<body([^>]*)>", html, flags=re.IGNORECASE)
    if match:
        insert_at = match.end()
        html = html[:insert_at] + "\n" + header + html[insert_at:]
    else:
        html = header + html
    if has_inf and job_id:
        export = export_menu_html(job_id, include_html=True, variant="inf")
        html = re.sub(
            r'(<div class="mc-inf">)(.*?)(</div>)',
            rf"\1\2{export}\3",
            html,
            count=1,
            flags=re.DOTALL,
        )
    if "fitWorkstationViewer" not in html:
        if "</body>" in html:
            html = html.replace("</body>", _VIEWER_FIT_SCRIPT + "</body>", 1)
        else:
            html = html + _VIEWER_FIT_SCRIPT
    return html


class WorkstationApp:
    """Shared mutable state for request handlers."""

    def __init__(
        self,
        catalog: Catalog,
        runner: JobRunner,
        docker_image: str,
        repo_root: Path,
    ):
        self.catalog = catalog
        self.runner = runner
        self.docker_image = docker_image
        self.repo_root = Path(repo_root)
        self.viewer_assets = self.repo_root / "data" / "viewer"


# Set by run_server before serve_forever(); read by WorkstationHandler.
_APP: Optional[WorkstationApp] = None


def _json_bytes(payload: Any, status: int = 200) -> Tuple[int, bytes, str]:
    """Encode a JSON response triple."""
    body = json.dumps(payload).encode("utf-8")
    return status, body, "application/json; charset=utf-8"


def run_server(
    workspace: Path,
    repo_root: Path,
    host: str = "127.0.0.1",
    port: int = 8765,
    docker_image: str = "rnacentral/r2dt:latest",
) -> None:
    """Start the workstation HTTP server (blocking)."""
    global _APP  # pylint: disable=global-statement

    ok, message = require_runtime()
    if not ok:
        raise SystemExit(message)

    workspace = workspace.expanduser().resolve()
    workspace.mkdir(parents=True, exist_ok=True)
    catalog = Catalog(workspace)
    runner = JobRunner(catalog, repo_root=repo_root, docker_image=docker_image)
    _APP = WorkstationApp(catalog, runner, docker_image, repo_root)

    server = ThreadingHTTPServer((host, port), WorkstationHandler)
    server.allow_reuse_address = True

    url = f"http://127.0.0.1:{port}/"
    rprint("[green]Local workstation[/green]")
    # Bare URL for copy/paste. Clickable links are unreliable in macOS Terminal
    # when output is piped from Docker; `just workstation` opens the browser.
    print(url, flush=True)
    if host != "127.0.0.1":
        rprint(
            f"[yellow]Bound to {host}:{port} inside the process; "
            "publish only to 127.0.0.1 from Docker. "
            "Mutating API calls still require Host/Origin localhost.[/yellow]"
        )
    rprint(f"[green]Workspace[/green] {workspace}")
    rprint(
        f"[dim]{message}. No auth — trust is loopback only. "
        "Structures stay on this machine. Ctrl+C to stop.[/dim]"
    )

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        rprint("\n[yellow]Workstation stopped.[/yellow]")
    finally:
        server.server_close()
        _APP = None


class WorkstationHandler(BaseHTTPRequestHandler):
    """Route table for the local R2DT workstation."""

    server_version = "R2DTWorkstation/0.1"

    @property
    def app(self) -> WorkstationApp:
        """Return the process-wide app state."""
        if _APP is None:
            raise RuntimeError("workstation app is not initialised")
        return _APP

    # Match BaseHTTPRequestHandler signature (format is a reserved word).
    def log_message(self, format, *args):  # pylint: disable=redefined-builtin
        """Suppress per-request access logging."""

    def do_GET(self):  # pylint: disable=invalid-name
        """Handle GET."""
        self._dispatch("GET")

    def do_POST(self):  # pylint: disable=invalid-name
        """Handle POST."""
        self._dispatch("POST")

    def do_PUT(self):  # pylint: disable=invalid-name
        """Handle PUT."""
        self._dispatch("PUT")

    def do_DELETE(self):  # pylint: disable=invalid-name
        """Handle DELETE."""
        self._dispatch("DELETE")

    def _dispatch(self, method: str) -> None:
        parsed = urlparse(self.path)
        path = parsed.path
        try:
            blocked = local_mutating_request_error(method, self.headers)
            if blocked:
                self._send(*_json_bytes({"error": blocked}, 403))
                return
            if not self._route(method, path, parsed):
                self._send(*_json_bytes({"error": "not found"}, 404))
        except ValueError as exc:
            self._send(*_json_bytes({"error": str(exc)}, 400))
        except Exception as exc:  # pylint: disable=broad-exception-caught
            self._send(*_json_bytes({"error": str(exc)}, 500))

    def _route(self, method: str, path: str, parsed) -> bool:
        # pylint: disable=too-many-return-statements,too-many-branches,too-many-statements
        if method == "GET" and path in ("/", "/index.html"):
            self._serve_home()
            return True
        if method == "GET" and path in ("/compare", "/compare/"):
            self._serve_compare()
            return True
        if method == "GET" and path in ("/compare/new", "/compare/new/"):
            self._serve_compare()
            return True
        if method == "GET" and path in ("/pdb", "/pdb/"):
            self._serve_pdb()
            return True
        if method == "GET" and path in ("/pdb/new", "/pdb/new/"):
            self._serve_pdb()
            return True
        if method == "GET" and path in ("/2d", "/2d/"):
            self._serve_draw()
            return True
        if method == "GET" and path in ("/2d/new", "/2d/new/"):
            self._serve_draw()
            return True
        match = re.fullmatch(r"/align/?", path)
        if method == "GET" and match:
            self._serve_align()
            return True
        match = re.fullmatch(r"/align/new/?", path)
        if method == "GET" and match:
            self._serve_align()
            return True
        match = re.fullmatch(r"/jobs/([^/]+)/results(?:/(.*))?$", path)
        if method == "GET" and match:
            self._serve_results(match.group(1), match.group(2) or "")
            return True
        if method == "GET" and path.startswith("/static/"):
            self._serve_static(path[len("/static/") :])
            return True
        if method == "GET" and path == "/api/runtime":
            self._send(*_json_bytes(self._runtime_payload()))
            return True
        match = re.fullmatch(r"/api/examples/stockholm/([^/]+)", path)
        if method == "GET" and match:
            return self._serve_stockholm_example(match.group(1))
        if method == "GET" and path == "/api/jobs":
            qs = parse_qs(parsed.query)
            mode = (qs.get("mode") or [None])[0]
            jobs = self.app.catalog.list_jobs(mode=mode)
            self._send(*_json_bytes({"jobs": jobs}))
            return True
        match = re.fullmatch(r"/api/jobs/([^/]+)", path)
        if method == "GET" and match:
            return self._get_job(match.group(1))
        match = re.fullmatch(r"/api/jobs/([^/]+)/log", path)
        if method == "GET" and match:
            return self._get_log(match.group(1), parsed)
        match = re.fullmatch(r"/api/jobs/([^/]+)/export/html", path)
        if method == "GET" and match:
            return self._export_job_html(match.group(1))
        match = re.fullmatch(r"/api/jobs/([^/]+)/export", path)
        if method == "GET" and match:
            return self._export_job(match.group(1))
        match = re.fullmatch(r"/api/uploads/([^/]+)/chains", path)
        if method == "GET" and match:
            return self._handle_chains(match.group(1))
        match = re.fullmatch(r"/jobs/([^/]+)/viewer(?:/(.*))?$", path)
        if method == "GET" and match:
            self._serve_viewer(match.group(1), match.group(2) or "index.html")
            return True
        if method == "GET" and path == "/__edit-api/ping":
            self._send(*_json_bytes({"ok": True, "editing": True}))
            return True
        match = re.fullmatch(r"/__edit-api/jobs/([^/]+)/basepairs", path)
        if method == "GET" and match:
            return self._get_basepairs(match.group(1))
        if method == "PUT" and match:
            return self._put_basepairs(match.group(1))
        if method == "POST" and path == "/api/uploads":
            self._handle_upload()
            return True
        if method == "POST" and path == "/api/jobs/import":
            self._handle_import_job()
            return True
        if method == "POST" and path == "/api/jobs":
            self._handle_create_job()
            return True
        match = re.fullmatch(r"/api/jobs/([^/]+)", path)
        if method == "DELETE" and match:
            return self._delete_job(match.group(1))
        return False

    def _serve_home(self) -> None:
        raw = (STATIC_DIR / "home.html").read_text(encoding="utf-8")
        html = _fill_page(raw, active_path=None)
        html = html.replace("<!--MODE_CARDS-->", _mode_cards_html())
        self._send(200, html.encode("utf-8"), "text/html; charset=utf-8")

    def _serve_compare(self) -> None:
        raw = (STATIC_DIR / "compare.html").read_text(encoding="utf-8")
        html = _fill_page(raw, active_path="compare")
        self._send(200, html.encode("utf-8"), "text/html; charset=utf-8")

    def _serve_pdb(self) -> None:
        raw = (STATIC_DIR / "pdb.html").read_text(encoding="utf-8")
        html = _fill_page(raw, active_path="pdb")
        self._send(200, html.encode("utf-8"), "text/html; charset=utf-8")

    def _serve_draw(self) -> None:
        raw = (STATIC_DIR / "draw.html").read_text(encoding="utf-8")
        html = _fill_page(raw, active_path="2d")
        self._send(200, html.encode("utf-8"), "text/html; charset=utf-8")

    def _serve_align(self) -> None:
        raw = (STATIC_DIR / "align.html").read_text(encoding="utf-8")
        html = _fill_page(raw, active_path="align")
        self._send(200, html.encode("utf-8"), "text/html; charset=utf-8")

    def _serve_stockholm_example(self, example_id: str) -> bool:
        rel = EXAMPLE_STOCKHOLM.get(example_id)
        if not rel:
            self._send(*_json_bytes({"error": "unknown example"}, 404))
            return True
        path = (self.app.repo_root / rel).resolve()
        if not path.is_file() or not path_is_within(path, self.app.repo_root):
            self._send(*_json_bytes({"error": "example file missing"}, 404))
            return True
        text = path.read_text(encoding="utf-8")
        self._send(
            *_json_bytes(
                {
                    "id": example_id,
                    "filename": path.name,
                    "stockholm": text,
                }
            )
        )
        return True

    def _serve_coming_soon(self, path_key: str) -> None:
        mode = MODE_BY_PATH.get(path_key)
        if not mode:
            self._send(*_json_bytes({"error": "not found"}, 404))
            return
        raw = (STATIC_DIR / "coming-soon.html").read_text(encoding="utf-8")
        html = _fill_page(
            raw,
            active_path=path_key,
            TITLE=mode["title"],
            BLURB=mode["blurb"],
        )
        self._send(200, html.encode("utf-8"), "text/html; charset=utf-8")

    def _redirect(self, location: str, status: int = 302) -> None:
        self.send_response(status)
        self.send_header("Location", location)
        self.send_header("Content-Length", "0")
        self.send_header("Cache-Control", "no-store")
        self.end_headers()

    def _get_job(self, job_id: str) -> bool:
        job = self.app.catalog.get_job(job_id)
        if not job:
            self._send(*_json_bytes({"error": "not found"}, 404))
            return True
        self._send(*_json_bytes({"job": job}))
        return True

    def _get_log(self, job_id: str, parsed) -> bool:
        if not self.app.catalog.get_job(job_id):
            self._send(*_json_bytes({"error": "not found"}, 404))
            return True
        qs = parse_qs(parsed.query)
        tail = int((qs.get("tail") or ["200"])[0])
        self._send(
            *_json_bytes(
                {
                    "id": job_id,
                    "log": self.app.catalog.read_log(job_id, tail=tail),
                    "current": self.app.runner.current_job_id,
                }
            )
        )
        return True

    def _delete_job(self, job_id: str) -> bool:
        if not self.app.catalog.delete_job(job_id):
            self._send(*_json_bytes({"error": "not found"}, 404))
            return True
        self._send(*_json_bytes({"ok": True}))
        return True

    def _export_job(self, job_id: str) -> bool:
        if not self.app.catalog.get_job(job_id):
            self._send(*_json_bytes({"error": "not found"}, 404))
            return True
        data, filename = export_job_bytes(self.app.catalog, job_id)
        self._send(
            200,
            data,
            "application/zip",
            extra_headers={
                "Content-Disposition": f'attachment; filename="{filename}"',
            },
        )
        return True

    def _export_job_html(self, job_id: str) -> bool:
        if not self.app.catalog.get_job(job_id):
            self._send(*_json_bytes({"error": "not found"}, 404))
            return True
        try:
            data, filename = export_shareable_html_bytes(self.app.catalog, job_id)
        except ValueError as exc:
            self._send(*_json_bytes({"error": str(exc)}, 400))
            return True
        self._send(
            200,
            data,
            "application/zip",
            extra_headers={
                "Content-Disposition": f'attachment; filename="{filename}"',
            },
        )
        return True

    def _handle_import_job(self) -> None:
        length = int(self.headers.get("Content-Length") or "0")
        if length <= 0:
            raise ValueError("Empty package")
        data = self.rfile.read(length)
        label = self.headers.get("X-Label") or self.headers.get("x-label") or ""
        notes = self.headers.get("X-Notes") or self.headers.get("x-notes") or ""
        result = import_job_bytes(
            self.app.catalog,
            data,
            label=label,
            notes=notes,
        )
        self._send(*_json_bytes(result, 201))

    def _runtime_payload(self) -> Dict[str, Any]:
        ok, message = require_runtime()
        return {
            "ok": ok,
            "message": message,
            "docker_image": self.app.docker_image,
            "workspace": str(self.app.catalog.workspace),
            "current_job": self.app.runner.current_job_id,
        }

    def _handle_upload(self) -> None:
        length = int(self.headers.get("Content-Length") or "0")
        if length <= 0:
            raise ValueError("Empty upload")
        filename = self.headers.get("X-Filename") or self.headers.get("x-filename")
        if not filename:
            raise ValueError("Missing X-Filename header")
        filename = Path(filename).name
        if not is_structure_filename(filename):
            raise ValueError("Unsupported file type (use .pdb / .cif, optionally .gz)")
        data = self.rfile.read(length)
        upload_id = uuid.uuid4().hex
        dest_dir = self.app.catalog.uploads_dir / upload_id
        dest_dir.mkdir(parents=True, exist_ok=True)
        dest = dest_dir / filename
        dest.write_bytes(data)
        info = list_rna_chains(dest)
        self._send(
            *_json_bytes(
                {
                    "upload_id": upload_id,
                    "filename": filename,
                    "size": length,
                    **info,
                }
            )
        )

    def _handle_chains(self, upload_id: str) -> bool:
        upload_dir = self.app.catalog.uploads_dir / upload_id
        if not upload_dir.is_dir():
            self._send(*_json_bytes({"error": "not found"}, 404))
            return True
        files = [p for p in upload_dir.iterdir() if p.is_file()]
        if not files:
            self._send(*_json_bytes({"error": "not found"}, 404))
            return True
        info = list_rna_chains(files[0])
        self._send(*_json_bytes({"upload_id": upload_id, **info}))
        return True

    def _handle_create_job(self) -> None:
        length = int(self.headers.get("Content-Length") or "0")
        raw = self.rfile.read(length) if length else b"{}"
        try:
            payload = json.loads(raw.decode("utf-8"))
        except ValueError as exc:
            raise ValueError("Invalid JSON body") from exc
        job_mode = payload.get("job_mode") or "compare"
        advanced = payload.get("advanced")
        if job_mode == "pdb":
            result = create_pdb_job_from_upload(
                self.app.catalog,
                self.app.runner,
                upload_id=payload.get("upload_id") or "",
                chain=payload.get("chain") or "",
                mode=payload.get("layout_mode") or payload.get("mode") or "auto",
                basepairs=payload.get("basepairs") or "fr3d",
                label=payload.get("label") or "",
                notes=payload.get("notes") or "",
                force=bool(payload.get("force")),
                advanced=advanced,
            )
            status = 200 if result.get("dedup") else 201
        elif job_mode == "draw":
            result = create_draw_jobs_from_fasta(
                self.app.catalog,
                self.app.runner,
                fasta_text=payload.get("fasta") or "",
                layout=payload.get("layout") or payload.get("layout_mode") or "auto",
                label=payload.get("label") or "",
                notes=payload.get("notes") or "",
                force=bool(payload.get("force")),
                advanced=advanced,
            )
            status = 201 if result.get("created") else 200
        elif job_mode == "align":
            result = create_align_job_from_stockholm(
                self.app.catalog,
                self.app.runner,
                stockholm_text=payload.get("stockholm") or "",
                stitch=bool(payload.get("stitch", True)),
                label=payload.get("label") or "",
                notes=payload.get("notes") or "",
                force=bool(payload.get("force")),
                advanced=advanced,
            )
            status = 200 if result.get("dedup") else 201
        else:
            layout_mode = payload.get("layout_mode") or payload.get("mode") or "auto"
            result = create_job_from_uploads(
                self.app.catalog,
                self.app.runner,
                ref_upload_id=payload.get("ref_upload_id") or "",
                model_upload_id=payload.get("model_upload_id") or "",
                chains=list(payload.get("chains") or []),
                model_chains=list(payload.get("model_chains") or []),
                mode=layout_mode,
                basepairs=payload.get("basepairs") or "fr3d",
                label=payload.get("label") or "",
                notes=payload.get("notes") or "",
                force=bool(payload.get("force")),
                advanced=advanced,
            )
            status = 200 if result.get("dedup") else 201
        self._send(*_json_bytes(result, status))

    def _get_basepairs(self, job_id: str) -> bool:
        if not self.app.catalog.get_job(job_id):
            self._send(*_json_bytes({"error": "not found"}, 404))
            return True
        data = edits_mod.load_overrides(self.app.catalog.job_dir(job_id))
        self._send(*_json_bytes(data))
        return True

    def _put_basepairs(self, job_id: str) -> bool:
        if not self.app.catalog.get_job(job_id):
            self._send(*_json_bytes({"error": "not found"}, 404))
            return True
        length = int(self.headers.get("Content-Length") or "0")
        raw = self.rfile.read(length) if length else b"{}"
        try:
            payload = json.loads(raw.decode("utf-8"))
        except ValueError as exc:
            raise ValueError("Invalid JSON body") from exc
        saved = edits_mod.save_overrides(self.app.catalog.job_dir(job_id), payload)
        counts = edits_mod.edit_counts(self.app.catalog.job_dir(job_id))
        self.app.catalog.update_meta(
            job_id,
            edits={
                "ref_count": counts["ref_count"],
                "model_count": counts["model_count"],
            },
        )
        self._send(*_json_bytes({"ok": True, **saved, **counts}))
        return True

    def _serve_results(self, job_id: str, rel: str) -> None:
        """Serve draw-mode results page or a file under the job directory."""
        meta = self.app.catalog.read_meta(job_id)
        if not meta:
            self._send(*_json_bytes({"error": "not found"}, 404))
            return
        job_dir = self.app.catalog.job_dir(job_id).resolve()
        if rel in ("", "index.html"):
            self._send_results_page(job_id, meta, job_dir)
            return
        target = (job_dir / rel).resolve()
        if not path_is_within(target, job_dir):
            self._send(*_json_bytes({"error": "forbidden"}, 403))
            return
        if not target.is_file():
            self._send(*_json_bytes({"error": "not found"}, 404))
            return
        self._send_file(target)

    def _send_results_page(
        self, job_id: str, meta: Dict[str, Any], job_dir: Path
    ) -> None:
        label = str(meta.get("label") or job_id)
        job_mode = normalize_job_mode(meta.get("mode"))
        active = next((m["path"] for m in MODES if m["id"] == job_mode), "2d")
        back_href = f"/{active}"
        back_label = next((m["nav"] for m in MODES if m["id"] == job_mode), "jobs")
        status = meta.get("status") or ""
        raw = (STATIC_DIR / "results.html").read_text(encoding="utf-8")
        html = raw.replace("<!--HEAD-->", _HEAD_LINKS)
        html = html.replace(
            "<!--HEADER-->",
            chrome_header(
                active_path=active,
                job_label=label,
                job_id=job_id,
                show_export=True,
                export_html=False,
            ),
        )
        html = html.replace("<!--JOB_LABEL-->", html_lib.escape(label))
        html = html.replace("<!--JOB_ID-->", html_lib.escape(job_id))
        html = html.replace("<!--JOB_STATUS-->", html_lib.escape(str(status)))

        if job_mode == "align":
            preview = self._align_results_body(
                job_id, meta, job_dir, status, back_href, back_label
            )
        else:
            preview = self._draw_results_body(
                job_id, meta, job_dir, status, label, back_href, back_label
            )
        html = html.replace("<!--RESULTS_BODY-->", preview)
        self._send(200, html.encode("utf-8"), "text/html; charset=utf-8")

    def _draw_results_body(  # pylint: disable=too-many-arguments,too-many-positional-arguments
        self,
        job_id: str,
        meta: Dict[str, Any],
        job_dir: Path,
        status: str,
        label: str,
        back_href: str,
        back_label: str,
    ) -> str:
        variants = find_draw_svg_variants(job_dir)
        # Prefer meta.svg_path as colored when it matches a known file.
        if svg_rel := meta.get("svg_path"):
            candidate = job_dir / str(svg_rel)
            if candidate.is_file():
                rel = candidate.relative_to(job_dir).as_posix()
                if not any(v["path"] == rel for v in variants):
                    variants.insert(0, {"id": "primary", "label": "SVG", "path": rel})
        if variants:
            default_id = next(
                (v["id"] for v in variants if v["id"] == "colored"),
                variants[0]["id"],
            )
            tabs = []
            panels = []
            for variant in variants:
                vid = html_lib.escape(variant["id"])
                vlabel = html_lib.escape(variant["label"])
                href = f"/jobs/{job_id}/results/{html_lib.escape(variant['path'])}"
                selected = variant["id"] == default_id
                tabs.append(
                    f'<button type="button" class="ws-svg-tab'
                    f'{" is-active" if selected else ""}" '
                    f'data-tab="{vid}" role="tab" '
                    f'aria-selected="{"true" if selected else "false"}">'
                    f"{vlabel}</button>"
                )
                panels.append(
                    f'<div class="ws-svg-panel'
                    f'{" is-active" if selected else ""}" '
                    f'data-panel="{vid}" role="tabpanel"'
                    f'{" hidden" if not selected else ""}>'
                    f'<div class="ws-svg-wrap">'
                    f'<img class="ws-svg-preview" src="{href}" '
                    f'alt="{html_lib.escape(label)} ({vlabel})">'
                    f"</div></div>"
                )
            default_href = next(
                f"/jobs/{job_id}/results/{v['path']}"
                for v in variants
                if v["id"] == default_id
            )
            tablist = ""
            if len(variants) > 1:
                tablist = (
                    '<div class="ws-svg-tabs" role="tablist" '
                    'aria-label="SVG variants">' + "".join(tabs) + "</div>"
                )
            return (
                f'<div class="ws-svg-stage" data-default-tab="{html_lib.escape(default_id)}">'
                f"{tablist}"
                f'<div class="ws-svg-panels">{"".join(panels)}</div>'
                f"</div>"
                f'<p class="ws-results-actions">'
                f'<a class="cta" id="ws-svg-download" '
                f'href="{html_lib.escape(default_href)}" download>'
                f"Download SVG</a>"
                f'<a class="ws-stub-link" href="{html_lib.escape(back_href)}">'
                f"← Back to {html_lib.escape(back_label)} jobs</a>"
                f"</p>"
            )
        if status in ("queued", "running"):
            return (
                '<p class="sub">Job still running — refresh this page shortly.</p>'
                f'<p><a class="ws-stub-link" href="{html_lib.escape(back_href)}">'
                f"← Back to {html_lib.escape(back_label)} jobs</a></p>"
            )
        err = html_lib.escape(str(meta.get("error") or "No SVG produced"))
        return (
            f'<p class="ws-coming">{err}</p>'
            f'<p><a class="ws-stub-link" href="{html_lib.escape(back_href)}">'
            f"← Back to {html_lib.escape(back_label)} jobs</a></p>"
        )

    def _align_results_body(  # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals
        self,
        job_id: str,
        meta: Dict[str, Any],
        job_dir: Path,
        status: str,
        back_href: str,
        back_label: str,
    ) -> str:
        gallery = meta.get("svg_gallery") or []
        if not gallery:
            paths = list_align_result_svgs(job_dir)
            gallery = [svg_gallery_entry(p, job_dir) for p in paths]
        if gallery:
            tiles = []
            for entry in gallery:
                rel = html_lib.escape(str(entry.get("path") or ""))
                href = f"/jobs/{job_id}/results/{rel}"
                caption = html_lib.escape(
                    str(entry.get("caption") or entry.get("name"))
                )
                kind = html_lib.escape(str(entry.get("kind") or ""))
                tiles.append(
                    f'<figure class="ws-gallery-tile" data-kind="{kind}">'
                    f"<figcaption>{caption}</figcaption>"
                    f'<div class="ws-svg-wrap">'
                    f'<img class="ws-svg-preview" src="{href}" alt="{caption}">'
                    f"</div>"
                    f'<a class="ws-stub-link" href="{href}" download>Download</a>'
                    f"</figure>"
                )
            note = ""
            params = meta.get("params") or {}
            if params.get("ran_rscape"):
                note = (
                    '<p class="sub">R-scape was run automatically '
                    "(no cov_* annotations in the input).</p>"
                )
            elif params.get("has_covariation"):
                note = (
                    '<p class="sub">Input already included R-scape '
                    "covariation annotations.</p>"
                )
            return (
                note
                + f'<div class="ws-gallery">{"".join(tiles)}</div>'
                + '<p class="ws-results-actions">'
                + f'<a class="ws-stub-link" href="{html_lib.escape(back_href)}">'
                + f"← Back to {html_lib.escape(back_label)} jobs</a>"
                + "</p>"
            )
        if status in ("queued", "running"):
            return (
                '<p class="sub">Job still running — refresh this page shortly.</p>'
                f'<p><a class="ws-stub-link" href="{html_lib.escape(back_href)}">'
                f"← Back to {html_lib.escape(back_label)} jobs</a></p>"
            )
        err = html_lib.escape(str(meta.get("error") or "No SVG produced"))
        return (
            f'<p class="ws-coming">{err}</p>'
            f'<p><a class="ws-stub-link" href="{html_lib.escape(back_href)}">'
            f"← Back to {html_lib.escape(back_label)} jobs</a></p>"
        )

    def _serve_viewer(self, job_id: str, rel: str) -> None:
        # Always serve the latest R2DT glue/CSS from the repo so existing jobs
        # pick up editing hooks without regenerating.
        live_names = {
            "r2dt-2d-3d-viewer.js",
            "r2dt-2d-3d-viewer.css",
        }
        if (name := Path(rel).name) in live_names:
            live = self.app.viewer_assets / name
            if live.is_file():
                self._send_file(live)
                return
        base = (self.app.catalog.job_dir(job_id) / "viewer").resolve()
        if not base.is_dir():
            self._send(*_json_bytes({"error": "viewer not ready"}, 404))
            return
        target = (base / rel).resolve()
        if not path_is_within(target, base):
            self._send(*_json_bytes({"error": "forbidden"}, 403))
            return
        if target.is_dir():
            target = target / "index.html"
        if not target.is_file():
            self._send(*_json_bytes({"error": "not found"}, 404))
            return
        if target.name == "index.html":
            meta = self.app.catalog.read_meta(job_id) or {}
            label = str(meta.get("label") or job_id)
            job_mode = normalize_job_mode(meta.get("mode"))
            active = next((m["path"] for m in MODES if m["id"] == job_mode), "compare")
            html = _inject_viewer_chrome(
                target.read_text(encoding="utf-8"),
                job_id,
                label,
                active_path=active,
            )
            self._send(200, html.encode("utf-8"), "text/html; charset=utf-8")
            return
        self._send_file(target)

    def _serve_static(self, rel: str) -> None:
        base = STATIC_DIR.resolve()
        target = (STATIC_DIR / rel).resolve()
        if not path_is_within(target, base):
            self._send(*_json_bytes({"error": "forbidden"}, 403))
            return
        if not target.is_file():
            self._send(*_json_bytes({"error": "not found"}, 404))
            return
        self._send_file(target)

    def _send_file(self, path: Path) -> None:
        ctype = mimetypes.guess_type(str(path))[0] or "application/octet-stream"
        data = path.read_bytes()
        self._send(200, data, ctype)

    def _send(
        self,
        status: int,
        body: bytes,
        content_type: str,
        extra_headers: Optional[Dict[str, str]] = None,
    ) -> None:
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        if extra_headers:
            for key, value in extra_headers.items():
                self.send_header(key, value)
        self.end_headers()
        self.wfile.write(body)
