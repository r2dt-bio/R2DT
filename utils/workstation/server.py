"""Stdlib HTTP server for the local curator workstation."""

from __future__ import annotations

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
from utils.workstation.catalog import Catalog
from utils.workstation.chains import is_structure_filename, list_rna_chains
from utils.workstation.jobs import (
    JobRunner,
    create_job_from_uploads,
    require_runtime,
)

STATIC_DIR = Path(__file__).resolve().parent / "static"


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

    rprint(f"[green]Local workstation: http://127.0.0.1:{port}/[/green]")
    if host != "127.0.0.1":
        rprint(
            f"[yellow]Bound to {host}:{port} inside the process; "
            "publish only to 127.0.0.1 from Docker.[/yellow]"
        )
    rprint(f"[green]Workspace: {workspace}[/green]")
    rprint(f"[dim]{message}. Structures stay on this machine. Ctrl+C to stop.[/dim]")

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        rprint("\n[yellow]Workstation stopped.[/yellow]")
    finally:
        server.server_close()
        _APP = None


class WorkstationHandler(BaseHTTPRequestHandler):
    """Route table for the curator workstation."""

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
            if not self._route(method, path, parsed):
                self._send(*_json_bytes({"error": "not found"}, 404))
        except ValueError as exc:
            self._send(*_json_bytes({"error": str(exc)}, 400))
        except Exception as exc:  # pylint: disable=broad-exception-caught
            self._send(*_json_bytes({"error": str(exc)}, 500))

    def _route(self, method: str, path: str, parsed) -> bool:
        # pylint: disable=too-many-return-statements,too-many-branches
        if method == "GET" and path in ("/", "/index.html"):
            self._send_file(STATIC_DIR / "index.html")
            return True
        if method == "GET" and path.startswith("/static/"):
            self._serve_static(path[len("/static/") :])
            return True
        if method == "GET" and path == "/api/runtime":
            self._send(*_json_bytes(self._runtime_payload()))
            return True
        if method == "GET" and path == "/api/jobs":
            self._send(*_json_bytes({"jobs": self.app.catalog.list_jobs()}))
            return True
        match = re.fullmatch(r"/api/jobs/([^/]+)", path)
        if method == "GET" and match:
            return self._get_job(match.group(1))
        match = re.fullmatch(r"/api/jobs/([^/]+)/log", path)
        if method == "GET" and match:
            return self._get_log(match.group(1), parsed)
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
        if method == "POST" and path == "/api/jobs":
            self._handle_create_job()
            return True
        match = re.fullmatch(r"/api/jobs/([^/]+)", path)
        if method == "DELETE" and match:
            return self._delete_job(match.group(1))
        return False

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
        result = create_job_from_uploads(
            self.app.catalog,
            self.app.runner,
            ref_upload_id=payload.get("ref_upload_id") or "",
            model_upload_id=payload.get("model_upload_id") or "",
            chains=list(payload.get("chains") or []),
            model_chains=list(payload.get("model_chains") or []),
            mode=payload.get("mode") or "auto",
            basepairs=payload.get("basepairs") or "fr3d",
            label=payload.get("label") or "",
            notes=payload.get("notes") or "",
            force=bool(payload.get("force")),
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
        if not str(target).startswith(str(base)):
            self._send(*_json_bytes({"error": "forbidden"}, 403))
            return
        if target.is_dir():
            target = target / "index.html"
        if not target.is_file():
            self._send(*_json_bytes({"error": "not found"}, 404))
            return
        self._send_file(target)

    def _serve_static(self, rel: str) -> None:
        target = (STATIC_DIR / rel).resolve()
        if not str(target).startswith(str(STATIC_DIR.resolve())):
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

    def _send(self, status: int, body: bytes, content_type: str) -> None:
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)
