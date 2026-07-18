"""Docker / runtime checks and job execution for the workstation."""

from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
import threading
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from utils.workstation.catalog import Catalog, utc_now
from utils.workstation.chains import list_rna_chains, normalize_suffix


def _resolve_upload(catalog: Catalog, upload_id: str) -> Path:
    """Return the path to a previously stored upload, or raise."""
    # upload_id is the directory name under uploads/
    upload_dir = catalog.uploads_dir / upload_id
    if not upload_dir.is_dir():
        raise ValueError(f"Unknown upload id: {upload_id}")
    files = [p for p in upload_dir.iterdir() if p.is_file()]
    if len(files) != 1:
        raise ValueError(f"Upload {upload_id} is incomplete")
    return files[0]


def in_docker() -> bool:
    """True when running inside a container."""
    return Path("/.dockerenv").exists()


def docker_available() -> bool:
    """True when the docker CLI can talk to a daemon."""
    if shutil.which("docker") is None:
        return False
    try:
        result = subprocess.run(
            ["docker", "info"],
            capture_output=True,
            check=False,
            timeout=20,
        )
        return result.returncode == 0
    except (OSError, subprocess.TimeoutExpired):
        return False


def require_runtime() -> Tuple[bool, str]:
    """Return (ok, message). Docker is always required."""
    if in_docker():
        return True, "running inside Docker"
    if docker_available():
        return True, "Docker daemon available on host"
    return (
        False,
        "Docker is required. Start the workstation with: just workstation",
    )


def new_job_id() -> str:
    """Return a filesystem-safe unique job id."""
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
    return f"{stamp}-{uuid.uuid4().hex[:6]}"


def content_hash(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    ref_path: Path,
    model_path: Path,
    chains: str,
    model_chains: str,
    mode: str,
    basepairs: str,
) -> str:
    """Stable hash of inputs + generate params (edits are not included)."""
    digest = hashlib.sha256()
    for path in (ref_path, model_path):
        digest.update(path.name.encode("utf-8"))
        digest.update(b"\0")
        with path.open("rb") as handle:
            while True:
                chunk = handle.read(1024 * 1024)
                if not chunk:
                    break
                digest.update(chunk)
        digest.update(b"\0")
    for part in (chains, model_chains, mode, basepairs):
        digest.update(part.encode("utf-8"))
        digest.update(b"\0")
    return "sha256:" + digest.hexdigest()


class JobRunner:
    """Serial job queue: one compare generate at a time."""

    def __init__(
        self,
        catalog: Catalog,
        repo_root: Path,
        docker_image: str = "rnacentral/r2dt:latest",
    ):
        self.catalog = catalog
        self.repo_root = repo_root
        self.docker_image = docker_image
        self._lock = threading.Lock()
        self._thread: Optional[threading.Thread] = None
        self._current: Optional[str] = None
        self._queue: List[str] = []

    @property
    def current_job_id(self) -> Optional[str]:
        """Job id currently running, if any."""
        return self._current

    def enqueue(self, job_id: str) -> None:
        """Append a job and kick the worker if idle."""
        with self._lock:
            if job_id not in self._queue and job_id != self._current:
                self._queue.append(job_id)
            self._ensure_worker()

    def _ensure_worker(self) -> None:
        if self._thread and self._thread.is_alive():
            return
        self._thread = threading.Thread(target=self._worker_loop, daemon=True)
        self._thread.start()

    def _worker_loop(self) -> None:
        while True:
            with self._lock:
                if not self._queue:
                    self._current = None
                    return
                job_id = self._queue.pop(0)
                self._current = job_id
            try:
                self._run_job(job_id)
            except Exception as exc:  # pylint: disable=broad-exception-caught
                self.catalog.append_log(job_id, f"ERROR: {exc}")
                self.catalog.update_meta(
                    job_id, status="failed", error=str(exc), finished=utc_now()
                )

    def _run_job(self, job_id: str) -> None:
        meta = self.catalog.read_meta(job_id)
        if not meta:
            return
        self.catalog.update_meta(
            job_id, status="running", started=utc_now(), error=None
        )
        self.catalog.append_log(job_id, f"Starting job {job_id}")

        inputs_dir = self.catalog.inputs_job_dir(job_id)
        out_dir = self.catalog.job_dir(job_id)
        params = meta.get("params") or {}
        inputs = meta.get("inputs") or {}
        ref_name = inputs["ref_name"]
        model_name = inputs["model_name"]
        ref_path = inputs_dir / ref_name
        model_path = inputs_dir / model_name

        chains = params.get("chains") or ""
        model_chains = params.get("model_chains") or ""
        mode = params.get("mode") or "auto"
        basepairs = params.get("basepairs") or "fr3d"

        cmd = self._build_cmd(
            ref_path=ref_path,
            model_path=model_path,
            out_dir=out_dir,
            chains=chains,
            model_chains=model_chains,
            mode=mode,
            basepairs=basepairs,
        )
        self.catalog.append_log(job_id, "Command: " + " ".join(cmd))

        env = os.environ.copy()
        # Stream logs line-by-line; context-manager wait would block the pipe.
        process = subprocess.Popen(  # pylint: disable=consider-using-with
            cmd,
            cwd=str(self.repo_root),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            env=env,
        )
        assert process.stdout is not None
        for line in process.stdout:
            self.catalog.append_log(job_id, line.rstrip("\n"))
        code = process.wait()

        viewer = out_dir / "viewer" / "index.html"
        if code == 0 and viewer.is_file():
            self.catalog.refresh_metrics(job_id)
            self.catalog.update_meta(
                job_id,
                status="ready",
                finished=utc_now(),
                viewer_url=f"/jobs/{job_id}/viewer/",
                error=None,
            )
            self.catalog.append_log(job_id, "Job completed successfully")
            return

        self.catalog.update_meta(
            job_id,
            status="failed",
            finished=utc_now(),
            error=f"pipeline exited {code}" if code else "viewer/index.html missing",
        )
        self.catalog.append_log(job_id, f"Job failed (exit {code})")

    def _build_cmd(  # pylint: disable=too-many-arguments,too-many-positional-arguments
        self,
        ref_path: Path,
        model_path: Path,
        out_dir: Path,
        chains: str,
        model_chains: str,
        mode: str,
        basepairs: str,
    ) -> List[str]:
        """Build the r2dt.py pdb --compare --model command (in- or out-of-docker)."""
        if in_docker():
            return self._local_pdb_cmd(
                ref_path, model_path, out_dir, chains, model_chains, mode, basepairs
            )
        # Host: run the same command inside a one-shot container with mounts.
        workspace = self.catalog.workspace.resolve()
        repo = self.repo_root.resolve()
        rel_ref = ref_path.resolve().relative_to(workspace)
        rel_model = model_path.resolve().relative_to(workspace)
        rel_out = out_dir.resolve().relative_to(workspace)
        inner = [
            "python3",
            "r2dt.py",
            "pdb",
            f"/workspace/{rel_ref.as_posix()}",
            f"/workspace/{rel_out.as_posix()}",
            "--chains",
            chains,
            "--compare",
            "--model",
            f"/workspace/{rel_model.as_posix()}",
            "--mode",
            mode,
            "--basepairs",
            basepairs,
            "--quiet",
        ]
        if model_chains:
            inner.extend(["--model-chains", model_chains])
        return [
            "docker",
            "run",
            "--rm",
            "-v",
            f"{repo}:/rna/r2dt",
            "-v",
            f"{workspace}:/workspace",
            "-w",
            "/rna/r2dt",
            self.docker_image,
            *inner,
        ]

    @staticmethod
    def _local_pdb_cmd(  # pylint: disable=too-many-arguments,too-many-positional-arguments
        ref_path: Path,
        model_path: Path,
        out_dir: Path,
        chains: str,
        model_chains: str,
        mode: str,
        basepairs: str,
    ) -> List[str]:
        cmd = [
            "python3",
            "r2dt.py",
            "pdb",
            str(ref_path),
            str(out_dir),
            "--chains",
            chains,
            "--compare",
            "--model",
            str(model_path),
            "--mode",
            mode,
            "--basepairs",
            basepairs,
            "--quiet",
        ]
        if model_chains:
            cmd.extend(["--model-chains", model_chains])
        return cmd


def create_job_from_uploads(  # pylint: disable=too-many-arguments,too-many-locals
    catalog: Catalog,
    runner: JobRunner,
    *,
    ref_upload_id: str,
    model_upload_id: str,
    chains: List[str],
    model_chains: List[str],
    mode: str = "auto",
    basepairs: str = "fr3d",
    label: str = "",
    notes: str = "",
    force: bool = False,
) -> Dict[str, Any]:
    """Materialise inputs, write meta, enqueue. Returns job meta."""
    if not chains:
        raise ValueError("Select at least one reference chain")
    if not model_chains:
        raise ValueError("Select at least one model chain")
    if len(chains) != len(model_chains):
        raise ValueError(
            "Reference and model must have the same number of chains "
            f"(got {len(chains)} vs {len(model_chains)})"
        )
    if not ref_upload_id or not model_upload_id:
        raise ValueError("Missing upload id(s)")

    ref_src = _resolve_upload(catalog, ref_upload_id)
    model_src = _resolve_upload(catalog, model_upload_id)

    ref_info = list_rna_chains(ref_src)
    if not ref_info.get("compare_ready"):
        raise ValueError(
            "Compare mode needs an mmCIF reference (.cif). "
            "Convert the reference or download the mmCIF from RCSB."
        )

    chains_csv = ",".join(chains)
    model_chains_csv = ",".join(model_chains)
    digest = content_hash(
        ref_src, model_src, chains_csv, model_chains_csv, mode, basepairs
    )
    if not force:
        existing = catalog.find_by_content_hash(digest)
        if existing:
            return {
                "dedup": True,
                "job": existing,
                "message": "Identical comparison already cached",
            }

    job_id = new_job_id()
    inputs_dir = catalog.inputs_job_dir(job_id)
    inputs_dir.mkdir(parents=True, exist_ok=True)
    ref_dest = inputs_dir / ref_src.name
    model_dest = inputs_dir / model_src.name
    shutil.copy2(ref_src, ref_dest)
    shutil.copy2(model_src, model_dest)

    display = label.strip() or f"{ref_src.stem} vs {model_src.stem}"
    meta: Dict[str, Any] = {
        "id": job_id,
        "label": display,
        "notes": notes.strip(),
        "created": utc_now(),
        "status": "queued",
        "params": {
            "chains": chains_csv,
            "model_chains": model_chains_csv,
            "mode": mode,
            "basepairs": basepairs,
        },
        "inputs": {
            "ref_name": ref_dest.name,
            "model_name": model_dest.name,
            "ref_format": normalize_suffix(ref_dest),
            "model_format": normalize_suffix(model_dest),
            "content_hash": digest,
        },
        "metrics": {},
        "viewer_url": f"/jobs/{job_id}/viewer/",
    }
    catalog.write_meta(job_id, meta)
    catalog.append_log(job_id, "Queued")
    runner.enqueue(job_id)
    return {"dedup": False, "job": meta}
