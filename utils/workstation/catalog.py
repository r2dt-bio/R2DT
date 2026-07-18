"""Workspace catalog: scan jobs/, read/write catalog.json and meta.json."""

from __future__ import annotations

import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

from utils.workstation.chrome import normalize_job_mode


def utc_now() -> str:
    """Return an ISO-8601 UTC timestamp."""
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


class Catalog:
    """Filesystem-backed job catalog under a workspace root."""

    def __init__(self, workspace: Path):
        self.workspace = workspace
        self.catalog_path = workspace / "catalog.json"
        self.jobs_dir = workspace / "jobs"
        self.inputs_dir = workspace / "inputs"
        self.uploads_dir = workspace / "uploads"
        self.jobs_dir.mkdir(parents=True, exist_ok=True)
        self.inputs_dir.mkdir(parents=True, exist_ok=True)
        self.uploads_dir.mkdir(parents=True, exist_ok=True)

    def list_jobs(self, mode: Optional[str] = None) -> List[Dict[str, Any]]:
        """Return jobs, newest first. Optional ``mode`` filters by job type."""
        jobs = []
        for child in sorted(self.jobs_dir.iterdir(), reverse=True):
            if not child.is_dir():
                continue
            meta = self.read_meta(child.name)
            if meta:
                jobs.append(meta)
        self._write_catalog(jobs)
        if mode:
            want = normalize_job_mode(mode)
            return [j for j in jobs if normalize_job_mode(j.get("mode")) == want]
        return jobs

    def get_job(self, job_id: str) -> Optional[Dict[str, Any]]:
        """Return one job's meta, or None."""
        return self.read_meta(job_id)

    def job_dir(self, job_id: str) -> Path:
        """Return ``jobs/<id>/``."""
        return self.jobs_dir / job_id

    def inputs_job_dir(self, job_id: str) -> Path:
        """Return ``inputs/<id>/``."""
        return self.inputs_dir / job_id

    def read_meta(self, job_id: str) -> Optional[Dict[str, Any]]:
        """Load ``jobs/<id>/meta.json``; tag legacy jobs as compare."""
        path = self.job_dir(job_id) / "meta.json"
        if not path.is_file():
            return None
        try:
            meta = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            return None
        if not isinstance(meta, dict):
            return None
        normalized = normalize_job_mode(meta.get("mode"))
        if meta.get("mode") != normalized:
            meta = dict(meta)
            meta["mode"] = normalized
            meta["id"] = job_id
            self._persist_meta_file(job_id, meta)
        return meta

    def _persist_meta_file(self, job_id: str, meta: Dict[str, Any]) -> None:
        """Write meta.json without refreshing the catalog (avoids recursion)."""
        job_path = self.job_dir(job_id)
        job_path.mkdir(parents=True, exist_ok=True)
        (job_path / "meta.json").write_text(
            json.dumps(meta, indent=2) + "\n", encoding="utf-8"
        )

    def write_meta(self, job_id: str, meta: Dict[str, Any]) -> None:
        """Persist meta.json and refresh catalog.json."""
        meta = dict(meta)
        meta["id"] = job_id
        meta["mode"] = normalize_job_mode(meta.get("mode"))
        self._persist_meta_file(job_id, meta)
        self.list_jobs()

    def update_meta(self, job_id: str, **fields: Any) -> Optional[Dict[str, Any]]:
        """Merge fields into meta.json."""
        meta = self.read_meta(job_id)
        if meta is None:
            return None
        meta.update(fields)
        self.write_meta(job_id, meta)
        return meta

    def delete_job(self, job_id: str) -> bool:
        """Remove job and input dirs. Returns False if unknown."""
        job_path = self.job_dir(job_id)
        if not job_path.exists():
            return False
        shutil.rmtree(job_path, ignore_errors=True)
        inputs = self.inputs_job_dir(job_id)
        if inputs.exists():
            shutil.rmtree(inputs, ignore_errors=True)
        self.list_jobs()
        return True

    def find_by_content_hash(self, content_hash: str) -> Optional[Dict[str, Any]]:
        """Return a ready job with the same content hash, if any."""
        for meta in self.list_jobs():
            if meta.get("status") != "ready":
                continue
            inputs = meta.get("inputs") or {}
            if inputs.get("content_hash") == content_hash:
                return meta
        return None

    def refresh_metrics(self, job_id: str) -> Optional[Dict[str, Any]]:
        """Copy INF / BP-diff from viewer/metrics.json into meta."""
        metrics_path = self.job_dir(job_id) / "viewer" / "metrics.json"
        if not metrics_path.is_file():
            return self.read_meta(job_id)
        try:
            raw = json.loads(metrics_path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            return self.read_meta(job_id)

        def _inf_val(key: str) -> Any:
            block = (raw.get("inf") or {}).get(key) or {}
            if isinstance(block, dict):
                return block.get("inf")
            return block

        diff = raw.get("diff") or {}
        metrics = {
            "inf": {
                "wc": _inf_val("wc"),
                "nwc": _inf_val("nwc"),
                "all": _inf_val("all"),
            },
            "matched": diff.get("matched"),
            "lost": diff.get("lost"),
            "added": diff.get("added"),
            "superpose_rmsd": raw.get("superpose_rmsd"),
        }
        return self.update_meta(job_id, metrics=metrics)

    def append_log(self, job_id: str, text: str) -> None:
        """Append to ``jobs/<id>/log.txt``."""
        log_path = self.job_dir(job_id) / "log.txt"
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with log_path.open("a", encoding="utf-8") as handle:
            handle.write(text)
            if text and not text.endswith("\n"):
                handle.write("\n")

    def read_log(self, job_id: str, tail: int = 200) -> str:
        """Return the last ``tail`` lines of the job log."""
        log_path = self.job_dir(job_id) / "log.txt"
        if not log_path.is_file():
            return ""
        lines = log_path.read_text(encoding="utf-8", errors="replace").splitlines()
        return "\n".join(lines[-tail:])

    def _write_catalog(self, jobs: List[Dict[str, Any]]) -> None:
        payload = {"updated": utc_now(), "jobs": jobs}
        self.catalog_path.write_text(
            json.dumps(payload, indent=2) + "\n", encoding="utf-8"
        )
