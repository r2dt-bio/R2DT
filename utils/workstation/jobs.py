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

from utils.workstation.align import (
    find_cacofold_r2r_sto,
    has_covariation_annotations,
    list_align_result_svgs,
    svg_gallery_entry,
)
from utils.workstation.catalog import Catalog, utc_now
from utils.workstation.chains import ensure_mmcif, normalize_suffix, structure_stem
from utils.workstation.fasta import FastaRecord, parse_fasta_records


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
    """Stable hash of compare inputs + generate params (edits are not included)."""
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


def pdb_content_hash(
    structure_path: Path,
    chain: str,
    mode: str,
    basepairs: str,
) -> str:
    """Stable hash for a single-structure 2D+3D job."""
    digest = hashlib.sha256()
    digest.update(b"pdb\0")
    digest.update(structure_path.name.encode("utf-8"))
    digest.update(b"\0")
    with structure_path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    digest.update(b"\0")
    for part in (chain, mode, basepairs):
        digest.update(part.encode("utf-8"))
        digest.update(b"\0")
    return "sha256:" + digest.hexdigest()


class JobRunner:
    """Serial job queue: one generate at a time."""

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

        job_mode = meta.get("mode") or "compare"
        if job_mode == "pdb":
            self._run_pdb_job(job_id, meta)
        elif job_mode == "draw":
            self._run_draw_job(job_id, meta)
        elif job_mode == "align":
            self._run_align_job(job_id, meta)
        else:
            self._run_compare_job(job_id, meta)

    def _finish_viewer_job(self, job_id: str, out_dir: Path, code: int) -> None:
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

    def _stream_command(
        self, job_id: str, cmd: List[str], *, cwd: Optional[Path] = None
    ) -> int:
        self.catalog.append_log(job_id, "Command: " + " ".join(cmd))
        env = os.environ.copy()
        workdir = str(cwd) if cwd is not None else str(self.repo_root)
        process = subprocess.Popen(  # pylint: disable=consider-using-with
            cmd,
            cwd=workdir,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            env=env,
        )
        assert process.stdout is not None
        for line in process.stdout:
            self.catalog.append_log(job_id, line.rstrip("\n"))
        return process.wait()

    def _run_compare_job(self, job_id: str, meta: Dict[str, Any]) -> None:
        inputs_dir = self.catalog.inputs_job_dir(job_id)
        out_dir = self.catalog.job_dir(job_id)
        params = meta.get("params") or {}
        inputs = meta.get("inputs") or {}
        ref_path = inputs_dir / inputs["ref_name"]
        model_path = inputs_dir / inputs["model_name"]
        cmd = self._build_compare_cmd(
            ref_path=ref_path,
            model_path=model_path,
            out_dir=out_dir,
            chains=params.get("chains") or "",
            model_chains=params.get("model_chains") or "",
            mode=params.get("mode") or "auto",
            basepairs=params.get("basepairs") or "fr3d",
        )
        code = self._stream_command(job_id, cmd)
        self._finish_viewer_job(job_id, out_dir, code)

    def _run_pdb_job(self, job_id: str, meta: Dict[str, Any]) -> None:
        inputs_dir = self.catalog.inputs_job_dir(job_id)
        out_dir = self.catalog.job_dir(job_id)
        params = meta.get("params") or {}
        inputs = meta.get("inputs") or {}
        structure_path = inputs_dir / inputs["structure_name"]
        cmd = self._build_pdb_cmd(
            structure_path=structure_path,
            out_dir=out_dir,
            chain=params.get("chain") or "",
            mode=params.get("mode") or "auto",
            basepairs=params.get("basepairs") or "fr3d",
        )
        code = self._stream_command(job_id, cmd)
        self._finish_viewer_job(job_id, out_dir, code)

    def _run_draw_job(self, job_id: str, meta: Dict[str, Any]) -> None:
        inputs_dir = self.catalog.inputs_job_dir(job_id)
        out_dir = self.catalog.job_dir(job_id)
        params = meta.get("params") or {}
        inputs = meta.get("inputs") or {}
        fasta_path = inputs_dir / inputs["fasta_name"]
        layout = params.get("layout") or "auto"
        cmd = self._build_draw_cmd(
            fasta_path=fasta_path,
            out_dir=out_dir,
            layout=layout,
        )
        code = self._stream_command(job_id, cmd)
        self._finish_draw_job(job_id, out_dir, code)

    def _run_align_job(self, job_id: str, meta: Dict[str, Any]) -> None:
        inputs_dir = self.catalog.inputs_job_dir(job_id)
        out_dir = self.catalog.job_dir(job_id)
        params = meta.get("params") or {}
        inputs = meta.get("inputs") or {}
        stockholm_path = inputs_dir / inputs["stockholm_name"]
        stitch = bool(params.get("stitch", True))

        text = stockholm_path.read_text(encoding="utf-8", errors="replace")
        if has_covariation_annotations(text):
            self.catalog.append_log(
                job_id, "R-scape covariation annotations present — skipping R-scape"
            )
            sto_for_r2dt = stockholm_path
            ran_rscape = False
        else:
            self.catalog.append_log(
                job_id, "No cov_* annotations — running R-scape --cacofold"
            )
            sto_for_r2dt = self._run_rscape(job_id, stockholm_path, out_dir)
            if sto_for_r2dt is None:
                self.catalog.update_meta(
                    job_id,
                    status="failed",
                    finished=utc_now(),
                    error="R-scape did not produce a *.cacofold.R2R.sto file",
                )
                self.catalog.append_log(job_id, "Job failed (R-scape)")
                return
            ran_rscape = True

        self.catalog.update_meta(
            job_id,
            params={
                **params,
                "ran_rscape": ran_rscape,
                "stockholm_used": sto_for_r2dt.name,
            },
        )
        cmd = self._build_stockholm_cmd(
            stockholm_path=sto_for_r2dt,
            out_dir=out_dir,
            stitch=stitch,
        )
        code = self._stream_command(job_id, cmd)
        self._finish_align_job(job_id, out_dir, code)

    def _run_rscape(
        self, job_id: str, stockholm_path: Path, out_dir: Path
    ) -> Optional[Path]:
        """Run R-scape --cacofold; return path to *.cacofold.R2R.sto or None."""
        rscape_dir = out_dir / "rscape"
        rscape_dir.mkdir(parents=True, exist_ok=True)
        local_sto = rscape_dir / "alignment.stk"
        shutil.copy2(stockholm_path, local_sto)
        cmd = self._build_rscape_cmd(rscape_dir)
        code = self._stream_command(
            job_id, cmd, cwd=rscape_dir if in_docker() else None
        )
        if code != 0:
            self.catalog.append_log(job_id, f"R-scape exited {code}")
            return None
        found = find_cacofold_r2r_sto(rscape_dir)
        if found is None:
            self.catalog.append_log(
                job_id, "R-scape finished but no *.cacofold.R2R.sto found"
            )
            return None
        self.catalog.append_log(job_id, f"Using R-scape output {found.name}")
        return found

    def _finish_align_job(self, job_id: str, out_dir: Path, code: int) -> None:
        svgs = list_align_result_svgs(out_dir)
        if code == 0 and svgs:
            gallery = [svg_gallery_entry(p, out_dir) for p in svgs]
            primary = gallery[0]["path"]
            self.catalog.update_meta(
                job_id,
                status="ready",
                finished=utc_now(),
                results_url=f"/jobs/{job_id}/results/",
                svg_path=primary,
                svg_gallery=gallery,
                error=None,
            )
            self.catalog.append_log(
                job_id, f"Job completed successfully ({len(gallery)} SVG(s))"
            )
            return
        self.catalog.update_meta(
            job_id,
            status="failed",
            finished=utc_now(),
            error=(
                f"pipeline exited {code}"
                if code
                else "no SVG outputs under results/svg/"
            ),
        )
        self.catalog.append_log(job_id, f"Job failed (exit {code})")

    def _finish_draw_job(self, job_id: str, out_dir: Path, code: int) -> None:
        svg = find_draw_svg(out_dir)
        if code == 0 and svg is not None:
            rel = svg.relative_to(out_dir).as_posix()
            self.catalog.update_meta(
                job_id,
                status="ready",
                finished=utc_now(),
                results_url=f"/jobs/{job_id}/results/",
                svg_path=rel,
                error=None,
            )
            self.catalog.append_log(job_id, f"Job completed successfully ({rel})")
            return
        self.catalog.update_meta(
            job_id,
            status="failed",
            finished=utc_now(),
            error=(
                f"pipeline exited {code}" if code else "no .colored.svg found in output"
            ),
        )
        self.catalog.append_log(job_id, f"Job failed (exit {code})")

    def _build_draw_cmd(
        self,
        fasta_path: Path,
        out_dir: Path,
        layout: str,
    ) -> List[str]:
        """Build r2dt.py draw or templatefree (in- or out-of-docker)."""
        if in_docker():
            return self._local_draw_cmd(fasta_path, out_dir, layout)
        workspace = self.catalog.workspace.resolve()
        repo = self.repo_root.resolve()
        rel_fa = fasta_path.resolve().relative_to(workspace)
        rel_out = out_dir.resolve().relative_to(workspace)
        subcmd = "templatefree" if layout == "templatefree" else "draw"
        inner = [
            "python3",
            "r2dt.py",
            subcmd,
            f"/workspace/{rel_fa.as_posix()}",
            f"/workspace/{rel_out.as_posix()}",
            "--quiet",
        ]
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
    def _local_draw_cmd(fasta_path: Path, out_dir: Path, layout: str) -> List[str]:
        subcmd = "templatefree" if layout == "templatefree" else "draw"
        return [
            "python3",
            "r2dt.py",
            subcmd,
            str(fasta_path),
            str(out_dir),
            "--quiet",
        ]

    def _build_rscape_cmd(self, rscape_dir: Path) -> List[str]:
        """Build R-scape --cacofold with cwd = rscape_dir (in- or out-of-docker)."""
        if in_docker():
            return ["R-scape", "--cacofold", "alignment.stk"]
        workspace = self.catalog.workspace.resolve()
        rel = rscape_dir.resolve().relative_to(workspace)
        return [
            "docker",
            "run",
            "--rm",
            "-v",
            f"{workspace}:/workspace",
            "-w",
            f"/workspace/{rel.as_posix()}",
            self.docker_image,
            "R-scape",
            "--cacofold",
            "alignment.stk",
        ]

    def _build_stockholm_cmd(
        self,
        stockholm_path: Path,
        out_dir: Path,
        stitch: bool,
    ) -> List[str]:
        """Build r2dt.py stockholm (in- or out-of-docker)."""
        if in_docker():
            return self._local_stockholm_cmd(stockholm_path, out_dir, stitch)
        workspace = self.catalog.workspace.resolve()
        repo = self.repo_root.resolve()
        rel_sto = stockholm_path.resolve().relative_to(workspace)
        rel_out = out_dir.resolve().relative_to(workspace)
        inner = [
            "python3",
            "r2dt.py",
            "stockholm",
            f"/workspace/{rel_sto.as_posix()}",
            f"/workspace/{rel_out.as_posix()}",
            "--quiet",
            "--stitch" if stitch else "--no-stitch",
        ]
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
    def _local_stockholm_cmd(
        stockholm_path: Path, out_dir: Path, stitch: bool
    ) -> List[str]:
        return [
            "python3",
            "r2dt.py",
            "stockholm",
            str(stockholm_path),
            str(out_dir),
            "--quiet",
            "--stitch" if stitch else "--no-stitch",
        ]

    def _build_compare_cmd(  # pylint: disable=too-many-arguments,too-many-positional-arguments
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
            return self._local_compare_cmd(
                ref_path, model_path, out_dir, chains, model_chains, mode, basepairs
            )
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

    def _build_pdb_cmd(  # pylint: disable=too-many-arguments,too-many-positional-arguments
        self,
        structure_path: Path,
        out_dir: Path,
        chain: str,
        mode: str,
        basepairs: str,
    ) -> List[str]:
        """Build r2dt.py pdb_2d_3d for a single-structure interactive viewer."""
        if in_docker():
            return self._local_pdb_2d3d_cmd(
                structure_path, out_dir, chain, mode, basepairs
            )
        workspace = self.catalog.workspace.resolve()
        repo = self.repo_root.resolve()
        rel_struct = structure_path.resolve().relative_to(workspace)
        rel_out = out_dir.resolve().relative_to(workspace)
        inner = [
            "python3",
            "r2dt.py",
            "pdb_2d_3d",
            f"/workspace/{rel_struct.as_posix()}",
            f"/workspace/{rel_out.as_posix()}",
            "--mode",
            mode,
            "--basepairs",
            basepairs,
            "--quiet",
        ]
        if chain:
            inner.extend(["--chain", chain])
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
    def _local_compare_cmd(  # pylint: disable=too-many-arguments,too-many-positional-arguments
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

    @staticmethod
    def _local_pdb_2d3d_cmd(
        structure_path: Path,
        out_dir: Path,
        chain: str,
        mode: str,
        basepairs: str,
    ) -> List[str]:
        cmd = [
            "python3",
            "r2dt.py",
            "pdb_2d_3d",
            str(structure_path),
            str(out_dir),
            "--mode",
            mode,
            "--basepairs",
            basepairs,
            "--quiet",
        ]
        if chain:
            cmd.extend(["--chain", chain])
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
    """Materialise compare inputs, write meta, enqueue. Returns job meta."""
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

    chains_csv = ",".join(chains)
    model_chains_csv = ",".join(model_chains)
    digest = content_hash(
        ref_src, model_src, chains_csv, model_chains_csv, mode, basepairs
    )
    if not force:
        existing = catalog.find_by_content_hash(digest, mode="compare")
        if existing:
            return {
                "dedup": True,
                "job": existing,
                "message": "Identical comparison already cached",
            }

    job_id = new_job_id()
    inputs_dir = catalog.inputs_job_dir(job_id)
    inputs_dir.mkdir(parents=True, exist_ok=True)
    # Compare path needs mmCIF reference (FR3D multichain reader) — same as
    # CASP official PDB refs, which are converted before --compare.
    ref_original = ref_src.name
    ref_dest = ensure_mmcif(ref_src, inputs_dir, label=structure_stem(ref_src))
    model_dest = inputs_dir / model_src.name
    if model_src.resolve() != model_dest.resolve():
        shutil.copy2(model_src, model_dest)

    display = (
        label.strip() or f"{structure_stem(ref_src)} vs {structure_stem(model_src)}"
    )
    meta: Dict[str, Any] = {
        "id": job_id,
        "mode": "compare",
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
            "ref_original_name": ref_original,
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


def create_pdb_job_from_upload(  # pylint: disable=too-many-arguments,too-many-locals
    catalog: Catalog,
    runner: JobRunner,
    *,
    upload_id: str,
    chain: str,
    mode: str = "auto",
    basepairs: str = "fr3d",
    label: str = "",
    notes: str = "",
    force: bool = False,
) -> Dict[str, Any]:
    """Materialise a single-structure 2D+3D job and enqueue it."""
    chain = (chain or "").strip()
    if not chain:
        raise ValueError("Select one RNA chain")
    if not upload_id:
        raise ValueError("Missing upload id")

    src = _resolve_upload(catalog, upload_id)
    digest = pdb_content_hash(src, chain, mode, basepairs)
    if not force:
        existing = catalog.find_by_content_hash(digest, mode="pdb")
        if existing:
            return {
                "dedup": True,
                "job": existing,
                "message": "Identical 2D+3D job already cached",
            }

    job_id = new_job_id()
    inputs_dir = catalog.inputs_job_dir(job_id)
    inputs_dir.mkdir(parents=True, exist_ok=True)
    dest = inputs_dir / src.name
    if src.resolve() != dest.resolve():
        shutil.copy2(src, dest)

    stem = structure_stem(src)
    display = label.strip() or f"{stem} chain {chain}"
    meta: Dict[str, Any] = {
        "id": job_id,
        "mode": "pdb",
        "label": display,
        "notes": notes.strip(),
        "created": utc_now(),
        "status": "queued",
        "params": {
            "chain": chain,
            "mode": mode,
            "basepairs": basepairs,
        },
        "inputs": {
            "structure_name": dest.name,
            "structure_format": normalize_suffix(dest),
            "content_hash": digest,
        },
        "metrics": {},
        "viewer_url": f"/jobs/{job_id}/viewer/",
    }
    catalog.write_meta(job_id, meta)
    catalog.append_log(job_id, "Queued")
    runner.enqueue(job_id)
    return {"dedup": False, "job": meta}


def find_draw_svg(job_dir: Path) -> Optional[Path]:
    """Return the best ``*.colored.svg`` under a draw job directory."""
    job_dir = Path(job_dir)
    if preferred := sorted(job_dir.glob("**/results/svg/*.colored.svg")):
        return preferred[0]
    fallback = sorted(job_dir.glob("**/*.colored.svg"))
    return fallback[0] if fallback else None


def draw_content_hash(record: FastaRecord, layout: str) -> str:
    """Stable hash for one draw sequence + layout."""
    digest = hashlib.sha256()
    digest.update(b"draw\0")
    digest.update(layout.encode("utf-8"))
    digest.update(b"\0")
    digest.update(record.seq_id.encode("utf-8"))
    digest.update(b"\0")
    digest.update(record.sequence.encode("utf-8"))
    digest.update(b"\0")
    if record.structure:
        digest.update(record.structure.encode("utf-8"))
    return "sha256:" + digest.hexdigest()


def create_draw_jobs_from_fasta(  # pylint: disable=too-many-arguments,too-many-locals
    catalog: Catalog,
    runner: JobRunner,
    *,
    fasta_text: str,
    layout: str = "auto",
    label: str = "",
    notes: str = "",
    force: bool = False,
) -> Dict[str, Any]:
    """Split FASTA into N draw jobs (one sequence each) and enqueue them."""
    layout = (layout or "auto").strip().lower()
    if layout not in ("auto", "templatefree"):
        raise ValueError("layout must be 'auto' or 'templatefree'")

    records = parse_fasta_records(fasta_text)
    if layout == "templatefree":
        missing = [r.seq_id for r in records if not r.structure]
        if missing:
            raise ValueError(
                "templatefree requires a third-line structure for each sequence "
                f"(missing for: {', '.join(missing[:5])}"
                + ("…" if len(missing) > 5 else "")
                + ")"
            )

    batch_id = uuid.uuid4().hex[:10]
    label_prefix = label.strip()
    created: List[Dict[str, Any]] = []
    deduped: List[Dict[str, Any]] = []

    for record in records:
        digest = draw_content_hash(record, layout)
        if not force:
            existing = catalog.find_by_content_hash(digest, mode="draw")
            if existing:
                deduped.append(existing)
                continue

        job_id = new_job_id()
        inputs_dir = catalog.inputs_job_dir(job_id)
        inputs_dir.mkdir(parents=True, exist_ok=True)
        fasta_name = f"{record.seq_id}.fasta"
        fasta_path = inputs_dir / fasta_name
        fasta_path.write_text(record.to_fasta_text(), encoding="utf-8")

        display = label_prefix or record.seq_id
        if label_prefix and len(records) > 1:
            display = f"{label_prefix} · {record.seq_id}"

        meta: Dict[str, Any] = {
            "id": job_id,
            "mode": "draw",
            "label": display,
            "notes": notes.strip(),
            "created": utc_now(),
            "status": "queued",
            "batch_id": batch_id,
            "params": {
                "layout": layout,
                "seq_id": record.seq_id,
                "length": record.length,
                "has_structure": bool(record.structure),
            },
            "inputs": {
                "fasta_name": fasta_name,
                "content_hash": digest,
            },
            "metrics": {},
            "results_url": f"/jobs/{job_id}/results/",
        }
        catalog.write_meta(job_id, meta)
        catalog.append_log(job_id, "Queued")
        runner.enqueue(job_id)
        created.append(meta)

    if not created and not deduped:
        raise ValueError("No jobs created")

    return {
        "batch_id": batch_id,
        "jobs": created + deduped,
        "created": created,
        "deduped": deduped,
        "message": (
            f"Queued {len(created)} job(s)"
            + (f", {len(deduped)} already cached" if deduped else "")
        ),
    }
