#!/usr/bin/env python3
"""Batch-generate CASP reference/model compare pages and aggregate metrics.

Season-agnostic: reads a ``targets.json`` manifest (as produced by
``casp_fetch.py`` for any season), runs the (Dockerised)
``r2dt.py pdb --compare --model`` pipeline for each target's top-N models, and
collects each run's ``viewer/metrics.json`` into a single ``<site>/results.json``
that drives the dashboard (``casp_dashboard.py``).

The pipeline needs the r2dt Docker image (Infernal, Traveler, FR3D, …), so each
pair is generated with ``docker run -v <repo>:/rna/r2dt``.  Paths in the manifest
are relative to the repo root (the mount point).

Manifest shape (JSON)::

    {
      "targets": {
        "R1108": {
          "reference": "temp/7qr3.cif",
          "ref_chains": "C",
          "models": [
            {"model": "R1108TS128_3", "file": "temp/R1108TS128_3.pdb",
             "tm": 0.62, "rank": 1, "model_chains": null}
          ]
        }
      }
    }

Idempotent (skips pairs whose ``viewer/index.html`` exists unless ``--force``)
and fault-isolated (a failing pair is recorded and the batch continues).
"""
import argparse
import json
import subprocess
from pathlib import Path

DOCKER_IMAGE = "rnacentral/r2dt:latest"


# pylint: disable=too-many-arguments,too-many-positional-arguments
def run_pair(repo, image, out_rel, reference, ref_chains, model_file, model_chains):
    """Run the compare pipeline for one pair inside Docker. Returns CompletedProcess."""
    cmd = [
        "docker",
        "run",
        "--rm",
        "-v",
        f"{repo}:/rna/r2dt",
        image,
        "python3",
        "r2dt.py",
        "pdb",
        reference,
        out_rel,
        "--compare",
        "--model",
        model_file,
        "--quiet",
    ]
    if ref_chains:
        cmd += ["--chains", str(ref_chains)]
    if model_chains:
        cmd += ["--model-chains", str(model_chains)]
    return subprocess.run(
        cmd, cwd=str(repo), capture_output=True, text=True, check=False
    )


def collect_metrics(viewer_dir):
    """Read a pair's viewer/metrics.json into flat dashboard fields (or {})."""
    mfile = viewer_dir / "metrics.json"
    if not mfile.exists():
        return {}
    md = json.loads(mfile.read_text())
    inf = md.get("inf", {})
    return {
        "inf": {k: (inf.get(k) or {}).get("inf") for k in ("wc", "nwc", "all")},
        "rmsd": md.get("superpose_rmsd"),
        "n_superposed": md.get("superpose_n_atoms"),
        "diff": md.get("diff"),
        "chains": md.get("chains"),
    }


def main():
    """Run the manifest-driven batch and write ``<site>/results.json``."""
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("manifest", help="targets.json manifest")
    ap.add_argument("--site", default="site", help="output site directory")
    ap.add_argument("--image", default=DOCKER_IMAGE, help="r2dt docker image")
    ap.add_argument("--top", type=int, default=5, help="models per target")
    ap.add_argument("--force", action="store_true", help="regenerate cached pairs")
    args = ap.parse_args()

    repo = Path.cwd().resolve()
    manifest = json.loads(Path(args.manifest).read_text())
    site = Path(args.site)
    site.mkdir(parents=True, exist_ok=True)

    results = []
    for target, cfg in manifest["targets"].items():
        default_reference = cfg["reference"]
        default_ref_chains = cfg.get("ref_chains")
        models = sorted(cfg["models"], key=lambda m: m.get("rank") or 999)[: args.top]
        for idx, m in enumerate(models, 1):
            rank = m.get("rank") or idx
            model_id = m["model"]
            # A model may carry its own reference/ref_chains, overriding the
            # target-level default — e.g. an ensemble state (R1149/R1156)
            # where different models best match different ensemble members.
            reference = m.get("reference") or default_reference
            ref_chains = m.get("ref_chains") or default_ref_chains
            out_rel = f"{args.site}/{target}/{rank:02d}-{model_id}"
            viewer = repo / out_rel / "viewer"
            index = viewer / "index.html"

            status, err = "ok", None
            if index.exists() and not args.force:
                status = "cached"
            else:
                proc = run_pair(
                    repo,
                    args.image,
                    out_rel,
                    reference,
                    ref_chains,
                    m["file"],
                    m.get("model_chains"),
                )
                if proc.returncode != 0 or not index.exists():
                    status = "failed"
                    err = ((proc.stderr or "") + (proc.stdout or "")).strip()[-800:]

            entry = {
                "target": target,
                "reference": reference,
                "model": model_id,
                "rank": rank,
                "tm": m.get("tm"),
                "status": status,
                "page": f"{target}/{rank:02d}-{model_id}/viewer/index.html",
            }
            entry.update(collect_metrics(viewer))
            if err:
                entry["error"] = err
            results.append(entry)
            print(f"[{status:6}] {target} rank {rank:>2} {model_id}")

    out = site / "results.json"
    out.write_text(json.dumps({"results": results}, indent=2) + "\n")
    n_ok = sum(1 for r in results if r["status"] in ("ok", "cached"))
    print(f"\nWrote {out} — {n_ok}/{len(results)} pairs available")


if __name__ == "__main__":
    main()
