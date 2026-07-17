#!/usr/bin/env python3
"""Fetch a CASP season's RNA references + top-N predicted models and emit a
manifest for ``casp_batch.py``.

For each RNA target with a public experimental reference (per
``casp_config.SEASONS[season]["target_pdb"]``), this:
  1. resolves the reference — either downloaded from RCSB, or (CASP15) taken
     from CASP's own official reference archive, one per *state* for targets
     solved in multiple conformations. A state may itself be an *ensemble*
     (R1149/R1156: alternate models fit to a low-resolution map) — every
     submitted model is checked against every ensemble member and the
     best-fitting one wins, per-model (see ``_rank_reference``),
  2. downloads predictions/RNA/<target>.tar.gz and extracts the top-N models
     (ranked by TM-score, from ``casp_rank.py``'s output),
  3. resolves which reference/model RNA chains co-index (equal sequence, same
     order) via BioPython sequence matching,
  4. writes ``targets.json`` for ``casp_batch.py`` — one manifest entry per
     target, or per state for multi-state targets (``<target>-<state>``, all
     states sharing the same ranked model list). Each model entry carries its
     own resolved ``reference``/``ref_chains`` (which can differ per model for
     ensemble states); the manifest entry also carries a target-level default
     (the most commonly picked reference) for display.

Targets without a public reference, or whose chains don't co-index, are skipped
with a logged reason. Downloads are cached.
"""
import argparse
import json
import sys
import tarfile
import warnings
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
# pylint: disable=wrong-import-position
from casp_config import download, get_season  # noqa: E402

RCSB = "https://files.rcsb.org/download/{pdb}.cif"
_STD = {"A": "A", "U": "U", "G": "G", "C": "C"}


def download_and_extract_archive(url, cache_tgz, extract_dir):
    """Idempotent download + flat extraction of an official CASP reference
    archive, including any nested per-target sub-tarballs (e.g. R1149.tgz,
    an ensemble shipped as a tarball-within-the-tarball)."""
    if extract_dir.exists() and any(extract_dir.iterdir()):
        return extract_dir
    download(url, cache_tgz)
    extract_dir.mkdir(parents=True, exist_ok=True)

    def _flat_extract(tar_path):
        with tarfile.open(tar_path) as tf:
            for member in tf.getmembers():
                if not member.isfile():
                    continue
                member.name = Path(member.name).name  # drop any archive-internal dirs
                tf.extract(member, extract_dir)

    _flat_extract(cache_tgz)
    for nested in sorted(extract_dir.glob("*.tgz")) + sorted(
        extract_dir.glob("*.tar.gz")
    ):
        _flat_extract(nested)
        nested.unlink()
    return extract_dir


def ensure_cif(path, cif_dir, label=None):
    """Return a ``.cif`` version of ``path`` (BioPython conversion, cached by
    ``label`` — or ``path``'s own stem if not given — in ``cif_dir``).
    ``r2dt.py``'s ``--chains``/``--compare`` multichain path requires mmCIF
    input; CASP's official reference archive ships plain PDB files, unlike
    RCSB (downloaded as ``.cif`` directly).

    ``label`` also becomes the viewer's displayed structure id (r2dt.py
    derives it from the reference filename), so callers pass something
    short and recognisable (e.g. the target id) instead of the archive's
    often long, underscore-heavy original filename.
    """
    path = Path(path)
    if path.suffix.lower() == ".cif" and label is None:
        return path
    out = cif_dir / f"{label or path.stem}.cif"
    if out.exists() and out.stat().st_size > 0:
        return out
    warnings.simplefilter("ignore")
    from Bio.PDB import (  # pylint: disable=import-outside-toplevel
        MMCIFIO,
        MMCIFParser,
        PDBParser,
    )

    cif_dir.mkdir(parents=True, exist_ok=True)
    parser = (
        MMCIFParser(QUIET=True)
        if path.suffix.lower() == ".cif"
        else PDBParser(QUIET=True)
    )
    structure = parser.get_structure(label or path.stem, str(path))
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(str(out))
    return out


def rna_chain_seqs(structure_file):
    """Return {chain_id: sequence} for RNA chains (standard bases; modified→'N')."""
    warnings.simplefilter("ignore")
    from Bio.PDB import (  # pylint: disable=import-outside-toplevel
        MMCIFParser,
        PDBParser,
    )

    parser = (
        MMCIFParser(QUIET=True)
        if str(structure_file).lower().endswith(".cif")
        else PDBParser(QUIET=True)
    )
    st = parser.get_structure("s", str(structure_file))
    model = next(iter(st))
    out = {}
    for ch in model:
        seq, std = [], 0
        for res in ch:
            if res.id[0] != " ":
                continue
            nm = res.resname.strip()
            if nm in _STD:
                seq.append(_STD[nm])
                std += 1
            elif res.has_id("P") or res.has_id("C1'"):
                seq.append("N")  # modified / non-standard nucleotide
        if std >= 3:  # an RNA chain (excludes protein/DNA/ligand)
            out[ch.id] = "".join(seq)
    return out


def _ident(a, b):
    if len(a) != len(b):
        return 0.0
    return sum(x == y or x == "N" or y == "N" for x, y in zip(a, b)) / len(a)


def _window(ref_seq, model_seq, thresh=0.9):
    """If ``ref_seq`` matches a same-length contiguous window of ``model_seq``,
    return that window's start index; else None. Handles references with
    unresolved *terminal* residues (model longer than the resolved reference)."""
    n, m = len(ref_seq), len(model_seq)
    if n == 0 or n > m:
        return None
    best_off, best_id = None, thresh
    for off in range(0, m - n + 1):
        score = _ident(ref_seq, model_seq[off : off + n])
        if score >= best_id:
            best_id, best_off = score, off
    return best_off


def _align_indices(ref_seq, model_seq, thresh=0.9):
    """Global-align the resolved reference chain to the model and return the list
    of model indices aligned to reference residues, in order (length ==
    len(ref_seq)), or None. Handles references with *internal* unresolved residues
    that the contiguous-window match can't. Every reference residue must map to a
    model residue and match at >= ``thresh`` identity (N is a wildcard)."""
    from Bio.Align import PairwiseAligner  # pylint: disable=import-outside-toplevel

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -6
    aligner.extend_gap_score = -0.5
    try:
        aln = aligner.align(ref_seq, model_seq)[0]
    except (IndexError, ValueError):
        return None
    ref_blocks, mod_blocks = aln.aligned  # equal-length aligned segments
    mapping = {}  # ref_idx -> model_idx
    matches = 0
    for (r0, r1), (m0, _m1) in zip(ref_blocks, mod_blocks):
        for d in range(r1 - r0):
            ri, mi = r0 + d, m0 + d
            mapping[ri] = mi
            a, b = ref_seq[ri], model_seq[mi]
            if a == b or a == "N" or b == "N":
                matches += 1
    if len(mapping) != len(ref_seq):  # some reference residue unmatched in model
        return None
    if matches / len(ref_seq) < thresh:
        return None
    return [mapping[i] for i in range(len(ref_seq))]


def resolve_chains(ref_file, model_file):
    """Resolve co-indexing chains, or None.

    Returns ``{"ref_chains": [...], "model_chains": [...], "trim": {chain: (start, count)}}``.
    The model reproduces the target sequence (usually one RNA chain); match it to
    the reference's RNA chain(s). ``trim`` (present only when needed) says to keep
    residues ``[start:start+count]`` of a model chain so it co-indexes with a
    reference that has unresolved terminal residues."""
    ref = rna_chain_seqs(ref_file)
    mod = rna_chain_seqs(model_file)
    if not ref or not mod:
        return None
    model_chains = list(mod.keys())
    model_seq = "".join(mod[c] for c in model_chains)

    # 1) single reference chain == whole model (exact length).
    for c, s in ref.items():
        if _ident(s, model_seq) >= 0.9:
            return {"ref_chains": [c], "model_chains": model_chains, "trim": {}}

    # 2) single-chain model, reference is a subset of it → trim the model to the
    #    reference's resolved residues. Try a contiguous window first (unresolved
    #    termini), then a gapped alignment (internal unresolved residues).
    if len(model_chains) == 1:
        mc = model_chains[0]
        for c, s in ref.items():
            off = _window(s, mod[mc])
            if off is not None:
                return {
                    "ref_chains": [c],
                    "model_chains": model_chains,
                    "trim": {mc: list(range(off, off + len(s)))},
                }
        for c, s in ref.items():
            idx = _align_indices(s, mod[mc])
            if idx is not None:
                return {
                    "ref_chains": [c],
                    "model_chains": model_chains,
                    "trim": {mc: idx},
                }

    # 3) per-model-chain greedy exact match (multimers).
    chosen, used = [], set()
    for mc in model_chains:
        pick = next(
            (c for c, s in ref.items() if c not in used and _ident(s, mod[mc]) >= 0.9),
            None,
        )
        if pick is None:
            return None
        chosen.append(pick)
        used.add(pick)
    return {"ref_chains": chosen, "model_chains": model_chains, "trim": {}}


def _c1_coords(structure_file, chain_id, indices=None):
    """C1' coordinates for one chain's RNA residues, in order. ``indices``
    (0-based, in the same sense as ``resolve_chains``' ``trim``) restricts to
    a subset — used to line up a trimmed model against the reference."""
    from Bio.PDB import (  # pylint: disable=import-outside-toplevel
        MMCIFParser,
        PDBParser,
    )

    parser = (
        MMCIFParser(QUIET=True)
        if str(structure_file).lower().endswith(".cif")
        else PDBParser(QUIET=True)
    )
    st = parser.get_structure("s", str(structure_file))
    model = next(iter(st))
    if chain_id not in model:
        return []
    ordered = [r for r in model[chain_id] if r.id[0] == " " and r.has_id("C1'")]
    if indices is not None:
        ordered = [ordered[i] for i in indices if 0 <= i < len(ordered)]
    return [r["C1'"].coord for r in ordered]


def _quick_rmsd(ref_file, model_file, resolved):
    """Cheap coordinate RMSD (Kabsch fit on co-indexed C1' atoms), used only to
    rank *candidate* references — e.g. which ensemble member (R1149/R1156)
    best fits a given model. Not the official score: that comes from the real
    R2DT/FR3D Docker pipeline (``casp_batch.py``) for whichever candidate
    wins here. Returns None if the atom counts don't line up 1:1."""
    warnings.simplefilter("ignore")
    import numpy as np  # pylint: disable=import-outside-toplevel
    from Bio.SVDSuperimposer import (  # pylint: disable=import-outside-toplevel
        SVDSuperimposer,
    )

    ref_coords = [
        c for cid in resolved["ref_chains"] for c in _c1_coords(ref_file, cid)
    ]
    model_coords = [
        c
        for cid in resolved["model_chains"]
        for c in _c1_coords(model_file, cid, indices=resolved["trim"].get(cid))
    ]
    if len(ref_coords) != len(model_coords) or len(ref_coords) < 3:
        return None
    sup = SVDSuperimposer()
    sup.set(np.array(ref_coords), np.array(model_coords))
    sup.run()
    return sup.get_rms()


def _rank_reference(ref_candidates, model_file):
    """Resolve ``model_file`` against every candidate reference and return the
    best-fitting ``(ref_file, resolved)``, or None if nothing resolves.

    A single-candidate list (every non-ensemble target) just resolves that one
    candidate directly — no ranking, no quick-RMSD cost. Multiple candidates
    (an ensemble) are ranked by :func:`_quick_rmsd`, lowest wins."""
    if len(ref_candidates) == 1:
        ref = ref_candidates[0]
        if not ref.exists():
            return None
        resolved = resolve_chains(str(ref), str(model_file))
        return (ref, resolved) if resolved is not None else None

    best = None
    for ref in ref_candidates:
        if not ref.exists():
            continue
        resolved = resolve_chains(str(ref), str(model_file))
        if resolved is None:
            continue
        score = _quick_rmsd(str(ref), str(model_file), resolved)
        if score is None:
            continue
        if best is None or score < best[0]:  # pylint: disable=unsubscriptable-object
            best = (score, ref, resolved)
    return (best[1], best[2]) if best else None


def write_trimmed_model(model_file, trim, out_file):
    """Write a model PDB keeping only the RNA residues at the given indices of each
    trimmed chain (``trim[chain] = [indices]``, 0-based within that chain's ordered
    RNA residues); other chains are dropped. Returns out_file."""
    warnings.simplefilter("ignore")
    from Bio.PDB import (  # pylint: disable=import-outside-toplevel
        PDBIO,
        PDBParser,
        Select,
    )

    st = PDBParser(QUIET=True).get_structure("m", str(model_file))
    model = next(iter(st))
    keep = set()
    for ch in model:
        if ch.id not in trim:
            continue
        ordered = [r for r in ch if r.id[0] == " "]
        for i in trim[ch.id]:
            if 0 <= i < len(ordered):
                keep.add((ch.id, ordered[i].id))

    class _Keep(Select):
        """Select only the (chain, residue id) pairs named in ``keep``."""

        def accept_residue(self, residue):  # noqa: N802
            """Bio.PDB.Select hook: keep only residues named in ``keep``."""
            return (residue.get_parent().id, residue.id) in keep

    io = PDBIO()
    io.set_structure(st)
    io.save(str(out_file), _Keep())
    return out_file


def extract_models(tar_path, model_names, out_dir):
    """Extract the named CASP model files from the tarball as <name>.pdb."""
    out_dir.mkdir(parents=True, exist_ok=True)
    wanted = {n: out_dir / f"{n}.pdb" for n in model_names}
    have = {n: p for n, p in wanted.items() if p.exists() and p.stat().st_size > 0}
    if len(have) == len(wanted):
        return wanted
    with tarfile.open(tar_path) as tf:
        members = {Path(m.name).name: m for m in tf.getmembers() if m.isfile()}
        for name, dest in wanted.items():
            m = members.get(name)
            if m is None:
                continue
            with tf.extractfile(m) as fh:
                dest.write_bytes(fh.read())
    return {n: p for n, p in wanted.items() if p.exists()}


def _resolve_ref_candidates(source, manifest_key, data, official_refs_dir):
    """Return (ref_candidates, ref_label) for a target-state's ``source``
    config, or (None, None) on an unrecoverable failure (logged here)."""
    if source["kind"] == "rcsb":
        pdb = source["pdb"]
        try:
            ref = download(RCSB.format(pdb=pdb), data / "refs" / f"{pdb}.cif")
        except Exception as e:  # pylint: disable=broad-except
            print(f"[skip ] {manifest_key}: reference download failed: {e}")
            return None, None
        return [ref], pdb

    # official archive
    if official_refs_dir is None:
        print(f"[skip ] {manifest_key}: no official reference archive available")
        return None, None
    candidates = [official_refs_dir / fn for fn in source["files"]]
    label = (
        candidates[0].name
        if len(candidates) == 1
        else f"{len(candidates)}-model ensemble"
    )
    return candidates, label


# pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals
def resolve_target(
    base_target, states, top, data, rankings, results, pred_url, official_refs_dir=None
):
    """Resolve every state of one base target against its ranked models. Appends
    ``(manifest_key, entry)`` pairs to ``results``; logs and returns on any
    unrecoverable failure (missing rankings/tarball/models)."""
    if base_target not in rankings:
        print(f"[skip ] {base_target}: no rankings")
        return
    try:
        tar = download(
            pred_url.format(target=base_target), data / "tars" / f"{base_target}.tar.gz"
        )
    except Exception as e:  # pylint: disable=broad-except
        print(f"[skip ] {base_target}: predictions download failed: {e}")
        return
    top_models = rankings[base_target][:top]
    models = extract_models(
        tar, [m["model"] for m in top_models], data / "models" / base_target
    )
    if not models:
        print(f"[skip ] {base_target}: no model files extracted")
        return

    for state, source in states:
        manifest_key = f"{base_target}-{state}" if state else base_target
        ref_candidates, ref_label = _resolve_ref_candidates(
            source, manifest_key, data, official_refs_dir
        )
        if ref_candidates is None:
            continue
        is_ensemble = len(ref_candidates) > 1

        entry_models = []
        chosen_refs = []
        n_trim = 0
        for m in top_models:
            if m["model"] not in models:
                continue
            mfile = models[m["model"]]
            picked = _rank_reference(ref_candidates, mfile)
            if picked is None:
                continue
            ref, resolved = picked
            # r2dt.py's --chains/--compare path requires mmCIF; RCSB refs are
            # already .cif (and keep their short PDB-code id unchanged).
            # Official-archive refs are plain .pdb with often long,
            # underscore-heavy filenames -- relabel to something short and
            # recognisable, since r2dt.py derives the viewer's displayed
            # structure id from this filename.
            if source["kind"] == "official":
                cif_label = (
                    manifest_key
                    if not is_ensemble
                    else f"{manifest_key}-{Path(ref).stem}"
                )
                ref = ensure_cif(ref, data / "refs_cif", label=cif_label)
            trim = resolved["trim"]
            if trim:
                # Keep the same basename (in a subdir) so the model id stays clean
                # — a ".trim" in the stem would leak a dot into the viewer's pdbId.
                tdir = mfile.parent / "trimmed" / (state or ".")
                tdir.mkdir(parents=True, exist_ok=True)
                tfile = tdir / mfile.name
                try:
                    write_trimmed_model(mfile, trim, tfile)
                    mfile = tfile
                    n_trim += 1
                except Exception:  # pylint: disable=broad-except
                    continue  # can't trim → skip this model
            chosen_refs.append(ref)
            entry_models.append(
                {
                    "model": m["model"],
                    "file": str(mfile),
                    "tm": float(m["tm"]) if m["tm"] not in ("-", "") else None,
                    "rank": m["rank"],
                    "model_chains": ",".join(resolved["model_chains"]),
                    # Authoritative for scoring (casp_batch.py); can differ per
                    # model when the state is an ensemble (each model may best
                    # match a different member).
                    "reference": str(ref),
                    "ref_chains": ",".join(resolved["ref_chains"]),
                }
            )

        if not entry_models:
            print(
                f"[skip ] {manifest_key}: no model co-indexes with any reference candidate"
            )
            continue

        # Target-level reference/ref_chains = the most commonly picked
        # candidate, for display only — each model's own fields above are
        # what casp_batch.py actually scores against.
        default_ref = str(Counter(chosen_refs).most_common(1)[0][0])
        default_ref_chains = next(
            m["ref_chains"] for m in entry_models if m["reference"] == default_ref
        )
        results.append(
            (
                manifest_key,
                {
                    "reference": default_ref,
                    "ref_chains": default_ref_chains,
                    "reference_pdb": ref_label,
                    "models": entry_models,
                },
            )
        )
        ensemble_note = (
            f", {len(set(chosen_refs))}/{len(ref_candidates)} ensemble members used"
            if is_ensemble
            else ""
        )
        print(
            f"[ok   ] {manifest_key}: ref {ref_label} · "
            f"{len(entry_models)} models ({n_trim} trimmed){ensemble_note}"
        )


def main():
    """CLI entry point: fetch one season's references + top-N models, write targets.json."""
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--season", required=True, choices=["casp15", "casp16"])
    ap.add_argument("--rankings", help="default temp/<season>/top5_by_target.json")
    ap.add_argument("--data", help="default temp/<season>/data")
    ap.add_argument("--out", help="default temp/<season>/targets.json")
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("--targets", nargs="*", help="subset of base target ids")
    args = ap.parse_args()

    cfg = get_season(args.season)
    pred_url = cfg["pred_url"]

    rankings_path = (
        Path(args.rankings)
        if args.rankings
        else Path("temp") / args.season / "top5_by_target.json"
    )
    data = Path(args.data) if args.data else Path("temp") / args.season / "data"
    out = Path(args.out) if args.out else Path("temp") / args.season / "targets.json"

    rankings = json.loads(rankings_path.read_text())
    target_pdb = cfg["target_pdb"]
    base_targets = args.targets or list(target_pdb)

    official_refs_dir = None
    if cfg.get("official_archive_url"):
        print(f"Fetching official reference archive ({cfg['official_archive_url']})...")
        official_refs_dir = download_and_extract_archive(
            cfg["official_archive_url"],
            data / "official_archive.tgz",
            data / "official_refs",
        )

    results = []
    for base_target in base_targets:
        states = target_pdb.get(base_target)
        if not states:
            print(f"[skip ] {base_target}: no reference PDB configured")
            continue
        resolve_target(
            base_target,
            states,
            args.top,
            data,
            rankings,
            results,
            pred_url,
            official_refs_dir=official_refs_dir,
        )

    manifest = {"targets": dict(results)}
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"\nWrote {out} — {len(manifest['targets'])} manifest entries")


if __name__ == "__main__":
    main()
