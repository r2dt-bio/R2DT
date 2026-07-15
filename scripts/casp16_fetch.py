#!/usr/bin/env python3
"""Fetch CASP16 RNA references + top-5 predicted models and emit a manifest.

For each RNA target that has a public experimental reference (PDB code in the
CASP16 target-list descriptions), this:
  1. downloads the reference mmCIF from RCSB,
  2. downloads predictions/RNA/<target>.tar.gz and extracts the top-5 models
     (ranked by TM-score in the CASP results table),
  3. resolves which reference/model RNA chains co-index (equal sequence, same
     order) via BioPython sequence matching,
  4. writes ``targets.json`` for ``casp16_batch.py``.

Targets without a public reference, or whose chains don't co-index, are skipped
with a logged reason. Downloads are cached.
"""
import argparse
import json
import tarfile
import urllib.request
import warnings
from pathlib import Path

from Bio.Align import PairwiseAligner
from Bio.PDB import PDBIO, MMCIFParser, PDBParser, Select

warnings.simplefilter("ignore")

# Target -> reference PDB code, from the CASP16 target-list description column.
# Only targets with a public deposited structure are listed (the rest have no
# reference, so reference-vs-model is impossible).
TARGET_PDB = {
    "R1203": "8uo6",
    "R1205": "9cfn",
    "R1209": "9c2k",
    "R1211": "9dcf",
    "R1212": "9b0l",
    "R1242": "9ely",
    "R1251": "9mee",
    "R1260": "9cbu",
    "R1261": "9bzc",
    "R1262": "9bz1",
    "R1263": "8vqv",
    "R1264": "8vvj",
    "R1285": "9mcw",
    "R1286": "9j6y",
    "D1273": "9hio",
}
RCSB = "https://files.rcsb.org/download/{pdb}.cif"
PRED = (
    "https://predictioncenter.org/download_area/CASP16/predictions/RNA/{target}.tar.gz"
)
_STD = {"A": "A", "U": "U", "G": "G", "C": "C"}


def download(url, dest):
    """Fetch ``url`` to ``dest`` if not already cached, and return ``dest``."""
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(url, dest)
    return dest


def rna_chain_seqs(structure_file):
    """Return {chain_id: sequence} for RNA chains (standard bases; modified→'N')."""
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


def write_trimmed_model(model_file, trim, out_file):
    """Write a model PDB keeping only the RNA residues at the given indices of each
    trimmed chain (``trim[chain] = [indices]``, 0-based within that chain's ordered
    RNA residues); other chains are dropped. Returns out_file."""
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
        def accept_residue(self, residue):  # noqa: N802
            """Keep only the trimmed residues recorded in ``keep``."""
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


def main():  # pylint: disable=too-many-statements
    """Fetch references + top-N models for each target and write the manifest."""
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rankings", default="temp/casp16/top5_by_target.json")
    ap.add_argument("--data", default="temp/casp16/data")
    ap.add_argument("--out", default="temp/casp16/targets.json")
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("--targets", nargs="*", help="subset of target ids")
    args = ap.parse_args()

    rankings = json.loads(Path(args.rankings).read_text())
    data = Path(args.data)
    targets = args.targets or list(TARGET_PDB)

    manifest = {"targets": {}}
    for target in targets:
        pdb = TARGET_PDB.get(target)
        if not pdb or target not in rankings:
            print(f"[skip ] {target}: no reference PDB or no rankings")
            continue
        try:
            ref = download(RCSB.format(pdb=pdb), data / "refs" / f"{pdb}.cif")
            tar = download(
                PRED.format(target=target), data / "tars" / f"{target}.tar.gz"
            )
        except Exception as e:  # pylint: disable=broad-except
            print(f"[skip ] {target}: download failed: {e}")
            continue
        top = rankings[target][: args.top]
        models = extract_models(
            tar, [m["model"] for m in top], data / "models" / target
        )
        if not models:
            print(f"[skip ] {target}: no model files extracted")
            continue

        # Resolve each available model independently (chains + per-model trim), so
        # a missing/oddly-shaped model doesn't sink the whole target.
        entry_models = []
        ref_chains = None
        n_trim = 0
        for m in top:
            if m["model"] not in models:
                continue
            mfile = models[m["model"]]
            resolved = resolve_chains(str(ref), str(mfile))
            if resolved is None:
                continue
            if ref_chains is None:
                ref_chains = resolved["ref_chains"]
            trim = resolved["trim"]
            if trim:
                # Keep the same basename (in a subdir) so the model id stays clean
                # — a ".trim" in the stem would leak a dot into the viewer's pdbId.
                tdir = mfile.parent / "trimmed"
                tdir.mkdir(exist_ok=True)
                tfile = tdir / mfile.name
                try:
                    write_trimmed_model(mfile, trim, tfile)
                    mfile = tfile
                    n_trim += 1
                except Exception:  # pylint: disable=broad-except
                    continue  # can't trim → skip this model
            entry_models.append(
                {
                    "model": m["model"],
                    "file": str(mfile),
                    "tm": float(m["tm"]) if m["tm"] not in ("-", "") else None,
                    "rank": m["rank"],
                    "model_chains": ",".join(resolved["model_chains"]),
                }
            )

        if not entry_models:
            print(f"[skip ] {target}: no model co-indexes with the reference")
            continue
        manifest["targets"][target] = {
            "reference": str(ref),
            "ref_chains": ",".join(ref_chains),
            "reference_pdb": pdb,
            "models": entry_models,
        }
        print(
            f"[ok   ] {target}: ref {pdb} chains {ref_chains} · "
            f"{len(entry_models)} models ({n_trim} trimmed)"
        )

    Path(args.out).write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"\nWrote {args.out} — {len(manifest['targets'])} targets")


if __name__ == "__main__":
    main()
