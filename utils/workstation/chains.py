"""RNA chain auto-detection for uploaded structures."""

from __future__ import annotations

import gzip
import shutil
import warnings
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from utils import fr3d as fr3d_utils
from utils.pdb_fetch import DecompressedStructureFile, open_structure_file

ALLOWED_SUFFIXES = {".cif", ".pdb", ".cif.gz", ".pdb.gz"}
IDENT_THRESHOLD = 0.9


def normalize_suffix(path: Path) -> str:
    """Return a lowercase structure suffix, including ``.gz`` when present."""
    name = path.name.lower()
    for suffix in (".cif.gz", ".pdb.gz", ".cif", ".pdb"):
        if name.endswith(suffix):
            return suffix
    return path.suffix.lower()


def structure_stem(path: Path) -> str:
    """Basename without ``.pdb`` / ``.cif`` / ``.gz`` suffixes."""
    name = path.name
    lower = name.lower()
    for suffix in (".cif.gz", ".pdb.gz", ".cif", ".pdb"):
        if lower.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def is_structure_filename(filename: str) -> bool:
    """True if filename looks like a PDB/mmCIF (optionally gzipped)."""
    return normalize_suffix(Path(filename)) in ALLOWED_SUFFIXES


def ensure_mmcif(path: Path, out_dir: Path, label: Optional[str] = None) -> Path:
    """Return an mmCIF path for ``path``, converting PDB when needed.

    Compare / multi-chain mode reads structures through FR3D's mmCIF reader
    (same constraint as ``r2dt.py pdb --compare`` and CASP ``ensure_cif``).
    Plain ``.cif`` is copied as-is; ``.pdb`` / gzipped inputs are written as
    ``{label}.cif`` under ``out_dir``.
    """
    # pylint: disable=import-outside-toplevel
    path = Path(path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    suffix = normalize_suffix(path)
    stem = label or structure_stem(path)
    out = out_dir / f"{stem}.cif"

    if suffix == ".cif":
        if path.resolve() != out.resolve():
            shutil.copy2(path, out)
        return out

    if suffix == ".cif.gz":
        with gzip.open(path, "rb") as src, out.open("wb") as dest:
            shutil.copyfileobj(src, dest)
        return out

    from Bio.PDB import MMCIFIO, PDBParser

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        parser = PDBParser(QUIET=True)
        if suffix == ".pdb.gz":
            with gzip.open(path, "rt") as handle:
                structure = parser.get_structure(stem, handle)
        else:
            structure = parser.get_structure(stem, str(path))
        io = MMCIFIO()
        io.set_structure(structure)
        io.save(str(out))
    return out


def sequence_preview(seq: str, head: int = 10, tail: int = 6) -> str:
    """Short 5′…3′ snippet for chain-picker labels."""
    text = seq or ""
    if len(text) <= head + tail + 1:
        return text
    return f"{text[:head]}…{text[-tail:]}"


def sequence_identity(left: str, right: str) -> float:
    """Fraction of aligned letters that match; ``N`` is a wildcard.

    Returns 0 when lengths differ or either string is empty.
    """
    if not left or not right or len(left) != len(right):
        return 0.0
    hits = sum(
        left_nt == right_nt or left_nt == "N" or right_nt == "N"
        for left_nt, right_nt in zip(left, right)
    )
    return hits / len(left)


def concatenated_length(
    details: Sequence[Mapping[str, Any]], chain_ids: Sequence[str]
) -> int:
    """Sum of ``length`` for ``chain_ids`` in that order."""
    by_id = {str(item["id"]): int(item.get("length") or 0) for item in details}
    total = 0
    for chain_id in chain_ids:
        key = str(chain_id)
        if key not in by_id:
            raise ValueError(f"Unknown chain id in details: {key}")
        total += by_id[key]
    return total


def list_rna_chains(structure_path: Path) -> Dict[str, Any]:
    """Return RNA chain ids and per-chain length/sequence for the first model.

    Prefers a one-pass residue scan (FR3D for mmCIF, BioPython for PDB). Falls
    back to ``fr3d.get_structure_info`` when that scan finds nothing.
    """
    path = Path(structure_path)
    suffix = normalize_suffix(path)
    details = rna_chain_details(path)
    fmt = "cif" if suffix.startswith(".cif") else "pdb"
    if details:
        chains = [str(item["id"]) for item in details]
        model_id = details[0].get("model")
    else:
        info = fr3d_utils.get_structure_info(str(path))
        models = info.get("models") or []
        chains_by_model = info.get("chains") or {}
        model_id = models[0] if models else None
        chains = (
            list(chains_by_model.get(model_id) or []) if model_id is not None else []
        )
    return {
        "filename": path.name,
        "format": fmt,
        "suffix": suffix,
        "model": model_id,
        "chains": chains,
        "chain_details": details,
        # PDB refs are auto-converted to mmCIF at job create time.
        "compare_ready": True,
        "needs_cif_conversion": fmt == "pdb",
    }


def rna_chain_details(structure_path: Path) -> List[Dict[str, Any]]:
    """Per-chain id, length, sequence, and author residue range (first model)."""
    path = Path(structure_path)
    suffix = normalize_suffix(path)
    try:
        if suffix.startswith(".cif"):
            details = _details_from_cif(path)
        else:
            details = _details_from_pdb(path)
    except Exception:  # pylint: disable=broad-exception-caught
        details = []
    if details:
        _mark_duplicate_sequences(details)
    return details


def assert_chains_known(
    requested: Sequence[str],
    available: Sequence[str],
    *,
    side: str,
) -> List[str]:
    """Return stripped chain ids, or raise ValueError if any are unknown."""
    cleaned = [str(chain).strip() for chain in requested if str(chain).strip()]
    available_list = [str(chain) for chain in available]
    available_set = set(available_list)
    unknown = [chain for chain in cleaned if chain not in available_set]
    if unknown:
        avail = ", ".join(available_list) if available_list else "(none)"
        raise ValueError(
            f"Unknown {side} chain(s): {', '.join(unknown)}. Available: {avail}"
        )
    return cleaned


def require_rna_chains(
    structure_path: Path,
    requested: Sequence[str],
    *,
    side: str,
) -> List[str]:
    """Validate ``requested`` against ``list_rna_chains`` for ``structure_path``."""
    info = list_rna_chains(structure_path)
    return assert_chains_known(requested, info.get("chains") or [], side=side)


def suggest_chain_mapping(
    ref_details: Sequence[Mapping[str, Any]],
    model_details: Sequence[Mapping[str, Any]],
) -> Optional[Dict[str, Any]]:
    """Suggest reference/model chain ids that co-index by sequence (or length).

    Keep in sync with ``suggestChainMapping`` in ``static/app.js``.
    """
    if not ref_details or not model_details:
        return None
    found = _mapping_by_identity(ref_details, model_details)
    if found:
        return found
    return _mapping_by_length(ref_details, model_details)


def assert_concatenated_lengths(
    ref_details: Sequence[Mapping[str, Any]],
    ref_ids: Sequence[str],
    model_details: Sequence[Mapping[str, Any]],
    model_ids: Sequence[str],
) -> None:
    """Raise ValueError when selected chains do not add up to the same length."""
    if not ref_details or not model_details:
        return
    try:
        ref_n = concatenated_length(ref_details, ref_ids)
        model_n = concatenated_length(model_details, model_ids)
    except ValueError:
        return
    if ref_n == model_n:
        return
    ref_label = "+".join(str(item) for item in ref_ids)
    model_label = "+".join(str(item) for item in model_ids)
    message = (
        f"Selected chains differ in length (reference {ref_label} = {ref_n} nt, "
        f"model {model_label} = {model_n} nt). Compare needs the same number of "
        "nucleotides in the same order."
    )
    suggestion = suggest_chain_mapping(ref_details, model_details)
    if suggestion:
        summary = suggestion.get("summary") or ""
        reason = suggestion.get("reason") or ""
        if summary:
            message += f" Suggested mapping: {summary}."
        if reason:
            message += f" {reason}"
    raise ValueError(message)


def _nt_letter(resname: str) -> str:
    """Map a residue name to A/C/G/U/T/N."""
    name = (resname or "").strip()
    if name in {"A", "C", "G", "U"}:
        return name
    if name in {"DA", "DC", "DG", "DT"}:
        return "T" if name == "DT" else name[1]
    parent = fr3d_utils._get_parent_base(name)  # pylint: disable=protected-access
    return parent if parent else "N"


def _is_polymer_nucleotide(resname: str) -> bool:
    """True for standard RNA/DNA letters or a known modified nucleotide."""
    name = (resname or "").strip()
    if name in {"A", "C", "G", "U", "DA", "DC", "DG", "DT"}:
        return True
    # FR3D keeps the modified-nucleotide table private.
    parent = fr3d_utils._get_parent_base(name)  # pylint: disable=protected-access
    return parent is not None


def _auth_from_unit_id(unit_id: str) -> str:
    """Author residue number from an FR3D unit id (field 5)."""
    parts = str(unit_id).split("|")
    if len(parts) >= 5:
        return parts[4]
    return ""


def _pdb_auth(residue: Any) -> str:
    """Author residue number plus insertion code from a BioPython residue."""
    _hetero, seqid, icode = residue.id
    extra = (icode or "").strip()
    return f"{seqid}{extra}" if extra else str(seqid)


def _chain_detail(
    chain_id: str,
    sequence: str,
    auth_start: str,
    auth_end: str,
    model_id: Any,
) -> Dict[str, Any]:
    """One picker row: ids, length, sequence, preview, residue range."""
    seq = sequence or ""
    return {
        "id": str(chain_id),
        "length": len(seq),
        "sequence": seq,
        "preview": sequence_preview(seq),
        "auth_start": auth_start,
        "auth_end": auth_end,
        "model": model_id,
        "same_sequence_as": [],
    }


def _mark_duplicate_sequences(details: List[Dict[str, Any]]) -> None:
    """Annotate chains that share an identical non-empty sequence."""
    for left in details:
        seq = left.get("sequence") or ""
        if not seq:
            left["same_sequence_as"] = []
            continue
        twins = [
            str(right["id"])
            for right in details
            if right["id"] != left["id"] and (right.get("sequence") or "") == seq
        ]
        left["same_sequence_as"] = twins


def _details_from_cif(path: Path) -> List[Dict[str, Any]]:
    """First-model RNA/DNA chains from mmCIF via FR3D."""
    # pylint: disable=import-outside-toplevel
    from fr3d.cif.reader import Cif

    with open_structure_file(path, "r") as handle:
        structure = Cif(handle).structure()
    bases = list(structure.residues(type=["RNA linking"]))
    if not bases:
        bases = list(structure.residues(type=["DNA linking"]))
    if not bases:
        return []
    first_model = bases[0].model
    by_chain: Dict[str, List[Any]] = defaultdict(list)
    for base in bases:
        if base.model != first_model:
            continue
        if fr3d_utils.is_symmetry_mate(base.unit_id()):
            continue
        by_chain[str(base.chain)].append(base)
    details = []
    for chain_id in sorted(by_chain):
        residues = sorted(by_chain[chain_id], key=lambda item: item.index)
        letters = []
        auths = []
        for residue in residues:
            letters.append(_nt_letter(residue.sequence))
            auths.append(_auth_from_unit_id(residue.unit_id()))
        if not letters:
            continue
        details.append(
            _chain_detail(
                chain_id,
                "".join(letters),
                auths[0] if auths else "",
                auths[-1] if auths else "",
                first_model,
            )
        )
    return details


def _details_from_pdb(path: Path) -> List[Dict[str, Any]]:
    """First-model RNA/DNA chains from PDB via BioPython."""
    # pylint: disable=import-outside-toplevel
    from Bio.PDB import PDBParser

    parser = PDBParser(QUIET=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with DecompressedStructureFile(path) as decompressed:
            structure = parser.get_structure("RNA", str(decompressed))
    model = next(iter(structure), None)
    if model is None:
        return []
    by_chain: Dict[str, List[Any]] = {}
    for chain in model:
        residues = []
        for residue in chain:
            if residue.id[0] != " ":
                continue
            resname = residue.resname.strip()
            if _is_polymer_nucleotide(resname):
                residues.append(residue)
        if residues:
            by_chain[str(chain.id)] = residues
    details = []
    for chain_id in sorted(by_chain):
        residues = by_chain[chain_id]
        letters = [_nt_letter(item.resname) for item in residues]
        auths = [_pdb_auth(item) for item in residues]
        details.append(
            _chain_detail(
                chain_id,
                "".join(letters),
                auths[0] if auths else "",
                auths[-1] if auths else "",
                model.id,
            )
        )
    return details


def _ids(details: Sequence[Mapping[str, Any]]) -> List[str]:
    return [str(item["id"]) for item in details]


def _seq_by_id(details: Sequence[Mapping[str, Any]]) -> Dict[str, str]:
    return {str(item["id"]): str(item.get("sequence") or "") for item in details}


def _concat_seq(details: Sequence[Mapping[str, Any]], chain_ids: Sequence[str]) -> str:
    lookup = _seq_by_id(details)
    return "".join(lookup.get(str(chain_id), "") for chain_id in chain_ids)


def _phrase(chain_ids: Sequence[str]) -> str:
    return "+".join(str(item) for item in chain_ids)


def _chain_word(chain_ids: Sequence[str]) -> str:
    if len(chain_ids) == 1:
        return f"chain {chain_ids[0]}"
    return f"chains {_phrase(chain_ids)}"


def _mapping(  # pylint: disable=too-many-arguments
    ref_ids: Sequence[str],
    model_ids: Sequence[str],
    reason: str,
    *,
    identity: Optional[float] = None,
    alternatives: Optional[List[Dict[str, Any]]] = None,
    by: str = "sequence",
) -> Dict[str, Any]:
    return {
        "ref_chains": [str(item) for item in ref_ids],
        "model_chains": [str(item) for item in model_ids],
        "reason": reason,
        "summary": f"{_phrase(ref_ids)} → {_phrase(model_ids)}",
        "identity": identity,
        "alternatives": alternatives or [],
        "by": by,
    }


def _alt_mapping(
    ref_id: str,
    model_ids: Sequence[str],
    *,
    identity: Optional[float],
    by: str,
    primary_id: str,
) -> Dict[str, Any]:
    kind = "sequence" if by == "sequence" else "length"
    return _mapping(
        [ref_id],
        model_ids,
        f"Same {kind} as chain {primary_id}.",
        identity=identity,
        by=by,
    )


def _model_seq_and_ids(
    model_details: Sequence[Mapping[str, Any]],
) -> Tuple[str, List[str], int]:
    model_ids = _ids(model_details)
    model_seq = _concat_seq(model_details, model_ids)
    model_len = len(model_seq) or concatenated_length(model_details, model_ids)
    return model_seq, model_ids, model_len


def _pairwise_identity_mapping(
    ref_details: Sequence[Mapping[str, Any]],
    model_details: Sequence[Mapping[str, Any]],
) -> Optional[Dict[str, Any]]:
    if len(ref_details) != len(model_details) or not ref_details:
        return None
    worst = 1.0
    for ref_item, model_item in zip(ref_details, model_details):
        score = sequence_identity(
            str(ref_item.get("sequence") or ""),
            str(model_item.get("sequence") or ""),
        )
        if score < IDENT_THRESHOLD:
            return None
        worst = min(worst, score)
    ref_ids = _ids(ref_details)
    model_ids = _ids(model_details)
    return _mapping(
        ref_ids,
        model_ids,
        "Each reference chain matches the model chain in the same order.",
        identity=worst,
    )


def _single_ref_matches(
    ref_details: Sequence[Mapping[str, Any]],
    model_seq: str,
) -> List[Tuple[str, float, int]]:
    hits: List[Tuple[str, float, int]] = []
    for item in ref_details:
        seq = str(item.get("sequence") or "")
        score = sequence_identity(seq, model_seq)
        if score >= IDENT_THRESHOLD:
            length = int(item.get("length") or len(seq))
            hits.append((str(item["id"]), score, length))
    if not hits:
        return []
    best = max(item[1] for item in hits)
    return [item for item in hits if item[1] == best]


def _single_ref_mapping(
    ref_details: Sequence[Mapping[str, Any]],
    model_details: Sequence[Mapping[str, Any]],
) -> Optional[Dict[str, Any]]:
    model_seq, model_ids, model_len = _model_seq_and_ids(model_details)
    hits = _single_ref_matches(ref_details, model_seq)
    if not hits:
        return None
    primary_id, identity, length = hits[0]
    alternatives = [
        _alt_mapping(
            item[0],
            model_ids,
            identity=item[1],
            by="sequence",
            primary_id=primary_id,
        )
        for item in hits[1:]
    ]
    reason = (
        f"Model {_chain_word(model_ids)} ({length or model_len} nt) matches "
        f"reference chain {primary_id}."
    )
    if alternatives:
        alt_ids = [item["ref_chains"][0] for item in alternatives]
        if len(alt_ids) == 1:
            reason += (
                f" Chain {alt_ids[0]} has the same sequence — pick that chain "
                "to score the other monomer. Interacting partners can still "
                "appear in the display."
            )
        else:
            listed = ", ".join(alt_ids)
            reason += (
                f" Chains {listed} have the same sequence — pick one of those "
                "to score the other monomer. Interacting partners can still "
                "appear in the display."
            )
    return _mapping(
        [primary_id],
        model_ids,
        reason,
        identity=identity,
        alternatives=alternatives,
    )


def _contiguous_identity_mapping(
    ref_details: Sequence[Mapping[str, Any]],
    model_details: Sequence[Mapping[str, Any]],
) -> Optional[Dict[str, Any]]:
    model_seq, model_ids, model_len = _model_seq_and_ids(model_details)
    ref_ids = _ids(ref_details)
    for width in range(2, len(ref_ids) + 1):
        for start in range(0, len(ref_ids) - width + 1):
            chunk = ref_ids[start : start + width]
            score = sequence_identity(_concat_seq(ref_details, chunk), model_seq)
            if score >= IDENT_THRESHOLD:
                return _mapping(
                    chunk,
                    model_ids,
                    (
                        f"Model {_chain_word(model_ids)} ({model_len} nt) matches "
                        f"reference {_phrase(chunk)} concatenated."
                    ),
                    identity=score,
                )
    return None


def _greedy_identity_mapping(
    ref_details: Sequence[Mapping[str, Any]],
    model_details: Sequence[Mapping[str, Any]],
) -> Optional[Dict[str, Any]]:
    if len(model_details) < 2:
        return None
    chosen: List[str] = []
    used = set()
    worst = 1.0
    for model_item in model_details:
        pick = None
        pick_score = -1.0
        model_seq = str(model_item.get("sequence") or "")
        for ref_item in ref_details:
            ref_id = str(ref_item["id"])
            if ref_id in used:
                continue
            score = sequence_identity(str(ref_item.get("sequence") or ""), model_seq)
            if score >= IDENT_THRESHOLD and score > pick_score:
                pick = ref_id
                pick_score = score
        if pick is None:
            return None
        used.add(pick)
        chosen.append(pick)
        worst = min(worst, pick_score)
    return _mapping(
        chosen,
        _ids(model_details),
        "Each model chain matches a distinct reference chain.",
        identity=worst,
    )


def _mapping_by_identity(
    ref_details: Sequence[Mapping[str, Any]],
    model_details: Sequence[Mapping[str, Any]],
) -> Optional[Dict[str, Any]]:
    return (
        _pairwise_identity_mapping(ref_details, model_details)
        or _single_ref_mapping(ref_details, model_details)
        or _contiguous_identity_mapping(ref_details, model_details)
        or _greedy_identity_mapping(ref_details, model_details)
    )


def _single_length_mapping(
    ref_details: Sequence[Mapping[str, Any]],
    model_details: Sequence[Mapping[str, Any]],
) -> Optional[Dict[str, Any]]:
    _model_seq, model_ids, model_len = _model_seq_and_ids(model_details)
    if model_len <= 0:
        return None
    hits = [
        str(item["id"])
        for item in ref_details
        if int(item.get("length") or 0) == model_len
    ]
    if not hits:
        return None
    primary = hits[0]
    alternatives = [
        _alt_mapping(item, model_ids, identity=None, by="length", primary_id=primary)
        for item in hits[1:]
    ]
    reason = (
        f"Reference chain {primary} and the model are both {model_len} nt "
        "(sequences differ). Confirm this is the intended pair."
    )
    if alternatives:
        alt_ids = [item["ref_chains"][0] for item in alternatives]
        reason += (
            " Chain"
            + ("s " if len(alt_ids) > 1 else " ")
            + ", ".join(alt_ids)
            + " have the same length."
        )
    return _mapping(
        [primary],
        model_ids,
        reason,
        alternatives=alternatives,
        by="length",
    )


def _contiguous_length_mapping(
    ref_details: Sequence[Mapping[str, Any]],
    model_details: Sequence[Mapping[str, Any]],
) -> Optional[Dict[str, Any]]:
    _model_seq, model_ids, model_len = _model_seq_and_ids(model_details)
    if model_len <= 0:
        return None
    ref_ids = _ids(ref_details)
    for width in range(2, len(ref_ids) + 1):
        for start in range(0, len(ref_ids) - width + 1):
            chunk = ref_ids[start : start + width]
            if concatenated_length(ref_details, chunk) == model_len:
                return _mapping(
                    chunk,
                    model_ids,
                    (
                        f"Reference {_phrase(chunk)} and the model are both "
                        f"{model_len} nt concatenated. Sequences were not an "
                        "automatic match — check the order."
                    ),
                    by="length",
                )
    return None


def _mapping_by_length(
    ref_details: Sequence[Mapping[str, Any]],
    model_details: Sequence[Mapping[str, Any]],
) -> Optional[Dict[str, Any]]:
    return _single_length_mapping(
        ref_details, model_details
    ) or _contiguous_length_mapping(ref_details, model_details)
