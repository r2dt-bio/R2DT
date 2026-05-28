"""
Convert R2DT outputs into the JSON payload consumed by ``pdb-rna-viewer``.

The viewer expects two blobs (see ``pdb-rna-viewer/src/app/data.ts``)::

    apiData = {
        svg_paths: string[],          # length N + 2
        dimensions: {width, height},
        sequence: string,             # length N
        label_seq_ids: number[],      # [None, 1, 2, ..., N, None]
        auth_seq_ids: number[],       # same shape as label_seq_ids
        unobserved_label_seq_ids: number[],
        pdb_ins_codes: string[],      # length N
    }

    FR3DData = {
        pdb_id, chain_id, modified: [], annotations: [
            {seq_id1, nt1, bp, seq_id2, nt2, crossing, ...}
        ]
    }

Both blobs use 1-based label indexing aligned to R2DT's full (mask-expanded)
sequence.  FR3D ``seq_id1``/``seq_id2`` are remapped from the structure's
author numbering into the same 1-based label space so the viewer can
position base-pair glyphs via ``locations.get(seq_id)``.
"""

import json
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def _real_nucleotides(colored_json: dict) -> List[dict]:
    """Return the nucleotide records from a Traveler-style JSON, dropping
    the synthetic 5'/3' label entries that R2DT emits at the ends."""
    mol = colored_json["rnaComplexes"][0]["rnaMolecules"][0]
    nucs = []
    for n in mol["sequence"]:
        name = n.get("residueName", "")
        if name in ("5'", "3'"):
            continue
        if len(name) != 1:
            continue
        nucs.append(n)
    nucs.sort(key=lambda n: n["residueIndex"])
    return nucs


def _svg_dimensions(svg_path: Path) -> Optional[Tuple[float, float]]:
    """Pull width/height from the colored SVG's root attributes."""
    if not svg_path.exists():
        return None
    head = svg_path.read_text(errors="ignore")[:2048]
    w = re.search(r'width="([\d.]+)"', head)
    h = re.search(r'height="([\d.]+)"', head)
    if w and h:
        return float(w.group(1)), float(h.group(1))
    return None


def _build_svg_paths(positions: List[Tuple[float, float]]) -> List[str]:
    """
    Encode nucleotide positions in the form the viewer parses.

    The viewer splits each ``svg_paths[i]`` on ``M`` and ``,`` then reads
    indices ``[1..4]`` as ``x1, y1, x2, y2``.  The midpoint becomes the
    nucleotide centre; consecutive midpoints are joined to draw the
    backbone.  Index ``0`` and ``1`` are dummy entries.
    """
    n = len(positions)
    paths: List[str] = []
    # Two dummy entries placing the start 5 px before nucleotide 1.
    x1, y1 = positions[0]
    dummy = f"M{x1 - 5},{y1 - 5},{x1},{y1}"
    paths.append(dummy)
    paths.append(dummy)
    # One path per nucleotide: previous → current.  For the very first
    # nucleotide use the dummy start as "previous".
    prev_x, prev_y = x1 - 5, y1 - 5
    for i, (x, y) in enumerate(positions):
        paths.append(f"M{prev_x},{prev_y},{x},{y}")
        prev_x, prev_y = x, y
        if i == n - 1:
            # Extra trailing segment mirroring the last step (the viewer
            # extrapolates here on its own, but the length must be N+2).
            pass
    return paths


def _parse_unit_id(unit_id: str) -> Tuple[str, str, str, int, str]:
    """Split ``pdb|model|chain|nt|seq_id[|ins_code]`` into typed parts.

    Returns ``(pdb_id, chain, nt, auth_seq_id, ins_code)``.  Raises
    ``ValueError`` if the field cannot be parsed.
    """
    parts = unit_id.split("|")
    if len(parts) < 5:
        raise ValueError(f"unexpected unit_id: {unit_id!r}")
    pdb_id = parts[0]
    chain = parts[2]
    nt = parts[3]
    auth_seq_id = int(parts[4])
    ins_code = parts[5] if len(parts) > 5 else ""
    return pdb_id, chain, nt, auth_seq_id, ins_code


def _resolved_to_full(mask: Optional[List[bool]], n_full: int) -> List[int]:
    """Return a list where ``out[j] = full_position`` (0-based) for the
    ``j``-th resolved nucleotide.  Without a mask every position is
    considered resolved."""
    if mask is None:
        return list(range(n_full))
    out = [i for i, ok in enumerate(mask) if ok]
    return out


# pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-branches
def build_api_data(
    colored_json_path: Path,
    structure_id: str,
    chain_id: Optional[str],
    resolved_mask: Optional[List[bool]],
    unit_id_to_position: Optional[Dict[str, int]],
    colored_svg_path: Optional[Path] = None,
) -> dict:
    """Build the ``apiData`` payload for ``pdb-rna-viewer``.

    Args:
        colored_json_path: Path to the Traveler-style ``*.colored.json``
            produced by R2DT's ``templatefree`` pipeline.  Nucleotide
            positions and the rendered sequence come from here.
        structure_id: Structure identifier (informational only).
        chain_id: Chain used during extraction (informational only).
        resolved_mask: Boolean mask aligned to the full deposited
            sequence (``True`` = resolved).  ``None`` if no expansion was
            performed.
        unit_id_to_position: Mapping returned by
            :func:`utils.fr3d.extract_sequence_from_pdb`/``..._cif``.
            Used to recover author residue numbers.  ``None`` falls back
            to label numbering for ``auth_seq_ids``.
        colored_svg_path: Optional path to the colored SVG, used to
            recover canvas dimensions when present.
    """
    colored = json.loads(colored_json_path.read_text())
    nucs = _real_nucleotides(colored)
    n = len(nucs)

    sequence = "".join(n_["residueName"] for n_ in nucs)
    positions = [(n_["x"], n_["y"]) for n_ in nucs]

    if colored_svg_path is not None:
        dim = _svg_dimensions(colored_svg_path)
    else:
        dim = None
    if dim is None:
        xs = [p[0] for p in positions]
        ys = [p[1] for p in positions]
        dim = (max(xs) + 20, max(ys) + 20)

    label_seq_ids: List[Optional[int]] = [None] + list(range(1, n + 1)) + [None]

    # auth_seq_ids: map each full (1-based) position back to the author
    # residue number via the resolved→full index.
    auth_per_full: List[Optional[int]] = [None] * n
    if unit_id_to_position:
        resolved_to_full = _resolved_to_full(resolved_mask, n)
        pos_to_unit = {v: k for k, v in unit_id_to_position.items()}
        for resolved_idx, full_idx in enumerate(resolved_to_full):
            unit = pos_to_unit.get(resolved_idx)
            if unit is None:
                continue
            try:
                _, _, _, auth, _ = _parse_unit_id(unit)
            except ValueError:
                continue
            if 0 <= full_idx < n:
                auth_per_full[full_idx] = auth
    # Fill unresolved positions by interpolating between neighbours so
    # molstar still gets a number it can click through.
    for i, v in enumerate(auth_per_full):
        if v is None:
            auth_per_full[i] = i + 1  # fallback: label numbering
    auth_seq_ids = [None] + auth_per_full + [None]

    unobserved: List[int] = []
    if resolved_mask is not None:
        for i, ok in enumerate(resolved_mask):
            if not ok:
                unobserved.append(i + 1)

    pdb_ins_codes = [""] * n

    return {
        "svg_paths": _build_svg_paths(positions),
        "dimensions": {"width": dim[0], "height": dim[1]},
        "sequence": sequence,
        "label_seq_ids": label_seq_ids,
        "auth_seq_ids": auth_seq_ids,
        "unobserved_label_seq_ids": unobserved,
        "pdb_ins_codes": pdb_ins_codes,
        "_meta": {
            "structure_id": structure_id,
            "chain_id": chain_id,
        },
    }


# pylint: disable=too-many-arguments,too-many-positional-arguments
def build_fr3d_data(
    basepair_txt_path: Path,
    structure_id: str,
    chain_id: Optional[str],
    unit_id_to_position: Dict[str, int],
    resolved_mask: Optional[List[bool]],
    n_full: int,
) -> dict:
    """Convert the FR3D tab-delimited basepair file into the JSON shape
    the viewer expects.

    Each line in the file looks like::

        1Y26|1|X|C|13<TAB>cWW<TAB>1Y26|1|X|G|83<TAB>0

    ``seq_id1``/``seq_id2`` in the output are 1-based positions in the
    full (mask-expanded) sequence -- the same space as ``label_seq_ids``
    in ``apiData`` -- because the viewer indexes nucleotide locations by
    that value.
    """
    resolved_to_full = _resolved_to_full(resolved_mask, n_full)

    annotations = []
    if not basepair_txt_path.exists():
        return {
            "pdb_id": structure_id,
            "chain_id": chain_id or "",
            "modified": [],
            "annotations": annotations,
        }

    # FR3D reports every pair twice -- once in each direction (e.g. both
    # "19 cSS 22" and "22 cSS 19"). The viewer draws one glyph per
    # annotation, so a second, direction-reversed copy renders on top of
    # the first. For symmetric symbols (cWW circle, cHH square) the copies
    # coincide and it's invisible, but the Sugar-edge triangle is oriented
    # along the pair axis, so the reversed copy points the opposite way and
    # the two triangles visibly overlap. Keep only the first occurrence of
    # each unordered pair.
    seen_pairs = set()

    for raw in basepair_txt_path.read_text().splitlines():
        row = raw.strip()
        if not row:
            continue
        parts = row.split("\t")
        if len(parts) < 4:
            continue
        unit1, bp, unit2, crossing = parts[0], parts[1], parts[2], parts[3]
        try:
            _, _, nt1, auth1, _ = _parse_unit_id(unit1)
            _, _, nt2, auth2, _ = _parse_unit_id(unit2)
        except ValueError:
            continue
        pos1 = unit_id_to_position.get(unit1)
        pos2 = unit_id_to_position.get(unit2)
        if pos1 is None or pos2 is None:
            continue
        if pos1 >= len(resolved_to_full) or pos2 >= len(resolved_to_full):
            continue
        pair_key = frozenset((pos1, pos2))
        if pair_key in seen_pairs:
            continue
        seen_pairs.add(pair_key)
        seq_id1 = resolved_to_full[pos1] + 1
        seq_id2 = resolved_to_full[pos2] + 1
        annotations.append(
            {
                "seq_id1": str(seq_id1),
                "3d_id1": str(auth1),
                "nt1": nt1,
                "unit1": nt1,
                "bp": bp,
                "seq_id2": str(seq_id2),
                "nt2": nt2,
                "unit2": nt2,
                "3d_id2": str(auth2),
                "crossing": crossing,
            }
        )

    return {
        "pdb_id": structure_id,
        "chain_id": chain_id or "",
        "modified": [],
        "annotations": annotations,
    }
