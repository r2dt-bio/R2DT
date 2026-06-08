"""Read base pairs straight from an mmCIF's own DNATCO/NDB annotation
(``_ndb_base_pair_list`` + ``_ndb_base_pair_annotation``) and write them in the
FR3D ``_basepair.txt`` format the rest of the pipeline already consumes -- so no
FR3D run is needed.

The CIF reading itself is done by the vendored layered-bp-notation scripts
in ``utils/layered_bp_notation/`` (upstream project, not pip-installable yet;
refreshed with ``just sync-lbn`` — see that directory's README). This module
is the thin adapter: it filters to one chain, derives a pseudoknot ``crossing``
level for each cWW pair from the same nesting logic the script's notation uses,
and emits the tab-delimited file keyed by FR3D-style unit ids.
"""

import sys
from pathlib import Path
from typing import Dict, FrozenSet, List, Optional, Tuple

_LBN_DIR = Path(__file__).resolve().parent / "layered_bp_notation"
if str(_LBN_DIR) not in sys.path:
    sys.path.insert(0, str(_LBN_DIR))

import common  # noqa: E402  pylint: disable=wrong-import-position,import-error
import standalone_lbn_script as lbn  # noqa: E402  pylint: disable=wrong-import-position,import-error

Residue = Tuple[str, str, int]  # (chain, symmetry_operator, number)


def has_annotation(cif_path) -> bool:
    """True when the CIF carries the base-pair annotation categories this
    source reads. Standard PDBe mmCIFs do not; DNATCO/NDB-annotated ones do."""
    lines = Path(cif_path).read_text().splitlines()
    return any(l.startswith("_ndb_base_pair_list.") for l in lines) and any(
        l.startswith("_ndb_base_pair_annotation.") for l in lines
    )


def _data_block_id(cif_path) -> Optional[str]:
    """The mmCIF ``data_XXXX`` block name. fr3d-python uses this as the unit-id
    prefix, so synthesised unit ids must match it for the overlay's exact
    lookup (``build_fr3d_data``) to find them."""
    for line in Path(cif_path).read_text().splitlines():
        if line.startswith("data_"):
            return line[len("data_") :].strip()
    return None


def resolve_chain(cif_path, chain: Optional[str] = None) -> Optional[str]:
    """The chain to render: the explicit one, else the first chain that carries
    annotated base pairs. Single source of truth so the basepair file and the
    position map always agree on the chain."""
    if chain:
        return chain
    avail = lbn.list_chains(Path(cif_path))
    return avail[0] if avail else None


def _resolved_residues(cif_path, chain: str) -> List[Tuple[int, str]]:
    """``[(auth_seq_num, base_letter), ...]`` for the chain's resolved residues
    (model 1, asymmetric unit) from ``_atom_site`` via the vendored reader --
    no fr3d-python, and no symmetry-mate expansion."""
    # _read_atom_site is the vendored reader's stable helper.
    # pylint: disable=protected-access
    lines = Path(cif_path).read_text().splitlines()
    return common._read_atom_site(lines, [chain]).get(chain, [])


def _unit_id(prefix: str, chain: str, letter: str, num: int) -> str:
    """FR3D-style unit id. Built identically for the basepair file and the
    position map so the overlay's exact lookup (``build_fr3d_data``) matches."""
    return f"{prefix}|1|{chain}|{letter.upper()}|{num}"


def _normalize_lw(lw: Optional[str]) -> Optional[str]:
    """Canonicalise the Leontis-Westhof edge letters (cWw -> cWW) while keeping
    the cis/trans prefix, so families match the viewer's checkboxes and the 2D
    plugin's path classes."""
    if lw and len(lw) == 3 and lw[0] in "ct":
        return lw[0] + lw[1].upper() + lw[2].upper()
    return lw


def _crossing_levels(
    cww_pairs: List[Tuple[Residue, Residue, str]],
    residues: List[Tuple[int, str]],
    chain: str,
) -> Dict[FrozenSet, int]:
    """Pseudoknot level per cWW pair: 0 for a nested pair, 1+ for each crossing
    layer. This mirrors the greedy class assignment in the vendored script's
    ``render()`` (``common._cross``), so the levels match the notation the
    script prints -- a nested pair is ``()``, a crossing pair ``[]``/``{}``/...
    """
    # common._cross is the vendored crossing predicate.
    # pylint: disable=protected-access
    col_of = {(chain, common.IDENTITY, num): k for k, (num, _) in enumerate(residues)}
    colpairs: List[Tuple[Tuple[int, int], FrozenSet]] = []
    for res_a, res_b, _ in cww_pairs:
        # A paired residue may be absent from the resolved set (e.g. modified
        # residue dropped by the atom_site reader); skip -- it just won't render.
        if res_a in col_of and res_b in col_of:
            i, j = sorted((col_of[res_a], col_of[res_b]))
            colpairs.append(((i, j), frozenset((res_a, res_b))))

    classed: List[Tuple[Tuple[int, int], int]] = []
    levels: Dict[FrozenSet, int] = {}
    for (i, j), key in sorted(colpairs):
        cls = 0
        while any(c == cls and common._cross((i, j), pq) for pq, c in classed):
            cls += 1
        classed.append(((i, j), cls))
        levels[key] = cls
    return levels


def write_basepair_txt(
    cif_path,
    output_txt_path,
    chain: Optional[str] = None,
    quiet: bool = False,
) -> Optional[str]:
    """Write an FR3D-format ``{id}_basepair.txt`` from the CIF annotation.

    Each line is ``unit1<TAB>lw<TAB>unit2<TAB>crossing`` where a unit id is
    ``{prefix}|1|{chain}|{nt}|{auth_seq_id}`` -- the same shape FR3D emits,
    matched downstream by chain + residue number. Only intra-chain,
    identity-symmetry pairs of the selected chain are written (the 2D pipeline
    is single-chain). Returns the chain used, or None if none could be resolved.
    """
    cif = Path(cif_path)
    chain = resolve_chain(cif, chain)
    if chain is None:
        return None

    residues = _resolved_residues(cif, chain)
    letter_of = dict(residues)
    prefix = _data_block_id(cif) or cif.stem

    # Keep only intra-chain pairs in the asymmetric unit (identity operator).
    pairs = [
        (res_a, res_b, _normalize_lw(lw))
        for res_a, res_b, lw in lbn.read_pairs(cif, {chain})
        if res_a[0] == chain
        and res_b[0] == chain
        and res_a[1] == common.IDENTITY
        and res_b[1] == common.IDENTITY
    ]

    cww_pairs = [(a, b, lw) for a, b, lw in pairs if lw == "cWW"]
    levels = _crossing_levels(cww_pairs, residues, chain)

    lines = []
    for res_a, res_b, lw in pairs:
        unit1 = _unit_id(prefix, chain, letter_of.get(res_a[2], "N"), res_a[2])
        unit2 = _unit_id(prefix, chain, letter_of.get(res_b[2], "N"), res_b[2])
        crossing = levels.get(frozenset((res_a, res_b)), 0)
        lines.append(f"{unit1}\t{lw}\t{unit2}\t{crossing}")

    out = Path(output_txt_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + ("\n" if lines else ""))

    if not quiet:
        print(
            f"CIF annotation: {len(pairs)} pairs on chain {chain} "
            f"({len(cww_pairs)} cWW) -> {out}"
        )
    return chain


def read_sequence_and_positions(
    cif_path,
    chain: Optional[str] = None,
    quiet: bool = False,
) -> Tuple[str, Dict[str, int], Optional[str]]:
    """Resolved RNA sequence + ``{unit_id: 0-based position}`` for one chain,
    read from ``_atom_site`` via the vendored script -- the fr3d-python-free
    replacement for ``fr3d.extract_sequence_from_cif``. Unit ids use the same
    convention as :func:`write_basepair_txt`, so the basepair file and this map
    line up exactly. Returns ``(sequence, positions, chain)``.
    """
    cif = Path(cif_path)
    chain = resolve_chain(cif, chain)
    if chain is None:
        return "", {}, None
    residues = _resolved_residues(cif, chain)
    prefix = _data_block_id(cif) or cif.stem
    sequence = "".join(letter.upper() for _, letter in residues)
    positions = {
        _unit_id(prefix, chain, letter, num): i
        for i, (num, letter) in enumerate(residues)
    }
    if not quiet:
        print(f"CIF sequence: {len(sequence)} nt on chain {chain}")
    return sequence, positions, chain


# pylint: disable=too-many-arguments,too-many-positional-arguments
def get_secondary_structure_cif(
    structure_file: str,
    output_dir: str,
    structure_id: str,
    chain_id: Optional[str] = None,
    include_pseudoknots: bool = False,
    quiet: bool = False,
) -> Tuple[Optional[str], Optional[str]]:
    """CIF-annotation counterpart to ``fr3d.get_secondary_structure_fr3d``.

    Reads the sequence, residue positions and base pairs entirely from the
    vendored script (no fr3d-python), writes the ``_basepair.txt``, and folds it
    to a dot-bracket with R2DT's pure-Python parser. Returns
    ``(sequence, dot_bracket)`` or ``(None, None)``.
    """
    # parse_fr3d_basepairs is pure Python; importing utils.fr3d does not pull in
    # the fr3d-python package (that import is lazy, inside the reader functions).
    from utils import fr3d as fr3d_utils  # pylint: disable=import-outside-toplevel

    sequence, unit_id_to_position, used_chain = read_sequence_and_positions(
        structure_file, chain_id, quiet=quiet
    )
    if not sequence or used_chain is None:
        return None, None

    txt = Path(output_dir) / f"{structure_id}_basepair.txt"
    write_basepair_txt(structure_file, txt, chain=used_chain, quiet=quiet)

    dot_bracket = fr3d_utils.parse_fr3d_basepairs(
        str(txt), unit_id_to_position, len(sequence), include_pseudoknots
    )
    return sequence, dot_bracket
