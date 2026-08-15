"""
Build layered dot-bracket notation (LBN) from fr3d_data + api_data for the
interactive viewer panel.

Works entirely in "label space" (1-based residue positions matching
api_data.label_seq_ids), so no CIF re-reading or coordinate remapping is
needed.  Pairs come from fr3d_data.annotations with seq_id1/seq_id2 already
in that space.

Output format::

    {
        "sequence": "AAGAAGAGU...",   # length N
        "rows": [
            {
                "label": "WC",        # LW family label
                "chars": ".((...",    # length N, dot-bracket string
                "partners": {         # 1-based pos (str) -> partner pos (int)
                    "2": 59,
                    "59": 2,
                    ...
                }
            },
            ...                       # one entry per non-empty LW layer
        ]
    }
"""

from typing import Dict, List, Tuple

# Bracket levels for pseudoknots within one layer.
_BRACKETS = ["()", "[]", "{}", "<>", "Aa", "Bb", "Cc", "Dd"]

# Canonical Watson-Crick base-pair combinations.
_WC_PAIRS = (frozenset({"A", "U"}), frozenset({"G", "C"}), frozenset({"A", "T"}))

# Canonical LW family order (excluding cWW, which is split into WC / cWW rows).
_DIRECTED_ORDER = [
    "cWH", "cWS", "cHW", "cHH", "cHS", "cSW", "cSH", "cSS",
    "tWW", "tWH", "tWS", "tHW", "tHH", "tHS", "tSW", "tSH", "tSS",
]


def _is_wc(bp: str, nt1: str, nt2: str) -> bool:
    return bp == "cWW" and frozenset({nt1.upper(), nt2.upper()}) in _WC_PAIRS


def _cross(i: int, j: int, k: int, l: int) -> bool:
    """True iff pairs (i,j) and (k,l) form a pseudoknot."""
    return i < k < j < l or k < i < l < j


def _pack(
    pairs: List[Tuple[int, int]],
) -> List[List[Tuple[int, int]]]:
    """Group pairs into the minimum number of buckets so each residue position
    appears at most once per bucket (handles multi-pairing residues)."""
    buckets: List[List[Tuple[int, int]]] = []
    for p in sorted(pairs):
        a, b = p
        for bucket in buckets:
            used = {x for x, y in bucket} | {y for x, y in bucket}
            if a not in used and b not in used:
                bucket.append(p)
                break
        else:
            buckets.append([p])
    return buckets


def _render_layer(
    pairs: List[Tuple[int, int]],
    n: int,
) -> Tuple[str, Dict[int, int]]:
    """Assign bracket characters to a set of pairs.

    Crossing pairs (pseudoknots) are placed on higher bracket levels
    (``[]``, ``{}``, …).

    Args:
        pairs: (pos1, pos2) with pos1 < pos2, 1-based.
        n: total sequence length.

    Returns:
        chars: dot-bracket string of length *n*.
        partners: mapping 1-based pos -> 1-based paired partner pos.
    """
    chars = ["."] * n
    partners: Dict[int, int] = {}
    placed: List[Tuple[int, int, int]] = []  # (i, j, bracket_level)

    for i, j in sorted(pairs):
        if not (1 <= i < j <= n):
            continue
        level = 0
        while level < len(_BRACKETS):
            if not any(
                lv == level and _cross(i, j, pi, pj)
                for pi, pj, lv in placed
            ):
                break
            level += 1
        level = min(level, len(_BRACKETS) - 1)

        placed.append((i, j, level))
        chars[i - 1] = _BRACKETS[level][0]
        chars[j - 1] = _BRACKETS[level][1]
        partners[i] = j
        partners[j] = i

    return "".join(chars), partners


def build_lbn_data(api_data: dict, fr3d_data: dict) -> dict:
    """Build the LBN dataset for the viewer's notation panel.

    Parameters
    ----------
    api_data:
        The dict written to ``api.json`` by :func:`viewer_export.build_api_data`.
    fr3d_data:
        The dict written to ``fr3d.json`` by :func:`viewer_export.build_fr3d_data`.

    Returns
    -------
    dict
        ``{sequence, rows}`` ready to be serialised as ``lbn.json``.
    """
    sequence: str = api_data.get("sequence", "")
    n = len(sequence)
    if n == 0:
        return {"sequence": "", "rows": []}

    wc_pairs: List[Tuple[int, int]] = []
    cww_pairs: List[Tuple[int, int]] = []
    by_type: Dict[str, List[Tuple[int, int]]] = {}

    for ann in fr3d_data.get("annotations", []):
        try:
            p1 = int(ann["seq_id1"])
            p2 = int(ann["seq_id2"])
        except (KeyError, ValueError, TypeError):
            continue
        if p1 == p2 or not (1 <= p1 <= n) or not (1 <= p2 <= n):
            continue
        bp = ann.get("bp", "")
        if not bp:
            continue
        if p1 > p2:
            p1, p2 = p2, p1
            nt1 = ann.get("nt2", "")
            nt2 = ann.get("nt1", "")
        else:
            nt1 = ann.get("nt1", "")
            nt2 = ann.get("nt2", "")

        if _is_wc(bp, nt1, nt2):
            wc_pairs.append((p1, p2))
        elif bp == "cWW":
            cww_pairs.append((p1, p2))
        else:
            by_type.setdefault(bp, []).append((p1, p2))

    # Deduplicate (FR3D occasionally reports a pair twice).
    wc_pairs = list(set(wc_pairs))
    cww_pairs = list(set(cww_pairs))
    for k in by_type:
        by_type[k] = list(set(by_type[k]))

    rows: List[dict] = []

    def _emit(label: str, pairs: List[Tuple[int, int]]) -> None:
        for bucket_idx, bucket in enumerate(_pack(pairs)):
            lbl = label if bucket_idx == 0 else f"{label}({bucket_idx + 1})"
            chars, partners = _render_layer(bucket, n)
            if partners:
                rows.append({
                    "label": lbl,
                    "chars": chars,
                    "partners": {str(k): v for k, v in partners.items()},
                })

    if wc_pairs:
        _emit("WC", wc_pairs)
    if cww_pairs:
        _emit("cWW", cww_pairs)
    for bp_type in _DIRECTED_ORDER:
        if bp_type in by_type:
            _emit(bp_type, by_type[bp_type])
    # Any unrecognised family names that weren't in DIRECTED_ORDER
    for bp_type, pairs in sorted(by_type.items()):
        if bp_type not in _DIRECTED_ORDER:
            _emit(bp_type, pairs)

    return {"sequence": sequence, "rows": rows}
