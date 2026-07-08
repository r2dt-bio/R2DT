"""
Multi-chain RNA extraction and concatenation for combined 2D diagrams.

Given a structure with two or more base-pairing RNA chains, this module
extracts each selected chain, concatenates them into a single label space,
and maps the FR3D base pairs — including *inter-chain* pairs, which the
single-chain path drops — onto that concatenation.

Base pairs are then partitioned into a maximum nested (pseudoknot-free) set,
which the layout engine can draw natively, and a remaining crossing set, which
is destined for an overlay.  The nested/crossing split is recomputed on the
concatenation itself rather than read from FR3D's own crossing column: FR3D
computes crossing on the 3D graph, but a loop-loop / kissing helix that is
"crossing=0" there becomes a pseudoknot once the two chains are laid out on one
line (see docs/multichain-2d-plan.md, 2D19).

Chain order is chosen automatically to maximize the nested set (minimize the
overlay); a caller may pin an explicit order instead.

Copyright [2009-present] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
     http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import itertools
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from utils import fr3d as fr3d_utils

# Colour-blind-safe palette used to tag each chain's 5'/3' termini.
_CHAIN_PALETTE = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02"]
# Overlay style for crossing (pseudoknot / inter-chain kissing) pairs.
_OVERLAY_STROKE = "#e6194b"

# Above this many chains the N! order search is skipped in favour of a greedy
# adjacency ordering (see _greedy_order).  Real RNA complexes sit well below
# this; the cap only guards against pathological inputs.
_MAX_BRUTEFORCE_CHAINS = 6

_CWW_INTERACTIONS = ("cWW", "cWw", "cwW")


@dataclass
class ChainSeq:
    """One extracted chain: its sequence and unit-id → local-position map."""

    chain_id: str
    sequence: str
    unit_id_to_local: Dict[str, int]  # 0-based position within this chain


@dataclass
class MultiChainStructure:  # pylint: disable=too-many-instance-attributes
    """Result of concatenating selected chains into one label space."""

    order: List[str]  # chain ids in concatenation order
    sequence: str  # concatenated sequence
    dot_bracket: str  # nested pairs only (crossing pairs excluded)
    nested_pairs: List[Tuple[int, int]]  # 0-based, i < j
    crossing_pairs: List[Tuple[int, int]]  # 0-based, i < j — for the overlay
    chain_of: List[str]  # per-position chain id, len == len(sequence)
    boundaries: List[Tuple[str, int, int]]  # (chain_id, start, end_exclusive)
    inter_pairs: List[Tuple[int, int]] = field(default_factory=list)

    def counts(self) -> Dict[str, int]:
        """Small summary useful for logging and the order scorer."""
        inter = set(self.inter_pairs)
        nested = set(self.nested_pairs)
        crossing = set(self.crossing_pairs)
        return {
            "length": len(self.sequence),
            "chains": len(self.order),
            "pairs": len(nested) + len(crossing),
            "nested": len(nested),
            "crossing": len(crossing),
            "inter": len(inter),
            "inter_nested": len(inter & nested),
            "inter_crossing": len(inter & crossing),
        }


# ---------------------------------------------------------------------------
# Chain discovery and extraction
# ---------------------------------------------------------------------------


def list_rna_chains(cif_file: str, model_id: Optional[int] = None) -> List[str]:
    """Return the RNA chain ids present in the first (or given) model, in
    structure order."""
    from fr3d.cif.reader import Cif  # pylint: disable=import-outside-toplevel

    with fr3d_utils.open_structure_file(Path(cif_file), "r") as handle:
        structure = Cif(handle).structure()

    bases = list(structure.residues(type=["RNA linking"]))
    if not bases:
        return []

    target_model = model_id if model_id is not None else bases[0].model
    seen: List[str] = []
    for base in bases:
        if base.model != target_model:
            continue
        if base.chain not in seen:
            seen.append(base.chain)
    return seen


def extract_chains(
    cif_file: str,
    chains: List[str],
    model_id: Optional[int] = None,
    quiet: bool = False,
) -> List[ChainSeq]:
    """Extract each named chain's sequence and unit-id map (one model)."""
    out: List[ChainSeq] = []
    for chain_id in chains:
        seq, unit_map = fr3d_utils.extract_sequence_from_cif(
            cif_file, chain_id=chain_id, model_id=model_id, quiet=quiet
        )
        if not seq:
            if not quiet:
                print(f"[multichain] chain {chain_id}: no RNA residues, skipping")
            continue
        out.append(ChainSeq(chain_id, seq, unit_map))
    return out


# ---------------------------------------------------------------------------
# Base-pair reading and concatenation
# ---------------------------------------------------------------------------


def read_cww_pairs(
    basepair_file: str, target_model: Optional[str] = None
) -> List[Tuple[str, str]]:
    """Read (unit_id1, unit_id2) cWW-type pairs from an FR3D basepair file.

    Optionally restrict to a single model (needed for NMR ensembles, where the
    file carries every model).  Pairs are deduplicated (unordered).
    """
    seen = set()
    pairs: List[Tuple[str, str]] = []
    with open(basepair_file, "r") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            unit1, interaction, unit2 = parts[0], parts[1], parts[2]
            if interaction not in _CWW_INTERACTIONS:
                continue
            if target_model is not None:
                m1 = unit1.split("|")[1] if len(unit1.split("|")) > 1 else None
                m2 = unit2.split("|")[1] if len(unit2.split("|")) > 1 else None
                if m1 != target_model or m2 != target_model:
                    continue
            key = frozenset((unit1, unit2))
            if key in seen:
                continue
            seen.add(key)
            pairs.append((unit1, unit2))
    return pairs


def _combined_map(
    chain_data: List[ChainSeq], order: List[str]
) -> Tuple[Dict[str, int], List[str], List[Tuple[str, int, int]], str]:
    """Build the concatenated sequence, unit-id → global-position map,
    per-position chain list, and chain boundaries for a given order."""
    by_id = {c.chain_id: c for c in chain_data}
    combined: Dict[str, int] = {}
    chain_of: List[str] = []
    boundaries: List[Tuple[str, int, int]] = []
    sequence_parts: List[str] = []
    offset = 0
    for chain_id in order:
        chain = by_id[chain_id]
        for unit_id, local in chain.unit_id_to_local.items():
            combined[unit_id] = offset + local
        chain_of.extend([chain_id] * len(chain.sequence))
        boundaries.append((chain_id, offset, offset + len(chain.sequence)))
        sequence_parts.append(chain.sequence)
        offset += len(chain.sequence)
    return combined, chain_of, boundaries, "".join(sequence_parts)


def _map_pairs(
    cww_pairs: List[Tuple[str, str]], combined: Dict[str, int]
) -> List[Tuple[int, int]]:
    """Map unit-id pairs to 0-based (i, j) positions in the concatenation.

    Keeps a pair only when *both* endpoints resolve — which now includes
    inter-chain pairs, since the combined map spans every selected chain.
    """
    out: List[Tuple[int, int]] = []
    for unit1, unit2 in cww_pairs:
        pos1 = (
            fr3d_utils._get_position_from_unit_id(  # pylint: disable=protected-access
                unit1, combined
            )
        )
        pos2 = (
            fr3d_utils._get_position_from_unit_id(  # pylint: disable=protected-access
                unit2, combined
            )
        )
        if pos1 is None or pos2 is None or pos1 == pos2:
            continue
        out.append((min(pos1, pos2), max(pos1, pos2)))
    # Enforce a matching (each position in one pair); drop later duplicates.
    used = set()
    matched: List[Tuple[int, int]] = []
    for i, j in sorted(set(out)):
        if i in used or j in used:
            continue
        used.add(i)
        used.add(j)
        matched.append((i, j))
    return matched


# ---------------------------------------------------------------------------
# Nested / crossing partition
# ---------------------------------------------------------------------------


def max_noncrossing(
    pairs: List[Tuple[int, int]], n: int
) -> Tuple[List[Tuple[int, int]], List[Tuple[int, int]]]:
    """Split a matching into a maximum nested subset and the crossing remainder.

    Interval DP over positions: ``dp[i][j]`` = max arcs fully inside ``[i, j]``.
    ``O(n^2)`` states, ``O(1)`` transition (each position opens at most one arc).
    Returns ``(nested, crossing)`` with the nested set as large as possible.
    """
    if not pairs:
        return [], []
    opener: Dict[int, int] = dict(pairs)

    # dp[i][j] via a flat dict keyed on (i, j); build by increasing length.
    dp: Dict[Tuple[int, int], int] = {}

    def get(i: int, j: int) -> int:
        if i > j:
            return 0
        return dp[(i, j)]

    for length in range(1, n + 1):
        for i in range(0, n - length + 1):
            j = i + length - 1
            best = get(i + 1, j)  # skip position i
            k = opener.get(i)
            if k is not None and k <= j:
                use = 1 + get(i + 1, k - 1) + get(k + 1, j)
                best = max(best, use)
            dp[(i, j)] = best

    # Reconstruct the chosen arcs with an explicit stack (avoid deep recursion).
    nested: List[Tuple[int, int]] = []
    stack: List[Tuple[int, int]] = [(0, n - 1)]
    while stack:
        i, j = stack.pop()
        if i > j:
            continue
        k = opener.get(i)
        if (
            k is not None
            and k <= j
            and dp[(i, j)] == 1 + get(i + 1, k - 1) + get(k + 1, j)
        ):
            nested.append((i, k))
            stack.append((i + 1, k - 1))
            stack.append((k + 1, j))
        else:
            stack.append((i + 1, j))

    nested_set = set(nested)
    crossing = [p for p in pairs if p not in nested_set]
    return sorted(nested_set), sorted(crossing)


def dot_bracket_from_pairs(pairs: List[Tuple[int, int]], n: int) -> str:
    """Render a nested (pseudoknot-free) pair set as a dot-bracket string."""
    chars = ["."] * n
    for i, j in pairs:
        chars[i] = "("
        chars[j] = ")"
    return "".join(chars)


# ---------------------------------------------------------------------------
# Order selection and assembly
# ---------------------------------------------------------------------------


def _is_inter(pair: Tuple[int, int], boundaries: List[Tuple[str, int, int]]) -> bool:
    """True if the two positions fall in different chain blocks."""
    i, j = pair

    def chain_of(pos: int) -> str:
        for cid, start, end in boundaries:
            if start <= pos < end:
                return cid
        return ""

    return chain_of(i) != chain_of(j)


def _build_for_order(
    chain_data: List[ChainSeq],
    order: List[str],
    cww_pairs: List[Tuple[str, str]],
) -> MultiChainStructure:
    combined, chain_of, boundaries, sequence = _combined_map(chain_data, order)
    mapped = _map_pairs(cww_pairs, combined)
    nested, crossing = max_noncrossing(mapped, len(sequence))
    inter = [p for p in mapped if _is_inter(p, boundaries)]
    return MultiChainStructure(
        order=list(order),
        sequence=sequence,
        dot_bracket=dot_bracket_from_pairs(nested, len(sequence)),
        nested_pairs=nested,
        crossing_pairs=crossing,
        chain_of=chain_of,
        boundaries=boundaries,
        inter_pairs=inter,
    )


def _greedy_order(
    chain_data: List[ChainSeq], cww_pairs: List[Tuple[str, str]]
) -> List[str]:
    """Fallback ordering for many chains: place the most strongly base-paired
    chain pairs adjacent.  Only used above _MAX_BRUTEFORCE_CHAINS."""
    ids = [c.chain_id for c in chain_data]
    # Count inter-chain pairs per unordered chain pair.
    unit_chain: Dict[str, str] = {}
    for chain in chain_data:
        for unit_id in chain.unit_id_to_local:
            unit_chain[unit_id] = chain.chain_id
    weight: Dict[frozenset, int] = {}
    for unit1, unit2 in cww_pairs:
        c1 = unit_chain.get(unit1)
        c2 = unit_chain.get(unit2)
        if not c1 or not c2 or c1 == c2:
            continue
        key = frozenset((c1, c2))
        weight[key] = weight.get(key, 0) + 1
    order = [sorted(ids)[0]]
    remaining = set(ids) - set(order)
    while remaining:
        last = order[-1]
        nxt = max(
            sorted(remaining),
            key=lambda c: weight.get(frozenset((last, c)), 0),
        )
        order.append(nxt)
        remaining.discard(nxt)
    return order


def choose_order(
    chain_data: List[ChainSeq], cww_pairs: List[Tuple[str, str]]
) -> MultiChainStructure:
    """Pick the concatenation order maximizing the nested set (min overlay).

    Brute-forces permutations for small chain counts; ties break
    deterministically (lexicographic by order) for reproducible diagrams.
    """
    ids = [c.chain_id for c in chain_data]
    if len(ids) > _MAX_BRUTEFORCE_CHAINS:
        candidates = [_greedy_order(chain_data, cww_pairs)]
    else:
        candidates = [list(p) for p in itertools.permutations(sorted(ids))]

    best: Optional[MultiChainStructure] = None
    for order in candidates:
        built = _build_for_order(chain_data, order, cww_pairs)
        if best is None:
            best = built
            continue
        # Maximize nested (== minimize crossing); tie-break lexicographic order.
        if len(built.nested_pairs) > len(best.nested_pairs):
            best = built
        elif len(built.nested_pairs) == len(best.nested_pairs) and (
            built.order < best.order
        ):
            best = built
    assert best is not None
    return best


def assemble(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    cif_file: str,
    output_dir: str,
    chains: Optional[List[str]] = None,
    auto_order: bool = True,
    model_id: Optional[int] = None,
    quiet: bool = False,
) -> Optional[MultiChainStructure]:
    """Extract, concatenate, and partition the selected RNA chains.

    Args:
        cif_file: Path to an mmCIF structure.
        output_dir: Directory for the FR3D basepair output.
        chains: Chain ids to include.  ``None`` selects all RNA chains.
        auto_order: If True, search for the best concatenation order among the
            selected chains.  If False, use ``chains`` as the literal order.
        model_id: Model to extract (default: first model).
        quiet: Suppress FR3D and progress output.

    Returns:
        A MultiChainStructure, or None if extraction failed.
    """
    if chains is None:
        chains = list_rna_chains(cif_file, model_id=model_id)
    if not chains:
        if not quiet:
            print("[multichain] no RNA chains found")
        return None

    chain_data = extract_chains(cif_file, chains, model_id=model_id, quiet=quiet)
    if len(chain_data) < 1:
        return None

    basepair_file = fr3d_utils.run_fr3d(cif_file, output_dir, quiet=quiet)
    if not basepair_file:
        if not quiet:
            print("[multichain] FR3D failed")
        return None

    # Derive the extracted model from the unit ids so NMR pair lines from other
    # models are filtered out.
    sample = next((uid for c in chain_data for uid in c.unit_id_to_local), None)
    target_model = sample.split("|")[1] if sample and "|" in sample else None
    cww_pairs = read_cww_pairs(basepair_file, target_model=target_model)

    if auto_order:
        result = choose_order(chain_data, cww_pairs)
    else:
        order = [c.chain_id for c in chain_data]  # preserve requested order
        result = _build_for_order(chain_data, order, cww_pairs)

    if not quiet:
        print(f"[multichain] order={'+'.join(result.order)} {result.counts()}")
    return result


# ---------------------------------------------------------------------------
# Combined-SVG post-processing: junction break, per-chain 5'/3', overlay arcs
# ---------------------------------------------------------------------------

_NT_G_RE = re.compile(
    r'<g><title>(\d+)[^<]*</title><text\s+x="([-\d.]+)"\s+y="([-\d.]+)"\s+'
    r'class="\w+"\s*>([^<]+)</text></g>'
)
_LINE_RE = re.compile(
    r'<line\s+x1="([-\d.]+)"\s+y1="([-\d.]+)"\s+x2="([-\d.]+)"\s+y2="([-\d.]+)"'
    r'\s+class="([^"]*)"\s*/>'
)
# Generic 5'/3' terminus markers R2DT emits for the whole (concatenated) strand.
_MARKER_RE = re.compile(
    r"<g><title>\d+[^<]*</title><text[^>]*>[53]&#8242;?</text></g>"
    r"|<g><title>\d+[^<]*</title><text[^>]*>[53]'</text></g>"
)


def _parse_nt_centres(content: str) -> Dict[int, Tuple[float, float]]:
    """Map 1-based nucleotide index → (x, y) centre from an R2DT SVG.

    Also captures the synthetic 5'/3' markers (index 0 and N+1); callers index
    only the real nucleotides (1..N)."""
    out: Dict[int, Tuple[float, float]] = {}
    for m in _NT_G_RE.finditer(content):
        out[int(m.group(1))] = (float(m.group(2)), float(m.group(3)))
    return out


def _nearest_index(x: float, y: float, centres: Dict[int, Tuple[float, float]]) -> int:
    best_idx, best_d = -1, float("inf")
    for idx, (cx, cy) in centres.items():
        d = math.hypot(cx - x, cy - y)
        if d < best_d:
            best_idx, best_d = idx, d
    return best_idx


def _terminus_label(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    centres: Dict[int, Tuple[float, float]],
    at_idx: int,
    neighbour_idx: int,
    tag: str,
    chain_id: str,
    colour: str,
) -> str:
    """A chain-tagged 5'/3' label, offset outward from the terminal nucleotide
    along the direction away from its sequence neighbour."""
    x, y = centres[at_idx]
    nx, ny = centres.get(neighbour_idx, (x, y))
    dx, dy = x - nx, y - ny
    length = math.hypot(dx, dy) or 1.0
    ox = x + dx / length * 14.0
    oy = y + dy / length * 14.0
    return (
        f'<text x="{ox:.2f}" y="{oy:.2f}" class="mc-terminus" fill="{colour}" '
        f'font-size="10" font-weight="bold" text-anchor="middle">'
        f"{chain_id}{tag}</text>"
    )


def _overlay_arc(a: Tuple[float, float], b: Tuple[float, float]) -> str:
    """A dashed quadratic-Bézier arc for a crossing (overlay) base pair."""
    x0, y0 = a
    x1, y1 = b
    mx, my = (x0 + x1) / 2, (y0 + y1) / 2
    chord = math.hypot(x1 - x0, y1 - y0) or 1.0
    bulge = max(chord * 0.3, 15.0)
    cx = mx - (y1 - y0) / chord * bulge
    cy = my + (x1 - x0) / chord * bulge
    return (
        f'<path d="M {x0:.1f},{y0:.1f} Q {cx:.1f},{cy:.1f} {x1:.1f},{y1:.1f}" '
        f'fill="none" stroke="{_OVERLAY_STROKE}" stroke-width="1" '
        f'stroke-dasharray="4,3" opacity="0.85" class="mc-crossing-bp"/>'
    )


def postprocess_combined_svg(  # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals
    svg_path: str,
    boundaries: List[Tuple[str, int, int]],
    nested_pairs: List[Tuple[int, int]],
    crossing_pairs: List[Tuple[int, int]],
    quiet: bool = False,
) -> None:
    """Turn a concatenated-strand SVG into a multi-chain diagram in place.

    - breaks the phantom backbone bond at each chain junction,
    - replaces the single 5'/3' pair with chain-tagged termini per chain,
    - draws crossing inter-chain pairs as dashed overlay arcs.
    """
    path = Path(svg_path)
    content = path.read_text()
    centres = _parse_nt_centres(content)
    if not centres:
        if not quiet:
            print(f"[multichain] no nucleotides parsed from {svg_path}; skipping")
        return

    # 1. Break the phantom backbone bond at each internal chain boundary.
    #    Boundary end `b` (0-based, exclusive) ⇒ 1-based junction nts (b, b+1).
    #    When the junction nucleotides are *also* a base pair (an antiparallel
    #    duplex pairs chain i's 3' with chain i+1's 5'), the backbone line and
    #    the base-pair line coincide — remove only the backbone (prefer the
    #    de-emphasised "gray" line) and keep the pair.
    nested_1based = {frozenset((i + 1, j + 1)) for i, j in nested_pairs}
    removed = 0
    for _, _, b in boundaries[:-1]:
        jkey = frozenset((b, b + 1))
        bp_here = 1 if jkey in nested_1based else 0
        matches: List[Tuple[str, bool]] = []
        for m in _LINE_RE.finditer(content):
            x1, y1, x2, y2 = map(float, m.groups()[:4])
            i = _nearest_index(x1, y1, centres)
            j = _nearest_index(x2, y2, centres)
            if frozenset((i, j)) == jkey:
                matches.append((m.group(0), m.group(5) == "gray"))
        # Keep bp_here line(s); remove the rest, dropping "gray" backbone first.
        matches.sort(key=lambda t: 0 if t[1] else 1)
        for text, _is_gray in matches[: max(0, len(matches) - bp_here)]:
            content = content.replace(text, "", 1)
            removed += 1

    # 2. Drop the generic strand 5'/3' markers, add chain-tagged termini.
    content = _MARKER_RE.sub("", content)
    labels: List[str] = []
    for ci, (chain_id, start, end) in enumerate(boundaries):
        colour = _CHAIN_PALETTE[ci % len(_CHAIN_PALETTE)]
        first, last = start + 1, end  # 1-based first / last nucleotide
        labels.append(
            _terminus_label(
                centres, first, min(first + 1, last), "5′", chain_id, colour
            )
        )
        labels.append(
            _terminus_label(centres, last, max(last - 1, first), "3′", chain_id, colour)
        )

    # 3. Overlay arcs for crossing pairs (0-based ⇒ 1-based nt indices).
    arcs: List[str] = []
    for i, j in crossing_pairs:
        a = centres.get(i + 1)
        b = centres.get(j + 1)
        if a and b:
            arcs.append(_overlay_arc(a, b))

    inject = (
        f'<g class="mc-overlay">{"".join(arcs)}</g>'
        f'<g class="mc-termini">{"".join(labels)}</g>'
    )
    content = content.replace("</svg>", inject + "</svg>", 1)
    path.write_text(content)

    if not quiet:
        print(
            f"[multichain] post-processed {path.name}: "
            f"{removed} junction bond(s) broken, {len(labels)} termini, "
            f"{len(arcs)} overlay arc(s)"
        )
