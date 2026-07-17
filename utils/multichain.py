# pylint: disable=too-many-lines
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
import random
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from utils import fr3d as fr3d_utils

# Colour-blind-safe palette used to tag each chain's 5'/3' termini.
_CHAIN_PALETTE = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02"]
# Overlay style for crossing (pseudoknot / inter-chain kissing) pairs.
_OVERLAY_STROKE = "#e6194b"

# Reference/model base-pair diff colours (approach B: draw the model's pairs on
# the reference layout).  Matched pairs are de-emphasised so the eye lands on
# the differences.
_DIFF_MATCHED = "#999999"  # agreement
_DIFF_LOST = "#4363d8"  # in reference, missing from model (false negative)
_DIFF_ADDED = "#e6194b"  # in model only (false positive)

# Above this many chains the N! order search is skipped in favour of a greedy
# adjacency ordering (see _greedy_order).  Real RNA complexes sit well below
# this; the cap only guards against pathological inputs.
_MAX_BRUTEFORCE_CHAINS = 6

_CWW_INTERACTIONS = ("cWW", "cWw", "cwW")

# A Leontis–Westhof base-pair family: cis/trans + two edges from {W,H,S}.
# Excludes FR3D "near" pairs (n-prefixed), stacking (s33/s55…), bifurcated, etc.
_LW_FAMILY_RE = re.compile(r"^[ct][whs][whs]$", re.IGNORECASE)


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
    auth_of: List[Optional[int]] = field(default_factory=list)  # author resnum
    inter_pairs: List[Tuple[int, int]] = field(default_factory=list)
    # Every base pair (all Leontis–Westhof families) as (i, j, family), for INF.
    all_pairs: List[Tuple[int, int, str]] = field(default_factory=list)

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
            # Skip crystal-symmetry contacts (a real residue paired with a
            # symmetry mate): they are inter-copy packing, not intramolecular
            # pairs, and would otherwise render as spurious lines once the
            # symmetry endpoint resolves back to its asymmetric-unit position.
            if fr3d_utils.is_symmetry_mate(unit1) or fr3d_utils.is_symmetry_mate(unit2):
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


def read_all_bp(
    basepair_file: str, target_model: Optional[str] = None
) -> List[Tuple[str, str, str]]:
    """Read every base pair (all LW families) as (unit_id1, family, unit_id2).

    Unlike :func:`read_cww_pairs` this keeps non-canonical families too, which
    the INF metric needs.  Near/stacking/other annotations are skipped, and
    pairs are deduplicated by unordered unit pair (keeping the first family).
    """
    seen = set()
    pairs: List[Tuple[str, str, str]] = []
    with open(basepair_file, "r") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            unit1, family, unit2 = parts[0], parts[1].strip(), parts[2]
            if not _LW_FAMILY_RE.match(family):
                continue
            # Skip crystal-symmetry contacts (see read_cww_pairs).
            if fr3d_utils.is_symmetry_mate(unit1) or fr3d_utils.is_symmetry_mate(unit2):
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
            pairs.append((unit1, family, unit2))
    return pairs


def _map_all_bp(
    all_bp: List[Tuple[str, str, str]], combined: Dict[str, int]
) -> List[Tuple[int, int, str]]:
    """Map (unit1, family, unit2) triples to 0-based (i, j, family) positions."""
    out: List[Tuple[int, int, str]] = []
    used = set()
    for unit1, family, unit2 in all_bp:
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
        key = (min(pos1, pos2), max(pos1, pos2))
        if key in used:
            continue
        used.add(key)
        out.append((key[0], key[1], family))
    return out


def _combined_map(chain_data: List[ChainSeq], order: List[str]) -> Tuple[
    Dict[str, int],
    List[str],
    List[Tuple[str, int, int]],
    str,
    List[Optional[int]],
]:
    """Build the concatenated sequence, unit-id → global-position map,
    per-position chain list, chain boundaries, and per-position author residue
    number (parsed from the unit id) for a given order."""
    by_id = {c.chain_id: c for c in chain_data}
    combined: Dict[str, int] = {}
    chain_of: List[str] = []
    boundaries: List[Tuple[str, int, int]] = []
    sequence_parts: List[str] = []
    auth_of: List[Optional[int]] = []
    offset = 0
    for chain_id in order:
        chain = by_id[chain_id]
        local_to_auth: Dict[int, Optional[int]] = {}
        for unit_id, local in chain.unit_id_to_local.items():
            combined[unit_id] = offset + local
            parts = unit_id.split("|")
            try:
                local_to_auth[local] = int(parts[4]) if len(parts) > 4 else None
            except ValueError:
                local_to_auth[local] = None
        chain_of.extend([chain_id] * len(chain.sequence))
        auth_of.extend(local_to_auth.get(i) for i in range(len(chain.sequence)))
        boundaries.append((chain_id, offset, offset + len(chain.sequence)))
        sequence_parts.append(chain.sequence)
        offset += len(chain.sequence)
    return combined, chain_of, boundaries, "".join(sequence_parts), auth_of


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
    combined, chain_of, boundaries, sequence, auth_of = _combined_map(chain_data, order)
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
        auth_of=auth_of,
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

    # Map the full base-pair set (all LW families) onto the chosen order's
    # label space so the INF metric can compare canonical vs non-canonical.
    combined, _, _, _, _ = _combined_map(chain_data, result.order)
    all_bp = read_all_bp(basepair_file, target_model=target_model)
    result.all_pairs = _map_all_bp(all_bp, combined)

    if not quiet:
        print(f"[multichain] order={'+'.join(result.order)} {result.counts()}")
    return result


def partition_components(
    cif_file: str,
    output_dir: str,
    model_id: Optional[int] = None,
    quiet: bool = False,
) -> Optional[List[List[str]]]:
    """Partition every RNA chain in the structure into connected components by
    inter-chain base pairing.

    Two chains land in the same component iff they share at least one
    inter-chain pair, transitively; a chain with no inter-chain pairs to
    anything is its own singleton component. This is the basis for deciding
    what a reference structure's viewer should show together: chains that
    base-pair with each other belong in one combined view, chains that don't
    interact with anything are independent and browsable one at a time.

    Returns a list of components (each a list of chain ids), or ``None`` if
    extraction/FR3D failed. An empty structure (no RNA chains) returns ``[]``.
    """
    chains = list_rna_chains(cif_file, model_id=model_id)
    if not chains:
        return []
    if len(chains) == 1:
        return [chains]

    result = assemble(
        cif_file,
        output_dir,
        chains=None,
        auto_order=True,
        model_id=model_id,
        quiet=quiet,
    )
    if result is None:
        return None

    parent = {c: c for c in result.order}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for i, j in result.inter_pairs:
        union(result.chain_of[i], result.chain_of[j])

    groups: Dict[str, List[str]] = {}
    for c in result.order:
        groups.setdefault(find(c), []).append(c)
    return list(groups.values())


# ---------------------------------------------------------------------------
# INF (Interaction Network Fidelity, Parisien et al. 2009)
# ---------------------------------------------------------------------------


def _is_canonical(family: str) -> bool:
    """Canonical = cis Watson-Crick (cWW).  Everything else is non-canonical."""
    return family.upper() == "CWW"


def _inf_score(ref: set, model: set) -> Dict[str, Optional[float]]:
    """INF for one interaction set: sqrt(PPV * STY).

    TP = ref ∩ model, FP = model \\ ref, FN = ref \\ model.  ``inf`` is None
    when neither structure has any interaction of this type.
    """
    tp = len(ref & model)
    fp = len(model - ref)
    fn = len(ref - model)
    if tp == 0:
        inf: Optional[float] = None if (not ref and not model) else 0.0
        ppv = sty = None if (not ref and not model) else 0.0
    else:
        ppv = tp / (tp + fp)
        sty = tp / (tp + fn)
        inf = math.sqrt(ppv * sty)
    return {"inf": inf, "ppv": ppv, "sty": sty, "tp": tp, "fp": fp, "fn": fn}


def _split_families(pairs: List[Tuple[int, int, str]]) -> Tuple[set, set, set]:
    """Return (canonical, non_canonical, all) sets of (i, j) position pairs."""
    canonical, non_canonical, everything = set(), set(), set()
    for i, j, family in pairs:
        key = (i, j)
        everything.add(key)
        (canonical if _is_canonical(family) else non_canonical).add(key)
    return canonical, non_canonical, everything


def compute_inf(
    ref_pairs: List[Tuple[int, int, str]],
    model_pairs: List[Tuple[int, int, str]],
) -> Dict[str, Dict[str, Optional[float]]]:
    """Interaction Network Fidelity of ``model`` against ``ref``.

    Returns INF for three interaction subsets:
    ``wc`` (canonical/cWW), ``nwc`` (non-canonical), and ``all`` (both).
    """
    r_wc, r_nwc, r_all = _split_families(ref_pairs)
    m_wc, m_nwc, m_all = _split_families(model_pairs)
    return {
        "wc": _inf_score(r_wc, m_wc),
        "nwc": _inf_score(r_nwc, m_nwc),
        "all": _inf_score(r_all, m_all),
    }


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
    # Use an inline ``style`` (not a ``fill`` attribute): R2DT's SVG carries a
    # bare ``text { fill: black; ... }`` CSS rule that overrides presentation
    # attributes, so only ``style`` wins.
    return (
        f'<text x="{ox:.2f}" y="{oy:.2f}" class="mc-terminus" '
        f'style="fill:{colour};font-size:10px;font-weight:bold;'
        f'text-anchor:middle">{chain_id}{tag}</text>'
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


# ---------------------------------------------------------------------------
# Reference/model base-pair comparison (approach B: model pairs on ref layout)
# ---------------------------------------------------------------------------

_CROSSING_PATH_RE = re.compile(r'<path[^>]*class="mc-crossing-bp"[^>]*/>')


def diff_pairs(
    ref_pairs: List[Tuple[int, int]], model_pairs: List[Tuple[int, int]]
) -> Tuple[List[Tuple[int, int]], List[Tuple[int, int]], List[Tuple[int, int]]]:
    """Compare two pair sets over the same coordinates.

    Returns ``(matched, lost, added)``: pairs in both, pairs only in the
    reference (missed by the model), and pairs only in the model (spurious).
    """
    ref = {tuple(p) for p in ref_pairs}
    model = {tuple(p) for p in model_pairs}
    return sorted(ref & model), sorted(ref - model), sorted(model - ref)


# ---------------------------------------------------------------------------
# 3D superposition (approach C: pre-align the predicted model onto the
# reference so the viewer can overlay both structures without an in-browser
# transform — the pinned pdbe-molstar CDN bundle exposes no transform API).
# ---------------------------------------------------------------------------


def _load_structure(path: str):
    """Parse a ``.cif``/``.pdb`` file into a BioPython structure.

    Uses author chain ids and residue numbers (MMCIFParser ``auth_*`` defaults,
    which PDBParser already does) so ``(chain, auth)`` lookups line up with the
    FR3D unit ids that produced ``chain_of``/``auth_of``.
    """
    # pylint: disable=import-outside-toplevel
    import warnings

    from Bio.PDB import MMCIFParser, PDBParser

    stem = Path(path).stem
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if Path(path).suffix.lower() == ".cif":
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)
        return parser.get_structure(stem, path)


def _atom_coords_by_residue(structure, atom_names=("C1'", "P")):
    """Map ``(chain_id, auth_seq_id) -> coord`` for the first model.

    Picks the first present atom in ``atom_names`` per residue (C1' for RNA,
    falling back to the phosphate). Heteroatoms/waters are skipped.
    """
    # pylint: disable=import-outside-toplevel
    import numpy as np

    model = next(iter(structure))
    out: Dict[Tuple[str, int], "np.ndarray"] = {}
    for chain in model:
        for residue in chain:
            if residue.id[0] != " ":  # HETATM / water
                continue
            key = (chain.id, residue.id[1])
            if key in out:
                continue
            for name in atom_names:
                if residue.has_id(name):
                    out[key] = np.asarray(residue[name].coord, dtype=float)
                    break
    return out


# pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals
def superpose_model_onto_reference(
    ref_file: str,
    model_file: str,
    ref_index: List[Optional[Tuple[str, int]]],
    model_index: List[Optional[Tuple[str, int]]],
    out_path: str,
    *,
    atom_names: Tuple[str, ...] = ("C1'", "P"),
    min_atoms: int = 3,
    quiet: bool = False,
) -> Optional[Tuple[float, int]]:
    """Rigid-superpose the model onto the reference and write the aligned model.

    ``ref_index[i]`` / ``model_index[i]`` give ``(chain, auth)`` for label
    position *i* in each structure (``None`` where unresolved). Both are
    co-indexed (approach B: same sequence, same chain order), so the fit runs on
    the positions resolved in *both* — an exact atom correspondence, no alignment.
    The resulting rotation+translation is applied to **every** model atom and the
    structure is written to ``out_path``.

    Returns ``(rmsd, n_atoms_used)``, or ``None`` if fewer than ``min_atoms``
    correspondences exist (caller should load the model unaligned instead).
    """
    # pylint: disable=import-outside-toplevel
    import numpy as np
    from Bio.PDB import MMCIFIO, PDBIO
    from Bio.SVDSuperimposer import SVDSuperimposer

    model_struct = _load_structure(model_file)
    ref_coords = _atom_coords_by_residue(_load_structure(ref_file), atom_names)
    model_coords = _atom_coords_by_residue(model_struct, atom_names)

    fixed, moving = [], []
    for rk, mk in zip(ref_index, model_index):
        if rk is None or mk is None:
            continue
        rc = ref_coords.get(rk)
        mc = model_coords.get(mk)
        if rc is None or mc is None:
            continue
        fixed.append(rc)
        moving.append(mc)

    if len(fixed) < min_atoms:
        if not quiet:
            print(
                f"[multichain] superpose: only {len(fixed)} matched atoms "
                f"(< {min_atoms}); loading model unaligned"
            )
        return None

    sup = SVDSuperimposer()
    sup.set(np.asarray(fixed), np.asarray(moving))
    sup.run()
    rot, tran = sup.get_rotran()
    rms = sup.get_rms()

    for atom in model_struct.get_atoms():
        atom.coord = np.dot(atom.coord, rot) + tran

    io = MMCIFIO() if str(out_path).lower().endswith(".cif") else PDBIO()
    io.set_structure(model_struct)
    io.save(str(out_path))
    if not quiet:
        print(
            f"[multichain] superposed model onto reference: {len(fixed)} atoms, "
            f"RMSD {rms:.2f} Å -> {Path(out_path).name}"
        )
    return rms, len(fixed)


def write_structure_chains(
    src_path: str,
    chains: List[str],
    out_path: str,
    *,
    quiet: bool = False,
) -> bool:
    """Write a copy of ``src_path`` keeping only the given author chains.

    A deposited reference entry often carries more than the RNA under analysis
    (extra RNA copies, ligands, whole protein chains). The 3D pane colours every
    non-selected atom with the reference base colour, so those extras render as
    additional structures beside the one the model is superimposed on. Keeping
    only the analysed chain(s) leaves a single reference structure in the view.

    Returns ``True`` if a filtered file was written, ``False`` on failure
    (caller should fall back to copying the original).
    """
    # pylint: disable=import-outside-toplevel
    from Bio.PDB import MMCIFIO, PDBIO, Select

    keep = {str(c) for c in chains}
    if not keep:
        return False
    try:
        structure = _load_structure(src_path)
    except Exception:  # pylint: disable=broad-except
        return False

    class _ChainSelect(Select):  # type: ignore[misc]
        def accept_chain(self, chain):  # noqa: N802
            """Keep only chains in ``keep``."""
            return 1 if chain.id in keep else 0

    io = MMCIFIO() if str(out_path).lower().endswith(".cif") else PDBIO()
    io.set_structure(structure)
    try:
        io.save(str(out_path), _ChainSelect())
    except Exception:  # pylint: disable=broad-except
        return False
    if not quiet:
        print(
            f"[multichain] wrote reference with chains {sorted(keep)} "
            f"-> {Path(out_path).name}"
        )
    return True


def simulate_model_pairs(
    ref_pairs: List[Tuple[int, int]],
    n_positions: int,
    *,
    seed: int = 0,
    remove_frac: float = 0.2,
    n_add: int = 2,
) -> List[Tuple[int, int]]:
    """Perturb a reference pair set into a synthetic "model" pair set.

    TEMPORARY testing aid for when no real predicted structure is available:
    drops a fraction of the reference pairs (simulated false negatives) and
    invents a few pairs between currently-unpaired positions (false positives).
    Deterministic for a given seed.
    """
    rng = random.Random(seed)
    kept = [tuple(p) for p in ref_pairs if rng.random() > remove_frac]
    used = {x for pair in kept for x in pair}
    free = [i for i in range(n_positions) if i not in used]
    rng.shuffle(free)
    added: List[Tuple[int, int]] = []
    while len(free) >= 2 and len(added) < n_add:
        a, b = free.pop(), free.pop()
        added.append((min(a, b), max(a, b)))
    return sorted(kept + added)


def _pair_glyph(
    centres: Dict[int, Tuple[float, float]],
    idx1: int,
    idx2: int,
    colour: str,
    *,
    dashed: bool,
    width: float,
) -> str:
    """A base-pair glyph (gentle Bézier arc) between two 1-based nucleotides."""
    if idx1 not in centres or idx2 not in centres:
        return ""
    x0, y0 = centres[idx1]
    x1, y1 = centres[idx2]
    chord = math.hypot(x1 - x0, y1 - y0) or 1.0
    bulge = chord * 0.15  # near-straight for helical pairs, arced for long ones
    cx = (x0 + x1) / 2 - (y1 - y0) / chord * bulge
    cy = (y0 + y1) / 2 + (x1 - x0) / chord * bulge
    dash = ' stroke-dasharray="4,3"' if dashed else ""
    return (
        f'<path d="M {x0:.1f},{y0:.1f} Q {cx:.1f},{cy:.1f} {x1:.1f},{y1:.1f}" '
        f'fill="none" stroke="{colour}" stroke-width="{width:.1f}"{dash} '
        f'opacity="0.9" class="mc-diff-bp"/>'
    )


def _remove_ref_pair_line(
    content: str, centres: Dict[int, Tuple[float, float]], idx1: int, idx2: int
) -> str:
    """Remove the engine-drawn base-pair line between two nucleotides.

    Guards against touching backbone lines by only removing lines whose
    endpoints join non-consecutive nucleotides.
    """
    if abs(idx1 - idx2) <= 1:
        return content  # adjacent: could be backbone — leave it
    target = frozenset((idx1, idx2))
    for m in _LINE_RE.finditer(content):
        x1, y1, x2, y2 = map(float, m.groups()[:4])
        i = _nearest_index(x1, y1, centres)
        j = _nearest_index(x2, y2, centres)
        if frozenset((i, j)) == target:
            return content.replace(m.group(0), "", 1)
    return content


def _diff_legend(n_matched: int, n_lost: int, n_added: int) -> str:
    rows = [
        (_DIFF_MATCHED, f"matched: {n_matched}"),
        (_DIFF_LOST, f"missing in model: {n_lost}"),
        (_DIFF_ADDED, f"model-only: {n_added}"),
    ]
    # Inline ``style`` beats R2DT's bare ``text { ... text-anchor: middle }``
    # rule; without it the labels centre on x and spill off the left edge.
    texts = "".join(
        f'<text x="6" y="{10 + i * 10}" style="fill:{colour};font-size:7px;'
        f'font-weight:bold;text-anchor:start">{label}</text>'
        for i, (colour, label) in enumerate(rows)
    )
    return f'<g class="mc-diff-legend">{texts}</g>'


def render_model_panel(
    ref_svg_path: str,
    out_svg_path: str,
    ref_pairs: List[Tuple[int, int]],
    model_pairs: List[Tuple[int, int]],
    quiet: bool = False,
) -> None:
    """Render a model panel on the reference layout, coloured by pair diff.

    Reuses the reference SVG's coordinates, backbone, nucleotides and termini
    verbatim (approach B), strips the reference's own base-pair rendering, and
    redraws every pair by diff category: matched (grey), missing-in-model
    (blue, dashed), model-only (red).
    """
    content = Path(ref_svg_path).read_text()
    centres = _parse_nt_centres(content)
    if not centres:
        if not quiet:
            print(f"[multichain] no nucleotides in {ref_svg_path}; skipping model")
        return

    matched, lost, added = diff_pairs(ref_pairs, model_pairs)

    # Strip the reference's base-pair rendering (crossing arcs + nested lines).
    content = _CROSSING_PATH_RE.sub("", content)
    for i, j in ref_pairs:
        content = _remove_ref_pair_line(content, centres, i + 1, j + 1)

    glyphs = [
        _pair_glyph(centres, i + 1, j + 1, _DIFF_MATCHED, dashed=False, width=1.0)
        for i, j in matched
    ]
    glyphs += [
        _pair_glyph(centres, i + 1, j + 1, _DIFF_LOST, dashed=True, width=1.2)
        for i, j in lost
    ]
    glyphs += [
        _pair_glyph(centres, i + 1, j + 1, _DIFF_ADDED, dashed=False, width=1.6)
        for i, j in added
    ]

    inject = (
        f'<g class="mc-diff">{"".join(glyphs)}</g>'
        f"{_diff_legend(len(matched), len(lost), len(added))}"
    )
    content = content.replace("</svg>", inject + "</svg>", 1)
    Path(out_svg_path).write_text(content)

    if not quiet:
        print(
            f"[multichain] model panel {Path(out_svg_path).name}: "
            f"{len(matched)} matched, {len(lost)} missing, {len(added)} model-only"
        )
