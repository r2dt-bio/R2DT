"""Shared building blocks for the layered base-pair notation scripts.

Used by both layered_basepairs.py (pairs read from an FR3D TSV) and
standalone_lbn_script.py (pairs read from a DNATCO-extended mmCIF, plus an
iCn3D 3D overlay). Each script reads its own pair list, then calls
build_notation(per_chain, pairs, chains, ...) from here to produce the
notation. Anything that is truly shared lives in this module; anything specific
to the input source or to 3D drawing stays in the caller.

Contents:
  - constants:  BRACKETS, OPENERS, CLOSE_TO_OPEN, DIRECTED, SEP, BLOCK,
                IDENTITY, NT, WC_PAIRS
  - LW utilities:                        _flip_lw, _is_canonical
  - CIF readers (sequence + numbering):  _loop_columns, _read_loop_rows,
                                         _read_poly_seq_scheme,
                                         _read_atom_site, read_residues
  - layout / render / parse:             _cross, _pack, render, parse
  - row-label helpers:                   _label, _span, _strand_token,
                                         _layer_label, _family_of,
                                         _res_id, _compact_line
  - top-level build:                     build_notation, _roundtrip, _wrap
"""

import shlex
from pathlib import Path


# --------------------------------------------------------------------------- #
#  Constants                                                                   #
# --------------------------------------------------------------------------- #

# Bracket classes for crossing pairs within one layer (pseudoknot levels).
BRACKETS = ["()", "[]", "{}", "<>", "Aa", "Bb", "Cc", "Dd"]
OPENERS = set("([{<ABCD")
CLOSE_TO_OPEN = {")": "(", "]": "[", "}": "{", ">": "<",
                 "a": "A", "b": "B", "c": "C", "d": "D"}

# Fixed layer order: the 18 directed LW types, cWW first.
DIRECTED = ["cWW", "cWH", "cWS", "cHW", "cHH", "cHS", "cSW", "cSH", "cSS",
            "tWW", "tWH", "tWS", "tHW", "tHH", "tHS", "tSW", "tSH", "tSS"]

SEP = "&"           # chain boundary marker
BLOCK = 100         # default columns per block when --block is given with no value
IDENTITY = "1_555"  # crystallographic identity operator (asymmetric-unit copy)

# Residue name -> one sequence letter (RNA and DNA); modified residues fall back
# to their last character, lowercased, so columns stay aligned.
NT = {"A": "A", "C": "C", "G": "G", "U": "U", "T": "T",
      "DA": "A", "DC": "C", "DG": "G", "DT": "T", "DU": "U"}

# Canonical Watson-Crick base pairs (used to split cWW into WC vs cWW layers).
WC_PAIRS = ({"A", "U"}, {"G", "C"}, {"A", "T"})


# --------------------------------------------------------------------------- #
#  LW utilities                                                                #
# --------------------------------------------------------------------------- #

def _flip_lw(lw: str) -> str:
    """Flip an LW code direction (cWH <-> cHW) when a pair is reordered."""
    if lw and len(lw) == 3 and lw[0] in "ct":
        return lw[0] + lw[2] + lw[1]
    return lw


def _is_canonical(lw: str, base_a: str, base_b: str) -> bool:
    """A canonical Watson-Crick pair: cWW geometry AND A-U/G-C/A-T bases. cWW is
    a *family*, not 'canonical' -- a cWW U-U or U-G wobble is non-canonical."""
    return lw == "cWW" and {base_a.upper(), base_b.upper()} in WC_PAIRS


# --------------------------------------------------------------------------- #
#  CIF readers (sequence and residue numbering)                                #
# --------------------------------------------------------------------------- #

def _loop_columns(lines: list[str], prefix: str, start: int = 0) -> tuple[dict, int]:
    """Find a `prefix` loop, return ({column_name: index}, index of first data row)."""
    i = start
    while i < len(lines) and not lines[i].startswith(prefix):
        i += 1
    cols: list[str] = []
    while i < len(lines) and lines[i].startswith(prefix):
        cols.append(lines[i].strip().split(".")[1])
        i += 1
    return {n: k for k, n in enumerate(cols)}, i


def _read_loop_rows(lines: list[str], i: int, n_cols: int) -> list[list[str]]:
    """Walk forward from `i` collecting tokens until each row has `n_cols`
    tokens. mmCIF allows long rows to wrap across several text lines, so a
    naive line-by-line read mis-counts. Stops at the next loop_, category
    header, comment, or blank line. Quote-aware via shlex."""
    rows: list[list[str]] = []
    buf: list[str] = []
    while i < len(lines):
        s = lines[i].strip()
        if s == "" or s.startswith(("#", "loop_", "_", "data_")):
            if buf:
                rows.append(buf)
            break
        try:
            buf.extend(shlex.split(s))
        except ValueError:
            buf.extend(s.split())
        while len(buf) >= n_cols:
            rows.append(buf[:n_cols])
            buf = buf[n_cols:]
        i += 1
    return rows


def _read_poly_seq_scheme(lines: list[str], chains: list[str]) -> dict[str, list]:
    """Per-chain [(number, letter), ...] from pdbx_poly_seq_scheme."""
    idx, i = _loop_columns(lines, "_pdbx_poly_seq_scheme.")
    c_mon, c_num, c_strand = idx["mon_id"], idx["pdb_seq_num"], idx["pdb_strand_id"]
    per_chain: dict[str, list] = {c: [] for c in chains}
    while i < len(lines):
        row = lines[i].strip()
        if row.startswith(("#", "_", "loop_")):
            break
        f = row.split()
        if len(f) >= len(idx) and f[c_strand] in per_chain:
            mon = f[c_mon]
            per_chain[f[c_strand]].append((int(f[c_num]), NT.get(mon, mon[-1].lower())))
        i += 1
    return per_chain


def _read_atom_site(lines: list[str], chains: list[str]) -> dict[str, list]:
    """Per-chain [(number, letter), ...] from _atom_site (model 1, polymer only).
    Fallback when pdbx_poly_seq_scheme is absent (e.g. minimal modelling CIFs)."""
    idx, i = _loop_columns(lines, "_atom_site.")
    c_ch, c_num, c_comp = idx["auth_asym_id"], idx["auth_seq_id"], idx["label_comp_id"]
    c_lab = idx["label_seq_id"]
    c_model = idx.get("pdbx_PDB_model_num")
    per_chain: dict[str, list] = {c: [] for c in chains}
    seen: set[tuple[str, str]] = set()
    while i < len(lines):
        row = lines[i].strip()
        if row.startswith(("#", "loop_", "_")):
            break
        f = row.split()
        if len(f) >= len(idx):
            if (c_model is not None and f[c_model] != "1") or f[c_lab] == ".":
                i += 1
                continue
            ch = f[c_ch]
            key = (ch, f[c_num])
            if ch in per_chain and key not in seen:
                seen.add(key)
                comp = f[c_comp]
                per_chain[ch].append((int(f[c_num]), NT.get(comp, comp[-1].lower())))
        i += 1
    return per_chain


def read_residues(cif: Path, chains: list[str]) -> dict[str, list]:
    """Per-chain [(number, letter), ...] from pdbx_poly_seq_scheme, falling back
    to _atom_site when that category is absent."""
    lines = cif.read_text().splitlines()
    has_scheme = any(l.startswith("_pdbx_poly_seq_scheme.") for l in lines)
    return (_read_poly_seq_scheme(lines, chains) if has_scheme
            else _read_atom_site(lines, chains))


# --------------------------------------------------------------------------- #
#  Layout, render and parse                                                    #
# --------------------------------------------------------------------------- #

def _cross(p: tuple[int, int], q: tuple[int, int]) -> bool:
    (i, j), (k, l) = p, q
    return i < k < j < l or k < i < l < j


def _pack(pairs: list) -> list[list]:
    """Pack same-type overflow pairs into as few layers as possible, each
    residue used at most once per layer."""
    buckets: list[list] = []
    for p in sorted(pairs):
        Ra, Rb, _ = p
        for bucket in buckets:
            keys = {k for x, y, _ in bucket for k in (x, y)}
            if Ra not in keys and Rb not in keys:
                bucket.append(p)
                break
        else:
            buckets.append([p])
    return buckets


def render(layer: list, key_to_col: dict, width: int, sep_cols: set[int]) -> str:
    """Render one layer to a dot-bracket line over the multi-strand ruler."""
    chars = [SEP if i in sep_cols else "." for i in range(width)]
    colpairs = [tuple(sorted((key_to_col[Ra], key_to_col[Rb]))) for Ra, Rb, _ in layer]
    classed: list[tuple[tuple[int, int], int]] = []
    for (i, j) in sorted(colpairs):
        cls = 0
        while any(c == cls and _cross((i, j), pq) for pq, c in classed):
            cls += 1
        classed.append(((i, j), cls))
        chars[i], chars[j] = BRACKETS[cls][0], BRACKETS[cls][1]
    return "".join(chars)


def parse(s: str, col_to_key: dict) -> set[frozenset]:
    """Read one rendered line back into {frozenset(Ra, Rb), ...}. '&' is skipped."""
    stacks: dict[str, list[int]] = {}
    out: set[frozenset] = set()
    for col, ch in enumerate(s):
        if ch in OPENERS:
            stacks.setdefault(ch, []).append(col)
        elif ch in CLOSE_TO_OPEN:
            o = stacks[CLOSE_TO_OPEN[ch]].pop()
            out.add(frozenset((col_to_key[o], col_to_key[col])))
    return out


# --------------------------------------------------------------------------- #
#  Row-label helpers                                                           #
# --------------------------------------------------------------------------- #

def _label(strand: tuple[str, str]) -> str:
    """'A[1_555]' -- compact strand label used inside per-block tags / metadata."""
    ch, sym = strand
    return f"{ch}[{sym}]"


def _span(residues: list, strand: tuple[str, str]) -> str:
    """First..last residue number of `strand`, e.g. '1-59'. Used in --metadata
    and in per-block tags for --block."""
    nums = [R[2] for R in residues if (R[0], R[1]) == strand]
    return f"{nums[0]}-{nums[-1]}" if nums else "absent"


def _strand_token(strand: tuple[str, str]) -> str:
    """ ('A', '1_555') -> 'A';   ('A', '2_755') -> 'A:2_755'.
    Identity operator is the default and stays implicit."""
    ch, sym = strand
    return ch if sym == IDENTITY else f"{ch}:{sym}"


def _layer_label(family: str, slot_num: int, show_layer: bool) -> str:
    """Row label for one layer. `family` is 'WC' for canonical cWW, or an LW
    code (cWW for non-canonical cWW, tWW, cWH, ...). When show_layer is True,
    prepend the slot number 'Ln '."""
    return f"L{slot_num} {family}" if show_layer else family


def _family_of(label: str) -> str:
    """Pull the family out of a layer label, regardless of --layer mode.
    'WC' -> 'WC'   'cWW' -> 'cWW'   'L0 WC' -> 'WC'   'L10 tWW' -> 'tWW'."""
    parts = label.split()
    return parts[1] if len(parts) >= 2 and parts[0].startswith("L") else parts[0]


def _format_unpaired(unpaired: list, strands: list[tuple[str, str]]) -> str:
    """ 'A8, A11, A29-32, A44, ...' -- compact comma-separated list of unpaired
    residues, with consecutive numbers within one strand collapsed to ranges.
    Symmetry mates use the same residue-with-bracketed-operator form as
    _res_id (e.g. 'A2[2_565]' for a single residue, 'A2-5[2_565]' for a range).
    Returns '(none)' when there are no unpaired residues."""
    if not unpaired:
        return "(none)"
    by_strand: dict[tuple, list[int]] = {}
    for ch, sym, num in unpaired:
        by_strand.setdefault((ch, sym), []).append(num)
    out: list[str] = []
    for s in strands:
        nums = sorted(by_strand.get(s, []))
        if not nums:
            continue
        ch, sym = s
        op_suffix = "" if sym == IDENTITY else f"[{sym}]"
        i = 0
        while i < len(nums):
            j = i
            while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
                j += 1
            out.append(f"{ch}{nums[i]}{op_suffix}" if i == j
                       else f"{ch}{nums[i]}-{nums[j]}{op_suffix}")
            i = j + 1
    return ", ".join(out)


def _res_id(R: tuple) -> str:
    """Compact residue label: chain+number, with the operator only when it is a
    symmetry mate (e.g. 'A24', or 'A1012[2_755]')."""
    ch, sym, num = R
    return f"{ch}{num}" if sym == IDENTITY else f"{ch}{num}[{sym}]"


def _compact_line(pairs: list) -> str:
    """A sparse layer as an explicit pair list, e.g. 'A24,A31 A25,A29'."""
    return " ".join(f"{_res_id(Ra)},{_res_id(Rb)}" for Ra, Rb, _ in sorted(pairs))


# --------------------------------------------------------------------------- #
#  Top-level notation builder                                                  #
# --------------------------------------------------------------------------- #

def build_notation(per_chain: dict[str, list], pairs: list, chains: list[str],
                   name: str = "RNA",
                   block: int | None = None, compact: bool = False,
                   noncanonical: bool = False, show_layer: bool = False,
                   show_metadata: bool = False,
                   show_unpaired: bool = False) -> tuple[str, bool]:
    """Build the layered notation and return (text, round-trip ok).

    Args:
        per_chain: {chain: [(residue_number, sequence_letter), ...]} -- caller
            obtains this via read_residues().
        pairs: [(Ra, Rb, lw)] where each R is (chain, symmetry, number) and lw
            is the Leontis-Westhof family. Caller obtains this from its own
            source (CIF annotation, FR3D TSV, ...).
        chains: chain ids in ruler order.
        name: label for the header.
        block: if set, wrap notation lines into blocks of this many columns.
        compact: only the canonical WC layer stays as full-width dot-bracket;
            every non-canonical layer is printed as an explicit pair list.
        noncanonical: drop true Watson-Crick pairs entirely (the WC layer
            becomes absent). cWW U-U / U-G wobbles still appear on the cWW row.
        show_layer: prepend slot numbers ('L0 WC', 'L1 cWW', 'L10 tWW', ...).
        show_unpaired: add a '# unpaired (N): A8, A11, A29-32, ...' comment
            line below the header listing residues that have no hydrogen-bond
            partner in the structure. Consecutive numbers within a strand are
            collapsed to ranges. The count is a structural property and stays
            the same regardless of --noncanonical: a stem residue that pairs
            only via WC is still paired in the structure, its pair is just
            hidden from the display.
        show_metadata: add a '# chains: ...' comment line below the header
            with per-strand chain/symmetry/range info.

    Row labels:
        WC                                       canonical Watson-Crick pairs
        cWW                                      non-canonical cWW (U-U, U-G wobbles)
        tWW cWH tWH cWS tWS cHW tHW cHH tHH
        cHS tHS cSW tSW cSH tSH cSS tSS          the other 17 directed LW families

    The header is NCBI/FASTA-style:
        >name|<strand1>|<strand2>|...
    where each strand is a chain id, with ':operator' appended only when the
    strand is a symmetry mate (not the identity 1_555)."""
    base_of = {(ch, num): letter
               for ch, lst in per_chain.items() for num, letter in lst}

    # Structural pairing snapshot, taken BEFORE the noncanonical display filter.
    # Used by show_unpaired so that residues paired only via WC are not falsely
    # reported as unpaired just because --noncanonical hid their pair from view.
    full_paired_keys = {R for Ra, Rb, _ in pairs for R in (Ra, Rb)}

    if noncanonical:
        pairs = [p for p in pairs
                 if not _is_canonical(p[2], base_of.get((p[0][0], p[0][2]), "?"),
                                            base_of.get((p[1][0], p[1][2]), "?"))]

    # One strand per (chain, symmetry): identity copy of each requested chain
    # plus any symmetry copies that appear in the pairs (deterministic order).
    sym_seen = sorted({R[:2] for Ra, Rb, _ in pairs for R in (Ra, Rb) if R[1] != IDENTITY})
    strands: list[tuple[str, str]] = []
    for c in chains:
        strands.append((c, IDENTITY))
        strands.extend(s for s in sym_seen if s[0] == c)

    columns: list[tuple] = []          # ("res", R, letter) | ("sep",)
    key_to_col: dict[tuple, int] = {}
    residues: list[tuple] = []
    for n, (ch, sym) in enumerate(strands):
        if n > 0:
            columns.append(("sep",))
        for num, letter in per_chain.get(ch, []):
            R = (ch, sym, num)
            key_to_col[R] = len(columns)
            columns.append(("res", R, letter))
            residues.append(R)
    width = len(columns)
    sep_cols = {i for i, col in enumerate(columns) if col[0] == "sep"}
    col_to_key = {i: col[1] for i, col in enumerate(columns) if col[0] == "res"}

    # Slot numbers for the 17 non-cWW LW families (L2..L18). cWW is split below.
    slot = {t: i + 1 for i, t in enumerate(DIRECTED)}

    cww_can: list = []                  # canonical cWW -> WC (slot 0)
    cww_non: list = []                  # non-canonical cWW -> cWW (slot 1)
    by_type: dict[str, list] = {}       # other 17 LW families -> slot 2..18
    for Ra, Rb, lw in pairs:
        if lw == "cWW":
            ba = base_of.get((Ra[0], Ra[2]), "?")
            bb = base_of.get((Rb[0], Rb[2]), "?")
            (cww_can if _is_canonical(lw, ba, bb) else cww_non).append((Ra, Rb, lw))
        else:
            by_type.setdefault(lw, []).append((Ra, Rb, lw))

    layers, overflow = [], []           # overflow: (type, bucket), packed per type

    def _place(label: str, group: list, lw_for_overflow: str) -> None:
        """Greedy single-layer placement; spill conflicts to overflow."""
        placed, used, extra = [], set(), []
        for Ra, Rb, lw in sorted(group):
            if Ra in used or Rb in used:
                extra.append((Ra, Rb, lw))
            else:
                placed.append((Ra, Rb, lw))
                used |= {Ra, Rb}
        layers.append((label, placed))
        overflow.extend((lw_for_overflow, b) for b in _pack(extra))

    if cww_can:
        _place(_layer_label("WC", 0, show_layer), cww_can, "cWW")
    if cww_non:
        _place(_layer_label("cWW", 1, show_layer), cww_non, "cWW")
    for t in DIRECTED:
        if t == "cWW" or t not in by_type:
            continue
        _place(_layer_label(t, slot[t], show_layer), by_type[t], t)
    for k, (t, bucket) in enumerate(overflow):
        layers.append((_layer_label(t, 19 + k, show_layer), bucket))

    seq_line = "".join(col[2] if col[0] == "res" else SEP for col in columns)
    rendered = [(label, render(p, key_to_col, width, sep_cols)) for label, p in layers]

    head = ">" + "|".join([name] + [_strand_token(s) for s in strands])
    if show_metadata:
        spans = ", ".join(f"{_label(s)}({_span(residues, s)})" for s in strands)
        head += f"\n# chains: {spans}; '{SEP}' separates chains"
    if show_unpaired:
        # Residues that have no hydrogen-bond partner in the FULL structure.
        # We use the pre-filter pair set so the count is a structural property
        # of the molecule, independent of --noncanonical (residues paired only
        # via Watson-Crick are still paired -- their pair is just hidden).
        unpaired = [R for R in residues if R not in full_paired_keys]
        head += f"\n# unpaired ({len(unpaired)}): {_format_unpaired(unpaired, strands)}"

    if compact:        # only the canonical WC layer stays as dot-bracket; every
        rows = [f"{'seq':12}: " + seq_line]   # non-canonical layer becomes a
        for (label, layer_pairs), (_, dbstr) in zip(layers, rendered):  # pair list
            line = dbstr if _family_of(label) == "WC" else _compact_line(layer_pairs)
            rows.append(f"{label:12}: " + line)
        text = head + "\n" + "\n".join(rows)
    elif block:
        body = "\n\n".join(_wrap(seq_line, rendered, columns, block))
        text = head + "\n\n" + body
    else:
        rows = [f"{'seq':12}: " + seq_line]
        rows += [f"{label:12}: " + string for label, string in rendered]
        text = head + "\n" + "\n".join(rows)

    ok = _roundtrip(rendered, col_to_key, pairs)
    return text, ok


def _roundtrip(rendered: list, col_to_key: dict, pairs: list) -> bool:
    """Parse the written lines back to pairs and require the exact input set.
    'WC' (canonical) maps back to LW family 'cWW' since both rows carry cWW
    geometry; whether a pair lands on the WC row or the cWW row is decided at
    re-render time from the bases."""
    recovered = set()
    for label, s in rendered:
        family = _family_of(label)
        lw = "cWW" if family == "WC" else family
        for keyset in parse(s, col_to_key):
            Ra, Rb = sorted(keyset)
            recovered.add((Ra, Rb, lw))
    return recovered == set(pairs)


def _wrap(seq_line: str, rendered: list, columns: list, block: int) -> list[str]:
    """Split the seq + layer lines into fixed-width blocks of `block` columns."""
    blocks = []
    for s in range(0, len(seq_line), block):
        e = min(s + block, len(seq_line))
        res = [col for col in columns[s:e] if col[0] == "res"]
        if res:
            r0, r1 = res[0][1], res[-1][1]
            tag = f"{_label(r0[:2])}:{r0[2]} .. {_label(r1[:2])}:{r1[2]}"
        else:
            tag = "separator"
        rows = [f"# cols {s}-{e - 1}  ({tag})",
                f"{'seq':12}: " + seq_line[s:e]]
        rows += [f"{label:12}: " + string[s:e] for label, string in rendered]
        blocks.append("\n".join(rows))
    return blocks
