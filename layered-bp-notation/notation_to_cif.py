"""Inverse of standalone_lbn_script: layered notation -> mmCIF.

Reads a layered-notation text file (what the standalone prints) plus any mmCIF
of the same structure, and writes a new CIF that carries the three categories
the forward script reads:

  _ndb_base_pair_list        - who pairs with whom (auth numbering + operator id)
  _ndb_base_pair_annotation  - the Leontis-Westhof family per pair
  _pdbx_struct_oper_list     - maps each operator id to its 1_555-style name

The two annotation loops are regenerated straight from the notation; any
existing copies are replaced. The operator list is left untouched in the common
case, but it is rewritten (with old operators preserved by name and any new
ones appended) when the notation references operators the input CIF does not
already list -- e.g. when notation built from a symmetry-mate structure is
applied to a CIF that only carries the identity operator. Everything else in
the CIF is left untouched.

The CIF only has to give the script the chain layout and per-residue metadata
(asym_id, entity_id, label_seq_id, comp_id, etc.). It does NOT have to be
DNATCO-extended -- a vanilla PDB-downloaded CIF works, and so does a stripped
CIF that only carries _atom_site. The script tries _ndb_base_pair_list rows
first (to mirror DNATCO's exact spellings when present), then
_pdbx_poly_seq_scheme, and finally _atom_site.

Perfect round-trip is checked by re-reading the rebuilt CIF with the standard
reader: the set of pairs (chain, num, chain, num, LW) must equal the set the
notation encoded.

Usage:
    python3 notation_to_cif.py <notation.txt> <any.cif> -o <new.cif>

    # produce the notation, then invert it, then verify
    python3 standalone_lbn_script.py 9cfn.cif > 9cfn.lbn
    python3 notation_to_cif.py 9cfn.lbn 9cfn.cif -o 9cfn_rebuilt.cif
"""

import argparse
import re
import sys
from pathlib import Path

# Shared parsers/utilities live in common.py; the CIF-specific pair-list and
# operator readers (read_pairs, list_chains) live in standalone_lbn_script.py.
sys.path.insert(0, str(Path(__file__).parent))
from common import (                                                 # noqa: E402
    CLOSE_TO_OPEN, IDENTITY, NT, OPENERS, SEP,
    _flip_lw, _loop_columns, read_residues,
)
from standalone_lbn_script import (                                  # noqa: E402
    list_chains, read_pairs,
)

# LW -> family number (DNATCO undirected numbering, 1..12).
LW_NUM = {
    "cWW": 1, "tWW": 2, "cWH": 3, "tWH": 4, "cWS": 5, "tWS": 6,
    "cHW": 3, "tHW": 4, "cHH": 7, "tHH": 8, "cHS": 9, "tHS": 10,
    "cSW": 5, "tSW": 6, "cSH": 9, "tSH": 10, "cSS": 11, "tSS": 12,
}
EDGE = {"W": "Watson-Crick", "H": "Hoogsteen", "S": "Sugar"}


# --------------------------------------------------------------------------- #
#  Read the notation file                                                      #
# --------------------------------------------------------------------------- #

_FAMILIES = {"WC", "cWW", "cWH", "cWS", "cHW", "cHH", "cHS", "cSW", "cSH", "cSS",
             "tWW", "tWH", "tWS", "tHW", "tHH", "tHS", "tSW", "tSH", "tSS"}


def parse_header_strands(header: str) -> list[tuple[str, str]] | None:
    """ Header parser, accepts both new and old formats.
        New (NCBI-style):  '>9cfn|A|B'  or  '>2q1r|A|A:2_755'
        Old (verbose):     '>9CFN complex: chains A[1_555](2-58), B[1_555](2-58)'
        Returns [(chain, operator), ...], or None if no header is found."""
    s = header.lstrip(">").strip()

    # New format: pipe-separated and never includes the word 'chains'.
    if "|" in s and "chains" not in s:
        parts = [p.strip() for p in s.split("|")]
        out: list[tuple[str, str]] = []
        for p in parts[1:]:                           # parts[0] is the name
            if not p:
                continue
            ch, _, op = p.partition(":")
            out.append((ch.strip(), op.strip() or IDENTITY))
        return out or None

    # Old format: 'chains A[1_555](..), B[1_555](..)  (...)'
    m = re.search(r"chains\s+(.+)", s)
    if not m:
        return None
    raw = re.split(r"\s{2,}\(", m.group(1))[0]
    out = []
    for part in raw.split(","):
        part = part.strip()
        mm = re.match(r"([^\[\s]+)\[([^\]]+)\]", part)
        if mm:
            out.append((mm.group(1), mm.group(2)))
    return out or None


def _label_to_family(label: str) -> str | None:
    """Pull the family out of a layer label. Accepts both modes:
        'WC'        -> 'WC'
        'cWW'       -> 'cWW'
        'L0 WC'     -> 'WC'
        'L10 tWW'   -> 'tWW'
    Returns None if the label is not a recognised family."""
    parts = label.split()
    if not parts:
        return None
    if parts[0].startswith("L") and parts[0][1:].isdigit() and len(parts) >= 2:
        cand = parts[1]
    elif len(parts) == 1:
        cand = parts[0]
    else:
        return None
    return cand if cand in _FAMILIES else None


def parse_notation(path: Path):
    """Return (strands, layers) where layers = [(lw, kind, value), ...] preserving
    file order; kind is 'db' (dot-bracket string) or 'compact' (pair list string).
    'WC' rows are reported with lw='cWW' since canonical and non-canonical cWW
    share the cWW family -- which row a pair lands on is decided at re-render
    time from the bases."""
    lines = path.read_text().splitlines()
    strands = None
    layers: list[tuple[str, str, str]] = []
    for ln in lines:
        if ln.startswith(">") and strands is None:
            strands = parse_header_strands(ln)
            continue
        if ":" not in ln or ln.lstrip().startswith("#"):
            continue
        label, _, val = ln.partition(":")
        label = label.strip()
        val = val.strip()
        if label == "seq":
            continue
        family = _label_to_family(label)
        if family is None:
            continue
        lw = "cWW" if family == "WC" else family
        # Dot-bracket lines contain only .&([{<ABCDabcd)]}>; compact lines have
        # digits and commas. Easy disambiguator: presence of a comma.
        kind = "compact" if "," in val else "db"
        layers.append((lw, kind, val))
    return strands, layers


def build_columns(strands, per_chain):
    """Same column layout the forward script builds. Returns col_to_key."""
    col_to_key: dict[int, tuple] = {}
    col = 0
    for n, (ch, sym) in enumerate(strands):
        if n > 0:
            col += 1                                  # separator column
        for num, _letter in per_chain.get(ch, []):
            col_to_key[col] = (ch, sym, num)
            col += 1
    return col_to_key


def parse_db(s: str, col_to_key: dict[int, tuple]) -> list[tuple[tuple, tuple]]:
    pairs: list[tuple[tuple, tuple]] = []
    stacks: dict[str, list[int]] = {}
    for col, ch in enumerate(s):
        if ch in OPENERS:
            stacks.setdefault(ch, []).append(col)
        elif ch in CLOSE_TO_OPEN:
            o = stacks[CLOSE_TO_OPEN[ch]].pop()
            Ra, Rb = sorted((col_to_key[o], col_to_key[col]))
            pairs.append((Ra, Rb))
    return pairs


_RES_RE = re.compile(r"([A-Za-z]+)(\d+)(?:\[([^\]]+)\])?")


def parse_compact(s: str) -> list[tuple[tuple, tuple]]:
    """ 'A24,A31 A25[2_755],A29' -> list of (Ra, Rb)."""
    out: list[tuple[tuple, tuple]] = []
    for token in s.split():
        a, _, b = token.partition(",")
        def res(x):
            m = _RES_RE.fullmatch(x)
            return (m.group(1), m.group(3) or IDENTITY, int(m.group(2)))
        Ra, Rb = sorted((res(a), res(b)))
        out.append((Ra, Rb))
    return out


def recover_pairs(notation_path: Path, cif_path: Path):
    strands, layers = parse_notation(notation_path)
    if strands is None:                              # no header -> default chain order
        strands = [(c, IDENTITY) for c in list_chains(cif_path)]
    chains = list(dict.fromkeys(s[0] for s in strands))
    per_chain = read_residues(cif_path, chains)
    col_to_key = build_columns(strands, per_chain)

    pairs: dict[frozenset, tuple] = {}                # de-dup: (Ra, Rb) -> (Ra, Rb, lw)
    for lw, kind, val in layers:
        bricks = parse_db(val, col_to_key) if kind == "db" else parse_compact(val)
        for Ra, Rb in bricks:
            lw_dir = lw if Ra <= Rb else _flip_lw(lw)
            Ra2, Rb2 = sorted((Ra, Rb))
            pairs[frozenset((Ra2, Rb2))] = (Ra2, Rb2, lw_dir)
    return sorted(pairs.values())


# --------------------------------------------------------------------------- #
#  Pull per-residue metadata + operator name<->id from the original CIF        #
# --------------------------------------------------------------------------- #

def read_meta(cif_path: Path, chains: list[str]):
    """ {(auth_asym_id, auth_seq_id): {asym_id, entity_id, label_seq_id,
                                       comp_id, ins, alt}}.
    Sources, in order of preference:
      1. existing _ndb_base_pair_list rows (when the CIF is DNATCO-extended, so
         we mirror DNATCO's exact spellings),
      2. _pdbx_poly_seq_scheme (every PDB-deposited CIF carries this),
      3. _atom_site (the final fallback that any CIF is guaranteed to have)."""
    lines = cif_path.read_text().splitlines()
    meta: dict[tuple, dict] = {}

    idx, i = _loop_columns(lines, "_ndb_base_pair_list.")
    if idx:
        while i < len(lines):
            row = lines[i].strip()
            if row.startswith(("#", "loop_", "_")):
                break
            f = row.split()
            if len(f) >= len(idx):
                for s in ("1", "2"):
                    key = (f[idx[f"auth_asym_id_{s}"]],
                           int(f[idx[f"auth_seq_id_{s}"]]))
                    meta.setdefault(key, {
                        "asym_id":      f[idx[f"asym_id_{s}"]],
                        "entity_id":    f[idx[f"entity_id_{s}"]],
                        "label_seq_id": f[idx[f"seq_id_{s}"]],
                        "comp_id":      f[idx[f"comp_id_{s}"]],
                        "ins":          f[idx[f"PDB_ins_code_{s}"]],
                        "alt":          f[idx[f"alt_id_{s}"]],
                    })
            i += 1

    idx2, j = _loop_columns(lines, "_pdbx_poly_seq_scheme.")
    if idx2:
        c_ins = idx2.get("pdb_ins_code")
        while j < len(lines):
            row = lines[j].strip()
            if row.startswith(("#", "loop_", "_")):
                break
            f = row.split()
            if len(f) >= len(idx2):
                key = (f[idx2["pdb_strand_id"]], int(f[idx2["pdb_seq_num"]]))
                if key not in meta:
                    meta[key] = {
                        "asym_id":      f[idx2["asym_id"]],
                        "entity_id":    f[idx2["entity_id"]],
                        "label_seq_id": f[idx2["seq_id"]],
                        "comp_id":      f[idx2["mon_id"]],
                        "ins":          f[c_ins] if c_ins is not None else "?",
                        "alt":          ".",
                    }
            j += 1

    # Final fallback: pull what we can from _atom_site (the only category every
    # CIF is guaranteed to have). label_seq_id may be '.' for non-polymer rows.
    idx3, k = _loop_columns(lines, "_atom_site.")
    if idx3 and ("auth_asym_id" in idx3 and "auth_seq_id" in idx3
                 and "label_comp_id" in idx3):
        c_asym  = idx3.get("label_asym_id")
        c_ent   = idx3.get("label_entity_id")
        c_label = idx3.get("label_seq_id")
        c_ins   = idx3.get("pdbx_PDB_ins_code")
        c_alt   = idx3.get("label_alt_id")
        c_model = idx3.get("pdbx_PDB_model_num")
        while k < len(lines):
            row = lines[k].strip()
            if row.startswith(("#", "loop_", "_")):
                break
            f = row.split()
            if len(f) >= len(idx3):
                if c_model is not None and f[c_model] != "1":
                    k += 1
                    continue
                key = (f[idx3["auth_asym_id"]], int(f[idx3["auth_seq_id"]]))
                if key not in meta:
                    meta[key] = {
                        "asym_id":      f[c_asym]  if c_asym  is not None else "?",
                        "entity_id":    f[c_ent]   if c_ent   is not None else "1",
                        "label_seq_id": f[c_label] if c_label is not None else "?",
                        "comp_id":      f[idx3["label_comp_id"]],
                        "ins":          f[c_ins]   if c_ins   is not None else "?",
                        "alt":          f[c_alt]   if c_alt   is not None else ".",
                    }
            k += 1
    return meta


def _read_loop_rows(lines: list[str], i: int, n_cols: int) -> list[list[str]]:
    """Walk forward from `i` in a CIF, accumulating tokens until each row has
    `n_cols` tokens. mmCIF allows long rows to wrap across several text lines,
    so a naive line-by-line read mis-counts. Stops at the next loop_, category
    header, comment, or blank line. Quote-aware via shlex."""
    import shlex
    rows: list[list[str]] = []
    buf: list[str] = []
    while i < len(lines):
        s = lines[i].strip()
        if s == "" or s.startswith(("#", "loop_", "_", "data_")):
            if buf:                                    # incomplete trailing row
                rows.append(buf)
            break
        try:
            buf.extend(shlex.split(s))
        except ValueError:                             # unbalanced quote; treat as raw
            buf.extend(s.split())
        while len(buf) >= n_cols:
            rows.append(buf[:n_cols])
            buf = buf[n_cols:]
        i += 1
    return rows


def read_oper_name_to_id(cif_path: Path) -> dict[str, str]:
    """ {operator_name: id} -- e.g. '1_555' -> '1'. Handles both the
    single-row key-value form (one operator) and the loop_ form (multiple
    operators, possibly with rows wrapped across several text lines, as in
    2Q1R)."""
    lines = cif_path.read_text().splitlines()
    out: dict[str, str] = {}
    kv = [l for l in lines if l.strip().startswith("_pdbx_struct_oper_list.")
          and len(l.split()) >= 2]
    if kv:                                            # single-row key-value form
        d = {l.split()[0].split(".")[1]: l.split()[1] for l in kv}
        if "name" in d and "id" in d:
            out[d["name"]] = d["id"]
            return out
    idx, i = _loop_columns(lines, "_pdbx_struct_oper_list.")
    if "name" in idx and "id" in idx:
        for row in _read_loop_rows(lines, i, len(idx)):
            out[row[idx["name"]]] = row[idx["id"]]
    if IDENTITY not in out:
        out[IDENTITY] = "1"
    return out


# --------------------------------------------------------------------------- #
#  Emit the two CIF blocks                                                     #
# --------------------------------------------------------------------------- #

LIST_HEADER = """loop_
_ndb_base_pair_list.base_pair_id
_ndb_base_pair_list.PDB_model_number
_ndb_base_pair_list.asym_id_1
_ndb_base_pair_list.entity_id_1
_ndb_base_pair_list.seq_id_1
_ndb_base_pair_list.comp_id_1
_ndb_base_pair_list.PDB_ins_code_1
_ndb_base_pair_list.alt_id_1
_ndb_base_pair_list.struct_oper_id_1
_ndb_base_pair_list.asym_id_2
_ndb_base_pair_list.entity_id_2
_ndb_base_pair_list.seq_id_2
_ndb_base_pair_list.comp_id_2
_ndb_base_pair_list.PDB_ins_code_2
_ndb_base_pair_list.alt_id_2
_ndb_base_pair_list.struct_oper_id_2
_ndb_base_pair_list.auth_asym_id_1
_ndb_base_pair_list.auth_seq_id_1
_ndb_base_pair_list.auth_asym_id_2
_ndb_base_pair_list.auth_seq_id_2
"""

ANN_HEADER = """loop_
_ndb_base_pair_annotation.id
_ndb_base_pair_annotation.base_pair_id
_ndb_base_pair_annotation.orientation
_ndb_base_pair_annotation.base_1_edge
_ndb_base_pair_annotation.base_2_edge
_ndb_base_pair_annotation.l-w_family_num
_ndb_base_pair_annotation.l-w_family
_ndb_base_pair_annotation.class
_ndb_base_pair_annotation.subclass
"""


def emit_blocks(pairs, meta, oper_name_to_id) -> tuple[str, str]:
    list_rows, ann_rows = [], []
    for i, (Ra, Rb, lw) in enumerate(pairs, 1):
        ch1, s1, n1 = Ra
        ch2, s2, n2 = Rb
        m1 = meta.get((ch1, n1), {})
        m2 = meta.get((ch2, n2), {})
        op1 = oper_name_to_id.get(s1, "1")
        op2 = oper_name_to_id.get(s2, "1")
        c1, c2 = m1.get("comp_id", "?"), m2.get("comp_id", "?")
        list_rows.append(" ".join([
            str(i), "1",
            m1.get("asym_id", "?"), m1.get("entity_id", "1"),
            m1.get("label_seq_id", "?"), c1,
            m1.get("ins", "?"), m1.get("alt", "."), op1,
            m2.get("asym_id", "?"), m2.get("entity_id", "1"),
            m2.get("label_seq_id", "?"), c2,
            m2.get("ins", "?"), m2.get("alt", "."), op2,
            ch1, str(n1), ch2, str(n2),
        ]))
        orient = "cis" if lw and lw[0] == "c" else "trans"
        e1 = EDGE.get(lw[1] if len(lw) >= 3 else "?", "?")
        e2 = EDGE.get(lw[2] if len(lw) >= 3 else "?", "?")
        fam_num = LW_NUM.get(lw, 0)
        b1 = NT.get(c1, c1[-1] if c1 else "?").upper()
        b2 = NT.get(c2, c2[-1] if c2 else "?").upper()
        klass = f"{lw}_{b1}-{b2}"
        sub = f"{klass}_1"
        ann_rows.append(" ".join([
            str(i), str(i), orient, e1, e2, str(fam_num), lw, klass, sub,
        ]))
    return LIST_HEADER + "\n".join(list_rows) + "\n", ANN_HEADER + "\n".join(ann_rows) + "\n"


# --------------------------------------------------------------------------- #
#  Splice the new blocks into the CIF                                          #
# --------------------------------------------------------------------------- #

def _drop_loop_for(category: str, lines: list[str]) -> list[str]:
    """Remove the loop_ + headers + data rows of a single category from `lines`.
    A 'loop_' line followed by header lines all sharing the category prefix is
    treated as that category's block."""
    out: list[str] = []
    i = 0
    prefix = f"_{category}."
    while i < len(lines):
        if lines[i].strip() == "loop_":
            # peek for the first header line
            j = i + 1
            while j < len(lines) and lines[j].strip() == "":
                j += 1
            if j < len(lines) and lines[j].lstrip().startswith(prefix):
                # gobble loop_ + all category headers
                k = j
                while k < len(lines) and lines[k].lstrip().startswith(prefix):
                    k += 1
                # gobble data rows until next loop_/_section/#/data_/blank
                while k < len(lines):
                    s = lines[k].strip()
                    if s == "" or s.startswith(("loop_", "_", "#", "data_")):
                        break
                    k += 1
                i = k
                continue
        out.append(lines[i])
        i += 1
    return out


def _drop_kv_for(category: str, lines: list[str]) -> list[str]:
    """Remove single-row key-value form lines for a category (e.g.
    '_pdbx_struct_oper_list.name 1_555' on its own line)."""
    prefix = f"_{category}."
    return [l for l in lines if not l.lstrip().startswith(prefix)]


# --------------------------------------------------------------------------- #
#  Manage _pdbx_struct_oper_list: extend it when the notation needs operators  #
#  the input CIF does not already provide.                                     #
# --------------------------------------------------------------------------- #

OPER_HEADER = ("loop_\n"
               "_pdbx_struct_oper_list.id\n"
               "_pdbx_struct_oper_list.type\n"
               "_pdbx_struct_oper_list.name\n"
               "_pdbx_struct_oper_list.symmetry_operation\n")


def ensure_operators(pairs, oper_name_to_id):
    """Make sure every operator name the notation pairs reference has an id in
    `oper_name_to_id`. Returns (updated_map, rebuilt_block_or_None). When the
    block is not None, the splice step drops the original _pdbx_struct_oper_list
    (in whichever form it had) and writes the rebuilt one.

    The rebuilt block uses minimal columns (id/type/name/symmetry_operation);
    matrix data for already-known operators is dropped because there is no way
    to synthesise matrices for the newly-added operators. This trade-off only
    fires when extension is actually needed (rare)."""
    needed = {R[1] for p in pairs for R in (p[0], p[1])}
    needed.add(IDENTITY)
    missing = needed - set(oper_name_to_id.keys())
    if not missing:
        return oper_name_to_id, None
    used_ids = {int(v) for v in oper_name_to_id.values() if v.isdigit()}
    next_id = max(used_ids, default=0) + 1
    updated = dict(oper_name_to_id)
    for name in sorted(missing):
        updated[name] = str(next_id)
        next_id += 1
    rows = []
    for name, oid in updated.items():
        if name == IDENTITY:
            rows.append((oid, "'identity operation'", name, "x,y,z"))
        else:
            rows.append((oid, "'crystal symmetry operation'", name, "?"))
    rows.sort(key=lambda r: int(r[0]))
    block = OPER_HEADER + "\n".join(" ".join(r) for r in rows) + "\n"
    return updated, block


def splice(cif_text: str, list_block: str, ann_block: str,
           oper_block: str | None) -> str:
    """Drop any existing copies of the loops being rebuilt, then append the new
    ones at the end (just before any trailing data_ or EOF). CIFs are unordered,
    so appending is semantically equivalent and avoids fragile in-place edits.
    When the CIF has no copy of a loop to start with (vanilla / minimal CIFs),
    the drop step is a no-op and we simply add the loop.

    `oper_block` is None in the common case (input already lists every operator
    the notation needs) -- the original `_pdbx_struct_oper_list` is then left
    completely untouched, preserving its matrix data."""
    lines = cif_text.splitlines(keepends=True)
    lines = _drop_loop_for("ndb_base_pair_list", lines)
    lines = _drop_loop_for("ndb_base_pair_annotation", lines)
    if oper_block is not None:
        lines = _drop_loop_for("pdbx_struct_oper_list", lines)
        lines = _drop_kv_for("pdbx_struct_oper_list", lines)
    body = "".join(lines).rstrip() + "\n"
    parts = [body, list_block, ann_block]
    if oper_block is not None:
        parts.append(oper_block)
    return "\n".join(parts).rstrip() + "\n"


# --------------------------------------------------------------------------- #
#  Main                                                                        #
# --------------------------------------------------------------------------- #

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Rebuild a CIF's base-pair annotation from layered notation.")
    ap.add_argument("notation", type=Path, help="layered notation text file")
    ap.add_argument("cif", type=Path,
                    help="any mmCIF of the same structure (DNATCO-extended, "
                         "vanilla PDB download, or even atom-site-only)")
    ap.add_argument("-o", "--out", type=Path, required=True, help="output CIF path")
    a = ap.parse_args()

    pairs = recover_pairs(a.notation, a.cif)
    print(f"# recovered {len(pairs)} pairs from notation")

    strands, _ = parse_notation(a.notation)
    if strands is None:
        strands = [(c, IDENTITY) for c in list_chains(a.cif)]
    chains = list(dict.fromkeys(s[0] for s in strands))

    meta = read_meta(a.cif, chains)
    oper_name_to_id = read_oper_name_to_id(a.cif)

    # If the notation references operators the input CIF does not list, extend
    # _pdbx_struct_oper_list so the rebuilt CIF is self-consistent.
    prev_opers = set(oper_name_to_id)
    oper_name_to_id, oper_block = ensure_operators(pairs, oper_name_to_id)
    if oper_block is not None:
        added = sorted(set(oper_name_to_id) - prev_opers)
        print(f"# extended _pdbx_struct_oper_list with operators: {added} "
              f"(matrix data not synthesised; pair round-trip is unaffected)")

    list_block, ann_block = emit_blocks(pairs, meta, oper_name_to_id)
    new_text = splice(a.cif.read_text(), list_block, ann_block, oper_block)
    a.out.write_text(new_text)
    print(f"# wrote {a.out}")

    # Verify: re-read the rebuilt CIF and require its pair set to equal exactly
    # what the notation encoded (this is the perfect round-trip; the original CIF
    # may have more pairs if the notation was --noncanonical).
    chains_set = set(chains)
    from_notation = set(pairs)
    from_rebuilt = set(read_pairs(a.out, chains_set))
    ok = from_notation == from_rebuilt
    print(f"# notation -> CIF -> pairs round-trip: {ok}  ({len(from_rebuilt)} pairs)")
    if not ok:
        miss = from_notation - from_rebuilt
        extra = from_rebuilt - from_notation
        if miss:
            print(f"#   missing from rebuilt ({len(miss)}):")
            for p in sorted(miss)[:10]:
                print(f"#     {p}")
        if extra:
            print(f"#   extra in rebuilt ({len(extra)}):")
            for p in sorted(extra)[:10]:
                print(f"#     {p}")
    # Informational only: when the original CIF already carries an annotation
    # block, report how the rebuilt set relates to it. Vanilla / minimal CIFs
    # have no such block, so this section quietly contributes nothing.
    try:
        orig_full = set(read_pairs(a.cif, chains_set))
    except (KeyError, ValueError):
        orig_full = None
    if orig_full is not None:
        if from_rebuilt == orig_full:
            print("# rebuilt CIF == original CIF (semantically; all pairs preserved).")
        else:
            dropped = orig_full - from_rebuilt
            if dropped:
                print(f"# note: {len(dropped)} pairs were absent from the notation "
                      f"(e.g. canonical WC pairs dropped by --noncanonical).")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
