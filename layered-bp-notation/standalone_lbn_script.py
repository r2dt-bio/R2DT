"""Layered base-pair notation + 3D view, from a single DNATCO-extended mmCIF.

One command does everything:
  1. reads the base pairs straight from the CIF's own annotation,
  2. prints the layered dot-bracket notation (lossless; round-trip checked),
  3. prints iCn3D commands you paste into iCn3D to draw every pair on the real
     3D structure, each pair its own colour (canonical Watson-Crick pairs grey).

You load the structure in iCn3D, paste the commands into its command box, press
Enter. No fragile long URLs. Use --script FILE to load them via File > Open
instead.

This file owns the CIF-specific reading (operator list, annotation, pair list,
3D coordinates) and the iCn3D overlay. The rest -- constants, sequence/numbering
reading, the layered-notation builder, parser helpers -- lives in common.py.

Usage:
    python3 standalone_lbn_script.py <cif> [chains...] [--name NAME] [--id PDBID]
                                     [--noncanonical] [--compact] [--block N]
                                     [--layer] [--metadata] [--script FILE]

    python3 standalone_lbn_script.py 9cfn.cif --name 9cfn          # auto chains
    python3 standalone_lbn_script.py 1XPE.cif A B --name 1XPE      # pick chains
    python3 standalone_lbn_script.py 9cfn.cif --noncanonical       # only non-canonical
    python3 standalone_lbn_script.py 9cfn.cif --script 9cfn.icn3d  # write a script file

The notation prints to stdout; the round-trip result and the iCn3D commands (or
a --script file path) follow.
"""

import argparse
import sys
from pathlib import Path

from common import (
    BLOCK, DIRECTED, IDENTITY,
    _flip_lw, _is_canonical, _loop_columns, _read_loop_rows,
    _res_id, build_notation, read_residues,
)


# --------------------------------------------------------------------------- #
#  iCn3D constants                                                             #
# --------------------------------------------------------------------------- #

ICN3D_URL = "https://www.ncbi.nlm.nih.gov/Structure/icn3d/full.html"

# Bold, well-separated colours so non-canonical pairs read from a distance.
# Canonical Watson-Crick is muted grey; every non-canonical pair gets a vivid hue.
CWW_COLOR = "888888"
BOLD = ["FF00FF", "FF8000", "00C800", "00E5E5", "FFD000", "FF0000",
        "0000FF", "8000FF", "FF0080", "00FF80", "A0FF00", "0080FF",
        "FF4060", "8B4513", "FF80C0", "00B0FF", "9400D3"]


# --------------------------------------------------------------------------- #
#  Read the base pairs from the CIF's own annotation                           #
# --------------------------------------------------------------------------- #

def _read_oper_map(lines: list[str]) -> dict[str, str]:
    """Map a struct_oper id to its symmetry name (e.g. '1' -> '1_555').
    Handles both the single-row key-value form and the loop_ form; loop_ rows
    may wrap across several text lines (e.g. 2Q1R), so tokens are accumulated
    until each row has the expected column count."""
    kv = [l for l in lines if l.strip().startswith("_pdbx_struct_oper_list.")]
    if kv and all(len(l.split()) >= 2 for l in kv):       # single-row key-value form
        d = {l.split()[0].split(".")[1]: l.split()[1] for l in kv}
        return {d.get("id", "1"): d.get("name", IDENTITY)}
    idx, i = _loop_columns(lines, "_pdbx_struct_oper_list.")
    m: dict[str, str] = {}
    if not idx:
        return m
    for row in _read_loop_rows(lines, i, len(idx)):
        m[row[idx["id"]]] = row[idx["name"]]
    return m


def _read_annotation(lines: list[str]) -> dict[str, str]:
    """Map base_pair_id -> Leontis-Westhof family string (cWW, tWH, ...) from
    _ndb_base_pair_annotation."""
    idx, i = _loop_columns(lines, "_ndb_base_pair_annotation.")
    c_id, c_lw = idx["base_pair_id"], idx["l-w_family"]
    lw_by_id: dict[str, str] = {}
    while i < len(lines):
        row = lines[i].strip()
        if row.startswith(("#", "loop_", "_")):
            break
        f = row.split()
        if len(f) >= len(idx):
            lw_by_id[f[c_id]] = f[c_lw]
        i += 1
    return lw_by_id


def read_pairs(cif: Path, chains: set[str]) -> list[tuple[tuple, tuple, str]]:
    """De-duplicated pairs (Ra, Rb, lw) read from the CIF annotation, where each
    residue R is (chain, symmetry, number). _ndb_base_pair_list gives the
    residues and operators; _ndb_base_pair_annotation gives the LW family,
    joined on base_pair_id. Ordering is canonical so the LW direction stays
    consistent."""
    lines = cif.read_text().splitlines()
    oper = _read_oper_map(lines)
    lw_by_id = _read_annotation(lines)
    idx, i = _loop_columns(lines, "_ndb_base_pair_list.")
    seen: dict[frozenset, tuple] = {}
    while i < len(lines):
        row = lines[i].strip()
        if row.startswith(("#", "loop_", "_")):
            break
        f = row.split()
        if len(f) >= len(idx):
            ch1, ch2 = f[idx["auth_asym_id_1"]], f[idx["auth_asym_id_2"]]
            if ch1 in chains and ch2 in chains:
                s1 = oper.get(f[idx["struct_oper_id_1"]], IDENTITY)
                s2 = oper.get(f[idx["struct_oper_id_2"]], IDENTITY)
                Ra = (ch1, s1, int(f[idx["auth_seq_id_1"]]))
                Rb = (ch2, s2, int(f[idx["auth_seq_id_2"]]))
                lw = lw_by_id.get(f[idx["base_pair_id"]])
                if Ra > Rb:                # canonical order; flip the label to match
                    Ra, Rb, lw = Rb, Ra, _flip_lw(lw)
                seen[frozenset((Ra, Rb))] = (Ra, Rb, lw)
        i += 1
    return sorted(seen.values())


def list_chains(cif: Path) -> list[str]:
    """Chains that appear in the base-pair list, sorted. Used when no chains
    are given; these are exactly the nucleic-acid chains that form annotated
    pairs."""
    lines = cif.read_text().splitlines()
    idx, i = _loop_columns(lines, "_ndb_base_pair_list.")
    chains: set[str] = set()
    while i < len(lines):
        row = lines[i].strip()
        if row.startswith(("#", "loop_", "_")):
            break
        f = row.split()
        if len(f) >= len(idx):
            chains.add(f[idx["auth_asym_id_1"]])
            chains.add(f[idx["auth_asym_id_2"]])
        i += 1
    return sorted(chains)


# --------------------------------------------------------------------------- #
#  Build the iCn3D 3D overlay (same pairs, plus coordinates from the CIF)      #
# --------------------------------------------------------------------------- #

def read_asu_c1(cif: Path, chains: set[str]) -> dict[tuple, tuple[float, float, float]]:
    """{(chain, num): (x, y, z)} for the asymmetric unit, using each residue's
    C1' atom (falls back to the residue's first atom if C1' is missing)."""
    lines = cif.read_text().splitlines()
    idx, i = _loop_columns(lines, "_atom_site.")
    c_ch, c_num, c_atom = idx["auth_asym_id"], idx["auth_seq_id"], idx["label_atom_id"]
    c_x, c_y, c_z = idx["Cartn_x"], idx["Cartn_y"], idx["Cartn_z"]
    c_model = idx.get("pdbx_PDB_model_num")
    coords: dict[tuple, tuple] = {}
    fallback: dict[tuple, tuple] = {}
    while i < len(lines):
        row = lines[i].strip()
        if row.startswith(("#", "loop_", "_")):
            break
        f = row.split()
        if len(f) >= len(idx):
            if c_model is not None and f[c_model] != "1":
                i += 1
                continue
            ch = f[c_ch]
            if ch in chains:
                key = (ch, int(f[c_num]))
                xyz = (float(f[c_x]), float(f[c_y]), float(f[c_z]))
                if f[c_atom].strip('"') == "C1'":
                    coords[key] = xyz
                else:
                    fallback.setdefault(key, xyz)
        i += 1
    for k, v in fallback.items():
        coords.setdefault(k, v)
    return coords


def read_oper_matrices(cif: Path) -> dict[str, tuple[tuple, tuple]]:
    """{operator_name: (3x3 row-major R, 3-vector t)} from _pdbx_struct_oper_list.
    Matrices are already Cartesian, so a symmetry mate of point p is R*p + t.
    Loop rows may wrap across text lines (e.g. 2Q1R)."""
    lines = cif.read_text().splitlines()
    out: dict[str, tuple] = {IDENTITY: ((1., 0., 0., 0., 1., 0., 0., 0., 1.),
                                        (0., 0., 0.))}
    if not any(l.strip().startswith("_pdbx_struct_oper_list.matrix[1][1]") for l in lines):
        return out
    idx, i = _loop_columns(lines, "_pdbx_struct_oper_list.")
    rk = ("matrix[1][1]", "matrix[1][2]", "matrix[1][3]",
          "matrix[2][1]", "matrix[2][2]", "matrix[2][3]",
          "matrix[3][1]", "matrix[3][2]", "matrix[3][3]")
    for row in _read_loop_rows(lines, i, len(idx)):
        out[row[idx["name"]]] = (tuple(float(row[idx[k]]) for k in rk),
                                 tuple(float(row[idx[k]]) for k in
                                       ("vector[1]", "vector[2]", "vector[3]")))
    return out


def _apply(R: tuple, asu: dict, opers: dict) -> tuple | None:
    """Cartesian C1' coordinate of residue R=(chain, symmetry, num), applying
    the symmetry operator to the asymmetric-unit position. None if unavailable."""
    chain, sym, num = R
    base = asu.get((chain, num))
    if base is None:
        return None
    if sym == IDENTITY:
        return base
    M = opers.get(sym)
    if M is None:
        return None
    r, t = M
    x, y, z = base
    return (r[0] * x + r[1] * y + r[2] * z + t[0],
            r[3] * x + r[4] * y + r[5] * z + t[1],
            r[6] * x + r[7] * y + r[8] * z + t[2])


def build_icn3d_lines(cif: Path, chains: list[str], noncanonical: bool = False
                      ) -> tuple[list[str], bool, int, list[tuple[str, str]]]:
    """iCn3D 'add line' commands for the same pairs the notation uses, drawn
    between C1' atoms. Every canonical Watson-Crick pair is the same thin grey;
    each non-canonical pair (including a cWW U-U) gets its own bold thick colour.
    Returns (commands, used a symmetry mate?, skipped count, legend list)."""
    pairs = read_pairs(cif, set(chains))
    asu = read_asu_c1(cif, set(chains))
    opers = read_oper_matrices(cif)
    base_of = {(ch, num): letter
               for ch, lst in read_residues(cif, list(chains)).items()
               for num, letter in lst}

    cmds: list[str] = []
    used_sym = False
    skipped = 0
    legend: list[tuple[str, str]] = []   # (colour, "pair family") per non-canonical pair
    nc = 0                               # running index -> a distinct colour per pair
    for Ra, Rb, lw in pairs:
        ba = base_of.get((Ra[0], Ra[2]), "?")
        bb = base_of.get((Rb[0], Rb[2]), "?")
        canonical = _is_canonical(lw, ba, bb)
        if noncanonical and canonical:
            continue
        pa, pb = _apply(Ra, asu, opers), _apply(Rb, asu, opers)
        if pa is None or pb is None:
            skipped += 1
            continue
        if Ra[1] != IDENTITY or Rb[1] != IDENTITY:
            used_sym = True
        if canonical:
            color, radius = CWW_COLOR, "0.3"             # canonical WC pair: grey, thin
        else:
            color, radius = BOLD[nc % len(BOLD)], "0.8"  # non-canonical pair: own colour
            nc += 1
            legend.append((color, f"{lw} {ba}-{bb} {_res_id(Ra)}-{_res_id(Rb)}"))
        cmds.append(f"add line | x {pa[0]:.3f} y {pa[1]:.3f} z {pa[2]:.3f} "
                    f"| x {pb[0]:.3f} y {pb[1]:.3f} z {pb[2]:.3f} "
                    f"| color {color} | dashed false | type {lw} | radius {radius}")
    return cmds, used_sym, skipped, legend


# --------------------------------------------------------------------------- #
#  Main                                                                        #
# --------------------------------------------------------------------------- #

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Layered base-pair notation + iCn3D commands to draw the "
                    "pairs in 3D, from a single annotated mmCIF.")
    ap.add_argument("cif", type=Path,
                    help="mmCIF with _ndb_base_pair_list/_annotation categories")
    ap.add_argument("chains", nargs="*",
                    help="chain ids in ruler order; default: all chains that pair")
    ap.add_argument("--name", default="RNA", help="label for the notation header")
    ap.add_argument("--id", default=None,
                    help="PDB id iCn3D should load (default: the cif file name); "
                         "must be an id iCn3D can load")
    ap.add_argument("--noncanonical", action="store_true",
                    help="draw only non-canonical pairs (cWW U-U/G-U included)")
    ap.add_argument("--block", type=int, nargs="?", const=BLOCK, default=None,
                    help=f"wrap notation lines into blocks of N columns "
                         f"(default {BLOCK})")
    ap.add_argument("--compact", action="store_true",
                    help="keep only the canonical WC layer as dot-bracket; "
                         "print every non-canonical layer as an explicit pair "
                         "list (e.g. A24,A31); saves space on large RNA")
    ap.add_argument("--layer", action="store_true",
                    help="prepend slot numbers to each layer label "
                         "(L0 WC, L1 cWW, L10 tWW, ...). Off by default: "
                         "labels are plain family names (WC, cWW, tWW, ...)")
    ap.add_argument("--metadata", action="store_true",
                    help="add a '# chains: ...' comment line below the header "
                         "with per-strand chain/symmetry/range info "
                         "(e.g. 'A[1_555](1-59)'); off by default")
    ap.add_argument("--script", metavar="FILE",
                    help="write the iCn3D commands to FILE (one per line) to "
                         "load via File > Open File > State/Script File, "
                         "instead of a URL")
    a = ap.parse_args()

    chains = a.chains or list_chains(a.cif)

    # 1) notation -- caller does the CIF reads, common.py does the layering
    per_chain = read_residues(a.cif, chains)
    pairs = read_pairs(a.cif, set(chains))
    text, ok = build_notation(per_chain, pairs, chains, name=a.name, block=a.block,
                              compact=a.compact, noncanonical=a.noncanonical,
                              show_layer=a.layer, show_metadata=a.metadata)
    print(text)
    print(f"\n# round-trip recovers all pairs exactly: {ok}", file=sys.stderr)

    # 2) 3D view
    cmds, used_sym, skipped, legend = build_icn3d_lines(a.cif, chains,
                                                        noncanonical=a.noncanonical)
    pdb_id = a.id or a.cif.stem.lower()   # pdbids load lowercase in iCn3D
    # 'set assembly on' must come first so the symmetry mate exists before the lines.
    draw_cmds = (["set assembly on"] if used_sym else []) + cmds

    if not cmds:
        msg = "non-canonical " if a.noncanonical else ""
        hint = " (drop --noncanonical to see the canonical pairs)" if a.noncanonical else ""
        print(f"\n# Nothing to draw: this structure has no {msg}pairs{hint}.",
              file=sys.stderr)
        sys.exit(0 if ok else 1)

    print(f"\n# iCn3D 3D view: {len(cmds)} pairs"
          f"{', non-canonical only' if a.noncanonical else ''}"
          f"{'; symmetry assembly included' if used_sym else ''}.")
    if legend:
        print("# non-canonical pairs (each drawn in its own colour):")
        for color, label in legend:
            print(f"#   {label}  #{color}")

    if a.script:
        # self-contained script file: load the structure, then draw the pairs.
        Path(a.script).write_text("\n".join([f"load pdb {pdb_id}"] + draw_cmds) + "\n")
        print(f"# Script written to {a.script} -> in iCn3D: File > Open File > State/Script File")
    else:
        print("#")
        print(f"# STEP 1 - open the structure in iCn3D (any browser):")
        print(f"#   {ICN3D_URL}?pdbid={pdb_id}")
        print("# STEP 2 - paste the lines below into iCn3D's command box (the '>' log at the")
        print("#   bottom), click right after the '>', and press Enter:")
        print("\n".join(draw_cmds))

    if skipped:
        print(f"# {skipped} pairs skipped (no coordinates)", file=sys.stderr)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
