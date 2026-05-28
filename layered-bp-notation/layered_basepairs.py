"""Layered base-pair notation from an FR3D base-pair TSV.

Chains are laid end to end on one ruler, separated by '&' (the ViennaRNA cofold
convention), and a pair between two chains is a bracket that opens in one chain
and closes in another. Each base is keyed by (chain, symmetry, number): the
symmetry field distinguishes a residue from a symmetry-generated copy of the
same chain, so self-complementary duplexes and crystal-packing contacts (where
FR3D reports a pair as chain A to chain A with a symmetry operator on one
partner) are represented as two strands rather than a spurious self-pairing.

The sequence and residue numbers come from pdbx_poly_seq_scheme when present,
otherwise from _atom_site (the coordinates), so structures lacking that
category (e.g. modelling results) still work.

This file owns the FR3D-specific reading (unit-id parsing, pair list); the rest
-- constants, sequence/numbering reading, the layered-notation builder, parser
helpers -- lives in common.py.

Usage:
    python3 layered_basepairs.py <cif> <tsv> <chains...> [--name NAME]
                                  [--block N] [--compact] [--noncanonical]
                                  [--layer] [--metadata]

Examples (mmCIF from files.rcsb.org/download/<ID>.cif, basepairs TSV from FR3D):
    python3 layered_basepairs.py 9CFN.cif 9cfn_fr3d_basepairs.tsv A --name 9CFN
    python3 layered_basepairs.py 1XPE.cif 1xpe_fr3d_basepairs.tsv A B --name 1XPE
    python3 layered_basepairs.py 2Q1R.cif 2q1r_fr3d_basepairs.tsv A --name 2Q1R

The notation prints to stdout; the round-trip result (True = lossless) to stderr.
"""

import argparse
import sys
from pathlib import Path

from common import (
    BLOCK, IDENTITY,
    _flip_lw, build_notation, read_residues,
)


# --------------------------------------------------------------------------- #
#  FR3D-specific reading                                                       #
# --------------------------------------------------------------------------- #

def _sym(unit: list[str]) -> str:
    """Symmetry operator from an FR3D unit id. A missing operator denotes the
    asymmetric-unit copy, which is the crystallographic identity, so it is
    reported explicitly as 1_555. FR3D unit id format is
    pdb|model|chain|comp|num|atom|alt|ins|symmetry."""
    return unit[8] if len(unit) > 8 and unit[8] else IDENTITY


def read_pairs(tsv: Path, chains: set[str]) -> list[tuple[tuple, tuple, str]]:
    """De-duplicated pairs (Ra, Rb, lw), where each residue R is
    (chain, symmetry, number). Both ends must be among `chains`; ordering is
    canonical so the LW direction is consistent."""
    seen: dict[frozenset, tuple] = {}
    for line in tsv.read_text().splitlines():
        c = line.strip().split("\t")
        if len(c) < 3:
            continue
        a, b = c[0].split("|"), c[2].split("|")
        if a[2] not in chains or b[2] not in chains:
            continue
        Ra, Rb, lw = (a[2], _sym(a), int(a[4])), (b[2], _sym(b), int(b[4])), c[1]
        if Ra > Rb:                       # canonical order; flip the label to match
            Ra, Rb, lw = Rb, Ra, _flip_lw(lw)
        seen[frozenset((Ra, Rb))] = (Ra, Rb, lw)
    return sorted(seen.values())


# --------------------------------------------------------------------------- #
#  Main                                                                        #
# --------------------------------------------------------------------------- #

if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Layered base-pair notation for an RNA complex "
                    "(symmetry-aware), driven by an FR3D base-pair TSV.")
    ap.add_argument("cif", type=Path, help="mmCIF file (sequence + residue numbers)")
    ap.add_argument("tsv", type=Path, help="FR3D basepairs TSV")
    ap.add_argument("chains", nargs="+",
                    help="chain ids, in the order to lay them on the ruler")
    ap.add_argument("--name", default="RNA", help="label for the header line")
    ap.add_argument("--block", type=int, nargs="?", const=BLOCK, default=None,
                    help=f"wrap lines into blocks of this many columns "
                         f"(default {BLOCK} when given with no value)")
    ap.add_argument("--compact", action="store_true",
                    help="keep only the canonical WC layer as dot-bracket; "
                         "print every non-canonical layer as an explicit pair "
                         "list (e.g. A24,A31); saves space on large RNA")
    ap.add_argument("--noncanonical", action="store_true",
                    help="show only non-canonical pairs (drop true Watson-Crick "
                         "A-U/G-C/A-T); a cWW U-U or U-G wobble is kept on the "
                         "cWW row")
    ap.add_argument("--layer", action="store_true",
                    help="prepend slot numbers to each layer label "
                         "(L0 WC, L1 cWW, L10 tWW, ...); off by default")
    ap.add_argument("--metadata", action="store_true",
                    help="add a '# chains: ...' comment line below the header "
                         "with per-strand chain/symmetry/range info; off by default")
    a = ap.parse_args()

    per_chain = read_residues(a.cif, a.chains)
    pairs = read_pairs(a.tsv, set(a.chains))
    text, ok = build_notation(per_chain, pairs, a.chains, name=a.name,
                              block=a.block, compact=a.compact,
                              noncanonical=a.noncanonical,
                              show_layer=a.layer, show_metadata=a.metadata)
    print(text)
    print(f"\n# round-trip recovers all pairs exactly: {ok}", file=sys.stderr)
    sys.exit(0 if ok else 1)
