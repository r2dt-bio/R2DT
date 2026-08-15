#!/usr/bin/env python3
"""Post-process already-generated compare viewers: rewrite each reference cif to
keep only the analysed RNA chain(s), so a deposited entry's extra chains don't
render as extra structures in 3D beside the model overlay.

This mirrors the fix now in ``r2dt.py`` (multichain.write_structure_chains) but
applies it to existing ``site/`` viewers without re-running the Docker pipeline.
The analysed chains are read from the reference panel's inlined ``chain_ids`` in
each ``index.html``.

Usage: python3 scripts/filter_reference_chains.py site
"""
import json
import re
import sys
import warnings
from pathlib import Path

from Bio.PDB import MMCIFIO, MMCIFParser, Select


def analysed_chains_and_ref(index_html: str):
    """Return (ref_structure_id, [chains]) from the first panel's inlined data."""
    txt = Path(index_html).read_text()
    # The first panel is the reference; grab its structureId and chain_ids.
    sid = re.search(r'"structureId"\s*:\s*"([^"]+)"', txt)
    ci = re.search(r'"chain_ids"\s*:\s*\[(.*?)\]', txt, re.S)
    if not sid or not ci:
        return None, []
    raw = "[" + ci.group(1) + "]"
    try:
        vals = json.loads(raw)
    except json.JSONDecodeError:
        return None, []
    chains = sorted({v for v in vals if v})
    return sid.group(1), chains


def filter_cif(cif_path: Path, keep):
    """Rewrite ``cif_path`` in place to keep only the chains in ``keep``."""
    keep = {str(c) for c in keep}

    class _Sel(Select):
        def accept_chain(self, chain):
            """Keep only chains in ``keep``."""
            return 1 if chain.id in keep else 0

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        structure = MMCIFParser(QUIET=True).get_structure(cif_path.stem, str(cif_path))
        present = {c.id for c in next(iter(structure))}
        if not present - keep:
            return "already-single"  # nothing extra to drop
        io = MMCIFIO()
        io.set_structure(structure)
        io.save(str(cif_path), _Sel())
    return f"kept {sorted(keep)} (dropped {sorted(present - keep)})"


def main(site_dir: str):
    """Filter every reference cif under ``site_dir`` to its analysed chains."""
    root = Path(site_dir)
    for index_html in root.rglob("viewer/index.html"):
        sid, chains = analysed_chains_and_ref(str(index_html))
        if not sid or not chains:
            continue
        cif = index_html.parent / f"{sid}.cif"
        if not cif.is_file():
            continue
        try:
            result = filter_cif(cif, chains)
        except Exception as exc:  # pylint: disable=broad-except
            result = f"ERROR: {exc}"
        print(f"{cif.relative_to(root)} -> {result}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "site")
