#!/usr/bin/env python3
"""Refresh LBN data and vendored viewer assets in published viewer folders.

Use when ``output/site/*/`` was built before the 1D (LBN) panel or 2D
selection styling landed, without re-running the full ``pdb_2d_3d`` pipeline.
"""

import argparse
import json
import re
import shutil
import sys
from pathlib import Path
from typing import Optional, Tuple

from utils import lbn_export, viewer_html

_REPO = Path(__file__).resolve().parent.parent
_VIEWER_SRC = _REPO / "data" / "viewer"
_CONFIG_RE = re.compile(r"window\.R2DT_CONFIG\s*=\s*(\{.*?\});", re.DOTALL)


def _manifest_field(folder: Path, key: str) -> Optional[str]:
    manifest_path = folder / "manifest.json"
    if not manifest_path.is_file():
        return None
    try:
        data = json.loads(manifest_path.read_text())
        val = data.get(key)
        return val if val is not None else None
    except (ValueError, OSError):
        return None


def _parse_config(index_html: str) -> Optional[dict]:
    match = _CONFIG_RE.search(index_html)
    if not match:
        return None
    try:
        return json.loads(match.group(1))
    except ValueError:
        return None


def _structure_from_folder(folder: Path) -> Tuple[str, str]:
    for pattern, fmt in (("*.cif", "cif"), ("*.pdb", "pdb")):
        hits = sorted(folder.glob(pattern))
        if hits:
            return hits[0].name, fmt
    raise FileNotFoundError(f"no .cif/.pdb structure file in {folder}")


def _annotation_source(index_html: str) -> Optional[str]:
    if "mmCIF's own" in index_html or "DNATCO/NDB" in index_html:
        return viewer_html.ANNOTATION_SOURCE_HTML["cif"]
    return None


def upgrade_viewer_folder(folder: Path) -> bool:
    """Write ``lbn.json``, copy viewer assets, and re-render ``index.html``."""
    api_path = folder / "api.json"
    fr3d_path = folder / "fr3d.json"
    if not api_path.is_file() or not fr3d_path.is_file():
        return False

    api_data = json.loads(api_path.read_text())
    fr3d_data = json.loads(fr3d_path.read_text())
    meta = api_data.get("_meta", {})

    index_path = folder / "index.html"
    index_html = index_path.read_text() if index_path.is_file() else ""
    config = _parse_config(index_html)

    structure_id = (
        (config or {}).get("structureId")
        or _manifest_field(folder, "structureId")
        or meta.get("structure_id")
        or folder.name.upper()
    )
    chain_id = (
        (config or {}).get("chainId")
        or _manifest_field(folder, "chainId")
        or meta.get("chain_id")
        or ""
    )

    if config and config.get("structureUrl"):
        structure_filename = config["structureUrl"].lstrip("./")
        structure_format = config.get("structureFormat") or "pdb"
    else:
        structure_filename, structure_format = _structure_from_folder(folder)

    lbn_data = lbn_export.build_lbn_data(api_data, fr3d_data)
    (folder / "lbn.json").write_text(json.dumps(lbn_data))
    for name in (
        viewer_html.VIEWER_JS_FILENAME,
        viewer_html.R2DT_CSS_FILENAME,
    ):
        shutil.copyfile(_VIEWER_SRC / name, folder / name)

    viewer_html.render(
        folder,
        structure_id=structure_id,
        chain_id=chain_id or None,
        structure_filename=structure_filename,
        structure_format=structure_format,
        annotation_source=_annotation_source(index_html),
    )
    return True


def upgrade_site(site_dir: Path) -> int:
    """Upgrade every viewer subfolder under ``site_dir``; return count upgraded."""
    count = 0
    for child in sorted(site_dir.iterdir()):
        if not child.is_dir():
            continue
        if upgrade_viewer_folder(child):
            count += 1
            print(f"  upgraded {child.name}/")
    return count


def main(argv=None):
    """CLI entry point for ``python -m utils.upgrade_viewers``."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "site_dirs",
        nargs="*",
        default=["output/site"],
        help="Gallery root(s) containing viewer subfolders (default: output/site)",
    )
    args = parser.parse_args(argv)

    total = 0
    for raw in args.site_dirs:
        site_dir = Path(raw)
        if not site_dir.is_dir():
            print(f"skip: {site_dir} (not a directory)", file=sys.stderr)
            continue
        print(f"Upgrading viewers under {site_dir}/")
        total += upgrade_site(site_dir)

    print(f"Done — {total} viewer folder(s) refreshed.")
    return 0 if total else 1


if __name__ == "__main__":
    sys.exit(main())
