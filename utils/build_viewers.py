#!/usr/bin/env python3
"""Build a gallery ``index.html`` for a folder of R2DT 2D+3D viewers.

Scans every immediate subfolder of the site directory that looks like a
published viewer (has both ``index.html`` and ``2d.svg``), reads the PDB
id / chain from ``api.json``, and writes a gallery ``index.html`` linking
to each structure's viewer with its non-interactive 2D preview.

Usage:

    python3 utils/build_viewers.py [SITE_DIR]

``SITE_DIR`` defaults to ``output/site``.  Safe to re-run; it only
rewrites the top-level ``index.html``.
"""

import argparse
import html
import json
import sys
from pathlib import Path


def discover(site_dir: Path):
    """Yield (folder_name, pdb_id, chain_id) for each viewer subfolder."""
    entries = []
    for child in sorted(site_dir.iterdir()):
        if not child.is_dir():
            continue
        if not (child / "index.html").is_file() or not (child / "2d.svg").is_file():
            continue
        pdb_id, chain_id = child.name, ""
        api = child / "api.json"
        if api.is_file():
            try:
                meta = json.loads(api.read_text()).get("_meta", {})
                pdb_id = meta.get("structure_id") or pdb_id
                chain_id = meta.get("chain_id") or ""
            except (ValueError, OSError):
                pass
        entries.append((child.name, pdb_id, chain_id))
    return entries


def render(entries):
    """Return the gallery HTML for the given (folder, pdb_id, chain) list."""
    cards = []
    for folder, pdb_id, chain_id in entries:
        chain = f" · chain {html.escape(chain_id)}" if chain_id else ""
        f = html.escape(folder)
        pid = html.escape(pdb_id)
        img = f'<img src="{f}/2d.svg" alt="{pid} 2D diagram" loading="lazy">'
        cards.append(
            f'    <a class="card" href="{f}/index.html">\n'
            f'      <div class="thumb">{img}</div>\n'
            f'      <div class="label"><span class="pdb">{pid}</span>{chain}</div>\n'
            f"    </a>"
        )
    cards_html = (
        "\n".join(cards) if cards else "    <p>No structures published yet.</p>"
    )
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>R2DT 2D+3D viewers</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
         margin: 0; padding: 32px; color: #222; background: #fafafa; }}
  h1 {{ font-size: 22px; margin: 0 0 4px; }}
  p.sub {{ color: #666; margin: 0 0 28px; font-size: 14px; }}
  .grid {{ display: grid; gap: 20px;
          grid-template-columns: repeat(auto-fill, minmax(220px, 1fr)); }}
  .card {{ display: block; text-decoration: none; color: inherit;
          background: #fff; border: 1px solid #e4e4e4; border-radius: 10px;
          overflow: hidden; transition: box-shadow .15s, transform .15s; }}
  .card:hover {{ box-shadow: 0 6px 18px rgba(0,0,0,.10); transform: translateY(-2px); }}
  .thumb {{ height: 200px; display: flex; align-items: center; justify-content: center;
           background: #fff; padding: 12px; box-sizing: border-box; }}
  .thumb img {{ max-width: 100%; max-height: 100%; }}
  .label {{ padding: 10px 14px; border-top: 1px solid #f0f0f0; font-size: 14px; }}
  .pdb {{ font-weight: 600; letter-spacing: .02em; }}
</style>
</head>
<body>
<h1>R2DT 2D + 3D structure viewers</h1>
<p class="sub">{len(entries)} structure(s). Click a diagram to open its
   interactive 2D + 3D viewer.</p>
<div class="grid">
{cards_html}
</div>
</body>
</html>
"""


def main(argv=None):
    """Parse arguments and write the gallery index.html for the site dir."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "site_dir",
        nargs="?",
        default="output/site",
        help="Folder of published viewers (default: output/site)",
    )
    args = parser.parse_args(argv)

    site_dir = Path(args.site_dir)
    if not site_dir.is_dir():
        parser.error(f"site directory not found: {site_dir}")

    entries = discover(site_dir)
    (site_dir / "index.html").write_text(render(entries))
    print(f"Wrote {site_dir / 'index.html'} with {len(entries)} structure(s):")
    for folder, pdb_id, chain_id in entries:
        suffix = f" chain {chain_id}" if chain_id else ""
        print(f"  {pdb_id}{suffix}  ({folder}/)")
    if not entries:
        print("  (none found -- did you publish any viewers?)", file=sys.stderr)


if __name__ == "__main__":
    main()
