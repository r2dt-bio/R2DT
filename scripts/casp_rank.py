#!/usr/bin/env python3
"""Download a CASP season's official RNA results table and emit
``top<N>_by_target.json`` — the per-target, TM-score-ranked model list that
``casp_fetch.py`` reads to pick which models to fetch/compare.

Two on-the-wire formats exist (both official predictioncenter.org tables, no
local scoring):
- ``casp_table`` (CASP16): whitespace-column blocks, one per target, starting
  with a ``Target: <id>`` line; rows already ordered by TM_score.
- ``casp15_csv`` (CASP15): one flat CSV for every target/group/model, columns
  ``target,gr_code,model,...,tm_score,...``; the on-disk model filename isn't
  given directly and is reconstructed as ``<target>TS<gr_code:03d>_<model>``
  (verified against the actual CASP15 predictions tarballs).

Output shape (unchanged from the original hand-derived CASP16 file this
replaces)::

    {"<target>": [{"rank": 1, "model": "R1107TS029_1", "tm": "0.368"}, ...]}
"""
import argparse
import csv
import io
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
# pylint: disable=wrong-import-position
from casp_config import download, get_season  # noqa: E402


def parse_casp_table(text):
    """Parse CASP16-style ``Target: <id>`` whitespace-column blocks."""
    out = {}
    blocks = re.split(r"^Target:\s*(\S+)\s*$", text, flags=re.MULTILINE)[1:]
    for target, body in zip(blocks[0::2], blocks[1::2]):
        rows = []
        for line in body.splitlines():
            fields = line.split()
            if len(fields) < 7 or not fields[0].isdigit():
                continue
            model, tm = fields[1], fields[6]
            rows.append((model, None if tm == "-" else float(tm)))
        out[target] = rows
    return out


def parse_casp15_csv(text):
    """Parse the CASP15 flat scores CSV (one row per target/group/model)."""
    out = {}
    for row in csv.DictReader(io.StringIO(text)):
        target = row["target"]
        model = f"{target}TS{int(row['gr_code']):03d}_{row['model']}"
        tm = row.get("tm_score") or ""
        out.setdefault(target, []).append(
            (model, None if tm in ("", "NA") else float(tm))
        )
    return out


PARSERS = {"casp_table": parse_casp_table, "casp15_csv": parse_casp15_csv}


def top_n(rows, n):
    """Sort ``[(model, tm)]`` by TM-score desc (missing scores last), rank 1..n."""
    ranked = sorted(rows, key=lambda r: (r[1] is None, -(r[1] or 0.0)))[:n]
    return [
        {"rank": i, "model": model, "tm": f"{tm:.3f}" if tm is not None else "-"}
        for i, (model, tm) in enumerate(ranked, 1)
    ]


def main():
    """CLI entry point: rank one season's official results table, write top-N JSON."""
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--season", required=True, choices=["casp15", "casp16"])
    ap.add_argument(
        "--cache",
        help="raw table cache path (default temp/<season>/<rank file basename>)",
    )
    ap.add_argument("--out", help="default temp/<season>/top5_by_target.json")
    ap.add_argument("--top", type=int, default=5)
    args = ap.parse_args()

    cfg = get_season(args.season)
    cache = (
        Path(args.cache)
        if args.cache
        else Path("temp") / args.season / Path(cfg["rank_url"]).name
    )
    out = (
        Path(args.out)
        if args.out
        else Path("temp") / args.season / f"top{args.top}_by_target.json"
    )

    download(cfg["rank_url"], cache)
    text = cache.read_text(errors="replace")
    rows_by_target = PARSERS[cfg["rank_format"]](text)

    manifest = {t: top_n(rows, args.top) for t, rows in rows_by_target.items()}
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(__import__("json").dumps(manifest, indent=2) + "\n")
    print(f"Wrote {out} — {len(manifest)} targets ranked from {cache}")


if __name__ == "__main__":
    main()
