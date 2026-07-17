#!/usr/bin/env python3
"""Render a static CASP dashboard (``<site>/index.html``) from ``results.json``.

Season-agnostic (pass ``--title``, e.g. "CASP16" or "CASP15"). Self-contained:
the results are embedded in the page and rendered client-side, so it works both
over an HTTP server and from ``file://``.  The table is sortable (click a
header) and filterable (type in the box); each row links to that
reference/model pair's compare viewer.
"""
import argparse
import json
from pathlib import Path

_PAGE = """<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title} — reference vs model dashboard</title>
<style>
  :root {{ color-scheme: light dark; }}
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
          margin: 0; padding: 24px; color: #1f2937; background: #fff; }}
  h1 {{ font-size: 20px; margin: 0 0 4px; }}
  .sub {{ color: #6b7280; font-size: 13px; margin: 0 0 12px; }}
  .legend {{ color: #6b7280; font-size: 12px; margin: 0 0 16px; line-height: 1.5; }}
  .legend b {{ color: inherit; font-weight: 600; }}
  .toolbar {{ display: flex; gap: 12px; align-items: center; margin: 0 0 12px; flex-wrap: wrap; }}
  input[type=search] {{ padding: 7px 10px; font-size: 13px; border: 1px solid #d1d5db;
                        border-radius: 6px; min-width: 220px; background: #fff; color: inherit; }}
  .count {{ color: #6b7280; font-size: 12px; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 13px; }}
  th, td {{ text-align: left; padding: 7px 10px; border-bottom: 1px solid #e5e7eb; white-space: nowrap; }}
  th {{ position: sticky; top: 0; background: #f9fafb; cursor: pointer; user-select: none;
       font-weight: 600; }}
  th.num, td.num {{ text-align: right; font-variant-numeric: tabular-nums; }}
  td.tid {{ color: #6b7280; font-variant-numeric: tabular-nums; }}
  tr:hover td {{ background: rgba(0,0,0,0.04); }}
  a {{ color: inherit; text-decoration: underline; text-underline-offset: 2px; }}
  a:hover {{ text-decoration: underline; opacity: 0.7; }}
  .metric {{ font-weight: 600; font-variant-numeric: tabular-nums; }}
  .badge {{ font-size: 11px; padding: 1px 6px; border-radius: 10px;
            background: rgba(0,0,0,0.08); color: inherit; }}
  .diff {{ font-variant-numeric: tabular-nums; }}
  @media (prefers-color-scheme: dark) {{
    body {{ background: #0b0f16; color: #e5e7eb; }}
    th {{ background: #111827; }}
    th, td {{ border-color: #1f2937; }}
    input[type=search] {{ background: #111827; border-color: #374151; }}
    tr:hover td {{ background: rgba(255,255,255,0.05); }}
    .badge {{ background: rgba(255,255,255,0.12); }}
  }}
</style>
</head>
<body>
<h1>{title} — reference vs model</h1>
<p class="sub">{summary}</p>
<p class="legend">
  <b>BP diff (M/L/A)</b> — base-pair differences between the predicted model and the
  experimental reference:
  <b>M</b> = matched (pairs present in both),
  <b>L</b> = lost (in the reference but missing from the model — false negatives),
  <b>A</b> = added (only in the model — false positives).
</p>
<div class="toolbar">
  <input id="q" type="search" placeholder="Filter by target or model…" autocomplete="off">
  <span class="count" id="count"></span>
</div>
<table id="tbl">
  <thead><tr>
    <th data-k="target" title="Sequential target number">#</th>
    <th data-k="target">Target</th>
    <th data-k="rank" class="num">Rank</th>
    <th data-k="model">Model</th>
    <th data-k="tm" class="num">TM-score</th>
    <th data-k="inf_wc" class="num">INF WC</th>
    <th data-k="inf_nwc" class="num">INF non-WC</th>
    <th data-k="inf_all" class="num">INF all</th>
    <th data-k="rmsd" class="num">RMSD (Å)</th>
    <th>BP diff (M/L/A)</th>
    <th></th>
  </tr></thead>
  <tbody id="rows"></tbody>
</table>
<script id="data" type="application/json">{data_json}</script>
<script>
(function () {{
  var RAW = JSON.parse(document.getElementById('data').textContent).results || [];
  // Sequential target id (1..N) from the sorted unique targets, so each target
  // has a stable short number regardless of the current sort order.
  var TNUM = {{}};
  Array.prototype.slice.call(RAW.map(function (r) {{ return r.target; }}))
    .filter(function (t, i, a) {{ return a.indexOf(t) === i; }})
    .sort()
    .forEach(function (t, i) {{ TNUM[t] = i + 1; }});
  var rows = RAW.map(function (r) {{
    var inf = r.inf || {{}};
    return {{
      tid: TNUM[r.target], target: r.target, rank: r.rank, model: r.model, tm: r.tm,
      inf_wc: inf.wc, inf_nwc: inf.nwc, inf_all: inf.all,
      rmsd: r.rmsd, diff: r.diff, status: r.status, page: r.page, error: r.error,
    }};
  }});
  function num(v, d) {{ return v == null ? '' : Number(v).toFixed(d); }}
  function metricCell(v, d) {{
    return v == null ? '<td class="num"></td>'
      : '<td class="num"><span class="metric">' + num(v, d) + '</span></td>';
  }}
  var state = {{ key: 'target', dir: 1, q: '' }};
  function render() {{
    var q = state.q.toLowerCase();
    var view = rows.filter(function (r) {{
      return !q || (r.target + ' ' + r.model).toLowerCase().indexOf(q) >= 0;
    }});
    view.sort(function (a, b) {{
      var x = a[state.key], y = b[state.key];
      if (x == null) return 1; if (y == null) return -1;
      if (typeof x === 'string') return state.dir * x.localeCompare(y);
      return state.dir * (x - y);
    }});
    document.getElementById('rows').innerHTML = view.map(function (r) {{
      var d = r.diff || {{}};
      var diffCell = r.diff
        ? '<td class="diff">' + d.matched + '/' + d.lost + '/' + d.added + '</td>'
        : '<td></td>';
      var canView = r.status !== 'failed' && r.page;
      // Target links to its own compare viewer page (same as the View link).
      var targetCell = canView
        ? '<td><a href="' + r.page + '">' + r.target + '</a></td>'
        : '<td>' + r.target + '</td>';
      var last = canView
        ? '<td><a href="' + r.page + '">View →</a></td>'
        : '<td><span class="badge" title="' + (r.error || '').replace(/"/g, '&quot;') + '">failed</span></td>';
      return '<tr>'
        + '<td class="tid">' + (r.tid == null ? '' : r.tid) + '</td>'
        + targetCell
        + '<td class="num">' + (r.rank == null ? '' : r.rank) + '</td>'
        + '<td>' + r.model + '</td>'
        + metricCell(r.tm, 3)
        + metricCell(r.inf_wc, 3) + metricCell(r.inf_nwc, 3) + metricCell(r.inf_all, 3)
        + '<td class="num">' + num(r.rmsd, 2) + '</td>'
        + diffCell + last + '</tr>';
    }}).join('');
    document.getElementById('count').textContent = view.length + ' of ' + rows.length + ' pairs';
  }}
  document.querySelectorAll('th[data-k]').forEach(function (th) {{
    th.addEventListener('click', function () {{
      var k = th.getAttribute('data-k');
      state.dir = state.key === k ? -state.dir : 1;
      state.key = k; render();
    }});
  }});
  document.getElementById('q').addEventListener('input', function (e) {{
    state.q = e.target.value; render();
  }});
  render();
}})();
</script>
</body>
</html>
"""


def main():
    """Render ``<site>/results.json`` into a self-contained ``<site>/index.html``."""
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--site", default="site", help="site dir containing results.json")
    ap.add_argument("--title", default="CASP", help='page title/heading, e.g. "CASP16"')
    args = ap.parse_args()

    site = Path(args.site)
    data = json.loads((site / "results.json").read_text())
    results = data.get("results", [])
    targets = sorted({r["target"] for r in results})
    n_ok = sum(1 for r in results if r.get("status") in ("ok", "cached"))
    summary = (
        f"{len(targets)} targets · {len(results)} reference/model pairs "
        f"({n_ok} available). Click a column header to sort."
    )
    html = _PAGE.format(
        title=args.title, summary=summary, data_json=json.dumps({"results": results})
    )
    out = site / "index.html"
    out.write_text(html)
    print(f"Wrote {out} — {len(targets)} targets, {len(results)} pairs")


if __name__ == "__main__":
    main()
