# Vendored 2D+3D viewer assets

These files are copied verbatim into every `viewer/` folder produced by
`r2dt.py pdb_2d_3d`. They are checked in here because the upstream build
artefacts are not published on a CDN or npm.

| File | Source | License |
| --- | --- | --- |
| `pdb-rna-viewer-plugin-0.3.0.js` | [PDBeurope/pdb-rna-viewer](https://github.com/PDBeurope/pdb-rna-viewer) v0.3.0 release | Apache-2.0 (see `*.LICENSE.txt`) |
| `pdb-rna-viewer-0.3.0.css` | same | Apache-2.0 |
| `viewer.js` | R2DT (this repository) | Apache-2.0 |

`viewer.js` is the glue that wires pdb-rna-viewer (2D) to pdbe-molstar
(3D); it is plain JavaScript and reads per-structure configuration from
`window.R2DT_CONFIG`, injected by the generated `index.html`.

`pdbe-molstar` itself is loaded at runtime from a pinned jsDelivr URL (see
`utils/viewer_html.py`), not vendored here.

To upgrade pdb-rna-viewer, download the new release's `build/` files,
replace the two filenames above (updating the version in the names and in
`utils/viewer_html.py`), and re-test the viewer.

## Local modification to the plugin

`pdb-rna-viewer-plugin-0.3.0.js` carries one local patch that must be
re-applied after any upgrade. Upstream draws a filled Leontis–Westhof circle
only for cWW **G-U/U-G** wobbles and a plain line for every other cWW pair --
which makes non-canonical cWW pairs (e.g. U-U) look identical to canonical
Watson–Crick pairs. We broadened the condition so the filled circle is drawn
for *any* cWW pair that is **not** a canonical WC pair (A-U, U-A, G-C, C-G,
A-T, T-A); canonical pairs keep the plain line. Find `if("cWW"==s)` and the
condition immediately after it:

```js
// upstream:
"G"==l&&"U"==u||"U"==l&&"G"==u
// patched:
!("A"==l&&"U"==u||"U"==l&&"A"==u||"G"==l&&"C"==u||"C"==l&&"G"==u||"A"==l&&"T"==u||"T"==l&&"A"==u)
```
