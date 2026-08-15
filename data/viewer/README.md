# Vendored 2D+3D viewer assets

These files are copied verbatim into every `viewer/` folder produced by
`r2dt.py pdb_2d_3d`. They are checked in here because the upstream build
artefacts are not published on a CDN or npm.

| File | Source | License |
| --- | --- | --- |
| `pdb-rna-viewer-plugin-0.3.0.js` | [PDBeurope/pdb-rna-viewer](https://github.com/PDBeurope/pdb-rna-viewer) v0.3.0 release | Apache-2.0 (see `pdb-rna-viewer.LICENSE`) |
| `pdb-rna-viewer-0.3.0.css` | same | Apache-2.0 |
| `r2dt-2d-3d-viewer.js` | R2DT (this repository) | Apache-2.0 |
| `r2dt-2d-3d-viewer.css` | R2DT (this repository) | Apache-2.0 |

`pdb-rna-viewer.LICENSE` is the upstream license verbatim: Apache-2.0 plus an
EMBL-EBI acknowledgement clause requiring that reuse acknowledge the original
source — which this README and the generated viewer pages do. The bundled
`pdb-rna-viewer-plugin-0.3.0.js.LICENSE.txt` is the tslib banner webpack
extracted from the upstream build, kept alongside as shipped.

`r2dt-2d-3d-viewer.js` exposes `R2DTViewer.create()` — the glue that wires
pdb-rna-viewer (2D) to pdbe-molstar (3D). Generated `index.html` calls it
inline; third-party pages can load the same scripts and call it with a
`baseUrl` pointing at a published artifact folder.

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
