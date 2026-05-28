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
