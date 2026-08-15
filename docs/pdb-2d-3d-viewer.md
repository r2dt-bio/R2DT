# Interactive 2D + 3D viewer

The `pdb_2d_3d` command runs the standard [`pdb`](./pdb.md) pipeline and additionally assembles a `viewer/` folder. Served over HTTP, `viewer/index.html` shows the R2DT 2D diagram and a 3D molstar view side by side, with clicks linked between them.

For an interactive local UI that generates these viewers (and compare / edit / export workflows), see the [local workstation](./workstation.md).

## Basic usage

```bash
# Using a PDB ID (downloaded from RCSB)
r2dt.py pdb_2d_3d <PDB-ID> <output_folder>

# Using a local PDB or mmCIF file
r2dt.py pdb_2d_3d <structure_file> <output_folder>
```

For example:

```bash
r2dt.py pdb_2d_3d 1Y26 output/
```

The command accepts the same options as [`pdb`](./pdb.md): `--mode`, `--basepairs`, `--format`, `--chain`, `--pseudoknots/--no-pseudoknots`, `--rnapuzzler`, and `-q/--quiet`. By default the 2D viewer renders base-pair annotations computed by FR3D. You can also read pairs straight from a DNATCO/NDB-annotated mmCIF with [`--basepairs cif`](./pdb.md#choosing-a-base-pair-extractor) (no FR3D run); the viewer header credits whichever source was used.

## Layout modes

The `--mode` option controls how the 2D layout is produced:

- **`--mode auto`** (default) — try templated first, fall back to templatefree. R2DT runs the template search; if a template matches it is used, otherwise R2DT automatically falls back to the templatefree layout. This always produces a diagram and gives the best layout available for the structure.
- **`--mode templated`** — templated only. R2DT runs the full template search (CRW, RiboVision SSU/LSU, Rfam, GtRNAdb, RNase P, tmRNA, Rfam-tRNA) and renders the diagram on the matched template via Traveler. Biologically meaningful and recommended for rRNAs, tRNAs, and anything with a curated template, but **fails** (with a hint) if no template matches.
- **`--mode templatefree`** — templatefree only. R2DT uses the FR3D-derived dot-bracket and lays out the diagram with R2R, RNApuzzler, or RNArtist (auto-picked). Fast, works on any RNA, but gives an ad-hoc layout that may not match conventional textbook diagrams.

The mode only affects the 2D *layout* — the base-pair overlay (from whichever `--basepairs` source) is identical in all three.

```bash
# Default: templated if a template matches, otherwise templatefree
r2dt.py pdb_2d_3d 9RJA output/

# Force the templated layout (errors if no template matches)
r2dt.py pdb_2d_3d 9RJA output/ --mode templated

# Force the quick templatefree layout
r2dt.py pdb_2d_3d 9RJA output/ --mode templatefree
```

The `--rnapuzzler` flag only applies to the templatefree path; in templated mode the layout comes from the matched template.

## Output layout

In addition to everything the `pdb` command writes, the output folder gets a `viewer/` directory:

```
<output_folder>/viewer/
├── index.html                          # viewer page (calls R2DTViewer.create)
├── manifest.json                       # structure id, chain, file URL (for embeds)
├── r2dt-2d-3d-viewer.js                # R2DTViewer.create() glue (vendored in R2DT)
├── r2dt-2d-3d-viewer.css               # toolbar / panel chrome (vendored in R2DT)
├── api.json                            # nucleotide positions, sequence, seq-id maps
├── fr3d.json                           # base-pair annotations
├── lbn.json                            # layered dot-bracket notation rows
├── <structure_id>.cif (or .pdb)        # local copy of the structure
├── pdb-rna-viewer-plugin-0.3.0.js      # 2D viewer plugin (vendored in R2DT)
└── pdb-rna-viewer-0.3.0.css            # 2D viewer stylesheet (vendored in R2DT)
```

`api.json` and `fr3d.json` are derived from R2DT's `*.colored.json`, `*.colored.svg`, and `*_basepair.txt` outputs and use 1-based positions aligned to the full (mask-expanded) deposited sequence, so unresolved residues are preserved.

R2DT ships the [pdb-rna-viewer](https://github.com/PDBeurope/pdb-rna-viewer) v0.3.0 build files (plus the `r2dt-2d-3d-viewer.js` glue) under `data/viewer/` and copies them into each output folder automatically — no separate setup is needed. The 3D pane loads a pinned `pdbe-molstar` build from jsDelivr, so opening the viewer requires a network connection for that script.

The folder is self-contained relative to its own location: every data file is fetched via a relative URL, so it works under any same-origin static host — a local `http.server`, GitHub Pages, or Cloudflare Pages — with no configuration.

## Comparing a reference structure against a predicted model

`pdb` also has a compare mode that produces a 3-panel viewer instead of the single-structure one above: the reference's 2D diagram, a predicted model's 2D diagram (drawn on the *same* combined layout, so equivalent residues line up), and one shared Mol* pane overlaying both 3D structures. It's the same `R2DTViewer.createCompare()` API used in [Comparing two 2D diagrams with one 3D view](#comparing-two-2d-diagrams-with-one-3d-view) below — this command is simply R2DT's own generator for that pattern, and it's what CASP-style reference/prediction dashboards are built from.

Compare mode needs an mmCIF reference and an explicit chain selection (multi-chain diagrams only work from mmCIF):

```bash
# Reference vs a real predicted model (same sequence, same chain order)
r2dt.py pdb reference.cif output/ --chains A --compare --model model.pdb

# --model implies --compare; --model-chains maps the model's own chain ids
# when they differ from the reference's (default: same order as --chains)
r2dt.py pdb reference.cif output/ --chains A,B --model model.cif --model-chains X,Y

# Without --model, --compare alone shows a randomly perturbed copy of the
# reference standing in for the model (useful for previewing the diff UI)
r2dt.py pdb reference.cif output/ --chains A --compare
```

Output layout adds to `viewer/`:

```
<output_folder>/viewer/
├── index.html                          # calls R2DTViewer.createCompare()
├── ref/                                 # reference panel's data
│   ├── api.json
│   ├── fr3d.json
│   ├── lbn.json
│   └── bp-compare.json                 # base-pair keys for the model panel's TP/FN badges
├── model/                               # model panel's data (same shape as ref/)
│   ├── api.json
│   ├── fr3d.json
│   ├── lbn.json
│   ├── bp-compare.json                 # base-pair keys for the reference panel's TP/FP badges
│   └── label-maps.json                 # 2D-label → model author residue/chain, for 3D click-through
├── api.json, fr3d.json                 # root-level copy of ref/'s data, used
│                                        # by the shared Mol* pane
├── metrics.json                        # INF + matched/lost/added base pairs
├── inf-pairs.json                      # INF scores + reference/model pair lists (download)
├── inf-pairs.csv                       # same data, spreadsheet-friendly
├── <reference_id>.cif                  # reference structure
├── <model_id>.aligned.cif              # model superposed onto the reference
└── … the same vendored viewer assets as above
```

As with the single-structure viewer, every file is fetched via a relative URL — `index.html` itself carries no structural data, so hosting a compare page is the same "upload the folder" deal as [Embedding on your own site](#embedding-on-your-own-site), just with `ref/`/`model/` added alongside the structure files.

## Running via Docker

`pdb_2d_3d` is not in a released image yet, so use the development build that has the code baked in (`rnacentral/r2dt:develop`; switch to `rnacentral/r2dt:latest` once a release includes it). Mount an output directory to get the results back on your host:

```bash
docker run --rm \
    -v $(pwd)/output:/rna/r2dt/output \
    -w /rna/r2dt \
    rnacentral/r2dt:develop \
    ./r2dt.py pdb_2d_3d 9RJA output/9RJA_2d3d
```

The output appears under `output/9RJA_2d3d/viewer/` on your host.

:::{note}
Don't use `rnacentral/r2dt:latest` yet — it predates this command. Once the feature is merged and released, `:latest` will include it. If you have a local checkout of this feature's branch, you can instead mount it over the image (`-v $(pwd):/rna/r2dt`) and run against any recent tag, since the mount supplies the code.
:::

## Opening the viewer

The viewer fetches its data files via relative URLs, which browsers block when a page is opened directly from disk (`file://`). Serve the folder over HTTP instead:

```bash
python3 -m http.server -d <output_folder>/viewer 8000
# then open http://localhost:8000/
```

To publish it, upload the `viewer/` folder to any static host (GitHub Pages, Cloudflare Pages, an S3 bucket, etc.) — no server-side configuration is required.

## Embedding on your own site

Generated `index.html` is itself an embed page: it loads the vendored scripts and calls `R2DTViewer.create()` to mount the viewer into a `<div>`. To put the same viewer on **your** website, host the whole `viewer/` folder on HTTPS and replicate that bootstrap on your page.

### 1. Host the artifact folder

Upload every file from `viewer/` to a static URL, for example:

```
https://cdn.example.com/rna/9SFQ/
├── index.html
├── manifest.json
├── r2dt-2d-3d-viewer.js
├── r2dt-2d-3d-viewer.css
├── api.json
├── fr3d.json
├── lbn.json
├── 9SFQ.cif
├── pdb-rna-viewer-plugin-0.3.0.js
└── pdb-rna-viewer-0.3.0.css
```

All paths are relative, so the folder can live at any prefix as long as the files stay together. Visitors' browsers also need network access to the pinned **pdbe-molstar** script on jsDelivr (same as the standalone page).

### 2. Add a mount point and load the scripts

On your page, include the same stylesheets and scripts as `index.html`, then call the API:

```html
<link rel="stylesheet" href="https://cdn.example.com/rna/9SFQ/pdb-rna-viewer-0.3.0.css">
<link rel="stylesheet" href="https://cdn.example.com/rna/9SFQ/r2dt-2d-3d-viewer.css">
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/pdbe-molstar@3.12.0/build/pdbe-molstar-light.css">

<div id="my-rna-viewer"></div>

<script src="https://cdn.example.com/rna/9SFQ/pdb-rna-viewer-plugin-0.3.0.js"></script>
<script src="https://cdn.jsdelivr.net/npm/pdbe-molstar@3.12.0/build/pdbe-molstar-plugin.js"></script>
<script src="https://cdn.example.com/rna/9SFQ/r2dt-2d-3d-viewer.js"></script>
<script>
R2DTViewer.create({
  mount: '#my-rna-viewer',
  baseUrl: 'https://cdn.example.com/rna/9SFQ/',
  structureId: '9SFQ',
  chainId: 'T',
  structureUrl: './9SFQ.cif',
  structureFormat: 'cif',
}).catch(function (err) {
  console.error(err);
});
</script>
```

Copy the `<script>` block verbatim from a generated `index.html` if you prefer — it already contains the correct `structureId`, `chainId`, and file URLs for that structure.

When `baseUrl` points at a hosted artifact folder, `structureId` / `chainId` / `structureUrl` can be omitted if `manifest.json` is present; the API reads them from there.

### API reference

`R2DTViewer.create(options)` returns a Promise that resolves to a handle:

| Option | Required | Description |
| --- | --- | --- |
| `mount` | yes | CSS selector or DOM element to mount into |
| `baseUrl` | yes | URL of the artifact folder (trailing slash optional) |
| `structureId` | if no manifest | PDB / mmCIF identifier |
| `chainId` | no | Author chain id (default: from manifest or `""`) |
| `structureUrl` | no | Relative path to structure file (default: from manifest or `<id>.cif`) |
| `structureFormat` | no | `"cif"` or `"pdb"` (default: from manifest or `"cif"`) |
| `layout` | no | `"side-by-side"` (default), `"stacked"`, `"2d-only"`, `"3d-only"` |
| `height` | no | Panel height, e.g. `640` or `"100%"` |
| `panelWidth` | no | Width of each pane in px (default `600`) |
| `showLbn` | no | Show layered dot-bracket panel (default `true`) |
| `showLegend` | no | Show Leontis–Westhof legend in the base-pairs filter (default `true`) |
| `onReady` | no | Callback `(handle) => {}` when init completes |

Handle methods:

- `handle.selectResidue(label)` — highlight a 1-based sequence position in 2D and 3D
- `handle.selectBasePair(a, b)` — highlight a base pair
- `handle.destroy()` — tear down the viewer and clear the mount point

:::{note}
**v1 limitation:** only one `R2DTViewer.create()` or `R2DTViewer.createCompare()` call per page. To show multiple independent structures on a dashboard, link to separate viewer pages or use `<iframe src="…/index.html">` for each structure.
:::

### Comparing two 2D diagrams with one 3D view

`R2DTViewer.createCompare()` mounts two (or more) 2D panels side by side plus a shared molstar pane. Clicks in one designated 2D panel are linked to the 3D view — useful when comparing two layouts of related RNAs while keeping a single 3D structure in frame.

Host one folder per structure (each with `api.json`, `fr3d.json`, and the structure file), then bootstrap from the parent page:

```html
<div id="rna-compare"></div>
<script>
R2DTViewer.createCompare({
  mount: '#rna-compare',
  panels: [
    {
      title: '8SH5',
      subtitle: '— 2D (chain R)',
      baseUrl: './8sh5/',
      structureId: '8SH5',
      chainId: 'R',
      structureUrl: './8SH5.pdb',
      structureFormat: 'pdb',
    },
    {
      title: 'RNApolis_R1293_m2',
      subtitle: '— 2D (chain 1)',
      baseUrl: './rnapolis/',
      structureId: 'RNApolis_R1293_m2',
      chainId: '1',
      structureUrl: './RNApolis_R1293_m2.pdb',
      structureFormat: 'pdb',
    },
  ],
  molstar: {
    panelIndex: 0,
    title: '8SH5',
    subtitle: '— 3D (molstar)',
    baseUrl: './8sh5/',
    structureId: '8SH5',
    chainId: 'R',
    structureUrl: './8SH5.pdb',
    structureFormat: 'pdb',
  },
});
</script>
```

| Option | Required | Description |
| --- | --- | --- |
| `mount` | yes | CSS selector or DOM element |
| `panels` | yes | Array of panel configs (`title`, `subtitle`, `baseUrl`, `structureId`, `chainId`, …) |
| `molstar` | no | Shared 3D pane; `panelIndex` (default `0`) selects which 2D panel links to 3D |
| `panelHeight` | no | Height of each pane, e.g. `480` or `"50vh"` |
| `fetchShim` | no | Route pdb-rna-viewer's EBI fetches to local `api.json` / `fr3d.json` (default `true`) |
| `onReady` | no | Callback `(handle) => {}` when init completes |

`utils/viewer_html.render_compare()` writes an `index.html` that calls `createCompare()` — see `output/compare/` for a hand-built working example, or generate a real one with [`pdb --compare --model`](#comparing-a-reference-structure-against-a-predicted-model) above, whose `viewer/ref/` and `viewer/model/` folders are exactly the per-panel `baseUrl` folders this API expects.

### iframe fallback

If you do not want to load the scripts on your own page, embed the hosted viewer page directly:

```html
<iframe
  src="https://cdn.example.com/rna/9SFQ/index.html"
  width="1232"
  height="700"
  style="border:0"
  allow="fullscreen"
></iframe>
```

This works without JavaScript integration but gives you a fixed-size page-in-a-page rather than a native widget.

## Galleries

Several `viewer/` folders can be combined into a browsable gallery with `utils/build_viewers.py`, which writes an `index.html` of 2D-diagram thumbnails linking to each interactive viewer. Two galleries are wired up as `just` recipes:

- **Main gallery** — `just viewers` regenerates a gallery of curated PDB structures (FR3D base pairs) under `output/site/`.
- **Workstream 1 dashboard** — `just ws1-viewers` regenerates a gallery that uses [`--basepairs cif`](./pdb.md#choosing-a-base-pair-extractor) on the FR3D-converted mmCIF inputs from the [na-hackathon](https://github.com/na-hackathon/na-hackathon-2026) repo. It downloads the inputs, runs `pdb_2d_3d --basepairs cif` on each, and publishes under `output/site/workstream1/` — beside, but separate from, the main gallery.

`just viewers-deploy` publishes all of `output/site/` (both galleries) to Cloudflare Pages — for example [https://2d-3d-viewer.na-hackathon-2026.pages.dev](https://2d-3d-viewer.na-hackathon-2026.pages.dev) — so the workstream1 dashboard lands at `/workstream1/`. Run `just viewers` and/or `just ws1-viewers` first to refresh the contents. Set `CLOUDFLARE_PROJECT` in a gitignored `.env`.

## Interaction model

- **2D → 3D.** Clicking a nucleotide in the 2D diagram selects and focuses the corresponding residue in the 3D view.
- **2D → 3D (base pairs).** Clicking a base-pair line selects both partner residues in 3D. The line keeps its TP/FP/FN (or default) colour and is marked selected with a thicker stroke plus a brighter same-hue glow; selected letters stay black with a pale structure-tinted badge (reference green / model blue in compare mode). The matching row in the base-pair list is outlined. Only one base pair stays highlighted at a time.
- **Compare TP/FP/FN colours.** In compare mode, true-positive pairs are green, false negatives blue, and false positives red — the same colours in the base-pair list, 2D diagram, and LBN. Selection does not recolour those strokes. Reference green / model blue on selection badges mark which panel is active; pair strokes keep TP/FN meaning.
- **Multi-chain pair list.** When more than one chain is present, the base-pair list is grouped into per-chain and between-chain sections; the INF bar can expand a by-chain score breakdown. **Download scores** (JSON / CSV) exports `inf-pairs.json` / `inf-pairs.csv` with the scores and the underlying pair lists.
- **3D → 2D.** Clicking or hovering a residue in molstar highlights the matching nucleotide in the 2D diagram.
- **Hover does nothing in 2D.** Hovering nucleotides or base pairs in the 2D diagram has no effect — only clicks drive a 3D response.
- **Backbone path overlay.** A faint line tracing the backbone is drawn under the nucleotide letters, on by default. A "Backbone" toggle in the toolbar shows or hides it.
- **Base-pair family filter.** Only families actually present in this structure are listed in the filter, and they are all enabled by default. In compare mode each panel's filter is scoped to that panel (duplicate plugin checkbox ids are not shared across panels).
- **Layered dot-bracket notation (LBN).** Below the 2D+3D panes, a scrollable panel shows the sequence and per-family dot-bracket rows; clicks highlight the corresponding residues in 2D and 3D.
- **Unresolved residues are dimmed.** Nucleotides missing from the 3D coordinates are shown in grey, consistent with the `pdb` command's behaviour (see [missing nucleotides in PDB structures](./pdb.md#missing-nucleotides-in-pdb-structures)).
- **Pseudoknot Watson–Crick pairs are lightened** so they don't dominate the nested cWW ladder.
- **Short-span LW glyphs.** Base-pair symbols between nearby nucleotides (e.g. adjacent stacking contacts) are shrunk and nudged off the letters so the sequence stays readable.

## Examples

```bash
# tRNA-Phe
r2dt.py pdb_2d_3d 1EHZ output/

# Local mmCIF
r2dt.py pdb_2d_3d ./my_rna.cif output/ --basepairs fr3d

# DNATCO/NDB-annotated mmCIF (pairs read from the file, no FR3D run)
r2dt.py pdb_2d_3d ./my_rna_dnatco.cif output/ --basepairs cif

# Specific chain
r2dt.py pdb_2d_3d 1S72 output/ --chain 9
```
