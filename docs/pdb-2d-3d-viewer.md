# Interactive 2D + 3D viewer

The `pdb_2d_3d` command runs the standard [`pdb`](./pdb.md) pipeline and additionally assembles a `viewer/` folder. Served over HTTP, `viewer/index.html` shows the R2DT 2D diagram and a 3D molstar view side by side, with clicks linked between them.

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
├── index.html                          # the viewer page
├── viewer.js                           # interaction glue (vendored in R2DT)
├── api.json                            # nucleotide positions, sequence, seq-id maps
├── fr3d.json                           # base-pair annotations
├── <structure_id>.cif (or .pdb)        # local copy of the structure
├── pdb-rna-viewer-plugin-0.3.0.js      # 2D viewer plugin (vendored in R2DT)
└── pdb-rna-viewer-0.3.0.css            # 2D viewer stylesheet (vendored in R2DT)
```

`api.json` and `fr3d.json` are derived from R2DT's `*.colored.json`, `*.colored.svg`, and `*_basepair.txt` outputs and use 1-based positions aligned to the full (mask-expanded) deposited sequence, so unresolved residues are preserved.

R2DT ships the [pdb-rna-viewer](https://github.com/PDBeurope/pdb-rna-viewer) v0.3.0 build files (plus the `viewer.js` glue) under `data/viewer/` and copies them into each output folder automatically — no separate setup is needed. The 3D pane loads a pinned `pdbe-molstar` build from jsDelivr, so opening the viewer requires a network connection for that script.

The folder is self-contained relative to its own location: every data file is fetched via a relative URL, so it works under any same-origin static host — a local `http.server`, GitHub Pages, or Cloudflare Pages — with no configuration.

## Running via Docker

`pdb_2d_3d` is not in a released image yet, so use the pull-request build that has the code baked in (`rnacentral/r2dt:pr-219`). Mount an output directory to get the results back on your host:

```bash
docker run --rm \
    -v $(pwd)/output:/rna/r2dt/output \
    -w /rna/r2dt \
    rnacentral/r2dt:pr-219 \
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

## Galleries

Several `viewer/` folders can be combined into a browsable gallery with `utils/build_viewers.py`, which writes an `index.html` of 2D-diagram thumbnails linking to each interactive viewer. Two galleries are wired up as `just` recipes:

- **Main gallery** — `just viewers` regenerates a gallery of curated PDB structures (FR3D base pairs) under `output/site/`.
- **Workstream 1 dashboard** — `just ws1-viewers` regenerates a gallery that uses [`--basepairs cif`](./pdb.md#choosing-a-base-pair-extractor) on the FR3D-converted mmCIF inputs from the [na-hackathon](https://github.com/na-hackathon/na-hackathon-2026) repo. It downloads the inputs, runs `pdb_2d_3d --basepairs cif` on each, and publishes under `output/site/workstream1/` — beside, but separate from, the main gallery.

`just viewers-deploy` publishes all of `output/site/` (both galleries) to Cloudflare Pages, so the workstream1 dashboard lands at `/workstream1/`. Run `just viewers` and/or `just ws1-viewers` first to refresh the contents. Set `CLOUDFLARE_PROJECT` in a gitignored `.env`.

## Interaction model

- **2D → 3D.** Clicking a nucleotide in the 2D diagram selects and focuses the corresponding residue in the 3D view.
- **2D → 3D (base pairs).** Clicking a base-pair line selects both partner residues in 3D and highlights the line in orange (the same colour the plugin uses for a clicked nucleotide). Only one base pair stays highlighted at a time.
- **3D → 2D.** Clicking or hovering a residue in molstar highlights the matching nucleotide in the 2D diagram.
- **Hover does nothing in 2D.** Hovering nucleotides or base pairs in the 2D diagram has no effect — only clicks drive a 3D response.
- **Backbone path overlay.** A faint line tracing the backbone is drawn under the nucleotide letters, on by default. A "Show backbone path" checkbox (in place of the plugin's "View as Path" dropdown) toggles it.
- **Base-pair family filter.** Only families actually present in this structure are listed in the filter, and they are all enabled by default.
- **Unresolved residues are dimmed.** Nucleotides missing from the 3D coordinates are shown in grey, consistent with the `pdb` command's behaviour (see [missing nucleotides in PDB structures](./pdb.md#missing-nucleotides-in-pdb-structures)).
- **Pseudoknot Watson–Crick pairs are lightened** so they don't dominate the nested cWW ladder.

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
