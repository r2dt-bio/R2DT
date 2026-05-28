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

The command accepts the same options as [`pdb`](./pdb.md): `--mode`, `--basepairs`, `--format`, `--chain`, `--pseudoknots/--no-pseudoknots`, `--rnapuzzler`, and `-q/--quiet`. The 2D viewer renders base-pair annotations from FR3D, so leaving `--basepairs` at its default (or setting `--basepairs fr3d`) is recommended.

## Layout modes

The 2D layout can come from one of two pipelines:

- **`--mode templatefree`** (default, via `--mode auto`) — fast. R2DT uses the FR3D-derived dot-bracket and lays out the diagram with R2R, RNApuzzler, or RNArtist (auto-picked). Works on any RNA but gives an ad-hoc layout that may not match conventional textbook diagrams.
- **`--mode templated`** — slower but biologically meaningful. R2DT runs the full template search (CRW, RiboVision SSU/LSU, Rfam, GtRNAdb, RNase P, tmRNA, Rfam-tRNA) and renders the diagram on the matched template via Traveler. Recommended for rRNAs, tRNAs, and any RNA that has a curated template.

In both modes the FR3D base-pair overlay is identical — only the 2D *layout* differs. If `--mode templated` finds no matching template, the command fails with a hint to retry with `--mode templatefree` or `--mode auto`.

```bash
# Default (templatefree, ~seconds)
r2dt.py pdb_2d_3d 9RJA output/

# Templated layout (matches an SSU rRNA template)
r2dt.py pdb_2d_3d 9RJA output/ --mode templated
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

The same command works inside the published Docker image:

```bash
docker run --rm \
    -v $(pwd):/rna/r2dt \
    -w /rna/r2dt \
    rnacentral/r2dt:latest \
    ./r2dt.py pdb_2d_3d 9RJA output/9RJA_2d3d
```

The output appears under `output/9RJA_2d3d/viewer/` on your host. Pull-request builds (for example `rnacentral/r2dt:pr-219`) work the same way.

## Opening the viewer

The viewer fetches its data files via relative URLs, which browsers block when a page is opened directly from disk (`file://`). Serve the folder over HTTP instead:

```bash
python3 -m http.server -d <output_folder>/viewer 8000
# then open http://localhost:8000/
```

To publish it, upload the `viewer/` folder to any static host (GitHub Pages, Cloudflare Pages, an S3 bucket, etc.) — no server-side configuration is required.

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

# Specific chain
r2dt.py pdb_2d_3d 1S72 output/ --chain 9
```
