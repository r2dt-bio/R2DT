# Processing Stockholm alignments

R2DT can process Stockholm-format multiple sequence alignments that contain named secondary structure regions. This is particularly useful for viral genome alignments where different RNA structures have been manually annotated with names.

## Overview

Many viral RNA databases use Stockholm alignments with special `#=GC` annotation lines to mark named structural elements. R2DT uses the following annotation lines:

- **#=GC SS_cons** — Consensus secondary structure in dot-bracket notation
- **#=GC structureID** — Individual structure names (e.g. SLI, IRES), separated by `|`
- **#=GC regionID** — Parent genomic region names (e.g. 5'UTR, NS5B), separated by `|`

R2DT extracts each named region, computes an RF-style consensus sequence using IUPAC ambiguity codes, and generates template-free visualizations for each structure. The outputs can be automatically stitched into a single combined diagram.

## Quick start

```bash
r2dt.py stockholm examples/hcv-alignment.stk output/hcv/
```

This command:
1. Parses the Stockholm alignment
2. Extracts named structural elements from `structureID`
3. Assigns parent regions from `regionID`
4. Computes RF consensus sequences with IUPAC codes (R for purine, Y for pyrimidine, etc.)
5. Generates template-free 2D diagrams for each valid region
6. Stitches all diagrams into a combined view

## Input format

The Stockholm file must contain `#=GC SS_cons` and at least one of:
- `#=GC structureID` — for multi-region alignments (e.g. viral genomes)
- `#=GC knownSS_names` — legacy format equivalent
- Neither — the entire alignment is treated as a single structure (e.g. Rfam seed alignments)

### structureID + regionID

The format uses two `#=GC` annotation lines:

- **`structureID`** labels individual structural elements (stem-loops, junctions, etc.)
- **`regionID`** labels the broader genomic region each structure belongs to

```text
# STOCKHOLM 1.0

seq1/1-200      AUCGAUCG...AUCGAUCG...AUCGAUCG...
seq2/1-200      AUCGAUCG...AUCGAUCG...AUCGAUCG...

#=GC SS_cons         ...(((...)))......(((....)))......(((.....)))...
#=GC structureID     |...SLI.........|......SLII......|.....SLIII...|
#=GC regionID        |.........5'UTR...................|..core_protein.|
//
```

Both lines use the same pipe-delimited format:

- `|` characters mark boundaries between named segments
- Text between pipes (after stripping dots and whitespace) becomes the name
- Dots (`.`) are filler characters — only the text content matters
- Columns in `structureID` / `regionID` correspond 1:1 with alignment and `SS_cons` columns

In this example, `SLI` and `SLII` both fall within the `5’UTR` region, while `SLIII` belongs to `core_protein`. R2DT uses the midpoint of each structure to determine its parent region.


### Format rules

The `|` characters in the annotation lines mark the boundaries between regions:

| Feature | Description |
|---------|-------------|
| Delimiters | Pipe `\|` characters separate named segments |
| Filler | Dots `.` pad segments to the correct column width |
| Names | Non-dot, non-whitespace text between pipes is the name |
| Empty segments | Segments with only dots/whitespace are unnamed (skipped) |
| Column alignment | Each character position corresponds to an alignment column |

### Simple alignment (e.g. Rfam seed)

When neither `structureID` nor `knownSS_names` is present, R2DT treats the entire alignment as a single structure. This works out of the box for Rfam seed alignments and any other Stockholm file with `SS_cons`.

```text
# STOCKHOLM 1.0
#=GF ID SAM
#=GF AC RF00162

seq1/1-108      CUCUUAUCAAGAG...
seq2/1-108      ACCUUAUUUUGAG...

#=GC SS_cons    (.(((((((,,,,<.<<<.<.--......
#=GC RF         c.ucUuAUcaAGAG.gGG.c.gG......
//
```

R2DT uses these annotations:

| Annotation | Purpose |
|---|---|
| `#=GF ID` | Region name (e.g. "SAM") |
| `#=GF AC` | Fallback name if ID is absent (e.g. "RF00162") |
| `#=GC RF` | Reference annotation — uppercase/lowercase positions define match columns; dots mark inserts to remove |
| `#=GC SS_cons` | Consensus structure, may use WUSS notation |

If `#=GC RF` is present, match columns are determined by non-dot RF positions. Otherwise, R2DT computes an IUPAC consensus and strips all-gap columns.

Since only one region is produced, stitching is automatically skipped — no need to pass `--no-stitch`.

### WUSS notation

Rfam and Infernal use [WUSS (Washington University Secondary Structure)](http://eddylab.org/infernal/Userguide-Infernal-1.1.5.pdf) notation in `SS_cons`. R2DT automatically converts WUSS to standard dot-bracket:

| WUSS character | Meaning | Converted to |
|---|---|---|
| `(` `)` | Base pair (depth 1) | `(` `)` |
| `<` `>` | Base pair (depth 2) | `(` `)` |
| `[` `]` | Base pair (depth 3) | `(` `)` |
| `{` `}` | Base pair (depth 4) | `(` `)` |
| `A`…`Z` / `a`…`z` | Pseudoknot (letter pair) | Kept as-is |
| `.` | Unpaired | `.` |
| `,` `-` `_` `:` `~` | WUSS unpaired variants | `.` |

## Example: HCV structural elements

The Hepatitis C virus genome contains numerous conserved RNA secondary structures. An example alignment of 57 HCV sequences with named structures is included at `examples/hcv-alignment.stk`. This file uses the `structureID` + `regionID` format.

### Step 1: Run stockholm command

```bash
r2dt.py stockholm examples/hcv-alignment.stk output/hcv/
```

Output:

```
Processing Stockholm alignment: HCV Nucleotide Structures.stk
  Found 57 sequences
  SS_cons length: 9831
  Found 32 named regions
  Processing region: SLI
    ✓ Generated SVG: SLI_13-27.svg
  Processing region: SLII
    ✓ Generated SVG: SLII_52-125.svg
  ...
  Processing region: X-tail
    Skipping: Unmatched ')' at position 23 (no matching '(')

Summary:
  Processed: 28 regions
  Skipped: 4 regions

Stitching SVG outputs...
✓ Stitched SVG written to: output/hcv/stitched.svg
```

### Step 2: View the stitched output

The stitched SVG shows all valid structural elements in alignment order, with each panel labelled by its parent genomic region (from `regionID`) and structure name:

```{figure} images/hcv-stockholm-stitched.svg
:alt: HCV RNA secondary structures generated from a Stockholm alignment
:width: 100%

All named RNA structures from an HCV alignment, automatically stitched into a single diagram. Each panel shows the consensus secondary structure for one named region, with IUPAC ambiguity codes reflecting conservation across 57 sequences.
```

### Step 3: Use the thumbnail view

The `stockholm` command also produces a **thumbnail** version of the stitched diagram. The thumbnail strips away text, numbering, base-pair lines, and pseudoknots, leaving only the backbone outline of each structure — a compact silhouette that is ideal for embedding in web pages, gallery views, or overview figures.

```{figure} images/hcv-stockholm-thumbnail.svg
:alt: HCV RNA secondary structures — thumbnail view
:width: 100%

Thumbnail view of the same HCV structures. Each panel is reduced to its backbone outline, making it easy to compare the shapes and relative sizes of all structural elements at a glance.
```

The thumbnail is written to `stitched-thumbnail.svg` alongside the full diagram. It is generated automatically — no extra flags are needed.

#### What the thumbnail removes

| Element | Reason |
|---------|--------|
| Nucleotide letters | Clutter at small sizes |
| Numbering lines and tick marks | Not needed for shape overview |
| Base-pair lines | Simplify the silhouette |
| Pseudoknot arcs | Remove visual noise |
| Circle markers | Keep only the backbone path |

The connecting outline between panels and the break glyphs are preserved so the overall genome layout remains clear.

## Command options

```bash
r2dt.py stockholm [OPTIONS] STOCKHOLM_INPUT OUTPUT_FOLDER
```

| Option | Description |
|--------|-------------|
| `--stitch/--no-stitch` | Enable/disable automatic stitching (default: enabled) |
| `--stitch-output PATH` | Custom path for stitched SVG |
| `--monochrome/--color` | Monochrome output (default) or preserve colors |
| `--color-by MODE` | Auto-colour panels: `none` (default), `structure`, or `region` |
| `--color-config PATH` | TSV file mapping structure/region names to SVG colours |
| `--quiet` | Suppress progress output |

## Output files

The `stockholm` command produces:

```
output/
├── processing_summary.txt    # Summary of processed/skipped regions
├── stitched.svg              # Combined diagram (if --stitch)
├── stitched-outline.svg      # Outline-only version (no fill colors)
├── stitched-thumbnail.svg    # Backbone silhouette for galleries/embeds
├── results/
│   ├── svg/                  # Individual SVG diagrams
│   │   ├── SLI_13-27.svg
│   │   ├── SLII_52-125.svg
│   │   └── ...
│   └── fasta/                # Consensus sequences with structures
│       ├── SLI.fasta
│       ├── SLII.fasta
│       └── ...
└── regions/                  # Full working files per region
    ├── SLI/
    │   ├── SLI.fasta
    │   └── r2r/              # R2R intermediate files
    └── ...
```

## RF consensus computation

R2DT computes consensus sequences using IUPAC ambiguity codes, similar to Infernal's RF line:

| Code | Nucleotides | Description |
|------|-------------|-------------|
| A, C, G, U | Single | Conserved nucleotide (>50% frequency) |
| R | A, G | Purine |
| Y | C, U | Pyrimidine |
| S | C, G | Strong |
| W | A, U | Weak |
| K | G, U | Keto |
| M | A, C | Amino |
| N | A, C, G, U | Any nucleotide |
| - | Gap | Majority gaps |

Lowercase letters indicate positions where the consensus nucleotide is present in 50-80% of sequences.

## Validation and skipping

Regions are skipped if they have:

- **Unbalanced brackets** — Mismatched `(` and `)` pairs
- **No base pairs** — Structure is all dots
- **Length mismatch** — Sequence and structure have different lengths

## Coloring panels

By default, stitched diagrams use monochrome styling. R2DT can colour each panel in the stitched output according to its structure name, genomic region, or a custom colour palette. The colour is applied to nucleotide letters, backbone lines, and the connecting outline between panels.

Three coloring modes are available:

### Auto-colour by structure name

Assign a deterministic colour to each panel based on its `structureID` name. Every unique structure gets a different colour:

```bash
r2dt.py stockholm examples/hcv-alignment.stk output/hcv-colored/ --color-by structure
```

### Auto-colour by genomic region

All structures that share a `regionID` get the same colour. This groups related structures visually (for example, all stem-loops within the 5′ UTR share one colour):

```bash
r2dt.py stockholm examples/hcv-alignment.stk output/hcv-colored/ --color-by region
```

### Custom colour palette from a TSV file

For full control, provide a tab-separated file that maps structure or region names to specific SVG colours:

```bash
r2dt.py stockholm examples/hcv-alignment.stk output/hcv-colored/ \
    --color-config examples/color-config.tsv
```

The colour config file has two tab-separated columns — name and colour — with no header row:

```text
# Example colour configuration for HCV genome regions.
# Lines starting with # are comments.
5'UTR	steelblue
core_protein	#e07a5f
E1_protein	#81b29a
E2_protein	#f2cc8f
NS3_protease/helicase	#457b9d
NS5B_RNA-dependent_RNA_polymerase	#2a9d8f
3'UTR	#e76f51
*	gray
```

Colour config rules:

| Feature | Description |
|---------|-------------|
| Columns | Two tab-separated columns: name and SVG colour |
| Comments | Lines starting with `#` are ignored |
| Name matching | Each panel is matched first by its `structureID`, then by its `regionID` |
| Default colour | The special key `*` sets the fallback colour for unmatched panels |
| Colour values | Any SVG-valid colour: named (`steelblue`), hex (`#e07a5f`), `rgb()`, etc. |

An example configuration file is provided at `examples/color-config.tsv`.

### What gets coloured

When panel colours are active, the following elements are styled:

- **Nucleotide letters** — text colour matches the panel accent
- **Backbone lines** — the "gray" backbone strokes use the panel colour
- **Connecting outline** — the path tracing through nucleotide positions uses per-panel colours instead of a single flat gray
- **Thumbnail** — the backbone silhouette in `stitched-thumbnail.svg` also reflects panel colours

:::{note}
`--color-by` and `--color-config` both override `--monochrome`. You do not need to pass `--color` explicitly.
:::
