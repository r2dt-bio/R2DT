# Covariation visualisation

R2DT can display **covariation annotations** from [R-scape](http://eddylab.org/R-scape/) analysis on top of RNA secondary structure diagrams. This is useful for assessing the statistical support for predicted base pairs: positions that covary across an alignment provide strong evidence that the inferred secondary structure is real.

## Motivation

[R-scape](http://eddylab.org/R-scape/) (RNA Structural Covariation Above Phylogenetic Expectation) identifies statistically significant covariation in a multiple sequence alignment and produces annotated Stockholm files with `#=GC cov_SS_cons` and `#=GC cov_h_SS_cons` lines that mark covarying base pairs and helices.

R-scape includes [R2R](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-3) for drawing secondary structure diagrams, but the R2R layout algorithm frequently produces **overlapping nucleotides**, especially for larger or more complex RNAs. This makes the resulting diagrams difficult to read.

R2DT solves this by using [RNApuzzler](https://doi.org/10.1371/journal.pcbi.1006768) to compute **overlap-free layouts**, while preserving the covariation colour annotations from R-scape. Like R-scape/R2R, R2DT draws coloured **rectangles behind base pairs** to show covariation.

## Example: GOLLD RNA (PDB [9LEL](https://www.ebi.ac.uk/pdbe/entry/pdb/9lel))

[GOLLD](https://rfam.org/family/RF01998) (Giant, Ornate, Lake- and Lactobacillus-Derived) RNA is a large bacterial non-coding RNA (~600 nt) related to the [ROOL](https://rfam.org/family/RF03087) RNA motif. Both are unusually large, structurally complex RNAs found in prophages and lactic acid bacteria.

### R-scape / R2R output (with overlaps)

The diagram below is the direct R2R output from R-scape. Green circles mark statistically significant covarying base pairs. Note the **extensive overlapping** of nucleotide labels, making the covariation pattern nearly impossible to interpret:

```{figure} images/9LEL-rscape-r2r.svg
:alt: R-scape R2R diagram of GOLLD RNA showing overlapping nucleotides

R-scape/R2R diagram of GOLLD RNA (9LEL). Green circles indicate covarying positions, but the R2R layout produces significant overlaps that obscure the structure.
```

### R2DT output (overlap-free, with covariation)

R2DT processes the same R-scape Stockholm alignment and produces an **overlap-free layout** with covariation annotations displayed as coloured rectangles behind base pairs:

```{figure} images/9LEL-r2dt-covariation.svg
:alt: R2DT covariation diagram of GOLLD RNA with overlap-free layout

R2DT diagram of GOLLD RNA (9LEL) with covariation annotations. Dark green rectangles highlight base pairs with statistically significant covariation signal. Light green rectangles mark base pairs within covarying helices. The RNApuzzler layout eliminates all overlaps.
```

For comparison, here is the same R2DT layout without covariation colouring:

```{figure} images/9LEL-r2dt-normal.svg
:alt: R2DT diagram of GOLLD RNA without covariation colouring

R2DT diagram of GOLLD RNA (9LEL) without covariation annotations — the default output showing the consensus sequence on an overlap-free layout.
```

## Colour scheme

R2DT uses two colours to distinguish covariation signals, drawn as rectangles that span both nucleotides of a base pair:

| Colour | Rendered | Meaning |
|--------|----------|--------|
| <span style="background-color: rgb(49,163,84); opacity: 0.7; padding: 4px 12px;">Dark green</span> | `rgb(49, 163, 84)` at 70% opacity | **Covariation signal** — the base pair has a statistically significant covariation score (`2` in `cov_SS_cons`) |
| <span style="background-color: rgb(175,240,168); opacity: 0.7; padding: 4px 12px;">Light green</span> | `rgb(175, 240, 168)` at 70% opacity | **Helix context** — the base pair is part of a helix marked in `cov_h_SS_cons` but does not itself have a significant pair-level covariation score |
| No rectangle | — | No covariation signal — unpaired positions and non-covarying pairs have no rectangle |

Each rectangle is oriented along the base-pair line connecting the two nucleotides. Dark green rectangles are compact; light green rectangles are slightly larger to provide visual context.

## How it works

R2DT automatically detects R-scape covariation annotations in Stockholm alignments and generates **two output files**:

1. **`{name}.svg`** — the standard diagram without covariation colouring
2. **`{name}.covariation.svg`** — the same layout with covariation rectangles

### Detection

R2DT looks for these `#=GC` annotation lines in the Stockholm file:

| Annotation | Source | Meaning |
|---|---|---|
| `#=GC cov_SS_cons` | R-scape | Per-column pair-level covariation. `2` = significant, `.` = none |
| `#=GC cov_h_SS_cons` | R-scape | Per-column helix-level covariation. `3` = covarying helix, `.` = none |

If either annotation is present, R2DT will automatically generate the covariation SVG.

### Pipeline

The covariation pipeline runs inside the `stockholm` command:

1. **Parse** the R-scape Stockholm file and extract `cov_SS_cons` and `cov_h_SS_cons`
2. **Map** annotation columns to gap-free consensus positions
3. **Classify** each position as `dark_green` (pair or helix covariation), `light_green` (helix context only), or `none`
4. **Symmetrise** — propagate the highest-priority colour to both sides of each base pair
5. **Inject** coloured `<rect>` elements into the SVG, oriented along each base-pair line

## Usage

### From an R-scape Stockholm file

If you have an R-scape `.R2R.sto` file (produced by `R-scape --cacofold`):

```bash
r2dt.py stockholm alignment.cacofold.R2R.sto output/
```

R2DT will automatically:
- Detect the `cov_SS_cons` / `cov_h_SS_cons` annotations
- Generate the overlap-free layout using RNApuzzler (with R2R fallback)
- Produce both `{name}.svg` and `{name}.covariation.svg`

### Running R-scape first

To generate the R-scape Stockholm file from a seed alignment:

```bash
# Run R-scape with CACOfold to predict structure with covariation
R-scape --cacofold alignment.sto

# This produces alignment.cacofold.R2R.sto with covariation annotations
# Then run R2DT
r2dt.py stockholm alignment.cacofold.R2R.sto output/
```

## Related

- [Annotations and data layers](annotations.md) — alignment accuracy and custom annotations
- [Processing Stockholm alignments](stockholm-alignments.md) — multi-region Stockholm input format
- [R-scape documentation](http://eddylab.org/R-scape/) — covariation analysis software
- [ROOL RNA (RF03087)](https://rfam.org/family/RF03087) — Rumen-Originating, Ornate, Large RNA
- [GOLLD RNA (RF01998)](https://rfam.org/family/RF01998) — Giant, Ornate, Lake- and Lactobacillus-Derived RNA
