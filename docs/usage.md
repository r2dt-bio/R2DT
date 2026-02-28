# Command line reference

Specify the input file in [FASTA format](https://en.wikipedia.org/wiki/FASTA_format) containing one or more RNA sequences as well as the path where the output files will be created (the folder will be created if it does not exist).

```bash
r2dt.py draw <input.fasta> <output_folder>
```

For example:

```bash
r2dt.py draw examples/examples.fasta temp/examples
```

R2DT will automatically select the best matching template and visualise the secondary structures.

## Output files

`r2dt.py draw` produces a folder called `results` with the following subfolders:

- `svg`: RNA secondary structure diagrams in SVG format
- `fasta`: input sequences and their secondary structure in dot-bracket notation
- `tsv`: a file `metadata.tsv` listing sequence ids, matching templates, and template sources
- `thumbnail`: secondary structure diagrams displayed as outlines in SVG format
- `json`: RNA secondary structure and its layout described using [RNA 2D JSON Schema](https://github.com/LDWLab/RNA2D-data-schema/)

## Manually selecting template type

If the RNA type of the input sequences is known in advance, it is possible to bypass the classification steps and achieve faster performance.

* CRW templates (5S and SSU rRNA)
    ```bash
    r2dt.py crw draw examples/crw-examples.fasta temp/crw-examples
    ```

* RiboVision LSU and SSU rRNA templates
    ```bash
    r2dt.py ribovision draw_lsu examples/lsu-examples.fasta temp/lsu-examples
    r2dt.py ribovision draw_ssu examples/ribovision-ssu-examples.fasta temp/ssu-examples
    ```

* Rfam families
    ```bash
    r2dt.py rfam draw RF00162 examples/RF00162.example.fasta temp/rfam-example
    ```

    By default, R2DT uses the template with minimum number of overlaps (R-scape or RNArtist). However, it is possible to specify the template type manually.

    Specify `--rnartist` to use the RNArtist template:

    ```bash
    r2dt.py rfam draw RF00162 examples/RF00162.example.fasta temp/rfam-example --rnartist
    ```

    Specify `--rscape` to use the R-scape template:

    ```bash
    r2dt.py rfam draw RF00162 examples/RF00162.example.fasta temp/rfam-example --rnartist
    ```

    These options could be useful to avoid overlaps in the diagrams, for example:

    ```{figure} images/rnartist-example.png
    :alt: R2R vs RNArtist

    A comparison of the cobalamin riboswitch ([RF00174](https://rfam.org/family/RF00174)) visualised using R2R (left) and RNArtist (right). The RNArtist template is more legible and avoids overlaps.
    ```

* RNAse P
    ```bash
    r2dt.py rnasep draw examples/rnasep.fasta temp/rnasep-example
    ```

* tRNAs (using GtRNAdb templates)
    ```bash
    # for tRNAs, provide domain and isotype (if known), or use tRNAScan-SE to classify
    r2dt.py gtrnadb draw examples/gtrnadb.E_Thr.fasta temp/gtrnadb
    r2dt.py gtrnadb draw examples/gtrnadb.E_Thr.fasta temp/gtrnadb --domain E --isotype Thr
    ```

## Manually selecting specific template

It is possible to select a specific template and skip the classification step altogether.

1. Get a list of all available templates and copy the template id:
    ```bash
    r2dt.py list-models
    ```

In addition, all models are listed in the file [models.json](https://github.com/r2dt-bio/R2DT/blob/main/data/models.json).

1. Specify the template (for example, `RNAseP_a_P_furiosus_JB`):
    ```bash
    r2dt.py draw --force_template <template_id> <input_fasta> <output_folder>
    ```

    For example:

    ```bash
    r2dt.py draw --force_template RNAseP_a_P_furiosus_JB examples/force/URS0001BC2932_272844.fasta temp/example
    ```

## Constraint-based folding for insertions

If a structure contains insertions relative to the R2DT template files, they are shown as large unstructured loops. Insertions larger than 100 nucleotides are replaced with a placeholder element (`XXXX` characters) to keep the rest of the diagram legible.

The `--constraint` flag allows de-novo prediction of the insertion structure using the RNAfold algorithm. There are currently three constraint folding modes available. R2DT will automatically predict which folding mode is best for a given molecule, but the mode can also be manually overridden using the --fold_type parameter. There are three options for fold_type.

* Let R2DT pick a fold_type
    ```bash
    r2dt.py draw --constraint <input_fasta> <output_folder>
    ```
* Fold insertions (along with adjacent unpaired nucleotides) one at a time. Recommended for large RNAs.
    ```bash
    r2dt.py draw --constraint --fold_type insertions_only <input_fasta> <output_folder>
    ```
* Run entire molecule through RNAfold at once. Base pairs predicted from the template are used as constraints for prediction.
    ```bash
    r2dt.py draw --constraint --fold_type full_molecule <input_fasta> <output_folder>
    ```
* Run entire molecule through RNAfold at once. Both conserved single-stranded regions and base pairs predicted from the template are used as constraints for prediction.
    ```bash
    r2dt.py draw --constraint --fold_type all_constraints_enforced <input_fasta> <output_folder>
    ```
* Prevent certain nucleotides from base pairing. This will only work for base pairs that are de-novo predicted.
The exclusion file should contain a string the same length as the input sequence composed of '.'s and 'x's. Positions with '.'s are allowed to base pair,
positions with 'x's are not.
Example string: `xxxx..............xx..............x............xx`
    ```bash
    r2dt.py draw --constraint --exclusion <exclusion_file> <input_fasta> <output_folder>
    ```

## Template-free visualisation

It is possible to visualise a sequence and its secondary structure using layouts generated by [R2R](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-3), [RNArtist](https://github.com/fjossinet/RNArtist), or [RNApuzzler](https://doi.org/10.1371/journal.pcbi.1006768) (via ViennaRNA). This functionality is useful as a starting point when [generating new templates](./templates.md) or in cases when the R2DT template library does not yet have a template for a certain RNA.

Example input found in `examples/template-free.fasta` (pseudoknots can be specified in the Aa, Bb notation):

```
>3SKZ_B
GGCCUUAUACAGGGUAGCAUAAUGGGCUACUGACCCCGCCUUCAAACCUAUUUGGAGACUAUAAGGUC
.((((((((A..((((((.....BB))))))(.....a)(((((((bb..)))))))..)))))))).
```

R2DT automatically detects if the secondary structure is present in the input FASTA file and will use it to generate a diagram:

```bash
r2dt.py draw examples/template-free.fasta temp/template-free-example
```

Alternatively, one can explicit specify the `templatefree` mode:

```bash
r2dt.py templatefree examples/template-free.fasta temp/template-free-example
```

Output diagram:

![Template free visualisation example](./images/template-free-example.png)

The output files are organised in the standard way and include SVG and RNA 2D JSON Schema files. The JSON file can be uploaded to [RNAcanvas](https://rnacanvas.app/) for further manual editing, if necessary.

## Processing Stockholm alignments

R2DT can process Stockholm-format multiple sequence alignments with named secondary structure regions. This is useful for viral genome alignments where structural elements have been manually annotated.

```bash
r2dt.py stockholm examples/hcv-alignment.stk output/hcv/
```

The command extracts named regions from `#=GC structureID` (or the legacy `#=GC knownSS_names`), computes RF-style consensus sequences with IUPAC codes, generates template-free visualizations, and optionally stitches them into a combined diagram. When `#=GC regionID` is present, each structure is labelled with its parent genomic region.

The stitched output includes three variants:
- **stitched.svg** — full diagram with nucleotide letters and annotations
- **stitched-outline.svg** — outline-only version (no fill colors)
- **stitched-thumbnail.svg** — compact backbone silhouette with no text, numbering, or base-pair lines (ideal for galleries and web embedding)

Use `--auto-repair` to attempt fixing unbalanced bracket structures in the annotations:

```bash
r2dt.py stockholm --auto-repair alignment.stk output/
```

### Coloring panels

Panels in the stitched output can be coloured by structure name, genomic region, or a custom palette:

```bash
# Auto-colour each panel by its structure name
r2dt.py stockholm examples/hcv-alignment.stk output/hcv/ --color-by structure

# Auto-colour by genomic region (structures sharing a region get the same colour)
r2dt.py stockholm examples/hcv-alignment.stk output/hcv/ --color-by region

# Use a custom TSV colour palette
r2dt.py stockholm examples/hcv-alignment.stk output/hcv/ \
    --color-config examples/color-config.tsv
```

The colour config file (`examples/color-config.tsv`) maps names to SVG colours:

```text
5'UTR	steelblue
core_protein	#e07a5f
NS3_protease/helicase	#457b9d
*	gray
```

The special key `*` sets the default colour for panels not listed explicitly. See [](stockholm-alignments.md) for full details.

For detailed documentation, see [](stockholm-alignments.md).

### Layout engines

Sometimes the output diagrams may contain overlaps:

![Bridge RNA visualised with R2R](./images/template-free-example-r2r.svg){w=500px}

Trying a different layout engine could help minimise the overlaps. R2DT supports three layout engines:

| Flag | Engine | Description |
|------|--------|-------------|
| `--rscape` | R2R (default) | Fast consensus layout |
| `--rnartist` | RNArtist | Java-based layout, good for complex structures |
| `--rnapuzzler` | RNApuzzler | Overlap-free layout via ViennaRNA |

```bash
# Use RNArtist
r2dt.py templatefree examples/bridge-rna.fasta temp/bridge-rna-rnartist --rnartist

# Use RNApuzzler for an overlap-free layout
r2dt.py templatefree examples/bridge-rna.fasta temp/bridge-rna-puzzler --rnapuzzler
```

Output diagram (RNArtist):

![Bridge RNA visualised with RNArtist](./images/template-free-example-rnartist.svg)

Use `--auto` to run all three engines and automatically pick the layout with the fewest overlaps:

```bash
r2dt.py templatefree examples/bridge-rna.fasta temp/bridge-rna-auto --auto
```

## Skipping ribovore filters

In some cases R2DT may not generate a diagram for a sequence because [ribovore](https://github.com/ncbi/ribovore) detects one or more [unexpected features](https://github.com/ncbi/ribovore/blob/main/documentation/ribotyper.md#unexpectedfeatures), such as having hits on both strands or having too many hits in the same sequence. You can use `--skip_ribovore_filters` to ignore these warnings and attempt to generate a secondary structure diagram anyway.

For example, the following command will produce no results because the sequence is close to a palindrome:

```bash
r2dt.py draw examples/ribovore-qc-example.fasta temp/examples
```

However, the following command generates a [valid diagram](https://github.com/r2dt-bio/R2DT/blob/main/tests/examples/skip-ribovore-filters/URS0000001EB3-RF00661.colored.svg):

```bash
r2dt.py draw --skip_ribovore_filters examples/ribovore-filters.fasta temp/examples
```

Please note that this option should be used with caution as sequences with unexpected features often result in poor diagrams.

## Stitching multiple diagrams

The `stitch` command combines multiple R2DT SVG diagrams into a single horizontal layout. This is useful for displaying multiple RNA structures from a viral genome, comparing related structures, or creating publication-ready figures.

### Basic usage

```bash
r2dt.py stitch diagram1.svg diagram2.svg diagram3.svg -o combined.svg
```

### Sorting by genomic coordinates

When SVG filenames include genomic coordinates (e.g., `RF00507_13469-13546.svg`), the `--sort` flag arranges panels in genomic order:

```bash
r2dt.py stitch output/rfam/*.svg -o stitched.svg --sort
```

### Adding captions

Add labels above each panel using `--captions` (specify once per panel):

```bash
r2dt.py stitch \
    5utr.svg fse.svg 3utr.svg \
    -o combined.svg \
    --captions "5′ UTR" --captions "FSE" --captions "3′ UTR"
```

### Available options

| Option | Default | Description |
|--------|---------|-------------|
| `-o, --output PATH` | Required | Output SVG file path |
| `--gap FLOAT` | 100 | Horizontal gap between panels (pixels) |
| `--sort` | False | Sort panels by genomic coordinates |
| `--glyph TYPE` | break | Join glyph: `none`, `bead`, `bar`, `break` |
| `--monochrome/--color` | monochrome | Monochrome (default) or preserve original colors |
| `--color-by MODE` | none | Auto-colour panels: `none`, `structure`, or `region` |
| `--color-config PATH` | None | TSV file mapping names to SVG colours |
| `--captions TEXT` | None | Caption for each panel (repeat for each) |
| `--caption-font-size FLOAT` | auto | Font size for captions |
| `--keep-intermediate-labels` | False | Show all 5′/3′ labels |
| `--no-gap-labels` | False | Hide nucleotide distance labels |
| `--gap-label-font-size FLOAT` | auto | Font size for gap labels |
| `--no-outline` | False | Disable connecting outline stroke |
| `--outline-color COLOR` | #cccccc | Outline stroke color |
| `--outline-width FLOAT` | 3.0 | Outline stroke width |
| `--outline-opacity FLOAT` | 0.6 | Outline opacity (0-1) |
| `--anchor-label-font-size FLOAT` | auto | Font size for 5′/3′ labels |
| `--normalize-font-size` | False | Scale panels to match nucleotide font size of first panel |
| `-q, --quiet` | False | Suppress output messages |

### Join glyphs

The `--glyph` option controls how panels are visually connected:

- `none` - No connecting element
- `bead` - Small circle at each join point
- `bar` - Vertical line connecting 3′ to 5′
- `break` - Double diagonal lines indicating a sequence break (default)

### Example with all options

```bash
r2dt.py stitch \
    panel1.svg panel2.svg panel3.svg \
    -o publication-figure.svg \
    --sort \
    --gap 150 \
    --glyph break \
    --captions "Region A" --captions "Region B" --captions "Region C" \
    --monochrome \
    --outline-color "#999999" \
    --outline-width 2.0
```

:::{note}
Font sizes for captions, gap labels, and anchor labels are automatically detected from the nucleotide font size in the input SVGs. Override with explicit values if needed.
:::

For a complete workflow example using `stitch` with viral genomes, see [Annotating viral genomes](viral-genomes.md).

## Other useful commands

* Print R2DT version
    ```bash
    r2dt.py version
    ```

* Classify example sequences using Ribotyper
    ```bash
    perl /rna/ribovore/ribotyper.pl -i data/cms/crw/modelinfo.txt -f examples/pdb.fasta temp/ribotyper-test
    ```

* Generate covariance models and modelinfo files
    ```bash
    python3 utils/generate_cm_library.py
    r2dt.py generatemodelinfo <path to covariance models>
    ```

* Precompute template library locally (may take several hours):
    ```bash
    r2dt.py setup
    ```

* Run R2DT with Singularity
    ```bash
    singularity exec --bind <path_to_cms>:/rna/r2dt/data/cms r2dt r2dt.py draw sequence.fasta output
    ```

* Convert a SVG diagram to a JSON file containing the paths per nucleotide and an ordinal numbering. Note that this *assumes* that the input pdb id is formatted like: `<PDB>_<Entity>_<chain>`, ie `1S72_1_0`.
    ```bash
    svg2json.py <pdb-id> diagram.svg <pdb-id>.json
    ```
