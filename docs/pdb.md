# Visualising RNA from PDB structures

The `pdb` command generates R2DT diagrams directly from 3D structures. It can work with either PDB IDs (downloading from RCSB) or local structure files.

## Basic usage

```bash
# Using a PDB ID (downloads from RCSB)
r2dt.py pdb <PDB-ID> <output_folder>

# Using a local file (PDB or mmCIF format)
r2dt.py pdb <structure_file> <output_folder>
```

For example, to visualise the tRNA-Phe from yeast:

```bash
r2dt.py pdb 1EHZ output/
```

## Using local structure files

R2DT supports local PDB and mmCIF files, including gzip-compressed files:

```bash
# Local PDB file
r2dt.py pdb structure.pdb output/

# Local mmCIF file
r2dt.py pdb structure.cif output/

# Gzip-compressed files (auto-detected)
r2dt.py pdb structure.pdb.gz output/
r2dt.py pdb structure.cif.gz output/
```

The file format is automatically detected from the file content, not the extension. This means renamed files will still work correctly.

:::{tip}
Local file support is useful when working with structures not yet deposited in the PDB, modified structures, or when you have already downloaded structure files.
:::

## Choosing a base pair extractor

R2DT supports two tools for extracting base pairs from 3D coordinates:

- **FR3D** (default) - Works with both PDB and mmCIF files, better for complex structures
- **RNAView** - Fast, but requires PDB format

Use the `--basepairs` option to override the default:

```bash
# Force RNAView (requires PDB format)
r2dt.py pdb 1S72 output/ --basepairs rnaview

# Force FR3D explicitly
r2dt.py pdb 9FN3 output/ --basepairs fr3d
```

:::{tip}
RNAView only works with PDB format files. Use `--basepairs rnaview` only when you specifically need RNAView output.
:::

## Specifying structure format

Control which file format to download from RCSB:

```bash
# Auto-select based on basepairs tool (default)
r2dt.py pdb 1S72 output/ --format auto

# Prefer PDB format
r2dt.py pdb 1S72 output/ --format pdb

# Prefer mmCIF format
r2dt.py pdb 1S72 output/ --format cif
```

:::{note}
RNAView only works with PDB format files. If you specify `--basepairs rnaview` with `--format cif`, R2DT will report an error.
:::

## Extracting a specific chain

By default, R2DT extracts the first RNA chain found in the structure. Use `--chain` to specify a particular chain:

```bash
# Extract chain A from tRNA structure
r2dt.py pdb 1EHZ output/ --chain A
```

## Output files

The `pdb` command creates the following in the output folder:

- `downloads/` - Downloaded structure files from RCSB
- `<PDB-ID>.fasta` - Extracted sequence and secondary structure
- `results/` - Standard R2DT output (SVG, JSON, thumbnails, etc.)

## Examples

```bash
# Visualise tRNA from PDB
r2dt.py pdb 1EHZ output/

# Visualise bridge RNA (mmCIF only, requires FR3D)
r2dt.py pdb 9FN3 output/ --basepairs fr3d

# Extract specific chain with FR3D
r2dt.py pdb 1EHZ output/ --basepairs fr3d --chain A
```

## Missing nucleotides in PDB structures

Many PDB structures have nucleotides that were not resolved experimentally (e.g. flexible loops or disordered regions). R2DT automatically detects these missing nucleotides and includes them in the 2D diagram, where they are shown in grey to distinguish them from resolved nucleotides shown in black.

For example, [9MME](https://www.rcsb.org/structure/9MME) has ~580 nucleotides in its deposited sequence but only ~521 are resolved in the 3D coordinates. R2DT generates a diagram based on the full deposited sequence, with the ~59 unresolved nucleotides greyed out:

```bash
r2dt.py pdb 9MME output/
```

```{figure} images/9mme-missing-nucleotides.svg
:alt: 9MME with missing nucleotides greyed out

2D diagram of PDB 9MME. Resolved nucleotides are shown in black, while unresolved nucleotides (missing from the 3D coordinates) are greyed out.
```

R2DT extracts the full deposited sequence from SEQRES records (PDB format) or `_entity_poly.pdbx_seq_one_letter_code_can` (mmCIF format) and compares it with the residues present in the 3D coordinates. The secondary structure is computed from the resolved nucleotides and then remapped onto the full sequence, with unresolved positions shown as unpaired.

:::{note}
If R2DT cannot determine the full sequence (e.g. if SEQRES records are missing), it falls back to using only the resolved nucleotides, and no greying is applied.
:::

## Animating RNA secondary structures

R2DT can be used to generate animations of RNA secondary structures. The `animate.py` script takes two R2DT SVG files and generates an animated SVG file that transitions between the two structures.

```bash
python3 utils/animate.py \
    examples/animate/PZ39_solution.svg \
    examples/animate/PZ39_Dfold_3.svg \
    animated.svg
```

:::{note}
Note that the two input SVG files must have the same number of nucleotides.
:::

### Animating from PDB files

R2DT can be used to generate secondary structure diagrams directly from PDB files. Additionally R2DT can generate animated SVG file of the transition between two structures between two 3D structures from PDB files.

Generation of a diagram from PDB file:

```bash
python3 utils/rnaview.py \
    examples/PZ1_Bujnicki_1.pdb
```
`PZ1_Bujnicki_1.colored.svg` will be generated.

:::{note}
Note - in cases that there are multiple interactions detected for one nucleotide, all of them will be omitted in the diagram.
:::

Animated transition between two structures from PDB:

```bash
python3 utils/animate3d.py \
    examples/PZ1_solution_0.pdb \
    examples/PZ1_Bujnicki_1.pdb
```

Performing animation on a reference PDB file, and a set of PDBs. User needs to specify the reference PDB file, and a directory with query PDB files:

```bash
python3 utils/animate3d.py -b \
    examples/PZ1_solution_0.pdb \
    examples/animate_bulk
```

:::{note}
Note - structures in the PDB files should be of the same length [nt].
:::
