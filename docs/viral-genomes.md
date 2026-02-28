# Annotating viral genomes

R2DT can automatically scan viral genomes to identify and visualise non-coding RNA elements using the [Rfam](https://rfam.org) database. This tutorial demonstrates the complete workflow from a viral genome FASTA file to a diagram showing all RNA structures.

## Overview

Many viruses contain conserved RNA secondary structures that play critical roles in their life cycles:

- **5′ and 3′ UTRs** - Untranslated regions with regulatory structures
- **Frameshift elements** - Programmed ribosomal frameshifting signals
- **Internal ribosome entry sites (IRES)** - Cap-independent translation initiation
- **Packaging signals** - RNA structures involved in genome packaging

R2DT uses Infernal's `cmscan` tool to search viral genomes against the Rfam covariance model library, then generates secondary structure diagrams for each identified RNA family.

## Quick start

```bash
r2dt.py viral-annotate genome.fasta output/
```

This command:
1. Scans the genome against Rfam using GA (gathering) thresholds
2. Filters and ranks overlapping hits
3. Generates 2D diagrams for each RNA family found
4. Optionally stitches all diagrams into a single combined view

## Example: SARS-CoV-2 coronavirus

The SARS-CoV-2 genome (~30,000 nt) contains several well-characterised RNA structures that are important for viral replication.

### Step 1: Prepare input

This example uses the SARS-CoV-2 genome with accession [OX309346.1](https://www.ebi.ac.uk/ena/browser/view/OX309346.1), which is included in the R2DT repository at `examples/viral/coronavirus.fasta`.

### Step 2: Run viral-annotate

```bash
r2dt.py viral-annotate examples/viral/coronavirus.fasta output/coronavirus/
```

Output:

```
# R2DT :: visualise RNA secondary structure using templates
# Version 2.2 (2026)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Step 1: Calculating genome size
  Genome size: 29,903 nt
  Database size (Z): 0.059806 Mb

Step 2: Running cmscan with Rfam GA thresholds
  Running: cmscan... (this may take several minutes)

Step 3: Parsing cmscan results
  Found 3 RNA family hits:
    RF03120 (Sarbecovirus-5UTR): 26-299 + score=309.9
    RF00507 (Corona_FSE): 13,469-13,546 + score=77.6
    RF03125 (Sarbecovirus-3UTR): 29,536-29,870 + score=406.9

Step 4: Generating 2D diagrams for each hit
  ✓ RF03120 (Sarbecovirus-5UTR)
  ✓ RF00507 (Corona_FSE)
  ✓ RF03125 (Sarbecovirus-3UTR)

Summary
  Genome: OX309346.1 (29,903 nt)
  RNA families found: 3
  Diagrams generated: 3
```

### Step 3: Examine the output

The `output/coronavirus/` directory contains:

```
output/coronavirus/
├── cmscan.tblout          # Tabular cmscan results
├── cmscan.out             # Full cmscan output
└── rfam/                  # Individual RNA diagrams
    ├── RF03120_26-299.colored.svg
    ├── RF00507_13469-13546.colored.svg
    └── RF03125_29536-29870.colored.svg
```

Each SVG file shows the secondary structure of one RNA element:

- **RF03120** (Sarbecovirus-5UTR): The 5′ untranslated region containing stem-loops SL1-SL5

```{figure} images/RF03120_26-299.svg
:alt: SARS-CoV-2 5′ UTR secondary structure
:width: 100%

5′ UTR of SARS-CoV-2 (nucleotides 26-299)
```

- **RF00507** (Corona_FSE): The programmed ribosomal frameshift element between ORF1a and ORF1b

```{figure} images/RF00507_13469-13546.svg
:alt: Coronavirus FSE (frameshift element)
:width: 50%

FSE - Frameshift element (nucleotides 13,469-13,546)
```

- **RF03125** (Sarbecovirus-3UTR): The 3′ untranslated region containing the s2m element

```{figure} images/RF03125_29536-29870.svg
:alt: SARS-CoV-2 3′ UTR secondary structure
:width: 100%

3′ UTR of SARS-CoV-2 (nucleotides 29,536-29,870)
```

### Step 4: Create a stitched diagram

To combine all diagrams into a single panoramic view:

```bash
r2dt.py stitch \
    output/coronavirus/rfam/*.colored.svg \
    -o coronavirus-stitched.svg \
    --sort \
    --captions "5′ UTR" --captions "FSE" --captions "3′ UTR"
```

The `--sort` flag arranges panels by genomic coordinates (extracted from filenames), and `--captions` adds labels above each panel.

The resulting stitched diagram shows all three RNA structures in genomic order:

```{raw} html
<script src="https://cdn.jsdelivr.net/npm/svg-pan-zoom@3.6.1/dist/svg-pan-zoom.min.js"></script>

<div style="border: 1px solid #ddd; border-radius: 4px; margin: 1em 0; position: relative;">
  <div style="background: #f5f5f5; padding: 8px 12px; border-bottom: 1px solid #ddd; font-size: 0.85em; color: #666;">
    🔍 Use mouse wheel to zoom, drag to pan. <button onclick="panZoomInstance.reset()" style="margin-left: 10px; cursor: pointer;">Reset view</button>
  </div>
  <div id="stitched-container" style="width: 100%; height: 400px; overflow: hidden;">
    <object id="stitched-svg" type="image/svg+xml" data="_images/coronavirus-stitched.svg" style="width: 100%; height: 100%;"></object>
  </div>
</div>

<script>
  var panZoomInstance;
  document.getElementById('stitched-svg').addEventListener('load', function() {
    var svgDoc = this.contentDocument;
    if (svgDoc) {
      var svgElement = svgDoc.querySelector('svg');
      if (svgElement) {
        panZoomInstance = svgPanZoom(svgElement, {
          zoomEnabled: true,
          controlIconsEnabled: true,
          fit: true,
          center: true,
          minZoom: 0.1,
          maxZoom: 10
        });
      }
    }
  });
</script>

<p style="text-align: center; font-style: italic; color: #666; margin-top: 0.5em;">
  Combined view of SARS-CoV-2 RNA structures: 5′ UTR, FSE (frameshift element), and 3′ UTR.
</p>
```

```{figure} images/coronavirus-stitched.svg
:alt: SARS-CoV-2 stitched RNA structures
:class: hidden

This figure is hidden but ensures Sphinx copies the image to _images/
```

```{raw} html
<style>
.hidden {
  display: none !important;
}
</style>
```

## Two approaches to viral RNA annotation

R2DT supports two complementary approaches for annotating viral genomes:

| Approach | Command | Best for |
|----------|---------|----------|
| **FASTA + Rfam** | `viral-annotate` | Automatic discovery — finds known RNA families in any genome |
| **Stockholm alignment** | `stockholm` | Expert curation — uses manually annotated structures and regions |

The FASTA approach (shown above with SARS-CoV-2) automatically scans the genome against Rfam and is fully automated. However, it can only find structures that already have Rfam models.

The Stockholm approach uses a curated multiple sequence alignment where structures have been manually annotated with `#=GC structureID` and `#=GC regionID` lines. This can capture structures that are not in Rfam, and groups them by genomic region (e.g. 5′UTR, NS5B).

### Example: HCV using Stockholm alignment

The HCV genome contains over 40 annotated RNA structures across the 5′UTR, coding regions, and 3′UTR. An alignment of 57 HCV sequences with named structures is included at `examples/hcv-alignment.stk`.

```bash
r2dt.py stockholm examples/hcv-alignment.stk output/hcv-stockholm/
```

Each structure is labelled with its parent genomic region, producing captions like _"SLI (5′UTR)"_ or _"5BSL3.1 (NS5B)"_ in the stitched output:

```{raw} html
<div style="border: 1px solid #ddd; border-radius: 4px; margin: 1em 0; position: relative;">
  <div style="background: #f5f5f5; padding: 8px 12px; border-bottom: 1px solid #ddd; font-size: 0.85em; color: #666;">
    🔍 Use mouse wheel to zoom, drag to pan. <button id="hcv-stk-reset" style="margin-left: 10px; cursor: pointer;">Reset view</button>
  </div>
  <div id="hcv-stk-container" style="width: 100%; height: 400px; overflow: hidden;">
    <object id="hcv-stk-svg" type="image/svg+xml" data="_images/hcv-stockholm-stitched.svg?v=2" style="width: 100%; height: 100%;"></object>
  </div>
</div>
```

```{figure} images/hcv-stockholm-stitched.svg
:alt: HCV RNA structures from Stockholm alignment with structureID and regionID annotations
:class: hidden

HCV RNA structures generated from a Stockholm alignment using `#=GC structureID` and `#=GC regionID` annotations. Structures are labelled with their parent genomic region.
```

Compare this with the [Rfam-based HCV diagram](#hepatitis-c-virus-hcv) in the gallery below — the Stockholm approach captures significantly more structures because it supports manually curated annotations, including families that have not yet been added to Rfam.

## Command reference

### viral-annotate

```bash
r2dt.py viral-annotate <genome.fasta> <output_folder> [OPTIONS]
```

**Arguments:**
- `genome.fasta` - Input viral genome in FASTA format (one sequence)
- `output_folder` - Directory for output files

**Options:**

| Option | Default | Description |
|--------|---------|-------------|
| `--stitch-output PATH` | None | Path for stitched SVG output |
| `--cm-library PATH` | data/rfam/cms/all.cm | Path to Rfam CM library |
| `--clanin PATH` | None | Path to Rfam.clanin for clan competition |
| `--cpu N` | 4 | Number of CPUs for cmscan |
| `--evalue`, `-E` | None | E-value threshold (default: use Rfam GA thresholds) |
| `--monochrome/--color` | monochrome | Monochrome (default) or preserve original colors |
| `--quiet` | False | Suppress progress messages |

### stitch

See the [stitch command reference](usage.md#stitching-multiple-diagrams) for all available options.

## Understanding the results

### cmscan.tblout format

The tabular output from cmscan contains detailed hit information:

```
#idx  target name         accession  query name  ...  seq from  seq to  strand  ...  score  E-value
1     Sarbecovirus-5UTR   RF03120    OX309346.1  ...  26        299     +       ...  309.9  3.2e-76
2     Corona_FSE          RF00507    OX309346.1  ...  13469     13546   +       ...  77.6   1.1e-15
3     Sarbecovirus-3UTR   RF03125    OX309346.1  ...  29536     29870   +       ...  406.9  1.5e-98
```

### Genomic coordinates in filenames

SVG filenames encode the genomic coordinates:

```
RF00507_13469-13546.colored.svg
   │      │     │
   │      │     └── End position (1-based)
   │      └── Start position (1-based)
   └── Rfam accession
```

The `stitch` command uses these coordinates to:
- Order panels by genomic position (`--sort`)
- Calculate nucleotide distances between panels
- Display gap labels showing the distance in nucleotides

## Advanced usage

### Using a custom CM library

To scan against a subset of Rfam families or custom models:

```bash
r2dt.py viral-annotate genome.fasta output/ --cm-library my-models.cm
```

The CM library will be automatically indexed with `cmpress` if needed.

### Customising the stitched output

The `stitch` command provides many options for customising the combined diagram:

```bash
r2dt.py stitch \
    output/rfam/*.svg \
    -o stitched.svg \
    --sort \
    --gap 150 \
    --glyph break \
    --captions "5′ UTR" --captions "FSE" --captions "3′ UTR" \
    --color \
    --no-outline
```

By default, stitched diagrams are monochrome (black and white), which is suitable for publications. Use `--color` to preserve the original nucleotide coloring:

```{raw} html
<div style="border: 1px solid #ddd; border-radius: 4px; margin: 1em 0; position: relative;">
  <div style="background: #f5f5f5; padding: 8px 12px; border-bottom: 1px solid #ddd; font-size: 0.85em; color: #666;">
    🔍 Use mouse wheel to zoom, drag to pan. <button onclick="panZoomColor.reset()" style="margin-left: 10px; cursor: pointer;">Reset view</button>
  </div>
  <div id="stitched-color-container" style="width: 100%; height: 400px; overflow: hidden;">
    <object id="stitched-color-svg" type="image/svg+xml" data="_images/coronavirus-stitched-color.svg" style="width: 100%; height: 100%;"></object>
  </div>
</div>

<script>
  var panZoomColor;
  document.getElementById('stitched-color-svg').addEventListener('load', function() {
    var svgDoc = this.contentDocument;
    if (svgDoc) {
      var svgElement = svgDoc.querySelector('svg');
      if (svgElement) {
        panZoomColor = svgPanZoom(svgElement, {
          zoomEnabled: true,
          controlIconsEnabled: true,
          fit: true,
          center: true,
          minZoom: 0.1,
          maxZoom: 10
        });
      }
    }
  });
</script>

<p style="text-align: center; font-style: italic; color: #666; margin-top: 0.5em;">
  Colored version of the stitched diagram using <code>--color</code> flag.
</p>
```

```{figure} images/coronavirus-stitched-color.svg
:alt: SARS-CoV-2 stitched RNA structures (colored)
:class: hidden

This figure is hidden but ensures Sphinx copies the image to _images/
```

See [Stitching multiple diagrams](usage.md#stitching-multiple-diagrams) for the full option reference.

### Integrating with other R2DT workflows

The `viral-annotate` command is a convenience wrapper. For more control, you can run the steps manually:

1. **Run cmscan separately:**
   ```bash
   cmscan -Z 0.06 --cut_ga --rfam --nohmmonly \
       --tblout hits.tblout --fmt 2 --cpu 4 \
       data/rfam/cms/all.cm genome.fasta > cmscan.out
   ```

2. **Extract hit regions and visualise:**
   ```bash
   # For each hit, extract the sequence and run:
   r2dt.py rfam draw RF00507 hit.fasta output/
   ```

3. **Stitch the results:**
   ```bash
   r2dt.py stitch output/results/svg/*.colored.svg -o combined.svg --sort
   ```

## Supported virus families

R2DT can identify RNA structures in many virus families. The Rfam database contains models for:

- **Coronaviridae** - 5′/3′ UTRs, frameshift elements, s2m
- **Flaviviridae** - 5′/3′ UTRs, pseudoknots
- **Picornaviridae** - IRES elements, cre
- **Retroviridae** - TAR, RRE, packaging signals
- **And many more...**

Search [Rfam](https://rfam.org/search) to find RNA families associated with your virus of interest.

## Troubleshooting

### No hits found

If cmscan finds no hits:
- The genome may not contain characterised RNA families in Rfam
- Try lowering the threshold with `--evalue 1` to see marginal hits

### cmscan is slow

For large genomes or many sequences:
- Increase `--cpu` to use more processors
- Consider searching with a smaller CM library targeting expected families

### Missing diagrams

If some hits don't produce diagrams:
- The Rfam family may not have a template in R2DT
- Check `output/rfam/` for error messages in the draw logs
- Some very divergent sequences may fail alignment

## Gallery

Examples of RNA structure annotations for different viral genomes, generated with `--normalize-font-size` for consistent visual appearance.

:::{note}
All gallery images are regenerated from example inputs with `just docs-images`. See [Updating documentation](docs.md#regenerating-doc-images) for details.
:::

### SARS-CoV-2 coronavirus

Genome: [OX309346.1](https://www.ebi.ac.uk/ena/browser/view/OX309346.1) (29,903 nt) | RNA structures: 3

```{raw} html
<div style="border: 1px solid #ddd; border-radius: 4px; margin: 1em 0; position: relative;">
  <div style="background: #f5f5f5; padding: 8px 12px; border-bottom: 1px solid #ddd; font-size: 0.85em; color: #666;">
    🔍 Use mouse wheel to zoom, drag to pan. <button id="covid-reset" style="margin-left: 10px; cursor: pointer;">Reset view</button>
  </div>
  <div id="covid-container" style="width: 100%; height: 400px; overflow: hidden;">
    <object id="covid-svg" type="image/svg+xml" data="_images/coronavirus-stitched.svg?v=2" style="width: 100%; height: 100%;"></object>
  </div>
</div>
```

### Hepatitis C virus (HCV)

Genome: [NC_038882.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_038882.1) (9,646 nt) | RNA structures: 12

```{raw} html
<div style="border: 1px solid #ddd; border-radius: 4px; margin: 1em 0; position: relative;">
  <div style="background: #f5f5f5; padding: 8px 12px; border-bottom: 1px solid #ddd; font-size: 0.85em; color: #666;">
    🔍 Use mouse wheel to zoom, drag to pan. <button id="hcv-reset" style="margin-left: 10px; cursor: pointer;">Reset view</button>
  </div>
  <div id="hcv-container" style="width: 100%; height: 400px; overflow: hidden;">
    <object id="hcv-svg" type="image/svg+xml" data="_images/hcv-stitched.svg?v=2" style="width: 100%; height: 100%;"></object>
  </div>
</div>
```

### Dengue virus serotype 2

Genome: [NC_001474.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_001474.2) (10,723 nt) | RNA structures: 4

```{raw} html
<div style="border: 1px solid #ddd; border-radius: 4px; margin: 1em 0; position: relative;">
  <div style="background: #f5f5f5; padding: 8px 12px; border-bottom: 1px solid #ddd; font-size: 0.85em; color: #666;">
    🔍 Use mouse wheel to zoom, drag to pan. <button id="dengue-reset" style="margin-left: 10px; cursor: pointer;">Reset view</button>
  </div>
  <div id="dengue-container" style="width: 100%; height: 400px; overflow: hidden;">
    <object id="dengue-svg" type="image/svg+xml" data="_images/dengue2-stitched.svg?v=2" style="width: 100%; height: 100%;"></object>
  </div>
</div>

<script src="https://cdn.jsdelivr.net/npm/svg-pan-zoom@3.6.1/dist/svg-pan-zoom.min.js"></script>
<script>
(function() {
  function initPanZoom(objectId, resetButtonId) {
    var obj = document.getElementById(objectId);
    var resetBtn = document.getElementById(resetButtonId);
    if (!obj) { console.log('Object not found:', objectId); return; }

    function setup() {
      try {
        var svgDoc = obj.contentDocument;
        if (svgDoc) {
          var svgElement = svgDoc.querySelector('svg');
          if (svgElement && typeof svgPanZoom !== 'undefined') {
            var pz = svgPanZoom(svgElement, {
              zoomEnabled: true,
              controlIconsEnabled: true,
              fit: true,
              center: true,
              minZoom: 0.1,
              maxZoom: 10
            });
            if (resetBtn) {
              resetBtn.onclick = function() { pz.reset(); };
            }
            console.log('Pan/zoom initialized for', objectId);
          }
        }
      } catch(e) {
        console.log('Error initializing', objectId, e);
      }
    }

    // Try immediately (in case already loaded)
    if (obj.contentDocument && obj.contentDocument.querySelector('svg')) {
      setup();
    } else {
      // Otherwise wait for load
      obj.addEventListener('load', setup);
    }
  }

  // Wait for DOM to be ready
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', function() {
      initPanZoom('covid-svg', 'covid-reset');
    initPanZoom('hcv-stk-svg', 'hcv-stk-reset');
    initPanZoom('hcv-svg', 'hcv-reset');
    initPanZoom('dengue-svg', 'dengue-reset');
    });
  } else {
    initPanZoom('covid-svg', 'covid-reset');
    initPanZoom('hcv-stk-svg', 'hcv-stk-reset');
    initPanZoom('hcv-svg', 'hcv-reset');
    initPanZoom('dengue-svg', 'dengue-reset');
  }
})();
</script>
```

```{figure} images/hcv-stitched.svg
:class: hidden
```

```{figure} images/dengue2-stitched.svg
:class: hidden
```

## See also

- [Stitching multiple diagrams](usage.md#stitching-multiple-diagrams) - Full stitch command reference
- [Rfam cmscan documentation](https://docs.rfam.org/en/latest/genome-annotation.html) - Rfam genome annotation guide
