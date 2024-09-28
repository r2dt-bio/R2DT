# Annotations and data layers

Starting with version 1.4, R2DT can visualise additional annotations as **data layers** on top of the secondary structure diagrams.

By default, R2DT displays the **expected accuracy** of the alignment between the query sequence and the template. The accuracy is calculated by the [Infernal software](http://eddylab.org/infernal) as **posterior probabilities** (see page 110 of the [Infernal User Guide](http://eddylab.org/infernal/Userguide.pdf) for more information).

## Alignment accuracy

For example, the following diagram shows the secondary structure of a [SAM riboswitch RNA sequence](https://rnacentral.org/rna/URS00002D29F6/224308) from *Bacillus subtilis* visualised using the [RF00162](https://rfam.org/family/RF00162) Rfam template:

```{figure} images/URS00002D29F6_224308-RF00162.enriched.svg
:width: 600px
:alt: SAM riboswitch

The colored circles highlight the hairpin loop region where Infernal is not 100% confident in the alignment.
```

Hovering over the nucleotides shows the alignment accuracy (posterior probability) values. Similar to Infernal `#=GR PP` lines, the posterior probability `p` is encoded as 11 possible characters from `0` to `10`:

- 0.0 ≤ p < 0.05 is coded as `0`
- 0.05 ≤ p < 0.15 is coded as `1` (... and so on ...)
- 0.85 ≤ p < 0.95 is coded as `9`
- 0.95 ≤ p ≤ 1.0 (shown as `*` in Infernal output) is coded as `10`.

The nucleotides are colored based on a color-blind friendly palette from a [paper by Bang Wong](https://www.nature.com/articles/nmeth.1618). The default colors can be found in the  [colorscheme.json](https://github.com/r2dt-bio/R2DT/blob/develop/utils/colorscheme.json) file and are shown below:

- <span style="background-color: rgb(213, 94, 0); padding: 4px;">0-5 or <55% accuracy</span>
- <span style="background-color: rgb(230, 159, 0); padding: 4px;">6 or ≥55-65% accuracy</span>
- <span style="background-color: rgb(240, 228, 66); padding: 4px;">7 or ≥65-75% accuracy</span>
- <span style="background-color: rgb(86, 180, 233); padding: 4px;">8 or ≥75-85% accuracy</span>
- <span style="background-color: rgb(0, 158, 115); padding: 4px;">9 or ≥95% accuracy</span>
- `10` and `*` are shown with white background.


## Example annotation file

The annotations can be provided in a tab-separated format where the first two columns are called `residue_index` and `residue_name`. The number of rows in the file (not counting the header line) and residue names must match the query sequence.

The third column (in this case `posterior_probability`) can have an arbitrary name matching the `colorscheme.json`.

The following example is shortened for display purposes but a [complete TSV file](./files/URS00002D29F6_224308_post_prob.txt) is available for reference:

```
residue_index	residue_name	posterior_probability
1	U	10
2	U	10
3	C	10
...
50	G	10
51	G	10
52	U	10
53	G	9
54	U	9
55	A	9
56	A	8
57	U	8
58	G	7
59	G	8
60	C	8
61	G	8
62	A	8
63	U	8
64	C	8
65	A	8
66	G	8
67	C	8
68	C	8
69	A	7
70	U	7
71	G	9
72	A	10
73	C	10
...
117	A	10
118	A	10
```

## Example color scheme file

The colors for each possible `posterior_probability` value from the TSV file are configured using a `colorscheme.json` file, for example:

```
{
    "coloring": {
      "posterior_probability": {
        "label": "Alignment accuracy",
        "values": {
          "0": "rgb(213, 94, 0)",
          "1": "rgb(213, 94, 0)",
          "2": "rgb(213, 94, 0)",
          "3": "rgb(213, 94, 0)",
          "4": "rgb(213, 94, 0)",
          "5": "rgb(213, 94, 0)",
          "6": "rgb(230, 159, 0)",
          "7": "rgb(240, 228, 66)",
          "8": "rgb(86, 180, 233)",
          "9": "rgb(0, 158, 115)",
          "10": "rgb(255, 255, 255)",
          "*": "rgb(255, 255, 255)"
        }
      }
    }
}
```

Note that the header of the third column of the annotation file must match the corresponding field of `colorscheme.json`.

## Example R2DT commands

You will need 3 files:

- `query.fasta` - a query sequence `rna_id` in fasta format,
- `annotations.tsv` - per-residue annotations in TSV format,
- `colorscheme.json` - RGB colors for each annotation category in JSON format.

1. Visualise secondary structure using R2DT:
    ```bash
    # predict and visualise secondary structure using templates
    r2dt.py draw <query.fasta> <output_folder>

    # alternatively, visualise provided secondary structure without a template
    r2dt.py templatefree <query.fasta> <output_folder>
    ```
2. Generate additional SVG file colored using custom annotations:

    ```bash
	python3 /rna/traveler/utils/enrich_json.py --input-json <output_folder>/json/<rna_id>.colored.json \
	--input-data <annotations.tsv> --output <output_folder>/json/<rna_id>.enriched.json

	python3 /rna/traveler/utils/json2svg.py -p <colorscheme.json> \
	-i <output_folder>/json/<rna_id>.enriched.json -o <output_folder>/svg/<rna_id>.enriched.svg
    ```

	The final result will be found in `<output_folder>/svg/<rna_id>.enriched.svg`.

These commands add metadata from the TSV file to an [RNA 2D JSON Schema](https://github.com/LDWLab/RNA2D-data-schema/) file using the [enrich_json.py](https://github.com/cusbg/traveler/tree/master/utils) script and then generate an SVG file using [json2svg.py](https://github.com/cusbg/traveler/tree/master/utils) (both scripts are part of the Traveler software and are available in the R2DT Docker container).

## Displaying other types of data

It is possible to visualise other residue annotations, such as the location of single nucleotide polymorphisms (SNPs), conservation, SHAPE reactivity, location of RNA protein binding sites, and other information.

To visualise custom data:

- the annotations should be provided in the third column of the TSV file,
- the `posterior_probability` field of the `colorscheme.json` file should match the name of the third column of the TSV file,
- the `values` field of `colorscheme.json` should match the possible annotation values.
