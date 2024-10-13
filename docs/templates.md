# Creating and updating templates

## Adding new templates

To submit a new template or update an existing one, please [create a GitHub issue](https://github.com/r2dt-bio/R2DT/issues/new) with the following information:

1. A RNA 2D JSON Schema file (recommended) - see [example](https://github.com/r2dt-bio/R2DT/blob/main/examples/RF02976.json).

    Alternatively, you can submit the following:
    1. A FASTA or BPSEQ file with a reference sequence and secondary structure - see [FASTA](https://github.com/r2dt-bio/R2DT/blob/main/data/rfam/RF00012/RF00012-traveler.fasta) and [BPSEQ](https://github.com/r2dt-bio/R2DT/blob/main/data/ribovision-ssu/bpseq/EC_SSU_3D.bpseq) examples
    1. A [Traveler XML file](https://github.com/cusbg/traveler#traveler-intermediate-format) - see [example](https://github.com/r2dt-bio/R2DT/blob/main/data/rfam/RF00003/traveler-template.xml)

1. Description of the new template and any relevant background information.

:::{note}
GitHub currently does not support attaching files with `.fasta`, `.bpseq`, or `.json` extensions so please attach the files as `.txt`.
:::

## Creating templates using RNA 2D JSON Schema files

It is possible to generate new templates using [RNA 2D JSON Schema](https://github.com/LDWLab/RNA2D-data-schema/) files as input.

### XRNA-React

1. Run your sequence through [R2DT](https://r2dt.bio) webserver and click `Edit in XRNA`
    1. Alternatively, generate a structure using R2DT on the command line and locate the JSON file in the `json` folder.
    1. Upload JSON file to an [interactive editor](https://ldwlab.github.io/XRNA-React)
1. Manually edit the layout
1. Download the edited structure as a new JSON file and save it in the `data/new` folder
1. Run the following command:

    ```bash
    r2dt.py generate-template data/new/<input.json>
    ```

1. The new Traveler template, covariance model, and a fasta file will be generated in the `data/local_data` folder.
1. Use the new template: `r2dt.py draw --force_template <template_name> <input.fasta> <output_folder>`

### RNAcanvas

1. Run your sequence through [R2DT](https://r2dt.bio) webserver and click Edit in `RNAcanvas`. Alternatively, open [RNAcanvas](https://rnacanvas.app) and click `Create new drawing`
1. Edit the drawing and click `RNA 2D` to export an RNA 2D JSON Schema file
1. Put the JSON file in any folder accessible to R2DT, for example, `data/new`
1. Run `r2dt.py generate-template data/new/<template_name.json>`
1. The new Traveler template, covariance model, and a fasta file can found in a new folder `data/local_data/template_name`, where `template_name` matches the name of the json file
1. Use the new template: `r2dt.py draw --force_template <template_name> <input.fasta> <output_folder>`

### Creating templates using FASTA/BPSEQ and Traveler XML files

1. Place new FASTA or BPSEQ file(s) in the `data/new` folder
1. Run `r2dt.py generatecm`. The command will generate new `.cm` file(s) with the covariance models
1. Move the new `.cm` and `.tr` files in the destination directory (for example, `ribovision-ssu` is where all SSU templates submitted by the RiboVision group are stored)
1. Run `r2dt.py generatemodelinfo <destination>/cms` to add new models to the list of searched models where `<destination>` is the same folder as in the step above
1. Update `metadata.tsv` file in the destination directory
1. Run `r2dt.py list-models` to update a list of all available models
1. Verify that the templates work as expected by testing on a fasta file with a sequence similar to the template
1. Run tests and update the model counts in `TestCovarianceModelDatabase` as needed

## Updating Rfam templates

The following procedure should be done after each Rfam release:

1. Recompute all Rfam templates (takes ~6h)
    ```bash
    r2dt.py setup-rfam
    ```
1. Run tests
1. Generate new precomputed library archive
    ```bash
    r2dt.py create-precomputed-library <release_number>
    ```
    The folder should contain 2 subfolders: `crw` and `rfam`.
1. Update the precomputed library link in Readme
1. Update a list of available models
    ```bash
    r2dt.py list-models
    ```
1. ⚠️ Note that the [tRNA Rfam Traveler template](https://github.com/r2dt-bio/R2DT/blob/main/data/rfam/RF00005/traveler-template.xml) has been manually edited to match the standard tRNA layout so the automatically generated `traveler-template.xml` file should be discarded and the current version should be kept.
