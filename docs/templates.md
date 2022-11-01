# Creating and updating templates

## Adding new templates

If you would like to submit a new template or replace an existing one, please [submit an issue](https://github.com/RNAcentral/R2DT/issues/new) on GitHub including:

1. A FASTA or BPSEQ file with a reference sequence and secondary structure - see [FASTA](https://github.com/RNAcentral/R2DT/blob/master/data/rfam/RF00012/RF00012-traveler.fasta) and [BPSEQ](https://github.com/RNAcentral/R2DT/blob/master/data/ribovision-ssu/bpseq/EC_SSU_3D.bpseq) examples
1. A [Traveler XML file](https://github.com/cusbg/traveler#traveler-intermediate-format) - see [example](https://github.com/RNAcentral/R2DT/blob/master/data/rfam/RF00003/traveler-template.xml)
1. Description of the new template and any relevant background information

One can create a new template locally using the [generate_cm_library.py](https://github.com/RNAcentral/R2DT/blob/master/utils/generate_cm_library.py) script with the FASTA and XML files described above. It is also possible to generate a new template using a special version of the XRNA software, [XRNA-GT](https://github.com/LDWLab/XRNA-GT).

:warning: GitHub currently does not support attaching files with `.fasta` or `.bpseq` extensions so please attach the files as `.txt`.

We will review the template and reply on GitHub as soon as possible.

## Manually creating templates

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
    tar -czvf cms.tar.gz <path/to/new/cms>
    ```
    The folder should contain 2 subfolders: `crw` and `rfam`.
1. Update the precomputed library link in Readme
1. Update a list of available models
    ```bash
    r2dt.py list-models
    ```
1. ⚠️ Note that the [tRNA Rfam Traveler template](https://github.com/RNAcentral/R2DT/blob/master/data/rfam/RF00005/traveler-template.xml) has been manually edited to match the standard tRNA layout so the automatically generated `traveler-template.xml` file should be discarded and the current version should be kept.