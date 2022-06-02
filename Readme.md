# R2DT

_Visualise RNA 2D structure in standard layouts_

**Version 1.2 (August 10, 2021)**

The R2DT software (RNA 2D Templates) automatically generates [RNA secondary structure](https://en.wikipedia.org/wiki/Nucleic_acid_secondary_structure) diagrams in standard layouts using a template library representing a wide range of RNAs:

 - 5S and SSU rRNA from [CRW](http://www.rna.ccbb.utexas.edu)
 - 3D-structure based SSU and LSU rRNA from [RiboVision](http://apollo.chemistry.gatech.edu/RiboVision/#)
 - tRNA from [GtRNAdb](http://gtrnadb.ucsc.edu)
 - RNAse P from [Ribonuclease P Database](https://academic.oup.com/nar/article/26/1/351/2379438)
 - over 4,000 RNA families from [Rfam](https://rfam.org) (release 14.8)

![R2DT method overview](./examples/method-overview.png)

R2DT is used by RNAcentral to visualise [>20 million RNA secondary structures](https://rnacentral.org/search?q=has_secondary_structure:%22True%22). See [method overview](#method-overview) for details or read the [R2DT paper](https://www.nature.com/articles/s41467-021-23555-5) in Nature Communications.

## Examples

The following example visualisations show LSU, SSU, and 5S rRNA, four tRNAs, two RNAse P, snoRNA, MoCo riboswitch, and U4 snRNA.

![R2DT examples](./examples/r2dt-examples.png)

## Getting started

R2DT can be used in a number of ways:

* [Web application](https://rnacentral.org/r2dt) hosted by RNAcentral
* [API](https://www.ebi.ac.uk/Tools/common/tools/help/index.html?tool=r2dt) powered by EMBL-EBI Web Services
* As a command line tool with [Docker](https://www.docker.com), [Singularity](https://sylabs.io/docs/#singularity), or in a bare metal installation

### Installation

![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/rnacentral/r2dt)

* Download the R2DT image from [Docker Hub](https://hub.docker.com/r/rnacentral/r2dt) and run it with Docker or Singularity.

    **Docker**

    ```
    docker pull rnacentral/r2dt
    docker run --entrypoint r2dt.py rnacentral/r2dt draw --help
    ```

    **Singularity**

    ```
    singularity build r2dt docker://rnacentral/r2dt    
    singularity exec r2dt r2dt.py draw --help
    ```

* :hammer_and_wrench: Development installation:
    ```
    # Get the code
    git clone https://github.com/RNAcentral/R2DT.git
    cd R2DT

    # Build and tag a Docker image
    docker build -t rnacentral/r2dt .
    docker-compose run cli
    ```
    The current directory is mounted inside the container so that all code and data changes are instantly reflected in the container.

* :hammer_and_wrench: Bare metal installation: if running R2DT using containers is not possible, follow instructions in the [Dockerfile](./Dockerfile).

### Initial setup

1. Download a [precomputed data library](https://ftp.ebi.ac.uk/pub/databases/RNAcentral/r2dt/1.2/cms.tar.gz) _(197 MB, last updated Aug 9, 2021)_ and uncompress it.

2. Enter an interactive Docker terminal session:

    ```
    docker run -it -v <path_to_cms>:/rna/r2dt/data/cms -v `pwd`:/rna/r2dt/temp rnacentral/r2dt
    ```

    - `-it` - start an interactive session
    - `-v <path_to_cms>:/rna/r2dt/data/cms` - mount the precomputed data library folder `<path_to_cms>` as `/rna/r2dt/data/cms` inside the container. :warning: Note that `<path_to_cms>` should be a full path.
    - make the current working directory available inside the container as `/rna/r2dt/temp`:

        ```
        -v `pwd`:/rna/r2dt/temp
        ```

Any file placed in `/rna/r2dt/temp` within the container will be available on the host machine after the Docker container exits.

## Usage

### Automatic template selection

Specify the input file in [FASTA format](https://en.wikipedia.org/wiki/FASTA_format) containing one or more RNA sequences as well as the path where the output files will be created (the folder will be created if it does not exist).

```
r2dt.py draw <input.fasta> <output_folder>
```

For example:

```
r2dt.py draw examples/examples.fasta temp/examples
```

R2DT will automatically select the best matching template and visualise the secondary structures.

### Manually selecting template category

If the RNA type of the input sequences is known in advance, it is possible to bypass the classification steps and achieve faster performance.

* CRW templates (5S and SSU rRNA)
    ```
    r2dt.py crw draw examples/crw-examples.fasta temp/crw-examples
    ```

* RiboVision LSU and SSU rRNA templates
    ```
    r2dt.py ribovision draw_lsu examples/lsu-examples.fasta temp/lsu-examples
    r2dt.py ribovision draw_ssu examples/ribovision-ssu-examples.fasta temp/ssu-examples
    ```

* Rfam families
    ```
    r2dt.py rfam draw RF00162 examples/RF00162.example.fasta temp/rfam-example
    ```

* RNAse P
    ```
    r2dt.py rnasep draw examples/rnasep.fasta temp/rnasep-example
    ```

* tRNAs (using GtRNAdb templates)
    ```
    # for tRNAs, provide domain and isotype (if known), or use tRNAScan-SE to classify
    r2dt.py gtrnadb draw examples/gtrnadb.E_Thr.fasta temp/gtrnadb
    r2dt.py gtrnadb draw examples/gtrnadb.E_Thr.fasta temp/gtrnadb --domain E --isotype Thr
    ```

### Manually selecting templates

It is possible to select a specific template and skip the classification step altogether.

1. Get a list of all available templates and copy the template id:
    ```
    r2dt.py list-models
    ```

In addition, all models are listed in the file [models.json](./data/models.json).

2. Specify the template (for example, `RNAseP_a_P_furiosus_JB`):
    ```
    r2dt.py draw --force_template <template_id> <input_fasta> <output_folder>
    ```

    For example:

    ```
    r2dt.py draw --force_template RNAseP_a_P_furiosus_JB examples/force/URS0001BC2932_272844.fasta temp/example
    ```
### Constraint-based folding for insertions

If a structure contains insertions that are not present in the R2DT template files, the --constraint flag will allow their folding to be de-novo predicted using the RNAfold algorithm.

There are currently three constraint folding modes available. R2DT will automatically predict which folding mode is best for a given molecule, but the mode can also be manually overridden using the --fold_type parameter. There are three options for fold_type.

* Let R2DT pick a fold_type  
    ```
    r2dt.py draw --constraint <input_fasta> <output_folder>
    ```
* Fold insertions (along with adjacent unpaired nucleotides) one at a time. Recommended for large RNAs.
    ```
    r2dt.py draw --constraint --fold_type insertions_only <input_fasta> <output_folder>
    ```
* Run entire molecule through RNAfold at once. Base pairs predicted from the template are used as constraints for prediction.
    ```
    r2dt.py draw --constraint --fold_type full_molecule <input_fasta> <output_folder>
    ```
* Run entire molecule through RNAfold at once. Both conserved single-stranded regions and base pairs predicted from the template are used as constraints for prediction.
    ```
    r2dt.py draw --constraint --fold_type all_constraints_enforced <input_fasta> <output_folder>
    ```
* Prevent certain nucleotides from base pairing. This will only work for base pairs that are de-novo predicted.
The exclusion file should contain a string the same length as the input sequence composed of '.'s and 'x's. Positions with '.'s are allowed to base pair,
positions with 'x's are not.
Example string: 'xxxx..............xx..............x............xx'
    ```
    r2dt.py draw --constraint --exclusion <exclusion_file> <input_fasta> <output_folder>
    ```

### Other useful commands

* Print R2DT version
    ```
    r2dt.py version
    ```

* Run [all tests](./tests/tests.py)
    ```
    r2dt.py test
    ```

    or

    ```
    python3 -m unittest
    ```

* Run a [single test](./tests/tests.py)
    ```
    python3 -m unittest tests.tests.TestRibovisionLSU
    ```

* Run tests and keep the results (useful when updating Traveler for example)
    ```
    R2DT_KEEP_TEST_RESULTS=1 r2dt.py test
    R2DT_KEEP_TEST_RESULTS=1 python3 -m unittest tests.tests.TestRnasep
    ```

* Classify example sequences using Ribotyper
    ```
    perl /rna/ribovore/ribotyper.pl -i data/cms/crw/modelinfo.txt -f examples/pdb.fasta temp/ribotyper-test
    ```

* Generate covariance models and modelinfo files
    ```    
    python3 utils/generate_cm_library.py
    r2dt.py generatemodelinfo <path to covariance models>
    ```

* Precompute template library locally (may take several hours):
    ```
    r2dt.py setup
    ```

* Run R2DT with Singularity
    ```
    singularity exec --bind <path_to_cms>:/rna/r2dt/data/cms r2dt r2dt.py draw sequence.fasta output
    ```

* Convert a SVG diagram to a JSON file containing the paths per nucleotide and an ordinal numbering. Note that this *assumes* that the input pdb id is formatted like: `<PDB>_<Entity>_<chain>`, ie `1S72_1_0`.
    ```
    svg2json.py <pdb-id> diagram.svg <pdb-id>.json
    ```

## Output files

`r2dt.py draw` produces a folder called `results` with the following subfolders:

- `svg`: RNA secondary structure diagrams in SVG format
- `fasta`: input sequences and their secondary structure in dot-bracket notation
- `tsv`: a file `metadata.tsv` listing sequence ids, matching templates, and template sources
- `thumbnail`: secondary structure diagrams displayed as outlines in SVG format
- `json`: RNA secondary structure and its layout described using [RNA 2D JSON Schema](https://github.com/LDWLab/RNA2D-data-schema/)

## How to add new templates

If you would like to submit a new template or replace an existing one, please [submit an issue](https://github.com/RNAcentral/R2DT/issues/new) including:

- A FASTA or BPSEQ file with a reference sequence and secondary structure - see [FASTA](./data/rfam/RF00002/RF00002-traveler.fasta) and [BPSEQ](./data/ribovision-ssu/bpseq/EC_SSU_3D.bpseq) examples
- A [Traveler XML file](https://github.com/davidhoksza/traveler#traveler-intermediate-format) - see [example](./data/rfam/RF00002/traveler-template.xml)
- Description of the new template and any relevant background information

One can create a new template locally using the [generate_cm_library.py](./utils/generate_cm_library.py) script with the FASTA and XML files described above. It is also possible to generate a new template using a special version of the XRNA software, [XRNA-GT](https://github.com/LDWLab/XRNA-GT).

:warning: GitHub currently does not support attaching files with `.fasta` or `.bpseq` extensions so please attach the files as `.txt`.

We will review the template and reply on GitHub as soon as possible.

### Manually creating templates

1. Place new FASTA or BPSEQ file(s) in the `data/new` folder
1. Run `r2dt.py generatecm`. The command will generate new `.cm` file(s) with the covariance models
1. Move the new `.cm` and `.tr` files in the destination directory (for example, `ribovision-ssu` is where all SSU templates submitted by the RiboVision group are stored)
1. Run `r2dt.py generatemodelinfo <destination>/cms` to add new models to the list of searched models where `<destination>` is the same folder as in the step above
1. Update `metadata.tsv` file in the destination directory
1. Run `r2dt.py list-models` to update a list of all available models
1. Verify that the templates work as expected by testing on a fasta file with a sequence similar to the template
1. Run tests and update the model counts in `TestCovarianceModelDatabase` as needed

### Updating Rfam library

The following procedure should be done after each Rfam release:

1. Recompute all Rfam templates (takes ~6h)
    ```
    r2dt.py setup-rfam
    ```

1. Run tests

1. Generate new precomputed library archive
    ```
    tar -czvf cms.tar.gz <path/to/new/cms>
    ```

    The folder should contain 3 subfolders: `crw`, `gtrnadb`, and `rfam`.

1. Update the precomputed library link in Readme

1. Update a list of available models
    ```
    r2dt.py list-models
    ```

1. :warning: Note that the [tRNA Rfam Traveler template](./data/rfam/RF00005/traveler-template.xml) has been manually edited to match the standard tRNA layout so the automatically generated `traveler-template.xml` file should be discarded and the current version should be kept.

## Method overview

The R2DT pipeline includes the following steps:

1. **Generate a library of covariance models** using BPSEQ files from [CRW](http://www.rna.icmb.utexas.edu/DAT/3C/Structure/index.php), RiboVision or another source with [Infernal](http://eddylab.org/infernal/). For best results, remove pseudoknots from the secondary structures using [RemovePseudoknots](https://rna.urmc.rochester.edu/Text/RemovePseudoknots.html) from the RNAStructure package.
1. **Select the best matching covariance model** for each input sequence
using [Ribovore](https://github.com/nawrockie/ribovore) or [tRNAScan-SE 2.0](http://lowelab.ucsc.edu/tRNAscan-SE/).
1. **Fold** input sequence into a secondary structure compatible with the template
using the top scoring covariance model.
1. **Generate secondary structure diagrams** using [Traveler](https://github.com/davidhoksza/traveler) and the secondary structure layouts.

See the [R2DT paper](https://www.nature.com/articles/s41467-021-23555-5) for more details.

## Release process

All R2DT releases are [available](https://github.com/RNAcentral/R2DT/releases) on GitHub. R2DT uses [git flow](https://github.com/RNAcentral/R2DT/wiki) for managing the release process.

## Contributors

- [David Hoksza](https://github.com/davidhoksza) (Traveler software)
- [Eric Nawrocki](https://github.com/nawrockie) (Ribovore and Infernal software)
- [Robin Gutell](https://scholar.google.com/citations?user=IdDGv6oAAAAJ&hl=en) and [Jamie Cannone](https://scholar.google.com/citations?user=PnqrMGAAAAAJ&hl=en) (CRW)
- [Anton S. Petrov](https://cool.gatech.edu/people/petrov-anton), [Loren D. Williams](https://cool.gatech.edu/people/williams-loren-dean), and the [RiboVision](http://apollo.chemistry.gatech.edu/RiboVision/#) team
- [Todd Lowe](https://users.soe.ucsc.edu/~lowe/) and [Patricia Chan](https://www.soe.ucsc.edu/people/pchan) from [GtRNAdb](http://gtrnadb.ucsc.edu)
- [Blake Sweeney](https://www.ebi.ac.uk/about/people/blake-sweeney), [Carlos Ribas](https://www.ebi.ac.uk/about/people/carlos-eduardo-ribas), [Fabio Madeira](https://www.ebi.ac.uk/about/people/fabio-madeira), [Rob Finn](https://www.ebi.ac.uk/about/people/rob-finn), [Anton I. Petrov](https://www.ebi.ac.uk/about/people/anton-petrov) from [EMBL-EBI](https://www.ebi.ac.uk)

:wave: We welcome additional contributions so please feel free to [raise an issue](https://github.com/RNAcentral/R2DT/issues) or submit a pull request.

## Acknowledgements

- [Sean Eddy](http://eddylab.org) - [Infernal software](http://eddylab.org/infernal)
- [David Mathews](http://rna.urmc.rochester.edu/RNAstructure.html) - [RNAstructure software](http://rna.urmc.rochester.edu/RNAstructure.html)
- [Elena Rivas](http://rivaslab.org/) - [R-scape software](http://eddylab.org/R-scape/)
- [Zasha Weinberg](https://zashaweinberglab.org) - [R2R software](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-3)
