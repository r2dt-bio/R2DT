# R2DT

_Visualise RNA 2D structure in standard layouts_

The R2DT software (RNA 2D Templates) automatically generates [RNA secondary structure](https://en.wikipedia.org/wiki/Nucleic_acid_secondary_structure) diagrams in standard layouts using a template library that contains >3,500 templates representing a wide range of RNAs from the following sources:

 - [CRW](http://www.rna.ccbb.utexas.edu) (5S and SSU rRNA)
 - [RiboVision](http://apollo.chemistry.gatech.edu/RiboVision/#) (3D-structure based SSU and LSU rRNA)
 - [GtRNAdb](http://gtrnadb.ucsc.edu) (tRNA)
 - [Ribonuclease P Database](https://academic.oup.com/nar/article/26/1/351/2379438) (RNAse P)
 - [Rfam](https://rfam.org) (>2,000 RNA families)

R2DT is used by RNAcentral to visualise [>14 million RNA secondary structures](https://rnacentral.org/search?q=has_secondary_structure:%22True%22). See [method overview](#method-overview) for details or read the [preprint](https://www.biorxiv.org/content/10.1101/2020.09.10.290924v1) on BioRxiv.

## Examples

Here are some example visualisations showing LSU, SSU, and 5S rRNA, several tRNAs, RNAse P and others.

![R2DT examples](./examples/r2dt-examples.png)

## Getting started

R2DT can be used in a number of ways:

* [Web application](https://rnacentral.org/r2dt) hosted by RNAcentral
* [API](https://www.ebi.ac.uk/Tools/common/tools/help/) powered by EMBL-EBI Web Services
* As a command line tool with [Docker](https://www.docker.com), [Singularity](https://sylabs.io/docs/#singularity), or in a bare metal installation

### Installation

![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/rnacentral/r2dt)

* Recommended installation: download the R2DT image from [Docker Hub](https://hub.docker.com/r/rnacentral/r2dt) and run it with Docker or Singularity.

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

1. Download a [precomputed data library](https://www.dropbox.com/s/3ie8kzb8ol658s0/cms.tar.gz?dl=1) _(190.1 MB, last updated Jan 7, 2021)_ and uncompress it.

2. Enter an interactive Docker terminal session:

```
docker run -it -v <path_to_cms>:/rna/r2dt/data/cms -v `pwd`:/rna/r2dt/temp rnacentral/r2dt
```

- `-it` - start an interactive session
- `-v <path_to_cms>:/rna/r2dt/data/cms` - mount the precomputed data library folder `<path_to_cms>` as `/rna/r2dt/data/cms` inside the container
- make the current working directory available inside the container as `/rna/r2dt/temp`:
    ```
    -v `pwd`:/rna/r2dt/temp
    ```

Any file placed in `/rna/r2dt/temp` will be available on the host machine after the Docker container exits.

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

### Specifying template category

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

### Manual template selection

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

### Other useful commands

* Perform one-time initial setup locally (may take up to several hours):
    ```
    r2dt.py setup
    ```

* Run the entire test suite
    ```
    python3 -m unittest
    ```

* Run individual tests
    ```
    python3 -m unittest tests.tests.TestRibovisionLSU
    ```

* Classify example sequences using Ribotyper
    ```
    perl ribotyper.pl -i data/cms/all.modelinfo.txt -f examples/pdb.fasta example-output
    ```

* Generate covariance models and modelinfo files
    ```    
    python3 utils/generate_cm_library.py
    r2dt.py generatemodelinfo <path to covariance models>
    ```

* Run R2DT with Singularity
    ```
    singularity exec --bind <path_to_cms>:/rna/r2dt/data/cms r2dt r2dt.py draw sequence.fasta output
    ```

## How to add new templates

If you would like to submit a new template or replace an existing one, please [submit an issue](https://github.com/RNAcentral/R2DT/issues/new) including:

- Description of the new template and any relevant background information
- A FASTA file with a reference sequence and secondary structure - see [example](./data/rfam/RF00002/RF00002-traveler.fasta)
- A [Traveler XML file](https://github.com/davidhoksza/traveler#traveler-intermediate-format) - see [example](./data/rfam/RF00002/traveler-template.xml)

:warning: GitHub currently does not support attaching files with `.fasta` or `.bpseq` extensions so please attach the files as `.txt`.

We will review the template and reply on GitHub as soon as possible.

## Method overview

1. **Generate a library of covariance models** using bpseq files from [CRW](http://www.rna.icmb.utexas.edu/DAT/3C/Structure/index.php), RiboVision or another source with [Infernal](http://eddylab.org/infernal/). For best results, remove pseudoknots from the secondary structures using [RemovePseudoknots](https://rna.urmc.rochester.edu/Text/RemovePseudoknots.html) from the RNAStructure package.
1. **Select the best matching covariance model** for each input sequence
using [Ribovore](https://github.com/nawrockie/ribovore) or [tRNAScan-SE 2.0](http://lowelab.ucsc.edu/tRNAscan-SE/).
1. **Fold** input sequence into a secondary structure compatible with the template
using the top scoring covariance model.
1. **Generate secondary structure diagrams** using [Traveler](https://github.com/davidhoksza/traveler) and the secondary structure layouts.

## Contributors

- [David Hoksza](https://github.com/davidhoksza) (Traveler software)
- [Eric Nawrocki](https://github.com/nawrockie) (Ribovore and Infernal software)
- [Robin Gutell]() and [Jamie Cannone]() (CRW)
- [Anton S. Petrov](https://cool.gatech.edu/people/petrov-anton), [Loren D. Williams](https://cool.gatech.edu/people/williams-loren-dean), and the [RiboVision](http://apollo.chemistry.gatech.edu/RiboVision/#) team
- [Todd Lowe](https://users.soe.ucsc.edu/~lowe/) and [Patricia Chan](https://www.soe.ucsc.edu/people/pchan) and the GtRNAdb team
- [Blake Sweeney](https://www.ebi.ac.uk/about/people/blake-sweeney), [Carlos Ribas](https://www.ebi.ac.uk/about/people/carlos-eduardo-ribas), [Fabio Madeira](https://www.ebi.ac.uk/about/people/fabio-madeira), [Rob Finn](https://www.ebi.ac.uk/about/people/rob-finn), [Anton I. Petrov](https://www.ebi.ac.uk/about/people/anton-petrov) ([EMBL-EBI](https://www.ebi.ac.uk))

## Acknowledgements

- [Sean Eddy](http://eddylab.org) - [Infernal software](http://eddylab.org/infernal)
- [David Mathews](http://rna.urmc.rochester.edu/RNAstructure.html) - [RNAstructure software](http://rna.urmc.rochester.edu/RNAstructure.html)
- [Elena Rivas](http://rivaslab.org/) - [R-scape software](http://eddylab.org/R-scape/)
- [Zasha Weinberg](https://zashaweinberglab.org) - [R2R software](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-3)
