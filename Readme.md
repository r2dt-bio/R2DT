
# R2DT

The R2DT software (RNA 2D Templates) automatically generates RNA secondary structure in standard layouts using templates from the following sources:

 - [CRW](http://www.rna.ccbb.utexas.edu) (5S and SSU rRNA)
 - [Rfam](https://rfam.org) (>2,000 RNA families)
 - [RiboVision](http://apollo.chemistry.gatech.edu/RiboVision/#) (LSU rRNA)
 - [GtRNAdb](http://gtrnadb.ucsc.edu) (tRNA)

**RNAcentral** uses R2DT to visualise RNA secondary structures. For more details see the [R2DT preprint](https://www.biorxiv.org/content/10.1101/2020.09.10.290924v1) or [browse all secondary  structures](https://rnacentral.org/search?q=has_secondary_structure:%22True%22).

## Method overview

1. **Generate a library of covariance models** using bpseq files from [CRW](http://www.rna.icmb.utexas.edu/DAT/3C/Structure/index.php), RiboVision or another source with [Infernal](http://eddylab.org/infernal/). For best results, remove pseudoknots from the secondary structures using [RemovePseudoknots](https://rna.urmc.rochester.edu/Text/RemovePseudoknots.html) from the RNAStructure package.
1. **Select the best matching covariance model** for each input sequence
using [Ribovore](https://github.com/nawrockie/ribovore) or [tRNAScan-SE 2.0](http://lowelab.ucsc.edu/tRNAscan-SE/).
1. **Fold** input sequence into a secondary structure compatible with the template
using the top scoring covariance model.
1. **Generate secondary structure diagrams** using [Traveler](https://github.com/davidhoksza/traveler) and the secondary structure layouts.

## Installation

![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/rnacentral/r2dt)

1. Pull from [Docker Hub](https://hub.docker.com/r/rnacentral/r2dt):

    ```
    docker pull rnacentral/r2dt
    ```

    or build your own Docker image:

    ```
    # Get the code:
    git clone https://github.com/RNAcentral/R2DT.git
    cd R2DT

    # Build and tag a Docker image:
    docker build -t rnacentral/r2dt .
    ```

2. Run the container:

    ```
    docker-compose run cli
    ```

    This command mounts the current directory so all code or data changes are instantly reflected in the container.

## Initial setup

Perform one-time initial setup:

```
r2dt.py setup
```

Alternatively, you can download a [precomputed data library](https://www.dropbox.com/s/q5l0s1nj5h4y6e4/cms.tar.gz?dl=0), uncompress and mount it in the container:

```
docker run -it -v `pwd`:/rna/r2dt -v <path_to_data_library>:/rna/r2dt/data/cms rnacentral/r2dt
```

Run tests to verify that the installation worked:
```
python3 -m unittest
```

## Usage

Run examples:

```
r2dt.py draw examples/examples.fasta temp/examples
```

To bypass classification steps, run the following commands:
```
r2dt.py crw draw examples/crw-examples.fasta temp/crw-examples
r2dt.py ribovision draw_lsu examples/lsu-examples.fasta temp/lsu-examples
r2dt.py ribovision draw_Ssu examples/ribovision-ssu-examples.fasta temp/ssu-examples
r2dt.py rfam draw RF00162 examples/RF00162.example.fasta temp/rfam-example

# for tRNAs, provide domain and isotype (if known), or use tRNAScan-SE to classify
r2dt.py gtrnadb draw examples/gtrnadb.E_Thr.fasta temp/gtrnadb
r2dt.py gtrnadb draw examples/gtrnadb.E_Thr.fasta temp/gtrnadb --domain E --isotype Thr
```

Additional commands:

```
# to run individual tests
python3 -m unittest tests.tests.TestRibovisionLSU

# classify example sequences using Ribotyper
perl ribotyper.pl -i data/cms/all.modelinfo.txt -f examples/pdb.fasta example-output

# to generate covariance models:
python3 utils/generate_cm_library.py
python3 utils/generate_lsu_cm_library.py

python3 utils/generate_model_info.py
python3 utils/generate_model_info.py --cm-library=data/ribovision/cms --rna-type=LSU
```

## How to add new templates

If you would like to submit a new template or replace an existing one, please [submit an issue](https://github.com/RNAcentral/R2DT/issues/new) including:

- Description of the new template and any relevant background information
- A FASTA file with a reference sequence and secondary structure - see [example](./data/rfam/RF00002/RF00002-traveler.fasta)
- A [Traveler XML file](https://github.com/davidhoksza/traveler#traveler-intermediate-format) - see [example](./data/rfam/RF00002/traveler-template.xml)

:warning: GitHub currently does not support attaching files with `.fasta` or `.bpseq` extensions so please attach the files as `.txt`.

We will review the template and reply on GitHub as soon as possible.

## Acknowledgements

- [David Hoksza](https://github.com/davidhoksza)
- [Eric Nawrocki](https://github.com/nawrockie)
- [Robin Gutell lab](http://www.rna.ccbb.utexas.edu)
- [Anton S. Petrov](https://cool.gatech.edu/people/petrov-anton) and the [RiboVision](http://apollo.chemistry.gatech.edu/RiboVision/#) team
- [Todd Lowe](https://users.soe.ucsc.edu/~lowe/) and [Patricia Chan](https://www.soe.ucsc.edu/people/pchan)
- [David Mathews lab](http://rna.urmc.rochester.edu/RNAstructure.html)
- [Elena Rivas](https://twitter.com/RivasElenaRivas)
