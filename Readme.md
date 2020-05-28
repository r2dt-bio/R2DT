
# Auto Traveler

Auto Traveler automatically generates RNA secondary structure in standard layouts using templates from the following sources:

 - [CRW](http://www.rna.ccbb.utexas.edu) (5S and SSU rRNA)
 - [Rfam](https://rfam.org) (>2,000 RNA families)
 - [RiboVision](http://apollo.chemistry.gatech.edu/RiboVision/#) (LSU rRNA)
 - [GtRNAdb](http://gtrnadb.ucsc.edu) (tRNA)

**RNAcentral** uses Auto Traveler to visualise RNA secondary structures. For more details see [RNAcentral help](https://rnacentral.org/help/secondary-structure) or [browse all secondary  structures](https://rnacentral.org/search?q=has_secondary_structure:%22True%22).

## Method overview

1. **Generate a library of covariance models** using bpseq files from [CRW](http://www.rna.icmb.utexas.edu/DAT/3C/Structure/index.php), RiboVision or another source with [Infernal](http://eddylab.org/infernal/). For best results, remove pseudoknots from the secondary structures using [RemovePseudoknots](https://rna.urmc.rochester.edu/Text/RemovePseudoknots.html) from the RNAStructure package.
1. **Select the best matching covariance model** for each input sequence
using [Ribovore](https://github.com/nawrockie/ribovore) or [tRNAScan-SE 2.0](http://lowelab.ucsc.edu/tRNAscan-SE/).
1. **Fold** input sequence into a secondary structure compatible with the template
using the top scoring covariance model.
1. **Generate secondary structure diagrams** using [Traveler](https://github.com/davidhoksza/traveler) and the secondary structure layouts.

## Installation

![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/rnacentral/auto-traveler)

1. Pull from [Docker Hub](https://hub.docker.com/r/rnacentral/auto-traveler):

    ```
    docker pull rnacentral/auto-traveler
    ```

    or build your own Docker image:

    ```
    # Get the code:
    git clone https://github.com/RNAcentral/auto-traveler.git
    cd auto-traveler

    # Build and tag a Docker image:
    docker build -t rnacentral/auto-traveler .
    ```

2. Run the container:

    ```
    docker-compose run cli
    ```

    This command mounts the current directory so all code or data changes are instantly reflected in the container.

## Initial setup

Perform one-time initial setup:

```
auto-traveler.py setup
```

Alternatively, you can download a [precomputed data library](https://www.dropbox.com/s/ypx64w3vfa8z219/cms.tar.gz?dl=0), uncompress and mount it in the container:

```
docker run -it -v `pwd`:/rna/auto-traveler -v <path_to_data_library>:/rna/auto-traveler/data/cms rnacentral/auto-traveler
```

Run tests to verify that the installation worked:
```
python3 -m unittest
```

## Usage

Run examples:

```
auto-traveler.py draw examples/examples.fasta temp/examples
```

To bypass classification steps, run the following commands:
```
auto-traveler.py crw draw examples/crw-examples.fasta temp/crw-examples
auto-traveler.py ribovision draw_lsu examples/lsu-examples.fasta temp/lsu-examples
auto-traveler.py ribovision draw_Ssu examples/ribovision-ssu-examples.fasta temp/ssu-examples
auto-traveler.py rfam draw RF00162 examples/RF00162.example.fasta temp/rfam-example

# for tRNAs, provide domain and isotype (if known), or use tRNAScan-SE to classify
auto-traveler.py gtrnadb draw examples/gtrnadb.E_Thr.fasta temp/gtrnadb
auto-traveler.py gtrnadb draw examples/gtrnadb.E_Thr.fasta temp/gtrnadb --domain E --isotype Thr
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

## Acknowledgements

- [David Hoksza](https://github.com/davidhoksza)
- [Eric Nawrocki](https://github.com/nawrockie)
- [Robin Gutell lab](http://www.rna.ccbb.utexas.edu)
- [Anton S. Petrov](https://cool.gatech.edu/people/petrov-anton) and the [RiboVision](http://apollo.chemistry.gatech.edu/RiboVision/#) team
- [Todd Lowe](https://users.soe.ucsc.edu/~lowe/) and [Patricia Chan](https://www.soe.ucsc.edu/people/pchan)
- [David Mathews lab](http://rna.urmc.rochester.edu/RNAstructure.html)
- [Elena Rivas](https://twitter.com/RivasElenaRivas)
