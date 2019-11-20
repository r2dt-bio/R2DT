
# Auto Traveler

Auto Traveler automatically generates RNA secondary structure in standard layouts using templates from the following sources:

 - [CRW](http://www.rna.ccbb.utexas.edu) (5S and SSU rRNA)
 - [Rfam](http://rfam.org) (>2,000 RNA families)
 - [RiboVision](http://apollo.chemistry.gatech.edu/RiboVision/#) (LSU rRNA)

**RNAcentral** uses Auto Traveler to visualise RNA secondary structures. For more details see [RNAcentral help](https://rnacentral.org/help/secondary-structure) or [browse all secondary  structures](https://rnacentral.org/search?q=has_secondary_structure:%22True%22).

## Method overview

1. **Generate a library of covariance models** using bpseq files from [CRW](http://www.rna.icmb.utexas.edu/DAT/3C/Structure/index.php), RiboVision or another source with [Infernal](http://eddylab.org/infernal/). For best results, remove pseudoknots from the secondary structures using [RemovePseudoknots](https://rna.urmc.rochester.edu/Text/RemovePseudoknots.html) from the RNAStructure package.
1. **Select the best matching covariance model** for each input sequence
using [Ribotyper](https://github.com/nawrockie/ribotyper-v1).
1. **Fold** input sequence into a secondary structure compatible with the template
using the top scoring covariance model.
1. **Generate secondary structure diagrams** using [Traveler](https://github.com/davidhoksza/traveler) and the secondary structure layouts.

## Installation

![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/rnacentral/auto-traveler)

Pull from [Docker Hub](https://hub.docker.com/r/rnacentral/auto-traveler):

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

## Usage

Run Docker container and mount the current directory inside the container:

```
docker run -it -v `pwd`:/rna/auto-traveler rnacentral/auto-traveler
```

Perform one-time initial setup:

```
cd auto-traveler
auto-traveler.py setup

# run tests to verify that the installation worked
python3 -m unittest
```

Run examples:

```
python3 auto-traveler.py draw examples/examples.fasta temp/examples
```

For performance reasons you can run the following commands:
```
auto-traveler.py crw draw examples/crw-examples.fasta temp/crw-examples
auto-traveler.py ribovision draw examples/lsu-examples.fasta temp/lsu-examples
auto-traveler.py rfam draw RF00162 examples/RF00162.example.fasta temp/rfam-example
```

Additional commands:

```
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
- [David Mathews lab](http://rna.urmc.rochester.edu/RNAstructure.html)
- [Elena Rivas](https://twitter.com/RivasElenaRivas)
