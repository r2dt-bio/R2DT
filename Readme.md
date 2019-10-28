
# Auto Traveler

This is a tool for automatic generation of RNA secondary structure in standard
[CRW](http://www.rna.ccbb.utexas.edu) layouts.

## Method overview

1. Generate a library of **covariance models** using [CRW bpseq files](http://www.rna.icmb.utexas.edu/DAT/3C/Structure/index.php)
and [Infernal](http://eddylab.org/infernal/). For best results, remove pseudoknots from CRW secondary structures
  using [RemovePseudoknots](https://rna.urmc.rochester.edu/Text/RemovePseudoknots.html) from the RNAStructure package.
1. **Select the best matching covariance model** for each input sequence
using [Ribotyper](https://github.com/nawrockie/ribotyper-v1)
1. **Fold** input sequence into a secondary structure compatible with the template
using the top scoring covariance model.
1. **Generate secondary structure diagrams** using [Traveler](https://github.com/davidhoksza/traveler) and CRW layouts.

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
docker build -t auto-traveler .

# Run Docker container and mount the current directory inside the container:
docker run -it -v `pwd`:/rna/auto-traveler auto-traveler

# generate extra files required by Ribotyper:
cd auto-traveler
python utils/generate_model_info.py
```

## Usage

Run Docker container as shown above, then inside the container:

```
cd auto-traveler

# run traveler on all sequences from an example file:
python auto-traveler.py --fasta-input examples/examples.fasta --output-folder example-output

# place your fasta file in a folder that is mounted in the container:
python auto-traveler.py --fasta-input /path/to/input.fasta --output-folder /path/to/output-folder
```

For Rfam families, the sequence classification step can be skipped
if the Rfam accession is provided:

```
# to process a specific Rfam family and store output in `rfam-output` folder
python auto-traveler.py --rfam-accession RF00162 --output-folder rfam-output

# to process all Rfam families
python auto-traveler.py --rfam-accession all --output-folder rfam-output

# to process sequences from a specific fasta file
python auto-traveler.py --rfam-accession RF00162 --output-folder rfam-output --fasta_input /path/to/fasta

# see help for more options
python auto-traveler.py --help
```

Additional commands:

```
# classify example sequences using Ribotyper
perl ribotyper.pl -i data/cms/all.modelinfo.txt -f examples/pdb.fasta example-output

# to generate covariance models:
perl utils/generate_cm_library.py
```

## Acknowledgements

- [David Hoksza](https://github.com/davidhoksza)
- [Eric Nawrocki](https://github.com/nawrockie)
- [Robin Gutell lab](http://www.rna.ccbb.utexas.edu)
- [David Mathews lab](http://rna.urmc.rochester.edu/RNAstructure.html)

Secondary structure information was downloaded from the [CRW](http://www.rna.ccbb.utexas.edu) website.
