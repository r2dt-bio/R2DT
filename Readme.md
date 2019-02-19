
# Auto Traveler

This is a tool for automatic generation of RNA secondary structure in standard
[CRW](http://www.rna.ccbb.utexas.edu) layouts.

## Method overview

1. Generate a library of **covariance models** using [CRW bpseq files](http://www.rna.icmb.utexas.edu/DAT/3C/Structure/index.php)
and [Infernal](http://eddylab.org/infernal/). For best results, use _pseudoknot-free_ CRW bpseq files
or remove pseudoknots using [RemovePseudoknots](https://rna.urmc.rochester.edu/Text/RemovePseudoknots.html) from the RNAStructure package.
1. **Select the best matching covariance model** for each input sequence
using [Ribotyper](https://github.com/nawrockie/ribotyper-v1)
1. **Fold** input sequence into a secondary structure compatible with the template
using the top scoring covariance model.
1. **Generate secondary structure diagrams** using [Traveler](https://github.com/davidhoksza/traveler) and CRW layouts.

## Installation

Download a precomputed library of covariance models:
https://www.dropbox.com/s/q5l0s1nj5h4y6e4/cms.tar.gz?dl=0

Uncompress it and place the folder in `auto-traveler/data/cms`.

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
python auto-traveler.py examples/examples.fasta example-output

# place your fasta file in a folder that is mounted in the container:
python auto-traveler.py /path/to/input.fasta /path/to/output-folder
```

Additional commands:

```
cd auto-traveler

# classify example sequences using Ribotyper
perl ribotyper.pl -i data/cms/all.modelinfo.txt -f examples/pdb.fasta example-output

# to generate covariance models:
python utils/generate_cm_library.py

# generate png version of images using Image Magick
for d in *.colored.svg; do convert $d ${d}.png; done
```

## Acknowledgements

- [David Hoksza](https://github.com/davidhoksza)
- [Eric Nawrocki](https://github.com/nawrockie)
- [Robin Gutell lab](http://www.rna.ccbb.utexas.edu)
- [Anton S. Petrov](https://scholar.google.com/citations?user=V9KP2IkAAAAJ&hl=en)
- [David Mathews lab](http://rna.urmc.rochester.edu/RNAstructure.html)

Secondary structure information was downloaded from the [CRW](http://www.rna.ccbb.utexas.edu) website.
