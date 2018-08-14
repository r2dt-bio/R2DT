
# Auto Traveler

This is a tool for automatic generation of RNA secondary structure in standard
[CRW](http://www.rna.ccbb.utexas.edu) layouts using

* [Traveler](https://github.com/davidhoksza/traveler),
* [Ribotyper](https://github.com/nawrockie/ribotyper-v1), and
* [Infernal](http://eddylab.org/infernal/).

## Installation and usage

```
docker build -t auto-traveler .
docker run -it -v `pwd`:/rna/auto-traveler

python auto-traveler/auto-traveler.py examples/examples.fasta auto-traveler/test-output

perl ribotyper.pl -i data/cms/all.modelinfo.txt -f data/pdb.fasta data/test-auto
```

## Acknowledgements

- [David Hoksza](https://github.com/davidhoksza)
- [Eric Nawrocki](https://github.com/nawrockie)
- [Robin Gutell](http://www.rna.ccbb.utexas.edu)
- [RNAStructure package](http://rna.urmc.rochester.edu/RNAstructure.html)

Secondary structure information was downloaded from the [CRW](http://www.rna.ccbb.utexas.edu) website.
