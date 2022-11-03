# R2DT documentation

[![RTD badge](https://readthedocs.org/projects/r2dt/badge/?version=latest)](https://r2dt.readthedocs.io/en/latest/?badge=latest)

## What is R2DT?

The [R2DT software](https://github.com/RNAcentral/R2DT) automatically generates [RNA secondary structure](https://en.wikipedia.org/wiki/Nucleic_acid_secondary_structure) diagrams in consistent, reproducible and recongnisable layouts using a library of templates representing a wide range of RNAs:

 * 3D-structure based SSU and LSU rRNA from [RiboVision](http://apollo.chemistry.gatech.edu/RiboVision/)
 * 5S and SSU rRNA from [CRW](http://www.rna.ccbb.utexas.edu)
 * tRNA from [GtRNAdb](http://gtrnadb.ucsc.edu)
 * RNAse P from [Ribonuclease P Database](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC148169/)
 * RNA families from [Rfam](https://rfam.org)

![Method overview](./method-overview.png)

Read the [R2DT paper](https://www.nature.com/articles/s41467-021-23555-5) for a detailed method description.

## Examples

![Examples](./r2dt-examples.png)

## Getting started

R2DT can be used in a number of ways:

* [Web application](https://rnacentral.org/r2dt) hosted by RNAcentral
* [API](https://www.ebi.ac.uk/Tools/common/tools/help/index.html?tool=r2dt) powered by EMBL-EBI Web Services (see [](api.md))
* As a command line tool with [Docker](https://www.docker.com), [Podman](https://podman.io), [Singularity](https://sylabs.io/docs/), or as bare metal installation

## Who uses R2DT?

* RNAcentral uses R2DT to visualise [>20 million RNA secondary structures](https://rnacentral.org/search?q=has_secondary_structure:%22True%22)
* Rfam displays R2DT diagrams in [sequence similarity search](https://rfam.org/search#tabview=tab1)
* PDBe uses R2DT to enable interactive navigation between sequence, 2D and 3D structure (for example, [1S72](https://www.ebi.ac.uk/pdbe/entry/pdb/1s72/RNA/1))
* [FlyBase](http://flybase.org/reports/FBgn0053537#gene_model_products) and [SGD](https://www.yeastgenome.org/locus/S000006550/sequence) show R2DT diagrams for RNA genes

```{eval-rst}
.. toctree::
   :maxdepth: 2
   :caption: Table of Contents

   installation
   usage
   api
   templates
   widget
   docs
   team
```

## Citation

If you use R2DT in your work, please consider citing the following paper:

> **R2DT is a framework for predicting and visualising RNA secondary structure using templates**
> [Nature Communications](https://www.nature.com/articles/s41467-021-23555-5)

## License

R2DT is available under the [Apache 2.0 license](https://github.com/RNAcentral/R2DT/blob/master/LICENSE).

## Get in touch

If you have any questions or feedback, feel free to [submit a GitHub issue](https://github.com/RNAcentral/r2dt/issues) or contact the [RNAcentral help desk](https://rnacentral.org/contact).
