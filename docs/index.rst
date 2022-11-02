Welcome to R2DT documentation!
================================

|docs|

What is R2DT?
-------------

The R2DT software (RNA 2D Templates) automatically generates `RNA secondary structure <https://en.wikipedia.org/wiki/Nucleic_acid_secondary_structure>`_ diagrams in standard layouts using a template library representing a wide range of RNAs:

 * 5S and SSU rRNA from `CRW <http://www.rna.ccbb.utexas.edu>`_
 * 3D-structure based SSU and LSU rRNA from `RiboVision <http://apollo.chemistry.gatech.edu/RiboVision/>`_
 * tRNA from `GtRNAdb <http://gtrnadb.ucsc.edu>`_
 * RNAse P from `Ribonuclease P Database <https://academic.oup.com/nar/article/26/1/351/2379438>`_
 * RNA families from `Rfam <https://rfam.org>`_

Read the `R2DT paper <https://www.nature.com/articles/s41467-021-23555-5>`_ in Nature Communications for description of the method.

Who uses R2DT?
--------------

* RNAcentral uses R2DT to visualise `>20 million RNA secondary structures <https://rnacentral.org/search?q=has_secondary_structure:%22True%22>`_
* `Rfam <https://rfam.org/search#tabview=tab1>`_ displays R2DT diagrams in sequence similarity search results
* `PDBe <https://www.ebi.ac.uk/pdbe/entry/pdb/1s72/RNA/1>`_ uses R2DT to enable interactive navigation between sequence, 2D and 3D structure
* `FlyBase <http://flybase.org/reports/FBgn0053537#gene_model_products>`_ and `SGD <https://www.yeastgenome.org/locus/S000006550/sequence>`_ show R2DT diagrams for RNA genes

.. toctree::
   :maxdepth: 2
   :caption: Table of Contents

   installation
   usage
   templates
   docs
   team

Citation
--------

If you use R2DT in your work, please consider citing the following paper:

**R2DT is a framework for predicting and visualising RNA secondary structure using templates**
*Nature Communications*
https://www.nature.com/articles/s41467-021-23555-5

License
-------
R2DT is available under the `Apache 2.0 license <https://github.com/RNAcentral/R2DT/blob/master/LICENSE>`_.

Get in touch
------------

If you have any questions or feedback, feel free to `submit a GitHub issue <https://github.com/RNAcentral/r2dt/issues>`_ or `contact the RNAcentral help desk <https://rnacentral.org/contact-us>`_.

.. |docs| image:: https://readthedocs.org/projects/r2dt/badge/?version=latest
    :alt: Documentation Status
    :scale: 100%
    :target: https://r2dt.readthedocs.io/en/latest/?badge=latest