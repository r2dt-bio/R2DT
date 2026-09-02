# Installation

1. Pull an image from [Docker Hub](https://hub.docker.com/r/rnacentral/r2dt):
    ```bash
    docker pull rnacentral/r2dt
    ```

    Alternatively, build a Docker image locally (requires [just](https://just.systems)):

    ```bash
    # Get the code
    git clone https://github.com/r2dt-bio/R2DT.git
    cd R2DT
    just build
    ```

    Or build a Singularity image:
    ```bash
    singularity build r2dt docker://rnacentral/r2dt
    ```

2. Enter an interactive terminal session:
    ```bash
    docker run -it -v `pwd`:/rna/r2dt/temp rnacentral/r2dt
    ```

    - `-it` - start an interactive session
    - make the current working directory available inside the container as `/rna/r2dt/temp`:
        ```bash
        -v `pwd`:/rna/r2dt/temp
        ```

    Any file placed in `/rna/r2dt/temp` within the container will be available on the host machine after the Docker container exits. The current directory is mounted inside the container so that all code and data changes are instantly reflected in the container.

:::{note}

Starting with version 2.2, downloading a precomputed library is no longer necessary. For older versions, however, you must download the library manually and mount it inside the container using Docker’s `-v` option. For example:

```bash
curl -O -L https://github.com/r2dt-bio/R2DT/releases/download/v2.0/cms.tar.gz
tar -xzf cms.tar.gz
export R2DT_LIBRARY=<path to precomputed library>
docker run -it -v $R2DT_LIBRARY:/rna/r2dt/data/cms -v `pwd`:/rna/r2dt/temp rnacentral/r2dt
```
:::

## Setup a development environment

To set up a development container, you can use [just](https://just.systems) by running the following commands:

```bash
# Display available commands
just

# Start the development container
just run
```

Alternatively, if you prefer not to use `just`, you can manually execute the commands listed in the `justfile` to achieve the same result.

## Manual installation

If it is not possible to use containers, follow instructions in the [base Dockerfile](https://github.com/r2dt-bio/R2DT/blob/main/base_image/Dockerfile) and [main Dockerfile](https://github.com/r2dt-bio/R2DT/blob/main/Dockerfile) to install all the requirements manually.

R2DT looks for the Traveler `utils` folder (providing `infernal2mapping.py`, `enrich_json.py`, and `json2svg.py`) in the following locations:

1. the `R2DT_TRAVELER_UTILS` environment variable,
2. `/rna/traveler/utils` (the location used in the Docker image),
3. the `utils` folder of the Traveler installation containing the `traveler` executable found in `PATH`.

If the folder cannot be found, R2DT prints a warning and falls back to `traveler --all`, which computes its own template mapping. The layouts may differ from the standard `traveler --draw` output and the enriched JSON/SVG files are not generated, so when installing manually make sure the Traveler `utils` folder is available, for example:

```bash
export R2DT_TRAVELER_UTILS=/path/to/traveler/utils
```

### Minimal installation for drawing with known templates

Running `r2dt.py draw --force_template <model_id>` (including [custom templates](./templates.md) in `local_data`) requires only a subset of the R2DT dependencies:

- [Traveler](https://github.com/cusbg/traveler) (including its `utils` folder, see above)
- [Infernal](http://eddylab.org/infernal/) (also provides the `esl-*` Easel programs)
- [Bio-Easel](https://github.com/nawrockie/Bio-Easel) (provides `esl-alidepair.pl` and the Perl modules used by the jiffy scripts)
- [jiffy-infernal-hmmer-scripts](https://github.com/nawrockie/jiffy-infernal-hmmer-scripts) (provides `ali-pfam-lowercase-rf-gap-columns.pl` and `ali-pfam-sindi2dot-bracket.pl`)
- the Python packages from [requirements-minimal.txt](https://github.com/r2dt-bio/R2DT/blob/main/requirements-minimal.txt):

    ```bash
    pip install -r requirements-minimal.txt
    ```

    These are all pure Python, so they install without a build toolchain on any
    supported Python version. In particular ViennaRNA, Biopython, NumPy, Pillow,
    scikit-image and CairoSVG are **not** needed for this code path — they are
    imported lazily and only when the features that use them are invoked.

tRNAscan-SE, Ribovore, RNAView, RNArtist, R-scape, and RNAstructure are not used by this code path. Note that template selection (running `draw` without `--force_template`) and template-free mode require the full installation.

The other requirements files are:

| File | Purpose |
| --- | --- |
| `requirements-minimal.txt` | drawing with known templates (`draw --force_template`) |
| `requirements.txt` | full runtime, adds ViennaRNA (constraint folding, RNApuzzler) and Biopython (`pdb` base-pair extraction) |
| `requirements-dev.txt` | everything needed to run `r2dt.py test` |
| `requirements-docs.txt` | building the documentation |
