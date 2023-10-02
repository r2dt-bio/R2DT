# Installation

1. Pull an image from [Docker Hub](https://hub.docker.com/r/rnacentral/r2dt):
    ```bash
    docker pull rnacentral/r2dt
    ```

    Alternatively, build a Docker image locally:

    ```bash
    # Get the code
    git clone https://github.com/RNAcentral/R2DT.git
    cd R2DT
    docker build -t rnacentral/r2dt .
    ```

    Or build a Singularity image:
    ```bash
    singularity build r2dt docker://rnacentral/r2dt
    ```

2. Setup a precomputed data library _(198 MB, last updated Aug 21, 2023)_:
    ```bash
    curl -O https://ftp.ebi.ac.uk/pub/databases/RNAcentral/r2dt/1.4/cms.tar.gz
    tar -xzf cms.tar.gz
    export R2DT_LIBRARY=<path to precomputed library>
    ```

3. Verify that the installation worked:
    ```bash
    docker run --entrypoint r2dt.py rnacentral/r2dt draw --help
    ```

    Or in Singularity:
    ```bash
    singularity exec r2dt r2dt.py draw --help
    ```

4. Enter an interactive terminal session:
    ```bash
    docker run -it -v $R2DT_LIBRARY:/rna/r2dt/data/cms -v `pwd`:/rna/r2dt/temp rnacentral/r2dt
    ```

    - `-it` - start an interactive session
    - `-v $R2DT_LIBRARY:/rna/r2dt/data/cms` - mount the precomputed data library folder as `/rna/r2dt/data/cms` inside the container. ⚠️ Note that `<path_to_cms>` should be a full path.
    - make the current working directory available inside the container as `/rna/r2dt/temp`:
        ```bash
        -v `pwd`:/rna/r2dt/temp
        ```

    Any file placed in `/rna/r2dt/temp` within the container will be available on the host machine after the Docker container exits. The current directory is mounted inside the container so that all code and data changes are instantly reflected in the container.

## Setup a development environment

Set up a precomputed library as described above and run the container using [just](https://just.systems):

```bash

# see what commands are available
just

# run the container
just run
```

## Manual installation

If it is not possible to use containers, follow instructions in the [base Dockerfile](https://github.com/RNAcentral/R2DT/blob/main/base_image/Dockerfile) and [main Dockerfile](https://github.com/RNAcentral/R2DT/blob/main/Dockerfile) to install all the requirements manually.
