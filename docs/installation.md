# Installation

![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/rnacentral/r2dt)

* Download the R2DT image from [Docker Hub](https://hub.docker.com/r/rnacentral/r2dt) and run it with Docker or Singularity.

    **Docker**

    ```
    docker pull rnacentral/r2dt
    docker run --entrypoint r2dt.py rnacentral/r2dt draw --help
    ```

    **Singularity**

    ```
    singularity build r2dt docker://rnacentral/r2dt
    singularity exec r2dt r2dt.py draw --help
    ```

* :hammer_and_wrench: Development installation:
    ```
    # Get the code
    git clone https://github.com/RNAcentral/R2DT.git
    cd R2DT

    # Build and tag a Docker image
    docker build -t rnacentral/r2dt .
    docker-compose run cli
    ```
    The current directory is mounted inside the container so that all code and data changes are instantly reflected in the container.

* :hammer_and_wrench: Bare metal installation: if running R2DT using containers is not possible, follow instructions in the [Dockerfile](./Dockerfile).

### Initial setup

1. Download a [precomputed data library](https://ftp.ebi.ac.uk/pub/databases/RNAcentral/r2dt/1.3/cms.tar.gz) _(198 MB, last updated Oct 14, 2022)_ and uncompress it.

2. Enter an interactive Docker terminal session:

    ```
    docker run -it -v <path_to_cms>:/rna/r2dt/data/cms -v `pwd`:/rna/r2dt/temp rnacentral/r2dt
    ```

    - `-it` - start an interactive session
    - `-v <path_to_cms>:/rna/r2dt/data/cms` - mount the precomputed data library folder `<path_to_cms>` as `/rna/r2dt/data/cms` inside the container. :warning: Note that `<path_to_cms>` should be a full path.
    - make the current working directory available inside the container as `/rna/r2dt/temp`:

        ```
        -v `pwd`:/rna/r2dt/temp
        ```

Any file placed in `/rna/r2dt/temp` within the container will be available on the host machine after the Docker container exits.
