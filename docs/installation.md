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
