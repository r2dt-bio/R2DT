# Base and R2DT Docker Images

## Overview

This project uses a combination of two Docker images: `base_image` and the R2DT image. The `base_image` contains all non-Python dependencies of R2DT, while the R2DT image includes the necessary Python packages.

## Base Image

The build process for the base image is detailed in [base_image/Dockerfile](https://github.com/RNAcentral/R2DT/blob/main/base_image/Dockerfile). It uses GCC to compile C/C++ code and installs additional dependencies, such as Perl-related tools.

Once built, the base image is uploaded to Docker Hub under the name `rnacentral/r2dt-base`.

Rebuilding the base image is rare and is necessary only when an R2DT dependency is updated and needs integration into the R2DT build.

View [rnacentral/r2dt-base](https://hub.docker.com/r/rnacentral/r2dt-base) on Docker Hub &rarr;

## R2DT Image

The R2DT image is constructed using the [Dockerfile](https://github.com/RNAcentral/R2DT/blob/main/Dockerfile) located in the repository's root folder. It inherits from `r2dt-base`.

This approach ensures that the R2DT image remains compact, containing only the Python and R2DT dependencies, excluding the extensive build tools used in the base image. Most changes in the Python code can be made using the prebuilt dependencies and do not require recompilation of C/C++ code.

View [rnacentral/r2dt](https://hub.docker.com/r/rnacentral/r2dt) on Docker Hub &rarr;

## Upgrade Process

When updating the base image, follow these steps:

1. Create a branch from `develop` and make necessary changes to `base_image/Dockerfile`.
2. Build the base image locally and test it. Update the `FROM` directive in the R2DT's image Dockerfile to reference the locally built base image. Use the `just bbuild` command to build the base image locally.
3. Submit a pull request, triggering image build via GitHub Actions. The resulting base image will be tagged with your pull request ID (e.g., `rnacentral/r2dt-base:pr-109`).
4. If you plan to use the newly built image in subsequent R2DT development, [notify](https://github.com/rnacentral/r2dt/issues/new) the repository maintainers. Ask them to tag your image as a new release version using the corresponding [workflow](https://github.com/RNAcentral/R2DT/actions/workflows/tag-base-image.yml). Only repository maintainers have permission to run this workflow. Update the `FROM` directive in the R2DT's image Dockerfile to reference the newly tagged base image.

_Note_: When changes to both base and R2DT images are required, it is recommended to submit separate pull requests. This allows the base image to be built with the new version on Dockerhub before the R2DT image is built with the updated base image.
