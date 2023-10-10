# Base and R2DT docker images

## Overview

This project uses combination of two Docker images: `base_image`, which contains all the non-Python dependencies of R2DT, and the R2DT image itself, which includes the necessary Python packages.

## Base Image

Base image's build process is detailed in `base_image/Dockerfile`. In essence, it uses GCC to compile C/C++ code and installs some other (e.g., Perl-related) dependencies to the `/install` folder.

Afterward, the image is uploaded to Docker Hub under the name `rnacentral/r2dt-base`.

It's anticipated that the base image will rarely need rebuilding. Such an action would only be needed when an R2DT dependency gets updated and subsequently needs to be integrated into the R2DT build.

## R2DT Image

This image is constructed using the `Dockerfile` located in the repository's root folder, inheriting from `r2dt-base`.

This approach allows to keep the R2DT compact, containing only the Python and R2DT dependencies, and excluding massive build tools used in the base image. Moreover, most of the changes in the Python code can be done using the prebuilt dependencies and doesn't require any recompilation of C/C++ code.

## Upgrade process

When there's a need to update the base image, one needs to:

* Create a branch off of `develop` and make changes to the `base_image/Dockerfile` as necessary.
* Build the base image locally and test it. Specifically, this requires updating the `FROM` directive in the R2DT's image Dockerfile to reference the locally built base image. You may use `just bbuild` command to build the base image locally.
* Once the base image is ready, submit a pull request. This will trigger image build via GitHub Actions. Resulting base image will have a tag after your branch name (e.g., `rnacentral/r2dt-base:my-feature-branch`).
* If you plan to use the newly built image in subsequent R2DT development, reach out to the repo maintainers and ask them to tag your image as new release version using corresponding [workflow](https://github.com/RNAcentral/R2DT/actions/workflows/tag-base-image.yml) (only repo maintainers have permissions run it) and update the `FROM` directive in the R2DT's image Dockerfile to reference the newly tagged base image.


_Note_: when changes to both base and R2DT images are required, it's recommended to put them in separate PRs, to allow the base image to be built with the new version on Dockerhub, so that the R2DT image can be built with the new base image.
