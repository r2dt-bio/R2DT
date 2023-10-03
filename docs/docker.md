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

* Create a branch from `develop` and make changes to the `base_image/Dockerfile` as necessary.
* Build the base image locally and test it. Specifically, this requires updating the `FROM` directive in the R2DT's image Dockerfile to reference the locally built base image.
* Once the base image is ready, push it to Docker Hub as `rnacentral/r2dt-base:<version>`. Ideally, retain the prior version's number, appending a distinct build number (e.g., `1.4.0-2`), unless changes to the R2DT are also required (for instance, when one of the dependencies changes its API).

_Note_: the base image is not tagged as `latest` to avoid accidental use of the latest version of the base image by R2DT.

_Note 2_: when changes to both base and R2DT images are required, it's recommended to put them in separate PRs, to allow the base image to be built with the new version on Dockerhub, so that the R2DT image can be built with the new base image.
