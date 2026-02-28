# Docker images

## Overview

This project uses a combination of two Docker images: `base_image` and the R2DT image. The `base_image` contains all non-Python dependencies of R2DT, while the R2DT image includes the necessary Python packages. The images are built and pushed to Docker Hub using [GitHub Actions](./github-actions.md).

## Base image

The build process for the base image is detailed in [base_image/Dockerfile](https://github.com/r2dt-bio/R2DT/blob/main/base_image/Dockerfile). It uses GCC to compile C/C++ code and installs additional dependencies, such as Perl-related tools.

Once built, the base image is uploaded to Docker Hub under the name `rnacentral/r2dt-base`.

Rebuilding the base image is rare and is necessary only when an R2DT dependency is updated and needs integration into the R2DT build.

View [rnacentral/r2dt-base](https://hub.docker.com/r/rnacentral/r2dt-base) on Docker Hub &rarr;

## Main image

The main R2DT image is constructed using the [Dockerfile](https://github.com/r2dt-bio/R2DT/blob/main/Dockerfile) located in the repository's root folder. It inherits from `r2dt-base`.

This approach ensures that the R2DT image remains compact, containing only the Python and R2DT dependencies, excluding the extensive build tools used in the base image. Most changes in the Python code can be made using the prebuilt dependencies and do not require recompilation of C/C++ code.

When a Git tag (e.g., `v2.1.3`) is pushed to the repository, the [main.yml](https://github.com/r2dt-bio/R2DT/blob/main/.github/workflows/main.yml) workflow automatically builds and tags the Docker image with the corresponding version tags (e.g., `2.1.3` and `2.1`). This ensures that stable version tags are available for the main R2DT image, making it easier to pin specific versions of R2DT as dependencies.

### Bundled covariance model libraries

Starting with version 2.2, the Rfam and CRW covariance models are bundled directly in the Docker image. The uncompressed `all.cm` files are too large for git (~275 MB for Rfam, ~508 MB for CRW), so they are stored as compressed `all.cm.tar.gz` archives in the repository. The Dockerfile extracts and indexes them at build time:

```dockerfile
ADD data/rfam/cms/all.cm.tar.gz data/rfam/cms/
RUN cmfetch --index data/rfam/cms/all.cm
```

When updating Rfam or CRW models, regenerate the tarballs before building:

```bash
python3 r2dt.py compress-rfam-crw
git add data/rfam/cms/all.cm.tar.gz data/crw/all.cm.tar.gz
git commit -m "Update CM library archives"
```

This eliminates the need for users to download and mount a precomputed library separately.

View [rnacentral/r2dt](https://hub.docker.com/r/rnacentral/r2dt) on Docker Hub &rarr;

## Upgrade process

When updating the base image, follow these steps:

1. Create a branch from `develop` and make necessary changes to `base_image/Dockerfile`.
2. Build the base image locally and test it. Use the `just bbuild` command to build the base image locally and `just full-build` to build both the base image and the main image locally.
3. Submit a pull request, triggering image build via GitHub Actions. The resulting base image will be tagged with your pull request ID (e.g., `rnacentral/r2dt-base:pr-109`).
4. If you plan to use the newly built image in subsequent R2DT development, [notify](https://github.com/r2dt-bio/r2dt/issues/new) the repository maintainers. Ask them to tag your image as a new release version using the corresponding [workflow](https://github.com/r2dt-bio/R2DT/actions/workflows/tag-base-image.yml). Only repository maintainers have permission to run this workflow. Update the `FROM` directive in the R2DT's image Dockerfile to reference the newly tagged base image.

_Note_: If you would like to build the R2DT image against a custom version of the base image, you can do so by using `just tag-build <TAG>` command:

```bash
# build r2dt image against the base image tagged as pr-109
just tag-build pr-109

# build r2dt image against the base image tagged as latest
just tag-build latest
```

_Note_: When changes to both base and R2DT images are required, it is recommended to submit separate pull requests. This allows the base image to be built with the new version on Dockerhub before the R2DT image is built with the updated base image.


## Building for non-default platforms/architectures

Should you need to build any of the images for a different platform (e.g., `arm64`), you can do so by specifying the `platform=PLATFORM_NAME` flag when building the image. For example:

```bash
# build base image for linux/arm64
just platform=linux/arm64 bbuild

# build r2dt image for linux/amd64
just platform=linux/amd64 build

# build both images for darwin/arm64
just platform=darwin/arm64 full-build
```
