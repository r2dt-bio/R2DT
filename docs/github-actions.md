# GitHub Actions

R2DT automates the building and deployment of [Docker images](./docker.md) using **GitHub Actions**. The project employs several workflows:

1.	Base Image Workflow – Creates a foundational Docker image.
2.	Main Workflow – Builds the final image using the base image.
3.	Tag Base Image Workflow – Manually tags the base image with a specific version.
4.	Tag Release Image Workflow – Automatically tags the main image with release versions.

This automation streamlines development, ensuring consistency, reducing build times, and minimising manual intervention.

## Base image

The base image serves as a starting point for the main Docker image. It includes essential dependencies and configurations that are common across different builds, ensuring consistency and reducing build times for the resulting images.

The base image is created via the workflow defined in [parallel-base-image.yml](https://github.com/r2dt-bio/R2DT/blob/main/.github/workflows/parallel-base-image.yml). This workflow provides a standardised environment for subsequent images. The image layers are built in parallel using [BuildJet](https://buildjet.com) and pushed to Docker Hub for caching.

## Main image

The main workflow for building the resulting Docker image is defined in [main.yml](https://github.com/r2dt-bio/R2DT/blob/main/.github/workflows/main.yml). This workflow is triggered on pushes and pull requests to the `main` and `develop` branches, as well as on tag pushes matching the pattern `v*` (e.g., `v1.0`, `v2.1.3`).

### Key steps

1. **Initial notification**: Sends a Slack notification indicating the start of the Docker image creation process.
2. **Create Docker tag**: Uses the Docker metadata action to generate tags for the Docker image based on the branch or tag being built.
3. **Build and push**: Builds the Docker image using `Dockerfile` and pushes it to the Docker registry. This step supports parallel builds for different platforms (e.g., `linux/amd64`, `linux/arm64`).
4. **Comment on pull request**: If the workflow is triggered by a pull request, it comments on the PR with the Docker image tags and labels.
5. **Final notification**: Sends a Slack notification indicating the completion of the Docker image creation process.

### Parallel builds in Dockerfile

The `Dockerfile` supports parallel builds by specifying multiple stages for different components. Each stage is responsible for building a specific part of the application, for example:

- **tRNAScan-SE**: Installs tRNAScan-SE, a tool for identifying tRNA genes.
- **Traveler**: Installs Traveler, a tool for RNA structure visualization.
- **Ribovore-Infernal-Easel**: Installs Ribovore and Infernal, tools for RNA sequence analysis.

These stages are combined in the final build stage to create a comprehensive Docker image that includes all the necessary tools and dependencies.

### Parallel builds in main workflow

The main workflow supports parallel builds for multiple platforms using Docker Buildx. This allows the resulting Docker image to be compatible with different architectures, enhancing its usability across various environments.

By automating these processes with GitHub Actions, the project ensures consistent and efficient Docker image builds, reducing manual effort and potential errors.

## Tagging workflows

### Tag base image

The workflow defined in [tag-base-image.yml](https://github.com/r2dt-bio/R2DT/blob/main/.github/workflows/tag-base-image.yml) allows repository maintainers to manually tag a specific version of the base image. This workflow is triggered manually and requires two inputs:
- **Source Tag**: The existing tag of the image to be retagged
- **Destination Tag**: The new tag to apply to the image

This workflow is typically used after a pull request that modifies the base image has been merged, to create a stable version tag for the base image.

### Tag release image

The workflow defined in [tag-release-image.yml](https://github.com/r2dt-bio/R2DT/blob/main/.github/workflows/tag-release-image.yml) automatically tags the main R2DT Docker image when a new release is published on GitHub. This workflow:

1. Is triggered when a new release is published
2. Extracts the Git tag that the release was created from
3. Checks if a Docker image with that tag exists
4. Tags the corresponding Docker image with both the full version number (e.g., `1.2.3`) and the major.minor version (e.g., `1.2`)
5. Sends a Slack notification about the new tagged release

If no Docker image exists with the tag corresponding to the Git release tag, the workflow will fail with an error message. This ensures that stable version tags are only created for Docker images that were specifically built for the release.

This workflow ensures that stable version tags are available for the main R2DT image, making it easier to pin specific versions as dependencies.
