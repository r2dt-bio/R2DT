# GitHub Actions

R2DT automates the building and deployment of [Docker images](./docker.md) using **GitHub Actions**. The project employs two primary workflows:

1.	Base Image Workflow – Creates a foundational Docker image.
2.	Main Workflow – Builds the final image using the base image.

This automation streamlines development, ensuring consistency, reducing build times, and minimising manual intervention.

## Base image

The base image serves as a starting point for the main Docker image. It includes essential dependencies and configurations that are common across different builds, ensuring consistency and reducing build times for the resulting images.

The base image is created via the workflow defined in [parallel-base-image.yml](https://github.com/r2dt-bio/R2DT/blob/main/.github/workflows/parallel-base-image.yml). This workflow provides a standardised environment for subsequent images. The image layers are built in parallel using [BuildJet](https://buildjet.com) and pushed to Docker Hub for caching.

## Main image

The main workflow for building the resulting Docker image is defined in [main.yml](https://github.com/r2dt-bio/R2DT/blob/main/.github/workflows/main.yml). This workflow is triggered on pushes and pull requests to the `main` and `develop` branches.

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
