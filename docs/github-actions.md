# GitHub Actions for Docker Image Builds

This document explains how GitHub Actions are used in this project to automate the building and pushing of Docker images. The project uses two main workflows: one for building the base image and another for building the resulting image.

## Building the Base Image

The base image is built using a separate workflow defined in `.github/workflows/base-image.yml`. This workflow is responsible for creating a foundational Docker image that other images can build upon. The base image is defined in `base_image/Dockerfile`.

### Key Steps in Base Image Workflow

1. **Checkout Code**: The workflow checks out the code from the repository to access the Dockerfile and any necessary scripts.
2. **Set Up Docker**: It sets up Docker Buildx, which is a Docker CLI plugin for extended build capabilities with BuildKit.
3. **Build and Push**: The base image is built using the instructions in `base_image/Dockerfile` and pushed to a Docker registry.

### Usage of Base Image

The base image serves as a starting point for the main Docker image. It includes essential dependencies and configurations that are common across different builds, ensuring consistency and reducing build times for the resulting images.

## Building the Resulting Image

The main workflow for building the resulting Docker image is defined in `.github/workflows/main.yml`. This workflow is triggered on pushes and pull requests to the `main` and `develop` branches.

### Key Steps in Main Workflow

1. **Initial Notification**: Sends a Slack notification indicating the start of the Docker image creation process.
2. **Create Docker Tag**: Uses the Docker metadata action to generate tags for the Docker image based on the branch or tag being built.
3. **Build and Push**: Builds the Docker image using `Dockerfile` and pushes it to the Docker registry. This step supports parallel builds for different platforms (e.g., `linux/amd64`, `linux/arm64`).
4. **Comment on Pull Request**: If the workflow is triggered by a pull request, it comments on the PR with the Docker image tags and labels.
5. **Final Notification**: Sends a Slack notification indicating the completion of the Docker image creation process.

### Parallel Builds in Dockerfile

The `Dockerfile` supports parallel builds by specifying multiple stages for different components. Each stage is responsible for building a specific part of the application, such as:

- **R-scape**: Installs R-scape, a tool for RNA sequence analysis.
- **tRNAScan-SE**: Installs tRNAScan-SE, a tool for identifying tRNA genes.
- **Bio-Easel**: Installs Bio-Easel, a library for biological sequence analysis.
- **Traveler**: Installs Traveler, a tool for RNA structure visualization.
- **Scripts**: Installs various scripts for RNA analysis.
- **Ribovore-Infernal-Easel**: Installs Ribovore and Infernal, tools for RNA sequence analysis.
- **RNArtist**: Installs RNArtist, a tool for RNA structure visualization.

These stages are combined in the final build stage to create a comprehensive Docker image that includes all necessary tools and dependencies.

### Parallel Builds in Main Workflow

The main workflow supports parallel builds for multiple platforms using Docker Buildx. This allows the resulting Docker image to be compatible with different architectures, enhancing its usability across various environments.

By automating these processes with GitHub Actions, the project ensures consistent and efficient Docker image builds, reducing manual effort and potential errors.
