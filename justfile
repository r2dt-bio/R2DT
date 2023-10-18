# Use values from .env automatically, if present

set dotenv-load := true

alias r := run
alias b := build
alias bb := bbuild

platform := ""
platform_arg := if platform == "" { "" } else { "--platform=" + platform }
base_image := "rnacentral/r2dt-base"
image := "rnacentral/r2dt"
default_tag := "latest"
data_version := "1.4"
data_dir := "./" + data_version

# Default recipe to display help information
default:
    @just --list

# Prepare and activate Python virtual environment
venv:
    python3 -m venv .venv
    pip3 install -r requirements.txt
    source .venv/bin/activate

# Download data from RNAcentral
download:
    curl -O https://ftp.ebi.ac.uk/pub/databases/RNAcentral/r2dt/{{ data_version }}/cms.tar.gz
    tar -xzf cms.tar.gz

# Run shell in docker
run:
    docker run {{ platform_arg }} -v $(pwd):/rna/r2dt -v {{ data_dir }}:/rna/r2dt/data/cms -it --rm {{ image }}

# Run all tests in docker
test-all:
    docker run {{ platform_arg }} --rm -it -v ./:/rna/r2dt/ -v {{ data_dir }}:/rna/r2dt/data/cms {{ image }} bash -c "./r2dt.py test"

# Run specific test in docker
test TEST:
    docker run {{ platform_arg }} --rm -it -v ./:/rna/r2dt/ -v {{ data_dir }}:/rna/r2dt/data/cms {{ image }} bash -c "./r2dt.py test Test{{ TEST }}"

# Build R2DT Docker image
build base_version="" tag=default_tag:
    #!/usr/bin/env bash
    set -euxo pipefail
    [[ "{{base_version}}" == "" ]] && build_arg="" || build_arg="--build-arg BASE_IMAGE_VERSION={{base_version}}"
    docker buildx build {{ platform_arg }} $build_arg -t {{ image }}:{{tag}}  .

# Shortcut to build the R2DT Docker image against custom base image
tag-build tag: (build tag tag)

# Build base image locally
bbuild:
    docker buildx build {{ platform_arg }} -t {{ base_image }} base_image

# Build base and then the r2dt images locally
full-build: bbuild (tag-build "latest")

# Start a development docs server
docs:
    sphinx-autobuild docs docs/_build/html

# Check links in docs
check-links:
    (cd docs && make linkcheck)

# Delete test rusults
clean:
    -rm -rf tests/results
    -rm tests/*.html
