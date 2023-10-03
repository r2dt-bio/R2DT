# Use values from .env automatically, if present
set dotenv-load

alias r := run

# platform := "--platform=linux/amd64"
platform := ""
image := "rnacentral/r2dt"
base_image := "rnacentral/r2dt-base"


# Default recipe to display help information
default:
  @just --list

# Prepare and activate Python virtual environment
venv:
    python3 -m venv .venv
    pip3 install -r requirements.txt
    source .venv/bin/activate


# Run shell in docker
run:
  docker run {{platform}} -v $(pwd):/rna/r2dt -v ./1.4:/rna/r2dt/data/cms -it --rm {{image}}

# Build R2DT Docker image
build:
  docker build -t rnacentral/r2dt .

# Start a development docs server
docs:
  sphinx-autobuild docs docs/_build/html

# Check links in docs
check-links:
  (cd docs && make linkcheck)
