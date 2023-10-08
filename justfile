# Use values from .env automatically, if present
set dotenv-load

alias r := run
alias b  := build
alias bb := bbuild

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

# Run all tests in docker
test-all:
  docker run {{platform}} --rm -it -v ./:/rna/r2dt/ -v ./1.4/:/rna/r2dt/data/cms {{image}} bash -c "cd r2dt && ./r2dt.py test"
# Run specific test in docker
test TEST:
  docker run {{platform}} --rm -it -v ./:/rna/r2dt/ -v ./1.4/:/rna/r2dt/data/cms {{image}} bash -c "cd r2dt && ./r2dt.py test Test{{TEST}}"

# Build R2DT Docker image
build:
  docker buildx build {{platform}} -t {{image}}  .

# Build base image
bbuild:
  docker buildx build {{platform}} -t {{base_image}} base_image

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
