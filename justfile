alias r := run


# platform := "--platform=linux/amd64"
platform := ""
image := "rnacentral/r2dt"
base_image := "rnacentral/r2dt-base"


# Default recipe to display help information
default:
  @just --list

# Run shell in docker
run:
  docker run {{platform}} -v $(pwd):/rna/r2dt -v ./1.4:/rna/r2dt/data/cms -it --rm {{image}}

# Start a development docs server
docs:
  sphinx-autobuild docs docs/_build/html

links:
  (cd docs && make linkcheck)
