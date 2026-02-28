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
port := "8000"

# Default recipe to display help information
default:
    @just --list

# Prepare and activate Python virtual environment
venv:
    python3 -m venv .venv
    pip3 install -r requirements.txt
    source .venv/bin/activate

# Download precomputed data from GitHub
download data_version="2.1":
    curl -O -L https://github.com/r2dt-bio/R2DT/releases/download/v{{ data_version }}/cms.tar.gz
    tar -xzf cms.tar.gz

# Run shell in docker
run tag=default_tag:
    docker run {{ platform_arg }} -v $(pwd):/rna/r2dt -it --rm {{ image }}:{{tag}}

# Run without mounting the current directory
run_no_mount tag=default_tag:
    docker run {{ platform_arg }} -it --rm {{ image }}:{{tag}}

# Run all tests in docker
test-all:
    docker run {{ platform_arg }} --rm -it -v ./:/rna/r2dt/ {{ image }} bash -c "./r2dt.py test"

# Run specific test in docker
test TEST:
    docker run {{ platform_arg }} --rm -it -v ./:/rna/r2dt/ {{ image }} bash -c "./r2dt.py test Test{{ TEST }}"

# Build R2DT Docker image
build base_version="" tag=default_tag:
    #!/usr/bin/env bash
    set -euxo pipefail
    [[ "{{base_version}}" == "" ]] && build_arg="" || build_arg="--build-arg BASE_IMAGE_VERSION={{base_version}}"
    docker buildx build --builder default --load {{ platform_arg }} $build_arg -t {{ image }}:{{tag}}  .

# Shortcut to build the R2DT Docker image against custom base image
tag-build tag: (build tag tag)

# Build base image locally
bbuild:
    docker buildx build --builder default --load {{ platform_arg }} -t {{ base_image }} base_image

# Build base and then the r2dt images locally
full-build: bbuild (tag-build "latest")

# Regenerate viral doc images from example inputs
docs-images tag=default_tag:
    #!/usr/bin/env bash
    set -euo pipefail
    img=docs/images
    tmp=temp/docs-images
    run="docker run {{ platform_arg }} --rm -v $(pwd):/rna/r2dt -w /rna/r2dt {{ image }}:{{tag}}"

    echo "==> SARS-CoV-2 (FASTA → viral-annotate + stitch)"
    $run python3 r2dt.py viral-annotate \
        examples/viral/coronavirus.fasta $tmp/coronavirus/ --quiet
    cp $tmp/coronavirus/rfam/*_26-299_*-RF03120*.colored.svg   $img/RF03120_26-299.svg
    cp $tmp/coronavirus/rfam/*_13469-13546_*-RF00507*.colored.svg $img/RF00507_13469-13546.svg
    cp $tmp/coronavirus/rfam/*_29536-29870_*-RF03125*.colored.svg $img/RF03125_29536-29870.svg
    $run python3 r2dt.py stitch \
        $tmp/coronavirus/rfam/*.colored.svg \
        -o $img/coronavirus-stitched.svg \
        --sort --normalize-font-size \
        --captions "5′ UTR" --captions "FSE" --captions "3′ UTR"
    $run python3 r2dt.py stitch \
        $tmp/coronavirus/rfam/*.colored.svg \
        -o $img/coronavirus-stitched-color.svg \
        --sort --color --normalize-font-size \
        --captions "5′ UTR" --captions "FSE" --captions "3′ UTR"

    echo "==> HCV (FASTA → viral-annotate + stitch)"
    $run python3 r2dt.py viral-annotate \
        examples/viral/hcv.fasta $tmp/hcv-rfam/ --quiet
    $run python3 r2dt.py stitch \
        $tmp/hcv-rfam/rfam/*.colored.svg \
        -o $img/hcv-stitched.svg \
        --sort --normalize-font-size

    echo "==> Dengue 2 (FASTA → viral-annotate + stitch)"
    $run python3 r2dt.py viral-annotate \
        examples/viral/dengue2.fasta $tmp/dengue/ --quiet
    $run python3 r2dt.py stitch \
        $tmp/dengue/rfam/*.colored.svg \
        -o $img/dengue2-stitched.svg \
        --sort --normalize-font-size

    echo "==> HCV (Stockholm → stockholm)"
    $run python3 r2dt.py stockholm \
        examples/hcv-alignment.stk $tmp/hcv-stockholm/ --quiet
    cp $tmp/hcv-stockholm/stitched.svg           $img/hcv-stockholm-stitched.svg
    cp $tmp/hcv-stockholm/stitched-thumbnail.svg $img/hcv-stockholm-thumbnail.svg

    rm -rf $tmp
    echo "✓ All viral doc images regenerated in $img/"

# Start a development docs server
docs:
    docker run {{platform}} -p {{port}}:{{port}} -v $(pwd):/rna/r2dt -it --rm {{image}} sphinx-autobuild --host 0.0.0.0 --port {{port}} docs docs/_build/html

# Check links in docs
check-links:
    (cd docs && make linkcheck)

# Delete test results
clean:
    -rm -rf tests/results
    -rm -rf tests/html/*.html
