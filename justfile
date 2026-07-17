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

# 2D+3D viewer gallery: where it is published and what it contains.
# Each entry is "PDB_ID:mode" (mode = auto | templated | templatefree).
viewers_dir := "output/site"
viewers_structures := "9RJA:auto 8SH5:auto 9CFN:auto 9SFQ:auto 9HRF:auto 9E5I:auto 8XZN:auto 8BWT:auto 8EYW:auto"

# Workstream 1 dashboard: base pairs read from the mmCIF's own annotation
# (--basepairs cif), using the FR3D-converted mmCIF inputs from the na-hackathon
# repo. Published under <viewers_dir>/workstream1/ so it sits beside, but
# clearly apart from, the main gallery.
ws1_dir := viewers_dir / "workstream1"
ws1_structures := "8bwt 8vjt 8xzn 9e5i 9hrf 9sfq"
ws1_base := "https://raw.githubusercontent.com/na-hackathon/na-hackathon-2026/main/data/tests/conversion_outputs/fr3d_mmcif_outputs"

# Cloudflare Pages project name. Set CLOUDFLARE_PROJECT in your (gitignored)
# .env or environment -- it is deliberately not hardcoded here.
cloudflare_project := env_var_or_default("CLOUDFLARE_PROJECT", "")

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

# Refresh the vendored layered-bp-notation script from upstream.
# The full project lives in the na-hackathon repo (not pip-installable yet);
# R2DT only vendors the two files needed for --basepairs cif.
# See utils/layered_bp_notation/README.md.
sync-lbn:
    #!/usr/bin/env bash
    set -euo pipefail
    base="https://raw.githubusercontent.com/na-hackathon/na-hackathon-2026/main/workstreams/ws2-prediction-non-Watson-Crick/layered-bp-notation"
    dest="utils/layered_bp_notation"
    for f in common.py standalone_lbn_script.py; do
        curl -fsSL "$base/$f" -o "$dest/$f"
        echo "✓ synced $dest/$f"
    done

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

# Regenerate every 2D+3D viewer and rebuild the gallery index
viewers tag=default_tag:
    #!/usr/bin/env bash
    set -euo pipefail
    run="docker run {{ platform_arg }} --rm -v $(pwd):/rna/r2dt -w /rna/r2dt {{ image }}:{{tag}}"
    work=output/viewers
    site={{ viewers_dir }}
    mkdir -p "$site"

    for entry in {{ viewers_structures }}; do
        id="${entry%%:*}"
        mode="${entry##*:}"
        echo "==> $id (--mode $mode)"
        rm -rf "$work/$id"
        $run python3 r2dt.py pdb_2d_3d "$id" "$work/$id" --mode "$mode" --quiet

        svg="$work/$id/results/results/svg/$id.colored.svg"
        if [[ ! -d "$work/$id/viewer" || ! -f "$svg" ]]; then
            echo "!! $id produced no viewer (no template match in templated mode?)" >&2
            exit 1
        fi
        rm -rf "$site/$id"
        mkdir -p "$site/$id"
        cp "$work/$id"/viewer/* "$site/$id/"
        cp "$svg" "$site/$id/2d.svg"
    done

    $run python3 utils/build_viewers.py "$site"
    echo "✓ Viewer gallery built in $site/ — preview with: just viewers-serve"

# Regenerate the Workstream 1 dashboard (--basepairs cif) under workstream1/
ws1-viewers tag=default_tag:
    #!/usr/bin/env bash
    set -euo pipefail
    run="docker run {{ platform_arg }} --rm -v $(pwd):/rna/r2dt -w /rna/r2dt {{ image }}:{{tag}}"
    work=output/workstream1
    site={{ ws1_dir }}
    mkdir -p "$work/inputs" "$site"

    for id in {{ ws1_structures }}; do
        cif="$work/inputs/${id}_fr3d.cif"
        [[ -f "$cif" ]] || curl -fsSL "{{ ws1_base }}/${id}_fr3d_basepairs.cif" -o "$cif"
        echo "==> $id (--basepairs cif)"
        rm -rf "$work/out/${id}_fr3d"
        # A structure may have no asymmetric-unit pairs (e.g. 8vjt: symmetry
        # contacts only) and produce no viewer; that is fine, skip it.
        $run python3 r2dt.py pdb_2d_3d "$work/inputs/${id}_fr3d.cif" \
            "$work/out/${id}_fr3d" --basepairs cif --mode auto --quiet || true

        svg="$work/out/${id}_fr3d/results/results/svg/${id}_fr3d.colored.svg"
        if [[ -d "$work/out/${id}_fr3d/viewer" && -f "$svg" ]]; then
            rm -rf "$site/${id}_fr3d"
            mkdir -p "$site/${id}_fr3d"
            cp "$work/out/${id}_fr3d"/viewer/* "$site/${id}_fr3d/"
            cp "$svg" "$site/${id}_fr3d/2d.svg"
            echo "   OK -> $site/${id}_fr3d"
        else
            echo "!! ${id}_fr3d: no viewer produced, skipping" >&2
        fi
    done

    $run python3 utils/build_viewers.py "$site"
    echo "✓ Workstream 1 dashboard built in $site/ — deploy with: just viewers-deploy"

# Serve the viewer gallery locally for preview
viewers-serve:
    python3 -m http.server -d {{ viewers_dir }} {{ port }}

# Publish the viewer gallery (incl. workstream1/) to Cloudflare Pages
# (needs wrangler + CLOUDFLARE_API_TOKEN). Deploys all of viewers_dir, so run
# `just viewers` and/or `just ws1-viewers` first to refresh its contents.
viewers-deploy:
    #!/usr/bin/env bash
    set -euo pipefail
    if [[ -z "{{ cloudflare_project }}" ]]; then
        echo "Set CLOUDFLARE_PROJECT (e.g. in .env) to your Cloudflare Pages project name." >&2
        exit 1
    fi
    npx wrangler@latest pages deploy {{ viewers_dir }} --project-name {{ cloudflare_project }}

# Fold a season's standalone CASP dashboard (output/casp<season>-deploy/,
# built by scripts/casp_rank.py -> casp_fetch.py -> casp_batch.py ->
# casp_dashboard.py) into the single-domain gallery at
# <viewers_dir>/casp<season>/. A wrangler deploy uploads its target directory
# as a full replacement snapshot, not a merge -- deploying casp<season>-deploy
# on its own would wipe out the rest of the gallery, so it has to be folded
# in as a subpath of viewers_dir first. Always rm -rf + fresh-copies the
# destination (never merges into whatever was there before), so a stale or
# partial previous sync can't linger and silently ship. season is "15" or
# "16". See docs/pdb-2d-3d-viewer.md for the full deploy writeup.
casp-sync season:
    #!/usr/bin/env bash
    set -euo pipefail
    src="output/casp{{ season }}-deploy"
    dst="{{ viewers_dir }}/casp{{ season }}"
    if [[ ! -d "$src" ]]; then
        echo "!! $src not found -- run the casp_rank/casp_fetch/casp_batch/casp_dashboard pipeline for casp{{ season }} first" >&2
        exit 1
    fi
    rm -rf "$dst"
    mkdir -p "$dst"
    cp -R "$src"/. "$dst"/
    echo "✓ Synced $src -> $dst"

# Sync a season's CASP dashboard into the gallery (see `casp-sync`), then
# deploy the whole gallery -- same Cloudflare Pages project/domain as
# `viewers-deploy`, casp15/casp16 landing as subpaths beside the rest.
# Nothing outside <viewers_dir>/casp<season>/ is touched by the sync step,
# but the deploy step still re-uploads all of viewers_dir (that's how Pages
# deploys work), so anything else you meant to update in the main gallery or
# workstream1/ should be refreshed (`just viewers` / `just ws1-viewers`)
# before running this too. season is "15" or "16".
casp-deploy season: (casp-sync season)
    just viewers-deploy

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
