"""
Use this script to build R2DT Docker images for amd64 and arm64 architectures
using the multi-Dockerfile approach.

Usage:
    # build both amd64 and arm64 images using the default Dockerfiles
    python3 docker-build.py

    # specify architecture and Dockerfiles
    python3 docker-build.py [arch] [base_dockerfile] [main_dockerfile]
"""

import argparse
import hashlib
import logging
import os
from pathlib import Path
from typing import List


def get_dockerfile_hash(dockerfile: Path) -> str:
    """Get md5 hash of the Dockerfile."""
    hasher = hashlib.md5()
    with dockerfile.open("rb") as f_dockerfile:
        buf = f_dockerfile.read()
        hasher.update(buf)
    return hasher.hexdigest()[:8]


def get_tag(dockerfile: Path, arch: str) -> str:
    """Get tag for the Docker image."""
    return f"{get_dockerfile_hash(dockerfile)}-{arch.replace('linux/', '')}"


def build_base_image(dockerfile: Path, arch: str) -> str:
    """Build the base image."""
    tag = get_tag(dockerfile, arch)
    cmd = (
        f"docker build -t rnacentral/r2dt-base:{tag} --platform={arch} "
        f"-f {str(dockerfile)} base_image"
    )
    logging.info("Building base image: %s", cmd)
    result = os.system(cmd)
    if result != 0:
        raise RuntimeError("Failed to build base image")
    return tag


def build_main_image(dockerfile: Path, base_tag: str, arch: str) -> str:
    """Build the main image."""
    tag = get_tag(dockerfile, arch)
    # create temporary Dockerfile referencing the base image tag
    with dockerfile.open("r") as f_dockerfile:
        lines = f_dockerfile.readlines()
    temp_dockerfile = Path("Dockerfile.tmp")
    with open(temp_dockerfile, "w") as f_dockerfile:
        for line in lines:
            if line.startswith("FROM"):
                f_dockerfile.write(f"FROM rnacentral/r2dt-base:{base_tag}\n")
            else:
                f_dockerfile.write(line)
    cmd = f"docker build -t rnacentral/r2dt:{tag} --platform={arch} -f {temp_dockerfile} ."
    logging.info("Building main image: %s", cmd)
    result = os.system(cmd)
    if result != 0:
        print(cmd)
        raise RuntimeError("Failed to build main image")
    temp_dockerfile.unlink()
    return tag


def build_images(base_dockerfile: Path, main_dockerfile: Path, arch: str) -> List[str]:
    """Build R2DT Docker images."""
    if arch not in ["linux/amd64", "linux/arm64", "amd64", "arm64", "all"]:
        raise ValueError(f"Unsupported architecture: {arch}")
    if arch == "all":
        archs = ["linux/amd64", "linux/arm64"]
    elif arch == "amd64":
        archs = ["linux/amd64"]
    elif arch == "arm64":
        archs = ["linux/arm64"]
    else:
        archs = [arch]

    for _arch in archs:
        base_tag = build_base_image(base_dockerfile, _arch)
        main_tag = build_main_image(main_dockerfile, base_tag, _arch)
        logging.info("Successfully built rnacentral/r2dt-base:%s", base_tag)
        logging.info("Successfully built rnacentral/r2dt:%s", main_tag)


def main():
    """Parse command-line arguments and build R2DT Docker images."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "arch",
        help="image architecture",
        choices=["linux/amd64", "linux/arm64", "amd64", "arm64", "all"],
        nargs="?",
        default="all",
    )
    parser.add_argument(
        "base_dockerfile",
        help="base Dockerfile",
        type=Path,
        nargs="?",
        default=Path("base_image/Dockerfile"),
    )
    parser.add_argument(
        "main_dockerfile",
        help="main Dockerfile",
        type=Path,
        nargs="?",
        default=Path("Dockerfile"),
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    build_images(args.base_dockerfile, args.main_dockerfile, args.arch)


if __name__ == "__main__":
    main()
