"""
This module contains a script to scale the coordinates in Traveler templates.
"""

import argparse
import glob
import os
import xml.etree.ElementTree as ET

SCALING_FACTOR = 4


def scale_coordinates(file_path, scaling_factor=SCALING_FACTOR):
    """Scale the coordinates in a Traveler template."""
    tree = ET.parse(file_path)
    root = tree.getroot()

    for point in root.findall("point"):
        x_coord = float(point.get("x"))
        y_coord = float(point.get("y"))
        point.set("x", str(x_coord * scaling_factor))
        point.set("y", str(y_coord * scaling_factor))
    tree.write(file_path)


def main():
    """Parse command-line arguments and scale the coordinates in Traveler templates."""
    parser = argparse.ArgumentParser(
        description="Scale the coordinates in Traveler templates."
    )
    parser.add_argument(
        "path",
        help="The path to the Traveler template or directory containing templates.",
    )
    parser.add_argument(
        "--scaling-factor",
        type=float,
        default=SCALING_FACTOR,
        help="The scaling factor to apply to the coordinates.",
    )
    args = parser.parse_args()

    path = args.path
    scaling_factor = args.scaling_factor

    if os.path.isdir(path):
        # Recursively find all files called traveler-template.xml
        for file_path in glob.glob(
            os.path.join(path, "**", "*traveler-template.xml"), recursive=True
        ):
            print(
                f"Scaling coordinates in {file_path} with scaling factor {scaling_factor}"
            )
            scale_coordinates(file_path, scaling_factor)
    elif os.path.isfile(path):
        print(f"Scaling coordinates in {path} with scaling factor {scaling_factor}")
        scale_coordinates(path, scaling_factor)
    else:
        print(f"Error: The path '{path}' is neither a file nor a directory.")


if __name__ == "__main__":
    main()
