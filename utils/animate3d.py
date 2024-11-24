import subprocess
import sys
import os
from pathlib import Path
import argparse

def process_bulk(ref_pdb, query_dir):
    ref_name = ref_pdb.stem
    for query_pdb in query_dir.iterdir():
        if query_pdb.suffix == '.pdb' and query_pdb.name != ref_pdb.name and not query_pdb.name.endswith('tmp.pdb'):
            query_name = query_pdb.stem
            ref_svg = ref_pdb.with_suffix('.colored.svg')
            query_svg = query_pdb.with_suffix('.colored.svg')
            output_svg = ref_pdb.parent / f"{ref_name}_to_{query_name}.animated.svg"
            print(f"Processing {ref_pdb}")
            run_rnaview(ref_pdb)
            print(f"Processing {query_pdb}")
            run_rnaview(query_pdb)
            print(f"Animating {output_svg}")
            run_animate(ref_svg, query_svg, output_svg)

def run_rnaview(pdb_file):
    command = f"python ./utils/rnaview.py {pdb_file}"
    print(f"Running command: {command}")
    subprocess.run(command, shell=True, check=True)

def run_animate(ref_svg, query_svg, output_svg):
    command = f"python ./utils/animate.py {ref_svg} {query_svg} {output_svg}"
    print(f"Running command: {command}")
    subprocess.run(command, shell=True, check=True)

def main(ref_pdb, query_pdb_or_dir, bulk=False):
    ref_pdb = Path(ref_pdb)
    query_pdb_or_dir = Path(query_pdb_or_dir)
    
    if bulk:
        process_bulk(ref_pdb, query_pdb_or_dir)
    else:
        ref_svg = ref_pdb.with_suffix('.colored.svg')
        query_svg = query_pdb_or_dir.with_suffix('.colored.svg')
        output_svg = ref_pdb.parent / f"{ref_pdb.stem}_to_{query_pdb_or_dir.stem}.animated.svg"

        run_rnaview(ref_pdb)
        run_rnaview(query_pdb_or_dir)
        run_animate(ref_svg, query_svg, output_svg)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Animate 3D structures.")
    parser.add_argument("ref_pdb", type=str, help="Reference PDB file")
    parser.add_argument("query_pdb_or_dir", type=str, help="Query PDB file or directory")
    parser.add_argument("-b", "--bulk", action="store_true", help="Process in bulk mode")
    
    args = parser.parse_args()
    main(args.ref_pdb, args.query_pdb_or_dir, args.bulk)