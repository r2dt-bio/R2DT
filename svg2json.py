#!/usr/bin/env python3

import json
from xml.dom import minidom

import click


@click.command()
@click.argument('pdb_id')
@click.argument('filename')
@click.argument('output', type=click.File('w'))
def main(pdb_id, filename, output):
    pdb, model, chain = pdb_id.split('_')
    doc = minidom.parse(filename)
    svg_lines = []

    residues = []
    prev = None
    index = 0
    width = float(doc.documentElement.getAttribute('width'))
    height = float(doc.documentElement.getAttribute('height'))
    for nt in doc.getElementsByTagName('text'):
        if 'numbering-label' in nt.getAttribute('class'):
            continue

        x2 = float(nt.getAttribute('x'))
        y2 = float(nt.getAttribute('y'))
        if prev is None:
            prev = (x2 - 5, y2 - 5)
        x1 = prev[0]
        y1 = prev[1]
        residues.append({
            'resnum': index + 1,
            'path': (x1, y1, x2, y2),
        })
        prev = (x2, y2)
        index += 1


    data = {
        pdb: {
            model: {
                chain: {
                    'rna_nucleotides': residues
                }
            }
        },
        'dimensions': {'width': width, 'height': height},
    }

    if not residues:
        raise ValueError("Did not extract any coordinates")

    json.dump(data, output)


if __name__ == '__main__':
    main()
