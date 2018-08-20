


import collections
import glob
import os
import json



pngs = 'png'


data = collections.defaultdict(list)
with open('pdb-rnacentral.json', 'r') as f:
    raw = json.load(f)
    for entry in raw:
        data[entry['upi'] + '_' + str(entry['taxid'])].append(entry)

# import pdb; pdb.set_trace()


head = """
<html>
<head>
<link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.1.3/css/bootstrap.min.css" integrity="sha384-MCw98/SFnGE8fJT3GXwEOngsV7Zt27NXFoaoApmYm81iuXoPkFOJwJ8ERdknLPMO" crossorigin="anonymous">
<script src="https://stackpath.bootstrapcdn.com/bootstrap/4.1.3/js/bootstrap.min.js" integrity="sha384-ChfqqxuZUCnJSK3+MXmPNIyE6ZbWh2IMqE241rYiqJxyMiZ6OW/JmZQ5stwEULTy" crossorigin="anonymous"></script>
</head>
<body>
<div class="container">
<div class="row">
<h1>RNA secondary structures for sequences from PDB generated using <a href="https://github.com/rnacentral/auto-traveler">auto-traveler</a></h1>

<p>
Currently only small ribosomal subunit and 5S rRNA have been analysed.
</p>
<hr>
</div>
"""

ordered_png = []
for png in glob.glob(os.path.join(pngs, '*.png')):
    rna_id, model_id = os.path.basename(png).replace('.colored.svg.png', '').split('-')
    try:
        length = data[rna_id][0]['len']
    except:
        continue
    ordered_png.append((png, length))
ordered_png = sorted(ordered_png, key=lambda k: k[1])


with open('index.html', 'w') as html_output:
    html_output.write(head)
    for entry in reversed(ordered_png): # sorted(glob.glob(os.path.join(pngs, '*.png'))):
        png = entry[0]
        # URS0000CBFF1F_5693-d.16.e.L.major.colored.svg.png
        rna_id, model_id = os.path.basename(png).replace('.colored.svg.png', '').split('-')
        try:
            first_entry = data[rna_id][0]
        except:
            continue
        html_output.write("<div class='row'><h2><a href='http://rnacentral.org/rna/{rna_id}'>{rna_id}</a></h2>".format(rna_id=rna_id, model_id=model_id))

        html = """
            <ul class="text-muted">
                <li>CRW template: <strong>{model_id}</strong></li>
                <li>3D structure species: <a href="https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id={taxid}">{classification}</a></li>
            </ul>
        """.format(classification=first_entry['classification'], taxid=first_entry['taxid'], model_id=model_id)

        html_output.write(html)
        img = '<div class="col-md-6"><a href="./svg-colored/{rna_id}-{model_id}.colored.svg"><img class="img-thumbnail rounded" title="Click to view in SVG" src=png/{rna_id}-{model_id}.colored.svg.png style="max-height: 500px;"></a></div>'.format(rna_id=rna_id, model_id=model_id)
        html_output.write(img)
        html_output.write('<div class="col-md-6" style="max-height: 500px; overflow: auto;"><table class="table table-hover">')
        html_output.write("<thead><th>PDB</th><th>Description</th></thead>")
        for pdb_entry in sorted(data[rna_id], key=lambda k: k['ac']):
            pdb_id, pdb_chain, pdb_entry_id = pdb_entry['ac'].split('_')
            pdb_url = '<a href="https://www.ebi.ac.uk/pdbe/entry/pdb/{}/RNA/{}">{}</a>'.format(pdb_id, pdb_entry_id, pdb_id)
            row = [pdb_url, pdb_entry['description']]
            line = "<tr><td>" + '</td><td>'.join(row) + '</td></tr>'
            html_output.write(line)
            # print pdb_entry
        html_output.write('</table></div></div><hr>')


    tail = """
    </div>
    </body>
    </html>
    """
    html_output.write(tail)
print 'Done'
