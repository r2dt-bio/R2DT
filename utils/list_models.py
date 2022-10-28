"""
Copyright [2009-present] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
     http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""


import glob
import json
import os
import re

from collections import defaultdict
from utils import config


def get_model_type(model_id):
    folder_mapping = {
        'GtRNAdb': 'gtrnadb',
        'Rfam': 'rfam',
        'CRW': 'crw',
        'RiboVision SSU': 'ribovision_ssu',
        'RiboVision LSU': 'ribovision_lsu',
        'RNAse P Database': 'rnasep',
    }
    model_type = None
    with open(os.path.join(config.DATA, 'models.json'), 'r') as models_json:
        data = json.load(models_json)
        for model in data:
            if model_id == model['model_id']:
                if model['source'] in folder_mapping:
                    model_type = folder_mapping[model['source']]
                break
    return model_type


def get_gtrnadb_models():
    data = []
    for cm_file in glob.glob(os.path.join(config.GTRNADB_CM_LIBRARY, '*.cm')):
        model_id = os.path.basename(cm_file.replace('.cm', ''))
        parts = model_id.split('_')
        short_domain = parts[0]
        if short_domain == 'E':
            domain = 'Eukaryotes'
        elif short_domain == 'B':
            domain = 'Bacteria'
        elif short_domain == 'A':
            domain = 'Archaea'
        elif short_domain == 'M':
            domain = 'Vertebrates mito'
        else:
            domain = ''
        data.append({
            'model_id': model_id,
            'source': 'GtRNAdb',
            'description': 'tRNA {} ({})'.format(' '.join(parts[1:]), domain),
        })
    return data


def parse_metadata(metadata_file):
    metadata = {}
    rna_types = {}
    with open(metadata_file, 'r') as f_metadata:
        for i, line in enumerate(f_metadata):
            if i == 0:
                continue
            fields = line.split('\t')
            metadata[fields[0]] = fields[1].strip()
            rna_types[fields[0]] = fields[3].strip() if len(fields) > 3 else ''
    return metadata, rna_types


def parse_modelinfo(modelinfo_file):
    model_ids = []
    with open(modelinfo_file, 'r') as f_modelinfo:
        for line in f_modelinfo:
            if 'all.cm' in line:
                continue
            fields = re.split(r'\s+', line.strip())
            model_ids.append(fields[0])
    return model_ids


def get_crw_models():
    data = []
    modelinfo_file = os.path.join(config.CRW_CM_LIBRARY, 'modelinfo.txt')
    metadata_file = os.path.join(config.DATA, 'crw-metadata.tsv')

    metadata, rna_types = parse_metadata(metadata_file)
    model_ids = parse_modelinfo(modelinfo_file)

    for model_id in model_ids:
        if model_id in metadata:
            species = metadata[model_id]
        else:
            species = ''
        rna_type = rna_types[model_id].replace('_', ' ')
        if 'intron' in rna_type:
            continue
        data.append({
            'model_id': model_id,
            'source': 'CRW',
            'description': '{} {} {}'.format(species, rna_type, model_id),
        })

    return data


def get_qualifier(model_id):
    """
    Return a manually curated string that helps distinguish similarly named
    templates.

    For example, HS_LSU_3D and mHS_LSU_3D are both Human LSUs but one is found
    in mitochondria.
    """
    if model_id in ['mt_TetT_LSU_3D', 'mHS_LSU_3D']:
        qualifier = ' (mitochondrial)'
    elif model_id in ['cSO_23S_3D']:
        qualifier = ' (chloroplast)'
    else:
        qualifier = ''
    return qualifier


def get_models(source, modelinfo_file, metadata_file):
    data = []
    metadata, _ = parse_metadata(metadata_file)
    model_ids = parse_modelinfo(modelinfo_file)

    if source == 'RiboVision SSU':
        rna_type = 'small subunit rRNA'
    elif source == 'RiboVision LSU':
        rna_type = 'large subunit rRNA'
    elif source == 'RNAse P Database':
        rna_type = 'RNAse P'
    else:
        rna_type = ''

    for model_id in model_ids:
        species = metadata.get(model_id, '')
        data.append({
            'model_id': model_id,
            'source': source,
            'description': '{} {}{}'.format(species, rna_type, get_qualifier(model_id)),
        })
    return data


def get_rfam_models():
    data = [{
        'model_id': 'RF00005',
        'source': 'Rfam',
        'description': 'tRNA RF00005 (Rfam)',
    }]
    modelinfo_file = os.path.join(config.DATA, 'cms', 'rfam', 'modelinfo.txt')
    model_ids = parse_modelinfo(modelinfo_file)
    accessions = {}

    descriptions = {}
    with open(os.path.join(config.RFAM_DATA, 'family.txt'), encoding='utf8', errors='ignore') as family_file:
        for line in family_file:
            fields = re.split(r'\t', line)
            descriptions[fields[1]] = fields[3]
            accessions[fields[1]] = fields[0]

    for model_id in model_ids:
        data.append({
            'model_id': model_id,
            'source': 'Rfam',
            'description': '{} ({})'.format(descriptions[model_id], accessions[model_id]),
        })
    return data


def check_unique_descriptions(data):
    """
    Check if any templates have identical descriptions. This is undesirable
    because the descriptions are used to manually select a template in the web
    interface.
    """
    counts = defaultdict(list)
    for entry in data:
        counts[entry['description']].append(entry['model_id'])
    for description in counts.keys():
        if len(counts[description]) > 1:
            msg = 'Please fix identical descriptions "{}": {}'.format(
                    description, '; '.join(counts[description]))
            print(msg)


def list_models():
    data = []
    data = data + get_gtrnadb_models()
    data = data + get_crw_models()
    data = data + get_models('RiboVision LSU', os.path.join(config.RIBOVISION_LSU_CM_LIBRARY, 'modelinfo.txt'), os.path.join(config.RIBOVISION_LSU, 'metadata.tsv'))
    data = data + get_models('RiboVision SSU', os.path.join(config.RIBOVISION_SSU_CM_LIBRARY, 'modelinfo.txt'), os.path.join(config.RIBOVISION_SSU, 'metadata.tsv'))
    data = data + get_models('RNAse P Database', os.path.join(config.RNASEP_CM_LIBRARY, 'modelinfo.txt'), os.path.join(config.RNASEP, 'metadata.tsv'))
    data = data + get_rfam_models()
    data.sort(key=lambda x: x['description'])
    return data
