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
import os
import re

from utils import config


def get_gtrnadb_models():
    data = []
    for cm_file in glob.glob(os.path.join(config.GTRNADB_CM_LIBRARY, '*.cm')):
        model_id = os.path.basename(cm_file.replace('.cm', ''))
        short_domain, isotype = model_id.split('_')
        if short_domain == 'E':
            domain = 'Eukaryotes'
        elif short_domain == 'B':
            domain = 'Bacteria'
        elif short_domain == 'A':
            domain = 'Archaea'
        else:
            domain = ''
        data.append({
            'model_id': model_id,
            'source': 'GtRNAdb',
            'description': 'tRNA {} ({})'.format(isotype, domain),
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
        if model_id in metadata:
            species = metadata[model_id]
        else:
            species = ''
        data.append({
            'model_id': model_id,
            'source': source,
            'description': '{} {}'.format(species, rna_type),
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
            'description': '{} {}'.format(accessions[model_id], descriptions[model_id]),
        })
    return data


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
