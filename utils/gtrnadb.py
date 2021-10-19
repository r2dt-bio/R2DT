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

import os
import re
import tempfile
import subprocess as sp
from pathlib import Path

from . import config
from . import shared


def setup():
    base = Path('/usr/local/lib/tRNAscan-SE/models')
    for path in base.glob('TRNAinf-*-iso'):
        domain = None
        if path.name == 'TRNAinf-arch-iso':
            domain = 'A'
        elif path.name == 'TRNAinf-bact-iso':
            domain = 'B'
        elif path.name == 'TRNAinf-euk-iso':
            domain = 'E'
        else:
            raise ValueError("Must be able to fetch from: " + path)

        with path.open('r') as raw:
            for line in raw:
                line = line.strip()
                if line.startswith('NAME'):
                    _, name = re.split('\s+', line, maxsplit=1)
                    _, isotype = name.split('-', 1)
                    get_trnascan_cm(domain, isotype)


def parse_trnascan_output(filename):
    """
    Sequence           		     tRNA	Bounds	tRNA	Anti	Intron  Bounds	Inf
    Name               	tRNA #	Begin	End	    Type	Codon	Begin	   End	Score	Note
    --------           	------	-----	------	----	-----	-----	----	------	------
    URS0000023412_9606 	1	1 	73	Thr	TGT	0	0	60.2
    """
    data = {}
    with open(filename, 'r') as f:
        for i, line in enumerate(f):
            if i in [0, 1, 2]:
                continue # skip 3 header lines
            parts = line.split('\t')
            data[parts[0].strip()] = {
                'score': float(parts[8].strip()),
                'isotype': parts[4].strip(),
                'note': parts[9].lower().strip(),
                'start': int(parts[2].strip()),
                'end': int(parts[3].strip()),
            }
    return data


def run_trnascan(fasta_input, output_folder, domain):
    output_file = os.path.join(output_folder, domain + '-' + os.path.basename(fasta_input).replace('.fasta', '.txt'))
    if not os.path.exists(output_file):
        cmd = 'tRNAscan-SE -{domain} -o {output_file} {input}'.format(
            domain=domain,
            input=fasta_input,
            output_file=output_file,
        )
        print(cmd)
        os.system(cmd)
    return parse_trnascan_output(output_file)


def skip_trna(entry):
    """
    Some tRNAs should not be drawn and need to be skipped.
    """
    if 'pseudo' in entry['note'] or entry['isotype'] in ['Undet', 'Sup']:
        return True
    return False


def classify_trna_sequences(fasta_input, output_folder):
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    bacteria = run_trnascan(fasta_input, output_folder, 'B')
    archaea = run_trnascan(fasta_input, output_folder, 'A')
    eukaryotes = run_trnascan(fasta_input, output_folder, 'E')

    rna_ids = set()
    rna_ids.update(list(bacteria.keys()), list(archaea.keys()), list(eukaryotes.keys()))

    data = []
    for rna_id in rna_ids:
        values = [
            bacteria[rna_id]['score'] if rna_id in bacteria else -1000,
            archaea[rna_id]['score'] if rna_id in archaea else -1000,
            eukaryotes[rna_id]['score'] if rna_id in eukaryotes else -1000,
        ]
        maximum = values.index(max(values))
        if maximum == 0:
            if skip_trna(bacteria[rna_id]):
                continue
            bacteria[rna_id]['domain'] = 'B'
            bacteria[rna_id]['id'] = rna_id
            data.append(bacteria[rna_id])
        elif maximum == 1:
            if skip_trna(archaea[rna_id]):
                continue
            archaea[rna_id]['domain'] = 'A'
            archaea[rna_id]['id'] = rna_id
            data.append(archaea[rna_id])
        else:
            if skip_trna(eukaryotes[rna_id]):
                continue
            eukaryotes[rna_id]['domain'] = 'E'
            eukaryotes[rna_id]['id'] = rna_id
            data.append(eukaryotes[rna_id])

    with open(os.path.join(output_folder, 'hits.txt'), 'w') as f_out:
        for entry in data:
            f_out.write('{}\t{}_{}\tPASS\n'.format(entry['id'], entry['domain'], entry['isotype']))
    return data


def visualise(domain, isotype, fasta_input, output_folder, test, constraint, exclusion):
    destination = '{}/{}'.format(output_folder, '_'.join([domain, isotype]))
    if not os.path.exists(destination):
        os.makedirs(destination)

    if not os.path.exists(fasta_input + '.ssi'):
        cmd = 'esl-sfetch --index {}'.format(fasta_input)
        os.system(cmd)

    cmd = "grep '>' {} > headers.txt"
    os.system(cmd.format(fasta_input))

    with open('headers.txt', 'r') as f:
        for i, line in enumerate(f):
            if test and i > 10:
                continue
            seq_id = line.split(' ', 1)[0].replace('>', '').strip()
            print(seq_id)
            generate_2d(domain, isotype, seq_id, None, None, fasta_input, destination, constraint, exclusion)
    os.system('rm headers.txt')


def get_trnascan_cm(domain, isotype):
    """
    Fetch a domain-specific isotype covariance model as a separate file.
    """
    if not os.path.exists(config.GTRNADB_CM_LIBRARY):
        os.mkdir(config.GTRNADB_CM_LIBRARY)
    cm_output = Path(config.GTRNADB_CM_LIBRARY) / '{}_{}.cm'.format(domain, isotype)
    if cm_output.exists():
        return str(cm_output)

    cm_library = Path('/usr/local/lib/tRNAscan-SE/models')
    if domain == 'A':
        cm_library = cm_library / 'TRNAinf-arch-iso'
        cm_name = 'arch-' + isotype
    elif domain == 'B':
        cm_library = cm_library / 'TRNAinf-bact-iso'
        cm_name = 'bact-' + isotype
    elif domain == 'E':
        cm_library = cm_library / 'TRNAinf-euk-iso'
        cm_name = 'euk-' + isotype
    else:
        raise ValueError("Unknown domain: %s" % domain)

    with cm_output.open('w') as out:
        cmd = ['cmfetch', str(cm_library), cm_name]
        sp.check_call(cmd, stdout=out)
    return cm_output


def get_traveler_template_xml(domain, isotype):
    if domain == 'A':
        return os.path.join(config.GTRNADB_ARCH, 'arch-{}-traveler-template.xml'.format(isotype))
    elif domain == 'B':
        return os.path.join(config.GTRNADB_BACT, 'bact-{}-traveler-template.xml'.format(isotype))
    else:
        return os.path.join(config.GTRNADB_EUK, 'euk-{}-traveler-template.xml'.format(isotype))


def get_traveler_fasta(domain, isotype):
    if domain == 'A':
        return os.path.join(config.GTRNADB_ARCH, 'arch-{}-traveler.fasta'.format(isotype))
    elif domain == 'B':
        return os.path.join(config.GTRNADB_BACT, 'bact-{}-traveler.fasta'.format(isotype))
    else:
        return os.path.join(config.GTRNADB_EUK, 'euk-{}-traveler.fasta'.format(isotype))


def generate_2d(domain, isotype, seq_id, start, end, fasta_input, output_folder, constraint, exclusion):
    temp_fasta = tempfile.NamedTemporaryFile()
    temp_sto = tempfile.NamedTemporaryFile()
    temp_depaired = tempfile.NamedTemporaryFile()
    temp_stk = tempfile.NamedTemporaryFile()
    temp_pfam_stk = tempfile.NamedTemporaryFile(delete=False)
    temp_afa = tempfile.NamedTemporaryFile()
    temp_map = tempfile.NamedTemporaryFile()

    if not os.path.exists(fasta_input + '.ssi'):
        cmd = 'esl-sfetch --index {}'.format(fasta_input)
        os.system(cmd)

    cmd = 'esl-sfetch {seq_range} {fasta_input} {seq_id} > {temp_fasta}'.format(
        fasta_input=fasta_input,
        seq_id=seq_id,
        temp_fasta=temp_fasta.name,
        seq_range='-c {}..{}'.format(start, end) if start and end else ''
    )
    os.system(cmd)

    cmd = "cmalign {trnascan_cm} {temp_fasta} > {temp_sto}".format(
        trnascan_cm=get_trnascan_cm(domain, isotype),
        temp_fasta=temp_fasta.name,
        temp_sto=temp_sto.name
    )
    result = os.system(cmd)
    if result:
        print("Failed cmalign of %s to %s" % (seq_id, isotype))
        return

    # remove non-canonical Watson-Crick basepairs (e.g. C:A in URS000008DB9C_7227)
    cmd = 'esl-alidepair.pl --nc 0.5 {} {}'.format(temp_sto.name, temp_depaired.name)
    result = os.system(cmd)
    if result:
        print("Failed esl-alidepair for {}".format(seq_id))

    cmd = 'esl-alimanip --rna --sindi --outformat pfam {} > {}'.format(temp_depaired.name, temp_stk.name)
    result = os.system(cmd)
    if result:
        print("Failed esl-alimanip for %s" % (seq_id))
        return

    cmd = 'ali-pfam-lowercase-rf-gap-columns.pl -s {} > {}'.format(temp_stk.name, temp_pfam_stk.name)
    result = os.system(cmd)
    if result:
        raise ValueError("Failed ali-pfam-lowercase-rf-gap-columns for %s" % (seq_id))

    os.system('cp {} {}/temp_pfam_stk.txt'.format(temp_pfam_stk.name, output_folder))
    shared.remove_large_insertions_pfam_stk(temp_pfam_stk.name)

    cmd = 'ali-pfam-sindi2dot-bracket.pl -l -n -w -a -c {} > {}'.format(temp_pfam_stk.name, temp_afa.name)
    result = os.system(cmd)
    if result:
        raise ValueError("Failed ali-pfam-sindi2dot-bracket for %s" % (seq_id))

    cmd = '/rna/python36/bin/python3.6 /rna/traveler/utils/infernal2mapping.py -i {} > {}'.format(temp_afa.name, temp_map.name)
    result = os.system(cmd)
    if result:
        raise ValueError("Failed infernal2mapping for %s" % (cmd))

    result_base = os.path.join(output_folder, seq_id.replace('/', '-') + '-' + domain + '_' + isotype)
    input_fasta = os.path.join(output_folder, seq_id + '.fasta')
    cmd = 'ali-pfam-sindi2dot-bracket.pl {} > {}'.format(temp_stk.name, input_fasta)
    os.system(cmd)

    if constraint:
        shared.fold_insertions(input_fasta, exclusion)
    elif exclusion:
        print('Exclusion ignored, enable --constraint to add exclusion file')

    log = result_base + '.log'
    cmd = ('traveler '
           '--verbose '
           '--target-structure {fasta} '
           '--template-structure --file-format traveler {traveler_template_xml} {traveler_fasta} '
           '--draw {map} {result_base} '
           '> {log}' ).format(
               fasta=input_fasta,
               result_base=result_base,
               traveler_template_xml=get_traveler_template_xml(domain, isotype),
               traveler_fasta=get_traveler_fasta(domain, isotype),
               log=log,
               map=temp_map.name
            )
    result = os.system(cmd)
    if result:
        print('Repeating using Traveler mapping')
        cmd = ('traveler '
               '--verbose '
               '--target-structure {fasta} '
               '--template-structure --file-format traveler {traveler_template_xml} {traveler_fasta} '
               '--all {result_base} '
               '> {log}' ).format(
                   fasta=input_fasta,
                   result_base=result_base,
                   traveler_template_xml=get_traveler_template_xml(domain, isotype),
                   traveler_fasta=get_traveler_fasta(domain, isotype),
                   log=log
                )
        print(cmd)
        result = os.system(cmd)

    print(cmd)

    temp_fasta.close()
    temp_sto.close()
    temp_depaired.close()
    temp_stk.close()
    temp_pfam_stk.close()
    temp_afa.close()
    temp_map.close()

    cmd = 'rm -f {0}/*.xml {0}/*.ps'.format(output_folder)
    os.system(cmd)

    overlaps = 0
    with open(log, 'r') as raw:
        for line in raw:
            match = re.search(r'Overlaps count: (\d+)', line)
            if match:
                if overlaps:
                    print('ERROR: Saw too many overlap counts')
                    break
                overlaps = int(match.group(1))

    with open(result_base + '.overlaps', 'w') as out:
        out.write(str(overlaps))
        out.write('\n')
