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


import csv
import os
import glob


BPSEQ_LOCATION = '/rna/auto-traveler/data/crw-bpseq'
CM_LIBRARY = '/rna/auto-traveler/data/cms'
CRW_FASTA_NO_PSEUDOKNOTS = '/rna/auto-traveler/data/crw-fasta-no-pseudoknots'


def get_crw_metadata(filename):
    """
    Example:
    structure	rna_type	rna_class	tax_id	accession(s)	cell_location
    a.I1.e.C.luteoviridis.B.C1.SSU.1052.ps	I	IC1	31301	X73998	Nucleus
    d.16.e.C.luteoviridis.ps	R	16S	31301	X73998	Nucleus
    """
    data = {}
    with open(filename, 'r') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            model_id = row['structure'].replace('.ps', '')
            data[model_id] = {
                'rna_type': row['rna_type'],
                'rna_class': row['rna_class'],
                'taxid': row['tax_id'],
                'accessions': row['accession(s)'].split(' '),
                'cell_location': row['cell_location'],
            }
    return data


def convert_bpseq_to_fasta(bpseq):
    fasta = bpseq.replace('.bpseq', '-with-knots.fasta')
    if not os.path.exists(fasta):
        cmd = 'python /rna/traveler/utils/bpseq2fasta.py -i {bpseq} -o {fasta}'.format(
            bpseq=bpseq,
            fasta=fasta
        )
        os.system(cmd)
    return fasta


def break_pseudoknots(fasta):
    fasta_no_knots = fasta.replace('-with-knots.fasta', '.fasta')
    if not os.path.exists(fasta_no_knots):
        cmd = 'RemovePseudoknots -b {fasta} {fasta_no_knots}'.format(
            fasta=fasta,
            fasta_no_knots=fasta_no_knots
        )
        os.system(cmd)
    return fasta_no_knots


def convert_fasta_to_stockholm(fasta):
    stockholm = fasta.replace('.fasta', '.sto')
    model_id = os.path.basename(stockholm).replace('.sto', '')
    if not os.path.exists(stockholm):
        with open(fasta, 'r') as f_input:
            with open(stockholm, 'w') as f_output:
                   lines = f_input.readlines()
                   f_output.write('# STOCKHOLM 1.0\n')
                   f_output.write('\n')
                   f_output.write('{0}{1}\n'.format(model_id.ljust(60), lines[1].strip()))
                   f_output.write('{0}{1}\n'.format('#=GC SS_cons'.ljust(60), lines[2].strip()))
                   f_output.write('//\n')
    return stockholm


def copy_cm_evalues(cm):
    """
    Update covariance files genenrated from CRW covariance models
    by copying E-values from Rfam CMs.
    """
    if not os.path.exists('RF00177.cm'):
        cmd = 'wget -O RF00177.cm http://rfam.org/family/RF00177/cm'
        os.system(cmd)
    cmd = 'perl /rna/jiffy-infernal-hmmer-scripts/cm-copy-evalue-parameters.pl RF00177.cm {cm}'.format(cm=cm)
    os.system(cmd)


def build_cm(stockholm):
    cm = os.path.join(CM_LIBRARY, os.path.basename(stockholm).replace('.sto', '.cm'))
    if not os.path.exists(cm):
        cmd = 'cmbuild {cm} {stockholm}'.format(
            cm=cm,
            stockholm=stockholm
        )
        os.system(cmd)
        copy_cm_evalues(cm)
    return cm


def main():

    for bpseq in glob.glob('%s/*.bpseq' % BPSEQ_LOCATION)[:2]:
        print os.path.basename(bpseq).replace('.bpseq', '')
        fasta = convert_bpseq_to_fasta(bpseq)
        fasta_no_knots = break_pseudoknots(fasta)
        stockholm = convert_fasta_to_stockholm(fasta_no_knots)
        build_cm(stockholm)
    crw_metadata = get_crw_metadata('data/rnacentral-nopbpseq-2019-02-15.tsv')

    print 'Done'


if __name__ == '__main__':
    main()
