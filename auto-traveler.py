#!/usr/bin/env python

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

import re
import os

import click


CM_LIBRARY = '/rna/auto-traveler/data/cms'
CRW_PS_LIBRARY = '/rna/auto-traveler/data/crw-ps'
CRW_FASTA_LIBRARY = '/rna/auto-traveler/data/rnacentral-nopbpseq-2019-02-15'


def get_ribotyper_output(fasta_input, output_folder, cm_library):
    ribotyper_long_out = os.path.join(output_folder, os.path.basename(output_folder) + '.ribotyper.long.out')
    if not os.path.exists(ribotyper_long_out):
        cmd = 'ribotyper.pl -i {CM_LIBRARY}/modelinfo.txt -f {fasta_input} {output_folder}'.format(
            CM_LIBRARY=cm_library,
            fasta_input=fasta_input,
            output_folder=output_folder
        )
        print(cmd)
        os.system(cmd)
    f_out = os.path.join(output_folder, 'hits.txt')
    cmd = "cat %s | grep -v '^#' | grep -v MultipleHits | grep PASS | awk -v OFS='\t' '{print $2, $8, $3}' > %s" % (ribotyper_long_out, f_out)
    os.system(cmd)
    return f_out


@click.command()
@click.option('--cm-library', default=CM_LIBRARY)
@click.option('--ps-library', default=CRW_PS_LIBRARY)
@click.option('--fasta-library', default=CRW_FASTA_LIBRARY)
@click.argument('fasta-input', type=click.Path())
@click.argument('output_folder', type=click.Path())
def main(fasta_input, output_folder, cm_library=None, ps_library=None, fasta_library=None):

    os.system('mkdir -p %s' % output_folder)

    with open(get_ribotyper_output(fasta_input, output_folder, cm_library), 'r') as f:
        for line in f.readlines():
            rnacentral_id, model_id, _ = line.split('\t')
            found = '%s-%s' % (rnacentral_id, model_id)

            cmd = 'esl-sfetch %s %s > temp.fasta' % (fasta_input, rnacentral_id)
            os.system(cmd)

            cmd = "cmalign %s.cm temp.fasta > temp.sto" % os.path.join(cm_library, model_id)
            os.system(cmd)

            cmd = 'esl-alimanip --sindi --outformat pfam temp.sto > temp.stk'
            os.system(cmd)

            cmd = 'ali-pfam-sindi2dot-bracket.pl temp.stk > %s/%s-%s.fasta' % (output_folder, rnacentral_id, model_id)
            os.system(cmd)

            result_base = os.path.join(output_folder, '{rnacentral_id}-{model_id}'.format(
                rnacentral_id=rnacentral_id,
                model_id=model_id,
            ))

            log = result_base + '.log'
            cmd = ('traveler '
                   '--verbose '
                   '--target-structure {result_base}.fasta '
                   '--template-structure {ps_library}/{model_id}.ps {fasta_library}/{model_id}.fasta '
                   '--all {result_base} > {log}').format(
                       result_base=result_base,
                       model_id=model_id,
                       ps_library=ps_library,
                       fasta_library=fasta_library,
                       log=log,
                   )
            print(cmd)
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



if __name__ == '__main__':
    main()
