
import os
from sys import argv


CM_LIBRARY = '/rna/data/cms'

CRW_PS_LIBRARY = '/rna/auto-traveler/data/crw-ps'
CRW_FASTA_LIBRARY = '/rna/auto-traveler/data/crw-fasta'


fasta_input = argv[1]  # /rna/examples/examples.fasta
output_folder = argv[2]  # output-test


def get_ribotyper_output():
    ribotyper_long_out = os.path.join(output_folder, os.path.basename(output_folder) + '.ribotyper.long.out')
    if not os.path.exists(ribotyper_long_out):
        cmd = 'perl /rna/ribotyper-v1/ribotyper.pl -i {CM_LIBRARY}/all.modelinfo.txt -f {fasta_input} {output_folder}'.format(
            CM_LIBRARY=CM_LIBRARY,
            fasta_input=fasta_input,
            output_folder=output_folder
        )
        print cmd
        os.system(cmd)
    f_out = os.path.join(output_folder, 'hits.txt')
    cmd = "cat %s | grep -v '^#' | grep PASS | awk -v OFS='\t' '{print $2, $8, $3}' > %s" % (ribotyper_long_out, f_out)
    os.system(cmd)
    return f_out


def main():

    with open(get_ribotyper_output(), 'r') as f:
        for line in f.readlines():
            rnacentral_id, model_id, _ = line.split('\t')
            print model_id

            cmd = 'esl-sfetch %s %s > temp.fasta' % (fasta_input, rnacentral_id)
            os.system(cmd)

            cmd = "cmalign %s.cm temp.fasta > temp.sto" % os.path.join(CM_LIBRARY, model_id)
            os.system(cmd)

            cmd = 'esl-alimanip --sindi --outformat pfam temp.sto > temp.stk'
            os.system(cmd)

            cmd = 'perl /rna/jiffy-infernal-hmmer-scripts/ali-pfam-sindi2dot-bracket.pl temp.stk > %s/%s-%s.fasta' % (output_folder, rnacentral_id, model_id)
            os.system(cmd)

            cmd = ('traveler '
            '--target-structure {output}/{rnacentral_id}-{model_id}.fasta '
            '--template-structure {CRW_PS_LIBRARY}/{model_id}.ps {CRW_FASTA_LIBRARY}/{model_id}.fasta '
            '--all {output}/{rnacentral_id}-{model_id}').format(
                rnacentral_id=rnacentral_id,
                model_id=model_id,
                output=output_folder,
                CRW_PS_LIBRARY=CRW_PS_LIBRARY,
                CRW_FASTA_LIBRARY=CRW_FASTA_LIBRARY
            )
            print cmd
            os.system(cmd)


if __name__ == '__main__':
    main()
