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

import click


def generate_traveler_fasta(rfam_acc):
    """
    Generate fasta format for Rfam consensus.

    Example:

    >RF00162
    NNCUUAUCNAGAGNGGYRGAGGGAYNGGCCCNRUGAARCCNCRGCAACCNNYNNNNNNNNNRNNANGGUGCYAANUCCNRCNRNNNNNNNNNNNYNGRRAGAURAGRR
    ((((((((......(((...(((.....)))......)))..(((.((((((.........)))..))))))........((((.........))))...))))))))
    """

    ss_cons = ''
    consensus = ''

    # get a list of alignments
    seeds = []
    for seed in glob.glob(os.path.join('temp', rfam_acc, '*.R2R.sto')):
        seeds.append(seed)
    if len(seeds) != 1:
        print("Error: unusual number of seed alignments")

    with open(seeds[0], 'r') as f:
        for line in f.readlines():
            if line.startswith('#=GC SS_cons '):
                parts = line.split()
                ss_cons += parts[2].replace('<', '(').replace('>', ')')
            elif line.startswith('#=GC cons') and 'conss' not in line:
                parts = line.split()
                consensus += parts[2]

    if not ss_cons:
        print('No SS_CONS found')
    elif len(ss_cons) == len(consensus):
        if '-' in consensus:
            new_ss_cons = []
            new_consensus = []
            for i, nt in enumerate(consensus):
                if nt == '-' and ss_cons[i] == '.':
                    pass
                elif nt == '-' and ss_cons[i] in '<>':
                    # example RF00016
                    new_ss_cons.append(ss_cons[i])
                    new_consensus.append('N')
                else:
                    new_ss_cons.append(ss_cons[i])
                    new_consensus.append(consensus[i])
            ss_cons = ''.join(new_ss_cons)
            consensus = ''.join(new_consensus)

        with open(os.path.join('temp', rfam_acc, '{}-traveler.fasta'.format(rfam_acc)), 'w') as f:
            f.write('>{}\n'.format(rfam_acc))
            f.write('{}\n'.format(consensus.upper()))
            f.write('{}\n'.format(ss_cons))
    else:
        print('Error: structure and consensus have different lengths')


def convert_path_to_text(line):
    """
    <!--
    <path
     fill="#d90000" stroke="#000000" stroke-width="0.72" stroke-linecap="butt" stroke-linejoin="miter"
     d="M 91.0717 207.452 A 2.5575,2.5575 0 0,1 85.9567,207.452A 2.5575,2.5575 0 0,1 91.0717,207.452Z"/>
     -->

    <text x="85.9567" y="210.352" id="text1002">
      <tspan fill="#807b88"  font-variant="normal" font-weight="normal" font-style="normal" font-family="Bitstream Vera Sans" font-size="7.5" id="tspan1003">N</tspan>
    </text>
    """
    try:
        match = re.search(r'd="M (\d+(\.\d+)?) (\d+(\.\d+)?) ', line)
        if match:
            x = float(match.group(1))
            y = float(match.group(3))
            # new_x = x - 3.75 - (0.72 * 2)
            new_x = x - 5.115
            new_y = y + 2.9

            text = """
            <text x="{}" y="{}" id="foobar">
              <tspan fill="#807b88"  font-variant="normal" font-weight="normal" font-style="normal" font-family="Bitstream Vera Sans" font-size="7.5">N</tspan>
            </text>"""

            xml = """<point x="{:.2f}" y="{:.2f}" b="N"/>\n"""

            return (text.format(new_x, new_y), xml.format(new_x, new_y))
        else:
            print line
            print 'no match text'
    except:
        # import pdb; pdb.set_trace()
        print line

def convert_text_to_xml(line):
    """
    # <text x="209.519" y="231.111" id="text1002"><tspan x="209.519" y="231.111" fill="#807b88"  font-variant="normal" font-weight="normal" font-style="normal" font-family="Bitstream Vera Sans" font-size="7.5" id="tspan1003">R</tspan></text>
    # <point x="209.52" y="0.52" b="231.111"/>

    <text x="85.8067" y="195.025" id="text1002"><tspan x="85.8067" y="195.025" fill="#807b88"  font-variant="normal" font-weight="normal" font-style="normal" font-family="Bitstream Vera Sans" font-size="7.5" id="tspan1003">C</tspan></text>
    """
    # match = re.search(r'<text x="(\d+\.\d+)" y="(\d+\.\d+)".+>(\w)</tspan></text>', line)
    match = re.search(r'<text x="(\d+(\.\d+)?)" y="(\d+(\.\d+)?)".+>(\w)</tspan></text>', line)

    if match:
        point = '<point x="{:.2f}" y="{:.2f}" b="{}"/>\n'
        return point.format(float(match.group(1)), float(match.group(3)), match.group(5))
    else:
        print line
        print 'no match'


def download_rfam_seed(rfam_acc):
    output = os.path.join('temp', rfam_acc, '{}.seed'.format(rfam_acc))
    if not os.path.exists(output):
        url = 'http://rfam.org/family/{}/alignment'.format(rfam_acc)
        cmd = 'wget -O {output} {url}'.format(output=output, url=url)
        os.system(cmd)
    return output


def get_all_rfam_acc():
    rfam_accs = []
    blacklist = ['RF02541', 'RF00177', 'RF01960', 'RF02540', 'RF02543', 'RF01959', 'RF02542', 'RF02546', 'RF02545', 'RF00009']
    family_file = os.path.join('temp', 'family.txt')
    if not os.path.exists(family_file):
        cmd = 'wget -O {0}.gz ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/family.txt.gz && gunzip {0}.gz'.format(family_file)
        os.system(cmd)
    with open(family_file, 'r') as f:
        for line in f.readlines():
            if line.startswith('RF'):
                rfam_acc = line[:7]
                if rfam_acc in blacklist:
                    continue
                rfam_accs.append(rfam_acc)
    print('Found {} Rfam accessions'.format(len(rfam_accs)))
    return rfam_accs


def run_rscape(rfam_acc, destination):
    """
    Run R-scape on Rfam seed alignment to get the R-scape/R2R layout.
    """
    rfam_seed = download_rfam_seed(rfam_acc)
    cmd = 'R-scape --outdir {folder} {rfam_seed}'.format(folder=destination, rfam_seed=rfam_seed)
    if not os.path.exists(os.path.join(destination, '{}.out'.format(rfam_acc))):
        os.system(cmd)

    rscape_svg = None
    for svg in glob.glob(os.path.join(destination, '*.svg')):
        if 'R2R.sto.svg' in svg:
            rscape_svg = svg
    if not rscape_svg:
        print('Error: R-scape SVG file not found')
    return rscape_svg


def convert_rscape_svg_to_one_line(rscape_svg, destination):
    """
    Convert R-scape SVG into SVG with 1 line per element.
    """
    output = '{}/temp.svg'.format(destination)
    cmd = r"perl -0777 -pe 's/\n +fill/ fill/g' {rscape_svg} | perl -0777 -pe 's/\n d=/ d=/g' | perl -0777 -pe 's/\n +<tspan/ <tspan/g' | perl -0777 -pe 's/\n<\/text>/<\/text>/g' > {output}".format(rscape_svg=rscape_svg, output=output)
    os.system(cmd)
    return output


def convert_rscape_svg_to_traveler(rscape_one_line_svg, destination):
    """
    Convert R-scape SVG into traveler xml and SVG.
    """
    header = """
    <svg
    xmlns:svg="http://www.w3.org/2000/svg"
    xmlns="http://www.w3.org/2000/svg"
    xmlns:xlink="http://www.w3.org/1999/xlink"
    version="1.0"
    width="1000"
    height="2000"
    xml:space="preserve">
    """
    footer = "</svg>"

    xml_header = "<structure>\n"
    xml_footer = "</structure>"

    with open(rscape_one_line_svg, 'r') as f_in:
        with open(os.path.join(destination, 'traveler-template.svg'), 'w') as f_out:
            with open(os.path.join(destination, 'traveler-template.xml'), 'w') as xml_out:
                f_out.write(header)
                xml_out.write(xml_header)

                for line in f_in.readlines():
                    if '#5c5c5c' in line:
                        continue
                    if 'text1000' in line:
                        continue
                    if '#d7efc5' in line:
                        continue
                    if '<path fill="none" stroke="#000000"' in line:
                        continue
                    if '&apos;' in line:
                        continue
                    if line.startswith('<path'):
                        text, xml = convert_path_to_text(line)
                        f_out.write(text)
                        xml_out.write(xml)
                    elif line.startswith('<text'):
                        xml = convert_text_to_xml(line)
                        if not xml:
                            continue
                        xml_out.write(xml)
                        f_out.write(line)
                    else:
                        continue
                f_out.write(footer)
                xml_out.write(xml_footer)


def rscape2traveler(rfam_acc):
    """
    """
    destination = os.path.join('temp', rfam_acc)
    if not os.path.exists(destination):
        os.makedirs(destination)

    rscape_svg = run_rscape(rfam_acc, destination)
    rscape_one_line_svg = convert_rscape_svg_to_one_line(rscape_svg, destination)
    convert_rscape_svg_to_traveler(rscape_one_line_svg, destination)
    generate_traveler_fasta(rfam_acc)



def generate_2d(rfam_acc, fasta, test):

    destination = 'output/{}'.format(rfam_acc)
    if not os.path.exists(destination):
        os.makedirs(destination)

    if not fasta:
        # use Rfam sequences
        fasta_input = 'rfam/{}.fa'.format(rfam_acc)
        if not os.path.exists(fasta_input):
            url = 'ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/{}.fa.gz'.format(rfam_acc)
            cmd = 'wget -O rfam/{rfam_acc}.fa.gz {url} > rfam/{rfam_acc}.fa.gz && gunzip rfam/{rfam_acc}.fa.gz'.format(url=url, rfam_acc=rfam_acc)
            os.system(cmd)
    else:
        fasta_input = fasta

    if not os.path.exists(fasta_input + '.ssi'):
        cmd = 'esl-sfetch --index {}'.format(fasta_input)
        os.system(cmd)

    # download Rfam covariance model
    if not os.path.exists(os.path.join('data/rfam_cms', rfam_acc + '.cm', )):
        url = 'http://rfam.org/family/{}/cm'.format(rfam_acc)
        cmd = 'wget {url} -O data/rfam_cms/{rfam_acc}.cm'.format(rfam_acc=rfam_acc, url=url)
        os.system(cmd)

    cmd = "grep '>' {} > headers.txt"
    os.system(cmd.format(fasta_input))

    with open('headers.txt', 'r') as f:
        for i, line in enumerate(f.readlines()):
            if test and i > 10:
                continue
            seq_id = line.split(' ', 1)[0].replace('>', '')
            print seq_id

            cmd = 'esl-sfetch %s %s > temp.fasta' % (fasta_input, seq_id)
            os.system(cmd)

            cmd = "cmalign %s temp.fasta > temp.sto" % 'data/rfam_cms/{}.cm'.format(rfam_acc)
            os.system(cmd)

            has_conserved_structure = False
            with open('temp.sto', 'r') as f:
                for line in f.readlines():
                    if line.startswith('#=GC SS_cons '):
                        if '<' in line:
                            has_conserved_structure = True
                        else:
                            print('This RNA does not have conserved structure')
                        break

            if not has_conserved_structure:
                continue

            cmd = 'esl-alimanip --sindi --outformat pfam temp.sto > temp.stk'
            os.system(cmd)

            cmd = 'ali-pfam-sindi2dot-bracket.pl temp.stk > traveler-input.fasta'
            os.system(cmd)

            cmd = ('traveler '
                   '--verbose '
                   '--target-structure traveler-input.fasta '
                   '--template-structure --file-format traveler temp/{rfam_acc}/traveler-template.xml temp/{rfam_acc}/{rfam_acc}-traveler.fasta '
                   '--all output/{rfam_acc}/{seq_id}').format(seq_id=seq_id.replace('/', '-'), rfam_acc=rfam_acc)
            print(cmd)
            os.system(cmd)

            cmd = 'rm -f output/{rfam_acc}/*.xml output/{rfam_acc}/*.ps'.format(rfam_acc=rfam_acc)
            os.system(cmd)



@click.command()
@click.argument('rfam_accession', default='RF00001')
@click.option('--fasta', default=None, help='Sequences to be analysed (by default Rfam hits are analysed)')
@click.option('--test', default=False, is_flag=True, help='Process only the first 10 sequences')
def main(rfam_accession, fasta, test):
    """
    Visualise sequences using the Rfam/R-scape consensus structure as template.

    RFAM_ACCESSION - Rfam family to process (RF00001, RF00002 etc)
    """

    print(rfam_accession)
    if rfam_accession == 'all':
        rfam_accs = get_all_rfam_acc()
    else:
        rfam_accs = [rfam_accession]

    for rfam_acc in rfam_accs:
        rscape2traveler(rfam_acc)
        generate_2d(rfam_acc, fasta, test)


if __name__ == '__main__':
    main()
