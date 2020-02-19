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
import io
import os
import re
import tempfile
import subprocess as sp

from . import config
from . import shared


# these RNAs are better handled by other methods
BLACKLIST = [
    'RF00001', # 5S
    'RF02541', # LSU_rRNA_bacteria
    'RF00177', # SSU_rRNA_bacteria
    'RF01960', # SSU_rRNA_eukarya
    'RF02540', # LSU_rRNA_archaea
    'RF02543', # LSU_rRNA_eukarya
    'RF01959', # SSU_rRNA_archaea
    'RF02542', # SSU_rRNA_microsporidia
    'RF02546', # LSU_trypano_mito
    'RF02545', # SSU_trypano_mito
    'RF00005', # tRNA
    'RF01852', # tRNA-Sec
]


def blacklisted():
    bad = set(BLACKLIST)
    with open(os.path.join(config.RFAM_DATA, 'no_structure.txt')) as raw:
        bad.update(l.strip() for l in raw)
    return bad


def get_traveler_template_xml(rfam_acc):
    filename = os.path.join(config.RFAM_DATA, rfam_acc, 'traveler-template.xml')
    return filename


def get_traveler_fasta(rfam_acc):
    filename = os.path.join(config.RFAM_DATA, rfam_acc,
                            '{}-traveler.fasta'.format(rfam_acc))
    return filename


def get_rfam_cm(rfam_acc):
    return os.path.join(config.RFAM_DATA, rfam_acc, rfam_acc + '.cm')


def get_rfam_cms():
    """
    Download non-blacklisted Rfam covariance models.
    """
    rfam_cm = os.path.join(config.RFAM_DATA, 'Rfam.cm')
    rfam_ids = os.path.join(config.RFAM_DATA, 'rfam_ids.txt')
    if not os.path.exists(rfam_cm):
        cmd = 'wget -O {0}.gz ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz && gunzip {0}.gz'.format(rfam_cm)
        os.system(cmd)
    if not os.path.exists(rfam_cm + '.ssi'):
        os.system('cmfetch --index {}'.format(rfam_cm))
    if not os.path.exists(rfam_ids):
        cmd = "awk '/ACC   RF/ {{print $2}}' {} > {}".format(rfam_cm, rfam_ids)
        os.system(cmd)
    with open(rfam_ids, 'r') as f_in:
        for line in f_in:
            rfam_acc = line.strip()
            if rfam_acc in blacklisted():
                continue
            cm_file = get_rfam_cm(rfam_acc)
            if not os.path.exists(cm_file):
                cmd = 'cmfetch {} {} > {}'.format(rfam_cm, rfam_acc, cm_file)
                os.system(cmd)
            target = os.path.join(os.path.abspath(config.CM_LIBRARY), os.path.basename(cm_file))
            if not os.path.exists(target):
                cmd = 'ln -s {} {}'.format(os.path.abspath(cm_file), target)
                print(cmd)
                os.system(cmd)


def setup(accessions=None):
    precomputed_archive = os.path.join(config.DATA, 'rfam.zip')
    if os.path.exists(precomputed_archive):
        print('Uncompressing precomputed Rfam template library')
        cmd = 'unzip -d {} {}'.format(config.RFAM_DATA, precomputed_archive)
        os.system(cmd)
    else:
        if not accessions:
            accessions = 'all'
        fetch_data(accessions)
    get_rfam_cms()


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
    for seed in glob.glob(os.path.join(config.RFAM_DATA, rfam_acc, '*.R2R.sto')):
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

        with open(get_traveler_fasta(rfam_acc), 'w') as f:
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
        print(line)
        print('convert_path_to_text did not find a match')


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
        print(line)
        print('convert_text_to_xml did not find a match')


def download_rfam_seed(rfam_acc):
    output = os.path.join(config.RFAM_DATA, rfam_acc, '{}.seed'.format(rfam_acc))
    if not os.path.exists(output):
        url = 'http://rfam.org/family/{}/alignment'.format(rfam_acc)
        cmd = 'wget -O {output} {url}'.format(output=output, url=url)
        os.system(cmd)
    return output


def get_all_rfam_acc():
    rfam_accs = []
    family_file = os.path.join(config.RFAM_DATA, 'family.txt')
    if not os.path.exists(family_file):
        cmd = 'wget -O {0}.gz ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/family.txt.gz && gunzip {0}.gz'.format(family_file)
        try:
            sp.check_output([cmd], shell=True)
        except sp.CalledProcessError as e:
            print(e.output)

    with open(family_file, encoding='utf8', errors='ignore') as f:
        for line in f:
            if line.startswith('RF'):
                rfam_acc = line[:7]
                if rfam_acc in BLACKLIST:
                    continue
                rfam_accs.append(rfam_acc)
    print('Found {} Rfam accessions'.format(len(rfam_accs)))
    return rfam_accs


def get_rfam_acc_by_id(rfam_id):
    family_file = os.path.join(config.RFAM_DATA, 'family.txt')
    if not os.path.exists(family_file):
        cmd = 'wget -O {0}.gz ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/family.txt.gz && gunzip {0}.gz'.format(family_file)
        os.system(cmd)
    os.system('cut -f 1,2 {0} > {0}.short'.format(family_file))
    with open(family_file + '.short') as f:
        for line in f:
            rfam_acc, rfam_identifier = line.split()
            if rfam_id == rfam_identifier:
                return rfam_acc


def remove_pseudoknot_from_ss_cons(rfam_seed):
    """
    Some Rfam alignments contain pseudoknot annotations in SS_cons. To avoid
    R-scape displaying them on diagrams, pseudoknots need to be removed before
    running R-scape.
    """
    seed_no_pk = os.path.join(os.path.dirname(rfam_seed), 'nopk-' + os.path.basename(rfam_seed))
    with io.open(rfam_seed, 'r', encoding='latin-1') as f_seed_in:
        with open(seed_no_pk, 'w') as f_seed_out:
            for line in f_seed_in.readlines():
                if line.startswith('#=GC SS_cons '):
                    match = re.match(r'(#=GC SS_cons)(\s+)(.+)', line)
                    no_pk = re.sub(r'\w', '.', match.group(3))
                    # import pdb; pdb.set_trace()
                    f_seed_out.write(''.join([match.group(1), match.group(2), no_pk, '\n']))
                else:
                    f_seed_out.write(line)
    return seed_no_pk


def run_rscape(rfam_acc, destination):
    """
    Run R-scape on Rfam seed alignment to get the R-scape/R2R layout.
    """
    rfam_seed = download_rfam_seed(rfam_acc)
    rfam_seed_no_pk = remove_pseudoknot_from_ss_cons(rfam_seed)
    if not os.path.exists(os.path.join(destination, 'rscape.done')):
        cmd = 'R-scape --outdir {folder} {rfam_seed} && touch {folder}/rscape.done'.format(folder=destination, rfam_seed=rfam_seed_no_pk)
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
    output = os.path.join(destination, 'rscape-one-line.svg')
    cmd = (r"perl -0777 -pe 's/\n +fill/ fill/g' {rscape_svg} | "
           r"perl -0777 -pe 's/\n d=/ d=/g' | "
           r"perl -0777 -pe 's/\n +<tspan/ <tspan/g' | "
           r"perl -0777 -pe 's/\n<\/text>/<\/text>/g' "
           r"> {output}").format(rscape_svg=rscape_svg, output=output)
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
                    if 'pk' in line:
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
    destination = os.path.join(config.RFAM_DATA, rfam_acc)
    if not os.path.exists(destination):
        os.makedirs(destination)

    if os.path.exists(get_traveler_fasta(rfam_acc)) and \
       os.path.exists(get_traveler_template_xml(rfam_acc)):
        return

    rscape_svg = run_rscape(rfam_acc, destination)
    rscape_one_line_svg = convert_rscape_svg_to_one_line(rscape_svg, destination)
    convert_rscape_svg_to_traveler(rscape_one_line_svg, destination)
    generate_traveler_fasta(rfam_acc)


def fetch_data(accessions):
    possible = get_all_rfam_acc()
    if accessions == 'all':
        accessions = possible

    for accession in accessions:
        rscape2traveler(accession)


def download_rfam_cm(rfam_acc):
    rfam_cm = get_rfam_cm(rfam_acc)
    if not os.path.exists(rfam_cm):
        os.system('mkdir -p {}'.format(os.path.join(config.RFAM_DATA, rfam_acc)))
        url = 'http://rfam.org/family/{}/cm'.format(rfam_acc)
        cmd = 'wget {url} -O {rfam_cm}'.format(rfam_cm=rfam_cm, url=url)
        os.system(cmd)
    return rfam_cm


def visualise_rfam(fasta_input, output_folder, seq_id, model_id):
    if not model_id.startswith('RF'):
        rfam_acc = get_rfam_acc_by_id(model_id)
    else:
        rfam_acc = model_id
    rfam_cm = download_rfam_cm(rfam_acc)
    rscape2traveler(rfam_acc)

    temp_fasta = tempfile.NamedTemporaryFile()
    temp_sto = tempfile.NamedTemporaryFile()
    temp_stk = tempfile.NamedTemporaryFile()

    cmd = 'esl-sfetch %s %s > %s' % (fasta_input, seq_id, temp_fasta.name)
    os.system(cmd)

    cmd = "cmalign {rfam_cm} {temp_fasta} > {temp_sto}".format(rfam_cm=rfam_cm,
        temp_fasta=temp_fasta.name, temp_sto=temp_sto.name)
    os.system(cmd)

    has_conserved_structure = False
    with open(temp_sto.name, 'r') as f:
        for line in f.readlines():
            if line.startswith('#=GC SS_cons '):
                if '<' in line:
                    has_conserved_structure = True
                else:
                    print('This RNA does not have a conserved structure')
                break

    if not has_conserved_structure:
        return

    cmd = 'esl-alimanip --sindi --outformat pfam {} > {}'.format(temp_sto.name, temp_stk.name)
    os.system(cmd)

    result_base = os.path.join(output_folder, seq_id.replace('/', '-'))
    input_fasta = os.path.join(output_folder, seq_id + '.fasta')
    cmd = 'ali-pfam-sindi2dot-bracket.pl {} > {}'.format(temp_stk.name, input_fasta)
    os.system(cmd)

    shared.remove_large_insertions(result_base + '.fasta')

    log = result_base + '.log'
    cmd = ('traveler '
           '--verbose '
           '--target-structure {fasta} '
           '--template-structure --file-format traveler {traveler_template_xml} {traveler_fasta} '
           '--all {result_base} '
           '> {log}' ).format(
               fasta=input_fasta,
               result_base=result_base,
               rfam_acc=rfam_acc,
               rfam_data=config.RFAM_DATA,
               traveler_template_xml=get_traveler_template_xml(rfam_acc),
               traveler_fasta=get_traveler_fasta(rfam_acc),
               log=log
            )
    os.system(cmd)

    temp_fasta.close()
    temp_sto.close()
    temp_stk.close()

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


def generate_2d(rfam_acc, output_folder, fasta, test):

    destination = '{}/{}'.format(output_folder, rfam_acc)
    if not os.path.exists(destination):
        os.makedirs(destination)

    if not fasta:
        # use Rfam sequences
        fasta_input = os.path.join(config.RFAM_DATA, rfam_acc, '{}.fa'.format(rfam_acc))
        if not os.path.exists(fasta_input):
            url = 'ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/{}.fa.gz'.format(rfam_acc)
            cmd = 'wget -O {fasta_input}.gz {url} && gunzip {fasta_input}.gz'.format(url=url, fasta_input=fasta_input)
            os.system(cmd)
    else:
        fasta_input = fasta

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
            visualise_rfam(fasta_input, destination, seq_id, rfam_acc)
    os.system('rm headers.txt')


def has_structure(rfam_acc):
    no_structure = []
    with open(os.path.join(config.RFAM_DATA, 'no_structure.txt'), 'r') as f:
        for line in f.readlines():
            no_structure.append(line.strip())
    return rfam_acc not in no_structure
