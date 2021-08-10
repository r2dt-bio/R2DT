import re


def get_r2dt_version_header():
    header = """# R2DT :: visualise RNA secondary structure using templates
# Version 1.2 (August 10, 2021)
# https://github.com/RNAcentral/R2DT
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"""
    return header


def remove_large_insertions_pfam_stk(filename):
    """
    The Pfam Stockholm files can contain 9 or 11 lines depending on whether
    the description line is present.
    """
    MAX_INSERTIONS = 100
    with open(filename, 'r') as f:
        lines = f.readlines()
        if len(lines) == 9:
            sequence = lines[3]
            gr_pp = lines[4]
            gr_ss = lines[5]
            gc_ss_cons = lines[6]
            gc_rf = lines[7]
        elif len(lines) == 11:
            sequence = lines[5]
            gr_pp = lines[6]
            gr_ss = lines[7]
            gc_ss_cons = lines[8]
            gc_rf = lines[9]
        else:
            print('Unexpected number of lines in pfam stk')
            return
        # the tilda and period characters represent insert states in WUSS notation
        match = re.finditer(r'([\.~]{' + str(MAX_INSERTIONS) + ',})', gc_ss_cons)
        if not match:
            return
        for span in match:
            sequence = sequence[:span.start()] + '@' * (span.end() - span.start()) + sequence[span.end():]
            gr_pp = gr_pp[:span.start()] + '@' * (span.end() - span.start()) + gr_pp[span.end():]
            gr_ss = gr_ss[:span.start()] + '@' * (span.end() - span.start()) + gr_ss[span.end():]
            gc_ss_cons = gc_ss_cons[:span.start()] + '@' * (span.end() - span.start()) + gc_ss_cons[span.end():]
            gc_rf = gc_rf[:span.start()] + '@' * (span.end() - span.start()) + gc_rf[span.end():]
        if len(lines) == 9:
            lines[3] = re.sub(r'@+', 'XXXX', sequence)
            lines[4] = re.sub(r'@+', 'XXXX', gr_pp)
            lines[5] = re.sub(r'@+', '~~~~', gr_ss)
            lines[6] = re.sub(r'@+', '~~~~', gc_ss_cons)
            lines[7] = re.sub(r'@+', 'xxxx', gc_rf)
        elif len(lines) == 11:
            lines[5] = re.sub(r'@+', 'XXXX', sequence)
            lines[6] = re.sub(r'@+', 'XXXX', gr_pp)
            lines[7] = re.sub(r'@+', '~~~~', gr_ss)
            lines[8] = re.sub(r'@+', '~~~~', gc_ss_cons)
            lines[9] = re.sub(r'@+', 'xxxx', gc_rf)
        with open(filename, 'w') as f:
            for line in lines:
                f.write(line)
