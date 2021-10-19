import re
import RNA

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

def fold_insertions(input_fasta, exclusion):
    with open(input_fasta, 'r') as f:
        orig = f.readlines()
    R2DT_constraint = orig[2].strip()
    orig_sequence = orig[1].strip()
    orig_constraint = ''
    if exclusion:
        orig_constraint = handle_exclusion(exclusion, R2DT_constraint)
    else:
        orig_constraint = R2DT_constraint
    md = RNA.md()
    md.min_loop_size = 0
    fc = RNA.fold_compound(orig_sequence, md)
    constraint_options = RNA.CONSTRAINT_DB | RNA.CONSTRAINT_DB_ENFORCE_BP | RNA.CONSTRAINT_DB_DEFAULT
    fc.hc_add_from_db(orig_constraint, constraint_options)
    (ss,mfe) = fc.mfe()
    if mfe < 99999:
        with open(input_fasta, 'w') as f:
            f.write(orig[0])
        constraint_differences = ''
        for i, val in enumerate(R2DT_constraint):
            if val == ss[i]:
                constraint_differences += '-'
            else:
                constraint_differences += '*'
        with open(input_fasta, 'a') as f:
            f.write(orig_sequence)
            f.write("\n" + ss)
            f.write("\n" + constraint_differences)
    else: 
        print('Structure exceeds limits of RNAfold, ignoring constraint')

def handle_exclusion(exclusion, R2DT_constraint):
    try:  
        with open(exclusion, 'r') as f:
            exclusion_string = f.read().strip()
            if len(exclusion_string) != len(R2DT_constraint):
                print('Exclusion ignored, not same length as sequence')
                return R2DT_constraint
            elif re.search('[^\.x]', exclusion_string):
                print('Invalid characters in exclusion string, should only be . and x')
                return R2DT_constraint
            else:
                temp_constraint = ''
                for i,val in enumerate(exclusion_string):
                    if exclusion_string[i] == 'x':
                        if R2DT_constraint[i] != '.':
                            print('Invalid exclusion in position ' + str(i) + ' conflicts with template, constraint ignored')
                            temp_constraint += R2DT_constraint[i]
                        else: 
                            temp_constraint += 'x'
                    else:
                        temp_constraint += R2DT_constraint[i]
                return temp_constraint
    except FileNotFoundError:
        print('Constraint file not found')
        return R2DT_constraint