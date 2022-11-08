import re
import RNA
import json
import requests


def get_r2dt_version_header():
    header = """# R2DT :: visualise RNA secondary structure using templates
# Version 1.3 (October 2022)
# https://github.com/RNAcentral/R2DT
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"""
    return header


def remove_large_insertions_pfam_stk(filename):
    """
    The Pfam Stockholm files can contain 9 or 11 lines depending on whether
    the description line is present.
    """
    MAX_INSERTIONS = 100
    with open(filename, "r") as f:
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
            print("Unexpected number of lines in pfam stk")
            return
        # the tilda and period characters represent insert states in WUSS notation
        match = re.finditer(r"([\.~]{" + str(MAX_INSERTIONS) + ",})", gc_ss_cons)
        if not match:
            return
        for span in match:
            sequence = (
                sequence[: span.start()]
                + "@" * (span.end() - span.start())
                + sequence[span.end() :]
            )
            gr_pp = (
                gr_pp[: span.start()]
                + "@" * (span.end() - span.start())
                + gr_pp[span.end() :]
            )
            gr_ss = (
                gr_ss[: span.start()]
                + "@" * (span.end() - span.start())
                + gr_ss[span.end() :]
            )
            gc_ss_cons = (
                gc_ss_cons[: span.start()]
                + "@" * (span.end() - span.start())
                + gc_ss_cons[span.end() :]
            )
            gc_rf = (
                gc_rf[: span.start()]
                + "@" * (span.end() - span.start())
                + gc_rf[span.end() :]
            )
        if len(lines) == 9:
            lines[3] = re.sub(r"@+", "XXXX", sequence)
            lines[4] = re.sub(r"@+", "XXXX", gr_pp)
            lines[5] = re.sub(r"@+", "~~~~", gr_ss)
            lines[6] = re.sub(r"@+", "~~~~", gc_ss_cons)
            lines[7] = re.sub(r"@+", "xxxx", gc_rf)
        elif len(lines) == 11:
            lines[5] = re.sub(r"@+", "XXXX", sequence)
            lines[6] = re.sub(r"@+", "XXXX", gr_pp)
            lines[7] = re.sub(r"@+", "~~~~", gr_ss)
            lines[8] = re.sub(r"@+", "~~~~", gc_ss_cons)
            lines[9] = re.sub(r"@+", "xxxx", gc_rf)
        with open(filename, "w") as f:
            for line in lines:
                f.write(line)


def get_insertions(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
        if len(lines) == 9:
            sequence = lines[3].split()[1]
        elif len(lines) == 11:
            sequence = lines[5].split()[1]
        trimmed = sequence.replace("-", "")
    match = re.finditer("[a-z]+", trimmed)
    return match


def get_full_constraint(filename):
    constraint = ""
    with open(filename, "r") as f:
        lines = f.readlines()
        if len(lines) == 9:
            gc_ss = lines[5].split()[3]
            sequence = lines[3].split()[1]
        elif len(lines) == 11:
            gc_ss = lines[7].split()[3]
            sequence = lines[5].split()[1]
    deletions = re.finditer("-+", sequence)
    sub_ss = list(gc_ss)
    trimmed = sequence.replace("-", "")
    match = re.finditer("[a-z]+", trimmed)
    adjust = 0
    for span in deletions:
        i = span.start() - adjust
        j = span.end() - adjust
        sub_ss = sub_ss[0:i] + sub_ss[j : len(sub_ss)]
        adjust = adjust + j - i
    match_arr = set()
    for span in match:
        i = span.start()
        j = span.end()
        while i < j:
            match_arr.add(i)
            i += 1
    for i, val in enumerate(sub_ss):
        if i in match_arr:
            constraint += "."
        elif val == "<" or val == "(" or val == "{":
            constraint += "("
        elif val == ">" or val == ")" or val == "}":
            constraint += ")"
        elif val == "," or val == "-" or val == ":":
            constraint += "x"
        else:
            constraint += "x"
    return constraint


def fold_insertions_only(sequence, constraint, filename):
    match = get_insertions(filename)
    list_seq = list(sequence)
    list_con = list(constraint)
    final_list = ""
    for span in match:
        i = span.start()
        j = span.end()
        while i > 1 and (list_con[i - 1] == "x" or list_con[i - 1] == "."):
            i -= 1
        while j < (len(list_con) - 2) and (list_con[j] == "x" or list_con[j] == "."):
            j += 1
        subsequence = "".join(list_seq[i:j])
        constart = list_con[0:i]
        subconstraint = "".join(list_con[i:j])
        formatted_constraint = subconstraint.replace("(", "x").replace(")", "x")
        conend = list_con[j : len(list_con) + 1]
        md = RNA.md()
        md.min_loop_size = 0
        fc = RNA.fold_compound(subsequence, md)
        constraint_options = (
            RNA.CONSTRAINT_DB | RNA.CONSTRAINT_DB_ENFORCE_BP | RNA.CONSTRAINT_DB_DEFAULT
        )
        fc.hc_add_from_db(formatted_constraint, constraint_options)
        (ss, mfe) = fc.mfe()
        if mfe < 99999:
            list_con = constart
            for i, val in enumerate(ss):
                if subconstraint[i] != ".":
                    list_con += subconstraint[i]
                else:
                    list_con += ss[i]
            list_con += conend
        else:
            print("Structure exceeds limits of RNAfold, ignoring constraint")
            return
    final_list = "".join(list_con).replace("x", ".")
    return final_list


def handle_exclusion(exclusion, R2DT_constraint):
    try:
        with open(exclusion, "r") as f:
            exclusion_string = f.read().strip()
            if len(exclusion_string) != len(R2DT_constraint):
                print("Exclusion ignored, not same length as sequence")
                return R2DT_constraint
            elif re.search("[^\.x]", exclusion_string):
                print("Invalid characters in exclusion string, should only be . and x")
                return R2DT_constraint
            else:
                temp_constraint = ""
                for i, val in enumerate(exclusion_string):
                    if exclusion_string[i] == "x":
                        if R2DT_constraint[i] != "." and R2DT_constraint[i] != "x":
                            print(
                                "Invalid exclusion in position "
                                + str(i)
                                + " conflicts with template, constraint ignored"
                            )
                            temp_constraint += R2DT_constraint[i]
                        else:
                            temp_constraint += "x"
                    else:
                        temp_constraint += R2DT_constraint[i]
                return temp_constraint
    except FileNotFoundError:
        print("Constraint file not found")
        return R2DT_constraint


def fold_insertions(input_fasta, exclusion, source, filename, model_id, fold_type):
    with open(input_fasta, "r") as f:
        orig = f.readlines()
    R2DT_constraint = orig[2].strip()
    orig_sequence = orig[1].strip()
    ss = ""
    if (
        fold_type != "insertions_only"
        and fold_type != "full_molecule"
        and fold_type != "all_constraints_enforced"
    ):
        if source == "rfam":
            insertions_only = [
                "rRNA",
                "lncRNA",
                "Cis-regulatory element",
                "ribozyme",
                "CRISPR",
                "antisense",
                "antitoxin",
                "Intron",
            ]
            full_molecule = ["snRNA", "snoRNA", "sRNA", "tRNA", "miRNA"]
            r = requests.get(
                "http://rfam.org/family/{}?content-type=application/json".format(
                    model_id
                )
            )
            if any(x in r.json()["rfam"]["curation"]["type"] for x in full_molecule):
                fold_type = "full_molecule"
            elif any(
                x in r.json()["rfam"]["curation"]["type"] for x in insertions_only
            ):
                fold_type = "insertions_only"
            else:
                fold_type = "insertions_only"
        elif source == "gtrnadb":
            fold_type = "full_molecule"
        else:
            fold_type = "insertions_only"
    if fold_type == "all_constraints_enforced":
        constraint = get_full_constraint(filename)
        final_constraint = ""
        if exclusion:
            final_constraint = handle_exclusion(exclusion, constraint)
        else:
            final_constraint = constraint
        md = RNA.md()
        md.min_loop_size = 0
        fc = RNA.fold_compound(orig_sequence, md)
        constraint_options = (
            RNA.CONSTRAINT_DB | RNA.CONSTRAINT_DB_ENFORCE_BP | RNA.CONSTRAINT_DB_DEFAULT
        )
        fc.hc_add_from_db(final_constraint, constraint_options)
        (ss, mfe) = fc.mfe()
        if mfe > 99999:
            print("Structure exceeds limits of RNAfold, ignoring constraint")
            return
    elif fold_type == "full_molecule":
        if exclusion:
            final_constraint = handle_exclusion(exclusion, R2DT_constraint)
        else:
            final_constraint = R2DT_constraint
        md = RNA.md()
        md.min_loop_size = 0
        fc = RNA.fold_compound(orig_sequence, md)
        constraint_options = (
            RNA.CONSTRAINT_DB | RNA.CONSTRAINT_DB_ENFORCE_BP | RNA.CONSTRAINT_DB_DEFAULT
        )
        fc.hc_add_from_db(final_constraint, constraint_options)
        (ss, mfe) = fc.mfe()
        if mfe > 99999:
            print("Structure exceeds limits of RNAfold, ignoring constraint")
            return
    # Default option
    elif fold_type == "insertions_only":
        constraint = ""
        if exclusion:
            constraint = handle_exclusion(exclusion, R2DT_constraint)
        else:
            constraint = R2DT_constraint
        ss = fold_insertions_only(orig_sequence, constraint, filename)
    if ss:
        with open(input_fasta, "w") as f:
            f.write(orig[0])
        constraint_differences = ""
        for i, val in enumerate(R2DT_constraint):
            if val == ss[i]:
                constraint_differences += "-"
            else:
                constraint_differences += "*"
        with open(input_fasta, "a") as f:
            f.write(orig_sequence)
            f.write("\n" + ss)
            f.write("\n" + constraint_differences)
