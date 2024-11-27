"""
Reads INFERNAL mapping which contains one line for sequence and one for the structure. The sequence contains
in upper case nucleotides which are mapped from template to target. Lower case letters denote nucleotides
which are inserted into target respective to template and "-" denote positions which were deleted with respect
to template. The structure then contains the corresponding secondary structure, i.e. secondary structure of template
with added .~ characters corresponding to insertions (these will match the lowercase positions). Removing the positions
corresponding to "-" nucleotides will result to target structure while removing ".~" characters results in the
template structure.

Parsing the INFERNAL input results in two structures that contain sequence and structure of both target and
template (the template sequence has only "-" at positions which were deleted since information about the template is not present in
the input file) and mapping list which gives mapping from target/template to the position in the infernal
mapping (since this contains both insertions and deletions).

The structures are then processed into the pre-order sorted list of nodes of the target and template trees. Each node
contains the INFERNAL-indexed index/indexes of the respective nucleotides from the tree-based structure representation
(i.e. base-pairs represented by inner nodes and non-paired nucleotides by leafs [except of some special cases when
hairpin ends with a base pair]).

The pre-ordered lists are then scanned to identify indexes of nodes in the template tree which need to be matched/inserted/deleted
to arrive to the target tree. This mapping is then sent to output and can be imported by Traveler.

Notes on the '.afa' files by Eric Nawrocki:

The aligned sequence and SS_cons strings will be the same length.

In the aligned sequences, there are three types of columns:
- lowercase nucleotides are insertions relative to the template

- uppercase nucleotides align to a template position

- dashes (-) align to a template position but indicate where a
  template position is missing (has been deleted)

In the aligned secondary structures, there are four types of columns:
- The four characters in the set :_,- represent positions that align
  to a single stranded position in the template (these will always
  match up with either an uppercase nucleotide or a dash (-) in the
  aligned sequence)

- The four characters in the set {<([ are left halves of basepairs
  (these will always match up with either an uppercase nucleotide or a
  dash (-) in the aligned sequence)

- The four characters in the set }>)] are right halves of basepairs
  (these will always match up with either an uppercase nucleotide or a
  dash (-) in the aligned sequence)

- The two characters in the set .~ are insertions relative to the
  template (these will always match up with a lowercase nucleotide in
  the aligned sequence)

If we strip out the characters in the set .~ from the aligned
secondary structures, we will always get the full length template
secondary structure.
"""
import argparse
import copy
import logging
import sys
from typing import List, Tuple


class SequenceStructureMapping:
    def __init__(self):
        self.sq = ""  # sequence
        self.str = ""  # structure
        self.m = (
            []
        )  # mapping of position in tgt or tmp to the infernal string (which comprises both tgt and tmp information)


def error_exit(message):
    logging.error(message)
    sys.exit()


def read_structures(f) -> List[SequenceStructureMapping]:
    f.readline()
    sq = f.readline().strip()
    f.readline()
    str_tgt = f.readline().strip()
    f.readline()
    str_temp = f.readline().strip()
    f.readline()
    str_temp_orig = (
        f.readline().strip()
    )  # the input can contain also original secondary structure before depairing carried out by the R2DT pipeline.
    # Depairing happens when there are non-canonical base pairs in the template - these are then removed from the template secondary structure
    # but then it can happen that the sequence and template structures are not valid with respect to the rules described in the header of this file.
    # However, to Traveler is then passed the original template, i.e. the one corresponding to the "orig", where the non-canonical base pairs exist.
    assert len(sq) == len(str_temp) == len(str_tgt)

    str_temp_used = str_temp_orig if str_temp_orig else str_temp

    tgt_s_s_m = SequenceStructureMapping()
    tmp_s_s_m = (
        SequenceStructureMapping()
    )  # as for sequence, template contains dashes at positions which were deleted in target
    # tmp_full_s_s_m = SequenceStructureMapping()

    # Non-bracket and non-alpha symbols correspond to gaps in sequence
    tgt_s_s_m.sq = sq.replace("-", "")
    # . and ~ symbols correspond to insertions with respect to template
    # tmp_full_s_s_m.str = copy.deepcopy(tmp_s_s_m.str)

    tmp_s_s_m.str = ""
    for nt, str_character in zip(sq, str_temp_used):
        # if nt is lowercase AND str_character is . or ~, then it is an insertion, so skip it
        if "a" <= nt <= "z" and str_character in ".~":
            continue
        tmp_s_s_m.str += str_character

    for i, letter in enumerate(sq):
        if letter != "-":
            tgt_s_s_m.str += str_temp[i]
            tgt_s_s_m.m.append(i)
        if not "a" <= letter <= "z" and letter != "X":
            tmp_s_s_m.sq += letter
            tmp_s_s_m.m.append(i)
            # tmp_full_s_s_m.m.append(i)

    # return [sq, tgt_s_s_m, tmp_s_s_m, tmp_full_s_s_m]
    return [sq, tgt_s_s_m, tmp_s_s_m]


def get_ix(ix_cap: int, str_ix: List[List[int]]) -> int:
    """
    Finds in the str_ix, which is a list of non-paired (list length = 1) or base-paired(list length = 2),
    first index with residue having first residue having value > ix_cap.
    :param ix_cap:
    :param str_ix:
    :return:
    """
    max_ix = 0
    for i, p in enumerate(str_ix):
        if p[0] > ix_cap:
            return i

    assert False


def read_str_ix(s_s_m: SequenceStructureMapping, affected_ix: List[int]):
    """
    Gives the list of paired positions
    :param s_s_m:
    :return:
    """
    str_ix = []
    stack_normal = []  # (
    stack_curly = []  # {
    stack_arrow = []  # <
    stack_square = []  # [
    str = s_s_m.str
    for ix in range(len(str)):
        if str[ix] not in "{<([}>)]" or s_s_m.m[ix] in affected_ix:
            str_ix.append([ix])
        elif str[ix] in "{<([":
            if str[ix] == "(":
                stack = stack_normal
            elif str[ix] == "{":
                stack = stack_curly
            elif str[ix] == "<":
                stack = stack_arrow
            elif str[ix] == "[":
                stack = stack_square
            stack.append(ix)
        else:
            if str[ix] == ")":
                stack = stack_normal
            elif str[ix] == "}":
                stack = stack_curly
            elif str[ix] == ">":
                stack = stack_arrow
            elif str[ix] == "]":
                stack = stack_square
            ix1 = stack.pop()

            # ix_str_ix = get_ix(ix_cap=ix1, str_ix=str_ix)
            # str_ix.insert(ix_str_ix, [ix1, ix])
            str_ix.append([ix1, ix])  # the list is sorted in the post order
            # str_ix = [[ix1, ix]] + str_ix

    return str_ix


def find_in_mapping(ixs: List[int], str_ix: List[List[int]]) -> int:
    """
    Find index of the node which contains the ixs (which can be either index of a residue or indexes of two
    residues if it corresponds to a node holding a base pair.
    :param ixs:
    :param str_ix:
    :return:
    """
    for i in range(len(str_ix)):
        ixs1 = str_ix[i]
        if len(ixs) == 1 and len(ixs1) == 1:
            if ixs[0] == ixs1[0]:
                return i
        elif len(ixs) == 2 and len(ixs1) == 2:
            if ixs[0] == ixs1[0] and ixs[1] == ixs1[1]:
                assert len(ixs) == len(ixs1)
                return i
    return -1


def get_mapping(
    tmp_str_ix: List[List[int]], tgt_str_ix: List[List[int]]
) -> List[Tuple[int, int]]:
    """
    Goes through the nodes of template tree and finds where each of the node is placed in the
    target tree. This uses the fact that residues are alignment-indexed, i.e. aligned template and target
    residues have the same index.

    :param tmp_str_ix: Template tree in the form of list of nodes.
    Each node is either a single-member or two-member array corresponding either to a leave or base pair.
    :param tgt_str_ix: Target tree.
    :return:
    """
    m = []
    tgt_mapped = []
    for i_tmp in range(len(tmp_str_ix)):
        i_tgt = find_in_mapping(tmp_str_ix[i_tmp], tgt_str_ix)
        assert i_tgt == -1 or len(tmp_str_ix[i_tmp]) == len(tgt_str_ix[i_tgt])
        m.append((i_tmp + 1, i_tgt + 1))
        tgt_mapped.append(i_tgt)

    for i_tgt in set(range(len(tgt_str_ix))).difference(tgt_mapped):
        m.append((0, i_tgt + 1))

    return m


def get_distance(m: List[Tuple[int]]) -> int:
    d = 0
    for i1, i2 in m:
        if i1 == 0 or i2 == 0:
            d += 1
    return d


def get_affected_positions_in_tmp(
    sq_aln, s_s_m: SequenceStructureMapping, str_ix: List[List[int]]
) -> List[int]:
    """
    Identifies list of positions in the aligned sequence (containing insertions to both template and target) corresponding
    to basepairs where at least one of the residues was removed in template and therefore the pairing fell apart.
    :param sq_aln:
    :param s_s_m:
    :param str_ix:
    :return:
    """

    affected: List[int] = []

    for p in str_ix:
        if len(p) == 1:
            continue

        ix_aln1 = s_s_m.m[p[0]]
        ix_aln2 = s_s_m.m[p[1]]
        if sq_aln[ix_aln1] == "-" or sq_aln[ix_aln2] == "-":
            affected.append(ix_aln1)
            affected.append(ix_aln2)

    return affected


def convert_to_aln_ix(
    str_ix: List[List[int]], s_s_m: SequenceStructureMapping
) -> List[List[int]]:
    """
    Convers strucutre-indexed positions to alignment-indexed position based on structure-alignment position
    mapping available in the provided sequence-structure mapping.

    :param str_ix:
    :param s_s_m:
    :return:
    """
    aln_ix: List[List[int]] = []
    for p in str_ix:
        pp = []
        for ix in p:
            pp.append(s_s_m.m[ix])
        aln_ix.append(pp)
    return aln_ix


def main():
    with open(args.input, "r") as fr:
        sq_aln, tgt_s_s_m, tmp_s_s_m = read_structures(fr)
        # tmp_str_ix = read_str_ix(tmp_full_s_s_m, [])
        tmp_str_ix = read_str_ix(tmp_s_s_m, [])
        # print(f"Template size: {len(tmp_str_ix)}")
        # tmp_affected = get_affected_positions_in_tmp(sq_aln=sq_aln, s_s_m=tmp_full_s_s_m, str_ix=tmp_str_ix)
        tmp_affected = get_affected_positions_in_tmp(
            sq_aln=sq_aln, s_s_m=tmp_s_s_m, str_ix=tmp_str_ix
        )
        tgt_str_ix = read_str_ix(tgt_s_s_m, tmp_affected)
        # print(f"Target size: {len([s for s in tgt_str_ix if s != '-'])}")

        # tmp_str_ix_aln = convert_to_aln_ix(tmp_str_ix, tmp_full_s_s_m)
        tmp_str_ix_aln = convert_to_aln_ix(tmp_str_ix, tmp_s_s_m)
        tgt_str_ix_aln = convert_to_aln_ix(tgt_str_ix, tgt_s_s_m)

        m = get_mapping(tmp_str_ix_aln, tgt_str_ix_aln)
        with (sys.stdout if args.output is None else open(args.output, "w")) as fw:
            fw.write("DISTANCE: {}\n".format(get_distance(m)))
            for i1, i2 in m:
                fw.write("{} {}\n".format(i1, i2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        metavar="FILE",
        help="Infernal FASTA-like .afa file with the alignment.",
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="FILE",
        help="Output file name for the Traveler-formatted mapping. "
        "If non entered, the standard output will be used.",
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s [%(levelname)s] %(module)s - %(message)s",
        datefmt="%H:%M:%S",
    )

    main()
