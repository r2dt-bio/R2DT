/* Example sequences for the 2D form (subset of https://r2dt.bio). */
window.R2DT_WS_EXAMPLES = {
  draw: [
    {
      id: "3skz",
      label: "3SKZ_B",
      note: "with structure → templatefree",
      layout: "templatefree",
      fasta:
        ">3SKZ_B\n" +
        "GGCCUUAUACAGGGUAGCAUAAUGGGCUACUGACCCCGCCUUCAAACCUAUUUGGAGACUAUAAGGUC\n" +
        ".((((((((A..((((((.....BB))))))(.....a)(((((((bb..)))))))..)))))))).\n",
    },
    {
      id: "sbox",
      label: "S box leader",
      note: "SAM riboswitch",
      layout: "auto",
      fasta:
        ">S_box_leader\n" +
        "ACCTTATTTTGAGAAGCTGAGGGATTTGGCCCATAGAAGCTTCAGCAACCGACTTTAAATAGCACGGTGCTAATACCAACGAGCAACTCGAATGATAAGTA\n",
    },
    {
      id: "trna",
      label: "TRT-TGT2-1",
      note: "tRNA",
      layout: "auto",
      fasta:
        ">TRT-TGT2-1\n" +
        "GGCTCCATAGCTCAGTGGTTAGAGCACTGGTCTTGTAAACCAGGGGTCGCGAGTTCGATCCTCGCTGGGGCCT\n",
    },
    {
      id: "u1",
      label: "RNVU1-1",
      note: "U1 snRNA",
      layout: "auto",
      fasta:
        ">RNVU1-1\n" +
        "AUACUUACCUGGCAGGGGAGAUACGAUGAUCACGAAGGUGGUUUUCCCAGGGAGAGGCUUAUCCAUUGCACUCCGGAUGUGCUGACCCCUGCCGUUUCCCCAAAUGUGGGAAACUCGACUGCAUAAUUUGUGGUAGUGGGGGACUGCGUUCGCGCUGUCCUCUG\n",
    },
    {
      id: "multi",
      label: "2× SAM riboswitch",
      note: "multi-FASTA → 2 jobs",
      layout: "auto",
      fasta:
        ">URS000053CEAC_SAM\n" +
        "CTCTTATCGAGAGTTGGGCGAGGGATTTGGCCTTTTGACCCCAAAAGCAACCGACCGTAATTCCATTGTGAAATGGGGCGCATTTTTTTCGCGCCGAGACGCTGGTCTCTTAAGGCACGGTGCTAATTCCATTCAGATCTGATCTGAGAGATAAGAG\n" +
        ">URS00001D0AD3_SAM\n" +
        "TCCTTATCAAGAGAGGTGGAGGGACTGGCCCTGCGATACCCGGCAACCGCTGTTTAACAGAATGGTGCTAAATCCTTTAGAGCAATGATTGCTCTTGAAGATAAGGT\n",
    },
  ],
  align: [
    {
      id: "RF00162",
      label: "RF00162 SAM",
      note: "seed without cov_* → R-scape then draw",
    },
  ],
};
