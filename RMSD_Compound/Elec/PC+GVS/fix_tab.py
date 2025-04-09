#!/usr/bin/env python3

import os

table = { "a1dexp": "Induction correction A",
          "bdexp": "Induction correction b",
          "alpha": "Polarizability $\\alpha$",
          "charge": "Charge",
          "chi": "Electronegativity $\\chi$",
          "delta_chi": "Bond electronegativity $\\Delta\\chi$",
          "delta_eta": "Bond hardness $\\Delta\\eta$",
          "eta": "Hardness $\\eta$",
          "zeta": "Gaussian distribution width $\\zeta$",
          "vs3sa": "Virtual site position $a$" }

outt = "output_table_PC+GVS.tex"
with open(outt, "w") as outf:
    intt = "output_table_SC-Water-Ion_PC+GVS.tex"
    with open(intt, "r") as inf:
        for line in inf:
            lll = line.strip()
            words = lll.split("&")
            if len(words) == 3 and words[0].strip() in table.keys():
                outf.write("%s & {\\verb %s} & %s\n" % ( table[words[0].strip()],
                                                        words[1].strip(), words[2] ))
            else:
                outf.write(line)

                
