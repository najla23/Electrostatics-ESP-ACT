#!/usr/bin/env python3

import os, sys

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

if len(sys.argv) != 3:
    sys.exit("Usage: %s input.tex output.tex" % sys.argv[0])
    
if os.path.exists(sys.argv[2]):
    sys.exit("Will not write over existing file %s please remove manually" % sys.argv[2])
    
with open(sys.argv[2], "w") as outf:
    with open(sys.argv[1], "r") as inf:
        for line in inf:
            lll = line.strip()
            words = lll.split("&")
            if len(words) == 3 and words[0].strip() in table.keys():
                outf.write("%s & \\verb^%s^ & %s\n" % ( table[words[0].strip()],
                                                        words[1].strip(), words[2] ))
            else:
                outf.write(line)

                
