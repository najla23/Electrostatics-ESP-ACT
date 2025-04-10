#!/usr/bin/env python3

import os

def get_train()->list:
    mylist = []
    with open("../Selection/ac-train.dat", "r") as inf:
        for line in inf:
            ww = line.strip().split("|")
            mylist.append(ww[0])
    return mylist
        
def run_one(filenm:str)->dict:
    tempf = "rmsd.txt"
    os.system("grep 'COULOMB RMSD' %s > %s" % ( filenm, tempf ) )
    
    mydict = {}
    with open(tempf, "r") as inf:
        for line in inf:
            words = line.strip().split()
            if len(words) == 14:
                mydict[words[13]] = { "RMSD": words[2], "N": words[12] }
    return mydict
    
if __name__ == "__main__":
    basefn = "train-Inter6.log"
    alldata = {}
    for mydir in [ "PC", "GC", "GC+PGV", "PC+GVS" ]:
        alldata[mydir] = run_one("Elec/" + mydir + "/" + basefn)
    
    train = get_train()
    texfn = "rmsdtable.tex"
    with open(texfn, "w") as outf:
        outf.write("\\begin{landscape}\n")
        outf.write("\\begin{longtable}{lccccc}\n")
        outf.write("\\caption{Root mean square deviation (kJ/mol) from SAPT electrostatics per compound dimer for the different ACT models. N is the number of conformations of each dimer used. Compound dimers used in training are printed in {\\bf bold font}.}\\\\\n")
        outf.write("\\hline\n")
        outf.write("Dimer & N ")
        for md in alldata.keys():
            outf.write(" & %s" % md)
        outf.write("\\\\\n")
        outf.write("\\hline\n")
        for dimer in sorted(alldata["PC"].keys()):
            ddd = dimer.replace("#", "-")
            if dimer in train:
                ddd = ("{\\bf %s}" % ddd)
            outf.write("%s & %s " % (ddd, alldata["PC"][dimer]["N"] ) )
            for md in alldata.keys():
                outf.write(" & %s" % alldata[md][dimer]["RMSD"])
            outf.write("\\\\\n")
        outf.write("\\hline\n")
        outf.write("\\end{longtable}\n")
        outf.write("\\end{landscape}\n")
    
        
