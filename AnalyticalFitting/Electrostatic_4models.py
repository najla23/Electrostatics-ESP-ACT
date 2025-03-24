#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import json
import math
import os
from potential_elec_functions import *

figs = "Figures"

def main(T:int):
    with open(f'params_4_{T}.json', 'r') as json_f:
        params = json.load(json_f)


    cations = ["Li", "Na", "K"]
    anions = ["F", "Cl", "Br"]
    dir = 'SAPT_Alkali_Halides'

    files = [
        os.path.join(dir, f'distances_Electrostatics-{cation.lower()}{anion.lower()}.txt')
        for cation in cations
        for anion in anions
    ]

    data2 = [np.loadtxt(file) for file in files]
    distances = np.arange(1, 4.6, 0.1)

    function_values = {}

    func_index_to_name = {
        0: "P+G",
        1: "P+1S",
        2: "P+G+G",
        3: "P+1S+2S"
    }

    func_index_to_function = {
        0: Point_core_gaussian_shell,
        1: Point_core_1slater_shell,
        2: Point_core_2gaussian_shell,
        3: Point_core_1slater_2slater_shell,
    }


    for file, dataset in zip(files, data2):
        fig, axes = plt.subplots(4, 1, figsize=(3, 9))

        for i in range(4):
            x, y = dataset[:, 0], dataset[:, 1]
            sorted_indices = np.argsort(x)
            x1 = x[sorted_indices]
            y1 = y[sorted_indices]

            for cation in cations:
                for anion in anions:
                    if file == f'SAPT_Alkali_Halides/distances_Electrostatics-{cation.lower()}{anion.lower()}.txt':
                        func_index = i
                        function_name = f"{func_index_to_name[func_index]}_{cation}{anion}"
                        function_values[function_name] = None

                        if func_index == 0 or func_index ==1:
                            if  func_index == 0:
                                CM= Point_core_gaussian_shell
                            if func_index == 1:
                                CM= Point_core_1slater_shell

                            function_values[function_name] = CM(
                                distances, params[f"{cation}"][f"q_c_{func_index}"],
                                params[f"{cation}"][f"q_s_{func_index}"],
                                params[f"{anion}"][f"q_c_{func_index}"],
                                params[f"{anion}"][f"q_s_{func_index}"],
                                params[f"{cation}"][f"z2_{func_index}"],
                                params[f"{anion}"][f"z2_{func_index}"]
                            )

                        elif func_index == 2 or func_index == 3:
                            if  func_index == 2:
                                CM= Point_core_2gaussian_shell
                            if  func_index == 3:
                                CM= Point_core_1slater_2slater_shell

                            function_values[function_name] = CM(
                                distances, 
                                params[f"{cation}"][f"q_c_{func_index}"],
                                params[f"{cation}"][f"q_s1_{func_index}"],
                                params[f"{cation}"][f"q_s2_{func_index}"],
                                params[f"{anion}"][f"q_c_{func_index}"],
                                params[f"{anion}"][f"q_s1_{func_index}"],
                                params[f"{anion}"][f"q_s2_{func_index}"],
                                params[f"{cation}"][f"z1_{func_index}"],
                                params[f"{anion}"][f"z1_{func_index}"],
                                params[f"{cation}"][f"z2_{func_index}"],
                                params[f"{anion}"][f"z2_{func_index}"]
                            )



                        axes[i].plot(x1, y1, label=f"SAPT_{cation}{anion}", color='r')
                        axes[i].plot(distances, function_values[function_name], color='black', label=function_name)
                        axes[i].set_xlabel('Distance ($\AA$)', fontsize=12)
                        axes[i].set_ylabel('Electrostatic energies (kJ/mol)', fontsize=8)
                        axes[i].tick_params(axis='x', labelsize=8)
                        axes[i].tick_params(axis='y', labelsize=8)
                        axes[i].legend(fontsize=12)
                        plt.tight_layout()
                        plt.subplots_adjust(hspace=0)
                        plt.savefig(f'{figs}/SAPT_{cation}{anion}_{T}.pdf')
#       plt.show()

def print_tex():
    with open("ionpair_energies.tex", "w") as outf:
        for cat in [ "Li", "Na", "K" ]:
            for an in [ "F", "Cl", "Br" ]:
                ionpair = cat+an
                outf.write("""
\\begin{figure}[htb!]
\\centering
\\begin{minipage}{0.5\\textwidth}
\\centering
\\includegraphics[width=0.9\\linewidth]{Figures/SAPT_%s_10.pdf}
\\end{minipage}%%
\\begin{minipage}{0.5\\textwidth}
\\centering
\\includegraphics[width=0.9\\linewidth]{Figures/SAPT_%s_100.pdf}
\\end{minipage}
\\caption{Electrostatic energies from SAPT0 with the aug-cc-pVTZ basis set and different charge models based on fitting the ESP from 2.0 to 4.5 {\\AA} (left) and 0.0 to 4.5 {\\AA} (right) for %s. Note different units on the y-axis.}
\\label{fig:pot_%s}
\\end{figure}
""" % ( ionpair, ionpair, ionpair, ionpair ) )

if __name__ == "__main__":
    print_tex()
    os.makedirs(figs, exist_ok=True)
    for T in [ 10, 100 ]:
        main(T)
