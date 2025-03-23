#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit
import json
from enum import Enum
from potential_elec_functions import *

def doit(T:int):

        def Point_core_1slater_shell_charge(distance, q_c, z2):
            q_s = charge - q_c
            fit_potential = []
            for i in range(len(distance)):
                if 0 == distance[i]:
                        S = q_s*z2
                else:
                        S =  (q_c / distance[i] + q_s * slater_charge(distance[i], z2))
                fit_potential.append(one_4pi_eps0 * S)
            return fit_potential


        def Point_core_1slater_2slater_shell_charge(distance, q_c, q_s2, z1, z2):
            q_s1 = charge - q_c - q_s2
            fit_potential = []
            for i in range(len(distance)):
                    if 0 == distance[i]:
                            S = q_s1 * z1 + q_s2 * z2/2
                    else:
                            S = (q_c/ distance[i] + q_s1*slater_charge(distance[i], z1)+q_s2*slater2_charge(distance[i], z2) )
                    fit_potential.append(one_4pi_eps0 * S)
            return fit_potential


        def point_core_gaussian_shell_charge(distance, q_c,  z2):
            q_s = charge - q_c
            fit_potential = []
            for i in range(len(distance)):
                    if 0 == distance[i]:
                            S = 2 * q_s * z2 / math.sqrt(math.pi)
                    else:
                            erf2 = math.erf(distance[i] * z2)
                            S = (q_c / distance[i] + q_s * erf2 / distance[i])
                    fit_potential.append(one_4pi_eps0 * S)
            return fit_potential


        def point_core_gaussian_gaussian_shell_charge(distance, q_c, q_s2, z1, z2):
            q_s1 = charge - q_c - q_s2
            fit_potential = []
            for i in range(len(distance)):
                if 0 == distance[i]:
                    S = 2 * (q_s1 * z1 + q_s2 * z2)/math.sqrt(math.pi)
                else:
                    erf1 = math.erf(distance[i] * z1)
                    erf2 = math.erf(distance[i] * z2)
                    #erf = math.erf(distance[i] * z1 * z2 / (z1**2 + z2**2)**0.5)
                    S = (q_c / distance[i] + q_s1 * erf1 / distance[i]+ q_s2 * erf2 / distance[i])
                fit_potential.append(one_4pi_eps0 * S)
            return fit_potential

        def point_core(distance, q_c):
            fit_potential = []
            for i in range(len(distance)):
                if 0 == distance[i]:
                    S = 0
                else:
                    S = (q_c / distance[i])
                fit_potential.append(one_4pi_eps0 * S)
            return fit_potential


        class ChargeModel(Enum):
            POINT_CORE_GAUSSIAN_SHELL = 0
            POINT_CORE_1_SLATER_SHELL = 1
            POINT_CORE_GAUSSIAN_GAUSSIAN_SHELL = 2
            POINT_CORE_1_SLATER_2_SLATER_SHELL = 3
            POINT_CORE = 4


        charge_models = {
            ChargeModel.POINT_CORE_GAUSSIAN_SHELL: "PC+G",
            ChargeModel.POINT_CORE_1_SLATER_SHELL: "PC+1S",
            ChargeModel.POINT_CORE_GAUSSIAN_GAUSSIAN_SHELL: "PC+G+G",
            ChargeModel.POINT_CORE_1_SLATER_2_SLATER_SHELL: "PC+1S+2S",
            ChargeModel.POINT_CORE: "PC",
        }


        parameter_names = {
            ChargeModel.POINT_CORE_GAUSSIAN_SHELL: ["q_c",  "z2"],
            ChargeModel.POINT_CORE_1_SLATER_SHELL: ["q_c", "z2"],
            ChargeModel.POINT_CORE_GAUSSIAN_GAUSSIAN_SHELL: ["q_c", "q_s2", "z1", "z2"],
            ChargeModel.POINT_CORE_1_SLATER_2_SLATER_SHELL: ["q_c",  "q_s2", "z1", "z2"],
            ChargeModel.POINT_CORE: ["q_c"],
        }

        # charge models
        functions = [

            point_core_gaussian_shell_charge,
            Point_core_1slater_shell_charge,
            point_core_gaussian_gaussian_shell_charge,
            Point_core_1slater_2slater_shell_charge,
            point_core
        ]


        with open(f'output_4_{T}.json', 'r') as json_f:
            output_data = json.load(json_f)




        #inputs from json file
        data   = output_data['data']
        charge = None
        initial_guesses = output_data['initial_guesses']
        bounds = output_data['bounds']


        fig, (axes1, axes2, axes3, axes4, axes5, axes6) = plt.subplots(len(data), 1, figsize=(6, 14))
        axes = [axes1, axes2, axes3, axes4, axes5, axes6]

        print(f'Compound & Charge model & charge_core & charge_shell (i) & charge_shell (ii) & zeta_shell (i) & zeta_shell (ii) & RMSE')
        all_params = {}

        for i, (compound, (distance_data, potential_data)) in enumerate(data.items()):

            params = {}
            popts_compound = []
            pcovs_compound = []

            if compound=='F' or compound=='Cl' or compound=='Br' or compound=='I':
                charge=-1
            elif compound=='Li' or compound=='Na' or compound=='K':
                charge=1
            else:
                sys.exit("Unknown compound %s" % compound)
            label = f' {compound}'
            if charge == -1:
                    label += "-"
            else:
                    label += "+"

            # Modify data to just fit to the difference, but not for the set starting at zero
            original_charge = 0
            if distance_data[0] > 0:
                    original_charge = charge
                    for d in range(len(distance_data)):
                            potential_data[d] -= one_4pi_eps0*charge/distance_data[d]
                    charge = 0
            for func_index, func in enumerate(functions):
                if func_index == 4:
                        continue
                try:
                        if func_index < len(initial_guesses[compound]):
                                myp0 = initial_guesses[compound][func_index]
                                mybounds = bounds[compound][func_index]
                        else:
                                myp0     = [ charge ]
                                mybounds = [ [ charge - 0.5], [ charge + 0.5 ] ]
                        popt, pcov = curve_fit(func, distance_data, potential_data,
                                               p0=myp0,
                                               bounds=mybounds,
                                               maxfev=10000)
                except:
                        print("Failed to do the curve_fit for compound %s func %d" % ( compound, func_index ))
                        continue

                charge_model = charge_models[ChargeModel(func_index)]
                param_names = parameter_names[ChargeModel(func_index)]
                popts_compound.append(popt)
                pcovs_compound.append(pcov)
                charge_model_compound = func(distance_data, *popt)


                if func_index == 0 or func_index == 1:
                    q_c_opt, z2_opt = popt
                    q_s_opt = charge - q_c_opt
                    params[f"q_c_{func_index}"] = q_c_opt
                    params[f"q_s_{func_index}"] = q_s_opt
                    params[f"z2_{func_index}"]  = z2_opt

                    rmse = np.sqrt(np.mean((np.array(charge_model_compound) - np.array(potential_data))**2))
                    print(f"{label} & {charge_model} &  {original_charge+q_c_opt:.3f} &  {q_s_opt:.3f} & - &  {z2_opt:.3f} & -  & {rmse:.5f} \\\\")

                elif func_index == 2 or func_index == 3:
                    q_c_opt, q_s2_opt, z1_opt, z2_opt = popt
                    q_s1_opt = charge - q_c_opt - q_s2_opt

                    params[f"q_c_{func_index}"]  = q_c_opt + original_charge
                    params[f"q_s1_{func_index}"] = q_s1_opt
                    params[f"q_s2_{func_index}"] = q_s2_opt
                    params[f"z1_{func_index}"]   = z1_opt
                    params[f"z2_{func_index}"]   = z2_opt

                    rmse = np.sqrt(np.mean((np.array(charge_model_compound) - np.array(potential_data))**2))
                    print(f"{label} & {charge_model} & {original_charge+q_c_opt:.3f} & {q_s1_opt:.3f} &  {q_s2_opt:.3f} & {z1_opt:.3f} & {z2_opt:.3f} & {rmse:.5f} \\\\")
                elif func_index == 4:
                    [q_c_opt] = popt
                    params[f"q_c_{func_index}"]  = q_c_opt + original_charge
                    rmse = np.sqrt(np.mean((np.array(charge_model_compound) - np.array(potential_data))**2))
                    print(f"{label} & {charge_model} & {original_charge+q_c_opt:.3f} &  &  & &  & {rmse:.3f} \\\\")

                axes[i].plot(np.array(distance_data), np.array(charge_model_compound)-np.array(potential_data), label=f'{charge_model}')
                axes[i].text(.82, .89, label, transform=axes[i].transAxes,  va='top', fontsize=18)
            all_params[compound] = params
        axes[2].set_ylabel('Residual electrostatic potential (kJ/mol e)', fontsize=18)
        axes[5].set_xlabel(f'Distance ($\AA$)', fontsize=18)
        axes[5].tick_params(axis='x', labelsize=14)
        for ix in range(6):
            axes[ix].tick_params(axis='y', labelsize=14)
        print()

        json_f_path = f"params_4_{T}.json"
        with open(json_f_path, "w") as json_f:
             json.dump(all_params, json_f, indent=4)
             json_f.write('\n')

        axes1.legend(bbox_to_anchor=(.05, .95), loc='upper left',fontsize=14)
        plt.tight_layout()
        plt.subplots_adjust(hspace=0)

        plt.savefig('Fit-Alkali_Halides-%d.pdf' % T)
        plt.show()

if __name__ == "__main__":
    for T in [ 100, 10 ]:
            doit(T)
