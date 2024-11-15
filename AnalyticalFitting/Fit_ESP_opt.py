#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit
import json
from enum import Enum




def slater2_charge(distance, alpha):
    return 1/distance - (6 + 9*distance*alpha + 6*pow(distance, 2)*pow(alpha, 2) + 2*pow(distance, 3)*pow(alpha, 3)) / (6*math.exp(2*distance*alpha)*distance)


def slater_charge(distance, alpha):
    return (1/distance) - (1 + distance * alpha) / (math.exp(2 * distance * alpha) * distance)


def Point_core_1slater_shell_charge(distance, q_c, z2):
    q_s = charge - q_c
    fit_potential = []
    for i in range(len(distance)):
            fit_potential.append(K * (q_c / distance[i] + q_s * slater_charge(distance[i], z2)))
    return fit_potential


def Point_core_2slater_shell_charge(distance, q_c, z2):
    q_s = charge - q_c
    fit_potential = []
    for i in range(len(distance)):
            erf2 = math.erf(distance[i] * z2)
            fit_potential.append(K * (q_c / distance[i] + q_s * slater2_charge(distance[i], z2)))
    return fit_potential


def Point_core_1slater_2slater_shell_charge(distance, q_c, q_s2, z1, z2):
    q_s1 = charge - q_c - q_s2
    fit_potential = []
    for i in range(len(distance)):
            fit_potential.append(K *(q_c/ distance[i] + q_s1*slater_charge(distance[i], z1)+q_s2*slater2_charge(distance[i], z2) ))
    return fit_potential


def slater_core_2slater_shell_charge(distance, q_c, z1, z2):
    q_s = charge - q_c
    fit_potential = []
    for i in range(len(distance)):
            fit_potential.append(K * (q_c * slater_charge(distance[i], z1) + q_s * slater2_charge(distance[i], z2)))
    return fit_potential


def slater1_core_slater1_shell_charge(distance, q_c, z1, z2):
    q_s = charge - q_c
    fit_potential = []
    for i in range(len(distance)):
            fit_potential.append(K * (q_c  * slater_charge(distance[i], z1) + q_s * slater_charge(distance[i], z2)))
    return fit_potential


def point_core_shell_charge(distance, q_c):
    q_c=charge
    fit_potential = []
    for i in range(len(distance)):
            fit_potential.append(K * q_c / distance[i])
    return fit_potential


def point_core_gaussian_shell_charge(distance, q_c,  z2):
    q_s = charge - q_c
    fit_potential = []
    for i in range(len(distance)):
            erf2 = math.erf(distance[i] * z2)
            fit_potential.append(K * (q_c / distance[i] + q_s * erf2 / distance[i]))
    return fit_potential


def gaussian_core_gaussian_shell_charge(distance, q_c,  z1, z2):
    q_s = charge - q_c
    fit_potential = []
    for i in range(len(distance)):
            erf1 = math.erf(distance[i] * z1)
            erf2 = math.erf(distance[i] * z2)
            fit_potential.append(K * (q_c * erf1 / distance[i]+ q_s * erf2 / distance[i]))
    return fit_potential


def gaussian_core_shell_charge(distance, q_c, z2):
    q_c = charge
    fit_potential = []
    for i in range(len(distance)):
            erf2 = math.erf(distance[i] * z2)
            fit_potential.append(K * (q_c* erf2 / distance[i]))
    return fit_potential


def point_core_gaussian_gaussian_shell_charge(distance, q_c, q_s2, z1, z2):
    q_s1 = charge - q_c - q_s2
    fit_potential = []
    for i in range(len(distance)):
            erf1 = math.erf(distance[i] * z1)
            erf2 = math.erf(distance[i] * z2)
            erf = math.erf(distance[i] * z1 * z2 / (z1**2 + z2**2)**0.5)
            fit_potential.append(K * (q_c / distance[i] + q_s1 * erf1 / distance[i]+ q_s2 * erf2 / distance[i]))
    return fit_potential



class ChargeModel(Enum):
    POINT_CORE_1_SLATER_SHELL = 0
    POINT_CORE_2_SLATER_SHELL = 1
    POINT_CORE_1_SLATER_2_SLATER_SHELL = 2
    SLATER_CORE_2_SLATER_SHELL = 3
    POINT_CORE_SHELL = 4
    POINT_CORE_GAUSSIAN_GAUSSIAN_SHELL = 5
    GAUSSIAN_CORE_SHELL = 6
    POINT_CORE_GAUSSIAN_SHELL = 7
    SLATER1_CORE_SLATER1_SHELL = 8
    GAUSSIAN_CORE_GAUSSIAN_SHELL = 9


charge_models = {
    ChargeModel.POINT_CORE_1_SLATER_SHELL: "point core 1 slater shell charge",
    ChargeModel.POINT_CORE_2_SLATER_SHELL: "point core 2 slater shell charge",
    ChargeModel.POINT_CORE_1_SLATER_2_SLATER_SHELL: "point core 1 slater 2 slater shell charge",
    ChargeModel.SLATER_CORE_2_SLATER_SHELL: "slater core 2 slater shell charge",
    ChargeModel.POINT_CORE_SHELL: "point core shell charge",
    ChargeModel.POINT_CORE_GAUSSIAN_GAUSSIAN_SHELL: "point core gaussian gaussian shell charge",
    ChargeModel.GAUSSIAN_CORE_SHELL: "gaussian core shell charge",
    ChargeModel.POINT_CORE_GAUSSIAN_SHELL: "point core gaussian shell charge",
    ChargeModel.SLATER1_CORE_SLATER1_SHELL: "slater1 core slater1 shell charge",
    ChargeModel.GAUSSIAN_CORE_GAUSSIAN_SHELL: "gaussian core gaussian shell charge"
}


parameter_names = {
    ChargeModel.POINT_CORE_1_SLATER_SHELL: ["q_c", "z2"],
    ChargeModel.POINT_CORE_2_SLATER_SHELL: ["q_c", "z2"],
    ChargeModel.POINT_CORE_1_SLATER_2_SLATER_SHELL: ["q_c",  "q_s2", "z1", "z2"],
    ChargeModel.SLATER_CORE_2_SLATER_SHELL: ["q_c",  "z1", "z2"],
    ChargeModel.POINT_CORE_SHELL: ["q_c"],
    ChargeModel.POINT_CORE_GAUSSIAN_GAUSSIAN_SHELL: ["q_c", "q_s2", "z1", "z2"],
    ChargeModel.GAUSSIAN_CORE_SHELL: ["q_c", "z2"],
    ChargeModel.POINT_CORE_GAUSSIAN_SHELL: ["q_c",  "z2"],
    ChargeModel.SLATER1_CORE_SLATER1_SHELL: ["q_c",  "z1", "z2"],
    ChargeModel.GAUSSIAN_CORE_GAUSSIAN_SHELL: ["q_c",  "z1", "z2"]
}

# charge models
functions = [
    Point_core_1slater_shell_charge,
    Point_core_2slater_shell_charge,
    Point_core_1slater_2slater_shell_charge,
    slater_core_2slater_shell_charge,
    point_core_shell_charge,
    point_core_gaussian_gaussian_shell_charge,
    gaussian_core_shell_charge,
    point_core_gaussian_shell_charge,
    slater1_core_slater1_shell_charge,
    gaussian_core_gaussian_shell_charge
]



with open('output.json', 'r') as json_f:
    output_data = json.load(json_f)


#inputs from json file
K=1389
data = output_data['data']
charge=output_data['charge']
initial_guesses = output_data['initial_guesses']
bounds = output_data['bounds']


fig, (axes1, axes2, axes3, axes4, axes5, axes6, axes7) = plt.subplots(len(data), 1, figsize=(6, 14))


all_params = {}

for i, (compound, (distance_data, potential_data)) in enumerate(data.items()):


    params = {}
    popts_compound = []
    pcovs_compound = []


    if compound=='F' or compound=='Cl' or compound=='Br' or compound=='I':
        charge=-1

    if compound=='Li' or compound=='Na' or compound=='K':
        charge=1

    for func_index, func in enumerate(functions):



        popt, pcov = curve_fit(func, distance_data, potential_data,
                                p0=initial_guesses[compound][func_index],
                                bounds=bounds[compound][func_index],
                                maxfev=2000)




        charge_model = charge_models[ChargeModel(func_index)]
        param_names = parameter_names[ChargeModel(func_index)]
        # params[charge_model] = {
        #     'popt': popt.tolist(),
        #     'pcov': pcov.tolist()
        # }

        popts_compound.append(popt)
        pcovs_compound.append(pcov)



        charge_model_compound = func(distance_data, *popt)


        if func_index in [0, 1, 7]:
            q_c_opt, z2_opt = popt
            q_s_opt = charge - q_c_opt
            params[f"q_c_{compound}_{func_index}"], params[f"q_s_{compound}_{func_index}"], \
            params[f"z2_{compound}_{func_index}"]= q_c_opt.tolist(), q_s_opt.tolist(), z2_opt.tolist()
            print(f"{compound} & {charge_model} & q_c_{compound}_{func_index}: {q_c_opt:.3f} & q_s_{compound}_{func_index}: {q_s_opt:.3f} & zeta_{compound}_{func_index}: {z2_opt*10:.2f} \\\\")

        elif func_index == 2 or func_index == 5:
            q_c_opt, q_s2_opt, z1_opt, z2_opt = popt
            q_s1_opt = charge - q_c_opt - q_s2_opt
            params[f"q_c_{compound}_{func_index}"], params[f"q_s1_{compound}_{func_index}"], params[f"q_s2_{compound}_{func_index}"], \
            params[f"z1_{compound}_{func_index}"], params[f"z2_{compound}_{func_index}"] = q_c_opt.tolist(), q_s1_opt.tolist(), \
            q_s2_opt.tolist(), z2_opt.tolist(), z1_opt.tolist()
            print(f"{compound} & {charge_model} & q_c_{compound}_{func_index}: {q_c_opt:.3f} & q_s1_{compound}_{func_index}: {q_s1_opt:.3f} & q_s2_{compound}_{func_index}: {q_s2_opt:.3f} & zeta_core_{compound}_{func_index}: {z1_opt*10:.2f} & zeta_shell_{compound}_{func_index}: {z2_opt*10:.2f} \\\\")

        elif func_index in [3, 8, 9]:
            q_c_opt,  z1_opt, z2_opt = popt
            q_s_opt = charge - q_c_opt
            params[f"q_c_{compound}_{func_index}"], params[f"q_s_{compound}_{func_index}"],  \
            params[f"z1_{compound}_{func_index}"], params[f"z2_{compound}_{func_index}"] = q_c_opt.tolist(), q_s_opt.tolist(), \
            z1_opt.tolist(), z2_opt.tolist()
            print(f"{compound} & {charge_model} & q_c_{compound}_{func_index}: {q_c_opt:.3f} & q_s_{compound}_{func_index}: {q_s_opt:.3f} & zeta_core_{compound}_{func_index}: {z1_opt*10:.2f} & zeta_shell_{compound}_{func_index}: {z2_opt*10:.2f} \\\\")

        elif func_index== 4:
            q_c_opt = popt
            q_c_opt=charge
            params[f"q_c_{compound}_{func_index}"] = q_c_opt
            print(f"{compound} & {charge_model} &  q_c_{compound}_{func_index}: {q_c_opt:.3f}  \\\\")

        elif func_index== 6:
            q_c_opt, z2_opt = popt
            params[f"q_c_{compound}_{func_index}"], params[f"z2_{compound}_{func_index}"] = q_c_opt.tolist(), z2_opt.tolist()
            print(f"{compound} & {charge_model} & q_c_{compound}_{func_index}: {q_c_opt:.3f} & zeta_core_{compound}_{func_index}: {z2_opt*10:.3f}  \\\\")



        rmse = np.sqrt(np.mean((np.array(charge_model_compound) - np.array(potential_data))**2))
        print(f"RMSE for {compound} & {charge_model}: {rmse:.2f} \\\\")




        if i == 0:
            label = 'A) Fluoride'
        elif i == 1:
            label = 'B) Chloride'
        elif i == 2:
            label = 'C) Bromide'
        elif i == 3:
            label = 'D) Iodide'

        elif i == 4:
            label = 'E) Lithium'

        elif i == 5:
            label = 'E) Sodium'

        elif i == 6:
            label = 'E) Potassium'

        axes = [axes1, axes2, axes3, axes4, axes5, axes6, axes7]
        axes[i].plot(np.array(distance_data), np.array(charge_model_compound), label=f'{charge_model}')
        axes[i].text(1.05, 1.05, label, transform=axes[i].transAxes, fontweight='bold', va='top', fontsize='xx-small')

    all_params[compound] = params
    print()




json_f_path = "params.json"
with open(json_f_path, "w") as json_f:
     json.dump(all_params, json_f, indent=4)
     json_f.write('\n')


axes1.legend(bbox_to_anchor=(1.05, 1), loc='upper left',fontsize='xx-small')
plt.tight_layout()


plt.savefig('Fit-Alkali_Halides.pdf')
plt.show()
