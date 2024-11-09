#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import json
import math
import os





#2S-2S- shell shell

def calculate_energy_shellshell2(r, xi, xj):
    rxi = r * xi
    rxj = r * xj

    S = (1 / r) * (
             (6 * np.exp(2 * (rxi + rxj)) * (np.power(np.power(rxi, 2) - np.power(rxj, 2), 7))) -
             (np.exp(2 * rxi) * np.power(rxi, 6) *
             (21 * np.power(rxi, 4) * np.power(rxj, 4) * (6 + 11 * rxj + 2 * np.power(rxj, 2)) -
              2 * np.power(rxj, 8) * (90 + 54 * rxj + 12 * np.power(rxj, 2) + np.power(rxj, 3)) +
              np.power(rxi, 8) * (6 + 9 * rxj + 6 * np.power(rxj, 2) + 2 * np.power(rxj, 3)) +
              np.power(rxi, 2) * np.power(rxj, 6) * (-390 - 69 * rxj + 18 * np.power(rxj, 2) + 4 * np.power(rxj, 3)) -
              np.power(rxi, 6) * np.power(rxj, 2) * (42 + 63 * rxj + 42 * np.power(rxj, 2) + 4 * np.power(rxj, 3)))) +
              (np.exp(2 * rxj) * np.power(rxj, 6) *
              (-24 * np.power(rxi, 10) - 2 * np.power(rxi, 11) - 69 * np.power(rxi, 7) * np.power(rxj, 2) +
              6 * np.power(rxj, 8) + 9 * rxi * np.power(rxj, 8) +
              4 * np.power(rxi, 9) * (-27 + np.power(rxj, 2)) +
              18 * np.power(rxi, 8) * (-10 + np.power(rxj, 2)) +
              6 * np.power(rxi, 2) * np.power(rxj, 6) * (-7 + np.power(rxj, 2)) -
              42 * np.power(rxi, 4) * np.power(rxj, 4) * (-3 + np.power(rxj, 2)) +
              np.power(rxi, 3) * np.power(rxj, 6) * (-63 + 2 * np.power(rxj, 2)) +
              6 * np.power(rxi, 6) * np.power(rxj, 2) * (-65 + 7 * np.power(rxj, 2)) +
              np.power(rxi, 5) * (231 * np.power(rxj, 4) - 4 * np.power(rxj, 6))))) / \
              (6 * np.exp(2 * (rxi + rxj)) * np.power(rxi - rxj, 7) * np.power(rxi + rxj, 7))

    return S


#core-2S
def slater2_core_shell_potential(r, xi):
    S = 1 / r - (6 + 9 * r * xi + 6 * pow(r, 2) * pow(xi, 2) + 2 * pow(r, 3) * pow(xi, 3)) / \
        (6 * math.exp(2 * r * xi) * r)

    return S



#core-1S
def slater_core_shell_potential(distance, alpha):
    return 1/distance - (1 + distance * alpha) / (math.exp(2 * distance * alpha) * distance)



#1S-1S- shell shell

def calculate_energy_shellshell(r, rxi, rxj):
    rxi2 = rxi * rxi
    rxj2 = rxj * rxj
    rxij = rxi + rxj
    exp2rxij = np.exp(2 * rxij)

    n = (exp2rxij * pow((rxi2 - rxj2), 3) +
                 np.exp(2 * rxj) * pow(rxj, 4) *
                 (-3 * rxi2 - pow(rxi, 3) + rxj2 + rxi * rxj2) -
                 np.exp(2 * rxi) * pow(rxi, 4) *
                 (rxi2 * (1 + rxj) - rxj2 * (3 + rxj)))

    d = exp2rxij * pow((rxi - rxj), 3) * pow(rxij, 3)
    energy = (1 / r) * (n / d)

    return energy


#1S-2S shell-shell




def double_Slater_1S_2S(r, xi, xj):
    rxi = r * xi
    rxj = r * xj
    n= (6 * np.exp(2*(rxi + rxj)) * np.power(np.power(rxi, 2) - np.power(rxj, 2), 5) +
         6 * np.exp(2*rxj) * np.power(rxj, 6) *
         (-4*np.power(rxi, 4) - np.power(rxi, 5) - 5*np.power(rxi, 2)*np.power(rxj, 2) +
          np.power(rxj, 4) + rxi*np.power(rxj, 4)) -
         np.exp(2*rxi) * np.power(rxi, 4) *
         (np.power(rxi, 6)*(6 + 9*rxj + 6*np.power(rxj, 2) + 2*np.power(rxj, 3)) -
          3*np.power(rxi, 4)*np.power(rxj, 2) *
          (10 + 15*rxj + 10*np.power(rxj, 2) + 2*np.power(rxj, 3)) +
          3*np.power(rxi, 2)*np.power(rxj, 4) *
          (20 + 33*rxj + 14*np.power(rxj, 2) + 2*np.power(rxj, 3)) -
          np.power(rxj, 6)*(84 + 63*rxj + 18*np.power(rxj, 2) + 2*np.power(rxj, 3))))

    #print(r)



    # d_denominator = rxi - rxj
    # d_numerator = rxi + rxj



    d=  (6 * np.exp(2*(rxi + rxj)) * np.power(rxi - rxj, 5) * np.power(rxi + rxj, 5))

    S = (1 / r) *  (n /d )


    return S



#
# (1/r)*(      (6*exp(2*(rxi + rxj)) *np.pow(np.pow(rxi, 2) - np.pow(rxj, 2), 5) +
#
#                           6*exp(2*rxj)*pow(rxj, 6)*
#
#                           (  -4*pow(rxi, 4) - pow(rxi, 5) - 5*pow(rxi, 2)*pow(rxj, 2) +
#
#                            pow(rxj, 4) + rxi*pow(rxj, 4)   ) -
#
#                           exp(2*rxi)*pow(rxi, 4)*
#
#                           (      pow(rxi, 6)*(6 + 9*rxj + 6*pow(rxj, 2) + 2*pow(rxj, 3)) -
#
#                            3*pow(rxi, 4)*pow(rxj, 2)*
#
#                            (  10 + 15*rxj + 10*pow(rxj, 2) + 2*pow(rxj, 3)) +
#
#                            3*pow(rxi, 2)*pow(rxj, 4)*
#
#                            (20 + 33*rxj + 14*pow(rxj, 2) + 2*pow(rxj, 3)) -
#
#                            pow(rxj, 6)*(84 + 63*rxj + 18*pow(rxj, 2) + 2*pow(rxj, 3)   )   )  )/
#
#                          (6*exp(2*(rxi + rxj))*pow(rxi - rxj, 5)*pow(rxi + rxj, 5))
#
#                          );




def gaussian_pointshell(distance, q_c, q_s,z2):
    erf2 = math.erf(distance * z2)
    energy=(q_c *q_s * erf2 / distance)
    return energy



def gaussian_shellshell(distance, q_c, q_s, z1, z2):
    erf = math.erf(distance * z1 * z2 / (z1**2 + z2**2)**0.5)
    energy= (q_c * q_s*erf / distance)
    return energy


########################################################################



def core_1slater_shell(distances, q_c_na, q_s_na, q_c_cl, q_s_cl, z_c_na, z_c_cl, z_s_na, z_s_cl):
    energy = []
    for i in range(len(distances)):

        energy.append(
            K * (

                (q_s_na*q_s_cl*(calculate_energy_shellshell(distances[i], z_s_cl*distances[i], z_s_na*distances[i])) ) +
                (q_c_na*q_c_cl*(calculate_energy_shellshell(distances[i], z_c_cl*distances[i], z_c_na*distances[i])) ) +
                (q_c_na*q_s_cl*(calculate_energy_shellshell(distances[i], z_s_cl*distances[i], z_c_na*distances[i])) ) +
                (q_s_na*q_c_cl*(calculate_energy_shellshell(distances[i], z_c_cl*distances[i], z_s_na*distances[i])) )

            )

        )

    distance_index = np.abs(distances - 1.).argmin()
    energy_at_0_9 = energy[distance_index]
    print(f"Energy at distance 1. Å from point-salter: {energy_at_0_9} kJ/mol")
    return energy




#"point core 1 slater shell charge"

def Point_core_1slater_shell(distances, q_c_na, q_s_na, q_c_cl, q_s_cl, z_na, z_cl):
    energy = []
    for i in range(len(distances)):

        energy.append(
            K * (
                (q_c_na * q_c_cl / distances[i] ) +
                (q_s_na*q_s_cl*(calculate_energy_shellshell(distances[i], z_cl*distances[i], z_na*distances[i])) ) +
                (q_s_cl*q_c_na * slater_core_shell_potential(distances[i], z_cl)) +
                (q_s_na*q_c_cl * slater_core_shell_potential(distances[i], z_na))
            )

        )

    distance_index = np.abs(distances - 1.).argmin()
    energy_at_0_9 = energy[distance_index]
    print(f"Energy at distance 1. Å from point-salter: {energy_at_0_9} kJ/mol")
    return energy



def Point_core_1slater_2slater_shell(distances, q_c_na, q_s1_na, q_s2_na, q_c_cl, q_s1_cl, q_s2_cl, z_c_na, z_c_cl, z_s_na, z_s_cl):
    energy = []
    for i in range(len(distances)):

        energy.append(
            K * (
                 (q_c_na * q_c_cl / distances[i] ) +
                 (q_s1_cl*q_c_na * slater_core_shell_potential(distances[i], z_c_cl)) +
                 (q_s1_na*q_c_cl * slater_core_shell_potential(distances[i], z_c_na))+
                 (q_s2_cl*q_c_na * slater2_core_shell_potential(distances[i], z_s_cl)) +
                 (q_s2_na*q_c_cl * slater2_core_shell_potential(distances[i], z_s_na))+
                 (q_s1_na*q_s1_cl* (calculate_energy_shellshell(distances[i], z_s_cl*distances[i], z_s_na*distances[i])) ) +
                 (q_s2_na*q_c_cl* (calculate_energy_shellshell2(distances[i], z_c_cl, z_s_na)) )  +
                 (q_s1_na*q_s2_cl* (double_Slater_1S_2S(distances[i], z_c_na, z_s_cl)))+
                 (q_s2_na*q_s1_cl* (double_Slater_1S_2S(distances[i], z_s_na, z_c_cl)))

            )

        )   #(q_s1_na*q_s2_cl**(double_Slater_1S_2S(distances[i], z_c_na, z_s_cl)

    # distance_index = np.abs(distances - 1.).argmin()
    # energy_at_0_9 = energy[distance_index]
    # print(f"Energy at distance 1. Å from point-salter: {energy_at_0_9} kJ/mol")
    return energy



# "slater core 2 slater shell charge"

def core_2slater_shell(distances, q_c_na, q_s_na, q_c_cl, q_s_cl, z_c_na, z_c_cl, z_s_na, z_s_cl):
    energy = []
    for i in range(len(distances)):

        energy.append(
            K * (

                (q_s_na*q_s_cl*(calculate_energy_shellshell2(distances[i], z_s_cl, z_s_na)) ) +
                (q_c_na*q_c_cl*(calculate_energy_shellshell2(distances[i], z_c_cl, z_c_na)) ) +
                (q_c_na*q_s_cl*(calculate_energy_shellshell2(distances[i], z_s_cl, z_c_na)) ) +
                (q_s_na*q_c_cl*(calculate_energy_shellshell2(distances[i], z_c_cl, z_s_na)) )

            )

        )

    distance_index = np.abs(distances - 1.).argmin()
    energy_at_0_9 = energy[distance_index]
    print(f"Energy at distance 1. Å from point-salter: {energy_at_0_9} kJ/mol")
    return energy


# "point core 2 slater shell charge"

def Point_core_2slater_shell(distances, q_c_na, q_s_na, q_c_cl, q_s_cl, z_na, z_cl):
    energy = []
    for i in range(len(distances)):

        energy.append(
            K * (
                (q_c_na * q_c_cl / distances[i] ) +
                (q_s_na*q_s_cl*(calculate_energy_shellshell2(distances[i], z_cl*distances[i], z_na*distances[i])) ) +
                (q_s_cl*q_c_na * slater2_core_shell_potential(distances[i], z_cl)) +
                (q_s_na*q_c_cl * slater2_core_shell_potential(distances[i], z_na))
            )

        )

    distance_index = np.abs(distances - 1.).argmin()
    energy_at_0_9 = energy[distance_index]
    print(f"Energy at distance 1. Å from point-salter: {energy_at_0_9} kJ/mol")
    return energy


#####################################################################

#"point core gaussian gaussian shell charge"
def Point_core_2gaussian_shell(distances, q_c_na, q_s1_na, q_s2_na, q_c_cl, q_s1_cl, q_s2_cl, z1_na, z1_cl, z2_na, z2_cl):
    energy = []
    for i in range(len(distances)):

        energy.append(
             (K *
                ((q_c_na * q_c_cl / distances[i])  +
                (gaussian_shellshell(distances[i], q_s1_cl, q_s1_na, z1_cl, z1_na)) +
                (gaussian_shellshell(distances[i], q_s2_cl, q_s2_na, z2_cl, z2_na)) +
                (gaussian_shellshell(distances[i], q_s1_cl, q_s2_na, z1_cl, z2_na)) +
                (gaussian_shellshell(distances[i], q_s2_cl, q_s1_na, z2_cl, z1_na)) +
                (gaussian_pointshell(distances[i], q_c_cl, q_s1_na, z1_na)) +
                (gaussian_pointshell(distances[i], q_c_cl, q_s2_na, z2_na)) +
                (gaussian_pointshell(distances[i], q_c_na, q_s1_cl, z1_cl)) +
                (gaussian_pointshell(distances[i], q_c_na, q_s2_cl, z2_cl)) )
            )

        )

    distance_index = np.abs(distances - 1.).argmin()
    energy_at_0_9 = energy[distance_index]
    print(f"Energy at distance 1. Å from point-2 Gaussians: {energy_at_0_9} kJ/mol")
    return energy


#"gaussian core shell charge"

def gaussian_core_shell(distances, q_c_na, q_c_cl, z1_na, z1_cl):
    energy = []
    for i in range(len(distances)):
        erf = math.erf(distances[i] * z1_na * z1_cl / (z1_na**2 + z1_cl**2)**0.5)
        energy.append(K * q_c_na* q_c_cl*erf / distances[i])



    distance_index = np.abs(distances - 1.).argmin()
    energy_at_0_9 = energy[distance_index]
    print(f"Energy at distance 1. Å from point-2 Gaussians: {energy_at_0_9} kJ/mol")
    return energy


#"point core shell charge"

def Point_core_shell(distances, q_c_na, q_c_cl):
    energy = []
    for i in range(len(distances)):

        energy.append(K*q_c_na * q_c_cl / distances[i])


    distance_index = np.abs(distances - 1.).argmin()
    energy_at_0_9 = energy[distance_index]
    print(f"Energy at distance 1. Å from point-2 Gaussians: {energy_at_0_9} kJ/mol")
    return energy


#"point core gaussian shell charge"

def Point_core_gaussian_shell(distances, q_c_na, q_s1_na, q_c_cl, q_s1_cl, z1_na, z1_cl):
    energy = []
    for i in range(len(distances)):

        energy.append(
             K *
                (  (q_c_na * q_c_cl / distances[i]) +
                (gaussian_shellshell(distances[i], q_s1_cl, q_s1_na, z1_cl, z1_na)) +

                (gaussian_pointshell(distances[i], q_c_cl, q_s1_na, z1_na)) +

                (gaussian_pointshell(distances[i], q_c_na, q_s1_cl, z1_cl))
                )
            )



    distance_index = np.abs(distances - 1.).argmin()
    energy_at_0_9 = energy[distance_index]
    print(f"Energy at distance 1. Å from point-2 Gaussians: {energy_at_0_9} kJ/mol")
    return energy



#"gaussian core gaussian shell charge"
def gaussian_core_gaussian_shell(distances, q_s1_na, q_s2_na, q_s1_cl, q_s2_cl, z1_na, z1_cl, z2_na, z2_cl):
    energy = []
    for i in range(len(distances)):

        energy.append(
             (K *
                (
                (gaussian_shellshell(distances[i], q_s1_cl, q_s1_na, z1_cl, z1_na))  +
                (gaussian_shellshell(distances[i], q_s2_cl, q_s2_na, z2_cl, z2_na))  +
                (gaussian_shellshell(distances[i], q_s1_cl, q_s2_na, z1_cl, z2_na))  +
                (gaussian_shellshell(distances[i], q_s2_cl, q_s1_na, z2_cl, z1_na))


                )
            )

        )

    distance_index = np.abs(distances - 1.).argmin()
    energy_at_0_9 = energy[distance_index]
    print(f"Energy at distance 1. Å from point-2 Gaussians: {energy_at_0_9} kJ/mol")
    return energy





##############################################################################
directory = 'SAPT_Alkali_Halides'
files = [os.path.join(directory, 'distances_Electrostatics-licl30.txt'),
         os.path.join(directory, 'distances_Electrostatics-nacl30.txt')]

data2 = [np.loadtxt(file) for file in files]
distances = np.arange(1, 4.6, 0.1)
K = 1389

with open('params.json', 'r') as json_f:
    params = json.load(json_f)

cations = "Li", "Na"
anion = "Cl"

function_values = {}

func_index_to_name = {
    0: "P+1S",
    1: "P+2S",
    2: "P+1S+2S",
    3: "2S+2S",
    4: "P",
    5: "P+2Gs",
    6: "G",
    7: "P+G",
    8: "1S+1S",
    9: "G+G"
}

func_index_to_function = {
    0: Point_core_1slater_shell,
    1: Point_core_2slater_shell,
    2: Point_core_1slater_2slater_shell,
    3: core_2slater_shell,
    4: Point_core_shell,
    5: Point_core_2gaussian_shell,
    6: gaussian_core_shell,
    7: Point_core_gaussian_shell,
    8: core_1slater_shell,
    9: gaussian_core_gaussian_shell
}

# for file, dataset in zip(files, data2):
#     fig, axes = plt.subplots(10, 1, figsize=(6, 14))
#
#     for i, cation in enumerate(cations):
#         x, y = dataset[:, 0], dataset[:, 1]
#         sorted_indices = np.argsort(x)
#         x1 = x[sorted_indices]
#         y1 = y[sorted_indices]
#
#         for func_index in range(10):
#             if func_index==2 or func_index==5:
#                 function_name = f"{func_index_to_name[func_index]}_{cation}{anion}"
#                 function_values[function_name] = func_index_to_function[func_index](
#                     distances,
#                     *[params[atom][f"q_{key}_{atom}_{func_index}"] for key in ["c", "s1", "s2"] for atom in [cation, anion]],
#                     *[params[atom][f"z{key}_{atom}_{func_index}"] for key in ["1", "2"] for atom in [cation, anion]]
#                 )
#
#
#             if func_index in [0, 1, 7]:
#                 function_name = f"{func_index_to_name[func_index]}_{cation}{anion}"
#                 function_values[function_name] = func_index_to_function[func_index](
#                     distances,
#                     *[params[atom][f"q_{key}_{atom}_{func_index}"] for key in ["c", "s"] for atom in [cation, anion]],
#                     *[params[atom][f"z{key}_{atom}_{func_index}"] for key in [ "2"] for atom in [cation, anion]]
#                 )
#
#             if func_index in [3, 8, 9]:
#                function_name = f"{func_index_to_name[func_index]}_{cation}{anion}"
#                function_values[function_name] = func_index_to_function[func_index](
#                    distances,
#                    *[params[atom][f"q_{key}_{atom}_{func_index}"] for key in ["c", "s"] for atom in [cation, anion]],
#                    *[params[atom][f"z{key}_{atom}_{func_index}"] for key in ["1", "2"] for atom in [cation, anion]]
#                )
#
#
#
#             if func_index ==4:
#                function_name = f"{func_index_to_name[func_index]}_{cation}{anion}"
#                function_values[function_name] = func_index_to_function[func_index](
#                    distances,
#                    *[params[atom][f"q_{key}_{atom}_{func_index}"] for key in ["c"] for atom in [cation, anion]]
#                )
#
#             if func_index==6:
#                function_name = f"{func_index_to_name[func_index]}_{cation}{anion}"
#                function_values[function_name] = func_index_to_function[func_index](
#                    distances,
#                    *[params[atom][f"q_{key}_{atom}_{func_index}"] for key in ["c"] for atom in [cation, anion]],
#                    *[params[atom][f"z{key}_{atom}_{func_index}"] for key in ["2"] for atom in [cation, anion]]
#                )
#
#
#             axes[func_index].plot(x1, y1, label=f"SAPT_{cation}{anion}", color='r')
#             axes[func_index].plot(distances, function_values[function_name], color='black', label=function_name)
#
#             axes[func_index].set_xlabel('Distance ($\AA$)', fontsize=12)
#             axes[func_index].set_ylabel('Electrostatic energies (kJ/mol)', fontsize=8)
#             axes[func_index].tick_params(axis='x', labelsize=8)
#             axes[func_index].tick_params(axis='y', labelsize=8)
#             axes[func_index].legend(fontsize=12)
#
#     plt.tight_layout()
#     plt.savefig(f'SAPT_{anion}_{cation}.pdf')
#     plt.show()





directory = 'SAPT_Alkali_Halides'
files = [os.path.join(directory, 'distances_Electrostatics-licl30.txt'),
         os.path.join(directory, 'distances_Electrostatics-nacl30.txt')]


data2 = [np.loadtxt(file) for file in files]
distances = np.arange(1, 4.6, 0.1)
K=1389

with open('params.json', 'r') as json_f:
    params = json.load(json_f)


cations = ["Li", "Na"]
anion = "Cl"

function_values = {}

func_index_to_name = {
    0: "P+1S",
    1: "P+2S",
    2: "P+1S+2S",
    3: "2S+2S",
    4: "P",
    5: "P+2Gs",
    6: "G",
    7: "P+G",
    8: "1S+1S",
    9: "G+G"
}

func_index_to_function = {
    0: Point_core_1slater_shell,
    1: Point_core_2slater_shell,
    2: Point_core_1slater_2slater_shell,
    3: core_2slater_shell,
    4: Point_core_shell,
    5: Point_core_2gaussian_shell,
    6: gaussian_core_shell,
    7: Point_core_gaussian_shell,
    8: core_1slater_shell,
    9: gaussian_core_gaussian_shell
}


for file, dataset in zip(files, data2):
    fig, axes = plt.subplots(10, 1, figsize=(6, 18))

    for i in range(10):
            x, y = dataset[:, 0], dataset[:, 1]
            sorted_indices = np.argsort(x)
            x1 = x[sorted_indices]
            y1 = y[sorted_indices]
            for k, cation in enumerate(cations):
                if file == f'SAPT_Alkali_Halides/distances_Electrostatics-{cation.lower()}cl30.txt':
                    distance_index = np.abs(x1 - 2.).argmin()
                    energy_at_0_9 = y1[distance_index]
                    print(f"Energy at distance 1. Å from SAPT: {energy_at_0_9} kJ/mol")

                    func_index = i
                    function_name = f"{func_index_to_name[func_index]}_{cation}{anion}"
                    function_values[function_name] = None

                    if func_index == 0:
                        function_name = f"P+1S_{cation}{anion}"
                        function_values[f"P+1S_{cation}{anion}"] = Point_core_1slater_shell(distances, params[f"{cation}"][f"q_c_{cation}_{func_index}"], params[f"{cation}"][f"q_s_{cation}_{func_index}"], params[f"{anion}"][f"q_c_{anion}_{func_index}"], params[f"{anion}"][f"q_s_{anion}_{func_index}"], params[f"{cation}"][f"z2_{cation}_{func_index}"], params[f"{anion}"][f"z2_{anion}_{func_index}"])

                    elif func_index == 1:
                        function_name = f"P+2S_{cation}{anion}"
                        function_values[f"P+2S_{cation}{anion}"] = Point_core_2slater_shell(distances, params[f"{cation}"][f"q_c_{cation}_{func_index}"], params[f"{cation}"][f"q_s_{cation}_{func_index}"], params[f"{anion}"][f"q_c_{anion}_{func_index}"], params[f"{anion}"][f"q_s_{anion}_{func_index}"], params[f"{cation}"][f"z2_{cation}_{func_index}"], params[f"{anion}"][f"z2_{anion}_{func_index}"])

                          #def Point_core_1slater_2slater_shell(distances, q_c_na, q_s1_na, q_s2_na, q_c_cl, q_s1_cl, q_s2_cl, z_c_na, z_c_cl, z_s_na, z_s_cl):
                    elif func_index == 2:
                       function_name = f"P+1S+2S_{cation}{anion}"
                       function_values[f"P+1S+2S_{cation}{anion}"] = Point_core_1slater_2slater_shell(distances, params[f"{cation}"][f"q_c_{cation}_{func_index}"], params[f"{cation}"][f"q_s1_{cation}_{func_index}"], params[f"{cation}"][f"q_s2_{cation}_{func_index}"], params[f"{anion}"][f"q_c_{anion}_{func_index}"], params[f"{anion}"][f"q_s1_{anion}_{func_index}"], params[f"{anion}"][f"q_s2_{anion}_{func_index}"], params[f"{cation}"][f"z1_{cation}_{func_index}"],params[f"{anion}"][f"z1_{anion}_{func_index}"], params[f"{cation}"][f"z2_{cation}_{func_index}"], params[f"{anion}"][f"z2_{anion}_{func_index}"])


                    elif func_index == 3:
                        function_name = f"2S+2S_{cation}{anion}"
                        function_values[f"2S+2S_{cation}{anion}"] = core_2slater_shell(distances, params[f"{cation}"][f"q_c_{cation}_{func_index}"], params[f"{cation}"][f"q_s_{cation}_{func_index}"], params[f"{anion}"][f"q_c_{anion}_{func_index}"], params[f"{anion}"][f"q_s_{anion}_{func_index}"], params[f"{cation}"][f"z1_{cation}_{func_index}"], params[f"{anion}"][f"z1_{anion}_{func_index}"], params[f"{cation}"][f"z2_{cation}_{func_index}"], params[f"{anion}"][f"z2_{anion}_{func_index}"])

                    elif func_index == 4:
                        function_values[function_name] = Point_core_shell(distances, params[f"{cation}"][f"q_c_{cation}_{func_index}"], params[f"{anion}"][f"q_c_{anion}_{func_index}"])
                    elif func_index == 5:
                        function_values[function_name] = Point_core_2gaussian_shell(distances, params[f"{cation}"][f"q_c_{cation}_{func_index}"], params[f"{cation}"][f"q_s1_{cation}_{func_index}"], params[f"{cation}"][f"q_s2_{cation}_{func_index}"], params[f"{anion}"][f"q_c_{anion}_{func_index}"], params[f"{anion}"][f"q_s1_{anion}_{func_index}"], params[f"{anion}"][f"q_s2_{anion}_{func_index}"], params[f"{cation}"][f"z1_{cation}_{func_index}"], params[f"{anion}"][f"z1_{anion}_{func_index}"], params[f"{cation}"][f"z2_{cation}_{func_index}"], params[f"{anion}"][f"z2_{anion}_{func_index}"])
                    elif func_index == 6:
                        function_values[function_name] = gaussian_core_shell(distances, params[f"{cation}"][f"q_c_{cation}_{func_index}"], params[f"{anion}"][f"q_c_{anion}_{func_index}"], params[f"{cation}"][f"z2_{cation}_{func_index}"], params[f"{anion}"][f"z2_{anion}_{func_index}"])
                    elif func_index == 7:
                        function_values[function_name] = Point_core_gaussian_shell(distances, params[f"{cation}"][f"q_c_{cation}_{func_index}"], params[f"{cation}"][f"q_s_{cation}_{func_index}"], params[f"{anion}"][f"q_c_{anion}_{func_index}"], params[f"{anion}"][f"q_s_{anion}_{func_index}"], params[f"{cation}"][f"z2_{cation}_{func_index}"], params[f"{anion}"][f"z2_{anion}_{func_index}"])
                    elif func_index == 8:
                        function_values[function_name] = core_1slater_shell(distances, params[f"{cation}"][f"q_c_{cation}_{func_index}"], params[f"{cation}"][f"q_s_{cation}_{func_index}"], params[f"{anion}"][f"q_c_{anion}_{func_index}"], params[f"{anion}"][f"q_s_{anion}_{func_index}"], params[f"{cation}"][f"z1_{cation}_{func_index}"], params[f"{anion}"][f"z1_{anion}_{func_index}"], params[f"{cation}"][f"z2_{cation}_{func_index}"], params[f"{anion}"][f"z2_{anion}_{func_index}"])
                    elif func_index == 9:
                        function_values[function_name] = gaussian_core_gaussian_shell(distances, params[f"{cation}"][f"q_c_{cation}_{func_index}"], params[f"{cation}"][f"q_s_{cation}_{func_index}"], params[f"{anion}"][f"q_c_{anion}_{func_index}"], params[f"{anion}"][f"q_s_{anion}_{func_index}"], params[f"{cation}"][f"z1_{cation}_{func_index}"], params[f"{anion}"][f"z1_{anion}_{func_index}"], params[f"{cation}"][f"z2_{cation}_{func_index}"], params[f"{anion}"][f"z2_{anion}_{func_index}"])


                    axes[i].plot(x1, y1, label=f"SAPT_{cation}{anion}", color='r')
                    axes[i].plot(distances, function_values[function_name], color='black', label=function_name)

    for i in range(10):
            axes[i].set_xlabel('Distance ($\AA$)', fontsize=12)
            axes[i].set_ylabel('Electrostatic energies (kJ/mol)', fontsize=8)
            axes[i].tick_params(axis='x', labelsize=8)
            axes[i].tick_params(axis='y', labelsize=8)
            axes[i].legend(fontsize=12)

    plt.tight_layout()
    plt.savefig(f'SAPT_{anion}_{cation}.pdf')
    plt.show()
