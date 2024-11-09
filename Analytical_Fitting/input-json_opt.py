#!/usr/bin/env python
import json
import numpy as np





def process_data(filename_prefix,filename_prefix_esp, skiprows):
    grid_data = np.loadtxt(f'ESP_Alkali_Halides/{filename_prefix}.dat', skiprows=skiprows)
    grid_esp_data = np.loadtxt(f'ESP_Alkali_Halides/grid_esp-{filename_prefix_esp}.dat', skiprows=skiprows)
    distances = np.sqrt(np.sum(grid_data ** 2, axis=1))
    esp_potentials = grid_esp_data * 2625.4996
    sorted_indices = np.argsort(distances)
    distance_data = distances[sorted_indices]
    potential_data = esp_potentials[sorted_indices]
    return distance_data, potential_data


distance_data_br, potential_data_br = process_data('grid-br2','hf-br', 0)
distance_data_f, potential_data_f = process_data('grid-f2','hf-f', 0)
distance_data_i, potential_data_i = process_data('grid-i2','hf-i', 0)
distance_data_cl, potential_data_cl = process_data('grid-cl2','hf-cl', 0)
distance_data_li, potential_data_li = process_data('grid-br2','hf-li', 0)
distance_data_na, potential_data_na = process_data('grid-f2','hf-na', 0)
distance_data_k, potential_data_k = process_data('grid-i2','hf-k', 0)



distance_data_br_list = distance_data_br.tolist()
potential_data_br_list = potential_data_br.tolist()
distance_data_f_list = distance_data_f.tolist()
potential_data_f_list = potential_data_f.tolist()
distance_data_i_list = distance_data_i.tolist()
potential_data_i_list = potential_data_i.tolist()
distance_data_cl_list = distance_data_cl.tolist()
potential_data_cl_list = potential_data_cl.tolist()
distance_data_li_list = distance_data_li.tolist()
potential_data_li_list = potential_data_li.tolist()
distance_data_na_list = distance_data_na.tolist()
potential_data_na_list = potential_data_na.tolist()
distance_data_k_list = distance_data_k.tolist()
potential_data_k_list = potential_data_k.tolist()



charge = -1


data = {
    'Br': [distance_data_br_list, potential_data_br_list],
    'F': [distance_data_f_list, potential_data_f_list],
    'Cl': [distance_data_cl_list, potential_data_cl_list],
    'I': [distance_data_i_list, potential_data_i_list],
    'Li': [distance_data_li_list, potential_data_li_list],
    'Na': [distance_data_na_list, potential_data_na_list],
    'K': [distance_data_k_list, potential_data_k_list]
}




initial_guesses = {
    'Br': [
        [.1, .001],          # Point_core_1slater_shell,               q_core, q_shell, zeta_shell
        [.1, .1],              # Point_core_2slater_shell,               q_core, q_shell, zeta_shell
        [.1, -1.1, 1, .1],      # Point_core_1slater_2slater_shell,       q_core, q_shell_1,q_shell_2, zeta_shell_1, zeta_shell_2
        [.1, .10, 1],         # slater_core_2slater_shell,              q_core, q_shell, zeta_core, zeta_shell
        [-1],                    # point_core_shell,                       q_core
        [.1, -1.10, .1, .1],  # point_core_gaussian_gaussian,           q_core, q_shell_1, q_shell_2, zeta_shell_1, zeta_shell_2
        [-1, .1],                 # gaussian_core_shell,                    q_core, zeta_core
        [.1, .1],             # point_core_gaussian_shell,              q_core, q_shell, zeta_shell
        [.1, 30, 1],          # slater1_core_slater1_shell,             q_core, q_shell, zeta_core, zeta_shell
        [.1, 30, 1]           # gaussian_core_gaussian_shell,           q_core, q_shell, zeta_core, zeta_shell
    ],
    'F': [
        [.1, .1],
        [.1, .1],
        [.1, -.10, 1, 1],
        [.1, 1.0, 1],
        [-1],
        [.1, -1.10, .1, .1],
        [-1, .1],
        [.1, .1],
        [.1, 30, 1],
        [.1, 30, 1]
    ],
    'Cl': [
        [.1, .3],
        [.1, .9],
        [.1, -.10, 1, 1],
        [.1, 3, 1],
        [-1],
        [.1, -1.10, 1, .1],
        [-1, .1],
        [.1,  .1],
        [.1, 30, 1],
        [.1,  30, 1]
    ],
    'I': [
        [.1,  .1],
        [.1,  .1],
        [.1,  -.10, 1, 1],
        [.1,  1.0, 1],
        [-1],
        [.1,  -1.10, .1, .1],
        [-1, .1],
        [.1,  .1],
        [.1,  30, 1],
        [.1,  30, 1]
    ],
       'Li': [
        [.1,  .1],
        [.1,  .1],
        [.1,  -.10, .1, .1],
        [.1,  1.0, .1],
        [1],
        [.1,  -.10, .1, .1],
        [.99999, .1],
        [.1, .1],
        [.1, 1, .1],
        [.1, 1, .1]
    ],
    'Na': [
        [.1, .1],
        [.1, .1],#1
        [.1, -.10, 1, 1],
        [.1, 1.0, 1],
        [1],
        [.1, -.10, 1, .1],
        [.99999, .1],
        [.1, .1],
        [.1, 10, 1],
        [.1, 30, 1]
    ],
    'K': [
        [.1, .1],
        [.1, .1],
        [.1, -.10, 1, 1],
        [.1, 1.0, 1],
        [1],
        [.1, -.10, .1, .1],
        [.99999, .1],
        [.1, .1],
        [.1, 30, 1],
        [.1, 30, 1]
    ]
}

bounds = {
    'Br': [
        ([.1, .001], [26, 10]),
        ([.1, .1], [26, 40]),
        ([.1, -10, 1.0, .1], [26, -0.01, 2.2, 40]),
        ([.1, .10, 1], [26, 40, 40]),
        ([-1], [-.99999999]),
        ([.1, -10, .1, .1], [26, -0.01, 40, 40]),
        ([-1, .1], [-.999, 40]),
        ([.1, .1], [16,  30]),
        ([.1,  1, 1], [46,  40, 19]),
        ([.1,  1, 1], [46,  40, 19])
    ],
    'F': [
        ([.1,  .1], [26,  40]),
        ([.1,  .1], [26,  40]),
        ([.01, -10, 1.0, 1], [26, -0.01,  2.2, 40]),
        ([.1, 1.0, 1], [26,  40, 40]),
        ([-1], [-.999999999]),
        ([.01,  -10, .1, .1], [5,  -0.01, 10, 10]),
        ([-1, .1], [-.999, 10]),
        ([.1,  .1], [16,  30]),
        ([.1,  1, 1], [46,  40,19]),
        ([.1,  1, 1], [46,  40,19])
    ],
    'Cl': [
        ([.1,  .1], [30, 4]),
        ([.1,  .1], [26,  4]),
        ([.1,  -10, 1.0, 1], [26,  -0.01, 2.2, 4]),
        ([.1,  1.0, 1], [26,  4, 2]),
        ([-1], [-.999999999]),
        ([.1, -10, 1.0, .1], [26, -0.01, 4, 2.2]),
        ([-1, .1], [-.999, 40]),
        ([.1,  .1], [16,  30]),
        ([.1,  1, 1], [46,  40, 20]),
        ([.1,  1, 1], [46,  40, 19])
    ],

    'I': [
        ([.1, .1], [26, 40]),
        ([.1,  .1], [26,  40]),
        ([.1,  -10, 1.0, 1], [26,  -0.01, 2.2, 40]),
        ([.1,  1.0, 1], [26,  40, 40]),
        ([-1], [-.999999999999]),
        ([.1,  -10, .1, .1], [26,  -0.01, 40, 40]),
        ([-1, .1], [-.9999999, 40]),
        ([.1,  .1], [16,  30]),
        ([.1,  1, 1], [46,  40, 19]),
        ([.1,  1, 1], [46,  40, 19])
    ],
    'Li': [
        ([.1,  .1], [26,  3]),
        ([.1,  .1], [26,  3]),
        ([.1,  -10, .1, .1], [26, -0.01,  2.2, 40]),
        ([.1,  1.0, .1], [26,  40, 1]),
        ([1], [1.000001]),
        ([.01,  -10, .1, .1], [56, -0.01, 100, 100]),
        ([.99999, .1], [1, 50]),#0.01 10
        ([.1,  .1], [16,  1]),
        ([.1,  1, .1], [46,  20,1]),
        ([.1,  1, .1], [46,  20,1])
    ],
    'Na': [
        ([.1, .1], [26,  40]),
        ([.1,  .1], [26,  40]),
        ([.1,  -10, 1.0, 1], [26,  -0.01, 2.2, 40]),
        ([.1,  1, 1], [26,  40, 20]),
        ([1], [1.0000001]),
        ([.1,  -10, 1, .1], [26,  -0.01, 40, 2.2]),
        ([.99999, .1], [1, 40]),
        ([.1,  .1], [26,  10]),
        ([.1,  1, 1], [46,  20, 20]),
        ([.1,  1, 1], [46,  40, 19])
    ],
    'K': [
        ([.1,  .1], [26,  40]),
        ([.1,  .1], [26,  40]),
        ([.1,  -10, 1.0, 1], [26,  -0.01, 2.2, 40]),
        ([.1,  1.0, 1], [26,  40, 40]),
        ([1], [1.000001]),
        ([.1,  -10, .1, .1], [26,  -0.01, 40, 40]),
        ([.99999, .1], [1, 40]),
        ([.1,  .1], [16,  30]),
        ([.1,  1, 1], [46,  40, 19]),
        ([.1,  1, 1], [46,  40, 19])
    ]
}



output_data = {
    "charge": charge,
    "data": data,
    "initial_guesses": initial_guesses,
    "bounds": bounds
}


with open('output.json', 'w') as json_file:
    json.dump(output_data, json_file, indent=4)
