#!/usr/bin/env python
import json
import numpy as np


T=100 #delta z= 0.01
S=000

def process_data(filename_prefix,filename_prefix_esp, skiprows):
    grid_data = np.loadtxt(f'ESP_Alkali_Halides/{filename_prefix}000.dat', skiprows=skiprows)
    grid_esp_data = np.loadtxt(f'ESP_Alkali_Halides/grid_esp-{filename_prefix_esp}000.dat', skiprows=skiprows)
    distances = np.sqrt(np.sum(grid_data ** 2, axis=1))
    esp_potentials = grid_esp_data * 2625.4996
    return distances, esp_potentials

N=200
data = {}
for elem in [ "Br", "F", "Cl", "Li", "Na", "K" ]:
    dist, pot = process_data('grid','hf-'+elem.lower(), N)
    data[elem] = [ dist.tolist(), pot.tolist() ]

bounds = {

    'F': [
        ([.1,  .1], [26,  30]),

        ([1,  .1], [26,  4]),

        ([1,  -10, .1, .1], [10,  -0.01, 30, 30]),

        ([1, -10, 1.0, 1], [20, -0.01,  30, 40]),


    ],

    'Cl': [
        ([1,  .1], [16,  30]),

        ([1,  .1], [30, 4]),

        ([4, -10, .10, .1], [30, -0.01, 30, 10]),

        ([4,  -10, 1., 1], [20,  -0.01, 30, 40])


    ],
    'Br': [
        ([1, .1], [50,  10]), #2

        ([1, .001], [100, 10]),

         ([5, -10, 1, .1], [40, -0.01, 40, 40]), #100

        ([5, -20, .1, .1], [40, -0.01, 30, 40]) #([.1, -30, 1, .1], [50, -9.01, 30, 40]) #1.0


    ],


    'Li': [
        ([1,  .1], [26,  5]),

        ([1,  .1], [26,  8]),

         ([2,  -10, .1, .1], [56, -0.01, 100, 100]),

        ([2,  -10, .1, .1], [16, -0.01,  40, 40])


    ],
    'Na': [
        ([1,  .1], [26,  10]),

        ([1, .1], [26,  40]),

        ([5,  -10, .1, .1], [26,  -0.01, 60, 60]),

        ([5,  -10, .1, .1], [26,  -0.01, 40, 40])

    ],
    'K': [
        ([1,  .1], [26,  10]),

        ([1,  .1], [26,  40]),

        ([5,  -10, .1, .1], [19,  -0.01, 40, 40]),

        ([5,  -10, .10, 1], [19,  -0.01, 40, 40]),

    ]
}

initial_guesses = {

    'F': [

        [.1, .1],

        [1, .1],

        [1, -1.10, .1, .1],

        [1, -.10, 1, 1]


    ],
    'Cl': [
         [1,  .1],

         [1, .1],

         [17, -1.10, .1, .1],

         [17, -.10, 1, 1]

    ],
        'Br': [

        [1, .1],              # point_core_gaussian_shell,              q_core, q_shell, zeta_shell

        [1, .001],          # Point_core_1slater_shell,               q_core, q_shell, zeta_shell

        [35, -1.10, 1, .1],  # point_core_gaussian_gaussian,           q_core, q_shell_1, q_shell_2, zeta_shell_1, zeta_shell_2

        [35, -10.1, 1, .1]      # Point_core_1slater_2slater_shell,       q_core, q_shell_1,q_shell_2, zeta_shell_1, zeta_shell_2

    ],

       'Li': [

        [1, .1],

        [1,  .1],

        [2,  -.10, .1, .1],

        [2,  -.10, .1, .1]

    ],
    'Na': [
        [1, .1],

        [1, .1],

        [10,  -.10, 1, .1],

        [10, -.10, .1, .1]

    ],
    'K': [
        [1, .1],

        [1, .1],

        [10,  -.10, .1, .1],

        [10, -.10, 1, 1]


    ]
}

output_data = {
    "charge": -1,
    "data": data,
    "initial_guesses": initial_guesses,
    "bounds": bounds
}


with open(f'output_4_{T}.json', 'w') as json_file:
    json.dump(output_data, json_file, indent=4)
