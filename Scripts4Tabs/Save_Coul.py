#!/usr/bin/env python
import os
import re
import json

log_file_path = "tune_ff-elec.log"

# Molecules + QM values
molecules = ['water#lithium', 'water#sodium', 'water#potassium', 'water#fluoride', 'water#chloride', 'water#bromide', 'water#water',
             'acetate#lithium', 'acetate#sodium', 'acetate#potassium', 'methylammonium#fluoride', 'methylammonium#chloride', 'methylammonium#bromide']
QM_values = [-123.68, -91.83, -67.7, -225.18, -71.91, -66.45, -31.59, -624.68, -659.87, -574.45, -556.74, -496.90, -472.40]

QM_values = [f'{value:.2f}' for value in QM_values]

coulomb_value_pattern = re.compile(r"^\s*(?:\d+\s+)?(?:-?\d+\.\d+\s+){3}(-?\d+\.\d+)\s+", re.MULTILINE)

results = {}

if os.path.exists('results.json'):
    with open('results.json', 'r') as json_file:
        results = json.load(json_file)

with open(log_file_path, 'r') as file:
    log_content = file.read()

for molecule, qm_value in zip(molecules, QM_values):
    molecule_match = re.search(fr"Name: {re.escape(molecule)}.*", log_content, re.DOTALL)
    if molecule_match:
        coulomb_section = molecule_match.group(0)
        qm_index = coulomb_section.split().index(str(qm_value))
        act_value = float(coulomb_section.split()[qm_index + 1])
    else:
        act_value = None

    results[molecule] = {'QM': qm_value, 'ACT': act_value}

with open('results.json', 'w') as json_file:
    json.dump(results, json_file, indent=4)

print("Data saved to results.json")
