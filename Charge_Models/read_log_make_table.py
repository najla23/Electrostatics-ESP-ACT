#!/usr/bin/env python
import os


# make latex table using log file from Alexandria 

a="PC+GS-model"
b="elec"

def extract(log_file):
    table = []
    with open(log_file, 'r') as file:
        start_reading = False
        for line in file:
            if "Here are the best parameters I found, together with some summary statistics of the last population:" in line:
                start_reading = True
                next(file)
                continue
            elif "Final best genome." in line:
                break
            elif start_reading:
                if not line.strip():
                    continue
                if line.startswith('|'):
                    row = line.strip().split('|')
                    row = [cell.strip() for cell in row]
                    table.append(row[1:4])
    return table



def save(table, output_file):
    with open(output_file, 'w') as file:
        file.write("\\begin{table}[ht]\n")
        file.write(f"\\caption{{best parameters for {a} using {b}.}}\n")
        file.write("\\begin{tabular}{lcccc}\n")
        file.write("\\hline\n")
        file.write("Parameter & Atom & Best (Train) \\\\ \n")
        file.write("\\hline\n")
        for row in table:
            file.write(" & ".join(row) + " \\\\ \n")
        file.write("\\hline\n")
        file.write("\\end{tabular}\n")
        file.write("\\end{table}")

log_fp = os.path.join(os.getcwd(), 'tune_ff-elec-p.log')
output_fp = os.path.join(os.getcwd(), f'output_table_{a}_{b}.tex')

table = extract(log_fp)

save(table, output_fp)
