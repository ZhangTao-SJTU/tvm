from toolbox.cSection import makeSampleCrossSection
from toolbox.patterns import Patterns
from toolbox.periodic import PeriodicTissue
from toolbox import stress
from pattern_trainer import train_random_cells
import matplotlib.pyplot as plt
from toolbox import stress
import os
import numpy as np
import pandas as pd
import glob

def copy_config(source_dir,destination_dir):
    # os.system("cp {}minimized.txt {}minimized.txt".format(source_dir, destination_dir))
    source_file = sorted(glob.glob("{}*.cellParameters.input".format(source_dir)))[-1]
    os.system("cp {} {}cellParameters.input".format(source_file, destination_dir))
    source_file = sorted(glob.glob("{}*.bulk.txt".format(source_dir)))[-1]
    os.system("cp {} {}minimized.txt".format(source_file, destination_dir))
def single_iteration (patternA, patternB):
    copy_config(source_dir = patternB._dir, destination_dir = patternA._dir)
    patternA._config.load_periodic_tissue_from_file("minimized.txt")
    patternA.load_cell_parameters()
    # patternA.run_to_max_iters(5)
    patternA.minimize_config()
    patternA.run()
    copy_config(source_dir = patternA._dir, destination_dir = patternB._dir)
    patternA._config.load_periodic_tissue_from_file("minimized.txt")
    patternB.load_cell_parameters()
    patternB.minimize_config()

    # patternB.run_to_max_iters(5)
    patternB.run()
# def initialize_pattern(dir, tolerance):
#     tissue = PeriodicTissue.from_config(dir,"minimized.txt")
#     pattern = Patterns.periodic_tissue(tissue)
#     pattern.set_tolerance(tolerance)
#     df = pd.read_csv("{}0.stresses.csv".format(pattern._dir))
#     cellID_to_final_stress = {cellID:stress for cellID,stress in zip(df["CellID"],df["Final"])}
#     pattern.set_target_cell_to_stress(cellID_to_final_stress)
#     pattern.set_tolerance(tolerance)
#     return pattern
def calculate_parameter_space_distance(patternA, patternB):
    cellParametersA = "{}cellParameters.input".format(patternA._dir)
    cellParametersB = "{}cellParameters.input".format(patternB._dir)
    df = pd.read_csv(cellParametersA, sep=" ",header=None)
    cellID_to_s0 = {int(i): [float(s0)] for i, s0 in zip(df[0].to_numpy(), df[2].to_numpy())}
    df = pd.read_csv(cellParametersB, sep=" ",header=None)
    for i, row in df.iterrows():
        cellID = row[0]
        s0 = row[2]
        if cellID in cellID_to_s0:
            cellID_to_s0[cellID].append(s0)
        else:
            raise ValueError("CellID {} not found in first pattern".format(cellID))
    distance = 0
    for cellID, s0s in cellID_to_s0.items():
        # print("CellID: {}, s0s: {}".format(cellID, s0s))
        distance += (s0s[0] - s0s[1])**2
    distance = distance**0.5
    # print(distance)
    return distance

def initialize_patterns(parent_dir):
    alpha = 0.025
    n_cells = 2
    tolerance = 1e-5
    patterns = {"patternA/":None, "patternB/":None}
    for dir in patterns:
        os.system("cp -r init/{} {}".format(parent_dir, dir))
        patterns[dir] = train_random_cells(parent_dir, dir, alpha, n_cells, tolerance)
    return patterns

def initialize_patterns_from_dirs(pattern_dirs):

    tolerance = 1e-5
    patterns = {}
    for dir in pattern_dirs:
        file = sorted(glob.glob("{}*.bulk.txt".format(dir)))[-1]
        file = file.split("/")[-1]
        print("Loading pattern from file:", file)
        tissue = PeriodicTissue.from_config(dir, file)
        patterns[dir] = Patterns.periodic_tissue(tissue)
        patterns[dir].set_tolerance(tolerance)
        df = pd.read_csv("{}initial_stress.csv".format(dir))
        cellID_to_stress = {int(i): float(stress) for i, stress in zip(df["cellID"].to_numpy(), df["Stress"].to_numpy())}
        patterns[dir].set_target_cell_to_stress(cellID_to_stress)
        patterns[dir]._cost_values = list(np.loadtxt("{}costs.txt".format(dir))) 
        patterns[dir].set_iter_counter(len(patterns[dir]._cost_values))
        cellParameters_file = sorted(glob.glob("{}*.cellParameters.input".format(dir)))[-1]
        os.system("cp {} {}cellParameters.input".format(cellParameters_file, dir))
    return patterns

def main():
    # dirA = "patternA/"
    # dirB = "patternB/"
    # for folder in [dirA, dirB]:
    #     if not os.path.isdir(folder):
    #         os.system("cp -r tests/{} {}".format(folder, folder))
    # tolerance = 1e-5
    # patternA = initialize_pattern(dirA, tolerance)
    # patternA.set_iter_counter(17)
    # patternB = initialize_pattern(dirB, tolerance)
    # patternB.set_iter_counter(3)
    out_dir = "patterns/"
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    n_iters = 20
    distances = []
    # patterns = initialize_patterns("7_1/")
    patterns = initialize_patterns_from_dirs(["patternA/", "patternB/"])
    patternA = patterns["patternA/"]
    patternB = patterns["patternB/"]

    distance = [calculate_parameter_space_distance(patternA, patternB)]
    np.savetxt("{}initial_distances.txt".format(out_dir), distance, fmt='%.5f')
    # print("Initial distance:", distance)
    for iteration in range(n_iters):
        single_iteration(patternA, patternB)
        distance = calculate_parameter_space_distance(patternA, patternB)
        distances.append(distance)
        np.savetxt("{}distances.txt".format(out_dir), distances, fmt='%.5f')
        for pattern in patterns:
            os.system("cp -r {} {}".format(pattern, out_dir))
        if distance < 1e-9:
            print("Converged at iteration", iteration)
            break

if __name__ == "__main__":
    main()