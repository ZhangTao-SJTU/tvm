from toolbox.cSection import makeSampleCrossSection
from toolbox.patterns import Patterns
from toolbox.periodic import PeriodicTissue
from toolbox import stress
from toolbox.pattern_trainer import train_random_cells
import matplotlib.pyplot as plt
from toolbox import stress
import os
import numpy as np
import pandas as pd
import glob
import sys


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
    patternA.minimize_config()
    patternA.run()
    copy_config(source_dir = patternA._dir, destination_dir = patternB._dir)
    patternA._config.load_periodic_tissue_from_file("minimized.txt")
    patternB.load_cell_parameters()
    patternB.minimize_config()
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

def initialize_patterns(init_dir,run_dir, **kwargs):
    #default parameters
    alpha = 0.025
    n_cells = 1
    tolerance = 1e-5
    if "alpha" in kwargs:
        alpha = kwargs["alpha"]
    if "n_cells" in kwargs:
        n_cells = kwargs["n_cells"]
    if "tolerance" in kwargs:
        tolerance = kwargs["tolerance"]
    patterns_dict = {"patternA/":None, "patternB/":None}
    for pattern_name in patterns_dict:

        dir = run_dir + pattern_name
        os.system("cp -r {} {}".format(init_dir, dir))
        patterns_dict[pattern_name] = train_random_cells(
            run_dir = dir, alpha = alpha, n_cells = n_cells, tolerance = tolerance)
    return patterns_dict

def initialize_patterns_from_run_dir(run_dir,**kwargs):
    #default parameters
    cpp_executable_dir = "/home/shabeeb/Projects/tvm-fire/build/"
    tolerance = 1e-5
    if "cpp_executable_dir" in kwargs:
        cpp_executable_dir = kwargs["cpp_executable_dir"]
    if "tolerance" in kwargs:
        tolerance = kwargs["tolerance"]
    patterns = {"patternA/":None, "patternB/":None}
    for pattern_name in patterns:
        dir = run_dir + pattern_name
        file = sorted(glob.glob("{}*.bulk.txt".format(dir)))[-1]
        file = file.split("/")[-1]
        print("Loading pattern from file:", file)
        tissue = PeriodicTissue.from_config(dir, file)
        patterns[pattern_name] = Patterns.periodic_tissue(tissue)
        patterns[pattern_name].set_cpp_executable_dir(cpp_executable_dir)
        patterns[pattern_name].set_tolerance(tolerance)
        df = pd.read_csv("{}initial_stress.csv".format(dir))
        cellID_to_stress = {int(i): float(stress) for i, stress in zip(df["cellID"].to_numpy(), df["Stress"].to_numpy())}
        patterns[pattern_name].set_target_cell_to_stress(cellID_to_stress)
        cost_values = list(np.loadtxt("{}costs.txt".format(dir)))
        patterns[pattern_name]._cost_values = cost_values
        patterns[pattern_name].set_iter_counter(len(cost_values))
        cellParameters_file = sorted(glob.glob("{}*.cellParameters.input".format(dir)))[-1]
        os.system("cp {} {}cellParameters.input".format(cellParameters_file, dir))
    return patterns

def main():
    n_iters = 1
    if not len(sys.argv) == 3:
        print("Usage: python3 multiple_patterns.py <input_dir> <run_dir>")
        sys.exit(1)
    input_dir = str(sys.argv[1])
    run_dir = str(sys.argv[2])
    print("Input dir:", input_dir)
    print("Run dir:", run_dir)
    
    if not os.path.isdir(run_dir):
        os.makedirs(run_dir)
    distances_array = []
    if os.path.isfile("{}distances.txt".format(run_dir)):
        distances_array = list(np.loadtxt("{}distances.txt".format(run_dir)))
    
    patterns = initialize_patterns(input_dir,run_dir)
    # patterns = initialize_patterns_from_run_dir(run_dir)
    patternA = patterns["patternA/"]
    patternB = patterns["patternB/"]

    distance = [calculate_parameter_space_distance(patternA, patternB)]
    if not os.path.isfile("{}initial_distances.txt".format(run_dir)):
        np.savetxt("{}initial_distances.txt".format(run_dir), distance, fmt='%.12e')
    # print("Initial distance:", distance)
    for iteration in range(n_iters):
        single_iteration(patternA, patternB)
        distance = calculate_parameter_space_distance(patternA, patternB)
        distances_array.append(distance)
        np.savetxt("{}distances.txt".format(run_dir), distances_array, fmt='%.12e')
        if distance < 1e-9:
            print("Converged at iteration", iteration)
            break
if __name__ == "__main__":
    main()