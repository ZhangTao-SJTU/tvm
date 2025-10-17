from toolbox.patterns import Patterns
from toolbox.periodic import PeriodicTissue
from toolbox.pattern_trainer import train_random_cells
import os
import numpy as np
import pandas as pd
import glob
import sys


def copy_config(source_dir,destination_dir):
    source_file = sorted(glob.glob("{}*.cellParameters.input".format(source_dir)))[-1]
    print("copying {} to {}cellParameters.input".format(source_file, destination_dir))
    os.system("cp {} {}cellParameters.input".format(source_file, destination_dir))
    source_file = sorted(glob.glob("{}*.bulk.txt".format(source_dir)))[-1]
    print("copying {} to {}minimized.txt".format(source_file, destination_dir))
    os.system("cp {} {}minimized.txt".format(source_file, destination_dir))

def single_iteration (patternA, patternB):
    "Retraining patternA using minimized config of patternB"
    copy_config(source_dir = patternB._dir, destination_dir = patternA._dir)
    patternA._config.load_periodic_tissue_from_file("minimized.txt")
    patternA.load_cell_parameters()
    patternA.minimize_config()
    # patternA.run()
    patternA.run_to_max_iters(max_iters=100)
    "Retraining patternB using minimized config of patternA"
    copy_config(source_dir = patternA._dir, destination_dir = patternB._dir)
    patternA._config.load_periodic_tissue_from_file("minimized.txt")
    patternB.load_cell_parameters()
    patternB.minimize_config()
    # patternB.run()
    patternB.run_to_max_iters(max_iters=100)

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
    cpp_executable_dir = "/home/shabeeb/Projects/tvm-fire/build/"
    alpha = 0.025
    n_cells = 1
    tolerance = 1e-5
    average_cells_only = False
    if "cpp_executable_dir" in kwargs:
        cpp_executable_dir = kwargs["cpp_executable_dir"]
    if "alpha" in kwargs:
        alpha = kwargs["alpha"]
    if "n_cells" in kwargs:
        n_cells = kwargs["n_cells"]
    if "tolerance" in kwargs:
        tolerance = kwargs["tolerance"]
    if "average_cells_only" in kwargs:
        average_cells_only = kwargs["average_cells_only"]
    patterns_dict = {"patternA/":None, "patternB/":None}
    print("Parameters for initialization:")
    print("cpp_executable_dir:", cpp_executable_dir)
    print("alpha:", alpha)
    print("n_cells:", n_cells)
    print("tolerance:", tolerance)
    for pattern_name in patterns_dict:
        dir = run_dir + pattern_name
        os.system("cp -r {} {}".format(init_dir, dir))
        if pattern_name == "patternA/":
            patterns_dict[pattern_name] = train_random_cells(
                run_dir = dir,
                alpha = alpha,
                n_cells = n_cells,
                tolerance = tolerance, 
                cpp_executable_dir = cpp_executable_dir,
                average_cells_only = average_cells_only)
        elif pattern_name == "patternB/":
            exclude_cells = list(patterns_dict["patternA/"]._target_cell_to_stress.keys())
            print("Excluding cells in patternA from patternB:", exclude_cells)
            patterns_dict[pattern_name] = train_random_cells(
                run_dir = dir,
                alpha = alpha,
                n_cells = n_cells,
                tolerance = tolerance, 
                cpp_executable_dir = cpp_executable_dir,
                average_cells_only = average_cells_only)
        else:
            raise ValueError("Unknown pattern name: {}".format(pattern_name))
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
    n_iters = 100
    n_cells = 2
    alpha = 0.05
    tolerance = 1e-5
    if os.path.isdir("/Users/shabeebameen/Projects/tvm-fire/build/"):
        cpp_executable_dir = "/Users/shabeebameen/Projects/tvm-fire/build/"
    elif os.path.isdir("/home/shabeeb/Projects/tvm-fire/build/"):
        cpp_executable_dir = "/home/shabeeb/Projects/tvm-fire/build/"
    elif os.path.isdir("/home/mameen/tvm/build/"):
        cpp_executable_dir = "/home/mameen/tvm/build/"
    
    if not len(sys.argv) == 3:
        print("Usage: python3 multiple_patterns.py <input_dir> <run_dir>")
        sys.exit(1)
    input_dir = str(sys.argv[1])
    run_dir = str(sys.argv[2])
    # print("Input dir:", input_dir)
    print("Run dir:", run_dir)
    
    os.makedirs(run_dir , exist_ok=True)
    distances_array = []
    if os.path.isfile("{}distances.txt".format(run_dir)):
        distances_array = list(np.loadtxt("{}distances.txt".format(run_dir)))
    
    patterns = initialize_patterns(input_dir, run_dir,
        cpp_executable_dir = cpp_executable_dir,
        n_cells = n_cells,
        average_cells_only = True,
        alpha = alpha,
        tolerance = tolerance)
    # patterns = initialize_patterns_from_run_dir(run_dir)
    patternA = patterns["patternA/"]
    patternB = patterns["patternB/"]

    init_distance = [calculate_parameter_space_distance(patternA, patternB)]
    if not os.path.isfile("{}initial_distances.txt".format(run_dir)):
        np.savetxt("{}initial_distances.txt".format(run_dir), init_distance, fmt='%.12e')
    for iteration in range(n_iters):
        print("Starting pattern-switching iteration:", iteration)
        single_iteration(patternA, patternB)
        current_distance = calculate_parameter_space_distance(patternA, patternB)
        distances_array.append(current_distance)
        np.savetxt("{}distances.txt".format(run_dir), distances_array, fmt='%.12e')
        if current_distance < 1e-9:
            print("Converged at iteration", iteration)
            break
if __name__ == "__main__":
    main()