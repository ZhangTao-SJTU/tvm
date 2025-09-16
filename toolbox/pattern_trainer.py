from toolbox.cSection import makeSampleCrossSection
from toolbox.patterns import Patterns
from toolbox.periodic import PeriodicTissue
from toolbox import stress
import matplotlib.pyplot as plt
from toolbox import stress
import os
import numpy as np
import pandas as pd

def train_random_cells(run_dir, **kwargs):
    print("Training random cells in directory:", run_dir)
    #default parameters
    cpp_executable_dir = "/home/shabeeb/Projects/tvm-fire/build/"
    alpha = 0.025
    n_cells = 1
    sign = 1
    tolerance = 1e-6
    max_iters = 100
    average_cells_only = False
    exclude_cells = []
    target_stress = None
    if "cpp_executable_dir" in kwargs:
        cpp_executable_dir = kwargs["cpp_executable_dir"]
    if "n_cells" in kwargs:
        n_cells = kwargs["n_cells"]
    if "tolerance" in kwargs:
        tolerance = kwargs["tolerance"]
    if "max_iters" in kwargs:
        max_iters = kwargs["max_iters"]
    if "alpha" in kwargs:
        alpha = kwargs["alpha"]
    if "sign" in kwargs:
        sign = kwargs["sign"]
    if "average_cells_only" in kwargs:
        average_cells_only = kwargs["average_cells_only"]
    if "exclude_cells" in kwargs:
        exclude_cells = kwargs["exclude_cells"]
    if "target_stress" in kwargs:
        target_stress = kwargs["target_stress"]
    print("Parameters for training:")
    print("cpp_executable_dir:", cpp_executable_dir)
    print("n_cells:", n_cells)
    print("tolerance:", tolerance)
    print("max_iters:", max_iters)
    print("target stress:", target_stress)

    file = "minimized.txt"
    if os.path.isfile("{}minimized.txt".format(run_dir)):
        print("File exists:", "{}minimized.txt".format(run_dir))
    tissue = PeriodicTissue.from_config(run_dir,file)
    training_instance = Patterns.periodic_tissue(tissue)
    training_instance.set_cpp_executable_dir(cpp_executable_dir)
    training_instance.minimize_config()
    training_instance.set_tolerance(tolerance)
    training_instance.set_random_target_cells(
        n_cells = n_cells, 
        average_cells_only = average_cells_only,
        exclude_cells = exclude_cells)

    for cellID in training_instance._target_cell_to_stress:
        # sign = np.random.choice([-1, 1])
        initial_stress = stress.calculate_max_shear_stress(training_instance._config, cellID)
        if target_stress is None:
            training_instance._target_cell_to_stress[cellID] = np.round((1+sign*alpha) * initial_stress,3)
        else:
            training_instance._target_cell_to_stress[cellID] = target_stress
        # training_instance._target_cell_to_stress[cellID] = alpha
        print("Initial stress for cell {}: {}".format(cellID, initial_stress))
    print(training_instance._target_cell_to_stress)
    training_instance.initialize()
    training_instance.run_to_max_iters(max_iters)
    return training_instance

def resume_run(run_dir, **kwargs):
    print("Resuming runs in directory:", run_dir)
    #default parameters
    cpp_executable_dir = "/home/shabeeb/Projects/tvm-fire/build/"
    tolerance = 1e-7
    max_iters = 10
    if "cpp_executable_dir" in kwargs:
        cpp_executable_dir = kwargs["cpp_executable_dir"]
    if "tolerance" in kwargs:
        tolerance = kwargs["tolerance"]
    if "max_iters" in kwargs:
        max_iters = kwargs["max_iters"]
    print("Parameters for training:")
    print("cpp_executable_dir:", cpp_executable_dir)
    print("tolerance:", tolerance)
    print("max iters:", max_iters)
    
    costs = np.loadtxt("{}costs.txt".format(run_dir))
    iteration = len(costs)
    print(iteration)
    config_file = "{:04d}.bulk.txt".format(iteration-1)
    print("Loading configuration from: ", config_file)
    sample = PeriodicTissue.from_config(run_dir,config_file)
    training_instance = Patterns.periodic_tissue(sample)
    training_instance._cost_values = list(costs)
    training_instance.set_cpp_executable_dir(cpp_executable_dir)
    training_instance.set_tolerance(tolerance)
    cell_parameters_file = "{:04d}.cellParameters.input".format(iteration-1)
    print("Loading cell parameters from: ", cell_parameters_file)
    training_instance.load_cell_parameters(cell_parameters_file)
    os.system("rm {}cellParameters.input".format(run_dir))
    training_instance.set_iter_counter(iteration)
    df = pd.read_csv("{}0000.stresses.csv".format(run_dir))
    training_instance.set_target_cell_to_stress(dict(zip(df['CellID'], df['Target'])))
    training_instance.run_to_max_iters(max_iters=100)