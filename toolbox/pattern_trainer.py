from toolbox.cSection import makeSampleCrossSection
from toolbox.patterns import Patterns
from toolbox.periodic import PeriodicTissue
from toolbox import stress
import matplotlib.pyplot as plt
from toolbox import stress
import os
import numpy as np

def train_random_cells(run_dir, **kwargs):
    #default parameters
    cpp_executable_dir = "/home/shabeeb/Projects/tvm-fire/build/"
    alpha = 0.025
    n_cells = 1
    tolerance = 1e-6
    if "cpp_executable_dir" in kwargs:
        cpp_executable_dir = kwargs["cpp_executable_dir"]
    if "n_cells" in kwargs:
        n_cells = kwargs["n_cells"]
    if "tolerance" in kwargs:
        tolerance = kwargs["tolerance"]
    if "alpha" in kwargs:
        alpha = kwargs["alpha"]
    file = "minimized.txt"
    tissue = PeriodicTissue.from_config(run_dir,file)
    training_instance = Patterns.periodic_tissue(tissue)
    training_instance.set_cpp_executable_dir(cpp_executable_dir)
    training_instance.minimize_config()
    training_instance.set_tolerance(tolerance)
    training_instance.set_random_target_cells(n_cells = n_cells)

    for cellID in training_instance._target_cell_to_stress:
        # sign = np.random.choice([-1, 1])
        sign = 1
        initial_stress = stress.calculate_max_shear_stress(training_instance._config, cellID)
        training_instance._target_cell_to_stress[cellID] = np.round((1+sign*alpha) * initial_stress,3)
        # training_instance._target_cell_to_stress[cellID] = alpha
        print("Initial stress for cell {}: {}".format(cellID, initial_stress))
    print(training_instance._target_cell_to_stress)
    training_instance.initialize()
    training_instance.run()
    return training_instance