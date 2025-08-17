from toolbox.cSection import makeSampleCrossSection
from toolbox.patterns import Patterns
from toolbox.periodic import PeriodicTissue
from toolbox import stress
import matplotlib.pyplot as plt
from toolbox import stress
import os
import numpy as np
import glob
import copy

####
# Algorithm:
# 1. Set a global target_cell_to_stress dictionary
# 2. (Repeat until) 
# Single iteration:
#   for each target cell:
#       refresh working directory (retain only the minimized file, etc)
#       create new patterns
#       cp working directory to cell_{}_iteration_{} (for future reference)
#   evaluate cost
def clear_directory(dir):
    for file in glob.glob("{}*.bulk.txt".format(dir)):
        os.system("rm {}".format(file))
    for file in glob.glob("{}*.bulk.vtk".format(dir)):
        os.system("rm {}".format(file))
    for file in glob.glob("{}*.stresses.csv".format(dir)):
        os.system("rm {}".format(file))
    for file in glob.glob("{}cellParameters.*.input".format(dir)):
        os.system("rm {}".format(file))
 
def create_random_target_cell_to_stress(training_instance, alpha, n_cells):
    # training_instance.set_tolerance(tolerance)
    training_instance.set_random_target_cells(n_cells=n_cells)
    cellID_to_target_stress = copy.deepcopy(training_instance._target_cell_to_stress)
    for cellID in training_instance._target_cell_to_stress:
        sign = np.random.choice([-1, 1])
        initial_stress = stress.calculate_max_shear_stress(training_instance._config, cellID)
        target_stress = np.round((1+sign*alpha) * initial_stress,3)
        cellID_to_target_stress[cellID] = target_stress
        # training_instance._target_cell_to_stress[cellID] = alpha
        print("Initial and target_stresses for cell {}: {} {}".format(cellID, initial_stress,target_stress))
    return cellID_to_target_stress

def single_cell_solver(**kwargs):
    if not "training_instance" in kwargs:
        raise ValueError("no training instance specified!")
    if not "cellID" in kwargs:
        raise ValueError("no cellID given!")
    if not "target_stress" in kwargs:
        raise ValueError("no target stress specified!")
    
    training_instance = kwargs["training_instance"]
    cellID = kwargs["cellID"]
    target_stress = kwargs["target_stress"]
    training_instance.set_target_cell_to_stress({cellID:target_stress})
    training_instance.run()
    clear_directory(training_instance._dir)

def single_iteration(**kwargs):
    if not "training_instance" in kwargs:
        raise ValueError("no training instance specified!")
    if not "cellID_to_target_stress" in kwargs:
        raise ValueError("no cellID_to_target_stress dict given!")
    
    training_instance = kwargs["training_instance"]
    cellID_to_target_stress = kwargs["cellID_to_target_stress"]
    for cellID, stress in cellID_to_target_stress.items():
        single_cell_solver(training_instance = training_instance, cellID = cellID, target_stress = stress)

def main():
    init_config = "7_0/"
    dir = "7_0/"
    if os.path.isdir(dir):
        os.system("rm -r {}".format(dir))
    os.system("cp -r init/{} {}".format(init_config,dir))
    alpha = 0.05
    n_cells = 2
    min_tolerance = 1e-6
    max_iters = 100
    file = "minimized.txt"
    tissue = PeriodicTissue.from_config(dir,file)
    training_instance = Patterns.periodic_tissue(tissue)
    # training_instance.set_tolerance(tolerance)
    training_instance.minimize_config()
    cellID_to_target_stress = create_random_target_cell_to_stress(training_instance, alpha, n_cells)
    costs = []
    for i in range(max_iters):
        training_instance.set_target_cell_to_stress(cellID_to_target_stress)
        cost = training_instance.evaluate_cost()
        costs.append(cost)
        np.savetxt("{}NetCosts.txt".format(dir), costs, fmt='%.2e')
        if cost < len(cellID_to_target_stress)*min_tolerance:
            print("Training complete!")
            break
        training_instance.set_tolerance(min_tolerance)
        # else:
            # training_instance.set_tolerance(min_tolerance)
        single_iteration(training_instance = training_instance,cellID_to_target_stress = cellID_to_target_stress)

if __name__ == "__main__":
    main()