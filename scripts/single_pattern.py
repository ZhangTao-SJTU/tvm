from toolbox.patterns import Patterns
from toolbox.periodic import PeriodicTissue
from toolbox.spheroid import Spheroid
import os
import sys
import random
import pandas as pd
from toolbox import stress
import numpy as np

def find_random_target_cells(sample, n_cells = 1, stress_limits = [],exclude_cells = []):
    for polygonID,polygon in sample.polygons_.items():
        polygon.vtk_scalar_ = 0
    target_cells = []
    while len(target_cells)<n_cells:
        cellID = random.choice(list(sample.cells_.keys()))
        cell = sample.cells_[cellID]
        if cell.crossBoundary_: 
            continue
        if sample.tissueType_ == "spheroid" and cell.is_surface_:
            continue
        if sample.tissueType_ == "spheroid" and not cell.type_:
            continue
        if cellID in target_cells:
            continue
        if len(stress_limits):
            cell.max_shear_stress_ = stress.calculate_max_shear_stress(sample,cellID)
            if (cell.max_shear_stress_ < stress_limits[0]):
                continue
            if (cell.max_shear_stress_ > stress_limits[1]):
                continue
        if len(exclude_cells) and cellID in exclude_cells:
                continue
        target_cells.append(cellID)
        for polygonID in cell.polygons_:
            polygon = sample.polygons_[polygonID]
            polygon.vtk_scalar_ = 1
        if len(target_cells) == n_cells:
            break
    np.savetxt("{}target_cells.txt".format(sample.config_dir_), target_cells, fmt='%d')
    sample.write_cell_collection_vtk(target_cells,"target_cells_isolated.vtk",use_scalar=False)
    

def find_target_cells_in_spheroid(sample, r_sphere = 2, n_cells = 1, stress_limits = [],exclude_cells = []):
    sample.calculate_periodic_sample_center()
    frozen_cells = [cellID for cellID,cell in sample.cells_.items() if cell.crossBoundary_ or np.linalg.norm(np.subtract(cell.center_,sample.periodic_sample_center_))>r_sphere]
    np.savetxt("{}frozen_cells.txt".format(sample.config_dir_), frozen_cells, fmt='%d')
    exclude_cells += frozen_cells
    find_random_target_cells(sample, n_cells=n_cells, stress_limits=stress_limits, exclude_cells=exclude_cells)
    # Also record spheroid cell IDs and an initial vtk
    spheroid_cells = [cellID for cellID in sample.cells_ if not cellID in frozen_cells]
    np.savetxt("{}spheroid_cells.txt".format(sample.config_dir_), spheroid_cells, fmt='%d')
    sample.write_cell_collection_vtk(spheroid_cells,"initial_spheroid.vtk",use_scalar=False)


def train_target_cell_to_stress(
        run_dir, 
        target_cell_to_stress,
        tissue_type = "periodic",
        cpp_executable_dir = "/home/shabeeb/Projects/tvm-fire/build/",
        tolerance = 1e-8,
        learning_rate = 10,
        clear_interval = 10,
        max_iters = 2000,
        frozen_cells = []):
    print("Training random cells in directory:", run_dir)
    print("     Parameters for training:")
    print("     cpp_executable_dir:", cpp_executable_dir)
    print("     tolerance:", tolerance)
    print("     learning_rate:", learning_rate)
    print("     clear_interval:", clear_interval)
    print("     max_iters:", max_iters)

    file = "minimized.txt"
    if os.path.isfile("{}minimized.txt".format(run_dir)):
        print("File exists:", "{}minimized.txt".format(run_dir))
    if tissue_type == "periodic":
        tissue = PeriodicTissue.from_config(run_dir,file)
    elif tissue_type == "spheroid":
        tissue = Spheroid.from_config(run_dir,file)
    
    training_instance = Patterns.from_sample(tissue)
    training_instance.set_cpp_executable_dir(cpp_executable_dir)
    training_instance.set_tolerance(tolerance)
    training_instance.set_learning_rate(learning_rate)
    training_instance.set_clear_interval(clear_interval)
    training_instance.set_frozen_cells(frozen_cells)
    training_instance.set_target_cell_to_stress(target_cell_to_stress)
    training_instance.initialize()
    training_instance.run_to_max_iters(max_iters)

def resume_run(run_dir, **kwargs):
    print("Resuming runs in directory:", run_dir)
    #default parameters
    cpp_executable_dir = "/home/shabeeb/Projects/tvm-fire/build/"
    tolerance = 1e-8
    max_iters = 2000
    learning_rate = 10
    target_cell_to_stress = None
    frozen_cells = []
    if "frozen_cells" in kwargs:    
        frozen_cells = kwargs["frozen_cells"]
    if "cpp_executable_dir" in kwargs:
        cpp_executable_dir = kwargs["cpp_executable_dir"]
    if "tolerance" in kwargs:
        tolerance = kwargs["tolerance"]
    if "learning_rate" in kwargs:
        learning_rate = kwargs["learning_rate"]
    if "max_iters" in kwargs:
        max_iters = kwargs["max_iters"]
    if "target_cell_to_stress" in kwargs:
        target_cell_to_stress = kwargs["target_cell_to_stress"]

    print("     Parameters for training:")
    print("     cpp_executable_dir:", cpp_executable_dir)
    print("     tolerance:", tolerance)
    print("     learning_rate:", learning_rate)
    print("     max iters:", max_iters)

    costs = np.loadtxt("{}costs.txt".format(run_dir))
    # if costs[-1]<tolerance:
    #     print("The last iteration already meets the tolerance requirement. No need to resume.")
    #     return
    q_values = np.loadtxt("{}q_values.txt".format(run_dir))
    last_iteration = len(costs) - 1
    config_file = "{:07d}.bulk.txt".format(last_iteration)
    print("Loading configuration from: ", config_file)
    sample = PeriodicTissue.from_config(run_dir,config_file)
    training_instance = Patterns.periodic_tissue(sample)
    training_instance._cost_values = list(costs)
    training_instance._q_values = list(q_values)
    training_instance.set_initial_config(PeriodicTissue.from_config(training_instance._dir,"init_config.txt".format(training_instance._dir)))
    training_instance.set_cpp_executable_dir(cpp_executable_dir)
    training_instance.set_tolerance(tolerance)
    training_instance.set_learning_rate(learning_rate)
    cell_parameters_file = "{:07d}.cellParameters.input".format(last_iteration)
    print("Loading cell parameters from: ", cell_parameters_file)
    training_instance.load_cell_parameters(cell_parameters_file)
    os.system("cp {}{} {}cellParameters.input".format(run_dir,cell_parameters_file,run_dir))
    training_instance.set_iter_counter(last_iteration+1)
    if target_cell_to_stress is None:
        df = pd.read_csv("{}{:07d}.stresses.csv".format(run_dir,0))
        target_cell_to_stress = dict(zip(df['CellID'], df['Target']))
    training_instance.set_target_cell_to_stress(target_cell_to_stress)
    training_instance.set_frozen_cells(frozen_cells)
    training_instance.run_to_max_iters(max_iters=max_iters)

def remove_last_iteration(run_dir):
    costs = np.loadtxt("{}costs.txt".format(run_dir))
    q_values = np.loadtxt("{}q_values.txt".format(run_dir))
    last_iteration = len(costs)-1
    costs = costs[:-1]
    q_values = q_values[:-1]
    np.savetxt("{}costs.txt".format(run_dir), costs,fmt='%.2e')
    np.savetxt("{}q_values.txt".format(run_dir), q_values,fmt='%.4f')
    if os.path.isfile("{}{:07d}.cellParameters.input".format(run_dir,last_iteration)):
        os.remove("{}{:07d}.cellParameters.input".format(run_dir,last_iteration))
    if os.path.isfile("{}{:07d}.stresses.csv".format(run_dir, last_iteration)):
        os.remove("{}{:07d}.stresses.csv".format(run_dir, last_iteration))
    if os.path.isfile("{}{:07d}.bulk.txt".format(run_dir, last_iteration)):
        os.remove("{}{:07d}.bulk.txt".format(run_dir, last_iteration))
    print("Removed iteration {}".format(last_iteration))

def single_pattern(run_dir):
    # Default parameters, can be overridden by files in run_dir
    max_iters = 100000
    tolerance = 1e-7
    learning_rate = 10
    clear_interval = 10
    if os.path.isdir("/Users/shabeebameen/Projects/tvm-fire/build/"):
        cpp_executable_dir = "/Users/shabeebameen/Projects/tvm-fire/build/"
    elif os.path.isdir("/home/shabeeb/Projects/tvm-fire/build/"):
        cpp_executable_dir = "/home/shabeeb/Projects/tvm-fire/build/"
    elif os.path.isdir("/home/mameen/tvm/build/"):
        cpp_executable_dir = "/home/mameen/tvm/build/"
    
    target_stress = np.loadtxt("{}target".format(run_dir))
    with open("{}target_cells.txt".format(run_dir), 'r') as f:
        target_cells = [int(line.strip()) for line in f.readlines()]
    frozen_cells = []
    if os.path.isfile("{}frozen_cells.txt".format(run_dir)):
        with open("{}frozen_cells.txt".format(run_dir), 'r') as f:
            frozen_cells = [int(line.strip()) for line in f.readlines()]
    if os.path.isfile("{}tolerance".format(run_dir)):
        tolerance = np.loadtxt("{}tolerance".format(run_dir))
    if os.path.isfile("{}learning_rate".format(run_dir)):
        learning_rate = np.loadtxt("{}learning_rate".format(run_dir))
    if os.path.isfile("{}max_iters".format(run_dir)):
        max_iters = int(np.loadtxt("{}max_iters".format(run_dir)))
    if os.path.isfile("{}clear_interval".format(run_dir)):
        clear_interval = int(np.loadtxt("{}clear_interval".format(run_dir)))
    if os.path.isfile(run_dir + "{:07d}.stresses.csv".format(0)):
        print("Run already started in dir: {}\n Resuming run from last completed iteration.\n".format(run_dir))
        resume_run(
            run_dir, 
            tolerance = tolerance,
            learning_rate = learning_rate, 
            cpp_executable_dir = cpp_executable_dir,
            max_iters = max_iters)
    else:
        print("Starting a new run in dir: {}".format(run_dir))

        train_target_cell_to_stress(
            run_dir,
            target_cell_to_stress = dict(zip(target_cells, [target_stress]*len(target_cells))),
            frozen_cells=frozen_cells,
            tolerance = tolerance, 
            learning_rate = learning_rate,
            clear_interval = clear_interval,
            cpp_executable_dir = cpp_executable_dir,
            max_iters = max_iters)
 
def main():
    run_dir = sys.argv[1]
    single_pattern(run_dir)
if __name__ == "__main__":
    main()