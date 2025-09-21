from toolbox.cSection import makeSampleCrossSection
from toolbox.patterns import Patterns
from toolbox.periodic import PeriodicTissue
from toolbox import stress
import matplotlib.pyplot as plt
from toolbox import stress
import os
import numpy as np
import pandas as pd
import random

# def set_central_target_cells(self, n_cells = 1):
#     r_lim = 0.9
#     self._target_cell_to_stress = {}
#     for cellID,cell in self._config.cells_.items():
#         if cell.crossBoundary_:
#             continue
#         r = np.subtract(cell.center_,self._config.sample_center_)
#         r = np.linalg.norm(r)
#         if r<r_lim:
#             self._target_cell_to_stress[cellID] = None
#         if len(self._target_cell_to_stress) == n_cells:
#             break
#     # For checking the above functionality with vtk:
#     for polygonID,polygon in self._config.polygons_.items():
#         polygon.vtk_scalar_ = 0
#     for cellID,_ in self._target_cell_to_stress.items():
#         cell = self._config.cells_[cellID]
#         for polygonID in cell.polygons_:
#             polygon = self._config.polygons_[polygonID]
#             polygon.vtk_scalar_ = 1
#     self._config.write_periodic_vtk(filename = "target_cells.vtk", use_scalar=True)

def set_random_target_cells(training_instance, n_cells = 1, target_stress =1, **kwargs):
    stress_limits = []
    exclude_cells = []
    if "stress_limits" in kwargs:
        stress_limits = kwargs["stress_limits"]
    if "exclude_cells" in kwargs:
        exclude_cells = kwargs["exclude_cells"]
    for polygonID,polygon in training_instance._config.polygons_.items():
        polygon.vtk_scalar_ = 0

    target_cell_to_stress = {}
    while len(target_cell_to_stress)<n_cells:
        cellID = random.choice(list(training_instance._config.cells_.keys()))
        cell = training_instance._config.cells_[cellID]
        if cell.crossBoundary_: 
            continue
        if cellID in target_cell_to_stress:
            continue
        if len(stress_limits):
            if cell.max_shear_stress_ is None:
                cell.max_shear_stress_ = stress.calculate_max_shear_stress(training_instance._config,cellID)
            lower_limit = stress_limits[0]
            upper_limit = stress_limits[1]
            if (cell.max_shear_stress_ < lower_limit):
                continue
            if (cell.max_shear_stress_ > upper_limit):
                continue
        if len(exclude_cells):
            if cellID in exclude_cells:
                continue
        target_cell_to_stress[cellID] = target_stress
        targets_share_polygons = False
        for polygonID in cell.polygons_:
            polygon = training_instance._config.polygons_[polygonID]
            if polygon.vtk_scalar_ == 1:
                targets_share_polygons = True
                break    
        if targets_share_polygons:
            continue

        for polygonID in cell.polygons_:
            polygon = training_instance._config.polygons_[polygonID]
            polygon.vtk_scalar_ = 1
        if len(target_cell_to_stress) == n_cells:
            break

    training_instance.set_target_cell_to_stress(target_cell_to_stress)
    training_instance._config.write_periodic_vtk(filename = "target_cells.vtk", use_scalar=True)
    training_instance._config.write_cell_collection_vtk(list(target_cell_to_stress.keys()),"target_cells_isolated.vtk",use_scalar=False)


def train_random_cells(run_dir, n_cells = 1, target_stress = 1, **kwargs):
    print("Training random cells in directory:", run_dir)
    #default parameters
    cpp_executable_dir = "/home/shabeeb/Projects/tvm-fire/build/"
    tolerance = 1e-8
    max_iters = 100
    stress_limits = []
    exclude_cells = []
    if "stress_limits" in kwargs:
        stress_limits = kwargs["stress_limits"]
    if "cpp_executable_dir" in kwargs:
        cpp_executable_dir = kwargs["cpp_executable_dir"]
    if "tolerance" in kwargs:
        tolerance = kwargs["tolerance"]
    if "max_iters" in kwargs:
        max_iters = kwargs["max_iters"]
    if "exclude_cells" in kwargs:
        exclude_cells = kwargs["exclude_cells"]

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
    set_random_target_cells(
        training_instance = training_instance,
        n_cells = n_cells,
        target_stress = target_stress, 
        exclude_cells = exclude_cells,
        stress_limits = stress_limits)

    training_instance.initialize()
    training_instance.run_to_max_iters(max_iters)
    return training_instance

def resume_run(run_dir, **kwargs):
    print("Resuming runs in directory:", run_dir)
    #default parameters
    cpp_executable_dir = "/home/shabeeb/Projects/tvm-fire/build/"
    tolerance = 1e-8
    max_iters = 100
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
    training_instance.run_to_max_iters(max_iters=max_iters)