from toolbox.cSection import makeSampleCrossSection
from toolbox.patterns import Patterns
from toolbox.periodic import PeriodicTissue
from toolbox.spheroid import Spheroid
from toolbox import stress
import os
import numpy as np
import pandas as pd
import random

# def find_random_target_cells(run_dir, filename = "minimized.txt", n_cells = 1, **kwargs):
#     tissue = PeriodicTissue.from_config(run_dir,filename)
#     training_instance = Patterns.periodic_tissue(tissue)
#     stress_limits = []
#     exclude_cells = []
#     output_vtk_file = "target_cells.vtk"
#     if "stress_limits" in kwargs:
#         stress_limits = kwargs["stress_limits"]
#     if "exclude_cells" in kwargs:
#         exclude_cells = kwargs["exclude_cells"]
#     if "output_vtk_file" in kwargs:
#         output_vtk_file = kwargs["output_vtk_file"]
#     for polygonID,polygon in training_instance._config.polygons_.items():
#         polygon.vtk_scalar_ = 0
#     target_cells = []
#     while len(target_cells)<n_cells:
#         cellID = random.choice(list(training_instance._config.cells_.keys()))
#         cell = training_instance._config.cells_[cellID]
#         if cell.crossBoundary_: 
#             continue
#         if cellID in target_cells:
#             continue
#         if len(stress_limits):
#             if cell.max_shear_stress_ is None:
#                 cell.max_shear_stress_ = stress.calculate_max_shear_stress(training_instance._config,cellID)
#             lower_limit = stress_limits[0]
#             upper_limit = stress_limits[1]
#             if (cell.max_shear_stress_ < lower_limit):
#                 continue
#             if (cell.max_shear_stress_ > upper_limit):
#                 continue
#         if len(exclude_cells):
#             if cellID in exclude_cells:
#                 continue
#         targets_share_polygons = False
#         for polygonID in cell.polygons_:
#             polygon = training_instance._config.polygons_[polygonID]
#             if polygon.vtk_scalar_ == 1:
#                 targets_share_polygons = True
#                 break    
#         if targets_share_polygons:
#             continue
#         for polygonID in cell.polygons_:
#             polygon = training_instance._config.polygons_[polygonID]
#             polygon.vtk_scalar_ = 1
#         target_cells.append(cellID)
#         if len(target_cells) == n_cells:
#             break
#     training_instance._config.write_cell_collection_vtk(cells_array= target_cells,filename = output_vtk_file)
#     return target_cells

    # def set_random_target_cells(training_instance, n_cells = 1, target_stress =1, **kwargs):
    #     stress_limits = []
    #     exclude_cells = []
    #     if "stress_limits" in kwargs:
    #         stress_limits = kwargs["stress_limits"]
    #     if "exclude_cells" in kwargs:
    #         exclude_cells = kwargs["exclude_cells"]
    #     for polygonID,polygon in training_instance._config.polygons_.items():
    #         polygon.vtk_scalar_ = 0

    #     target_cell_to_stress = {}
    #     while len(target_cell_to_stress)<n_cells:
    #         cellID = random.choice(list(training_instance._config.cells_.keys()))
    #         cell = training_instance._config.cells_[cellID]
    #         if cell.crossBoundary_: 
    #             continue
    #         if training_instance._config.tissueType_ == "spheroid" and cell.is_surface_:
    #             continue
    #         if training_instance._config.tissueType_ == "spheroid" and cell.type_ == 0:
    #             continue
    #         if cellID in target_cell_to_stress:
    #             continue
    #         if len(stress_limits):
    #             if cell.max_shear_stress_ is None:
    #                 cell.max_shear_stress_ = stress.calculate_max_shear_stress(training_instance._config,cellID)
    #             lower_limit = stress_limits[0]
    #             upper_limit = stress_limits[1]
    #             if (cell.max_shear_stress_ < lower_limit):
    #                 continue
    #             if (cell.max_shear_stress_ > upper_limit):
    #                 continue
    #         if len(exclude_cells):
    #             if cellID in exclude_cells:
    #                 continue
    #         target_cell_to_stress[cellID] = target_stress
    #         targets_share_polygons = False
    #         for polygonID in cell.polygons_:
    #             polygon = training_instance._config.polygons_[polygonID]
    #             if polygon.vtk_scalar_ == 1:
    #                 targets_share_polygons = True
    #                 break    
    #         if targets_share_polygons:
    #             continue
    #         for polygonID in cell.polygons_:
    #             polygon = training_instance._config.polygons_[polygonID]
    #             polygon.vtk_scalar_ = 1
    #         if len(target_cell_to_stress) == n_cells:
    #             break

    #     training_instance.set_target_cell_to_stress(target_cell_to_stress)
    #     if training_instance._config.tissueType_ == "periodic":
    #         training_instance._config.write_periodic_vtk(filename = "target_cells.vtk", use_scalar=True)
    #     # elif training_instance._config.tissueType_ == "spheroid":
    #     #     training_instance._config.write_spheroid_vtk(cells_array= list(target_cell_to_stress.keys()),filename = "target_cells.vtk",)
    #     training_instance._config.write_cell_collection_vtk(list(target_cell_to_stress.keys()),"target_cells_isolated.vtk",use_scalar=False)

# A target cell trainer. Input run directory and a dictionary of target cells to target stress values.
def train_target_cells(run_dir, target_cell_to_stress, **kwargs):
    print("Training target cells in directory:", run_dir)
    #default parameters
    cpp_executable_dir = "/home/shabeeb/Projects/tvm-fire/build/"
    tolerance = 1e-8
    max_iters = 2000
    learning_rate = 10
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


    print("Parameters for training:")
    print("cpp_executable_dir:", cpp_executable_dir)
    print("tolerance:", tolerance)
    print("learning_rate:", learning_rate)
    print("max_iters:", max_iters)
    print("target_cell_to_stress:", target_cell_to_stress)

    file = "minimized.txt"
    if os.path.isfile("{}minimized.txt".format(run_dir)):
        print("File exists:", "{}minimized.txt".format(run_dir))
    tissue = PeriodicTissue.from_config(run_dir,file)
    training_instance = Patterns.periodic_tissue(tissue)
    training_instance.set_cpp_executable_dir(cpp_executable_dir)
    training_instance.minimize_config()
    training_instance.set_tolerance(tolerance)
    training_instance.set_learning_rate(learning_rate)
    training_instance.set_target_cell_to_stress(target_cell_to_stress)
    training_instance.set_frozen_cells(frozen_cells)
    training_instance.initialize()
    training_instance.run_to_max_iters(max_iters)

# ## A shortcut function that picks and trains random cells.
# def train_random_cells(run_dir, n_cells = 1, target_stress = 1, tissue_type = "periodic", **kwargs):
#     print("Training random cells in directory:", run_dir)
#     #default parameters
#     cpp_executable_dir = "/home/shabeeb/Projects/tvm-fire/build/"
#     tolerance = 1e-8
#     max_iters = 2000
#     learning_rate = 10
#     stress_limits = []
#     exclude_cells = []
#     frozen_cells = []
#     if "stress_limits" in kwargs:
#         stress_limits = kwargs["stress_limits"]
#     if "cpp_executable_dir" in kwargs:
#         cpp_executable_dir = kwargs["cpp_executable_dir"]
#     if "tolerance" in kwargs:
#         tolerance = kwargs["tolerance"]
#     if "learning_rate" in kwargs:
#         learning_rate = kwargs["learning_rate"]
#     if "max_iters" in kwargs:
#         max_iters = kwargs["max_iters"]
#     if "exclude_cells" in kwargs:
#         exclude_cells = kwargs["exclude_cells"]
#     if "frozen_cells" in kwargs:
#         frozen_cells = kwargs["frozen_cells"]
    

#     print("Train Random Cells:\n \tParameters for training:")
#     print("cpp_executable_dir:", cpp_executable_dir)
#     print("n_cells:", n_cells)
#     print("tolerance:", tolerance)
#     print("learning_rate:", learning_rate)
#     print("max_iters:", max_iters)
#     print("target stress:", target_stress)
#     print("stress_limits:", stress_limits)

#     file = "minimized.txt"
#     if os.path.isfile("{}minimized.txt".format(run_dir)):
#         print("File exists:", "{}minimized.txt".format(run_dir))
#     if tissue_type == "periodic":
#         tissue = PeriodicTissue.from_config(run_dir,file)
#         training_instance = Patterns.periodic_tissue(tissue)
#     elif tissue_type == "spheroid":
#         tissue = Spheroid.from_config(run_dir,file)
#         training_instance = Patterns.spheroid(tissue)
#     training_instance.set_cpp_executable_dir(cpp_executable_dir)
#     training_instance.minimize_config()
#     training_instance.set_tolerance(tolerance)
#     training_instance.set_learning_rate(learning_rate)
#     training_instance.set_frozen_cells(frozen_cells)
#     training_instance.set_random_target_cells(
#         n_cells = n_cells,
#         target_stress = target_stress, 
#         exclude_cells = exclude_cells,
#         stress_limits = stress_limits)

#     training_instance.initialize()
#     training_instance.run_to_max_iters(max_iters)

# def resume_run(run_dir, **kwargs):
#     print("Resuming runs in directory:", run_dir)
#     #default parameters
#     cpp_executable_dir = "/home/shabeeb/Projects/tvm-fire/build/"
#     tolerance = 1e-8
#     max_iters = 2000
#     learning_rate = 10
#     target_cell_to_stress = None
#     frozen_cells = []
#     if "frozen_cells" in kwargs:    
#         frozen_cells = kwargs["frozen_cells"]
#     if "cpp_executable_dir" in kwargs:
#         cpp_executable_dir = kwargs["cpp_executable_dir"]
#     if "tolerance" in kwargs:
#         tolerance = kwargs["tolerance"]
#     if "learning_rate" in kwargs:
#         learning_rate = kwargs["learning_rate"]
#     if "max_iters" in kwargs:
#         max_iters = kwargs["max_iters"]
#     if "target_cell_to_stress" in kwargs:
#         target_cell_to_stress = kwargs["target_cell_to_stress"]


#     print("Parameters for training:")
#     print("cpp_executable_dir:", cpp_executable_dir)
#     print("tolerance:", tolerance)
#     print("learning_rate:", learning_rate)
#     print("max iters:", max_iters)
    
#     costs = np.loadtxt("{}costs.txt".format(run_dir))
#     # if costs[-1]<tolerance:
#     #     print("The last iteration already meets the tolerance requirement. No need to resume.")
#     #     return
#     q_values = np.loadtxt("{}q_values.txt".format(run_dir))
#     last_iteration = len(costs) - 1
#     config_file = "{:07d}.bulk.txt".format(last_iteration)
#     print("Loading configuration from: ", config_file)
#     sample = PeriodicTissue.from_config(run_dir,config_file)
#     training_instance = Patterns.periodic_tissue(sample)
#     training_instance._cost_values = list(costs)
#     training_instance._q_values = list(q_values)
#     training_instance.set_initial_config(PeriodicTissue.from_config(training_instance._dir,"init_config.txt".format(training_instance._dir)))
#     training_instance.set_cpp_executable_dir(cpp_executable_dir)
#     training_instance.set_tolerance(tolerance)
#     training_instance.set_learning_rate(learning_rate)
#     cell_parameters_file = "{:07d}.cellParameters.input".format(last_iteration)
#     print("Loading cell parameters from: ", cell_parameters_file)
#     training_instance.load_cell_parameters(cell_parameters_file)
#     os.system("cp {}{} {}cellParameters.input".format(run_dir,cell_parameters_file,run_dir))
#     training_instance.set_iter_counter(last_iteration+1)
#     if target_cell_to_stress is None:
#         df = pd.read_csv("{}{:07d}.stresses.csv".format(run_dir,0))
#         target_cell_to_stress = dict(zip(df['CellID'], df['Target']))
#     training_instance.set_target_cell_to_stress(target_cell_to_stress)
#     training_instance.set_frozen_cells(frozen_cells)
#     training_instance.run_to_max_iters(max_iters=max_iters)

# def remove_last_iteration(run_dir):
#     costs = np.loadtxt("{}costs.txt".format(run_dir))
#     q_values = np.loadtxt("{}q_values.txt".format(run_dir))
#     last_iteration = len(costs)-1
#     costs = costs[:-1]
#     q_values = q_values[:-1]
#     np.savetxt("{}costs.txt".format(run_dir), costs,fmt='%.2e')
#     np.savetxt("{}q_values.txt".format(run_dir), q_values,fmt='%.4f')
#     if os.path.isfile("{}{:07d}.cellParameters.input".format(run_dir,last_iteration)):
#         os.remove("{}{:07d}.cellParameters.input".format(run_dir,last_iteration))
#     if os.path.isfile("{}{:07d}.stresses.csv".format(run_dir, last_iteration)):
#         os.remove("{}{:07d}.stresses.csv".format(run_dir, last_iteration))
#     if os.path.isfile("{}{:07d}.bulk.txt".format(run_dir, last_iteration)):
#         os.remove("{}{:07d}.bulk.txt".format(run_dir, last_iteration))
#     print("Removed iteration {}".format(last_iteration))

