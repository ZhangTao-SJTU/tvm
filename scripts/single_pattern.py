from toolbox.patterns import Patterns
from toolbox.periodic import PeriodicTissue
from toolbox.spheroid import Spheroid
import os
import sys
import numpy as np
## A shortcut function that picks and trains random cells.
def train_random_cells(run_dir, n_cells = 1, target_stress = 1, tissue_type = "periodic", **kwargs):
    print("Training random cells in directory:", run_dir)
    #default parameters
    cpp_executable_dir = "/home/shabeeb/Projects/tvm-fire/build/"
    tolerance = 1e-8
    max_iters = 2000
    learning_rate = 10
    clear_interval = 10
    stress_limits = []
    exclude_cells = []
    frozen_cells = []
    if "stress_limits" in kwargs:
        stress_limits = kwargs["stress_limits"]
    if "cpp_executable_dir" in kwargs:
        cpp_executable_dir = kwargs["cpp_executable_dir"]
    if "tolerance" in kwargs:
        tolerance = kwargs["tolerance"]
    if "learning_rate" in kwargs:
        learning_rate = kwargs["learning_rate"]
    if "max_iters" in kwargs:
        max_iters = kwargs["max_iters"]
    if "clear_interval" in kwargs:
        clear_interval = kwargs["clear_interval"]
    if "exclude_cells" in kwargs:
        exclude_cells = kwargs["exclude_cells"]
    if "frozen_cells" in kwargs:
        frozen_cells = kwargs["frozen_cells"]
    

    print("Train Random Cells:\n \tParameters for training:")
    print("cpp_executable_dir:", cpp_executable_dir)
    print("n_cells:", n_cells)
    print("tolerance:", tolerance)
    print("learning_rate:", learning_rate)
    print("clear_interval:", clear_interval)
    print("max_iters:", max_iters)
    print("target stress:", target_stress)
    print("stress_limits:", stress_limits)

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
    training_instance.set_random_target_cells(
        n_cells = n_cells,
        target_stress = target_stress, 
        exclude_cells = exclude_cells,
        stress_limits = stress_limits)

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

def main():
    run_dir = sys.argv[1]
    # Default values
    max_iters = 100000
    n_cells = 1
    stress_limits = []
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
    if os.path.isfile("{}n_cells".format(run_dir)):
        n_cells = int(np.loadtxt("{}n_cells".format(run_dir)))
    if os.path.isfile("{}stress_limits".format(run_dir)):
        stress_limits = list(np.loadtxt("{}stress_limits".format(run_dir)))
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
        train_random_cells(
            run_dir,
            n_cells = n_cells, 
            tolerance = tolerance, 
            learning_rate = learning_rate,
            clear_interval = clear_interval,
            cpp_executable_dir = cpp_executable_dir,
            target_stress = target_stress,
            stress_limits = stress_limits,
            max_iters = max_iters)

if __name__ == "__main__":
    main()