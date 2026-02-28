import glob
from toolbox.patterns import Patterns
from toolbox.periodic import PeriodicTissue
from toolbox.spheroid import Spheroid
import os
import sys
import random
import pandas as pd
from toolbox import stress
import numpy as np

def train_target_cell_to_stress(
        run_dir, 
        target_cell_to_stress,
        tissue,
        cpp_executable_dir,
        tolerance,
        learning_rate,
        clear_interval,
        max_iters,
        frozen_cells):
    print("Training random cells in directory:", run_dir)
    print("     Parameters for training:")
    print("     cpp_executable_dir:", cpp_executable_dir)
    print("     tolerance:", tolerance)
    print("     learning_rate:", learning_rate)
    print("     clear_interval:", clear_interval)
    print("     max_iters:", max_iters)

    training_instance = Patterns.from_sample(tissue)
    training_instance.set_cpp_executable_dir(cpp_executable_dir)
    training_instance.set_tolerance(tolerance)
    training_instance.set_learning_rate(learning_rate)
    training_instance.set_clear_interval(clear_interval)
    training_instance.set_frozen_cells(frozen_cells)
    training_instance.set_target_cell_to_stress(target_cell_to_stress)
    training_instance.initialize()
    training_instance.run_to_max_iters(max_iters)


def single_pattern(run_dir):
    # Default parameters, can be overridden by files in run_dir
    
    tolerance = 1e-7
    learning_rate = 10
    clear_interval = 10
    max_iters = 100000
    frozen_cells = []
    
    if os.path.isdir("/Users/shabeebameen/Projects/tvm-fire/build/"):
        cpp_executable_dir = "/Users/shabeebameen/Projects/tvm-fire/build/"
    elif os.path.isdir("/home/shabeeb/Projects/tvm-fire/build/"):
        cpp_executable_dir = "/home/shabeeb/Projects/tvm-fire/build/"
    elif os.path.isdir("/home/mameen/tvm/build/"):
        cpp_executable_dir = "/home/mameen/tvm/build/"
    
    target_stress = np.loadtxt("{}target".format(run_dir))
    with open("{}target_cells.txt".format(run_dir), 'r') as f:
        target_cells = [int(line.strip()) for line in f.readlines()]
    
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

    # if os.path.isfile(run_dir + "{:07d}.stresses.csv".format(0)):
    #     print("Run already started in dir: {}\n Resuming run from last completed iteration.\n".format(run_dir))
    #     resume_run(
    #         run_dir, 
    #         tolerance = tolerance,
    #         learning_rate = learning_rate, 
    #         cpp_executable_dir = cpp_executable_dir,
    #         max_iters = max_iters)
    # else:
    #     print("Starting a new run in dir: {}".format(run_dir))
    #     train_target_cell_to_stress(
    #         run_dir,
    #         target_cell_to_stress = dict(zip(target_cells, [target_stress]*len(target_cells))),
    #         frozen_cells=frozen_cells,
    #         tolerance = tolerance, 
    #         learning_rate = learning_rate,
    #         clear_interval = clear_interval,
    #         cpp_executable_dir = cpp_executable_dir,
    #         max_iters = max_iters)

    print("Starting a new run in dir: {}".format(run_dir))
    train_target_cell_to_stress(
        run_dir,
        target_cell_to_stress = dict(zip(target_cells, [target_stress]*len(target_cells))),
        tissue = PeriodicTissue.from_config(run_dir,"minimized.txt"),
        cpp_executable_dir = cpp_executable_dir,
        tolerance = tolerance, 
        learning_rate = learning_rate,
        clear_interval = clear_interval,
        max_iters = max_iters,
        frozen_cells = frozen_cells)

def main():
    run_dir = sys.argv[1]
    single_pattern(run_dir)

if __name__ == "__main__":
    main()