from toolbox.pattern_trainer import train_random_cells, resume_run
import os
import sys
import numpy as np
def main():
    run_dir = sys.argv[1]
    target_stress = np.loadtxt("{}target".format(run_dir))
    n_cells = np.loadtxt("{}n_cells".format(run_dir))
    tolerance = np.loadtxt("{}tolerance".format(run_dir))
    average_cells_only = True
    # cpp_executable_dir = "/Users/shabeebameen/Projects/tvm-fire/build/"
    cpp_executable_dir = "/home/mameen/tvm/build/"
    train_random_cells(
        run_dir,
        n_cells = n_cells, 
        tolerance = tolerance, 
        average_cells_only = average_cells_only, 
        cpp_executable_dir = cpp_executable_dir,
        target_stress = target_stress)
    
    # resume_run(
    #     run_dir, 
    #     tolerance = tolerance, 
    #     cpp_executable_dir = cpp_executable_dir,
    #     max_iters = 100)


if __name__ == "__main__":
    main()
