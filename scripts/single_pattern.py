from toolbox.pattern_trainer import train_random_cells, resume_run
import os
import sys
import numpy as np

def main():
    run_dir = sys.argv[1]
    # Default values
    max_iters = 10000
    n_cells = 1
    stress_limits = []
    tolerance = 1e-7
    learning_rate = 10

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
    
    if os.path.isfile(run_dir + "0000.stresses.csv"):
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
            cpp_executable_dir = cpp_executable_dir,
            target_stress = target_stress,
            stress_limits = stress_limits,
            max_iters = max_iters)

if __name__ == "__main__":
    main()