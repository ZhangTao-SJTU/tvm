from toolbox.pattern_trainer import train_random_cells, resume_run
import os
import sys
def main():
    input_dir = sys.argv[1]  # e.g., "/home/mameen/init_homogeneous/7_0/"
    run_dir = sys.argv[2]  # e.g., "/home/mameen/mp_1_cell/run_0/"
    alpha = 0.75
    n_cells = 1
    tolerance = 1e-9
    sign = -1
    average_cells_only = True
    # cpp_executable_dir = "/Users/shabeebameen/Projects/tvm-fire/build/"
    cpp_executable_dir = "/home/mameen/tvm/build/"
    os.makedirs(run_dir, exist_ok=True)
    # os.system("cp {}conf {}".format(input_dir, run_dir))
    # os.system("cp {}minimized.txt {}".format(input_dir, run_dir))
    # train_random_cells(
    #     run_dir, 
    #     alpha = alpha,
    #     sign = sign,
    #     n_cells = n_cells, 
    #     tolerance = tolerance, 
    #     average_cells_only = average_cells_only, 
    #     cpp_executable_dir = cpp_executable_dir)
    resume_run(
        run_dir, 
        tolerance = tolerance, 
        cpp_executable_dir = cpp_executable_dir,
        max_iters = 100)
if __name__ == "__main__":
    main()