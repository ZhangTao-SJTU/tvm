from toolbox.pattern_trainer import train_random_cells
import os
import sys
def main():
    input_dir = sys.argv[1]  # e.g., "/home/mameen/init_homogeneous/7_0/"
    run_dir = sys.argv[2]  # e.g., "/home/mameen/mp_1_cell/run_0/"
    alpha = 0.05
    n_cells = 5
    tolerance = 1e-6
    sign = 1
    average_cells_only = True
    # cpp_executable_dir = "/Users/shabeebameen/Projects/tvm-fire/build/"
    cpp_executable_dir = "/home/mameen/tvm/build/"
    os.makedirs(run_dir, exist_ok=True)
    os.system("cp {}conf {}".format(input_dir, run_dir))
    os.system("cp {}minimized.txt {}".format(input_dir, run_dir))
    train_random_cells(
        run_dir, 
        alpha = alpha,
        sign = sign,
        n_cells = n_cells, 
        tolerance = tolerance, 
        average_cells_only = average_cells_only, 
        cpp_executable_dir = cpp_executable_dir)

if __name__ == "__main__":
    main()