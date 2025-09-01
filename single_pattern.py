from toolbox.pattern_trainer import train_random_cells
import os
def main():
    input_dir = "init/7_0/"
    run_dir = "single_cell_test/"
    alpha = 0.025
    n_cells = 1
    tolerance = 1e-6
    if os.path.isdir(run_dir):
        os.system("rm -r {}".format(run_dir))
    os.system("cp -r {} {}".format(input_dir, run_dir))
    train_random_cells(run_dir, alpha, n_cells, tolerance)

if __name__ == "__main__":
    main()