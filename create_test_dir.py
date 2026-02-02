#single pattern 
import os
import numpy as np
from toolbox.periodic import PeriodicTissue
from toolbox.spheroid import Spheroid
from toolbox.patterns import Patterns

def create_test_directory_spheroid():
    test_dir = "tests/1_cell_spheroid_l_6/"
    stresses = np.loadtxt("init/init_homogeneous_6/stresses.txt")
    os.system("rm -rf {}".format(test_dir))
    os.system("cp -r init/init_homogeneous_6/000 {}".format(test_dir))
    os.system("echo {} > {}target".format(np.mean(stresses), test_dir))
    os.system("echo 1e-6 > {}tolerance".format(test_dir))
    os.system("echo 10000 > {}max_iters".format(test_dir))
    os.system("echo 1 > {}clear_interval".format(test_dir))
    os.system("echo 10 > {}learning_rate".format(test_dir))

    sample = PeriodicTissue.from_config(test_dir, "minimized.txt")
    training_instance = Patterns.from_sample(sample)
    training_instance.edit_conf(kv = 10)
    training_instance.find_random_target_cells_in_spheroid(sample,r_sphere = 2, n_cells =1)

def create_test_directory():
    n_cells = 1
    test_dir = "tests/{}_cells_periodic/".format(n_cells)
    print("Creating test directory at:", test_dir)
    stresses = np.loadtxt("init/kv_10_l_4/stresses.txt")
    os.system("rm -rf {}".format(test_dir))
    os.system("cp -r init/kv_10_l_4/000 {}".format(test_dir))
    os.system("echo {} > {}target".format(np.mean(stresses), test_dir))
    os.system("echo 1e-6 > {}tolerance".format(test_dir))
    os.system("echo 10000 > {}max_iters".format(test_dir))
    os.system("echo 10 > {}clear_interval".format(test_dir))
    os.system("echo 10 > {}learning_rate".format(test_dir))

    sample = PeriodicTissue.from_config(test_dir, "minimized.txt")
    Patterns.find_random_target_cells(sample,n_cells=n_cells)

if __name__ == "__main__":
    create_test_directory()