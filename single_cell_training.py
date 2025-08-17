from toolbox.cSection import makeSampleCrossSection
from toolbox.patterns import Patterns
from toolbox.periodic import PeriodicTissue
from toolbox import stress
import matplotlib.pyplot as plt
from toolbox import stress
import os
import numpy as np

##

def train_random_cells(test_sample, alpha, n_cells, tolerance):
    if os.path.isdir(test_sample):
        os.system("rm -r {}".format(test_sample))
    os.system("cp -r init/{} {}".format(test_sample,test_sample))
    dir = test_sample
    file = "minimized.txt"
    tissue = PeriodicTissue.from_config(dir,file)
    training_instance = Patterns.periodic_tissue(tissue)
    training_instance.minimize_config()
    training_instance.set_tolerance(tolerance)
    # training_instance.set_random_target_cells(n_cells=n_cells)
    training_instance.set_central_target_cells(n_cells=n_cells)
    for cellID in training_instance._target_cell_to_stress:
        # sign = np.random.choice([-1, 1])
        sign = -1
        initial_stress = stress.calculate_max_shear_stress(training_instance._config, cellID)
        training_instance._target_cell_to_stress[cellID] = np.round((1+sign*alpha) * initial_stress,3)
        # training_instance._target_cell_to_stress[cellID] = alpha
        print("Initial stress for cell {}: {}".format(cellID, initial_stress))
    print(training_instance._target_cell_to_stress)
    training_instance.run()

def main():
    test_sample = "7_3/"
    alpha = 0.5
    n_cells = 1
    tolerance = 1e-6
    train_random_cells(test_sample, alpha, n_cells, tolerance)

if __name__ == "__main__":
    main()