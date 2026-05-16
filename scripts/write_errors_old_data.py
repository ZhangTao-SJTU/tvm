from toolbox.stress import calculate_max_shear_stress
from toolbox.periodic import PeriodicTissue
from toolbox.training import Training
import pandas as pd
import numpy as np
import math
import os
target_stress = 0.23
def write_error_in_dir(dir):
    if not os.path.isfile(dir + "costs.txt"):
        print("costs.txt not found in {}".format(dir))
        return
    costs = np.genfromtxt(dir + "costs.txt")

    if costs.shape == ():
        return
    if costs.shape == (0,):
        return
    if math.isnan(costs[-1]):
        return
    errors = []
    target_cells = pd.read_csv(f"{dir}0000.stresses.csv")["CellID"].to_list()
    for i in range(len(costs)):
        training_instance = Training.from_sample(PeriodicTissue.from_config(dir,f"{i:04d}.bulk.txt"))
        training_instance.load_cell_parameters(f"{i:04d}.cellParameters.input")
        stresses = [calculate_max_shear_stress(training_instance._config, cellID) for cellID in target_cells]
        errors.append(np.mean([abs(s-target_stress)/target_stress for s in stresses]))
    np.savetxt(f"{dir}errors.txt", errors)
    print("wrote errors for {}".format(dir))
    return

def main():
    for i in range(100):
        dir = f"4_cells_mean_l_6/{i:03d}/"
        write_error_in_dir(dir)

if __name__ == "__main__":    
    main()