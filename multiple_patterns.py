from toolbox.cSection import makeSampleCrossSection
from toolbox.patterns import Patterns
from toolbox.periodic import PeriodicTissue
from toolbox import stress
import matplotlib.pyplot as plt
from toolbox import stress
import os
import numpy as np
import pandas as pd

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
    training_instance.set_random_target_cells(n_cells=n_cells)

    for cellID in training_instance._target_cell_to_stress:
        # sign = np.random.choice([-1, 1])
        sign = 1
        initial_stress = stress.calculate_max_shear_stress(training_instance._config, cellID)
        training_instance._target_cell_to_stress[cellID] = np.round((1+sign*alpha) * initial_stress,3)
        # training_instance._target_cell_to_stress[cellID] = alpha
        print("Initial stress for cell {}: {}".format(cellID, initial_stress))
    print(training_instance._target_cell_to_stress)
    training_instance.run()

def copy_files(source_dir,destination_dir):
    os.system("cp {}minimized.txt {}minimized.txt".format(source_dir, destination_dir))
    os.system("cp {}cellParameters.input {}cellParameters.input".format(source_dir, destination_dir))
def single_iteration (patternA, patternB):
    copy_files(source_dir = patternB._dir, destination_dir = patternA._dir)
    patternA.run()
    copy_files(source_dir = patternA._dir, destination_dir = patternB._dir)
    patternB.run()
def initialize_pattern(dir, tolerance):
    tissue = PeriodicTissue.from_config(dir,"minimized.txt")
    pattern = Patterns.periodic_tissue(tissue)
    pattern.set_tolerance(tolerance)
    df = pd.read_csv("{}0.stresses.csv".format(pattern._dir))
    cellID_to_final_stress = {cellID:stress for cellID,stress in zip(df["CellID"],df["Final"])}
    pattern.set_target_cell_to_stress(cellID_to_final_stress)
    pattern.set_tolerance(tolerance)
    return pattern
def calculate_parameter_space_distance(patternA, patternB):
    cellParametersA = "{}cellParameters.input".format(patternA._dir)
    cellParametersB = "{}cellParameters.input".format(patternB._dir)
    df = pd.read_csv(cellParametersA, sep=" ",header=None)
    cellID_to_s0 = {int(i): [float(s0)] for i, s0 in zip(df[0].to_numpy(), df[2].to_numpy())}
    df = pd.read_csv(cellParametersB, sep=" ",header=None)
    for i, row in df.iterrows():
        cellID = row[0]
        s0 = row[2]
        if cellID in cellID_to_s0:
            cellID_to_s0[cellID].append(s0)
        else:
            raise ValueError("CellID {} not found in first pattern".format(cellID))
    distance = 0
    for cellID, s0s in cellID_to_s0.items():
        # print("CellID: {}, s0s: {}".format(cellID, s0s))
        distance += (s0s[0] - s0s[1])**2
    distance = distance**0.5
    # print(distance)
    return distance
def main():
    dirA = "patternA/"
    dirB = "patternB/"
    for folder in [dirA, dirB]:
        if not os.path.isdir(folder):
            os.system("cp -r tests/{} {}".format(folder, folder))
    tolerance = 1e-5
    patternA = initialize_pattern(dirA, tolerance)
    patternA.set_iter_counter(17)
    patternB = initialize_pattern(dirB, tolerance)
    patternB.set_iter_counter(3)
    print("Initial distance:", calculate_parameter_space_distance(patternA, patternB))
    single_iteration(patternA, patternB)
    print("Final distance:", calculate_parameter_space_distance(patternA, patternB))

if __name__ == "__main__":
    main()