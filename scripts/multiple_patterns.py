from toolbox.patterns import Patterns
from toolbox.periodic import PeriodicTissue
from toolbox.stress import calculate_max_shear_stress
from toolbox.pattern_trainer import find_random_target_cells, train_target_cells, resume_run
import os
import numpy as np
import pandas as pd
import sys

class multiple_patterns:
    def __init__(self):
        self._run_dir = None
        self._target_cell_to_stress_A = None
        self._target_cell_to_stress_B = None
        self._target_stress = None
        self._tolerance = 1e-8
        self._n_cells_A = 2
        self._n_cells_B = 2
        self._max_iters = 10
        self._learning_rate = 10

    @classmethod
    def from_dir(cls, dir):
        inst = cls()
        inst._run_dir = dir
        # set cpp_executable_dir based on existing directories
        if os.path.isdir("/Users/shabeebameen/Projects/tvm-fire/build/"):
            cls._cpp_executable_dir = "/Users/shabeebameen/Projects/tvm-fire/build/"
        elif os.path.isdir("/home/shabeeb/Projects/tvm-fire/build/"):
            cls._cpp_executable_dir = "/home/shabeeb/Projects/tvm-fire/build/"
        elif os.path.isdir("/home/mameen/tvm/build/"):
            cls._cpp_executable_dir = "/home/mameen/tvm/build/"

        if os.path.isfile("{}target".format(dir)):
            with open("{}target".format(dir), "r") as f:
                inst._target_stress = float(f.read().strip())
        if os.path.isfile("{}tolerance".format(dir)):
            with open("{}tolerance".format(dir), "r") as f:
                inst._tolerance = float(f.read().strip())
        if os.path.isfile("{}learning_rate".format(dir)):
            with open("{}learning_rate".format(dir), "r") as f:
                inst._learning_rate = float(f.read().strip())
        if os.path.isfile("{}max_iters".format(dir)):
            with open("{}max_iters".format(dir), "r") as f:
                inst._max_iters = int(f.read().strip())
        if os.path.isfile("{}n_cells_A".format(dir)):
            with open("{}n_cells_A".format(dir), "r") as f:
                inst._n_cells_A = int(f.read().strip())
        if os.path.isfile("{}n_cells_B".format(dir)):
            with open("{}n_cells_B".format(dir), "r") as f:
                inst._n_cells_B = int(f.read().strip())
        return inst
  
    # set self._target_cell_to_stress_A and self._target_cell_to_stress_B
    # requires initialization of self._target_stress (can be read from dir/target if using alternate constructor.
    # This will be the uniform target stress for all target cells in both patterns.
    # Save the corresponding vtks and initial_stresses.csv files in the run_dir
    def set_new_uniform_target_stress_patterns(self):
        target_cells_A = find_random_target_cells(self._run_dir, n_cells = self._n_cells_A,output_vtk_file = "target_cells_A.vtk")
        self._target_cell_to_stress_A ={i:self._target_stress for i in target_cells_A}
        target_cells_B = find_random_target_cells(self._run_dir, n_cells = self._n_cells_B, exclude_cells = target_cells_A, output_vtk_file = "target_cells_B.vtk")
        self._target_cell_to_stress_B ={i:self._target_stress for i in target_cells_B}
        # Save initial_stress_A/B.csv files
        tissue = PeriodicTissue.from_config(self._run_dir, "minimized.txt")
        initial_stress_A = {cellID: calculate_max_shear_stress(tissue,cellID)for cellID in self._target_cell_to_stress_A}
        df = pd.DataFrame(list(initial_stress_A.items()), columns=['CellID', 'Current'])        
        df.to_csv("{}initial_stress_A.csv".format(self._run_dir), index=False)
        initial_stress_B = {cellID: calculate_max_shear_stress(tissue,cellID)for cellID in self._target_cell_to_stress_B}
        df = pd.DataFrame(list(initial_stress_B.items()), columns=['CellID', 'Current'])        
        df.to_csv("{}initial_stress_B.csv".format(self._run_dir), index=False)
    def reload_uniform_target_stress_patterns(self):
        target_cells_A = pd.read_csv("{}initial_stress_A.csv".format(self._run_dir))["CellID"].to_numpy()
        self._target_cell_to_stress_A ={int(cellID):self._target_stress for cellID in target_cells_A}
        target_cells_B = pd.read_csv("{}initial_stress_B.csv".format(self._run_dir))["CellID"].to_numpy()
        self._target_cell_to_stress_B ={int(cellID):self._target_stress for cellID in target_cells_B}
    def single_iteration(self):
        # start a new run if no costs.txt file exists
        # Otherwise, use resume_run
        # Either way, first train pattern A for self._max_iters
        if not os.path.isfile("{}costs.txt".format(self._run_dir)):
            print("Starting new run: Training pattern A")
            train_target_cells(self._run_dir, self._target_cell_to_stress_A, learning_rate=self._learning_rate, max_iters=self._max_iters, cpp_executable_dir=self._cpp_executable_dir, tolerance=self._tolerance)
        else:
            print("Training pattern A")
            resume_run(self._run_dir, target_cell_to_stress = self._target_cell_to_stress_A, learning_rate=self._learning_rate, max_iters=self._max_iters, cpp_executable_dir=self._cpp_executable_dir, tolerance=self._tolerance)
        # Now train pattern B for self._max_iters
        print("Training pattern B")
        resume_run(self._run_dir, target_cell_to_stress = self._target_cell_to_stress_B, learning_rate=self._learning_rate, max_iters=self._max_iters, cpp_executable_dir=self._cpp_executable_dir, tolerance=self._tolerance)
    def run(self, iterations):
        print(self._target_cell_to_stress_A)
        print(self._target_cell_to_stress_B)
        print("Target stress:", self._target_stress)
        print("Tolerance:", self._tolerance)
        print("Max iterations:", self._max_iters)
        print("cpp_executable_dir:", self._cpp_executable_dir)
        print("Number of target cells in pattern A:", self._n_cells_A)
        print("Number of target cells in pattern B:", self._n_cells_B)
        for _ in range(iterations):
            self.single_iteration()
    
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

def run(trainer,iterations):
    print(trainer._target_cell_to_stress_A)
    print(trainer._target_cell_to_stress_B)
    print("Target stress:", trainer._target_stress)
    print("Tolerance:", trainer._tolerance)
    print("Max iterations:", trainer._max_iters)
    print("cpp_executable_dir:", trainer._cpp_executable_dir)
    print("Number of target cells in pattern A:", trainer._n_cells_A)
    print("Number of target cells in pattern B:", trainer._n_cells_B)
    for _ in range(iterations):
        trainer.single_iteration()
# def new_multiple_patterns_trainer(run_dir, **kwargs):
#     trainer = multiple_patterns.from_dir(run_dir)
#     trainer.set_uniform_target_stress_patterns()
#     iterations = 100
#     if "iterations" in kwargs:
#         iterations = kwargs["iterations"]
#     run(trainer, iterations)

# def resume_multiple_patterns_trainer(run_dir, **kwargs):
#     trainer = multiple_patterns.from_dir(run_dir)
#     trainer.reload_uniform_target_stress_patterns()
#     iterations = 100
#     if "iterations" in kwargs:
#         iterations = kwargs["iterations"]
#     run(trainer, iterations)
# def new_run(run_dir, **kwargs):
#     iterations = 100
#     if "iterations" in kwargs:
#         iterations = kwargs["iterations"]
#     trainer = multiple_patterns.from_dir(run_dir)
#     trainer.set_uniform_target_stress_patterns()

#     print(trainer._target_cell_to_stress_A)
#     print(trainer._target_cell_to_stress_B)
#     print("Target stress:", trainer._target_stress)
#     print("Tolerance:", trainer._tolerance)
#     print("Max iterations:", trainer._max_iters)
#     print("cpp_executable_dir:", trainer._cpp_executable_dir)
#     print("Number of target cells in pattern A:", trainer._n_cells_A)
#     print("Number of target cells in pattern B:", trainer._n_cells_B)
#     for _ in range(iterations):
#         trainer.single_iteration()

def main():
    if len(sys.argv) < 2:
        print("Usage: python multiple_patterns.py <run_dir>")
        sys.exit(1)
    run_dir = sys.argv[1]
    trainer = multiple_patterns.from_dir(run_dir)
    if not os.path.isfile("{}costs.txt".format(run_dir)):
        trainer.set_new_uniform_target_stress_patterns()
    else:
        trainer.reload_uniform_target_stress_patterns()
    trainer.run(iterations=1000)
if __name__ == "__main__":
    main()