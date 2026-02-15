from toolbox.patterns import Patterns
from toolbox.periodic import PeriodicTissue
from toolbox.stress import calculate_max_shear_stress
import os
import numpy as np
import pandas as pd
import glob
import sys

class multiple_patterns:
    def __init__(self):
        self._run_dir = None
        self._cpp_executable_dir = None
        self._target_cells = None
        self._frozen_cells = []
        self._target_stress = None
        self._tolerance = None
        self._max_iters = None
        self._learning_rate = None
        self._clear_interval = None
        self._convergence_check_interval = 10
        self._max_epochs = 1000
        self._net_error = None
        self._distance = None
        self._epoch = 0
        self._pattern_trainer:Patterns = None
        self._subpattern_n_cells = None

    @classmethod
    def from_dir(cls, dir):
        inst = cls()
        inst._run_dir = dir
        # set cpp_executable_dir based on existing directories
        if os.path.isdir("/Users/shabeebameen/Projects/tvm-fire/build/"):
            inst._cpp_executable_dir = "/Users/shabeebameen/Projects/tvm-fire/build/"
        elif os.path.isdir("/home/shabeeb/Projects/tvm-fire/build/"):
            inst._cpp_executable_dir = "/home/shabeeb/Projects/tvm-fire/build/"
        elif os.path.isdir("/home/mameen/tvm/build/"):
            inst._cpp_executable_dir = "/home/mameen/tvm/build/"

        if os.path.isfile("{}target".format(dir)):
            inst._target_stress = float(np.loadtxt("{}target".format(dir)))
        if os.path.isfile("{}target_cells.txt".format(dir)):
            with open("{}target_cells.txt".format(dir), "r") as f:
                inst._target_cells = [int(line.strip()) for line in f.readlines()]
            # print(inst._target_cells)
        if os.path.isfile("{}frozen_cells.txt".format(dir)):
            with open("{}frozen_cells.txt".format(dir), 'r') as f:
                inst._frozen_cells = [int(line.strip()) for line in f.readlines()]
        if os.path.isfile("{}tolerance".format(dir)):
            inst._tolerance = np.loadtxt("{}tolerance".format(dir))
        if os.path.isfile("{}learning_rate".format(dir)):
            inst._learning_rate = np.loadtxt("{}learning_rate".format(dir))
        if os.path.isfile("{}max_iters".format(dir)):
            inst._max_iters = int(np.loadtxt("{}max_iters".format(dir)))
        if os.path.isfile("{}clear_interval".format(dir)):
            inst._clear_interval = int(np.loadtxt("{}clear_interval".format(dir)))
        if os.path.isfile("{}subpattern_n_cells".format(dir)):
            with open("{}subpattern_n_cells".format(dir), 'r') as f:
                inst._subpattern_n_cells = [int(line.strip()) for line in f.readlines()]
            if not len(inst._target_cells) == sum(inst._subpattern_n_cells):
                raise ValueError("target cells and subpattern prescription do not match")
        os.makedirs(inst._run_dir + "files/", exist_ok=True)
        inst.initialize_pattern_trainer()
        return inst
  

    def initialize_pattern_trainer(self):
        self._pattern_trainer = Patterns.from_sample(PeriodicTissue.from_config(self._run_dir,"minimized.txt"))
        self._pattern_trainer.set_cpp_executable_dir(self._cpp_executable_dir)
        self._pattern_trainer.set_tolerance(self._tolerance)
        self._pattern_trainer.set_learning_rate(self._learning_rate)
        self._pattern_trainer.set_clear_interval(self._clear_interval)
        self._pattern_trainer.set_frozen_cells(self._frozen_cells)
        self._pattern_trainer.set_target_cell_to_stress({cellID:self._target_stress for cellID in self._target_cells})
        self._pattern_trainer.initialize()

    def single_epoch(self):
        for i, n_cells in enumerate(self._subpattern_n_cells):
            cells_before = sum(self._subpattern_n_cells[:i])
            subpattern = {i:self._target_stress for i in self._target_cells[cells_before:cells_before+n_cells]}
            print("Training subpattern {}: {}".format(i,subpattern))
            self._pattern_trainer.set_target_cell_to_stress(subpattern)
            self._pattern_trainer.run_to_max_iters()
            self.write_info(i)
        self._epoch += 1

    def run(self):
        for i in range(self._max_epochs):
            self.single_epoch()
            if self._net_error < self._tolerance:
                print("Converged with net error:", self._net_error)
                break
            # check if distance is not changing every convergence_check_interval epochs. 
            # ... But finish 2*convergence_check_interval epochs first.
            if i>2*self._convergence_check_interval:
                distances = pd.read_csv("{}info.csv".format(self._run_dir))["Distance"].to_numpy()
                if np.allclose(distances[-self._convergence_check_interval:], distances[-1]):
                    print("Parameter space distance did not change for the last {} epochs.".format(self._convergence_check_interval))
                    break

    def evaluate_net_error(self):
        self._pattern_trainer.set_target_cell_to_stress({cellID:self._target_stress for cellID in self._target_cells})
        self._net_error = self._pattern_trainer.evaluate_cost()

    def evaluate_parameter_space_distance(self):
        file_init = "{}files/0000000.cellParameters.input".format(self._run_dir)
        file_current = "{}cellParameters.input".format(self._run_dir)
        s0_init = pd.read_csv(file_init, sep = " ",header = None)[2].to_numpy()
        s0_current = pd.read_csv(file_current, sep = " ",header = None)[2].to_numpy()
        if not (len(s0_init) == len(s0_current)):
            raise ValueError("Distance cannot be evaluated: cell parameter files have different number of cells")
        self._distance = np.sqrt(np.mean((s0_init-s0_current)**2))
        
    # Moves files to files/ directory and appends info to info.csv file
    def write_info(self,subpattern):
        iter = self._pattern_trainer._iter_counter
        overlap = self._pattern_trainer._q_values[-1]
        self.evaluate_net_error()
        self.evaluate_parameter_space_distance()
        if not os.path.isfile("{}info.csv".format(self._run_dir)):
            with open("{}info.csv".format(self._run_dir), "w") as f:
                f.write("Epoch,Iter,Pattern,Error,Overlap,Distance\n")
        with open("{}info.csv".format(self._run_dir), "a") as f:
            f.write("{},{},{},{},{},{}\n".format(self._epoch, iter, subpattern, self._net_error, overlap, self._distance))

def main():
    if not len(sys.argv) == 2:
        raise ValueError("Usage: python multiple_patterns.py <run_dir>")
    run_dir = sys.argv[1]
    trainer = multiple_patterns.from_dir(run_dir)
    trainer.run()

if __name__ == "__main__":
    main()