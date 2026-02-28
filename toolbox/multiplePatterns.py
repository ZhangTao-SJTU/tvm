from toolbox.patterns import Patterns
import os
import numpy as np
import pandas as pd


class MultiplePatterns:
    def __init__(self):
        self._run_dir = None
        self._target_cells = None
        self._target_stress = None
        self._convergence_check_interval = 50
        self._max_epochs = 1000
        self._net_error = None
        self._distances = []
        self._epoch = 0
        self._pattern:Patterns = None
        self._subpattern_n_cells = None

    @classmethod
    def from_pattern(cls,pattern:Patterns,subpatterns:list[int]):
        inst = cls()
        inst._run_dir = pattern.get_dir()
        inst._pattern = pattern
        inst._subpattern_n_cells = subpatterns
        inst._target_cells = list(pattern._target_cell_to_stress.keys())
        values = set(pattern._target_cell_to_stress.values())
        if len(values) != 1:
            raise ValueError("MultiplePatterns requires uniform target stress.")
        inst._target_stress = values.pop()
        if sum(subpatterns) != len(inst._target_cells):
            raise ValueError("Subpattern sizes do not match number of target cells.")
        return inst
    def set_target_cells(self,target_cells):
        self._target_cells = target_cells
    def set_target_stress(self,target_stress):
        self._target_stress = target_stress
    def set_subpattern_n_cells(self,subpattern_n_cells):
        self._subpattern_n_cells = subpattern_n_cells

    def single_epoch(self):
        for i, n_cells in enumerate(self._subpattern_n_cells):
            cells_before = sum(self._subpattern_n_cells[:i])
            subpattern = {cellID:self._target_stress for cellID in self._target_cells[cells_before:cells_before+n_cells]}
            print("Training subpattern {}: {}".format(i,subpattern))
            self._pattern.set_target_cell_to_stress(subpattern)
            self._pattern.run_to_max_iters()
            self.write_info(i)
            if self._net_error < self._pattern.get_tolerance():
                print("Converged with net error:", self._net_error)
                break
            # else, we should run the next subpattern. so we set the iter counter for the next iteration...
            self._pattern.set_iter_counter(self._pattern.get_iter_counter()+1)
        self._epoch += 1

    def run(self):
        for _ in range(self._max_epochs):
            self.single_epoch()
            if self._net_error < self._pattern.get_tolerance():
                return
            # check if distance is not changing every convergence_check_interval epochs. 
            # ... But finish 2*convergence_check_interval epochs first.
            if self._epoch>2*self._convergence_check_interval and self._epoch%self._convergence_check_interval == 0:
                if np.allclose(self._distances[-self._convergence_check_interval:], self._distances[-1], rtol = 1e-7):
                    print("Parameter space distance did not change for the last {} subpatterns.".format(self._convergence_check_interval))
                    break

    def evaluate_net_error(self):
        self._pattern.set_target_cell_to_stress({cellID:self._target_stress for cellID in self._target_cells})
        self._net_error = self._pattern.evaluate_cost()

    def evaluate_parameter_space_distance(self):
        file_init = "{}files/0000000.cellParameters.input".format(self._run_dir)
        file_current = "{}cellParameters.input".format(self._run_dir)
        s0_init = pd.read_csv(file_init, sep = " ",header = None)[2].to_numpy()
        s0_current = pd.read_csv(file_current, sep = " ",header = None)[2].to_numpy()
        if not (len(s0_init) == len(s0_current)):
            raise ValueError("Distance cannot be evaluated: cell parameter files have different number of cells")
        self._distances.append(np.sqrt(np.mean((s0_init-s0_current)**2)))
        
    # Moves files to files/ directory and appends info to info.csv file
    def write_info(self,subpattern):
        iter = self._pattern.get_iter_counter()
        overlap = self._pattern.get_last_overlap()
        self.evaluate_net_error()
        self.evaluate_parameter_space_distance()
        info_file = os.path.join(self._run_dir,"info.csv")
        if not os.path.isfile(info_file):
            with open(info_file, "w") as f:
                f.write("Epoch,Iter,Pattern,Error,Overlap,Distance\n")
        with open(info_file, "a") as f:
            f.write("{},{},{},{},{},{}\n".format(self._epoch, iter, subpattern, self._net_error, overlap, self._distances[-1]))
