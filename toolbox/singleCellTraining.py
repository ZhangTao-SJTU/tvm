from toolbox.training import Training
from toolbox import stress
import numpy as np
import os
import random
import pandas as pd

# This single cell trainer repackages the patterns class.
# Use cases:
# 1. Training a single cell to reach a target stress.
# 2. Training a multiple cell pattern in parallel.

class singleCellTraining(Training):
    def __init__(self):
        super().__init__()
        self._target_cell_ID = None
        self._target_stress = None
        self._clamping_max_iters = 10
        self._clamping_correction_factor = 10
        self._clamping_FIRE_only = True
        self._clamping_tolerance = 1e-3
        
    @classmethod
    def periodic_tissue(cls,tissue):
        inst = super().periodic_tissue(tissue)
        inst._modified_cells = list(inst._config.cells_.keys())
        return inst
    def set_clamping_FIRE_only(self, clamping_FIRE_only):
        self._clamping_FIRE_only = clamping_FIRE_only
    def set_clamping_max_iters(self, clamping_max_iters):
        self._clamping_max_iters = clamping_max_iters
    def set_clamping_correction_factor(self, clamping_correction_factor):
        self._clamping_correction_factor = clamping_correction_factor
    def set_clamping_tolerance(self, tol):
        self._clamping_tolerance = tol
    def set_target_cell_ID(self,target_cell_ID):
        self._target_cell_ID = target_cell_ID
    def set_target_stress(self,target_stress):
        self._target_stress = target_stress

    def calculate_max_shear_stress(self):
        cell = self._config.cells_[self._target_cell_ID]
        cell.max_shear_stress_ = stress.calculate_max_shear_stress(self._config, self._target_cell_ID)

    def evaluate_cost(self):
        multiplier = 1
        self.calculate_max_shear_stress()
        cell = self._config.cells_[self._target_cell_ID]
        cost = multiplier * (cell.max_shear_stress_ - self._target_stress) ** 2
        return cost
    
    def single_iteration(self):
        # Starting with a minimized config with loaded cell properties...
        print("\n\n=====================================\n\n")  
        print("Starting iteration: {:d}".format(self._iter_counter))
        print("\n\n=====================================\n\n")
        print("Number of modifiable cells: {:d}".format(len(self._modified_cells)))
        print("\n\n-------------------------------------")
        print("Step 1: Evaluate and store the current (free state) areas of hidden (non-target) cells")
        print("-------------------------------------\n\n")
        # The stored cell areas will be used to calculate learning DOF changes
        free_state_areas = {}
        for cellID,cell in self._config.cells_.items():
            if cellID == self._target_cell_ID:
                continue
            free_state_areas[cellID] = cell.surface_area_

        print("\n\n-------------------------------------")
        print("Step 2: CLAMPING")
        print("-------------------------------------\n\n")
        self.clamp()

        print("\n\n-------------------------------------")
        print("Step 3: Use the clamped state areas to update all learning degrees of freedom.")
        print("-------------------------------------\n\n")
        for cellID in free_state_areas:
            cell = self._config.cells_[cellID]
            del_area = cell.surface_area_ - free_state_areas[cellID]
            s0_change = self._learning_rate * del_area
            cell.s0_ -= s0_change
        
        print("\n\n-------------------------------------")
        print("Step 4: UNCLAMP target cell; write cell parameters")
        print("-------------------------------------\n\n")
        self._config.cells_[self._target_cell_ID].s0_ = self._config.s0_
        self.write_cell_parameters()
        self.minimize_config()
        """
        self.write_configuration(filename = "{}.bulk.txt".format(self._iter_counter))
        for cellID,cell in self._config.cells_.items():
            for polygonID in cell.polygons_:
                self._config.polygons_[polygonID].vtk_scalar_ = cell.s0_
        self._config.write_periodic_vtk(filename = "{}.bulk.vtk".format(self._iter_counter),use_scalar=True)
        os.system("cp {}cellParameters.input {}cellParameters.{}.input".format(self._dir,self._dir,self._iter_counter))
        """
    # The goal of clamping is to ensure that the target cell reaches the final target stress as
    # in self._target__stress, upto clamping tolerance.
    # Holding the configuration fixed, we iteratively adjust s0 of the target cells 
    # such that stress overshoots/undershoots the target stress and energy minimize the config at each iteration.    
    # Minimizing changes the surface area and hence the stress, so process is repeated.
    def clamp(self):
        for iter in range(self._clamping_max_iters):
            print("Clamping iteration {}".format(iter))
            current_stress = stress.calculate_max_shear_stress(self._config, self._target_cell_ID)
            if abs(current_stress - self._target_stress) < self._clamping_tolerance:
                print("Clamping successful to tolerance")
                return
            current_stress = stress.calculate_max_shear_stress(self._config, self._target_cell_ID)
            # The temporary target stress to be solved for should be an overshoot/undershoot of self._target_stress
            temp_target_stress = self._target_stress * (1 + self._clamping_correction_factor * (self._target_stress - current_stress))
            # The temporary target stress should be positive, 
            # so we take the absolute value in case we hit a really small number
            temp_target_stress = abs(temp_target_stress)
            print("Final Target Stress: {} Temporary Target Stress: {}, Current stress: {}, s0: {}".format(
                 self._target_stress, temp_target_stress, current_stress, self._config.cells_[self._target_cell_ID].s0_))
            self.solve_cell_s0_for_target_stress(self._target_cell_ID, temp_target_stress)
            current_stress = stress.calculate_max_shear_stress(self._config, self._target_cell_ID)
            self.write_cell_parameters()
            self.minimize_config(FIRE_only = self._clamping_FIRE_only)

    # Binary search for s0 that produces the right stress        
    def solve_cell_s0_for_target_stress(self, cellID, target_stress):
        cell = self._config.cells_[cellID]
        upper_limit = 5.3
        lower_limit = 4.6
        s0_guesses = {i:None for i in np.linspace(lower_limit, upper_limit, 1000)}
        for s0 in s0_guesses:
            cell.s0_ = s0
            current_stress = stress.calculate_max_shear_stress(self._config, cellID)
            s0_guesses[s0] = current_stress 
        init_guess = min(s0_guesses, key=lambda x: abs(s0_guesses[x] - target_stress))
        cell.s0_ = init_guess
        return

    def run(self):
        #store initial cell parameters
        os.system("cp {}cellParameters.input {}cellParameters.init.input".format(self._dir,self._dir))
        os.system("cp {}minimized.txt {}init_config.txt".format(self._dir,self._dir))
        self._cost_values = []
        cost = self.evaluate_cost()
        self._cost_values.append(cost)
        initial_stress = stress.calculate_max_shear_stress(self._config, self._target_cell_ID)
        print("Initial Cost: {:.2e}".format(cost))
        while cost > self._tolerance:
            self.set_clamping_tolerance(cost * 1)
            self.single_iteration()
            cost = self.evaluate_cost()
            self._cost_values.append(cost)
            print("Iteration: {:d}, Cost: {:.2e}".format(self._iter_counter,cost))
            # Rewrite costs.txt
            np.savetxt("{}costs.txt".format(self._dir), self._cost_values, fmt='%.2e')
            # Rewrite stresses.csv
            current_stress = stress.calculate_max_shear_stress(self._config,self._target_cell_ID)
            results = {"CellID": self._target_cell_ID,
                    "Target": self._target_stress,
                    "Initial": initial_stress,
                    "Current": current_stress}
            df = pd.DataFrame(results)
            df.to_csv("{}{}.stresses.csv".format(self._dir, self._iter_counter), index=False)
            self._iter_counter += 1
        print("Optimization finished at iteration: {:d}".format(self._iter_counter-1))
        print("Final cost: {:.2e}".format(cost))