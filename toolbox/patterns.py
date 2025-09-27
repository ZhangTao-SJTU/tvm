from toolbox.periodic import PeriodicTissue
from toolbox.training import Training
from toolbox.overlap import calculate_Q
from toolbox import stress
import numpy as np
import os
import random
import pandas as pd

class Patterns(Training):
    def __init__(self):
        super().__init__()
        self._target_cell_to_stress = None
        self._clamping_max_iters = 20
        self._clamping_s0_lower_limit = 4.6
        self._clamping_s0_upper_limit = 5.3
        self._clamping_correction_factor = 0.1
        self._clamping_FIRE_only = True
        # self._clamping_tolerance = 1e-7
        
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
    # def set_clamping_tolerance(self, tol):
    #     self._clamping_tolerance = tol
    def set_target_cell_to_stress(self,target_cell_to_stress):
        self._target_cell_to_stress = target_cell_to_stress
        # For checking the above functionality with vtk:
        for polygonID,polygon in self._config.polygons_.items():
            polygon.vtk_scalar_ = 0
        for cellID in self._target_cell_to_stress:
            cell = self._config.cells_[cellID]
            for polygonID in cell.polygons_:
                polygon = self._config.polygons_[polygonID]
                polygon.vtk_scalar_ = 1
        self._config.write_periodic_vtk(filename = "target_cells.vtk", use_scalar=True)
                
    def calculate_max_shear_stresses(self):
        for cellID in self._target_cell_to_stress:
            cell = self._config.cells_[cellID]
            cell.max_shear_stress_ = stress.calculate_max_shear_stress(self._config, cellID)

    def evaluate_cost(self):
        multiplier = 1
        cost = 0
        self.calculate_max_shear_stresses()
        for cellID, target_stress in self._target_cell_to_stress.items():
            cell = self._config.cells_[cellID]
            cost += multiplier * (cell.max_shear_stress_ - target_stress) ** 2
        return cost

    def single_iteration(self):
        if not len(self._target_cell_to_stress):
            print("No target cells, iteration terminated.")
            return
        free_state_areas = {}
        # Starting with a minimized config with loaded cell properties...
        print("\n\n=====================================\n\n")  
        print("Starting iteration: {:d}".format(self._iter_counter))
        print("\n\n=====================================\n\n")
        print("Number of modifiable cells: {:d}".format(len(self._modified_cells)))
        
        print("\n\n-------------------------------------")
        print("Step 1: Evaluate and store the current (free state) areas of hidden (non-target) cells")
        print("-------------------------------------\n\n")
        
        # The stored cell areas will be used to calculate learning DOF changes
        for cellID,cell in self._config.cells_.items():
            if cellID in self._target_cell_to_stress:
                continue
            free_state_areas[cellID] = cell.surface_area_

        print("\n\n-------------------------------------")
        print("Step 2: CLAMPING")
        print("-------------------------------------\n\n")

        self.clamp_target_cells()

        print("\n\n-------------------------------------")
        print("Step 3: Use the clamped state areas to update all learning degrees of freedom.")
        print("-------------------------------------\n\n")

        for cellID in free_state_areas:
            cell = self._config.cells_[cellID]
            del_area = cell.surface_area_ - free_state_areas[cellID]
            s0_change = self._learning_rate * del_area
            ## TESTING
            cell.s0_ -= s0_change
        
        print("\n\n-------------------------------------")
        print("Step 4: UNCLAMP target cells; write cell parameters")
        print("-------------------------------------\n\n")
        for cellID in self._target_cell_to_stress:
            cell = self._config.cells_[cellID]
            cell.s0_ = self._config.s0_
        self.write_cell_parameters()
        self.minimize_config()
        # Step 5: Logging, etc
        self.write_configuration(filename = "{:04d}.bulk.txt".format(self._iter_counter))
        for cellID,cell in self._config.cells_.items():
            for polygonID in cell.polygons_:
                self._config.polygons_[polygonID].vtk_scalar_ = cell.s0_
        # self._config.write_periodic_vtk(filename = "{:04d}.bulk.vtk".format(self._iter_counter),use_scalar=True)
        os.system("cp {}cellParameters.input {}{:04d}.cellParameters.input".format(
            self._dir,self._dir,self._iter_counter))

    # The goal of clamping is to ensure that the target cells reach the final target stress as
    # in self._target_cell_to_stress, upto clamping tolerance.
    # Holding the configuration fixed, we iteratively adjust s0 of the target cells 
    # such that stress overshoots/undershoots the target stress and energy minimize the config at each iteration.    
    # Minimizing changes the surface area and hence the stress, so process is repeated.
    def clamp_target_cells(self):
        self.write_cell_parameters()
        cell_parameters_file = self._dir+"cellParameters.input"
        for iter in range(self._clamping_max_iters):
            print("Clamping iteration {}".format(iter))
            needs_clamping = []
            for cellID, final_target_stress in self._target_cell_to_stress.items():
                cell = self._config.cells_[cellID]
                current_stress = stress.calculate_max_shear_stress(self._config, cellID)
                # if abs(current_stress - final_target_stress) > self._clamping_tolerance:
                if abs(current_stress - final_target_stress) > self._tolerance:
                    needs_clamping.append(cellID)
            if not len(needs_clamping):
                print("Clamping successful to tolerance")
                return
            pre_clamping_s0 = pd.read_csv(cell_parameters_file,header=None,sep=" ")[2]
            for cellID in needs_clamping:
                cell = self._config.cells_[cellID]
                final_target_stress = self._target_cell_to_stress[cellID]
                current_stress = stress.calculate_max_shear_stress(self._config, cellID)
                # The temporary target stress to be solved for should be an overshoot/undershoot of self._target_stress
                temp_target_stress = final_target_stress + self._clamping_correction_factor * (final_target_stress - current_stress)
                # The temporary target stress should be positive,
                # so if the overshooting makes it negative, avoid it.
                if temp_target_stress<0:
                    temp_target_stress = final_target_stress
                print("Cell: {}, Final Target Stress: {} Temporary Target Stress: {}, Current stress: {}, s0: {}".format(
                    cellID, final_target_stress, temp_target_stress, current_stress, cell.s0_))
                self.solve_cell_s0_for_target_stress(cellID, temp_target_stress)
                current_stress = stress.calculate_max_shear_stress(self._config, cellID)
            self.write_cell_parameters()
            self.minimize_config(FIRE_only = self._clamping_FIRE_only)
            post_clamping_s0 = pd.read_csv(cell_parameters_file,header=None,sep=" ")[2]
            if post_clamping_s0.equals(pre_clamping_s0):
                print("Clamping terminated: doesn't change s0")
                return
            
    # Binary search for s0 that produces the right stress        
    # def solve_cell_s0_for_target_stress(self, cellID, target_stress):
    #     cell = self._config.cells_[cellID]
    #     upper_limit = self._clamping_s0_upper_limit
    #     lower_limit = self._clamping_s0_lower_limit
    #     steps = self._clamping_steps
    #     s0_guesses = {i:None for i in np.linspace(lower_limit, upper_limit, steps)}
    #     for s0 in s0_guesses:
    #         cell.s0_ = s0
    #         current_stress = stress.calculate_max_shear_stress(self._config, cellID)
    #         s0_guesses[s0] = current_stress 
    #     init_guess = min(s0_guesses, key=lambda x: abs(s0_guesses[x] - target_stress))
    #     cell.s0_ = init_guess
    #     return
    def solve_cell_s0_for_target_stress(self,cellID, target_stress):
        cell = self._config.cells_[cellID]
        def clamping_error(s0):
            cell.s0_ = s0
            return stress.calculate_max_shear_stress(self._config,cellID)-target_stress
        upper_limit = self._clamping_s0_upper_limit
        lower_limit = self._clamping_s0_lower_limit

        root_interval = None
        s0_to_clamping_error = {s0:clamping_error(s0) for s0 in np.linspace(lower_limit, upper_limit,15)}
        for i, s0 in enumerate(list(s0_to_clamping_error.keys())):
            if i == len(s0_to_clamping_error)-1:
                continue
            if np.sign(s0_to_clamping_error[s0]) == np.sign(list(s0_to_clamping_error.values())[i+1]):
                continue
            root_interval = [s0,list(s0_to_clamping_error.keys())[i+1]]
            break
        
        if root_interval is None:
            print("No root interval found, picking closest guess")
            min_guess = min(s0_to_clamping_error, key=lambda x: abs(s0_to_clamping_error[x]))
            cell.s0_ = min_guess
            return
        # Binary search within the root interval
        # while abs(root_interval[1]-root_interval[0])>self._clamping_tolerance:
        while abs(root_interval[1]-root_interval[0])>self._tolerance:
            mid = (root_interval[0]+root_interval[1])/2
            if np.sign(clamping_error(mid)) == np.sign(clamping_error(root_interval[0])):
                root_interval[0] = mid
            else:
                root_interval[1] = mid
        cell.s0_ = (root_interval[0]+root_interval[1])/2
        print("Binary search complete, found s0:", cell.s0_)

    def initialize(self):
        os.system("cp {}cellParameters.input {}cellParameters.init.input".format(self._dir,self._dir))
        os.system("cp {}minimized.txt {}init_config.txt".format(self._dir,self._dir))
        self.set_initial_config(PeriodicTissue.from_config(self._dir,"init_config.txt".format(self._dir)))
        self._initial_config.evaluate_cell_neighbors()
        np.savetxt("{}initial_cost.txt".format(self._dir), [self.evaluate_cost()], fmt='%.2e')

        initial_stresses = {}
        self.calculate_max_shear_stresses()
        for cellID in self._target_cell_to_stress:
            cell = self._config.cells_[cellID]
            initial_stresses[cellID] = cell.max_shear_stress_
        df = pd.DataFrame(list(initial_stresses.items()), columns=['CellID', 'Current'])        
        # df = pd.DataFrame(self._target_cell_to_stress.items(), columns=['cellID', 'target_stress'])
        df.to_csv("{}initial_stress.csv".format(self._dir), index=False)
        # print("Initial Cost: {:.2e}".format(cost))
        self._cost_values = []
        self._q_values = []

    def run_to_max_iters(self,max_iters=2000):
        cost = self.evaluate_cost()
        for _ in range(max_iters):
            self.single_iteration()
            cost = self.evaluate_cost()            
            self._cost_values.append(cost)
            self._config.evaluate_cell_neighbors()
            self._q_values.append(calculate_Q(self._initial_config.cell_neighbors_, self._config.cell_neighbors_))
            print("Iteration: {:d}, Cost: {:.2e}".format(self._iter_counter,cost))
            # Rewrite costs.txt
            np.savetxt("{}costs.txt".format(self._dir), self._cost_values, fmt='%.2e')
            # Rewrite q_values.txt
            np.savetxt("{}q_values.txt".format(self._dir), self._q_values,fmt = "%.4f")
            # Write stresses.csv
            current_stresses = []
            for cellID in self._target_cell_to_stress:
                cell = self._config.cells_[cellID]
                current_stresses.append(cell.max_shear_stress_)
            results = {"CellID": list(self._target_cell_to_stress.keys()),
                    "Target": list(self._target_cell_to_stress.values()),
                    "Current": current_stresses}
            df = pd.DataFrame(results)
            df.to_csv("{}{:04d}.stresses.csv".format(self._dir, self._iter_counter), index=False)
            self._iter_counter += 1
            if cost < self._tolerance:
                break
        print("Optimization finished at iteration: {:d}".format(self._iter_counter-1))
        print("Final cost: {:.2e}".format(cost))
    
    def run(self):
        # #store initial cell parameters
        # os.system("cp {}cellParameters.input {}cellParameters.init.input".format(self._dir,self._dir))
        # os.system("cp {}minimized.txt {}init_config.txt".format(self._dir,self._dir))
        # self._cost_values = []
        cost = self.evaluate_cost()
        # self._cost_values.append(cost)
        # initial_stresses = []
        # self.calculate_max_shear_stresses()
        # for cellID in self._target_cell_to_stress:
        #     cell = self._config.cells_[cellID]
        #     initial_stresses.append(cell.max_shear_stress_)
        
        # # df = pd.DataFrame(self._target_cell_to_stress.items(), columns=['cellID', 'target_stress'])
        # # df.to_csv("{}target_cells.csv".format(self._dir), index=False)
        # print("Initial Cost: {:.2e}".format(cost))
        while cost > self._tolerance:
            # self.set_clamping_tolerance(cost * 1)
            self.single_iteration()
            cost = self.evaluate_cost()
            self._cost_values.append(cost)
            print("Iteration: {:d}, Cost: {:.2e}".format(self._iter_counter,cost))
            # Rewrite costs.txt
            np.savetxt("{}costs.txt".format(self._dir), self._cost_values, fmt='%.2e')
            # Rewrite stresses.csv
            current_stresses = []
            for cellID in self._target_cell_to_stress:
                cell = self._config.cells_[cellID]
                current_stresses.append(cell.max_shear_stress_)
            results = {"CellID": list(self._target_cell_to_stress.keys()),
                    "Target": list(self._target_cell_to_stress.values()),
                    "Current": current_stresses}
            df = pd.DataFrame(results)
            df.to_csv("{}{:04d}.stresses.csv".format(self._dir, self._iter_counter), index=False)
            self._iter_counter += 1
        print("Optimization finished at iteration: {:d}".format(self._iter_counter-1))
        print("Final cost: {:.2e}".format(cost))