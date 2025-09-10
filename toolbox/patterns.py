from toolbox.training import Training
from toolbox import stress
import numpy as np
import os
import random
import pandas as pd

class Patterns(Training):
    def __init__(self):
        super().__init__()
        self._target_cell_to_stress = None
        self._clamping_max_iters = 10 #try 10
        self._clamping_correction_factor = 10 #try 10
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

    def set_central_target_cells(self, n_cells = 1):
        r_lim = 0.9
        self._target_cell_to_stress = {}
        for cellID,cell in self._config.cells_.items():
            if cell.crossBoundary_:
                continue
            r = np.subtract(cell.center_,self._config.sample_center_)
            r = np.linalg.norm(r)
            if r<r_lim:
                self._target_cell_to_stress[cellID] = None
            if len(self._target_cell_to_stress) == n_cells:
                break
        # For checking the above functionality with vtk:
        for polygonID,polygon in self._config.polygons_.items():
            polygon.vtk_scalar_ = 0
        for cellID,_ in self._target_cell_to_stress.items():
            cell = self._config.cells_[cellID]
            for polygonID in cell.polygons_:
                polygon = self._config.polygons_[polygonID]
                polygon.vtk_scalar_ = 1
        self._config.write_periodic_vtk(filename = "target_cells.vtk", use_scalar=True)

    def set_random_target_cells(self, n_cells = 1, average_cells_only = False, exclude_cells = []):
        self._target_cell_to_stress = {}
        for polygonID,polygon in self._config.polygons_.items():
            polygon.vtk_scalar_ = 0
        if average_cells_only:
            for cellID,cell in self._config.cells_.items():
                cell.max_shear_stress_ = stress.calculate_max_shear_stress(self._config, cellID)
            avg_stress = np.mean([cell.max_shear_stress_ for cellID,cell in self._config.cells_.items()])
            std_stress = np.std([cell.max_shear_stress_ for cellID,cell in self._config.cells_.items()])
        while len(self._target_cell_to_stress)<n_cells:
            cellID = random.choice(list(self._config.cells_.keys()))
            cell = self._config.cells_[cellID]
            if cell.crossBoundary_: 
                continue
            if cellID in self._target_cell_to_stress:
                continue
            if average_cells_only:
                lower_limit = avg_stress - std_stress
                if lower_limit < 0:
                    continue
                upper_limit = avg_stress + std_stress
                if (cell.max_shear_stress_ < lower_limit):
                    continue
                if (cell.max_shear_stress_ > upper_limit):
                    continue
            if len(exclude_cells):
                if cellID in exclude_cells:
                    continue
            self._target_cell_to_stress[cellID] = None
            targets_share_polygons = False
            for polygonID in cell.polygons_:
                polygon = self._config.polygons_[polygonID]
                if polygon.vtk_scalar_ == 1:
                    targets_share_polygons = True
                    break    
            if targets_share_polygons:
                continue
            for polygonID in cell.polygons_:
                polygon = self._config.polygons_[polygonID]
                polygon.vtk_scalar_ = 1
            if len(self._target_cell_to_stress) == n_cells:
                break
        # # For checking the above functionality with vtk:
        # for polygonID,polygon in self._config.polygons_.items():
        #     polygon.vtk_scalar_ = 0
        for cellID in self._target_cell_to_stress:
            cell = self._config.cells_[cellID]
            # for polygonID in cell.polygons_:
            #     polygon = self._config.polygons_[polygonID]
            #     polygon.vtk_scalar_ = 1
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
        # The stored cell areas will be used to calculate learning DOF changes
        for cellID,cell in self._config.cells_.items():
            if cellID in self._target_cell_to_stress:
                continue
            free_state_areas[cellID] = cell.surface_area_

        print("\n\n-------------------------------------")
        print("Step 2: CLAMPING")
        self.clamp_target_cells()

        print("\n\n-------------------------------------")
        print("Step 3: Use the clamped state areas to update all learning degrees of freedom.")

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
        for iter in range(self._clamping_max_iters):
            print("Clamping iteration {}".format(iter))
            needs_clamping = []
            for cellID, final_target_stress in self._target_cell_to_stress.items():
                cell = self._config.cells_[cellID]
                current_stress = stress.calculate_max_shear_stress(self._config, cellID)
                if abs(current_stress - final_target_stress) > self._clamping_tolerance:
                    needs_clamping.append(cellID)
            if not len(needs_clamping):
                print("Clamping successful to tolerance")
                return

            for cellID in needs_clamping:
                cell = self._config.cells_[cellID]
                final_target_stress = self._target_cell_to_stress[cellID]
                current_stress = stress.calculate_max_shear_stress(self._config, cellID)
                # The temporary target stress to be solved for should be an overshoot/undershoot of self._target_stress
                temp_target_stress = final_target_stress * (1 + self._clamping_correction_factor * (final_target_stress - current_stress))
                # The temporary target stress should be positive, 
                # so we take the absolute value in case we hit a really small number
                temp_target_stress = abs(temp_target_stress)
                print("Cell: {}, Final Target Stress: {} Temporary Target Stress: {}, Current stress: {}, s0: {}".format(
                    cellID, final_target_stress, temp_target_stress, current_stress, cell.s0_))
                self.solve_cell_s0_for_target_stress(cellID, temp_target_stress)
                current_stress = stress.calculate_max_shear_stress(self._config, cellID)

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
    def initialize(self):
        os.system("cp {}cellParameters.input {}cellParameters.init.input".format(self._dir,self._dir))
        os.system("cp {}minimized.txt {}init_config.txt".format(self._dir,self._dir))
        
        cost = [self.evaluate_cost()]
        np.savetxt("{}initial_cost.txt".format(self._dir), cost, fmt='%.2e')

        initial_stresses = {}
        self.calculate_max_shear_stresses()
        for cellID in self._target_cell_to_stress:
            cell = self._config.cells_[cellID]
            initial_stresses[cellID] = cell.max_shear_stress_
        df = pd.DataFrame(list(initial_stresses.items()), columns=['cellID', 'Stress'])        
        # df = pd.DataFrame(self._target_cell_to_stress.items(), columns=['cellID', 'target_stress'])
        df.to_csv("{}initial_stress.csv".format(self._dir), index=False)
        # print("Initial Cost: {:.2e}".format(cost))
        self._cost_values = []

    def run_to_max_iters(self,max_iters=100):
        cost = self.evaluate_cost()
        for _ in range(max_iters):
            self.set_clamping_tolerance(cost * 1)
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
            if cost <= self._tolerance:
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
            self.set_clamping_tolerance(cost * 1)
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