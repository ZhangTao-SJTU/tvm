from toolbox.periodic import PeriodicTissue
from toolbox.spheroid import Spheroid
from toolbox.training import Training
from toolbox.overlap import calculate_Q
from toolbox import stress
import numpy as np
import os
import glob
import random
import pandas as pd

class Patterns(Training):
    def __init__(self):
        super().__init__()
        self._target_cell_to_stress = None
        self._frozen_cells = []
        self._clamping_max_iters = 10
        self._clamping_s0_lower_limit = 4
        self._clamping_s0_upper_limit = 6
        self._clamping_correction_factor = 0.1
        self._clamping_FIRE_only = True
        self._clamping_tolerance = 1e-15
        self._learning_rate = 10
        self._clear_interval = 100
    @classmethod
    def from_sample(cls,tissue):
        inst = super().from_sample(tissue)
        return inst    

    
    def set_clamping_FIRE_only(self, clamping_FIRE_only):
        self._clamping_FIRE_only = clamping_FIRE_only
    def set_clamping_max_iters(self, clamping_max_iters):
        self._clamping_max_iters = clamping_max_iters
    def set_clamping_correction_factor(self, clamping_correction_factor):
        self._clamping_correction_factor = clamping_correction_factor
    def set_target_cell_to_stress(self,target_cell_to_stress):
        self._target_cell_to_stress = target_cell_to_stress
        self._config.write_cell_collection_vtk(self._target_cell_to_stress,"target_cells.vtk")
    def set_frozen_cells(self, frozen_cells):
        self._frozen_cells = frozen_cells            
    def set_learning_rate(self,learning_rate):
        self._learning_rate = learning_rate
    def set_clear_interval(self, clear_interval):
        self._clear_interval = clear_interval
        
    def calculate_max_shear_stresses(self):
        for cellID in self._target_cell_to_stress:
            self._config.cells_[cellID].max_shear_stress_ = stress.calculate_max_shear_stress(self._config, cellID)

    def evaluate_cost(self):
        self.calculate_max_shear_stresses()
        return np.mean([
            abs(self._config.cells_[cellID].max_shear_stress_ - target_stress)/target_stress 
            for cellID, target_stress in self._target_cell_to_stress.items()])
    
    def initialize(self):
        if not os.path.isfile("{}cellParameters.input".format(self._dir)):
            self.write_cell_parameters()
        self.minimize_config()
        self.write_configuration(filename = "{:07d}.bulk.txt".format(self._iter_counter))
        self.write_cell_parameters(filename = "{:07d}.cellParameters.input".format(self._iter_counter))
        if self._config.tissueType_ == "periodic":
            self.set_initial_config(PeriodicTissue.from_config(self._dir,"{:07d}.bulk.txt".format(self._iter_counter)))
        elif self._config.tissueType_ == "spheroid":
            self.set_initial_config(Spheroid.from_config(self._dir,"{:07d}.bulk.txt".format(self._iter_counter)))
        # np.savetxt("{}initial_cost.txt".format(self._dir), [self.evaluate_cost()], fmt='%.2e')
        self._cost_values = [self.evaluate_cost()]
        self._q_values = [1]
        np.savetxt("{}costs.txt".format(self._dir), self._cost_values, fmt='%.2e')
        np.savetxt("{}q_values.txt".format(self._dir), self._q_values,fmt = "%.4f")
        # Write stresses.csv
        results = {"CellID": list(self._target_cell_to_stress.keys()),
                "Target": list(self._target_cell_to_stress.values()),
                "Current": [self._config.cells_[cellID].max_shear_stress_ for cellID in self._target_cell_to_stress]}
        df = pd.DataFrame(results)
        df.to_csv("{}{:07d}.stresses.csv".format(self._dir, self._iter_counter), index=False)
        self.clear_directory()
        self._iter_counter += 1

    # Binary search for s0 that produces the right stress
    def solve_cell_s0_for_target_stress(self,cellID, target_stress):
        cell = self._config.cells_[cellID]
        def clamping_error(s0):
            cell.s0_ = s0
            return stress.calculate_max_shear_stress(self._config,cellID)-target_stress
        upper_limit = self._clamping_s0_upper_limit
        lower_limit = self._clamping_s0_lower_limit

        root_interval = None
        s0_to_clamping_error = {s0:clamping_error(s0) for s0 in np.linspace(lower_limit, upper_limit,5)}
        for i, s0 in enumerate(list(s0_to_clamping_error.keys())):
            if i == len(s0_to_clamping_error)-1:
                continue
            if np.sign(s0_to_clamping_error[s0]) == np.sign(list(s0_to_clamping_error.values())[i+1]):
                continue
            root_interval = [s0,list(s0_to_clamping_error.keys())[i+1]]
            break
        
        if root_interval is None:
            min_guess = min(s0_to_clamping_error, key=lambda x: abs(s0_to_clamping_error[x]))
            cell.s0_ = min_guess
            print("No root interval found, picking closest guess: {}".format(cell.s0_))
            return
        # Binary search within the root interval
        # while abs(root_interval[1]-root_interval[0])>self._tolerance:
        while abs(root_interval[1]-root_interval[0])>self._clamping_tolerance:
            mid = (root_interval[0]+root_interval[1])/2
            if np.sign(clamping_error(mid)) == np.sign(clamping_error(root_interval[0])):
                root_interval[0] = mid
            else:
                root_interval[1] = mid
        cell.s0_ = (root_interval[0]+root_interval[1])/2
        print("Binary search complete, found s0 = {}".format(cell.s0_))

    # The goal of clamping is to ensure that the target cells reach the final target stress as
    # in self._target_cell_to_stress, upto clamping tolerance.
    # Holding the configuration fixed, we iteratively adjust s0 of the target cells 
    # such that stress overshoots/undershoots the target stress and energy minimize the config at each iteration.    
    # Minimizing changes the surface area and hence the stress, so process is repeated.
    def clamp_target_cells(self):
        # self.write_cell_parameters() ## IS THIS REALLY NEEDED?
        for iter in range(self._clamping_max_iters):
            print("Clamping iteration {}".format(iter))
            # needs_clamping = []
            self.calculate_max_shear_stresses()
            # for cellID, final_target_stress in self._target_cell_to_stress.items():
            #     cell = self._config.cells_[cellID]
            #     if abs(cell.max_shear_stress_ - final_target_stress) > self._clamping_tolerance:
            #         needs_clamping.append(cellID)
            needs_clamping = [cellID for cellID, final_target_stress in self._target_cell_to_stress.items() 
                              if abs(self._config.cells_[cellID].max_shear_stress_ - final_target_stress) > self._clamping_tolerance]
            if not len(needs_clamping):
                print("Clamping successful to clamping tolerance")
                return
            pre_clamping_s0 = pd.read_csv(self._dir+"cellParameters.input",header=None,sep=" ")[2]
            for cellID in needs_clamping:
                final_target_stress = self._target_cell_to_stress[cellID]
                current_stress = stress.calculate_max_shear_stress(self._config, cellID)
                # The temporary target stress to be solved for should be an overshoot/undershoot of self._target_stress
                temp_target_stress = final_target_stress + self._clamping_correction_factor * (final_target_stress - current_stress)
                # The temporary target stress should be positive,
                # so if the overshooting makes it negative, avoid it.
                if temp_target_stress<0:
                    temp_target_stress = final_target_stress
                print("Cell: {}, Final Target Stress: {} Temporary Target Stress: {}, Current stress: {}, s0: {}".format(
                    cellID, final_target_stress, temp_target_stress, current_stress, self._config.cells_[cellID].s0_))
                self.solve_cell_s0_for_target_stress(cellID, temp_target_stress)
            # Write updated cell parameters in preparation for clamped state minimization
            self.write_cell_parameters()
            post_clamping_s0 = pd.read_csv(self._dir+"cellParameters.input",header=None,sep=" ")[2]
            if post_clamping_s0.equals(pre_clamping_s0):
                print("Clamping terminated: doesn't change s0")
                return
            self.minimize_config(FIRE_only = self._clamping_FIRE_only)
            
    
    def single_iteration(self):
        if not len(self._target_cell_to_stress):
            print("No target cells, iteration terminated.")
            return
        # Starting with a minimized config with loaded cell properties...
        print("\n\n=====================================\n\n")  
        print("Starting iteration: {:d}".format(self._iter_counter))
        print("\n\n=====================================\n\n")
        
        print("\n\n-------------------------------------")
        print("Step 1: Evaluate and store the current (free state) areas of hidden (non-target) cells")
        print("-------------------------------------\n\n")
        
        # The stored cell areas will be used to calculate learning DOF changes
        free_state_areas = {cellID:cell.surface_area_ for cellID,cell in self._config.cells_.items() if cellID not in self._target_cell_to_stress}
        print("\n\n-------------------------------------")
        print("Step 2: CLAMPING")
        print("-------------------------------------\n\n")

        self.clamp_target_cells()

        print("\n\n-------------------------------------")
        print("Step 3: Use the clamped state areas to update all learning degrees of freedom.")
        print("-------------------------------------\n\n")

        for cellID in free_state_areas:
            cell = self._config.cells_[cellID]
            if cellID in self._target_cell_to_stress:
                continue
            if cellID in self._frozen_cells:
                continue
            if self._config.tissueType_ == "spheroid" and not cell.type_:
                continue
            del_area = cell.surface_area_ - free_state_areas[cellID]
            s0_change = self._learning_rate * del_area
            ## TESTING
            cell.s0_ -= s0_change
        
        print("\n\n-------------------------------------")
        print("Step 4: UNCLAMP target cells; write unclamped cell parameters, minimize,log")
        print("-------------------------------------\n\n")
        for cellID in self._target_cell_to_stress:
            self._config.cells_[cellID].s0_ = self._config.s0_
        self.write_cell_parameters()
        self.minimize_config()
        # Logging
        self.write_cell_parameters(filename = "{:07d}.cellParameters.input".format(self._iter_counter))
        self.write_configuration(filename = "{:07d}.bulk.txt".format(self._iter_counter))

    def run_to_max_iters(self,max_iters=2000):
        cost = self.evaluate_cost()
        for _ in range(max_iters):
            if cost < self._tolerance:
                print("Iteration {:d} not commenced: cost < tolerance".format(self._iter_counter))
                break
            self.single_iteration()
            cost = self.evaluate_cost()
            self._config.evaluate_cell_neighbors()
            q_value = calculate_Q(self._initial_config.cell_neighbors_, self._config.cell_neighbors_)           
            self._cost_values.append(cost)
            self._q_values.append(q_value)
            print("Iteration: {:d}, Cost: {:.2e} Q2: {:.4f}".format(self._iter_counter,cost,q_value))
            np.savetxt("{}costs.txt".format(self._dir), self._cost_values, fmt='%.2e')
            np.savetxt("{}q_values.txt".format(self._dir), self._q_values,fmt = "%.4f")
            results = {"CellID": list(self._target_cell_to_stress.keys()),
                    "Target": list(self._target_cell_to_stress.values()),
                    "Current": [self._config.cells_[cellID].max_shear_stress_ for cellID in self._target_cell_to_stress]}
            df = pd.DataFrame(results)
            df.to_csv("{}{:07d}.stresses.csv".format(self._dir, self._iter_counter), index=False)
            if self._iter_counter % self._clear_interval == 0:
                self.clear_directory()
            self._iter_counter += 1
        print("Final cost: {:.2e}".format(cost))
    
    def clear_directory(self):
        os.makedirs("{}files".format(self._dir), exist_ok=True)
        for filename in ["bulk.txt","cellParameters.input","stresses.csv"]:
            os.system("cp {}{:07d}.{} {}files/".format(self._dir,self._iter_counter,filename,self._dir))
            # remove all other files except the latest one
            for file in glob.glob("{}*.{}".format(self._dir,filename)):
                os.remove(file)

    def set_random_target_cells(self, n_cells = 1, target_stress = 1, **kwargs):
        stress_limits = []
        exclude_cells = []
        if "stress_limits" in kwargs:
            stress_limits = kwargs["stress_limits"]
        if "exclude_cells" in kwargs:
            exclude_cells = kwargs["exclude_cells"]
        for polygonID,polygon in self._config.polygons_.items():
            polygon.vtk_scalar_ = 0

        target_cell_to_stress = {}
        while len(target_cell_to_stress)<n_cells:
            cellID = random.choice(list(self._config.cells_.keys()))
            cell = self._config.cells_[cellID]
            if cell.crossBoundary_: 
                continue
            if self._config.tissueType_ == "spheroid" and cell.is_surface_:
                continue
            if self._config.tissueType_ == "spheroid" and not cell.type_:
                continue
            if cellID in target_cell_to_stress:
                continue
            if len(stress_limits):
                cell.max_shear_stress_ = stress.calculate_max_shear_stress(self._config,cellID)
                if (cell.max_shear_stress_ < stress_limits[0]):
                    continue
                if (cell.max_shear_stress_ > stress_limits[1]):
                    continue
            if len(exclude_cells) and cellID in exclude_cells:
                    continue
            target_cell_to_stress[cellID] = target_stress
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
            if len(target_cell_to_stress) == n_cells:
                break
        self.set_target_cell_to_stress(target_cell_to_stress)
        self._config.write_cell_collection_vtk(list(target_cell_to_stress.keys()),"target_cells_isolated.vtk",use_scalar=False)
